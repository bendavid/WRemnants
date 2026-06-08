#!/usr/bin/env python3
"""Propagate the stage-1 FLOW (template) statistical uncertainty into the
binned-θ covariance — Hessian-free, exact column-by-column (no sketching).

Two-step (plug-in) M-estimation sandwich: stage 2 solves ∂L₂/∂w(ŵ, φ̂) = 0
with the flow weights φ̂ plugged in from stage 1 (∂L₁/∂φ(φ̂) = 0), so

    Cov_flow(ŵ) = H_ww⁻¹ · H_wφ · H₁⁻¹ · H_φw · H_ww⁻¹           (MLE form)

with w = (active θ_scale | active θ_smear | background MLP) — the background
is PROFILED inside H_ww, exactly as in the data-statistics covariance — and
H₁ the stage-1 Hessian over the flow weights. Everything is evaluated
matrix-free:

  • only the n_θ θ-columns of the chain are needed (the flow uncertainty
    enters θ through a rank-≤ n_θ map), so the cost is n_θ CG solves, never
    an n_φ² object;
  • H·v by double-backward (∇⟨∇L, v⟩) accumulated over in-memory event
    chunks; the mixed block H_φw·x the same way across the two parameter
    sets;
  • CG with ridge damping (H + λI) — λ is the regularisation of the
    near-null (early-stopping-flat) flow directions; scan it for a plateau.

The stage-2 bread H_ww and the stage-1 meat H₁ play DIFFERENT roles, so they
are treated differently (``--w-solver``):

  • H₁ (the central PSD factor / Cov(φ̂)) MUST be PD — its negative
    eigenvalues would make the covariance non-PSD — so it is ridged (CG).
  • H_ww enters only as the congruence bread M = H_ww⁻¹H_wφ; the θθ block
    Cov_θθ = Mᵀ H₁⁻¹ M is PSD for ANY invertible H_ww, so H_ww need NOT be
    PD. The OBSERVED H_ww is in fact genuinely indefinite (the nonlinear
    background MLP is non-convex even at a minimum), which a uniform ridge
    cannot fix. Default ``--w-solver fisher`` uses the empirical Fisher
    J_ww = Σ wᵢsᵢsᵢᵀ (PSD by construction, the SAME bread the data-stat
    covariance uses → coherent addition) via one dense Cholesky solve for all
    columns (no CG_w). ``observed`` keeps the double-backward Hessian (CG +
    ridge escalation); ``minres`` inverts the indefinite observed H_ww
    directly (valid by the congruence above) as a cross-check.

Per column j (a unit vector in the θ-active block):
    x_j = (H_ww + λ_w)⁻¹ e_j          (CG, stage-2 HVPs)
    u_j = H_φw x_j                     (one mixed double-backward)
    z_j = (H₁ + λ_φ)⁻¹ u_j            (CG, stage-1 HVPs)
    Cov_flow[j, k] = u_jᵀ z_k / α₁

α₁ = (full flow-sample Σw)/(seen Σw) corrects an event-capped stage-1 sample;
the stage-2 scale factors cancel exactly in the chain (H_ww⁻¹·H_wφ·…·H_ww⁻¹
carries net power 0 of the stage-2 normalisation), so only α₁ survives.

Scope: binned θ AND θ-mlp. For ``--theta-mlp`` the quantities of interest are
the net OUTPUTS T_j on a fixed (η-bin-centre × φ-bin-centre) grid — physical
(A, e, M) and O(1) effective (a, c), the SAME active table layout and units as
the output-space Fisher (``--output-fisher``) — propagated by the delta
method: Cov(T) = J_T·Cov(ŵ)·J_Tᵀ with w = (θ-net ⊕ background-MLP weights),
so the only change to the chain is the RHS column e_j → g_j = ∂T_j/∂w (one
cheap autograd row each). Exact, column per grid output, no sketching — the
cost scales linearly in n_η·n_φ·n_comp. The θ-net weight space is
over-parameterised (near-null weight directions), so ``--ridge-w`` is more
consequential there: scan it for a plateau like ``--ridge-flow``.

Flow archs with a plain NLL stage-1 objective (compact/dcb/ege: −Σw·log p₀ on
the FLOW window) and gf/nsf (truncated window-norm + gauge penalty,
replicating step1); nce's BCE objective is not supported.

Output: ``flow_uncertainty.pt``. Binned: the raw-θ covariance in the SAME
active layout as ``empirical_fisher.pt`` (combine by addition), the
physical-units extras via the shared ``_theta_cov_extras`` (softplus
delta-method included), and CG diagnostics; with ``--fisher
<empirical_fisher.pt>`` the combined (data + flow) covariance and a
per-parameter σ budget table are also written. θ-mlp: the grid-output
covariance plus the φ-collapsed per-η keys mirroring the output-space Fisher
file (``*_flow`` suffix); with ``--fisher`` (an ``--output-fisher`` file with
the same n_phi) the blockwise totals (``*_total``) and a per-component σ
budget summary.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from jpsi_mass_arrow_loader import JpsiMassArrowLoader  # noqa: E402
from jpsi_mass_fit_diagnostics import load_model_from_checkpoint  # noqa: E402
import train_jpsi_mass_fit as T  # noqa: E402
from train_jpsi_mass_fit import _move_batch, discover_shards  # noqa: E402


# ---------------------------------------------------------------------------
# Parameter packing: w = [θ_scale (all) | θ_smear (active cols) | MLP (all)]
# — the SAME active layout as compute_empirical_fisher_joint, so the saved
# covariance combines with empirical_fisher.pt by addition.
# ---------------------------------------------------------------------------


def _w_layout(model):
    """Return (params, slices, labels, n_theta_active, smear_cols): the
    stage-2 parameter tensors, the per-tensor (offset, active-flat-index)
    mapping into the packed active vector, and the θ-active count."""
    params, segs, labels = [], [], []
    smear_cols = T._smear_active_cols(model) if model.smearing_enabled else []
    off = 0
    if model.scale_enabled:
        n_eta, nc = model.theta_scale.shape
        idx = torch.arange(model.theta_scale.numel())
        params.append(model.theta_scale)
        segs.append((off, idx))
        comp = ("A", "e", "M")
        labels += [f"{comp[c]}[{b}](raw)" for b in range(n_eta) for c in range(nc)]
        off += idx.numel()
    if model.smearing_enabled and smear_cols:
        n_eta, n_comp = model.theta_smear.shape
        idx = torch.tensor([b * n_comp + c for b in range(n_eta)
                            for c in smear_cols], dtype=torch.long)
        params.append(model.theta_smear)
        segs.append((off, idx))
        comp = ("a", "c")
        labels += [f"{comp[c]}[{b}](raw)" for b in range(n_eta) for c in smear_cols]
        off += idx.numel()
    n_theta_active = off
    if model.background_enabled:
        for i, p in enumerate(model.mlp.parameters()):
            idx = torch.arange(p.numel())
            params.append(p)
            segs.append((off, idx))
            labels += [f"mlp[{i}][{k}]" for k in range(p.numel())]
            off += idx.numel()
    return params, segs, labels, n_theta_active, smear_cols


def _w_layout_mlp(model):
    """θ-mlp w-layout: w = [θ-net (all) | background MLP (all)]. The grid
    outputs T(η, φ) are FUNCTIONS of w, so every net weight is active and the
    chain's RHS columns are the output gradients ∂T_j/∂w (delta method)
    instead of unit vectors. Returns (params, segs, n_net)."""
    params, segs = [], []
    off = 0
    for p in model.theta_net.parameters():
        idx = torch.arange(p.numel())
        params.append(p)
        segs.append((off, idx))
        off += idx.numel()
    n_net = off
    if model.background_enabled:
        for p in model.mlp.parameters():
            idx = torch.arange(p.numel())
            params.append(p)
            segs.append((off, idx))
            off += idx.numel()
    return params, segs, n_net


def _grid_outputs(model, eta_edges, n_phi_g, scale_cols, smear_cols, device):
    """Evaluate the θ-net on the fixed (η-bin-centre × φ-bin-centre) grid and
    return the ACTIVE output vector (autograd graph attached) in the SAME
    active ordering and units as ``compute_output_fisher_2d`` / the
    ``--output-fisher`` file: scale block then smear block, each (η outer, φ
    middle, active-col inner); scale = PHYSICAL (A, e, M), smear = O(1)
    EFFECTIVE (a, c) so the positivity-reparam delta method is automatic
    through autograd. Returns (o_active, labels, eta_centres, phi_centres)."""
    mdt = next(model.theta_net.parameters()).dtype
    centres = 0.5 * (np.asarray(eta_edges[:-1]) + np.asarray(eta_edges[1:]))
    n_eta = int(len(centres))
    eta_c = torch.as_tensor(centres, dtype=mdt, device=device)
    phi_edges = torch.linspace(-float(np.pi), float(np.pi), n_phi_g + 1,
                               dtype=mdt, device=device)
    phi_c = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    eg = eta_c[:, None, None].expand(n_eta, n_phi_g, 2).reshape(-1, 2)
    pg = phi_c[None, :, None].expand(n_eta, n_phi_g, 2).reshape(-1, 2)
    AeM, ac = model.theta_net(eg, pg)                    # physical A,e,M; raw a,c
    AeM = AeM[:, 0, :]                                   # [n_eta*n_phi_g, 3]
    ac_eff = model._smear_raw_to_effective(ac[:, 0, :])  # [n_eta*n_phi_g, 2]
    rows, labels = [], []
    comp_s, comp_c = ("A", "e", "M"), ("a", "c")
    for ie in range(n_eta):
        for ip in range(n_phi_g):
            for c in scale_cols:
                rows.append(AeM[ie * n_phi_g + ip, c])
                labels.append(f"{comp_s[c]}[{ie},{ip}]")
    for ie in range(n_eta):
        for ip in range(n_phi_g):
            for c in smear_cols:
                rows.append(ac_eff[ie * n_phi_g + ip, c])
                labels.append(f"{comp_c[c]}[{ie},{ip}]")
    return (torch.stack(rows), labels, centres,
            phi_c.detach().cpu().numpy())


def _pack(grads, params, segs, n_w, device, dtype=torch.float64):
    """Per-tensor gradients → packed active vector [n_w] (float64)."""
    v = torch.zeros(n_w, device=device, dtype=dtype)
    for g, p, (off, idx) in zip(grads, params, segs):
        if g is None:
            continue
        v[off:off + idx.numel()] = g.reshape(-1)[idx.to(g.device)].to(dtype)
    return v


def _unpack(v, params, segs):
    """Packed active vector → list of param-shaped direction tensors (zeros
    at inactive entries), in each param's dtype/device."""
    outs = []
    for p, (off, idx) in zip(params, segs):
        d = torch.zeros(p.numel(), device=p.device, dtype=p.dtype)
        d[idx.to(p.device)] = v[off:off + idx.numel()].to(p.dtype)
        outs.append(d.view_as(p))
    return outs


# ---------------------------------------------------------------------------
# Chunked losses + HVPs
# ---------------------------------------------------------------------------


class ChunkedLoss:
    """SUM-form weighted NLL over materialised event tensors, evaluated in
    chunks; provides H·v (double backward) and the mixed H_other,this·x."""

    def __init__(self, loss_chunk_fn, n_events, chunk):
        self.loss_chunk_fn = loss_chunk_fn   # (i0, i1) -> scalar loss (graph)
        self.n = int(n_events)
        self.chunk = int(chunk)

    def hvp(self, v_list, params):
        """Σ_chunks ∇⟨∇L_chunk, v⟩ w.r.t. ``params`` (list of grads)."""
        acc = [torch.zeros_like(p, dtype=torch.float64) for p in params]
        for i0 in range(0, self.n, self.chunk):
            L = self.loss_chunk_fn(i0, min(i0 + self.chunk, self.n))
            g = torch.autograd.grad(L, params, create_graph=True,
                                    allow_unused=True)
            s = sum((gi * vi).sum() for gi, vi in zip(g, v_list)
                    if gi is not None)
            h = torch.autograd.grad(s, params, retain_graph=False,
                                    allow_unused=True)
            for a, hi in zip(acc, h):
                if hi is not None:
                    a += hi.detach().double()
        return acc

    def mixed(self, x_list, params_in, params_out):
        """Σ_chunks ∇_{params_out} ⟨∇_{params_in} L_chunk, x⟩."""
        acc = [torch.zeros_like(p, dtype=torch.float64) for p in params_out]
        for i0 in range(0, self.n, self.chunk):
            L = self.loss_chunk_fn(i0, min(i0 + self.chunk, self.n))
            g = torch.autograd.grad(L, params_in, create_graph=True,
                                    allow_unused=True)
            s = sum((gi * xi).sum() for gi, xi in zip(g, x_list)
                    if gi is not None)
            h = torch.autograd.grad(s, params_out, retain_graph=False,
                                    allow_unused=True)
            for a, hi in zip(acc, h):
                if hi is not None:
                    a += hi.detach().double()
        return acc

    def mixed_batch(self, X, params_in, segs_in, n_in, params_out, n_out, dev):
        """BATCHED mixed block: U = H_{out,in} X for X [n_in (active), K] →
        U [n_out, K], sharing the per-chunk forward + the single ∇_in L
        backward across all K columns (only the ∇_out vjp is K-fold, via
        is_grads_batched). Same result as K separate ``mixed`` calls."""
        K = X.shape[1]
        acc = torch.zeros(n_out, K, dtype=torch.float64, device=dev)
        eyeK = torch.eye(K, device=dev)
        for i0 in range(0, self.n, self.chunk):
            L = self.loss_chunk_fn(i0, min(i0 + self.chunk, self.n))
            g = torch.autograd.grad(L, params_in, create_graph=True,
                                    allow_unused=True)
            # inner_k = ⟨∇_in L, x_k⟩, assembled per active param block.
            inner = X.new_zeros(K)
            for gi, p, (off, idx) in zip(g, params_in, segs_in):
                if gi is None:
                    continue
                inner = inner + (gi.reshape(-1)[idx.to(gi.device)]
                                 @ X[off:off + idx.numel()])
            h = torch.autograd.grad(inner, params_out,
                                    grad_outputs=eyeK.to(inner.dtype),
                                    is_grads_batched=True, allow_unused=True)
            off = 0
            for hi, p in zip(h, params_out):
                if hi is not None:
                    acc[off:off + p.numel()] += hi.reshape(K, -1).t().double()
                off += p.numel()
        return acc


def _cg(apply_A, b, tol, max_iter, label="", progress=True, report_dt=20.0):
    """Standard CG on the (damped) SPD system; returns (x, iters, rel_res).
    Aborts with a clear message on negative curvature (raise the ridge).
    With ``progress`` prints a throttled per-iteration line (every
    ``report_dt`` seconds) — each CG iteration is one full HVP over all the
    chunked events, so a column can take minutes; this shows it is alive."""
    x = torch.zeros_like(b)
    r = b.clone()
    p = r.clone()
    rs = float(r @ r)
    b_norm = max(float(b.norm()), 1e-300)
    it = 0
    t_last = time.time()
    while it < max_iter and (rs ** 0.5) / b_norm > tol:
        Ap = apply_A(p)
        pAp = float(p @ Ap)
        if pAp <= 0.0:
            raise RuntimeError(
                f"CG[{label}]: non-positive curvature (pᵀAp = {pAp:.3e}) at "
                f"iteration {it} — the (damped) Hessian is not PSD here. "
                f"For a 'w:' solve the usual cause is too small "
                f"--max-events-fit: the observed H_ww (incl. the nonlinear "
                f"background block) is indefinite when estimated from few "
                f"events at the full-fit optimum — raise --max-events-fit "
                f"(≫ n_w) before reaching for --ridge-w. For a 'φ:' solve "
                f"raise --ridge-flow (the flow's flat early-stopping "
                f"directions need regularising).")
        alpha = rs / pAp
        x += alpha * p
        r -= alpha * Ap
        rs_new = float(r @ r)
        p = r + (rs_new / rs) * p
        rs = rs_new
        it += 1
        if progress and (time.time() - t_last) > report_dt:
            print(f"        CG[{label}] it {it}/{max_iter}  "
                  f"rel {(rs ** 0.5) / b_norm:.2e} (tol {tol:g})", flush=True)
            t_last = time.time()
    return x, it, (rs ** 0.5) / b_norm


def _minres(apply_A, b, tol, max_iter, label="", progress=True, report_dt=20.0):
    """MINRES for a SYMMETRIC (possibly INDEFINITE) but invertible operator —
    minimises ‖Ax−b‖ over the Krylov subspace, so it does NOT abort on
    negative curvature (unlike CG). Used for the observed H_ww 'minres' solver,
    where the bread is indefinite (nonlinear background block) but the
    propagated θθ covariance M_θᵀ H₁⁻¹ M_θ is PSD by congruence regardless.
    Canonical Paige–Saunders recurrence (matching scipy.sparse.linalg.minres,
    no preconditioner); returns (x, iters, rel_res)."""
    x = torch.zeros_like(b)
    beta1 = float(b.norm())
    b_norm = max(beta1, 1e-300)
    if beta1 == 0.0:
        return x, 0, 0.0
    oldb = 0.0
    beta = beta1
    dbar = 0.0
    epsln = 0.0
    phibar = beta1
    cs = -1.0
    sn = 0.0
    w = torch.zeros_like(b)
    w2 = torch.zeros_like(b)
    r1 = b.clone()
    r2 = b.clone()
    y = b.clone()
    it = 0
    t_last = time.time()
    while it < max_iter and phibar / b_norm > tol:
        it += 1
        s = 1.0 / beta
        v = s * y
        y = apply_A(v)
        if it >= 2:
            y = y - (beta / oldb) * r1
        alfa = float(v @ y)
        y = y - (alfa / beta) * r2
        r1 = r2
        r2 = y
        oldb = beta
        beta = float(r2.norm())
        # apply previous rotation and compute the new one
        oldeps = epsln
        delta = cs * dbar + sn * alfa
        gbar = sn * dbar - cs * alfa
        epsln = sn * beta
        dbar = -cs * beta
        gamma = max((gbar * gbar + beta * beta) ** 0.5, 1e-300)
        cs = gbar / gamma
        sn = beta / gamma
        phi = cs * phibar
        phibar = sn * phibar
        denom = 1.0 / gamma
        w1 = w2
        w2 = w
        w = (v - oldeps * w1 - delta * w2) * denom
        x = x + phi * w
        if progress and (time.time() - t_last) > report_dt:
            print(f"        MINRES[{label}] it {it}/{max_iter}  "
                  f"rel {phibar / b_norm:.2e} (tol {tol:g})", flush=True)
            t_last = time.time()
    res = float((apply_A(x) - b).norm()) / b_norm   # true residual
    return x, it, res


# ---------------------------------------------------------------------------
# Event materialisation
# ---------------------------------------------------------------------------

_KEYS = ("mll", "pt_pm", "eta_pm", "phi_pm", "q_pm", "b_pm", "cond_std", "w")


def _materialise(loader, device, dtype, take_data_branch, mc_as_data,
                 cap, label):
    """Stream the loader once; keep the branch rows (concatenated, up to
    ``cap`` events) on ``device`` and keep counting Σw past the cap for the
    rescale factor. Returns (tensors dict, w_seen, w_total)."""
    bufs = {k: [] for k in _KEYS}
    w_seen = 0.0
    w_total = 0.0
    n_kept = 0
    for batch in loader:
        if take_data_branch:
            sel = (~batch["is_data_mask"]) if mc_as_data else batch["is_data_mask"]
        else:
            sel = ~batch["is_data_mask"]
        sel = sel & (batch["w"] > 0)
        if not bool(sel.any()):
            continue
        w_b = float(batch["w"][sel].sum())
        w_total += w_b
        if cap and n_kept >= cap:
            continue
        idx = sel.nonzero(as_tuple=True)[0]
        if cap and n_kept + idx.numel() > cap:
            idx = idx[: cap - n_kept]
        for k in _KEYS:
            t = batch[k][idx]
            if t.is_floating_point():
                t = t.to(device, dtype)
            else:
                t = t.to(device)
            bufs[k].append(t)
        n_kept += int(idx.numel())
        w_seen += float(batch["w"][idx].sum())
    if n_kept == 0:
        raise RuntimeError(f"no events materialised for the {label} sample")
    out = {k: torch.cat(v) for k, v in bufs.items()}
    print(f"  {label}: {n_kept} events on {device} "
          f"(Σw seen {w_seen:.4g} / total {w_total:.4g}; "
          f"rescale α = {w_total / max(w_seen, 1e-300):.4f})")
    return out, w_seen, w_total


# ---------------------------------------------------------------------------
# θ-mlp grid post-processing (mirrors _run_output_fisher_mlp's saved keys)
# ---------------------------------------------------------------------------


def _phi_collapse_np(block, aw, n_eta, n_phi_g, cols):
    """Occupancy-weighted φ-collapse of an active-table covariance block:
    cov_φmean[η,c,η',c'] = Σ_{φ,φ'} a[η,φ] a[η',φ'] B[η,φ,c,η',φ',c']."""
    nc = len(cols)
    B = block.reshape(n_eta, n_phi_g, nc, n_eta, n_phi_g, nc)
    return np.einsum("ep,EP,epcEPC->ecEC", aw, aw, B)


def _mlp_grid_keys(C, n_eta, n_phi_g, scale_cols, smear_cols, aw, suffix=""):
    """Build the φ-collapsed per-η covariance keys + the full 2-D blocks from
    the active-table covariance ``C`` [n_act, n_act] (np.float64), in the SAME
    shapes/units as the ``--output-fisher`` empirical_fisher.pt (physical
    A,e,M; smear converted O(1) effective → PHYSICAL via SMEAR_VAR_SCALE)."""
    from jpsi_mass_model import SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C
    out = {}
    n_sa = n_eta * n_phi_g * len(scale_cols)
    if scale_cols:
        Cs = _phi_collapse_np(C[:n_sa, :n_sa], aw, n_eta, n_phi_g, scale_cols)
        cov_s = np.zeros((n_eta, 3, n_eta, 3), dtype=np.float64)
        for i, ci in enumerate(scale_cols):
            for j, cj in enumerate(scale_cols):
                cov_s[:, ci, :, cj] = Cs[:, i, :, j]
        cs2 = cov_s.reshape(n_eta * 3, n_eta * 3)     # symmetrise ULP asymmetry
        cov_s = (0.5 * (cs2 + cs2.T)).reshape(n_eta, 3, n_eta, 3)
        out[f"covariance_24_3_24_3{suffix}"] = torch.tensor(
            cov_s, dtype=torch.float32)
        out[f"sigma_scale_24_3{suffix}"] = torch.sqrt(torch.clamp(
            torch.tensor(np.einsum("icic->ic", cov_s)), min=0.0)).float()
        Cs2d = C[:n_sa, :n_sa].reshape(n_eta, n_phi_g, len(scale_cols),
                                       n_eta, n_phi_g, len(scale_cols))
        full_s = np.zeros((n_eta, n_phi_g, 3, n_eta, n_phi_g, 3),
                          dtype=np.float64)
        for i, ci in enumerate(scale_cols):
            for j, cj in enumerate(scale_cols):
                full_s[:, :, ci, :, :, cj] = Cs2d[:, :, i, :, :, j]
        out[f"covariance_scale_2d{suffix}"] = torch.tensor(
            full_s, dtype=torch.float32)
    if smear_cols:
        sv = [SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C]
        Cc = _phi_collapse_np(C[n_sa:, n_sa:], aw, n_eta, n_phi_g, smear_cols)
        cov_c = np.zeros((n_eta, 2, n_eta, 2), dtype=np.float64)
        for i, ci in enumerate(smear_cols):
            for j, cj in enumerate(smear_cols):
                cov_c[:, ci, :, cj] = Cc[:, i, :, j] * sv[ci] * sv[cj]
        cc2 = cov_c.reshape(n_eta * 2, n_eta * 2)     # symmetrise ULP asymmetry
        cov_c = (0.5 * (cc2 + cc2.T)).reshape(n_eta, 2, n_eta, 2)
        out[f"covariance_smear_24_2_24_2{suffix}"] = torch.tensor(
            cov_c, dtype=torch.float32)
        out[f"sigma_smear_eff_24_2{suffix}"] = torch.sqrt(torch.clamp(
            torch.tensor(np.einsum("icic->ic", cov_c)), min=0.0)).float()
        Cc2d = C[n_sa:, n_sa:].reshape(n_eta, n_phi_g, len(smear_cols),
                                       n_eta, n_phi_g, len(smear_cols))
        full_c = np.zeros((n_eta, n_phi_g, 2, n_eta, n_phi_g, 2),
                          dtype=np.float64)
        for i, ci in enumerate(smear_cols):
            for j, cj in enumerate(smear_cols):
                full_c[:, :, ci, :, :, cj] = Cc2d[:, :, i, :, :, j] \
                    * sv[ci] * sv[cj]
        out[f"covariance_smear_2d{suffix}"] = torch.tensor(
            full_c, dtype=torch.float32)
    return out


def _grid_phi_weights(fit_ev, n_eta, n_phi_g):
    """Per-(η,φ)-cell Σw occupancy of the fit sample (both muons) → the
    normalised per-η φ-weights for the collapse (uniform for empty η rows).
    Returns (aw [n_eta, n_phi_g], cell_w [n_eta, n_phi_g])."""
    bnd = torch.linspace(-float(np.pi), float(np.pi), n_phi_g + 1,
                         device=fit_ev["phi_pm"].device,
                         dtype=fit_ev["phi_pm"].dtype)[1:-1].contiguous()
    pb = torch.clamp(torch.bucketize(fit_ev["phi_pm"], bnd), 0, n_phi_g - 1)
    eb = fit_ev["b_pm"].long()
    flat = (eb * n_phi_g + pb).reshape(-1)
    cw = torch.zeros(n_eta * n_phi_g, dtype=torch.float64,
                     device=flat.device)
    cw.index_add_(0, flat,
                  fit_ev["w"].double().unsqueeze(1).expand(-1, 2).reshape(-1))
    cell_w = cw.view(n_eta, n_phi_g).cpu().numpy()
    aw = cell_w.copy()
    rs = aw.sum(axis=1, keepdims=True)
    unif = (rs <= 0).reshape(-1)
    aw[unif] = 1.0
    rs[unif.reshape(-1, 1)] = n_phi_g
    return aw / rs, cell_w


# ---------------------------------------------------------------------------
# Empirical Fisher (PSD bread) for H_ww
# ---------------------------------------------------------------------------


def _empirical_fisher_dense(per_event_fn, n_fit, score_chunk, params, segs,
                            n_w, dev, progress=True):
    """Dense empirical Fisher ``J_ww = Σ_i w_i s_i s_iᵀ`` over the active
    w-space, with per-event scores ``s_i = ∂(per-event NLL)/∂w`` (batched vjp,
    per-event-loop fallback). PSD by construction → a valid, well-conditioned
    bread even where the OBSERVED H_ww is indefinite (nonlinear background).
    Same J the data-stat covariance uses, so the two combine coherently.
    Returns J [n_w, n_w] (fp64, on ``dev``)."""
    J = torch.zeros((n_w, n_w), dtype=torch.float64, device=dev)
    batched = True
    t0 = time.time()
    n_chunks = (n_fit + score_chunk - 1) // score_chunk
    for ci, i0 in enumerate(range(0, n_fit, score_chunk)):
        i1 = min(i0 + score_chunk, n_fit)
        per, w = per_event_fn(i0, i1)            # per [c] (graph), w [c]
        c = int(per.shape[0])
        S = None
        if batched:
            try:
                eye = torch.eye(c, device=dev, dtype=per.dtype)
                g = torch.autograd.grad(per, params, grad_outputs=eye,
                                        is_grads_batched=True,
                                        retain_graph=False, allow_unused=True)
                S = torch.zeros((c, n_w), dtype=torch.float64, device=dev)
                for gi, p, (off, idx) in zip(g, params, segs):
                    if gi is None:
                        continue
                    S[:, off:off + idx.numel()] = \
                        gi.reshape(c, -1)[:, idx.to(gi.device)].double()
            except (RuntimeError, NotImplementedError) as e:
                batched = False
                print(f"  note: batched per-event score unavailable "
                      f"({type(e).__name__}: {str(e).splitlines()[0][:70]}); "
                      f"per-event loop", flush=True)
        if S is None:
            rows = []
            for jj in range(c):
                gj = torch.autograd.grad(per[jj], params,
                                         retain_graph=(jj < c - 1),
                                         allow_unused=True)
                rows.append(_pack(gj, params, segs, n_w, dev))
            S = torch.stack(rows)
        J += (S * w.double().unsqueeze(1)).t() @ S
        if progress and (ci % 50 == 0 or i1 == n_fit):
            print(f"  fisher build: {i1}/{n_fit} events "
                  f"({time.time() - t0:.0f}s, chunk {ci + 1}/{n_chunks})",
                  flush=True)
    return 0.5 * (J + J.T)


def _efvp(per_event_fn, model_params, n_events, chunk, v_list):
    """MATRIX-FREE empirical-Fisher-vector product J·v = Σ_e w_e g_e (g_eᵀv),
    g_e = ∂(per-event NLL)/∂θ, accumulated over event chunks — for the flow
    meat H₁ (58208-dim, too large to densify). PSD by construction, so CG on it
    never hits negative curvature (no escalation, only a small null-space ridge).

    Computed REVERSE-mode only (no forward-mode AD rules needed for erf/erfcx),
    via the double-VJP identity for the per-event directional derivative:
      A(u) = ∂⟨per, u⟩/∂θ = Σ_e u_e g_e         (1st backward, create_graph)
      s    = ∂⟨A(u), v⟩/∂u = (g_eᵀv)_e          (2nd backward, in event space)
      J·v  = ∂⟨per, w⊙s⟩/∂θ = Σ_e w_e s_e g_e   (3rd backward)
    One forward + three backwards per chunk (~1.5× an HVP). ``per_event_fn``
    returns ``(per [c] with graph, w [c])``; returns the param-shaped grads."""
    acc = [torch.zeros_like(p, dtype=torch.float64) for p in model_params]
    for i0 in range(0, n_events, chunk):
        per, w = per_event_fn(i0, min(i0 + chunk, n_events))
        c = int(per.shape[0])
        u = torch.zeros(c, device=per.device, dtype=per.dtype,
                        requires_grad=True)
        Au = torch.autograd.grad((per * u).sum(), model_params,
                                 create_graph=True, allow_unused=True)
        inner = sum((a * vi).sum() for a, vi in zip(Au, v_list)
                    if a is not None)
        s = torch.autograd.grad(inner, u, retain_graph=True)[0]      # [c]
        ws = (w * s).detach()
        g = torch.autograd.grad((ws * per).sum(), model_params,
                                allow_unused=True)
        for a, gi in zip(acc, g):
            if gi is not None:
                a += gi.detach().double()
    return acc


def _efvp_batch(per_event_fn, model_params, n_events, chunk, V, offs):
    """BATCHED matrix-free empirical-Fisher-vector product: J·V for V
    [n_φ, K] → [n_φ, K], sharing the per-chunk forward AND the
    column-independent A(u) backward across all K columns (only the two
    v-dependent vjps are K-fold, via is_grads_batched). Same as K separate
    ``_efvp`` calls. ``offs`` = [(offset, numel)] per param (flat layout)."""
    K = V.shape[1]
    dev = V.device
    acc = torch.zeros(V.shape[0], K, dtype=torch.float64, device=dev)
    eyeK = torch.eye(K, device=dev)
    for i0 in range(0, n_events, chunk):
        per, w = per_event_fn(i0, min(i0 + chunk, n_events))
        c = int(per.shape[0])
        u = torch.zeros(c, device=per.device, dtype=per.dtype,
                        requires_grad=True)
        Au = torch.autograd.grad((per * u).sum(), model_params,
                                 create_graph=True, allow_unused=True)
        Au_flat = torch.cat([
            (a.reshape(-1) if a is not None
             else torch.zeros(num, device=dev, dtype=per.dtype))
            for a, (o, num) in zip(Au, offs)])                      # [n_φ]
        inner = Au_flat @ V.to(per.dtype)                           # [K]
        s = torch.autograd.grad(inner, u, grad_outputs=eyeK.to(per.dtype),
                                is_grads_batched=True, retain_graph=True)[0]  # [K,c]
        obj = ((w.unsqueeze(0) * s).detach() * per.unsqueeze(0)).sum(1)   # [K]
        g = torch.autograd.grad(obj, model_params,
                                grad_outputs=eyeK.to(per.dtype),
                                is_grads_batched=True, allow_unused=True)
        for gi, (o, num) in zip(g, offs):
            if gi is not None:
                acc[o:o + num] += gi.reshape(K, -1).t().double()
    return acc


def _fisher_diag(per_event_fn, model_params, n_events, chunk, n_param, dev,
                 progress=True):
    """Exact diagonal of the empirical Fisher, diag_p = Σ_e w_e g_e[p]², via
    batched per-event scores (one chunk of [c, n_param] at a time — never the
    full [n_events, n_param]). The Jacobi preconditioner M⁻¹ = 1/(diag + λ)
    for the (ill-conditioned) flow-Fisher CG."""
    diag = torch.zeros(n_param, dtype=torch.float64, device=dev)
    t0 = time.time()
    for ci, i0 in enumerate(range(0, n_events, chunk)):
        per, w = per_event_fn(i0, min(i0 + chunk, n_events))
        c = int(per.shape[0])
        eye = torch.eye(c, device=dev, dtype=per.dtype)
        g = torch.autograd.grad(per, model_params, grad_outputs=eye,
                                is_grads_batched=True, retain_graph=False,
                                allow_unused=True)
        off = 0
        for gi, p in zip(g, model_params):
            if gi is not None:
                S = gi.reshape(c, -1).double()                      # [c, numel]
                diag[off:off + p.numel()] += (w.double().unsqueeze(1) * S * S).sum(0)
            off += p.numel()
        if progress and ci % 50 == 0:
            print(f"  preconditioner diag: chunk {ci + 1} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    return diag


def _pcg_batched(apply_A_batch, B, Minv, tol, max_iter, label="",
                 progress=True, report_dt=20.0):
    """Jacobi-PRECONDITIONED, BATCHED CG: solve A·Z = B for the K columns of
    B [n, K] simultaneously (per-column scalars; a column freezes once its
    relative residual < tol). ``apply_A_batch``: [n, K] → [n, K] (one shared
    matvec per iteration). ``Minv`` [n] = the diagonal preconditioner. A is
    SPD here (empirical Fisher + ridge), so no negative-curvature guard is
    needed. Returns (Z [n, K], iters, max relative residual)."""
    n, K = B.shape
    X = torch.zeros_like(B)
    R = B.clone()
    Z = Minv.unsqueeze(1) * R
    P = Z.clone()
    rz = (R * Z).sum(0)                                  # [K]
    bnorm = R.norm(dim=0).clamp_min(1e-300)              # [K]
    rel = R.norm(dim=0) / bnorm
    it = 0
    t_last = time.time()
    while it < max_iter and bool((rel > tol).any()):
        active = (rel > tol).to(B.dtype)                 # [K] 1/0
        AP = apply_A_batch(P)
        pAp = (P * AP).sum(0).clamp_min(1e-300)          # [K] (SPD → >0)
        alpha = active * rz / pAp
        X = X + alpha.unsqueeze(0) * P
        R = R - alpha.unsqueeze(0) * AP
        Znew = Minv.unsqueeze(1) * R
        rz_new = (R * Znew).sum(0)
        beta = active * rz_new / rz.clamp_min(1e-300)
        P = Znew + beta.unsqueeze(0) * P
        rz = rz_new
        rel = R.norm(dim=0) / bnorm
        it += 1
        if progress and (time.time() - t_last) > report_dt:
            print(f"        PCG[{label}] it {it}/{max_iter}  "
                  f"max rel {float(rel.max()):.2e} "
                  f"({int((rel > tol).sum())}/{K} active, tol {tol:g})",
                  flush=True)
            t_last = time.time()
    return X, it, float(rel.max())


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    # ArgumentDefaultsHelpFormatter appends "(default: …)" to every option's
    # help line, so -h prints the default for each argument.
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--checkpoint", required=True,
                   help="Stage-2 fit checkpoint (fit_best.pt). REQUIRED.")
    p.add_argument("--shards", required=True,
                   help="Arrow shard dir/file(s). REQUIRED.")
    p.add_argument("--device",
                   default="cuda:0" if torch.cuda.is_available() else "cpu",
                   help="Compute device (default tracks CUDA availability).")
    p.add_argument("--batch-size", type=int, default=262144,
                   help="Loader batch size for streaming the event samples.")
    p.add_argument("--eval-chunk", type=int, default=4096,
                   help="Events per HVP autograd chunk — the PEAK-MEMORY knob. "
                   "Each chunk holds a full double-backward graph (the gh_qop "
                   "quadrature × the flow forward, retained for the 2nd "
                   "derivative), so memory scales with this. LOWER it (e.g. "
                   "2048/1024) if you hit a CUDA OOM or the NVML_SUCCESS "
                   "allocator assert; the result is identical (pure "
                   "accumulation), only peak memory changes.")
    p.add_argument("--max-events-fit", type=int, default=0,
                   help="Cap the materialised stage-2 (fit) events; 0 = all. "
                   "VALID speedup: the stage-2 normalisation CANCELS exactly in "
                   "the chain (net power 0), so a subset gives an UNBIASED "
                   "full-fit covariance with NO rescale — only extra estimation "
                   "noise on H_ww/H_wφ. Keep N ≫ n_w and enough to condition "
                   "the profiled background block (else raise --ridge-w).")
    p.add_argument("--max-events-flow", type=int, default=0,
                   help="Cap the materialised stage-1 (flow) events; 0 = all. "
                   "VALID speedup (the biggest CG_φ one): the covariance is "
                   "rescaled by the full/seen Σw ratio (α₁) to the FULL flow "
                   "sample — unbiased, with extra noise (the ridge regularises "
                   "the under-determined flow Hessian; the chain only probes its "
                   "≤n_θ-dim subspace, so N need not exceed n_φ).")
    p.add_argument("--ridge-w", type=float, default=1e-6,
                   help="RELATIVE damping for H_ww: λ_w = ridge·tr(H_ww)/n_w "
                   "(trace from Hutchinson probes at startup), so the value "
                   "is sample-size and unit independent.")
    p.add_argument("--ridge-flow", type=float, default=1e-4,
                   help="RELATIVE damping for H₁: λ_φ = ridge·tr(H₁)/n_φ — "
                   "the regulariser of the near-null (early-stopping-flat) "
                   "flow directions. SCAN this (e.g. ×10 up/down) and quote "
                   "the plateau.")
    p.add_argument("--ridge-escalate", type=float, default=10.0,
                   help="On a non-PD CG abort, multiply that solve's ridge by "
                   "this factor and retry (the escalated ridge persists across "
                   "columns). Keeps an indefinite/under-sampled Hessian from "
                   "crashing the run.")
    p.add_argument("--ridge-escalate-max", type=int, default=4,
                   help="Max ridge escalations per solve before giving up "
                   "(0 = disable → hard-fail on non-PD, the strict behaviour). "
                   "NOTE: escalation OVER-regularises and UNDERESTIMATES the "
                   "flow σ in affected columns — it is a keep-alive/rough-"
                   "bound, not a faithful covariance; a trustworthy run needs "
                   "ZERO escalations (add events / scan --ridge-flow).")
    p.add_argument("--w-solver", choices=["fisher", "observed", "minres"],
                   default="fisher",
                   help="How to invert the stage-2 bread H_ww. 'fisher' "
                   "(default): the empirical Fisher J_ww = Σ wᵢsᵢsᵢᵀ from "
                   "per-event scores — PSD by construction (the nonlinear "
                   "background block's OBSERVED Hessian is genuinely "
                   "indefinite), consistent with the --fisher data-stat "
                   "covariance (same bread), and built once (dense Cholesky "
                   "solve for all θ columns, no CG_w). 'observed': the "
                   "double-backward Hessian via CG (+ ridge escalation) — the "
                   "exact estimating-equation Jacobian, but indefinite for the "
                   "background MLP. 'minres': observed H_ww via MINRES (handles "
                   "the indefinite-but-invertible matrix; the θθ covariance is "
                   "PSD by congruence regardless) with a small ridge for the "
                   "near-null background directions — a cross-check of 'fisher' "
                   "(they agree iff the negative directions decouple from θ).")
    p.add_argument("--score-chunk", type=int, default=64,
                   help="(--w-solver fisher) Per-event-score batch size for the "
                   "empirical-Fisher build (batched vjp). Smaller = less memory.")
    p.add_argument("--flow-solver", choices=["fisher", "observed"],
                   default="fisher",
                   help="How to form the flow meat H₁ (the CENTRAL PSD factor "
                   "Cov(φ̂), which MUST be PD). 'fisher' (default): the "
                   "empirical Fisher I_φ = Σ wᵢsᵢsᵢᵀ via a MATRIX-FREE "
                   "Fisher-vector product (58208-dim, can't densify) — PSD by "
                   "construction (the standard inverse-Fisher parameter "
                   "covariance), so CG_φ never hits negative curvature and only "
                   "a small null-space --ridge-flow is needed. 'observed': the "
                   "double-backward flow Hessian, which is genuinely indefinite "
                   "at the early-stopping point (not a minimum) → needs a large "
                   "--ridge-flow / escalation.")
    p.add_argument("--check-stationarity", action="store_true",
                   help="Before solving, report ‖∂L₂/∂w‖ (total, θ-block, "
                   "background-block) at the checkpoint. A non-PSD/indefinite "
                   "H_ww is only a meaningful covariance if the fit is "
                   "STATIONARY (gradient ≈ 0); a large θ-block gradient means "
                   "the fit did not converge and NO solver gives a valid σ.")
    p.add_argument("--grid-nphi", type=int, default=0,
                   help="θ-mlp only: uniform φ bins of the fixed (η,φ) output "
                   "grid (η = the stats η-bin centres); 0 = match the "
                   "checkpoint's output_fisher_nphi (itself default 4) so the "
                   "covariance combines with the --output-fisher file. The "
                   "columns are EXACT — cost is linear in n_η·n_φ·n_comp.")
    p.add_argument("--flow-precond", choices=["jacobi", "none"],
                   default="jacobi",
                   help="(--flow-solver fisher) Preconditioner for the "
                   "ill-conditioned flow-Fisher CG_φ. 'jacobi': M⁻¹ = "
                   "1/(diag(I_φ) + λ_φ), one extra pass to build the exact "
                   "diagonal — helps when the Fisher's scale varies widely "
                   "across the flow weights (the usual case), but can HURT a "
                   "well-scaled/correlated Fisher, so the per-block iteration "
                   "counts are printed: compare against 'none' on your problem.")
    p.add_argument("--col-block", type=int, default=8,
                   help="Solve the θ-columns in BLOCKS of this size — the "
                   "mixed step and (for --flow-solver fisher) the "
                   "preconditioned CG_φ apply their matvec to the whole block "
                   "at once, sharing the per-chunk forward across columns (big "
                   "speedup; collapses the per-column Python loop). Memory "
                   "scales with the block (the batched backward holds a "
                   "block-fold graph) — lower it on OOM, 1 = per-column.")
    p.add_argument("--cg-tol", type=float, default=1e-4,
                   help="CG relative-residual tolerance. VALID speedup: loosen "
                   "to ~1e-3 to roughly halve the iterations (covariance error "
                   "~O(tol) — fine for an uncertainty); the printed asymmetry "
                   "is the convergence check. (Unlike raising --ridge-*, which "
                   "biases the covariance, loosening tol does not.)")
    p.add_argument("--cg-max-iter-w", type=int, default=200,
                   help="Max CG iterations for the H_ww (stage-2) solve.")
    p.add_argument("--cg-max-iter-flow", type=int, default=200,
                   help="Max CG iterations for the H₁ (stage-1 flow) solve.")
    p.add_argument("--precision", choices=["fp32", "fp64"], default="fp64",
                   help="Compute precision (fp64 strongly recommended: CG "
                   "residuals live below fp32 HVP noise).")
    p.add_argument("--fisher", default=None,
                   help="Optional empirical_fisher.pt to combine with "
                   "(data-stat + flow totals written and printed); None = skip.")
    p.add_argument("--output", default=None,
                   help="Output .pt; None → <ckpt dir>/flow_uncertainty.pt.")
    p.add_argument("--progress", action="store_true", default=True,
                   help="Print per-column CG progress (on by default).")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    dev = args.device
    if dev.startswith("cuda") and not torch.cuda.is_available():
        print("CUDA unavailable; using CPU")
        dev = "cpu"
    dtype = torch.float64 if args.precision == "fp64" else torch.float32

    print(f"loading checkpoint: {args.checkpoint}")
    model, stats, targs, ckpt = load_model_from_checkpoint(args.checkpoint, dev)
    if str(ckpt.get("stage", "")) == "flow":
        print("error: this is a stage-1 flow checkpoint; run on fit_best.pt",
              file=sys.stderr)
        return 1
    is_mlp = getattr(model, "theta_mode", "binned") == "mlp"
    arch = getattr(model, "flow_arch", "gf")
    if arch == "nce":
        print("error: the nce stage-1 objective (paired BCE) is not an NLL; "
              "flow-uncertainty propagation is not defined here.",
              file=sys.stderr)
        return 1
    if dtype == torch.float64:
        model.double()
    model.eval()

    flow_params = [p for p in model.flow.parameters()]
    for p in flow_params:
        p.requires_grad_(True)
    n_phi = sum(p.numel() for p in flow_params)
    if is_mlp:
        params, segs, n_net = _w_layout_mlp(model)
        n_w_qoi = n_net          # θ-net weights = QoI block (vs background)
        scale_cols = ([c for c in range(3)
                       if float(model.scale_param_mask[c]) != 0.0]
                      if model.scale_enabled else [])
        smear_cols = (T._smear_active_cols(model)
                      if model.smearing_enabled else [])
    else:
        params, segs, labels, n_theta, smear_cols = _w_layout(model)
        n_w_qoi = n_theta        # active θ entries = QoI block (vs background)
    for p in params:
        p.requires_grad_(True)
    n_w = segs[-1][0] + segs[-1][1].numel()
    if is_mlp:
        # Fixed (η, φ) output grid: T_j = the θ-net outputs at the grid points
        # in the output-fisher table layout; the chain columns are the delta-
        # method RHS g_j = ∂T_j/∂w (zeros on the background-MLP block).
        n_phi_g = args.grid_nphi or int(targs.get("output_fisher_nphi", 4)
                                        or 4)
        o_active, labels, eta_centres, phi_centres = _grid_outputs(
            model, stats.eta_edges, n_phi_g, scale_cols, smear_cols, dev)
        n_theta = int(o_active.numel())
        n_eta_g = len(eta_centres)
        n_net_tensors = sum(1 for _ in model.theta_net.parameters())
        net_params = params[:n_net_tensors]
        net_segs = segs[:n_net_tensors]
        G = torch.zeros((n_theta, n_w), dtype=torch.float64, device=dev)
        for k in range(n_theta):
            g = torch.autograd.grad(o_active[k], net_params,
                                    retain_graph=(k < n_theta - 1),
                                    allow_unused=True)
            G[k] = _pack(g, net_params, net_segs, n_w, dev)
        del o_active
        print(f"w-space: {n_w} active ({n_net} θ-net + {n_w - n_net} "
              f"background); φ-space: {n_phi} flow weights; arch={arch}")
        print(f"output grid: {n_eta_g} η × {n_phi_g} φ × "
              f"{len(scale_cols) + len(smear_cols)} comp = {n_theta} columns "
              f"(physical A,e,M; O(1) effective a,c)")
    else:
        print(f"w-space: {n_w} active ({n_theta} θ + {n_w - n_theta} "
              f"background); φ-space: {n_phi} flow weights; arch={arch}")

    # ---- loaders (mirror the fit's and the flow's event selections) -------
    shard_files = discover_shards([args.shards])
    validation = bool(targs.get("validation", False))
    no_split = bool(targs.get("no_validation_split", False))
    half_fit = (1 if validation and not no_split else None)
    half_flow = (0 if validation and not no_split else None)
    n_eta = len(stats.eta_edges) - 1
    inj = inj_sm = None
    inj_bkg = None
    if validation:
        ia, ie, im = (float(targs.get(k, 0.0) or 0.0)
                      for k in ("inject_A", "inject_e", "inject_M"))
        if ia or ie or im:
            inj = np.zeros((n_eta, 3))
            inj[:, 0], inj[:, 1], inj[:, 2] = ia, ie, im
        isa, isc = (float(targs.get(k, 0.0) or 0.0)
                    for k in ("inject_a", "inject_c"))
        if isa or isc:
            inj_sm = np.zeros((n_eta, 2))
            inj_sm[:, 0], inj_sm[:, 1] = isa, isc
        bf0 = float(targs.get("inject_bkg_f0", 0.0) or 0.0)
        bf1 = float(targs.get("inject_bkg_f1", 0.0) or 0.0)
        if bf0 > 0.0 or bf1 > 0.0:
            inj_bkg = (bf0, bf1)
    inj_prod = None
    if validation:
        sp = (float(targs.get("inject_prod_ptll_slope", 0.0) or 0.0),
              float(targs.get("inject_prod_yll_slope", 0.0) or 0.0),
              float(targs.get("inject_prod_costheta_slope", 0.0) or 0.0))
        if any(v != 0.0 for v in sp):
            inj_prod = sp
    fl_lo = float(targs.get("flow_m_lo") or stats.m_lo)
    fl_hi = float(targs.get("flow_m_hi") or stats.m_hi)
    ft_lo = float(targs.get("fit_m_lo") or fl_lo)
    ft_hi = float(targs.get("fit_m_hi") or fl_hi)
    common = dict(batch_size=args.batch_size, split="train", val_fraction=0.0,
                  holdout_fraction=0.0, drop_last=False,
                  inject_seed=int(targs.get("inject_smear_seed", 12345)),
                  cond_basis=targs.get("cond_basis", "muon_kin"),
                  inject_nonuniform=bool(targs.get("inject_nonuniform", False)),
                  event_fraction=float(targs.get("event_fraction", 1.0) or 1.0))
    fit_loader = JpsiMassArrowLoader(
        shard_files, stats, half=half_fit, inject_theta_scale=inj,
        inject_theta_smear=inj_sm, inject_bkg=inj_bkg,
        m_window=(ft_lo, ft_hi), inject_prod=inj_prod, **common)
    fm = targs.get("flow_monitor", "train")
    flow_loader = JpsiMassArrowLoader(
        shard_files, stats, half=half_flow,
        m_window=(fl_lo, fl_hi), **{**common,
        "val_fraction": (0.0 if fm == "train"
                         else float(targs.get("val_fraction", 0.1))),
        "holdout_fraction": (0.0 if fm == "train"
                             else float(targs.get("holdout_fraction", 0.05)))})

    print("materialising event samples...")
    fit_ev, _, _ = _materialise(
        fit_loader, dev, dtype, take_data_branch=True, mc_as_data=validation,
        cap=args.max_events_fit, label="stage-2 (fit)")
    flow_ev, w1_seen, w1_total = _materialise(
        flow_loader, dev, dtype, take_data_branch=False, mc_as_data=True,
        cap=args.max_events_flow, label="stage-1 (flow)")
    alpha1 = w1_total / max(w1_seen, 1e-300)
    if dev.startswith("cuda"):
        torch.cuda.empty_cache()   # release the loader's transient buffers

    n_iter = int(targs.get("continuity_n_iter", 2))

    def fit_loss(i0, i1):
        sl = slice(i0, i1)
        per = model.data_nll_continuity(
            fit_ev["mll"][sl], fit_ev["pt_pm"][sl], fit_ev["eta_pm"][sl],
            fit_ev["phi_pm"][sl], fit_ev["q_pm"][sl], fit_ev["b_pm"][sl],
            fit_ev["cond_std"][sl],
            torch.ones(i1 - i0, dtype=torch.bool, device=dev), n_iter=n_iter)
        return (fit_ev["w"][sl] * per).sum()

    window_norm = (arch in ("gf", "nsf")
                   and not bool(targs.get("no_flow_window_norm", False)))
    gauge_lambda = float(targs.get("flow_gauge_penalty", 1e-3) or 0.0)

    def flow_loss(i0, i1):
        # Replicates the step1 objective: plain NLL for the window-defined
        # flows; truncated window-norm + gauge penalty for gf/nsf.
        sl = slice(i0, i1)
        m = flow_ev["mll"][sl]
        mk = flow_ev["cond_std"][sl]
        w = flow_ev["w"][sl]
        logp = model.log_p_nominal(m, mk)
        loss = -(w * logp).sum()
        if window_norm:
            log_Z = model._flow_log_window_Z(
                m.new_full(m.shape, model._flow_m_lo_f),
                m.new_full(m.shape, model._flow_m_hi_f), mk)
            loss = loss + (w * log_Z).sum()
            if gauge_lambda > 0.0:
                loss = loss + gauge_lambda * (w * log_Z ** 2).sum()
        return loss

    L2 = ChunkedLoss(fit_loss, fit_ev["mll"].shape[0], args.eval_chunk)
    L1 = ChunkedLoss(flow_loss, flow_ev["mll"].shape[0], args.eval_chunk)

    def _raw_hvp_w(v):
        vl = _unpack(v, params, segs)
        return _pack(L2.hvp(vl, params), params, segs, n_w, dev)

    def _unpack_phi(v):
        vl, off = [], 0
        for p in flow_params:
            vl.append(v[off:off + p.numel()].view_as(p).to(p.dtype))
            off += p.numel()
        return vl

    def _raw_hvp_phi(v):
        h = L1.hvp(_unpack_phi(v), flow_params)
        return torch.cat([x.reshape(-1) for x in h]).double()

    n_flow = flow_ev["mll"].shape[0]

    def flow_per_event(i0, i1):
        # Per-event flow NLL (= the empirical-Fisher score source): −log p₀
        # (+ the window-norm term for gf/nsf). The gauge penalty is a
        # regulariser, NOT a per-event likelihood term, so it is excluded —
        # the empirical Fisher I_φ = Σ wᵢsᵢsᵢᵀ is the LIKELIHOOD information.
        sl = slice(i0, i1)
        m = flow_ev["mll"][sl]
        mk = flow_ev["cond_std"][sl]
        per = -model.log_p_nominal(m, mk)
        if window_norm:
            per = per + model._flow_log_window_Z(
                m.new_full(m.shape, model._flow_m_lo_f),
                m.new_full(m.shape, model._flow_m_hi_f), mk)
        return per, flow_ev["w"][sl]

    def _raw_efvp_phi(v):
        h = _efvp(flow_per_event, flow_params, n_flow, args.eval_chunk,
                  _unpack_phi(v))
        return torch.cat([x.reshape(-1) for x in h]).double()

    # The flow meat operator: empirical Fisher (PSD, default) or observed H₁.
    raw_phi = _raw_efvp_phi if args.flow_solver == "fisher" else _raw_hvp_phi

    offs_phi, _o = [], 0
    for p in flow_params:
        offs_phi.append((_o, p.numel()))
        _o += p.numel()

    def mixed_u(x):
        xl = _unpack(x, params, segs)
        h = L2.mixed(xl, params, flow_params)
        return torch.cat([t.reshape(-1) for t in h]).double()

    def mixed_u_batch(X):                       # [n_w, K] → [n_φ, K]
        return L2.mixed_batch(X, params, segs, n_w, flow_params, n_phi, dev)

    def efvp_phi_batch(V):                       # [n_φ, K] → [n_φ, K]
        return _efvp_batch(flow_per_event, flow_params, n_flow,
                           args.eval_chunk, V, offs_phi)

    n_fit = fit_ev["mll"].shape[0]

    def fit_per_event(i0, i1):
        sl = slice(i0, i1)
        per = model.data_nll_continuity(
            fit_ev["mll"][sl], fit_ev["pt_pm"][sl], fit_ev["eta_pm"][sl],
            fit_ev["phi_pm"][sl], fit_ev["q_pm"][sl], fit_ev["b_pm"][sl],
            fit_ev["cond_std"][sl],
            torch.ones(i1 - i0, dtype=torch.bool, device=dev), n_iter=n_iter)
        return per, fit_ev["w"][sl]

    # ---- stationarity check: ‖∂L₂/∂w‖ at the checkpoint -------------------
    # The covariance is only meaningful if the fit is at a stationary point
    # (∂L₂/∂w ≈ 0); a large QoI-block gradient means the fit did not converge.
    if args.check_stationarity:
        print("checking stage-2 stationarity (‖∂L₂/∂w‖)...")
        acc = [torch.zeros_like(p, dtype=torch.float64) for p in params]
        for i0 in range(0, n_fit, args.eval_chunk):
            L = fit_loss(i0, min(i0 + args.eval_chunk, n_fit))
            g = torch.autograd.grad(L, params, allow_unused=True)
            for a, gi in zip(acc, g):
                if gi is not None:
                    a += gi.double()
        gvec = _pack(acc, params, segs, n_w, dev)
        g_qoi = float(gvec[:n_w_qoi].norm())
        g_bkg = float(gvec[n_w_qoi:].norm())
        # per-event Σw to gauge a "small" gradient scale (grad ∝ Σw).
        sw_fit = float(fit_ev["w"].sum())
        print(f"  ‖∂L₂/∂w‖ total = {float(gvec.norm()):.3e}  "
              f"(QoI {g_qoi:.3e}, background {g_bkg:.3e}); Σw = {sw_fit:.3e}; "
              f"QoI rel = {g_qoi / max(sw_fit, 1e-300):.2e}")
        if g_qoi / max(sw_fit, 1e-300) > 1e-3:
            print("  WARNING: large QoI-block gradient — the fit may not be "
                  "converged; NO solver gives a valid covariance at a "
                  "non-stationary point.", file=sys.stderr)
        if dev.startswith("cuda"):
            torch.cuda.empty_cache()

    # ---- ridge scales (Hutchinson for the directions still solved by CG) --
    gen = torch.Generator(device="cpu").manual_seed(7)

    def _trace_scale(raw_hvp, n, label):
        acc = 0.0
        for k in range(2):
            t = time.time()
            v = (torch.randint(0, 2, (n,), generator=gen,
                               dtype=torch.int64).double() * 2.0 - 1.0).to(dev)
            acc += float(v @ raw_hvp(v)) / n
            print(f"  trace[{label}] probe {k + 1}/2 done "
                  f"({time.time() - t:.0f}s)", flush=True)
        return acc / 2.0

    # H₁ (flow) ridge always via Hutchinson; it stays a CG solve (on the
    # empirical Fisher I_φ for --flow-solver fisher, the observed H₁ otherwise).
    _phi_tag = "I_φ" if args.flow_solver == "fisher" else "H₁(obs)"
    print(f"estimating {_phi_tag} trace scale (Hutchinson; --eval-chunk="
          f"{args.eval_chunk})...")
    sc_phi = _trace_scale(raw_phi, n_phi, _phi_tag)
    lam_phi = args.ridge_flow * abs(sc_phi)

    # ---- H_ww bread: build/ridge per --w-solver --------------------------
    L_chol = X_fisher = None
    if args.w_solver == "fisher":
        print(f"building empirical Fisher J_ww (PSD; per-event scores, "
              f"--score-chunk={args.score_chunk})...")
        J_ww = _empirical_fisher_dense(
            fit_per_event, n_fit, args.score_chunk, params, segs, n_w, dev,
            progress=args.progress)
        sc_w = float(torch.diagonal(J_ww).mean())          # tr(J_ww)/n
        lam_w = args.ridge_w * abs(sc_w)
        eye_w = torch.eye(n_w, dtype=torch.float64, device=dev)
        L_chol = torch.linalg.cholesky(J_ww + lam_w * eye_w)
        del J_ww
        print(f"  tr(J_ww)/n = {sc_w:.4e} → λ_w = {lam_w:.4e}; "
              f"tr({_phi_tag})/n = {sc_phi:.4e} → λ_φ = {lam_phi:.4e}")
    else:
        print(f"estimating H_ww trace scale (Hutchinson)...")
        sc_w = _trace_scale(_raw_hvp_w, n_w, "H_ww")
        lam_w = args.ridge_w * abs(sc_w)
        print(f"  tr(H_ww)/n = {sc_w:.4e} → λ_w = {lam_w:.4e}; "
              f"tr({_phi_tag})/n = {sc_phi:.4e} → λ_φ = {lam_phi:.4e}")
    if dev.startswith("cuda"):
        torch.cuda.empty_cache()

    # RHS matrix B_rhs [n_w, n_theta]: e_j (binned) or g_j = ∂T_j/∂w (mlp).
    if is_mlp:
        B_rhs = G.t().contiguous()
    else:
        B_rhs = torch.zeros((n_w, n_theta), dtype=torch.float64, device=dev)
        B_rhs[torch.arange(n_theta), torch.arange(n_theta)] = 1.0
    if args.w_solver == "fisher":
        X_fisher = torch.cholesky_solve(B_rhs, L_chol)     # all columns at once

    # ---- robust damped solve: on a non-PD CG abort, escalate the ridge
    # (×escalate, up to escalate_max times) and retry, so an indefinite
    # (e.g. under-sampled) Hessian yields a PSD result instead of crashing.
    # The escalated ridge PERSISTS across columns (ridge_state) so only the
    # first affected column pays the search. WARNING: a larger ridge SHRINKS
    # the covariance, so escalated columns UNDERESTIMATE the flow uncertainty
    # — the result there is over-regularised, not faithful (add events for a
    # 'w' solve / scan --ridge-flow for a 'φ' solve). All escalations are
    # logged and recorded in the output.
    ridge_state = {"w": lam_w, "phi": lam_phi}
    escal = {"w": {"n_cols": 0, "max_ridge": lam_w},
             "phi": {"n_cols": 0, "max_ridge": lam_phi}}

    def _solve(raw_hvp, b, key, label, max_iter):
        lam = ridge_state[key]
        n_esc = 0
        while True:
            try:
                x, it, res = _cg(lambda p, _l=lam: raw_hvp(p) + _l * p, b,
                                 args.cg_tol, max_iter, label=label,
                                 progress=args.progress)
                ridge_state[key] = lam            # keep working ridge as floor
                if n_esc:
                    escal[key]["n_cols"] += 1
                    escal[key]["max_ridge"] = max(escal[key]["max_ridge"], lam)
                return x, it, res, lam
            except RuntimeError as ex:
                if ("non-positive curvature" not in str(ex)
                        or args.ridge_escalate_max <= 0
                        or n_esc >= args.ridge_escalate_max):
                    if args.ridge_escalate_max > 0:
                        raise RuntimeError(
                            f"{label}: ridge escalation exhausted "
                            f"({args.ridge_escalate_max} × "
                            f"{args.ridge_escalate:g}); the (damped) Hessian "
                            f"is still not PSD at λ={lam:.3e}. The matrix is "
                            f"badly indefinite — add events ('w': "
                            f"--max-events-fit) or pre-scan the ridge; no "
                            f"regularisation recovers a faithful covariance "
                            f"here.") from ex
                    raise
                lam *= args.ridge_escalate
                n_esc += 1
                print(f"    WARNING [{label}]: non-PD → escalating ridge to "
                      f"{lam:.3e} (×{args.ridge_escalate:g}, "
                      f"attempt {n_esc}/{args.ridge_escalate_max}) — this "
                      f"column will be over-regularised (flow σ UNDERestimated)",
                      flush=True)

    def _solve_w(j):
        """x_j = (H_ww + λ_w)⁻¹ B_rhs[:, j] via the chosen --w-solver."""
        if args.w_solver == "fisher":
            return X_fisher[:, j], 0, 0.0, lam_w        # direct (Cholesky) solve
        b = B_rhs[:, j]
        if args.w_solver == "minres":
            x, it, res = _minres(lambda p: _raw_hvp_w(p) + lam_w * p, b,
                                 args.cg_tol, args.cg_max_iter_w,
                                 label=f"w:{labels[j]}", progress=args.progress)
            return x, it, res, lam_w
        return _solve(_raw_hvp_w, b, "w", f"w:{labels[j]}", args.cg_max_iter_w)

    # ---- Jacobi preconditioner for the (ill-conditioned) flow-Fisher CG ----
    Minv_phi = None
    if args.flow_solver == "fisher" and args.flow_precond == "jacobi":
        print(f"building Jacobi preconditioner diag(I_φ) "
              f"(--score-chunk={args.score_chunk})...")
        diag_phi = _fisher_diag(flow_per_event, flow_params, n_flow,
                                args.score_chunk, n_phi, dev,
                                progress=args.progress)
        Minv_phi = 1.0 / (diag_phi + lam_phi)
        if dev.startswith("cuda"):
            torch.cuda.empty_cache()
    elif args.flow_solver == "fisher":
        Minv_phi = torch.ones(n_phi, dtype=torch.float64, device=dev)

    # ---- the column loop (BLOCKED over --col-block columns) ----------------
    bsz = max(1, int(args.col_block))
    batched = args.flow_solver == "fisher"     # batched mixed + PCG_φ path
    print(f"running {n_theta} θ-columns in blocks of {bsz} "
          f"(w-solver={args.w_solver}, flow-solver={args.flow_solver}, "
          f"λ_w={lam_w:.3e}, λ_φ={lam_phi:.3e}, tol={args.cg_tol:g}, "
          f"α₁={alpha1:.4f})...")
    U = torch.zeros((n_theta, n_phi), dtype=torch.float64)
    Z = torch.zeros((n_theta, n_phi), dtype=torch.float64)
    t0 = time.time()
    cg_stats = []
    for bs in range(0, n_theta, bsz):
        cols = list(range(bs, min(bs + bsz, n_theta)))
        K = len(cols)
        # w-solve for the block (fisher: slice the precomputed solve).
        if args.w_solver == "fisher":
            Xb = X_fisher[:, cols]                      # [n_w, K]
            it_w, res_w, lw_used = 0, 0.0, lam_w
        else:
            xs, it_w, res_w, lw_used = [], 0, 0.0, lam_w
            for j in cols:
                print(f"  [{j + 1:3d}/{n_theta}] {labels[j]:14s} solving "
                      f"{args.w_solver}_w...", flush=True)
                xj, itj, resj, lw_used = _solve_w(j)
                xs.append(xj); it_w = max(it_w, itj); res_w = max(res_w, resj)
            Xb = torch.stack(xs, dim=1)
        # mixed: U = H_φw X  (batched when fisher-flow, else per column).
        if batched:
            Ub = mixed_u_batch(Xb)                      # [n_φ, K]
            Zb, it_f, res_f = _pcg_batched(
                lambda P: efvp_phi_batch(P) + lam_phi * P, Ub, Minv_phi,
                args.cg_tol, args.cg_max_iter_flow,
                label=f"φ:{labels[cols[0]]}..", progress=args.progress)
            lf_used = lam_phi
        else:
            Ucols, Zcols, it_f, res_f, lf_used = [], [], 0, 0.0, lam_phi
            for k, j in enumerate(cols):
                u = mixed_u(Xb[:, k])
                z, itf, resf, lf_used = _solve(
                    raw_phi, u, "phi", f"φ:{labels[j]}", args.cg_max_iter_flow)
                Ucols.append(u); Zcols.append(z)
                it_f = max(it_f, itf); res_f = max(res_f, resf)
            Ub = torch.stack(Ucols, dim=1)
            Zb = torch.stack(Zcols, dim=1)
        for k, j in enumerate(cols):
            U[j] = Ub[:, k].cpu()
            Z[j] = Zb[:, k].cpu()
            cg_stats.append((it_w, res_w, it_f, res_f, lw_used, lf_used))
        if dev.startswith("cuda"):
            torch.cuda.empty_cache()   # limit fragmentation over the long loop
        el = time.time() - t0
        w_tag = ("Fisher" if args.w_solver == "fisher"
                 else f"{args.w_solver}_w {it_w:3d} it")
        phi_tag = (f"PCG_φ {it_f:3d} it (max res {res_f:.1e})" if batched
                   else f"CG_φ {it_f:3d} it (res {res_f:.1e})")
        print(f"  [{cols[0] + 1:3d}-{cols[-1] + 1:3d}/{n_theta}] {w_tag}  "
              f"{phi_tag}  elapsed {el / 60:.1f} min", flush=True)

    cov_flow = (U @ Z.T) / alpha1
    cov_flow = 0.5 * (cov_flow + cov_flow.T)          # numerical symmetrise
    asym = float((U @ Z.T - Z @ U.T).abs().max()
                 / cov_flow.abs().max().clamp_min(1e-300))
    print(f"covariance asymmetry (CG-residual scale check): {asym:.2e}")
    if escal["w"]["n_cols"] or escal["phi"]["n_cols"]:
        print("\n*** WARNING: ridge auto-escalation was triggered — the result "
              "is OVER-REGULARISED and UNDERESTIMATES the flow uncertainty in "
              "the affected columns (NOT a faithful covariance):")
        if escal["w"]["n_cols"]:
            print(f"    H_ww: {escal['w']['n_cols']}/{n_theta} columns, "
                  f"λ_w up to {escal['w']['max_ridge']:.3e} "
                  f"(from {lam_w:.3e}) — H_ww is indefinite; "
                  f"INCREASE --max-events-fit (≫ n_w={n_w}).")
        if escal["phi"]["n_cols"]:
            print(f"    H₁:   {escal['phi']['n_cols']}/{n_theta} columns, "
                  f"λ_φ up to {escal['phi']['max_ridge']:.3e} "
                  f"(from {lam_phi:.3e}) — SCAN --ridge-flow to a plateau "
                  f"instead of relying on escalation.")
        print("    Use this only as a rough/bounding figure; for the quoted "
              "number, rerun with no escalation needed.\n")

    out = {
        "covariance_flow": cov_flow,
        "labels": labels if is_mlp else labels[:n_theta],
        "smear_cols": smear_cols,
        "alpha1": alpha1,
        "ridge_w": args.ridge_w,
        "ridge_flow": args.ridge_flow,
        "lambda_w_abs": lam_w,
        "lambda_flow_abs": lam_phi,
        "cg_tol": args.cg_tol,
        "cg_stats": cg_stats,
        "asymmetry": asym,
        "ridge_escalation": escal,   # {} entries non-zero ⇒ over-regularised
        "w_solver": args.w_solver,
        "flow_solver": args.flow_solver,
        "flow_precond": args.flow_precond,
        "col_block": int(args.col_block),
    }

    if is_mlp:
        # ``covariance_flow`` is the active-table grid covariance (physical
        # A,e,M; O(1) effective a,c) in the output-fisher ordering. Collapse
        # over φ (occupancy-weighted) and emit the same per-η keys as the
        # --output-fisher file, ``_flow`` suffix; if --fisher points at a
        # matching output-fisher file, also write the blockwise totals.
        aw, cell_w = _grid_phi_weights(fit_ev, n_eta_g, n_phi_g)
        ef = None
        if args.fisher and os.path.exists(args.fisher):
            ef = torch.load(args.fisher, map_location="cpu",
                            weights_only=False)
            if (ef.get("theta_mode") != "mlp"
                    or int(ef.get("n_phi", -1)) != n_phi_g):
                print(f"warning: --fisher file is not an --output-fisher "
                      f"θ-mlp file with n_phi={n_phi_g} "
                      f"(theta_mode={ef.get('theta_mode')!r}, "
                      f"n_phi={ef.get('n_phi')!r}); not combined",
                      file=sys.stderr)
                ef = None
            else:
                pw = ef.get("phi_weights")
                if pw is not None and tuple(pw.shape) == (n_eta_g, n_phi_g):
                    # Collapse with the SAME φ-weights the data file used so
                    # the collapsed totals are exactly additive.
                    aw = pw.double().numpy()
        out.update({
            "theta_mode": "mlp",
            "n_phi": n_phi_g,
            "scale_cols": scale_cols,
            "eta_centres": torch.tensor(eta_centres, dtype=torch.float32),
            "phi_centres": torch.tensor(phi_centres, dtype=torch.float32),
            "phi_weights": torch.tensor(aw, dtype=torch.float32),
            "cell_w": torch.tensor(cell_w, dtype=torch.float32),
            "param_space": ("2-D (η,φ) θ-table outputs "
                            "(physical A,e,M; O(1) a,c)"),
        })
        out.update(_mlp_grid_keys(cov_flow.numpy(), n_eta_g, n_phi_g,
                                  scale_cols, smear_cols, aw,
                                  suffix="_flow"))
        if ef is not None:
            for k in ("covariance_24_3_24_3", "covariance_scale_2d",
                      "covariance_smear_24_2_24_2", "covariance_smear_2d"):
                kf = f"{k}_flow"
                if (k in ef and kf in out
                        and tuple(ef[k].shape) == tuple(out[kf].shape)):
                    out[f"{k}_total"] = (ef[k].double()
                                         + out[kf].double()).float()
            comp_s, comp_c = ("A", "e", "M"), ("a", "c")
            if "covariance_24_3_24_3_total" in out:
                ct = out["covariance_24_3_24_3_total"].double().numpy()
                out["sigma_scale_24_3_total"] = torch.sqrt(torch.clamp(
                    torch.tensor(np.einsum("icic->ic", ct)), min=0.0)).float()
            if "covariance_smear_24_2_24_2_total" in out:
                ct = out["covariance_smear_24_2_24_2_total"].double().numpy()
                out["sigma_smear_eff_24_2_total"] = torch.sqrt(torch.clamp(
                    torch.tensor(np.einsum("icic->ic", ct)), min=0.0)).float()
            print("\nσ budget (per-η φ-mean, median over η): "
                  "data-stat | flow | total (inflation)")

            def _budget_row(name, sd, sf, st):
                d, f, t = (float(np.median(x)) for x in (sd, sf, st))
                print(f"  {name:3s} {d:.3e} | {f:.3e} | {t:.3e}  "
                      f"(x{t / max(d, 1e-300):.2f})")

            if "sigma_scale_24_3_total" in out and "sigma_scale_24_3" in ef:
                for c in scale_cols:
                    _budget_row(comp_s[c],
                                ef["sigma_scale_24_3"].numpy()[:, c],
                                out["sigma_scale_24_3_flow"].numpy()[:, c],
                                out["sigma_scale_24_3_total"].numpy()[:, c])
            if ("sigma_smear_eff_24_2_total" in out
                    and "sigma_smear_eff_24_2" in ef):
                for c in smear_cols:
                    _budget_row(comp_c[c],
                                ef["sigma_smear_eff_24_2"].numpy()[:, c],
                                out["sigma_smear_eff_24_2_flow"].numpy()[:, c],
                                out["sigma_smear_eff_24_2_total"].numpy()[:, c])
    else:
        n_scale = model.theta_scale.numel() if model.scale_enabled else 0
        out["n_scale"] = n_scale
        extras = T._theta_cov_extras(cov_flow, model, smear_cols, n_scale)
        for k, v in extras.items():
            out[f"{k}_flow"] = v

        if args.fisher and os.path.exists(args.fisher):
            ef = torch.load(args.fisher, map_location="cpu",
                            weights_only=False)
            cov_d = ef["covariance"].double()
            if cov_d.shape == cov_flow.shape:
                cov_tot = cov_d + cov_flow
                out["covariance_total"] = cov_tot
                ex_t = T._theta_cov_extras(cov_tot, model, smear_cols, n_scale)
                for k, v in ex_t.items():
                    out[f"{k}_total"] = v
                sd = cov_d.diagonal().clamp_min(0).sqrt()
                sf = cov_flow.diagonal().clamp_min(0).sqrt()
                st = cov_tot.diagonal().clamp_min(0).sqrt()
                print("\nσ budget (raw θ units): data-stat | flow | total "
                      "(inflation)")
                for j in range(n_theta):
                    infl = float(st[j] / sd[j].clamp_min(1e-300))
                    print(f"  {labels[j]:14s} {float(sd[j]):.3e} | "
                          f"{float(sf[j]):.3e} | {float(st[j]):.3e}  "
                          f"(x{infl:.2f})")
            else:
                print(f"warning: --fisher covariance shape "
                      f"{tuple(cov_d.shape)} does not match the θ-active "
                      f"block {tuple(cov_flow.shape)}; not combined",
                      file=sys.stderr)

    out_path = args.output or os.path.join(
        os.path.dirname(os.path.abspath(args.checkpoint)),
        "flow_uncertainty.pt")
    torch.save(out, out_path)
    print(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
