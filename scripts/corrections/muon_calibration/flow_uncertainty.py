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


def _cg(apply_A, b, tol, max_iter, label="", progress=True):
    """Standard CG on the (damped) SPD system; returns (x, iters, rel_res).
    Aborts with a clear message on negative curvature (raise the ridge)."""
    x = torch.zeros_like(b)
    r = b.clone()
    p = r.clone()
    rs = float(r @ r)
    b_norm = max(float(b.norm()), 1e-300)
    it = 0
    while it < max_iter and (rs ** 0.5) / b_norm > tol:
        Ap = apply_A(p)
        pAp = float(p @ Ap)
        if pAp <= 0.0:
            raise RuntimeError(
                f"CG[{label}]: non-positive curvature (pᵀAp = {pAp:.3e}) at "
                f"iteration {it} — the (damped) Hessian is not PSD here; "
                f"raise the corresponding --ridge-*.")
        alpha = rs / pAp
        x += alpha * p
        r -= alpha * Ap
        rs_new = float(r @ r)
        p = r + (rs_new / rs) * p
        rs = rs_new
        it += 1
    return x, it, (rs ** 0.5) / b_norm


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
# main
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--checkpoint", required=True,
                   help="Stage-2 fit checkpoint (fit_best.pt).")
    p.add_argument("--shards", required=True, help="Arrow shard dir/file(s).")
    p.add_argument("--device", default="cuda:0" if torch.cuda.is_available()
                   else "cpu")
    p.add_argument("--batch-size", type=int, default=262144)
    p.add_argument("--eval-chunk", type=int, default=16384,
                   help="Events per HVP autograd chunk.")
    p.add_argument("--max-events-fit", type=int, default=0,
                   help="Cap the materialised stage-2 (fit) events; 0 = all.")
    p.add_argument("--max-events-flow", type=int, default=0,
                   help="Cap the materialised stage-1 (flow) events; 0 = all. "
                   "The covariance is rescaled by the full/seen Σw ratio.")
    p.add_argument("--ridge-w", type=float, default=1e-6,
                   help="RELATIVE damping for H_ww: λ_w = ridge·tr(H_ww)/n_w "
                   "(trace from Hutchinson probes at startup), so the value "
                   "is sample-size and unit independent.")
    p.add_argument("--ridge-flow", type=float, default=1e-4,
                   help="RELATIVE damping for H₁: λ_φ = ridge·tr(H₁)/n_φ — "
                   "the regulariser of the near-null (early-stopping-flat) "
                   "flow directions. SCAN this (e.g. ×10 up/down) and quote "
                   "the plateau.")
    p.add_argument("--grid-nphi", type=int, default=0,
                   help="θ-mlp only: uniform φ bins of the fixed (η,φ) output "
                   "grid (η = the stats η-bin centres); 0 = match the "
                   "checkpoint's output_fisher_nphi (default 4) so the "
                   "covariance combines with the --output-fisher file. The "
                   "columns are EXACT — cost is linear in n_η·n_φ·n_comp.")
    p.add_argument("--cg-tol", type=float, default=1e-4,
                   help="CG relative-residual tolerance.")
    p.add_argument("--cg-max-iter-w", type=int, default=200)
    p.add_argument("--cg-max-iter-flow", type=int, default=200)
    p.add_argument("--precision", choices=["fp32", "fp64"], default="fp64",
                   help="Compute precision (fp64 strongly recommended: CG "
                   "residuals live below fp32 HVP noise).")
    p.add_argument("--fisher", default=None,
                   help="Optional empirical_fisher.pt to combine with "
                   "(data-stat + flow totals written and printed).")
    p.add_argument("--output", default=None,
                   help="Output .pt (default: <ckpt dir>/flow_uncertainty.pt)")
    p.add_argument("--progress", action="store_true", default=True)
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
        scale_cols = ([c for c in range(3)
                       if float(model.scale_param_mask[c]) != 0.0]
                      if model.scale_enabled else [])
        smear_cols = (T._smear_active_cols(model)
                      if model.smearing_enabled else [])
    else:
        params, segs, labels, n_theta, smear_cols = _w_layout(model)
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
        m_window=(ft_lo, ft_hi), **common)
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

    def _raw_hvp_phi(v):
        vl = []
        off = 0
        for p in flow_params:
            vl.append(v[off:off + p.numel()].view_as(p).to(p.dtype))
            off += p.numel()
        h = L1.hvp(vl, flow_params)
        return torch.cat([x.reshape(-1) for x in h]).double()

    # Trace-scale estimates (Hutchinson, 2 Rademacher probes each) so the
    # ridge inputs are RELATIVE — independent of sample size and units.
    gen = torch.Generator(device="cpu").manual_seed(7)

    def _trace_scale(raw_hvp, n):
        acc = 0.0
        for _ in range(2):
            v = (torch.randint(0, 2, (n,), generator=gen,
                               dtype=torch.int64).double() * 2.0 - 1.0).to(dev)
            acc += float(v @ raw_hvp(v)) / n
        return acc / 2.0

    print("estimating Hessian trace scales (Hutchinson)...")
    sc_w = _trace_scale(_raw_hvp_w, n_w)
    sc_phi = _trace_scale(_raw_hvp_phi, n_phi)
    lam_w = args.ridge_w * abs(sc_w)
    lam_phi = args.ridge_flow * abs(sc_phi)
    print(f"  tr(H_ww)/n = {sc_w:.4e} → λ_w = {lam_w:.4e}; "
          f"tr(H₁)/n = {sc_phi:.4e} → λ_φ = {lam_phi:.4e}")

    def hvp_w(v):
        return _raw_hvp_w(v) + lam_w * v

    def hvp_phi(v):
        return _raw_hvp_phi(v) + lam_phi * v

    def mixed_u(x):
        xl = _unpack(x, params, segs)
        h = L2.mixed(xl, params, flow_params)
        return torch.cat([t.reshape(-1) for t in h]).double()

    # ---- the column loop ---------------------------------------------------
    print(f"running {n_theta} θ-columns "
          f"(λ_w={lam_w:.3e}, λ_φ={lam_phi:.3e}, "
          f"tol={args.cg_tol:g}, α₁={alpha1:.4f})...")
    U = torch.zeros((n_theta, n_phi), dtype=torch.float64)
    Z = torch.zeros((n_theta, n_phi), dtype=torch.float64)
    t0 = time.time()
    cg_stats = []
    for j in range(n_theta):
        if is_mlp:
            e = G[j]                    # delta-method RHS g_j = ∂T_j/∂w
        else:
            e = torch.zeros(n_w, dtype=torch.float64, device=dev)
            e[j] = 1.0
        x, it_w, res_w = _cg(hvp_w, e, args.cg_tol, args.cg_max_iter_w,
                             label=f"w:{labels[j]}")
        u = mixed_u(x)
        z, it_f, res_f = _cg(hvp_phi, u, args.cg_tol, args.cg_max_iter_flow,
                             label=f"φ:{labels[j]}")
        U[j] = u.cpu()
        Z[j] = z.cpu()
        cg_stats.append((it_w, res_w, it_f, res_f))
        el = time.time() - t0
        print(f"  [{j + 1:3d}/{n_theta}] {labels[j]:14s} "
              f"CG_w {it_w:3d} it (res {res_w:.1e})  "
              f"CG_φ {it_f:3d} it (res {res_f:.1e})  "
              f"elapsed {el / 60:.1f} min", flush=True)

    cov_flow = (U @ Z.T) / alpha1
    cov_flow = 0.5 * (cov_flow + cov_flow.T)          # numerical symmetrise
    asym = float((U @ Z.T - Z @ U.T).abs().max()
                 / cov_flow.abs().max().clamp_min(1e-300))
    print(f"covariance asymmetry (CG-residual scale check): {asym:.2e}")

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
