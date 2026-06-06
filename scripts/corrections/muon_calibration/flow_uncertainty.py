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

Scope: binned θ (the θ-mlp output-grid case needs the sketched variant to be
affordable and is rejected with a message). Flow archs with a plain NLL
stage-1 objective (compact/dcb/ege: −Σw·log p₀ on the FLOW window) and
gf/nsf (truncated window-norm + gauge penalty, replicating step1); nce's
BCE objective is not supported.

Output: ``flow_uncertainty.pt`` with the raw-θ covariance in the SAME active
layout as ``empirical_fisher.pt`` (combine by addition), the physical-units
extras via the shared ``_theta_cov_extras`` (softplus delta-method included),
and CG diagnostics. With ``--fisher <empirical_fisher.pt>`` the combined
(data + flow) covariance and a per-parameter σ budget table are also written.
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
    if getattr(model, "theta_mode", "binned") == "mlp":
        print("error: θ-mlp output-space propagation needs the sketched "
              "variant (exact columns over the output grid are not "
              "affordable); binned θ only for now.", file=sys.stderr)
        return 1
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
    params, segs, labels, n_theta, smear_cols = _w_layout(model)
    for p in params:
        p.requires_grad_(True)
    n_w = segs[-1][0] + segs[-1][1].numel()
    print(f"w-space: {n_w} active ({n_theta} θ + {n_w - n_theta} background); "
          f"φ-space: {n_phi} flow weights; arch={arch}")

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

    n_scale = model.theta_scale.numel() if model.scale_enabled else 0
    extras = T._theta_cov_extras(cov_flow, model, smear_cols, n_scale)
    out = {
        "covariance_flow": cov_flow,
        "labels": labels[:n_theta],
        "n_scale": n_scale,
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
    for k, v in extras.items():
        out[f"{k}_flow"] = v

    if args.fisher and os.path.exists(args.fisher):
        ef = torch.load(args.fisher, map_location="cpu", weights_only=False)
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
            print(f"warning: --fisher covariance shape {tuple(cov_d.shape)} "
                  f"does not match the θ-active block "
                  f"{tuple(cov_flow.shape)}; not combined", file=sys.stderr)

    out_path = args.output or os.path.join(
        os.path.dirname(os.path.abspath(args.checkpoint)),
        "flow_uncertainty.pt")
    torch.save(out, out_path)
    print(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
