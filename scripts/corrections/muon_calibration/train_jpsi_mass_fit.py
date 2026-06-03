"""Train the unbinned per-event J/ψ mass-fit calibration.

End-to-end driver:
  1. Discover Arrow shards (MC + data) produced by ``jpsi_mass_fit_snapshot.py``.
  2. Compute per-column standardisation stats over the full dataset.
  3. Build the :class:`JpsiMassMixtureModel` — a θ-conditioned flow with
     per-η-bin scale (A, e, M) and smearing (a, c) nuisances.
  4. Iterate the joint MLE: the MC branch forward-folds smear+scale at
     sampled, detached θ̃ to train the flow's conditional shape; the data
     branch fits θ through the flow's conditioning. Until early stop / epoch
     limit.
  5. Optionally, compute the plug-in (observed) Fisher information w.r.t.
     ``theta_scale`` at the converged point and persist the covariance.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from dataclasses import replace
from typing import List

import numpy as np
import torch
from tqdm import tqdm

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from jpsi_mass_arrow_loader import (
    JpsiMassArrowLoader,
    JpsiMassPreprocStats,
    compute_jpsi_mass_stats,
    discover_shards,
)
from jpsi_mass_model import (
    JpsiMassMixtureModel, SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C, THETA_SCALE_REF)


# ---------------------------------------------------------------------------
# Stats persistence
# ---------------------------------------------------------------------------


def _stats_to_dict(s: JpsiMassPreprocStats) -> dict:
    return {
        "mll_mean": float(s.mll_mean),
        "mll_std": float(s.mll_std),
        "y_event_mean": s.y_event_mean.tolist(),
        "y_event_std": s.y_event_std.tolist(),
        "muon_kin_mean": s.muon_kin_mean.tolist(),
        "muon_kin_std": s.muon_kin_std.tolist(),
        "eta_edges": s.eta_edges.tolist(),
        "m_lo": float(s.m_lo),
        "m_hi": float(s.m_hi),
        "k_moments": (None if s.k_moments is None else s.k_moments.tolist()),
    }


def _stats_from_dict(d: dict) -> JpsiMassPreprocStats:
    return JpsiMassPreprocStats(
        mll_mean=float(d["mll_mean"]),
        mll_std=float(d["mll_std"]),
        y_event_mean=np.asarray(d["y_event_mean"], dtype=np.float32),
        y_event_std=np.asarray(d["y_event_std"], dtype=np.float32),
        muon_kin_mean=np.asarray(d["muon_kin_mean"], dtype=np.float32),
        muon_kin_std=np.asarray(d["muon_kin_std"], dtype=np.float32),
        eta_edges=np.asarray(d["eta_edges"], dtype=np.float64),
        m_lo=float(d["m_lo"]),
        m_hi=float(d["m_hi"]),
        k_moments=(None if d.get("k_moments") is None
                   else np.asarray(d["k_moments"], dtype=np.float64)),
    )


# ---------------------------------------------------------------------------
# Per-batch loss helpers
# ---------------------------------------------------------------------------


def _make_amp(precision: str, device: str):
    """Return ``(autocast_ctx_factory, scaler)`` for ``--precision``.

    Matches the convention in ``train_muon_response_flow.py``:
      * fp32 → autocast disabled, no GradScaler.
      * bf16 → bfloat16 autocast, GradScaler disabled (bf16 has the
        same exponent range as fp32; loss scaling is unnecessary).
      * fp16 → float16 autocast + enabled GradScaler for loss scaling.
    """
    if device.startswith("cuda"):
        amp_device_type = "cuda"
    elif device.startswith("xpu"):
        amp_device_type = "xpu"
    else:
        amp_device_type = "cpu"
    if precision == "fp32":
        amp_dtype = torch.float32
        amp_enabled = False
    elif precision == "bf16":
        amp_dtype = torch.bfloat16
        amp_enabled = True
    elif precision == "fp16":
        amp_dtype = torch.float16
        amp_enabled = True
    else:
        raise ValueError(f"unknown precision {precision!r}")
    amp_ctx = lambda: torch.amp.autocast(
        device_type=amp_device_type, dtype=amp_dtype, enabled=amp_enabled,
    )
    scaler = torch.amp.GradScaler(
        amp_device_type, enabled=(precision == "fp16"),
    )
    return amp_ctx, scaler


def _move_batch(batch: dict, device: str) -> dict:
    return {
        k: v.to(device, non_blocking=device.startswith("cuda"))
        for k, v in batch.items()
    }


def _lr_str(optim: torch.optim.Optimizer) -> str:
    """Compact current-lr string for the progress bar: a single value if all
    param groups share an lr, else the per-group lrs joined by '/'."""
    lrs = [g["lr"] for g in optim.param_groups]
    if len(set(lrs)) == 1:
        return f"{lrs[0]:.2g}"
    return "/".join(f"{x:.2g}" for x in lrs)


def _make_scheduler(args, optim, epochs):
    """Build the LR scheduler per ``--lr-schedule``. ``plateau`` reduces lr by
    ``--lr-reduce-factor`` when the val metric stalls for ``--lr-reduce-patience``
    epochs (down to ``--min-lr``); ``cosine`` decays lr→min_lr over ``epochs``;
    ``none`` keeps it fixed. Returns ``(scheduler | None, kind)``."""
    kind = getattr(args, "lr_schedule", "none")
    if kind == "plateau":
        # threshold_mode="abs": detect a plateau by the SAME absolute NLL
        # decrease the early-stop uses (val < best - patience_threshold).
        # The default "rel" mode tests val < best*(1-threshold), which with a
        # NEGATIVE NLL (a log-density) moves the bar toward zero — i.e. a flat
        # or slightly-worse epoch still "improves" — so the LR would never
        # reduce while the absolute-threshold early-stop fires anyway.
        return torch.optim.lr_scheduler.ReduceLROnPlateau(
            optim, mode="min", factor=args.lr_reduce_factor,
            patience=args.lr_reduce_patience, min_lr=args.min_lr,
            threshold=args.patience_threshold, threshold_mode="abs"), kind
    if kind == "cosine":
        return torch.optim.lr_scheduler.CosineAnnealingLR(
            optim, T_max=max(1, epochs), eta_min=args.min_lr), kind
    return None, "none"


def _make_fit_optimizer(args, groups, minibatch_loop=False, kind=None):
    """Build the stage-2 optimizer over ``groups`` (list of {params, lr} dicts)
    per ``--fit-optimizer``. 'adam' (default): the historical torch.optim.Adam.
    'soap': SOAP (Shampoo-in-the-Adam-eigenbasis) from pytorch_optimizer — a
    Kronecker-factored curvature preconditioner that whitens the per-tensor
    gradient, so correlated/ill-scaled directions (the θ_net / background MLP
    weights especially) converge more completely than diagonal Adam toward the
    local minimum. Per-group lrs, the plateau scheduler, _lr_str, and the
    bootstrap reset all work unchanged (SOAP keeps the param_groups interface).
    'lbfgs': torch.optim.LBFGS with a strong-Wolfe line search — a quasi-Newton
    full-batch minimiser purpose-built for EXHAUSTIVE descent on a smooth
    deterministic objective (the stage-2 NLL is deterministic: fixed pseudo-data,
    seeded smear, all events summed). It builds a curvature model from the
    gradient history and line-searches to a tiny gradient norm in few outer
    iterations, handling the A/e, a/c ill-conditioning natively. REQUIRES the
    closure / full-batch loop in _run_epochs (driven there); torch.optim.LBFGS
    also rejects per-parameter groups, so the groups are flattened to one param
    list (the per-group lrs don't apply — the line search auto-scales the step;
    the O(1) THETA_SCALE_REF / SMEAR_VAR_SCALE rescaling still conditions θ).

    ``minibatch_loop=True`` (the bootstrap's per-replica refit, which steps without
    a closure) → 'lbfgs' falls back to Adam with a note, since LBFGS needs the
    closure loop.

    Note: SOAP preconditions WITHIN each parameter tensor, so it complements but
    does not replace --theta-whiten, whose analytic eigenbasis couples the small
    CROSS-tensor degeneracies (A/e, a/c, scale/smear) that live across the tiny
    binned-θ tensors. The two can be combined.

    'adam+lbfgs' / 'soap+lbfgs': a two-phase HYBRID (orchestrated in train_stage2,
    not here) — run Adam/SOAP to its normal stopping for robust bulk descent
    (Adam/SOAP escape the near-flat θ≈0 start that stalls L-BFGS's line search,
    but, normalising by the gradient RMS, plateau the NLL with a non-negligible
    gradient), THEN warm-start L-BFGS for --lbfgs-final-epochs to polish to a
    tight gradient norm (which Adam/SOAP structurally cannot reach). ``kind``
    lets the caller request a specific phase ('adam'/'soap'/'lbfgs') overriding
    --fit-optimizer."""
    if kind is None:
        kind = getattr(args, "fit_optimizer", "adam")
    # Hybrids ('adam+lbfgs', 'soap+trust-krylov', …) and the closure/scipy-driven
    # optimisers need their own driver loop; the minibatch bootstrap refit can
    # only use a plain step-based optimiser, so fall back to the base (adam/soap)
    # or adam.
    if minibatch_loop and ("+" in kind or kind in ("lbfgs", "trust-krylov",
                                                    "trust-ncg", "trust-exact")):
        base = kind.split("+")[0]
        kind = base if base in ("adam", "soap") else "adam"
        print(f"  note: --fit-optimizer {getattr(args, 'fit_optimizer', '')} "
              f"needs its own driver loop; the bootstrap refit uses {kind}.")
    if kind == "adam":
        return torch.optim.Adam(groups)
    if kind == "lbfgs":
        flat = [p for g in groups for p in g["params"]]
        return torch.optim.LBFGS(
            flat, lr=float(args.lbfgs_lr), max_iter=int(args.lbfgs_max_iter),
            history_size=int(args.lbfgs_history_size),
            line_search_fn="strong_wolfe",
            tolerance_grad=float(args.lbfgs_tolerance_grad),
            tolerance_change=float(args.lbfgs_tolerance_change))
    if kind == "soap":
        try:
            from pytorch_optimizer import SOAP
        except ImportError as e:
            raise RuntimeError(
                "--fit-optimizer soap requires the pytorch_optimizer package "
                f"(import failed: {e})")
        # weight_decay=0.0 (SOAP defaults to 0.01!) — a calibration fit must NOT
        # be pulled toward θ=0; this matches the Adam(groups) path (wd=0).
        # precondition_1d=True (SOAP defaults to False!) — else 1-D parameters
        # (every nn.Linear bias, incl. the θ_net final-layer (A,e,M,a,c) bias,
        # and any 1-D θ) are left UN-preconditioned / pure-Adam; we want SOAP to
        # condition them too. precondition_frequency: steps between the (cheap,
        # tiny-tensor here) eigendecompositions. eps: the denominator floor in
        # the rotated space — acts as a ridge on the preconditioner (larger →
        # less aggressive whitening of the sloppy/near-degenerate directions).
        sb = args.soap_shampoo_beta
        return SOAP(groups, weight_decay=0.0, precondition_1d=True,
                    eps=float(args.soap_eps),
                    shampoo_beta=(float(sb) if sb is not None and sb >= 0 else None),
                    precondition_frequency=int(args.soap_precondition_frequency))
    raise ValueError(f"unknown --fit-optimizer {kind!r}")


class _StepProfiler:
    """Phase timer to diagnose whether the training step is data/host-bound or
    compute-bound (and, on CUDA, GPU-bound vs launch/sync-bound). Times three
    phases over the first ``n`` steps of an epoch:

      • data   — loader fetch + ``_move_batch`` host→device copy (the gap from the
                 end of the previous step to the start of forward);
      • fwd    — ``step_fn`` up to the loss;
      • bwd    — ``loss.backward()`` + ``optim.step()``.

    On CUDA each phase is bracketed by ``torch.cuda.synchronize()`` so the async
    kernels are attributed to the phase that launched them (else wall-clock would
    measure only Python launch time). Prints the mean split + GPU util/mem when a
    GPU is active. Reading: data≫compute → CPU/dataloader-bound; compute high +
    GPU util high → GPU-bound; compute high + GPU util low → launch/sync-bound
    (small batches, the per-event fixed-point / GH-quadrature Python overhead)."""

    def __init__(self, device, n):
        self.n = int(n)
        self.is_cuda = (str(device).startswith("cuda")
                        and torch.cuda.is_available())
        self.data = self.fwd = self.bwd = 0.0
        self.count = 0
        self._t = None
        self._pd = self._pf = 0.0   # pending (this-iter) data / fwd, committed at backward

    def _sync(self):
        if self.is_cuda:
            torch.cuda.synchronize()

    def active(self):
        return self.count < self.n

    def mark_iter_start(self):
        """Call at the top of the loop body (a batch has just been yielded)."""
        if self.active():
            self._pd = self._pf = 0.0
            self._sync(); self._t = time.time()

    def after_move(self):
        if self.active():
            self._sync(); now = time.time(); self._pd = now - self._t; self._t = now

    def after_forward(self):
        if self.active():
            self._sync(); now = time.time(); self._pf = now - self._t; self._t = now

    def after_backward(self):
        # Commit all three phases together so bailed steps (sw<=0 / NaN-skip,
        # which never reach here) don't leak a partial iter into the averages.
        if self.active():
            self._sync(); now = time.time()
            self.data += self._pd; self.fwd += self._pf; self.bwd += now - self._t
            self.count += 1

    def report(self, stage_name, epoch):
        if self.count == 0:
            return
        n = self.count
        d, f, b = self.data / n, self.fwd / n, self.bwd / n
        tot = d + f + b
        if tot <= 0:
            return
        msg = (f"[{stage_name}] epoch {epoch:>3} profile (mean of {n} steps): "
               f"data={d*1e3:.1f}ms ({100*d/tot:.0f}%)  "
               f"fwd={f*1e3:.1f}ms ({100*f/tot:.0f}%)  "
               f"bwd+step={b*1e3:.1f}ms ({100*b/tot:.0f}%)")
        if self.is_cuda:
            try:
                util = torch.cuda.utilization()
                mem = torch.cuda.max_memory_allocated() / 1024**3
                msg += f"  | GPU util≈{util}%  peak_mem={mem:.2f}GB"
            except Exception:
                pass
        else:
            msg += "  | device=cpu (no GPU in use)"
        print(msg)


def _validation_half(args, which: str) -> "int | None":
    """Which event half to use in validation mode for the named stage
    (``which`` ∈ {'flow', 'fit'}). Defaults to the disjoint half split
    (flow ← 0, fit ← 1) so the fit pseudo-data is statistically independent
    from the flow's training sample. With ``--no-validation-split`` returns
    ``None`` for both stages — all events are used for both. Returns
    ``None`` outside validation mode."""
    if not getattr(args, "validation", False):
        return None
    if getattr(args, "no_validation_split", False):
        return None
    return 0 if which == "flow" else 1


def _smear_active_cols(model) -> List[int]:
    """Column indices of ``theta_smear`` actually fit (smear_param_mask != 0).
    The non-fit column is zeroed post-softplus → identically zero gradient, so
    it must be excluded from the Fisher Hessian (else a singular row/col)."""
    return [c for c in range(model.theta_smear.shape[1])
            if float(model.smear_param_mask[c]) != 0.0]


def _active_param_labels(model, smear_cols) -> List[str]:
    """Per-active-parameter labels in the [θ_scale | active θ_smear] order
    shared by the Fisher and bootstrap covariances."""
    comp = ["A", "e", "M"]
    smear_comp = ["a", "c"]
    labels: List[str] = []
    if model.scale_enabled:
        for b in range(model.theta_scale.shape[0]):
            for c in range(model.theta_scale.shape[1]):
                labels.append(f"{comp[c]}[{b}]")
    for b in range(model.theta_smear.shape[0]):
        for c in smear_cols:
            labels.append(f"{smear_comp[c]}[{b}](raw)")
    return labels


def _record_active_theta(model, smear_cols) -> torch.Tensor:
    """Flat active parameter vector ``[θ_scale (linear) | active θ_smear (raw)]``
    on cpu, matching the Fisher's active ordering (bin-major, col-minor)."""
    parts = []
    if model.scale_enabled:
        parts.append(model.theta_scale.detach().reshape(-1))
    if smear_cols:
        parts.append(model.theta_smear.detach()[:, smear_cols].reshape(-1))
    return torch.cat(parts).detach().cpu()


def _hessian_block_loop(g_full, params, active_idx, n_act):
    """Per-batch active Hessian block ``[n_act, n_act]`` by one second backward
    per row (robust; works through the nested autograd in the continuity
    density). ``g_full`` is the flat first-order gradient with its graph kept."""
    H = torch.zeros((n_act, n_act), device=g_full.device, dtype=g_full.dtype)
    for r in range(n_act):
        i = int(active_idx[r])
        row = torch.autograd.grad(
            g_full[i], params, retain_graph=True, allow_unused=True)
        row_full = torch.cat([
            (ri if ri is not None else torch.zeros_like(p)).reshape(-1)
            for ri, p in zip(row, params)])
        H[r] = row_full[active_idx]
    return H


def _hessian_block_batched(g_active, params, active_idx, n_act, chunk=None):
    """Per-batch active Hessian block ``[n_act, n_act]`` via ``is_grads_batched``
    (vmaps the per-row vjp over the identity basis), in row CHUNKS of ``chunk``
    (default: all ``n_act`` at once). Vectorising ALL rows holds ~n_act copies of
    the backward graph — fine for the small binned/output table (n_act ~ tens)
    but catastrophic when the background MLP is marginalised in (n_act ~ 1000s),
    so the caller passes a chunk to bound memory. The engine's vmap may not
    support every op in this double-backward path (the inner ``autograd.grad`` of
    the change-of-variables Jacobian, the fixed-point clamps) — the caller falls
    back to the per-row loop on failure/OOM."""
    c = int(chunk) if chunk else n_act
    H = torch.zeros((n_act, n_act), device=g_active.device, dtype=g_active.dtype)
    for r0 in range(0, n_act, max(1, c)):
        r1 = min(r0 + c, n_act)
        nb = r1 - r0
        eye = torch.zeros((nb, n_act), device=g_active.device, dtype=g_active.dtype)
        eye[torch.arange(nb, device=g_active.device),
            torch.arange(r0, r1, device=g_active.device)] = 1.0
        rows = torch.autograd.grad(
            g_active, params, grad_outputs=eye,
            is_grads_batched=True, retain_graph=True, allow_unused=True)
        parts = []
        for ri, p in zip(rows, params):
            if ri is None:
                parts.append(torch.zeros(nb, p.numel(),
                                         device=g_active.device, dtype=g_active.dtype))
            else:
                parts.append(ri.reshape(nb, -1))
        row_full = torch.cat(parts, dim=1)        # [nb, n_full]
        H[r0:r1] = row_full[:, active_idx]        # [nb, n_act]
    return H


def compute_fisher_info_continuity(
    model: JpsiMassMixtureModel,
    loader: JpsiMassArrowLoader,
    device: str,
    *,
    mc_as_data: bool = False,
    n_iter: int = 2,
    progress: bool = True,
    vectorized: bool = True,
):
    """Observed (plug-in) Fisher information for the two-stage continuity fit,
    over ``theta_scale`` + the ACTIVE ``theta_smear`` columns jointly, with the
    flow and the background MLP held FIXED (option 1: conditional / fixed-φ).

    H = Σ_events w · ∂²(−ln p_mixture)/∂θ² at the fit point, on the data branch
    (``data_nll_continuity`` — the objective stage 2 actually minimises).
    ``theta_smear`` are the signed width coefficients.

    Accumulated per batch (build the batch gradient with ``create_graph`` then
    immediately take its ∂/∂θ rows and free the graph) so memory stays at one
    batch regardless of dataset size; cost is O(N_batches · n_active) backward
    passes. With ``mc_as_data`` the data branch is the MC (``~is_data_mask``)
    rows — the validation-mode pseudo-data.

    Returns ``(H [n_act, n_act] on cpu, layout dict)``.
    """
    model.eval()
    for p in model.flow.parameters():
        p.requires_grad_(False)
    for p in model.mlp.parameters():
        p.requires_grad_(False)
    model.theta_scale.requires_grad_(model.scale_enabled)
    model.theta_smear.requires_grad_(model.smearing_enabled)

    # Parameters to differentiate + the active flat layout. theta_scale: all
    # active; theta_smear: only the fitted column(s).
    params: list = []
    blocks: list = []          # (name, numel, active_local_indices)
    if model.scale_enabled:
        params.append(model.theta_scale)
        blocks.append(("scale", model.theta_scale.numel(),
                       list(range(model.theta_scale.numel()))))
    smear_cols = _smear_active_cols(model) if model.smearing_enabled else []
    if smear_cols:
        params.append(model.theta_smear)
        n_eta, n_comp = model.theta_smear.shape
        active = [b * n_comp + c for b in range(n_eta) for c in smear_cols]
        blocks.append(("smear", model.theta_smear.numel(), active))
    if not params:
        raise RuntimeError(
            "compute_fisher_info_continuity: no free parameters "
            "(--disable-scale and --disable-smearing / no active smear term).")

    # Active index into the concatenated flat [scale_flat | smear_flat] vector.
    active_idx, offset = [], 0
    for _name, numel, act in blocks:
        active_idx += [offset + a for a in act]
        offset += numel
    active_idx = torch.tensor(active_idx, dtype=torch.long, device=device)
    n_act = int(active_idx.numel())

    H = torch.zeros((n_act, n_act), device=device, dtype=torch.float32)
    grad = torch.zeros(n_act, device=device, dtype=torch.float32)  # Σ ∂(NLL)/∂θ
    sw = 0.0
    seen = 0
    use_batched = bool(vectorized)  # may flip to False after a fallback
    bar = tqdm(loader, desc="fisher", disable=not progress, unit="batch")
    for batch in bar:
        batch = _move_batch(batch, device)
        data_mask = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
        if not bool(data_mask.any()):
            continue
        per = model.data_nll_continuity(
            batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
            batch["q_pm"], batch["b_pm"], batch["cond_std"], data_mask,
            n_iter=n_iter)
        w = batch["w"] * data_mask.to(batch["w"].dtype)
        nll = (w * per).sum()
        if not torch.isfinite(nll):
            continue
        g = torch.autograd.grad(nll, params, create_graph=True)
        g_full = torch.cat([gi.reshape(-1) for gi in g])
        g_active = g_full[active_idx]
        grad += g_active.detach()
        # Vectorised second backward (one vmapped vjp over the identity basis);
        # fall back to the per-row loop on any engine failure / OOM, once.
        if use_batched:
            try:
                Hb = _hessian_block_batched(g_active, params, active_idx, n_act)
            except (RuntimeError, NotImplementedError) as e:
                use_batched = False
                if device.startswith("cuda"):
                    torch.cuda.empty_cache()
                bar.write(
                    f"  note: vectorised Hessian unavailable "
                    f"({type(e).__name__}: {str(e).splitlines()[0][:80]}); "
                    f"using the per-row loop")
                Hb = _hessian_block_loop(g_full, params, active_idx, n_act)
        else:
            Hb = _hessian_block_loop(g_full, params, active_idx, n_act)
        H += Hb.detach()
        sw += float(w.sum().item())
        seen += int(data_mask.sum().item())
        bar.set_postfix_str(f"events={seen:,}")
    bar.close()
    if seen == 0:
        raise RuntimeError(
            "compute_fisher_info_continuity: loader yielded zero events on the "
            "data branch (mc_as_data=%s). Check the split / --validation." % mc_as_data)
    H = 0.5 * (H + H.T)
    layout = {"blocks": blocks, "smear_cols": smear_cols,
              "n_scale": (model.theta_scale.numel() if model.scale_enabled else 0),
              "sw": sw, "seen": seen, "grad": grad.detach().cpu()}
    return H.detach().cpu(), layout


def _fisher_save_dict(H: torch.Tensor, layout: dict, model: JpsiMassMixtureModel) -> dict:
    """Invert H → covariance and package it (with labels, the θ_scale block in
    the 24×3×24×3 layout for the diagnostics, and delta-method effective
    σ for the raw θ_smear)."""
    n_act = H.shape[0]
    # Positive-definiteness check: at a true optimum the observed information is
    # PD. Negative/zero eigenvalues flag a non-converged fit or unidentified
    # (e.g. event-starved) η-bins; the corresponding variances are not
    # trustworthy. eigvalsh is exact for the symmetric H.
    try:
        eig = torch.linalg.eigvalsh(H)
        min_eig = float(eig.min())
        n_neg_eig = int((eig <= 0).sum())
    except RuntimeError:
        min_eig = float("nan")
        n_neg_eig = -1
    cov = None
    ok = False
    try:
        cov = torch.linalg.inv(H)
        ok = bool(torch.isfinite(cov).all()) and n_neg_eig == 0
    except RuntimeError:
        ok = False
    if cov is None or not bool(torch.isfinite(cov).all()):
        try:
            cov = torch.linalg.pinv(H)
        except Exception:
            cov = None
    # Same fix as in _empirical_cov_theta_block: any parameter with H_ii = 0
    # (frozen via mask → bit-zero per-event gradient → bit-zero diagonal) is
    # not fit, so its variance is 0 by definition. Zero the corresponding
    # rows/cols in `cov` to remove any inv/pinv-amplified garbage that would
    # otherwise poison the chi² tolerance and the correlation heatmap.
    if cov is not None:
        zero_diag = (torch.diag(H) == 0)
        if zero_diag.any():
            cov[zero_diag, :] = 0.0
            cov[:, zero_diag] = 0.0

    n_scale = layout["n_scale"]
    smear_cols = layout["smear_cols"]
    labels = _active_param_labels(model, smear_cols)

    out: dict = {
        "hessian": H,
        "covariance": cov,
        "ok": ok,
        "labels": labels,
        "n_scale": n_scale,
        "smear_cols": smear_cols,
        "smear_fit_params": model.smear_fit_params,
        "smear_param_form": getattr(model, "smear_param_form", "linear"),
        "param_space": "scale: linear (A,e,M); smear: " + getattr(model, "smear_param_form", "linear"),
        "n_events": layout["seen"],
        "sum_weight": layout["sw"],
        "min_eig": min_eig,
        "n_negative_eig": n_neg_eig,
    }
    # Estimated distance to minimum: EDM = ½ gᵀ V g (V = covariance = H⁻¹), the
    # predicted remaining decrease in the NLL to reach the optimum (MINUIT
    # convention). ≈0 at convergence; large ⇒ the fit hasn't reached a minimum.
    grad = layout.get("grad")
    if grad is not None:
        out["grad"] = grad
        out["grad_norm"] = float(grad.norm())
        if cov is not None:
            out["edm"] = 0.5 * float(grad @ (cov @ grad))
    # θ_scale block in the 24×3×24×3 layout (consumed by the diagnostics for ±1σ
    # bands + the correlation heatmap). This is the scale block of the JOINT
    # inverse, so it carries the θ_scale↔θ_smear correlation.
    n_eta_s = model.theta_scale.shape[0]          # generic η-bin count (was 24)
    if n_scale == 3 * n_eta_s:
        out["hessian_24_3_24_3"] = H[:n_scale, :n_scale].view(n_eta_s, 3, n_eta_s, 3)
        if cov is not None:
            # θ_scale is O(1); physical (A,e,M) = θ·THETA_SCALE_REF → rescale the
            # covariance block to PHYSICAL units (cov_phys = cov·ref⊗ref).
            refN = torch.tensor(list(THETA_SCALE_REF) * n_eta_s, dtype=cov.dtype)
            cov_scale = cov[:n_scale, :n_scale] * refN.unsqueeze(0) * refN.unsqueeze(1)
            out["covariance_24_3_24_3"] = cov_scale.view(n_eta_s, 3, n_eta_s, 3)
            out["sigma_scale_24_3"] = torch.sqrt(
                torch.clamp(torch.diag(cov_scale), min=0.0)).view(n_eta_s, 3)
    # PHYSICAL σ on the smear a/c per η-bin. The fit param θ_smear is O(1); the
    # physical qop-variance coefficient is `effective(θ)·SMEAR_VAR_SCALE`, where
    # `effective` is the identity ('linear'), `softplus` ('softplus'), or `θ²`
    # ('square'). Delta method: σ_phys = |d effective/dθ|·SMEAR_VAR_SCALE·σ_raw —
    # Jacobian = 1 (linear), sigmoid(θ̂) (softplus), |2·θ̂| (square; → 0 at θ̂=0,
    # the documented raw-θ singularity for a bin pinned at zero — use
    # --output-fisher there). Evaluated at the fitted θ. Inactive columns → 0.
    if smear_cols and cov is not None:
        n_eta, n_comp = model.theta_smear.shape
        cov_smear = cov[n_scale:, n_scale:]
        sig_raw = torch.sqrt(torch.clamp(torch.diag(cov_smear), min=0.0))
        smear_scale = (SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C)
        form = getattr(model, "smear_param_form", "linear")
        with torch.no_grad():
            th = model.theta_smear.detach().cpu()
            if form == "softplus":
                jac = torch.sigmoid(th)
            elif form == "square":
                jac = (2.0 * th).abs()
            else:
                jac = torch.ones_like(th)
        sig_eff = torch.zeros(n_eta, n_comp)
        k = 0
        for b in range(n_eta):
            for c in smear_cols:
                sig_eff[b, c] = jac[b, c] * smear_scale[c] * sig_raw[k]
                k += 1
        out["sigma_smear_eff_24_2"] = sig_eff
    return out


# ---------------------------------------------------------------------------
# Main training loop
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Two-stage continuity training (default)
#   stage 1 ("flow"): nominal flow p₀(m|muon_kin) at θ=0 on simulation, no θ
#                     conditioning, MC rows only.
#   stage 2 ("fit"):  freeze the flow, fit θ_scale / θ_smear / background-MLP on
#                     data via the analytic continuity tilt (model.data_nll_continuity).
# ---------------------------------------------------------------------------


def _setup_common(args, *, stats_override=None):
    """Shards + stats + train/val loaders (shared by both stages).

    ``stats_override`` (used by ``--stage fit``) supplies the standardisation
    stats from the flow checkpoint, taking precedence over ``--stats-in`` /
    recomputation so the data branch standardises exactly as the flow was
    trained.
    """
    shard_files = discover_shards(args.inputs)
    if not shard_files:
        print("error: no Arrow shards found under --inputs", file=sys.stderr)
        return None
    print(f"discovered {len(shard_files)} shard(s)")
    # Resolve the binned-θ η binning. The flags default to None so we can tell an
    # explicit request from the 24/±2.4 default: when stats come from a checkpoint
    # or --stats-in, an explicit --n-eta-bins/--eta-range rebuilds stats.eta_edges
    # (the binning is independent of the flow), otherwise the loaded edges are kept.
    eta_bins_req = args.n_eta_bins is not None or args.eta_range is not None
    n_eta_bins = int(args.n_eta_bins) if args.n_eta_bins is not None else 24
    eta_range = float(args.eta_range) if args.eta_range is not None else 2.4
    if stats_override is not None:
        stats = stats_override
        print("using preproc stats from the flow checkpoint")
    elif args.stats_in is not None and os.path.exists(args.stats_in):
        with open(args.stats_in) as f:
            stats = _stats_from_dict(json.load(f))
        print(f"loaded preproc stats from {args.stats_in}")
    else:
        t0 = time.time()
        eta_edges = np.linspace(-eta_range, eta_range, n_eta_bins + 1, dtype=np.float64)
        stats = compute_jpsi_mass_stats(shard_files, m_lo=args.m_lo, m_hi=args.m_hi,
                                        eta_edges=eta_edges)
        print(f"computed preproc stats in {time.time() - t0:.1f}s "
              f"({n_eta_bins} η-bins over ±{eta_range:g})")
        eta_bins_req = False  # freshly computed with the requested binning already
    if eta_bins_req:
        eta_edges = np.linspace(-eta_range, eta_range, n_eta_bins + 1, dtype=np.float64)
        stats = replace(stats, eta_edges=eta_edges)
        print(f"  overriding binned-θ η binning: {n_eta_bins} bins over ±{eta_range:g} "
              "(flow conditioning is unaffected)")
    print(f"  mll: μ={stats.mll_mean:.4f}  σ={stats.mll_std:.4f}  "
          f"window [{stats.m_lo}, {stats.m_hi}]")
    os.makedirs(args.output, exist_ok=True)
    stats_path = os.path.join(args.output, "preproc_stats.json")
    with open(stats_path, "w") as f:
        json.dump(_stats_to_dict(stats), f, indent=2)
    print(f"wrote {stats_path}")
    train_loader, val_loader = _make_loaders(args, shard_files, stats)
    return shard_files, stats, train_loader, val_loader


def _make_loaders(args, shard_files, stats, *, half=None, inject_theta=None,
                  inject_smear=None, val_fraction=None, holdout_fraction=None):
    """Build the ``(train, val)`` loaders for one stage. ``half`` selects a
    deterministic disjoint event half (0/1) — used by the MC-closure
    validation mode (stage 1 ← half 0, stage 2 ← half 1); ``None`` = all
    events. ``inject_theta`` ([n_eta, 3]) / ``inject_smear`` ([n_eta, 2]) inject
    a known θ_scale shift / per-muon qop smear into the (pseudo-)data m_ll
    (validation closure). ``val_fraction``/``holdout_fraction`` override the
    args defaults (the fit stage passes 0/0 to use ALL events — within the
    half — with no held-out split)."""
    seed = int(getattr(args, "inject_smear_seed", 12345))
    me = int(getattr(args, "max_events", 0) or 0)
    ef = float(getattr(args, "event_fraction", 1.0) or 1.0)
    nu = bool(getattr(args, "inject_nonuniform", False))
    vf = args.val_fraction if val_fraction is None else float(val_fraction)
    hf = args.holdout_fraction if holdout_fraction is None else float(holdout_fraction)
    train_loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=args.batch_size, split="train",
        val_fraction=vf, holdout_fraction=hf,
        drop_last=True, half=half, inject_theta_scale=inject_theta,
        inject_theta_smear=inject_smear, inject_seed=seed,
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=me, event_fraction=ef, inject_nonuniform=nu)
    val_loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=args.batch_size, split="val",
        val_fraction=vf, holdout_fraction=hf,
        drop_last=False, half=half, inject_theta_scale=inject_theta,
        inject_theta_smear=inject_smear, inject_seed=seed,
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=me, event_fraction=ef, inject_nonuniform=nu)
    return train_loader, val_loader


def _inject_theta_np(args, n_eta):
    """``[n_eta, 3]`` injected θ_scale (constant A, e, M across η) for the
    validation closure, or ``None`` if no shift was requested."""
    a = float(getattr(args, "inject_A", 0.0) or 0.0)
    e = float(getattr(args, "inject_e", 0.0) or 0.0)
    m = float(getattr(args, "inject_M", 0.0) or 0.0)
    if a == 0.0 and e == 0.0 and m == 0.0:
        return None
    t = np.zeros((int(n_eta), 3), dtype=np.float64)
    t[:, 0] = a; t[:, 1] = e; t[:, 2] = m
    return t


def _inject_smear_np(args, n_eta):
    """``[n_eta, 2]`` injected PHYSICAL qop-variance coefficients (a, c) for the
    validation closure (or ``None``): σ²_qop = a + c·k², the values the loader
    applies directly. ``--inject-a/-c`` are the physical coefficients (a ~ 1e-7,
    c ~ 1e-6 are typical); the O(1) optimizer rescaling (SMEAR_VAR_SCALE) is an
    internal detail applied to the fit parameter, not here."""
    a = float(getattr(args, "inject_a", 0.0) or 0.0)
    c = float(getattr(args, "inject_c", 0.0) or 0.0)
    if a == 0.0 and c == 0.0:
        return None
    t = np.zeros((int(n_eta), 2), dtype=np.float64)
    t[:, 0] = a; t[:, 1] = c
    return t


# Args that fix the flow's parameter shapes OR conditioning semantics — these
# must match the saved flow when reloading it for --stage fit (the model is
# rebuilt from args before the flow weights are loaded). The MLP size and
# θ-smear fit choice are NOT here: only the flow is loaded from the checkpoint,
# so those stay free stage-2 knobs. cond_basis is included because a flow
# trained on one basis is meaningless under the other (same shape, different
# meaning), so stage-2 must inherit the flow's basis.
_FLOW_ARCH_KEYS = (
    "flow_arch", "flow_n_transforms", "flow_hidden", "flow_n_hidden",
    "gf_components", "nsf_bins", "cond_basis",
)


def _apply_flow_arch_from_ckpt(args, ck_args: dict) -> None:
    """Override the flow-architecture args on ``args`` with the values stored
    in the flow checkpoint so the rebuilt flow matches the saved weights.
    Logs any field that changed; silently keeps the CLI value for keys the
    checkpoint doesn't carry (older checkpoints)."""
    changed = []
    for k in _FLOW_ARCH_KEYS:
        if k in ck_args:
            old = getattr(args, k, None)
            new = ck_args[k]
            if old != new:
                changed.append(f"{k}: {old}→{new}")
            setattr(args, k, new)
    if changed:
        print("  flow-architecture args set from checkpoint: " + ", ".join(changed))
    else:
        print("  flow-architecture args already match the checkpoint")


def _build_model(args, stats, device):
    return JpsiMassMixtureModel(
        m_lo=stats.m_lo, m_hi=stats.m_hi, mll_log_scale=stats.mll_log_scale,
        mll_mean=stats.mll_mean, mll_std=stats.mll_std,
        y_event_mean=torch.from_numpy(stats.y_event_mean),
        y_event_std_tensor=torch.from_numpy(stats.y_event_std),
        muon_kin_mean=torch.from_numpy(stats.muon_kin_mean),
        muon_kin_std_tensor=torch.from_numpy(stats.muon_kin_std),
        flow_arch=args.flow_arch, flow_n_transforms=args.flow_n_transforms,
        flow_hidden_features=args.flow_hidden, flow_n_hidden_layers=args.flow_n_hidden,
        flow_gf_components=args.gf_components, flow_nsf_bins=args.nsf_bins,
        mlp_hidden=args.mlp_hidden, mlp_n_layers=args.mlp_n_layers,
        smearing_enabled=not args.disable_smearing,
        scale_enabled=not args.disable_scale,
        qop_floor_frac=args.qop_floor_frac, smear_fit_params=args.smear_fit_params,
        scale_fit_params=getattr(args, "scale_fit_params", "AM"),
        smear_flow_steps=getattr(args, "smear_flow_steps", 1),
        smear_operator=getattr(args, "smear_operator", "pf_ode"),
        n_gh_nodes=getattr(args, "n_gh_nodes", 8),
        jacobian_form=getattr(args, "jacobian_form", "softlog"),
        smear_param_form=getattr(args, "smear_param_form", "linear"),
        norm_correction=getattr(args, "norm_correction", "none"),
        background_enabled=not getattr(args, "no_background", False),
        theta_mode=("mlp" if getattr(args, "theta_mlp", False) else "binned"),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        theta_mlp_hidden=getattr(args, "theta_mlp_hidden", 32),
        theta_mlp_layers=getattr(args, "theta_mlp_layers", 2),
        theta_whiten=getattr(args, "theta_whiten", False),
        theta_whiten_max_rho=getattr(args, "theta_whiten_max_rho", 0.99),
        k_moments=stats.k_moments,
        # binned θ tables must match the stats' η binning (b_pm is bucketized
        # against stats.eta_edges in the loader) — derive from the edges, not the
        # N_ETA_BINS default, so --n-eta-bins (or a --stats-in with different
        # binning) is honoured.
        n_eta_bins=len(stats.eta_edges) - 1,
    ).to(device)


def _run_epochs(args, model, optim, train_loader, val_loader, stats, *,
                step_fn, ckpt_prefix, stage_name, epochs, monitor="val"):
    """Generic weighted-NLL epoch loop with best/last checkpoints + early-stop.
    ``step_fn(model, batch) -> (loss, sum_w)`` returns the batch's weighted-mean
    NLL (scalar tensor) over the rows the stage uses and the corresponding Σw.

    ``monitor`` selects the metric the early-stop / plateau-LR / best-checkpoint
    act on: ``"val"`` (the held-out validation NLL, stage 1) or ``"train"`` (the
    training NLL itself — stage 2 uses ALL events with no held-out split, so the
    train-NLL plateau is the convergence signal; the val pass is skipped). The
    per-epoch line always reports the train NLL and its change from the previous
    epoch. Returns the best monitored metric."""
    best_val = float("inf"); no_improve = 0; prev_train_nll = None
    best_ckpt = os.path.join(args.output, f"{ckpt_prefix}_best.pt")
    last_ckpt = os.path.join(args.output, f"{ckpt_prefix}_last.pt")
    device = args.device
    # The streaming loader has no __len__; learn the batch count on epoch 1 so
    # epochs ≥2 show a true percentage-complete bar.
    n_batches_total = None
    # L-BFGS governs its own step via the line search; an external LR scheduler
    # would fight it, so disable scheduling for L-BFGS (early-stop still applies).
    if isinstance(optim, torch.optim.LBFGS):
        sched, sched_kind = None, "none"
    else:
        sched, sched_kind = _make_scheduler(args, optim, epochs)
    if sched_kind != "none":
        print(f"  lr schedule: {sched_kind}"
              + (f" (factor={args.lr_reduce_factor:g}, patience={args.lr_reduce_patience}, "
                 f"min_lr={args.min_lr:g}); lr reductions tick in parallel with "
                 f"early-stop — they do NOT reset the early-stop counter, so "
                 f"`--patience {args.patience}` is the total stalled epochs "
                 f"allowed since the last improvement regardless of lr drops"
                 if sched_kind == "plateau"
                 else f" (T_max={epochs}, min_lr={args.min_lr:g})"))

    def _ckpt(epoch, bv, vm):
        return {
            "epoch": epoch, "stage": stage_name,
            "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()},
            "theta_scale": model.theta_scale.detach().cpu(),
            "theta_smear": model.theta_smear.detach().cpu(),
            "stats": _stats_to_dict(stats), "best_val": bv, "val_metric": vm,
            "args": vars(args),
        }

    prof_steps = int(getattr(args, "profile_steps", 0) or 0)
    is_lbfgs = isinstance(optim, torch.optim.LBFGS)
    for epoch in range(1, epochs + 1):
        t0 = time.time(); model.train()
        tr_sum = 0.0; tr_w = 0.0; n_seen = 0
        lr_str = _lr_str(optim)
        epoch_gnorm = None   # Σw-weighted gradient norm of the mean NLL (convergence)
        prof = _StepProfiler(device, prof_steps) if prof_steps > 0 else None
        if is_lbfgs:
            # Full-batch closure: one optim.step(closure) per "epoch" runs up to
            # --lbfgs-max-iter quasi-Newton inner iterations, each re-evaluating
            # the WEIGHTED-MEAN NLL + grad over the WHOLE loader (the objective is
            # deterministic, so this is a fixed smooth function of θ). The strong-
            # Wolfe line search guarantees descent; convergence is to a tiny
            # gradient norm, not a plateau heuristic.
            closure_stats = {"call": 0}
            lbfgs_params = [p for g in optim.param_groups for p in g["params"]]

            def closure():
                # GRADIENT ACCUMULATION: backward EACH batch and free its graph
                # immediately, accumulating only the scalar Σw-weighted total —
                # never holding the whole-dataset autograd graph at once (doing so
                # exhausts GPU memory → the CUDA allocator NVML assert). Σw is
                # θ-independent (fixed weights × mask), so grad(Σw·nll/Σw) =
                # (Σ_batch grad of Σw_batch·meanNLL_batch)/Σw — i.e. sum the
                # per-batch grads then rescale by 1/Σw at the end.
                optim.zero_grad(set_to_none=True)
                closure_stats["call"] += 1
                call = closure_stats["call"]
                s = 0.0; w = 0.0; nb = 0
                # Per-pass progress: each closure call is one full-batch
                # eval/line-search probe (up to --lbfgs-max-iter per epoch); the
                # bar tracks batches within the pass + running mean NLL.
                cbar = tqdm(train_loader, total=n_batches_total,
                            desc=f"[{stage_name}] ep{epoch:>3} pass {call}",
                            leave=False, disable=not args.progress, unit="batch")
                for batch in cbar:
                    batch = _move_batch(batch, device)
                    loss, sw = step_fn(model, batch)   # weighted-MEAN NLL over batch
                    nb += 1
                    if sw <= 0 or not torch.isfinite(loss):
                        continue
                    (loss * sw).backward()             # accumulates into .grad; frees graph
                    s += float(loss.item()) * sw; w += sw
                    cbar.set_postfix_str(f"nll={s / max(w, 1e-30):+.4f}")
                cbar.close()
                if w <= 0:
                    raise RuntimeError(
                        f"L-BFGS closure: no usable batches [{stage_name}] "
                        f"epoch {epoch}")
                inv_w = 1.0 / w
                for p in lbfgs_params:                 # rescale Σw·Σ → weighted mean
                    if p.grad is not None:
                        p.grad.mul_(inv_w)
                gnorm = sum(float((p.grad.detach()**2).sum())
                            for p in lbfgs_params if p.grad is not None) ** 0.5
                closure_stats["s"] = s; closure_stats["w"] = w; closure_stats["nb"] = nb
                # Always-visible one-liner per pass (survives --no-progress): the
                # full-batch passes are slow, so report mean NLL + |g| as each lands.
                print(f"  [{stage_name}] ep{epoch:>3} pass {call:>2}: "
                      f"nll={s * inv_w:+.5f}  |g|={gnorm:.3e}", flush=True)
                # L-BFGS only needs the scalar objective value (it reads .grad from
                # the params); return it as a detached tensor on the right device.
                return torch.as_tensor(s * inv_w, dtype=torch.float32, device=device)

            optim.step(closure)
            tr_sum = closure_stats["s"]; tr_w = closure_stats["w"]; n_seen = closure_stats["nb"]
            # gradient-norm convergence readout (the L-BFGS stopping signal)
            gsq = sum(float((p.grad.detach()**2).sum()) for g in optim.param_groups
                      for p in g["params"] if p.grad is not None)
            lr_str = f"|g|={gsq**0.5:.2e}"
        else:
            # total=n_batches_total → tqdm renders a % bar (None on epoch 1).
            bar = tqdm(train_loader, total=n_batches_total,
                       desc=f"[{stage_name}] epoch {epoch:>3}/{epochs}",
                       leave=False, disable=not args.progress, unit="batch")
            # Σw-weighted gradient accumulator (float64) → the full-dataset
            # gradient of the MEAN NLL, ‖g‖ = ‖Σ_b sw_b·∇meanNLL_b‖ / Σw. Like
            # the reported train_nll it averages over the epoch's moving θ, so it
            # is a per-epoch convergence readout, not the gradient at a fixed θ
            # (it → the true ‖g‖ as θ stops moving near the optimum).
            fit_params = [p for grp in optim.param_groups for p in grp["params"]]
            g_acc = {}
            for batch in bar:
                if prof is not None:
                    prof.mark_iter_start()
                batch = _move_batch(batch, device)
                if prof is not None:
                    prof.after_move()
                optim.zero_grad(set_to_none=True)
                loss, sw = step_fn(model, batch)
                n_seen += 1
                if sw <= 0:
                    continue
                if not torch.isfinite(loss):
                    if args.nan_on_step == "raise":
                        raise RuntimeError(f"non-finite loss [{stage_name}] epoch {epoch}")
                    if args.nan_on_step == "skip":
                        bar.set_postfix_str("SKIP NaN"); continue
                if prof is not None:
                    prof.after_forward()
                loss.backward()
                for p in fit_params:                  # accumulate BEFORE the step
                    if p.grad is not None:
                        if p not in g_acc:
                            g_acc[p] = torch.zeros(p.numel(), dtype=torch.float64,
                                                   device=p.device)
                        g_acc[p].add_(p.grad.detach().reshape(-1).double(), alpha=sw)
                optim.step()
                if prof is not None:
                    prof.after_backward()
                tr_sum += float(loss.item()) * sw; tr_w += sw
                bar.set_postfix_str(f"nll={tr_sum / max(tr_w, 1e-30):+.4f} lr={lr_str}")
            bar.close()
            if g_acc:
                gsq = sum(float(v.pow(2).sum()) for v in g_acc.values())
                epoch_gnorm = (gsq ** 0.5) / max(tr_w, 1e-30)
            if prof is not None:
                prof.report(stage_name, epoch)
        if n_batches_total is None:
            n_batches_total = n_seen   # exact count for the % bar from epoch 2 on
        train_nll = tr_sum / max(tr_w, 1e-30)

        # Monitored metric: held-out val NLL (stage 1) or the training NLL itself
        # (stage 2 — all events, no held-out split; skip the val pass).
        if monitor == "val" and val_loader is not None:
            model.eval(); v_sum = 0.0; v_w = 0.0
            for batch in val_loader:
                batch = _move_batch(batch, device)
                with torch.enable_grad():
                    loss, sw = step_fn(model, batch)
                if sw <= 0:
                    continue
                v_sum += float(loss.item()) * sw; v_w += sw
            val_nll = v_sum / max(v_w, 1e-30); metric = val_nll
        else:
            val_nll = train_nll; v_w = tr_w; metric = train_nll

        d_train = None if prev_train_nll is None else (train_nll - prev_train_nll)
        # Fixed notation for visible changes; scientific once |Δ| is too small to
        # show at 4 decimals (e.g. near the train-NLL plateau), so it never reads
        # as a flat "+0.0000".
        d_str = ("n/a" if d_train is None
                 else (f"{d_train:+.4f}" if abs(d_train) >= 1e-3
                       else f"{d_train:+.2e}"))
        extra = ""
        if stage_name == "fit":
            if model.theta_mode == "mlp":
                wn = max((p.detach().abs().max().item()
                          for p in model.theta_net.parameters()), default=0.0)
                extra = f" θ_net‖w‖∞={wn:.3e}"
            else:
                extra = (f" θ_scale‖∞={model.theta_scale.abs().max().item():.3e}"
                         f" θ_smear‖∞={model.theta_smear.abs().max().item():.3e}")
        val_str = f"val_nll={val_nll:+.4f} " if monitor == "val" else ""
        g_str = f"|g|={epoch_gnorm:.2e} " if epoch_gnorm is not None else ""
        print(f"[{stage_name}] epoch {epoch:>3}: train_nll={train_nll:+.4f} "
              f"(Δ={d_str}) {val_str}(Σw={v_w:.2e}) lr={lr_str} {g_str}"
              f"dt={time.time()-t0:.1f}s{extra}")
        prev_train_nll = train_nll

        improved = metric < best_val - args.patience_threshold
        if improved:
            best_val = metric; no_improve = 0
        else:
            no_improve += 1
        torch.save(_ckpt(epoch, best_val, metric), last_ckpt)
        if improved:
            torch.save(_ckpt(epoch, best_val, metric), best_ckpt)
        # Step the LR schedule. The early-stop counter is NOT reset on lr
        # reduction — `patience` is the total number of stalled epochs allowed
        # since the last improvement, regardless of any lr drops that happen
        # within that span. The plateau schedule still ticks (lr_reduce_patience
        # → factor) but those reductions just run in parallel with the
        # early-stop accumulator; they don't grant additional credit. Tested at
        # `lr = fit_{scale,smear}_lr = 0.1` (well-conditioned, O(1) θ): the
        # protective effect of the old reset is no longer needed, and the new
        # semantics save the ~5–10 min/fit otherwise spent traversing the full
        # lr schedule on already-converged runs. Re-enable manually by raising
        # `--patience` or with `--no-early-stop` for diagnostic runs.
        if sched is not None:
            sched.step(metric) if sched_kind == "plateau" else sched.step()
        if not improved and not args.no_early_stop and no_improve >= args.patience:
            print(f"[{stage_name}] early-stop: no improvement for {no_improve} epochs")
            break

    if os.path.exists(best_ckpt):
        ckpt = torch.load(best_ckpt, map_location=device, weights_only=False)
        model.load_state_dict(ckpt["state_dict"])
        print(f"[{stage_name}] reloaded best ({best_ckpt}, "
              f"{monitor}_nll={ckpt['best_val']:+.4f})")
    return best_val


def _run_trust_region(args, model, params, train_loader, stats, step_fn=None, *,
                      stage_name="fit", mc_as_data=False):
    # ``step_fn`` is accepted for call-site parity with _run_epochs but UNUSED:
    # the trust driver builds its own per-event objective (``_per_event_term``)
    # so the Σw-weighted reduction + backward run in FLOAT64 (the model stays at
    # --precision). It reproduces step2's data-branch NLL (mc_as_data-aware).
    """SciPy trust-region (trust-krylov / trust-ncg) minimisation of the stage-2
    NLL — a genuinely second-order, line-search-free minimiser for the smooth,
    deterministic objective, robust to the indefinite/degenerate curvature
    (A/e, a/c) that stalls the L-BFGS line search.

    Cost design (see the discussion): the GRADIENT is exact on the FULL sample
    (one pass/step, detached — it sets the step and is the convergence test);
    the exact (true-Hessian) HVP is computed WITHOUT forming H, on a fixed event
    SUBSET, by one of:

      • scheme A (``--trust-hvp reuse``, default): build the differentiable
        gradient g_S once per step with create_graph=True and RETAIN its graph;
        each Krylov HVP is then ``grad(g_S, θ, grad_outputs=v, retain_graph=True)``
        — a single second backward, NO data pass / forward per HVP. Cheapest per
        HVP, but the whole |S|-event forward+1st-backward graph stays resident for
        the inner solve (memory ∝ |S|, not chunkable).
      • scheme B (``--trust-hvp recompute``): each HVP re-does forward+backward
        over S (double-backward), chunked over events so peak memory is bounded
        by the chunk — the only scheme that could scale to the full sample.

    The full-gradient + full-objective trust-region ratio test makes ANY Hessian
    subset safe: a noisy H_S only yields a sub-optimal step (ρ small → shrink Δ),
    never a wrong-direction one. The Hessian subset is FIXED across all HVPs of a
    step (resampled per step) — mandatory for the Krylov inner solve's operator
    consistency. Returns the best (lowest) full-sample NLL reached."""
    try:
        from scipy.optimize import minimize as _scipy_min
    except ImportError as e:
        raise RuntimeError(f"--fit-optimizer trust-* requires scipy ({e})")
    device = args.device
    method = getattr(args, "trust_method", "trust-krylov")
    hvp_mode = getattr(args, "trust_hvp", "reuse")
    # 0 (default) → no sub-batch chunking: a double-backward HVP is only ~2× a
    # plain gradient in peak memory (it holds the forward activations + the
    # retained create_graph backward graph), so if a full-batch GRADIENT fits —
    # which it must, the training step runs full batches — the full-batch HVP
    # fits too. Sub-chunking the batch is pure Python/kernel overhead for the
    # same total compute; only set >0 if genuinely memory-bound.
    hvp_chunk = int(getattr(args, "trust_hvp_chunk", 0) or 0)
    hsub = int(getattr(args, "hess_subsample_events", 0) or 0)   # 0 = full sample
    sub_str = (f"{hsub:,}" if hsub > 0 else "full")
    max_iter = int(args.fit_epochs or args.epochs)
    # Flat parameter vector ↔ the active θ (+ bkg) params.
    shapes = [p.shape for p in params]
    numels = [int(p.numel()) for p in params]
    n_par = int(sum(numels))
    verbose = bool(getattr(args, "progress", True))
    # Evaluation counters (objective/grad passes, g_S rebuilds, Krylov HVPs) — the
    # trust-region inner solve is otherwise silent, so a slow step looks like a
    # hang; report them so the per-iteration cost is visible.
    ctr = {"fg": 0, "gS": 0, "hvp": 0, "hvp_since": 0, "iter": 0}

    def _pbar(desc):
        return tqdm(total=None, desc=desc, leave=False,
                    disable=not verbose, unit="batch")

    def _set_flat(x_np):
        x = torch.as_tensor(x_np, dtype=torch.float32, device=device)
        off = 0
        with torch.no_grad():
            for p, n, sh in zip(params, numels, shapes):
                p.copy_(x[off:off + n].view(sh)); off += n
        return x

    def _flat_grad(grads):
        return torch.cat([
            (g if g is not None else torch.zeros_like(p)).reshape(-1)
            for g, p in zip(grads, params)])

    # Per-event NLL + weights for ONE batch, reduced in FLOAT64. The model
    # (flow + gh_qop) runs at --precision; we upcast the per-event log-density to
    # double BEFORE the Σw-weighted reduction + backward, so the in-batch sum and
    # the gradient accumulation along the backward path are float64 (the dominant
    # remaining floor — a 65k-event float32 sum carries ~√N·ε ≈ 3e-5 relative
    # cancellation). This does NOT recover the model's internal float32 round-off
    # (per-event log p is only float32-accurate), only stops adding to it.
    def _per_event_term(batch):
        dm = (~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"])
        per = model.data_nll_continuity(
            batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
            batch["q_pm"], batch["b_pm"], batch["cond_std"], dm,
            n_iter=int(getattr(args, "continuity_n_iter", 2)))
        w = (batch["w"] * dm.to(batch["w"].dtype)).double()      # float64 weights
        term = (w * per.double()).sum()                          # Σw·NLL in float64
        return term, float(w.sum())

    # ---- full-sample exact objective + gradient (one pass/step, detached) ----
    def _fun_and_grad(x_np):
        _set_flat(x_np)
        model.zero_grad(set_to_none=True)
        ctr["fg"] += 1
        s = 0.0; w = 0.0; nev = 0
        # Accumulate in float64: the per-batch grads are float32 (float32 model),
        # but summing O(N/batch) of them in float32 incurs ~√N·ε ≈ 3e-5 relative
        # cancellation noise — right at the ‖g‖ floor that stalled trust-krylov.
        # float64 accumulation drops that floor to ~1e-12 so a small gtol is
        # meaningful; the Σw-weighted loss sum is likewise summed in float64.
        gacc = torch.zeros(n_par, device=device, dtype=torch.float64)
        bar = _pbar(f"[{stage_name}] grad pass #{ctr['fg']} (full sample)")
        for batch in train_loader:
            batch = _move_batch(batch, device)
            term, sw = _per_event_term(batch)        # Σw·NLL (float64), Σw
            nb = int(batch["mll"].shape[0]); nev += nb; bar.update(1)
            if sw <= 0 or not torch.isfinite(term):
                continue
            g = torch.autograd.grad(term, params, allow_unused=True)
            gacc += _flat_grad(g).detach().double()
            s += float(term.item()); w += sw         # term is already Σw·NLL
            bar.set_postfix_str(f"events={nev:,} nll={s / max(w, 1e-30):+.5f}")
        bar.close()
        if w <= 0:
            raise RuntimeError(f"trust-region: no usable batches [{stage_name}]")
        inv_w = 1.0 / w
        _fun_and_grad.last_nll = s * inv_w
        gnorm = float((gacc * inv_w).norm())
        if verbose:
            print(f"  [{stage_name}] f/grad eval #{ctr['fg']}: nll={s*inv_w:+.6f} "
                  f"‖g‖={gnorm:.3e}  (HVPs since last iter: {ctr['hvp_since']})",
                  flush=True)
        ctr["hvp_since"] = 0
        return float(s * inv_w), (gacc * inv_w).double().cpu().numpy()
    _fun_and_grad.last_nll = float("inf")

    # ---- Hessian subset loader: a fixed, deterministic subset of S events ----
    def _hess_batches():
        seen = 0
        for batch in train_loader:
            yield batch
            seen += int(batch["mll"].shape[0])
            if hsub > 0 and seen >= hsub:
                return

    # Scheme A: per-step retained differentiable gradient over the subset.
    state = {"x_key": None, "ver": None, "gS": None, "wS": 1.0}

    def _param_versions():
        # In-place op counter per param leaf; bumps whenever _set_flat's copy_
        # mutates a param (e.g. scipy evaluating a trial step fun(x+p)).
        return tuple(p._version for p in params)

    def _build_reuse_grad():
        # Build g_S = ∂(Σ_S w·nll)/∂θ with the graph retained (create_graph).
        # NOTE: the whole |S|-event graph stays alive for all HVPs this step.
        model.zero_grad(set_to_none=True)
        ctr["gS"] += 1
        accs = None; w = 0.0; nev = 0
        bar = _pbar(f"[{stage_name}] HVP graph build #{ctr['gS']} "
                    f"({sub_str} events)")
        for batch in _hess_batches():
            batch = _move_batch(batch, device)
            term, sw = _per_event_term(batch)        # Σw·NLL (float64), Σw
            nev += int(batch["mll"].shape[0]); bar.update(1)
            if sw <= 0 or not torch.isfinite(term):
                continue
            accs = term if accs is None else accs + term
            w += sw
            bar.set_postfix_str(f"events={nev:,}")
        bar.close()
        if accs is None or w <= 0:
            raise RuntimeError("trust-region HVP(reuse): no usable subset events")
        gS = torch.autograd.grad(accs, params, create_graph=True)
        state["gS"] = _flat_grad(gS)        # differentiable, graph retained
        state["wS"] = w
        if verbose:
            print(f"  [{stage_name}] built retained HVP graph #{ctr['gS']} "
                  f"({nev:,} events) — subsequent HVPs reuse it (no data pass)",
                  flush=True)

    def _hvp_reuse(x_np, v_np):
        # Rebuild g_S when the retained graph is stale: either the iterate x
        # changed, OR the param leaves were modified in place since the build
        # (scipy mutates them via _set_flat when evaluating a trial step — on a
        # REJECTED step x is unchanged but the params were left perturbed, so an
        # x-only check would reuse a stale graph → "modified by an inplace op"
        # version error). Within one subproblem solve scipy makes no fun calls,
        # so versions are stable and the retained graph is reused across HVPs.
        if (state["gS"] is None or state["x_key"] != _arr_key(x_np)
                or state["ver"] != _param_versions()):
            _set_flat(x_np)
            _build_reuse_grad()
            state["x_key"] = _arr_key(x_np)
            state["ver"] = _param_versions()
        v = torch.as_tensor(v_np, dtype=torch.float32, device=device)
        hv = torch.autograd.grad(state["gS"], params, grad_outputs=v,
                                 retain_graph=True, allow_unused=True)
        hv = _flat_grad(hv) / state["wS"]
        ctr["hvp"] += 1; ctr["hvp_since"] += 1
        if verbose and ctr["hvp_since"] % 10 == 0:
            print(f"  [{stage_name}] Krylov HVP {ctr['hvp_since']} this "
                  f"subproblem ({ctr['hvp']} total)", flush=True)
        return hv.detach().double().cpu().numpy()

    # Scheme B: recompute forward+double-backward per HVP, chunked over events.
    def _hvp_recompute(x_np, v_np):
        _set_flat(x_np)
        ctr["hvp"] += 1; ctr["hvp_since"] += 1
        v = torch.as_tensor(v_np, dtype=torch.float32, device=device)
        hv = torch.zeros(n_par, device=device, dtype=torch.float64); w = 0.0
        # NB tqdm ticks once per CHUNK (hvp_chunk events), not per loader batch —
        # hvp_chunk<=0 → whole batch (no sub-batch chunking); tqdm then ticks once
        # per loader batch (a double-backward HVP fits at ~2× a full-batch
        # gradient). Only >0 sub-chunks the batch (more, smaller passes) for the
        # rare memory-bound case.
        chunk_str = (f"{hvp_chunk}/chunk" if hvp_chunk > 0 else "whole-batch")
        bar = _pbar(f"[{stage_name}] HVP #{ctr['hvp']} (recompute, "
                    f"{sub_str} events, {chunk_str})")
        for batch in _hess_batches():
            batch = _move_batch(batch, device)
            n = int(batch["mll"].shape[0])
            step = n if hvp_chunk <= 0 else hvp_chunk
            for c0 in range(0, n, step):
                sub = {k: (val[c0:c0 + step] if torch.is_tensor(val)
                           and val.shape[:1] == (n,) else val)
                       for k, val in batch.items()}
                term, sw = _per_event_term(sub)      # Σw·NLL (float64 reduction)
                bar.update(1)
                if sw <= 0 or not torch.isfinite(term):
                    continue
                g = torch.autograd.grad(term, params, create_graph=True)
                gflat = _flat_grad(g)                 # float32 (fp32 param leaves)
                hvc = torch.autograd.grad(gflat, params, grad_outputs=v,
                                          retain_graph=False, allow_unused=True)
                hv += _flat_grad(hvc).detach().double(); w += sw
        bar.close()
        if w <= 0:
            raise RuntimeError("trust-region HVP(recompute): no usable subset events")
        if verbose and ctr["hvp_since"] % 10 == 0:
            print(f"  [{stage_name}] Krylov HVP {ctr['hvp_since']} this "
                  f"subproblem ({ctr['hvp']} total)", flush=True)
        return (hv / w).double().cpu().numpy()

    # trust-exact wants the FULL Hessian matrix (n_par×n_par), accumulated over
    # the subset and assembled with the vectorised batched second-backward (same
    # _hessian_block_batched used by the observed output-fisher), with a per-row
    # loop fallback. The reduction is float64 (per-event term upcast). Feasible
    # only when n_par is small (binned θ); a guard warns/blocks for large n_par.
    _h_active_idx = torch.arange(n_par, device=device)
    _h_use_batched = bool(getattr(args, "fisher_vectorized", True))

    def _full_hessian(x_np):
        nonlocal _h_use_batched
        _set_flat(x_np)
        H = torch.zeros((n_par, n_par), device=device, dtype=torch.float64)
        w = 0.0
        bar = _pbar(f"[{stage_name}] Hessian build (exact, {sub_str} events)")
        for batch in _hess_batches():
            batch = _move_batch(batch, device)
            term, sw = _per_event_term(batch)        # Σw·NLL (float64 reduction)
            bar.update(1)
            if sw <= 0 or not torch.isfinite(term):
                continue
            g = torch.autograd.grad(term, params, create_graph=True)
            g_full = torch.cat([gi.reshape(-1) for gi in g])
            if _h_use_batched:
                try:
                    Hb = _hessian_block_batched(
                        g_full[_h_active_idx], params, _h_active_idx, n_par,
                        chunk=max(1, int(args.empirical_fisher_chunk)))
                except (RuntimeError, NotImplementedError) as e:
                    _h_use_batched = False
                    if str(device).startswith("cuda"):
                        torch.cuda.empty_cache()
                    bar.write(f"  note: vectorised Hessian unavailable "
                              f"({type(e).__name__}); using the per-row loop")
                    Hb = _hessian_block_loop(g_full, params, _h_active_idx, n_par)
            else:
                Hb = _hessian_block_loop(g_full, params, _h_active_idx, n_par)
            H += Hb.detach().double(); w += sw
        bar.close()
        if w <= 0:
            raise RuntimeError("trust-region: no usable subset events (Hessian)")
        ctr["hvp"] += 1                              # count as one curvature eval
        H = H / w                                    # per-unit-weight (mean) Hessian
        H = 0.5 * (H + H.T)                           # symmetrise ULP asymmetry
        if verbose:
            print(f"  [{stage_name}] built exact Hessian ({n_par}×{n_par}) "
                  f"this iter", flush=True)
        return H.cpu().numpy()

    use_exact_hess = (method == "trust-exact")
    if use_exact_hess and n_par > int(getattr(args, "trust_exact_max_par", 400)):
        raise RuntimeError(
            f"--fit-optimizer trust-exact builds the full {n_par}×{n_par} Hessian "
            f"per iteration, which is too large (>{getattr(args,'trust_exact_max_par',400)} "
            f"params; e.g. --theta-mlp). Use trust-krylov (Hessian-free HVP) or "
            f"raise --trust-exact-max-par if you really intend this.")
    hessp = _hvp_reuse if hvp_mode == "reuse" else _hvp_recompute
    x0 = torch.cat([p.detach().reshape(-1) for p in params]).double().cpu().numpy()
    if use_exact_hess:
        print(f"  optimizer: {method} (2nd-order trust region) — full-sample exact "
              f"gradient; FULL {n_par}×{n_par} Hessian (vectorised second-backward, "
              f"float64) on {sub_str} events; max {max_iter} iters, "
              f"gtol={args.trust_gtol:g}")
    else:
        print(f"  optimizer: {method} (2nd-order trust region) — full-sample exact "
              f"gradient; exact HVP via scheme {'A/reuse' if hvp_mode=='reuse' else 'B/recompute'} "
              f"on {sub_str} events; max {max_iter} iters, gtol={args.trust_gtol:g}")

    best = {"nll": float("inf"), "x": x0.copy()}

    def _callback(xk, *a):
        ctr["iter"] += 1
        nll = _fun_and_grad.last_nll
        if nll < best["nll"]:
            best["nll"] = nll; best["x"] = np.array(xk, copy=True)
        print(f"  [{stage_name}] trust iter {ctr['iter']}/{max_iter}: "
              f"nll={nll:+.6f}  (cum: {ctr['fg']} f/grad evals, "
              f"{ctr['gS']} graph builds, {ctr['hvp']} HVPs)", flush=True)

    if use_exact_hess:
        res = _scipy_min(
            _fun_and_grad, x0, method=method, jac=True, hess=_full_hessian,
            callback=_callback,
            options={"maxiter": max_iter, "gtol": float(args.trust_gtol)})
    else:
        res = _scipy_min(
            _fun_and_grad, x0, method=method, jac=True, hessp=hessp,
            callback=_callback,
            options={"maxiter": max_iter, "gtol": float(args.trust_gtol)})
    # Use the best iterate seen (scipy returns the last, which the ratio test
    # guarantees is ≤ start, but the callback-tracked best is safest).
    x_final = best["x"] if best["nll"] <= float(res.fun) else res.x
    _set_flat(x_final)
    final_nll = min(best["nll"], float(res.fun))
    print(f"  [{stage_name}] trust-region done: nll={final_nll:+.6f}, "
          f"‖g‖={np.linalg.norm(res.jac):.3e}, {res.nit} iters, "
          f"success={res.success} ({res.message})")
    ck = {
        "epoch": res.nit, "stage": stage_name,
        "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()},
        "theta_scale": model.theta_scale.detach().cpu(),
        "theta_smear": model.theta_smear.detach().cpu(),
        "stats": _stats_to_dict(stats), "best_val": final_nll,
        "val_metric": final_nll, "args": vars(args),
    }
    torch.save(ck, os.path.join(args.output, "fit_best.pt"))
    torch.save(ck, os.path.join(args.output, "fit_last.pt"))
    return final_nll


def _arr_key(x_np):
    """Cheap content key for a small parameter vector (to detect when scipy
    hands the same iterate across the HVPs of one trust-region step)."""
    return (float(x_np[0]), float(x_np[-1]), float(np.asarray(x_np).sum()),
            int(x_np.size))


def train_stage1(args, model, train_loader, val_loader, stats) -> float:
    """Stage 1: fit the nominal flow p₀(m|muon_kin) on simulation (MC rows)."""
    print("\n=== stage 1: nominal flow on simulation (θ=0, no θ conditioning) ===")
    optim = torch.optim.Adam(model.flow.parameters(), lr=args.lr,
                             weight_decay=args.weight_decay)
    print(f"  optimizer: flow ({sum(p.numel() for p in model.flow.parameters()):,} params), lr={args.lr:g}")

    def step1(model, batch):
        idx = (~batch["is_data_mask"]).nonzero(as_tuple=True)[0]
        if idx.numel() == 0:
            return torch.zeros((), dtype=torch.float64,
                               device=batch["mll"].device), 0.0
        logp = model.log_p_nominal(batch["mll"][idx], batch["cond_std"][idx])
        w = batch["w"][idx].double()                 # float64 reduction: the model
        sw = float(w.sum().clamp_min(1e-30))         # runs at --precision, but the
        return -(w * logp.double()).sum() / sw, sw   # Σw·NLL sum + its backward are
        # float64 (a ~65k-event float32 sum carries ~√N·ε cancellation; backward
        # accumulates per-event grad contributions in float64, cast to the fp32
        # leaf only at the end). Benefits adam/soap/lbfgs alike (all call step_fn).

    return _run_epochs(args, model, optim, train_loader, val_loader, stats,
                       step_fn=step1, ckpt_prefix="flow", stage_name="flow",
                       epochs=args.flow_epochs or args.epochs)


def train_stage2(args, model, train_loader, val_loader, stats,
                 *, mc_as_data: bool = False) -> float:
    """Stage 2: freeze the flow, fit θ + background MLP on data (continuity).

    With ``mc_as_data`` (MC-closure validation mode) the simulation rows are
    treated as the pseudo-data branch (``~is_data_mask``) so θ is fit against a
    disjoint half of simulation; the closure target is θ → 0.
    """
    src = "MC pseudo-data" if mc_as_data else "data"
    print(f"\n=== stage 2: θ + background fit on {src} (frozen flow, #2 direct-eval) ===")
    # Freeze the flow; it is the nominal template from stage 1.
    for p in model.flow.parameters():
        p.requires_grad_(False)
    # If the background mixture is disabled the data branch is pure signal —
    # the MLP is bypassed in data_nll_continuity and stays at its random init
    # with no gradient flow. Freeze its parameters and drop the group from
    # the optimiser so nothing depends on its (now-irrelevant) values.
    groups = []
    tags = []
    if model.background_enabled:
        groups.append({"params": model.mlp.parameters(), "lr": args.fit_mlp_lr})
        tags.append("mlp")
    else:
        for p in model.mlp.parameters():
            p.requires_grad_(False)
        print("  background DISABLED: f_data ≡ [0, 0, 1] (pure signal); MLP frozen + excluded from optimiser")
    if model.theta_mode == "mlp":
        # Continuous θ(η,φ): float the ThetaNet (zero-init → 0). One group;
        # the output reference scaling differentiates the A,e,M vs a,c magnitudes.
        groups.append({"params": model.theta_net.parameters(),
                       "lr": args.fit_theta_mlp_lr}); tags.append("θ_net")
        print(f"  optimizer groups: {', '.join(tags)}  "
              f"(lr mlp={args.fit_mlp_lr:g} θ_net={args.fit_theta_mlp_lr:g}); "
              f"θ = ThetaNet(η, φ) [continuous]")
    else:
        # θ_scale is the advective shift (init 0, signed). θ_smear are the
        # qop-variance coefficients. Their DEFAULT init comes from the model
        # (0 for linear/softplus; SMEAR_SQUARE_INIT_RAW for `square`, since raw=0
        # is a dead saddle there — ∂effective/∂raw=0 → the smear stays frozen at
        # 0). For `softplus`, `θ=0` is physical c ≈ softplus(0)·SCALE_C ≈ 0.69·
        # SCALE_C, near the saturation knee — use --init-theta-{a,c} to start at
        # a positive raw θ. --init-theta-{a,c} (when set) OVERRIDES the model
        # default here for the binned table.
        a0 = float(getattr(args, "init_theta_a", 0.0))
        c0 = float(getattr(args, "init_theta_c", 0.0))
        with torch.no_grad():
            model.theta_scale.zero_()
            if a0 != 0.0 or c0 != 0.0:
                model.theta_smear[:, 0].fill_(a0)
                model.theta_smear[:, 1].fill_(c0)
                print(f"  θ_smear init (raw, all η-bins): a={a0:g}  c={c0:g} "
                      f"(per-bin mask still applied in the forward pass; the "
                      f"frozen column is inert regardless of init)")
        if not args.disable_scale:
            groups.append({"params": [model.theta_scale], "lr": args.fit_scale_lr}); tags.append("θ_scale")
        if not args.disable_smearing:
            groups.append({"params": [model.theta_smear], "lr": args.fit_smear_lr}); tags.append("θ_smear")
        print(f"  optimizer groups: {', '.join(tags)}  "
              f"(lr mlp={args.fit_mlp_lr:g} scale={args.fit_scale_lr:g} smear={args.fit_smear_lr:g})")
    print(f"  signal density: #2 direct-eval (advection + probability-flow smear, "
          f"flow_steps={getattr(args, 'smear_flow_steps', 1)}, n_iter="
          f"{args.continuity_n_iter}); normalised by construction")

    def step2(model, batch):
        # In validation mode the simulation rows play the role of data.
        data_mask = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
        per = model.data_nll_continuity(
            batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
            batch["q_pm"], batch["b_pm"], batch["cond_std"], data_mask,
            n_iter=args.continuity_n_iter)
        # float64 Σw·NLL reduction + backward (model stays at --precision); same
        # rationale as step1 / the trust driver's _per_event_term, applied to all
        # optimisers that go through step_fn (adam/soap/lbfgs).
        w = (batch["w"] * data_mask.to(batch["w"].dtype)).double()
        sw = float(w.sum().clamp_min(1e-30))
        return (w * per.double()).sum() / sw, sw

    fit_opt = getattr(args, "fit_optimizer", "adam")
    max_epochs = args.fit_epochs or args.epochs
    if fit_opt in ("trust-krylov", "trust-ncg", "trust-exact"):
        # Second-order trust region (scipy) over the active θ (+ bkg) params.
        tr_params = [p for g in groups for p in g["params"]]
        args.trust_method = fit_opt
        return _run_trust_region(args, model, tr_params, train_loader, stats,
                                 step2, stage_name="fit", mc_as_data=mc_as_data)
    if "+" in fit_opt:
        # Two-phase hybrid base+polish: run the base optimiser (adam/soap) to its
        # normal stopping for robust bulk descent, then warm-start the polish
        # (lbfgs / trust-krylov / trust-ncg — both reach a tighter gradient norm
        # than Adam/SOAP can) from the reloaded best. Groups depend only on the
        # (unchanged) model state, so they're rebuilt fresh per phase.
        base, polish = fit_opt.split("+", 1)
        is_trust = polish in ("trust-krylov", "trust-ncg", "trust-exact")
        n2 = (max_epochs if is_trust          # trust uses --fit-epochs as its iter cap
              else int(getattr(args, "lbfgs_final_epochs", 1)))
        print(f"  optimizer: HYBRID {fit_opt} — phase 1 = {base} (≤{max_epochs} "
              f"epochs), phase 2 = {polish} (warm-started)")
        print(f"  --- phase 1: {base} ---")
        optim = _make_fit_optimizer(args, groups, kind=base)
        best1 = _run_epochs(args, model, optim, train_loader, val_loader, stats,
                            step_fn=step2, ckpt_prefix="fit", stage_name="fit",
                            epochs=max_epochs, monitor="train")
        # Snapshot phase-1's best so a stalled/regressing polish can't leave a
        # WORSE checkpoint (phase 2 may overwrite fit_best.pt regardless).
        fit_best = os.path.join(args.output, "fit_best.pt")
        ph1_ckpt = os.path.join(args.output, "fit_phase1_best.pt")
        if os.path.exists(fit_best):
            shutil.copyfile(fit_best, ph1_ckpt)
        print(f"  --- phase 2: {polish} polish (warm-started from {base}, "
              f"phase-1 best train_nll={best1:+.5f}) ---")
        if is_trust:
            tr_params = [p for g in groups for p in g["params"]]
            args.trust_method = polish
            best2 = _run_trust_region(args, model, tr_params, train_loader, stats,
                                      step2, stage_name="fit", mc_as_data=mc_as_data)
        else:
            optim = _make_fit_optimizer(args, groups, kind="lbfgs")
            best2 = _run_epochs(args, model, optim, train_loader, val_loader, stats,
                                step_fn=step2, ckpt_prefix="fit", stage_name="fit",
                                epochs=n2, monitor="train")
        if best1 < best2 and os.path.exists(ph1_ckpt):
            # Polish regressed — restore phase-1's checkpoint + weights.
            print(f"  {polish} polish did not improve (phase1={best1:+.5f} < "
                  f"phase2={best2:+.5f}); keeping phase-1 result.")
            shutil.copyfile(ph1_ckpt, fit_best)
            ck = torch.load(fit_best, map_location=args.device, weights_only=False)
            model.load_state_dict(ck["state_dict"])
            best2 = best1
        if os.path.exists(ph1_ckpt):
            os.remove(ph1_ckpt)
        return best2
    optim = _make_fit_optimizer(args, groups)
    print(f"  optimizer: {fit_opt}")
    return _run_epochs(args, model, optim, train_loader, val_loader, stats,
                       step_fn=step2, ckpt_prefix="fit", stage_name="fit",
                       epochs=max_epochs, monitor="train")


def train_loop(args: argparse.Namespace) -> int:
    """The two-stage continuity pipeline (stage 1 flow → stage 2 fit)."""
    device = args.device
    if device.startswith("cuda") and not torch.cuda.is_available():
        print("warning: CUDA requested but not available; falling back to CPU")
        device = args.device = "cpu"

    # --stage uncertainties: load an existing FULL fit (flow + MLP + θ) from
    # --checkpoint and run only the uncertainty estimation (Fisher / bootstrap),
    # no training. Handled in its own path.
    if args.stage == "uncertainties":
        return _run_uncertainties_stage(args, device)

    # --stage fit: load the flow checkpoint up front so its architecture and
    # stats drive the model build + standardisation (the model is rebuilt from
    # args before the flow weights are loaded, so they must match).
    flow_ck = None
    flow_ckpt = None
    stats_override = None
    if args.stage == "fit":
        flow_ckpt = args.checkpoint or os.path.join(args.output, "flow_best.pt")
        if not os.path.exists(flow_ckpt):
            print(f"error: --stage fit needs a stage-1 flow; {flow_ckpt!r} not found "
                  f"(run --stage flow first or pass --checkpoint)", file=sys.stderr)
            return 1
        print(f"loading stage-1 flow checkpoint: {flow_ckpt}")
        flow_ck = torch.load(flow_ckpt, map_location=device, weights_only=False)
        _apply_flow_arch_from_ckpt(args, flow_ck.get("args", {}) or {})
        # Reuse the flow's own standardisation unless the user forces --stats-in.
        if args.stats_in is None and flow_ck.get("stats") is not None:
            stats_override = _stats_from_dict(flow_ck["stats"])

    setup = _setup_common(args, stats_override=stats_override)
    if setup is None:
        return 1
    shard_files, stats, train_loader, val_loader = setup

    model = _build_model(args, stats, device)
    print(f"model: flow={model.flow_arch} (no θ conditioning), "
          f"scale={'on' if model.scale_enabled else 'off'} "
          f"smear={'on' if model.smearing_enabled else 'off'} "
          f"smear_fit={model.smear_fit_params}")

    if args.stage == "fit":
        # Load ONLY the flow weights; the background MLP and θ start fresh
        # (stage 1 leaves them at init anyway) so the MLP size / smear-fit
        # choice remain free stage-2 knobs.
        flow_sd = {k[len("flow."):]: v for k, v in flow_ck["state_dict"].items()
                   if k.startswith("flow.")}
        model.flow.load_state_dict(flow_sd)
        print(f"loaded flow weights from {flow_ckpt} ({len(flow_sd)} tensors; "
              f"mlp + θ start fresh)")

    # MC-closure validation: both stages run on simulation, with a deterministic
    # disjoint half each (stage 1 ← half 0, stage 2 ← half 1 treated as
    # pseudo-data). The disjoint halves keep stage 2 from fitting θ against the
    # very events stage 1's flow was trained on; the closure target is θ → 0.
    if args.validation:
        inj = _inject_theta_np(args, len(stats.eta_edges) - 1)
        inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1)
        tgt = "θ → 0"
        if inj is not None or inj_sm is not None:
            parts = []
            if inj is not None:
                parts.append(f"(A,e,M)=({args.inject_A:g},{args.inject_e:g},{args.inject_M:g})")
            if inj_sm is not None:
                parts.append(f"(a,c)=({args.inject_a:g},{args.inject_c:g})")
            tgt = "θ → injected " + " ".join(parts)
        h_flow, h_fit = _validation_half(args, "flow"), _validation_half(args, "fit")
        split_desc = ("stage 1 ← ALL events, stage 2 ← ALL events as pseudo-data "
                      "(no half split per --no-validation-split)" if h_flow is None
                      else "stage 1 ← half 0, stage 2 ← half 1 as pseudo-data")
        print(f"\n*** MC-closure validation mode: simulation for both stages "
              f"({split_desc}); target {tgt} ***")
        if inj is not None:
            print("    injecting the θ_scale shift into the stage-2 pseudo-data m_ll")
        if inj_sm is not None:
            print("    injecting the per-muon qop smear into the stage-2 pseudo-data m_ll")
        s1_train, s1_val = _make_loaders(args, shard_files, stats, half=h_flow)   # flow: NOT injected
        # Fit: ALL events of its half (no held-out val/holdout); stops on train NLL.
        s2_train, s2_val = _make_loaders(args, shard_files, stats, half=h_fit,
                                         inject_theta=inj, inject_smear=inj_sm,
                                         val_fraction=0.0, holdout_fraction=0.0)
    else:
        if (_inject_theta_np(args, len(stats.eta_edges) - 1) is not None
                or _inject_smear_np(args, len(stats.eta_edges) - 1) is not None):
            print("warning: --inject-A/e/M/a/c only apply in --validation mode; ignoring.",
                  file=sys.stderr)
        s1_train, s1_val = train_loader, val_loader
        # Fit: ALL events (no held-out val/holdout); stops on train NLL.
        s2_train, s2_val = _make_loaders(args, shard_files, stats,
                                         val_fraction=0.0, holdout_fraction=0.0)

    if args.stage in ("both", "flow"):
        train_stage1(args, model, s1_train, s1_val, stats)
    if args.stage in ("both", "fit"):
        train_stage2(args, model, s2_train, s2_val, stats, mc_as_data=args.validation)
        if args.fisher_info:
            _run_fisher_continuity(args, model, shard_files, stats, device)
        if args.empirical_fisher or getattr(args, "output_fisher", False):
            _run_empirical_fisher(args, model, shard_files, stats, device)
        if args.bootstrap > 0:
            run_bootstrap_continuity(args, model, shard_files, stats, device,
                                     mc_as_data=args.validation)
    return 0


def _run_fisher_continuity(args, model, shard_files, stats, device) -> None:
    """Observed Fisher info for the two-stage fit (θ_scale + active θ_smear,
    fixed flow + MLP) → ``<output>/fisher_info.pt``. Evaluated on the data the
    fit used (``--fisher-split``, default the train split; half 1 in
    validation mode), so the covariance has the right statistical scale."""
    if not (model.scale_enabled or model.smearing_enabled):
        print("skipping Fisher info: both --disable-scale and --disable-smearing.")
        return
    if model.theta_mode == "mlp":
        print("skipping observed Fisher info: not implemented for --theta-mlp "
              "(binned 24×3 θ layout); use --empirical-fisher over the net weights "
              "or re-run binned for the per-bin covariance.", file=sys.stderr)
        return
    half = _validation_half(args, "fit")
    inj = _inject_theta_np(args, len(stats.eta_edges) - 1) if args.validation else None
    inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1) if args.validation else None
    loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=args.batch_size, split=args.fisher_split,
        val_fraction=0.0, holdout_fraction=0.0,   # all events, matching the all-events fit
        drop_last=False, half=half, inject_theta_scale=inj,
        inject_theta_smear=inj_sm, inject_seed=int(args.inject_smear_seed),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=int(getattr(args, "max_events", 0) or 0),
        event_fraction=float(getattr(args, "event_fraction", 1.0) or 1.0))
    print(f"\ncomputing observed Fisher information (θ_scale + active θ_smear, "
          f"fixed flow + MLP) on split={args.fisher_split}"
          + ("  half=%s (MC pseudo-data)" % ('all' if half is None else half)
             if args.validation else "")
          + f"  [smear_fit={model.smear_fit_params}]")
    t0 = time.time()
    H, layout = compute_fisher_info_continuity(
        model, loader, device, mc_as_data=args.validation,
        n_iter=args.continuity_n_iter,
        progress=args.progress, vectorized=args.fisher_vectorized)
    out = _fisher_save_dict(H, layout, model)
    path = os.path.join(args.output, "fisher_info.pt")
    torch.save(out, path)
    print(f"  wrote {path}: {H.shape[0]}×{H.shape[0]} info matrix over "
          f"{len(out['labels'])} params ({layout['seen']:,} events, "
          f"Σw={layout['sw']:.2e}); "
          f"{'PD, inverted' if out['ok'] else 'NOT positive-definite → pinv (cov approximate)'} "
          f"in {time.time()-t0:.1f}s")
    if out.get("edm") is not None:
        print(f"  EDM (½ gᵀV g, est. NLL distance to minimum) = {out['edm']:.3e}"
              f"  (‖grad‖={out['grad_norm']:.3e})")
    if out["n_negative_eig"] != 0:
        print(f"  WARNING: observed information has {out['n_negative_eig']} "
              f"non-positive eigenvalue(s) (min={out['min_eig']:.3e}) — the fit is "
              f"not at a clean optimum or some η-bins are event-starved; the "
              f"corresponding variances are unreliable.")


def _bootstrap_save_dict(cov, mean, TH, eff_smears, conv_epochs, model, smear_cols):
    """Package the warm-start bootstrap covariance, diagnostics-compatible with
    fisher_info.pt (same θ_scale-block / σ keys; no Hessian/EDM keys)."""
    n_scale = model.theta_scale.numel() if model.scale_enabled else 0
    out = {
        "method": "warm-start Poisson bootstrap",
        "covariance": cov,
        "mean": mean,
        "replicas": TH,                       # [B, n_act] raw active vectors
        "labels": _active_param_labels(model, smear_cols),
        "n_scale": n_scale,
        "smear_cols": smear_cols,
        "smear_fit_params": model.smear_fit_params,
        "n_replicas": int(TH.shape[0]),
        "convergence_epochs": conv_epochs,
        "param_space": "scale: linear (A,e,M); smear: raw pre-softplus theta_smear",
    }
    n_eta_s = model.theta_scale.shape[0]
    if n_scale == 3 * n_eta_s:
        # θ_scale O(1) → physical (A,e,M) = θ·THETA_SCALE_REF; cov_phys = cov·ref⊗ref.
        refN = torch.tensor(list(THETA_SCALE_REF) * n_eta_s, dtype=cov.dtype)
        cov_scale = cov[:n_scale, :n_scale] * refN.unsqueeze(0) * refN.unsqueeze(1)
        out["covariance_24_3_24_3"] = cov_scale.view(n_eta_s, 3, n_eta_s, 3)
        out["sigma_scale_24_3"] = torch.sqrt(
            torch.clamp(torch.diag(cov_scale), min=0.0)).view(n_eta_s, 3)
    # Physical smear σ straight from the replicas' effective_theta_smear (already
    # physical) — the bootstrap gives the effective-space spread exactly.
    if smear_cols and eff_smears:
        EFF = torch.stack(eff_smears)                       # [B, n_smear_active]
        sig = EFF.std(0, unbiased=True) if EFF.shape[0] > 1 else torch.zeros(EFF.shape[1])
        n_eta, n_comp = model.theta_smear.shape
        sig_eff = torch.zeros(n_eta, n_comp)
        k = 0
        for b in range(n_eta):
            for c in smear_cols:
                sig_eff[b, c] = sig[k]
                k += 1
        out["sigma_smear_eff_24_2"] = sig_eff
    return out


def run_bootstrap_continuity(args, model, shard_files, stats, device, *,
                             mc_as_data: bool) -> None:
    """Warm-start Poisson bootstrap of the stage-2 fit → ``<output>/bootstrap_cov.pt``.

    Each replica restarts from the nominal (θ̂, φ̂), resets Adam to the initial
    fit LRs (re-raised, so the replica can actually relax), and refits θ AND the
    background MLP jointly on the SAME data Poisson(1)-reweighted (an independent
    per-event count, fixed across the replica's epochs), until its reweighted NLL
    plateaus. The covariance of {θ̂_b} folds in the background (and every other)
    uncertainty with no Hessian. Warm-starting only changes how each replica
    reaches its optimum, not where — provided it re-converges (hence the per-
    replica early stop on the reweighted NLL)."""
    B = int(args.bootstrap)
    if B <= 0:
        return
    if model.theta_mode == "mlp":
        print("skipping bootstrap: not implemented for --theta-mlp (binned 24×3 "
              "θ layout / per-bin covariance).", file=sys.stderr)
        return
    if not (model.scale_enabled or model.smearing_enabled):
        print("skipping bootstrap: both --disable-scale and --disable-smearing.")
        return
    smear_cols = _smear_active_cols(model) if model.smearing_enabled else []
    # mc_as_data is always == args.validation at the call site, so the helper
    # gives the right half (or None when --no-validation-split is set).
    half = _validation_half(args, "fit") if mc_as_data else None
    inj = _inject_theta_np(args, len(stats.eta_edges) - 1) if mc_as_data else None
    inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1) if mc_as_data else None
    loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=args.batch_size, split=args.fisher_split,
        val_fraction=0.0, holdout_fraction=0.0,   # all events, matching the all-events fit
        drop_last=False, half=half, inject_theta_scale=inj,
        inject_theta_smear=inj_sm, inject_seed=int(args.inject_smear_seed),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=int(getattr(args, "max_events", 0) or 0),
        event_fraction=float(getattr(args, "event_fraction", 1.0) or 1.0))
    nominal_sd = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    # Freeze the flow; float the MLP + active θ (the MLP must re-fit per replica
    # so its uncertainty enters the spread).
    for p in model.flow.parameters():
        p.requires_grad_(False)
    for p in model.mlp.parameters():
        p.requires_grad_(True)
    model.theta_scale.requires_grad_(model.scale_enabled)
    model.theta_smear.requires_grad_(model.smearing_enabled)

    # Convergence settings synchronised with the nominal stage-2 fit: the
    # per-replica early stop (patience + threshold), epoch cap, and LR schedule
    # all default to the nominal fit's, so each replica converges as thoroughly
    # as the nominal (monitored on the replica's reweighted NLL, since a
    # bootstrap replica has no separate validation set).
    patience = (args.bootstrap_patience if args.bootstrap_patience is not None
                else args.patience)
    max_epochs = (args.bootstrap_epochs if args.bootstrap_epochs is not None
                  else (args.fit_epochs or args.epochs))
    print(f"\nwarm-start Poisson bootstrap: {B} replicas × ≤{max_epochs} epochs "
          f"on split={args.fisher_split}"
          + ("  half=%s (MC pseudo-data)" % ('all' if half is None else half)
             if mc_as_data else "")
          + f"  [smear_fit={model.smear_fit_params}; "
          + ("no early stop" if args.no_early_stop
             else f"patience={patience}, threshold={args.patience_threshold:g}")
          + f", lr-schedule={args.lr_schedule}]")
    gen = torch.Generator(device=device)
    replicas, eff_smears, conv_epochs = [], [], []
    n_batches_total = None   # learned on the first epoch → % on the inner bar after
    t0 = time.time()
    rbar = tqdm(range(B), desc="bootstrap", disable=not args.progress, unit="replica")
    for b in rbar:
        model.load_state_dict(nominal_sd)          # warm-start from the nominal fit
        groups = []
        if model.background_enabled:
            groups.append({"params": model.mlp.parameters(), "lr": args.fit_mlp_lr})
        if model.scale_enabled:
            groups.append({"params": [model.theta_scale], "lr": args.fit_scale_lr})
        if model.smearing_enabled:
            groups.append({"params": [model.theta_smear], "lr": args.fit_smear_lr})
        optim = _make_fit_optimizer(args, groups, minibatch_loop=True)  # fresh state (Adam/SOAP; LBFGS→Adam here)
        sched, sched_kind = _make_scheduler(args, optim, max_epochs)  # same as nominal
        seed = args.bootstrap_seed + b
        best = float("inf"); no_improve = 0; used = 0
        model.train()
        for epoch in range(max_epochs):
            gen.manual_seed(seed)                  # same per-event Poisson each epoch
            tr_sum = 0.0; tr_w = 0.0; n_seen = 0
            # Inner per-epoch bar over batches so a (slow) epoch visibly advances;
            # total is learned on epoch 1 so later epochs render a % bar.
            ebar = tqdm(loader, total=n_batches_total, leave=False,
                        desc=f"  replica {b + 1}/{B} epoch {epoch + 1}/{args.bootstrap_epochs}",
                        disable=not args.progress, unit="batch")
            for batch in ebar:
                n_seen += 1
                batch = _move_batch(batch, device)
                data_mask = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
                if not bool(data_mask.any()):
                    continue
                pois = torch.poisson(
                    torch.ones(batch["mll"].shape[0], device=device), generator=gen)
                w = batch["w"] * pois * data_mask.to(batch["w"].dtype)
                sw = float(w.sum().clamp_min(1e-30))
                if sw <= 0:
                    continue
                per = model.data_nll_continuity(
                    batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
                    batch["q_pm"], batch["b_pm"], batch["cond_std"], data_mask,
                    n_iter=args.continuity_n_iter)
                loss = (w * per).sum() / sw
                if not torch.isfinite(loss):
                    continue
                optim.zero_grad(set_to_none=True)
                loss.backward()
                optim.step()
                tr_sum += float(loss.item()) * sw; tr_w += sw
                ebar.set_postfix_str(f"nll={tr_sum / max(tr_w, 1e-30):+.4f}")
            ebar.close()
            if n_batches_total is None:
                n_batches_total = n_seen
            used = epoch + 1
            nll = tr_sum / max(tr_w, 1e-30)
            improved = nll < best - args.patience_threshold
            if improved:
                best = nll; no_improve = 0
            else:
                no_improve += 1
            # Same LR schedule + early-stop coupling as the nominal fit: a LR
            # reduction resets the early-stop counter (so it only fires once the
            # reductions are exhausted), monitored on the replica's reweighted NLL.
            if sched is not None:
                lr_before = [g["lr"] for g in optim.param_groups]
                sched.step(nll) if sched_kind == "plateau" else sched.step()
                if any(g["lr"] < lb - 1e-12 for g, lb in zip(optim.param_groups, lr_before)):
                    no_improve = 0
            if not improved and not args.no_early_stop and no_improve >= patience:
                break
        replicas.append(_record_active_theta(model, smear_cols))
        if smear_cols:
            eff_smears.append(
                model.effective_theta_smear().detach()[:, smear_cols].reshape(-1).cpu())
        conv_epochs.append(used)
        rbar.set_postfix_str(f"ep={used} nll={best:+.4f}")
    rbar.close()
    model.load_state_dict(nominal_sd)              # leave the model at the nominal fit

    TH = torch.stack(replicas)                     # [B, n_act]
    mean = TH.mean(0)
    Xc = TH - mean
    cov = (Xc.t() @ Xc) / max(B - 1, 1)
    out = _bootstrap_save_dict(cov, mean, TH, eff_smears, conv_epochs, model, smear_cols)
    path = os.path.join(args.output, "bootstrap_cov.pt")
    torch.save(out, path)
    ce = torch.tensor(conv_epochs, dtype=torch.float32)
    n_hit_cap = int((ce >= max_epochs).sum()) if not args.no_early_stop else 0
    print(f"  wrote {path}: {B} replicas over {TH.shape[1]} params; "
          f"epochs/replica median={int(ce.median())} max={int(ce.max())} "
          f"in {time.time()-t0:.1f}s")
    if n_hit_cap:
        print(f"  WARNING: {n_hit_cap}/{B} replica(s) hit the {max_epochs}-epoch "
              f"cap without plateauing — they may be under-converged (variance "
              f"underestimated). Raise --bootstrap-epochs.")
    if "sigma_scale_24_3" in out:
        ss = out["sigma_scale_24_3"]
        print(f"  bootstrap σ(A,e,M) median over bins = "
              f"({float(ss[:,0].median()):.2e}, {float(ss[:,1].median()):.2e}, "
              f"{float(ss[:,2].median()):.2e})")


def _theta_cov_extras(cov_theta: torch.Tensor, model, smear_cols, n_scale: int) -> dict:
    """Diagnostics-compatible extras from a θ-block covariance (cpu): the
    θ_scale block in the 24×3×24×3 layout, its √diag σ, and the delta-method
    effective σ for the raw θ_smear. Shared by the Fisher / empirical builders."""
    out: dict = {}
    n_eta_s = model.theta_scale.shape[0]
    if n_scale == 3 * n_eta_s:
        # θ_scale O(1) → physical (A,e,M) = θ·THETA_SCALE_REF; cov_phys = cov·ref⊗ref.
        refN = torch.tensor(list(THETA_SCALE_REF) * n_eta_s, dtype=cov_theta.dtype)
        cs = cov_theta[:n_scale, :n_scale] * refN.unsqueeze(0) * refN.unsqueeze(1)
        out["covariance_24_3_24_3"] = cs.view(n_eta_s, 3, n_eta_s, 3)
        out["sigma_scale_24_3"] = torch.sqrt(
            torch.clamp(torch.diag(cs), min=0.0)).view(n_eta_s, 3)
    if smear_cols:
        n_eta, n_comp = model.theta_smear.shape
        cv = cov_theta[n_scale:, n_scale:]
        sig_raw = torch.sqrt(torch.clamp(torch.diag(cv), min=0.0))
        smear_scale = (SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C)  # θ→physical (linear)
        sig_eff = torch.zeros(n_eta, n_comp)
        k = 0
        for b in range(n_eta):
            for c in smear_cols:
                sig_eff[b, c] = smear_scale[c] * sig_raw[k]
                k += 1
        out["sigma_smear_eff_24_2"] = sig_eff
        # Full PHYSICAL (a,c) covariance (24×2×24×2) — needed for the whitened
        # (stiff/sloppy) band, which mixes a and c. Only when BOTH float (the
        # whitened smear plot requires smear_fit_params=='both'); cv is then
        # ordered (bin0_a, bin0_c, bin1_a, …).
        if smear_cols == [0, 1]:
            sv = torch.tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C] * n_eta,
                              dtype=cv.dtype)
            cv_phys = cv * sv.unsqueeze(0) * sv.unsqueeze(1)
            out["covariance_smear_24_2_24_2"] = cv_phys.view(n_eta, 2, n_eta, 2)
    return out


def compute_empirical_fisher_joint(
    model: JpsiMassMixtureModel,
    loader: JpsiMassArrowLoader,
    device: str,
    *,
    mc_as_data: bool = False,
    n_iter: int = 2,
    chunk_events: int = 64,
    max_events: int = 0,
    progress: bool = True,
    vectorized: bool = True,
):
    """Joint (θ, φ) empirical Fisher ``J = Σ_i w_i s_i s_iᵀ`` over θ_scale + the
    active θ_smear + ALL background-MLP params, from per-event scores
    ``s_i = ∇_(θ,φ) log p_{θ,φ}(x_i)`` (the data branch). PSD by construction
    (sum of outer products), so its pseudo-inverse never yields a negative
    variance; the θ-block of ``pinv(J)`` is the nuisance-marginalised θ
    covariance — background included — and uses only FIRST derivatives (no MLP
    Hessian). Per-event scores are taken in chunks via ``is_grads_batched`` with
    a per-event-loop fallback. Returns ``(J [n_act, n_act] cpu, layout)`` with
    the active vector ordered ``[θ-active | mlp]``."""
    model.eval()
    for p in model.flow.parameters():
        p.requires_grad_(False)
    for p in model.mlp.parameters():
        p.requires_grad_(True)
    model.theta_scale.requires_grad_(model.scale_enabled)
    model.theta_smear.requires_grad_(model.smearing_enabled)

    params, names, numels = [], [], []
    if model.scale_enabled:
        params.append(model.theta_scale); names.append("scale")
        numels.append(model.theta_scale.numel())
    smear_cols = _smear_active_cols(model) if model.smearing_enabled else []
    if model.smearing_enabled:
        params.append(model.theta_smear); names.append("smear")
        numels.append(model.theta_smear.numel())
    mlp_params = list(model.mlp.parameters())
    n_mlp = int(sum(p.numel() for p in mlp_params))
    for p in mlp_params:
        params.append(p); names.append("mlp"); numels.append(p.numel())
    if not params or (not model.scale_enabled and not smear_cols):
        raise RuntimeError("empirical Fisher: no active θ parameters.")

    # Active flat indices into the concatenated score: θ_scale (all) + active
    # θ_smear cols + all MLP — θ first, then φ.
    active_idx, off, n_theta_active = [], 0, 0
    for nm, ne in zip(names, numels):
        if nm == "scale":
            active_idx += list(range(off, off + ne)); n_theta_active += ne
        elif nm == "smear":
            n_eta, n_comp = model.theta_smear.shape
            a = [off + b * n_comp + c for b in range(n_eta) for c in smear_cols]
            active_idx += a; n_theta_active += len(a)
        else:
            active_idx += list(range(off, off + ne))
        off += ne
    active_idx = torch.tensor(active_idx, dtype=torch.long, device=device)
    n_act = int(active_idx.numel())

    J = torch.zeros((n_act, n_act), device=device, dtype=torch.float32)
    sw = 0.0; seen = 0; hit_cap = False
    use_batched = bool(vectorized)
    bar = tqdm(loader, desc="emp-fisher", disable=not progress, unit="batch")
    for batch in bar:
        if max_events > 0 and seen >= max_events:
            hit_cap = True; break
        batch = _move_batch(batch, device)
        data_mask = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
        di = data_mask.nonzero(as_tuple=True)[0]
        if di.numel() == 0:
            continue
        per = model.data_nll_continuity(
            batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
            batch["q_pm"], batch["b_pm"], batch["cond_std"], data_mask,
            n_iter=n_iter)
        per_d = per[di]
        w_d = batch["w"][di].detach()
        nd_total = int(per_d.shape[0])
        # Truncate to the remaining max_events budget if it falls mid-batch:
        # only the first `nd` data events of this batch get score-computed.
        if max_events > 0 and seen + nd_total > max_events:
            nd = max(0, max_events - seen)
        else:
            nd = nd_total
        if nd == 0:
            hit_cap = True; break
        for s0 in range(0, nd, max(1, chunk_events)):
            s1 = min(s0 + chunk_events, nd)
            c = s1 - s0
            S = None
            if use_batched:
                try:
                    eye = torch.zeros((c, nd), device=device, dtype=per_d.dtype)
                    eye[torch.arange(c, device=device),
                        torch.arange(s0, s1, device=device)] = 1.0
                    g = torch.autograd.grad(
                        per_d, params, grad_outputs=eye, is_grads_batched=True,
                        retain_graph=True, allow_unused=True)
                    S = torch.cat([
                        (gi if gi is not None else torch.zeros((c,) + p.shape,
                         device=device, dtype=per_d.dtype)).reshape(c, -1)
                        for gi, p in zip(g, params)], dim=1)
                except (RuntimeError, NotImplementedError) as e:
                    use_batched = False
                    bar.write(f"  note: vectorised per-event score unavailable "
                              f"({type(e).__name__}); using the per-event loop")
            if S is None:
                rows = []
                for j in range(s0, s1):
                    gj = torch.autograd.grad(per_d[j], params, retain_graph=True,
                                             allow_unused=True)
                    rows.append(torch.cat([
                        (x if x is not None else torch.zeros_like(p)).reshape(-1)
                        for x, p in zip(gj, params)]))
                S = torch.stack(rows)
            S = S[:, active_idx]                       # [c, n_act]
            wc = w_d[s0:s1]
            J += (S * wc.unsqueeze(1)).t() @ S         # Σ w_i s_i s_iᵀ
        # Account only the events we actually scored (after the mid-batch cap).
        seen += nd; sw += float(w_d[:nd].sum())
        bar.set_postfix_str(f"events={seen:,}")
    bar.close()
    if seen == 0:
        raise RuntimeError("empirical Fisher: zero data-branch events seen.")
    J = 0.5 * (J + J.T)                                # symmetrise (numerical)
    layout = {"n_theta_active": n_theta_active,
              "n_scale": (model.theta_scale.numel() if model.scale_enabled else 0),
              "smear_cols": smear_cols, "n_mlp": n_mlp,
              "sw": sw, "seen": seen, "hit_cap": hit_cap}
    return J.detach().cpu(), layout


def _empirical_cov_theta_block(J: torch.Tensor, n_theta: int, ridge: float) -> torch.Tensor:
    """θ-block covariance from the joint empirical Fisher ``J`` (PSD).

    ``ridge == 0`` → Moore–Penrose ``pinv(J)``: unconstrained / degenerate
    directions land in the null space and get **zero** variance (the misleading
    σ≈0 on a near-degenerate parameter, e.g. the A/e degeneracy over the J/ψ
    pt range).

    ``ridge > 0`` → scale-aware ridge ``cov = (J + ridge·diag(J))⁻¹``, computed
    as a plain inverse in the per-parameter standardised space
    ``D⁻¹(D⁻¹JD⁻¹ + ridge·I)⁻¹D⁻¹`` with ``D = √diag(J)``. A flat direction then
    reads as a LARGE but finite variance (~1/(ridge·J_ii)) rather than 0, and
    the regularisation is proportional to each parameter's own information
    (dimensionless ``ridge``), so it is well-behaved across the mixed-unit
    A/e/M/smear/MLP blocks."""
    J = 0.5 * (J + J.T)
    # Parameters whose per-event score is bit-zero (frozen via the
    # smear/scale param mask, or just an MLP weight that doesn't feed the
    # data branch) have J_ii = 0 exactly. They're not fit, so their variance
    # is 0 by definition — but the inversion can leak large garbage into
    # those rows/cols: pinv may amplify near-zero singular values via the
    # joint-cov coupling, and the scale-aware ridge inflates them by 1/d² as
    # d → its 1e-12·dmax floor, producing ~10²³ artefacts that pollute the
    # χ² eigenvalue tolerance and the correlation plot. Capture the frozen
    # mask up front to zero out those rows/cols at the end.
    zero_diag_mask = (torch.diag(J) == 0)
    if ridge <= 0.0:
        cov = torch.linalg.pinv(J)
    else:
        d = torch.sqrt(torch.clamp(torch.diag(J), min=0.0))
        dmax = float(d.max()) if d.numel() else 1.0
        dinv = 1.0 / torch.clamp(d, min=1e-12 * (dmax if dmax > 0 else 1.0))
        Jt = J * dinv.unsqueeze(0) * dinv.unsqueeze(1)            # D⁻¹ J D⁻¹ (PSD)
        n = Jt.shape[0]
        eye = torch.eye(n, dtype=J.dtype, device=J.device)
        cov = torch.linalg.inv(Jt + ridge * eye)                 # PD → plain inv
        cov = cov * dinv.unsqueeze(0) * dinv.unsqueeze(1)        # back to raw units
    if zero_diag_mask.any():
        cov[zero_diag_mask, :] = 0.0
        cov[:, zero_diag_mask] = 0.0
    return cov[:n_theta, :n_theta].contiguous()


def compute_empirical_fisher_net(
    model: JpsiMassMixtureModel,
    loader: JpsiMassArrowLoader,
    device: str,
    *,
    mc_as_data: bool = False,
    n_iter: int = 2,
    chunk_events: int = 64,
    max_events: int = 0,
    progress: bool = True,
    vectorized: bool = True,
):
    """Empirical Fisher ``J = Σ_i w_i s_i s_iᵀ`` over the θ-NET weights (+ the
    background MLP, for marginalisation) from per-event data-branch scores, for
    ``--theta-mlp``. There the per-η (A,e,M)/(a,c) come from ``theta_net``, so
    the statistical uncertainty lives in its weights; the theta_net block of
    ``(J + ridge)⁻¹`` is the weight covariance, which the caller propagates to
    the per-η outputs via the net's output Jacobian. θ_net params are ordered
    FIRST so that block is ``cov[:n_net, :n_net]``. Returns ``(J cpu, layout)``.
    """
    model.eval()
    for p in model.flow.parameters():
        p.requires_grad_(False)
    # theta_scale / theta_smear are inert in mlp mode — freeze them.
    model.theta_scale.requires_grad_(False)
    model.theta_smear.requires_grad_(False)
    net_params = list(model.theta_net.parameters())
    for p in net_params:
        p.requires_grad_(True)
    bg_params = []
    if getattr(model, "background_enabled", True):
        bg_params = list(model.mlp.parameters())
        for p in bg_params:
            p.requires_grad_(True)
    params = net_params + bg_params
    n_net = int(sum(p.numel() for p in net_params))
    n_bg = int(sum(p.numel() for p in bg_params))
    if n_net == 0:
        raise RuntimeError("empirical Fisher (mlp): theta_net has no parameters.")

    n_act = n_net + n_bg
    J = torch.zeros((n_act, n_act), device=device, dtype=torch.float32)
    sw = 0.0; seen = 0; hit_cap = False
    use_batched = bool(vectorized)
    bar = tqdm(loader, desc="emp-fisher(net)", disable=not progress, unit="batch")
    for batch in bar:
        if max_events > 0 and seen >= max_events:
            hit_cap = True; break
        batch = _move_batch(batch, device)
        data_mask = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
        di = data_mask.nonzero(as_tuple=True)[0]
        if di.numel() == 0:
            continue
        per = model.data_nll_continuity(
            batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
            batch["q_pm"], batch["b_pm"], batch["cond_std"], data_mask,
            n_iter=n_iter)
        per_d = per[di]
        w_d = batch["w"][di].detach()
        nd_total = int(per_d.shape[0])
        nd = (max(0, max_events - seen)
              if max_events > 0 and seen + nd_total > max_events else nd_total)
        if nd == 0:
            hit_cap = True; break
        for s0 in range(0, nd, max(1, chunk_events)):
            s1 = min(s0 + chunk_events, nd)
            c = s1 - s0
            S = None
            if use_batched:
                try:
                    eye = torch.zeros((c, nd), device=device, dtype=per_d.dtype)
                    eye[torch.arange(c, device=device),
                        torch.arange(s0, s1, device=device)] = 1.0
                    g = torch.autograd.grad(
                        per_d, params, grad_outputs=eye, is_grads_batched=True,
                        retain_graph=True, allow_unused=True)
                    S = torch.cat([
                        (gi if gi is not None else torch.zeros((c,) + p.shape,
                         device=device, dtype=per_d.dtype)).reshape(c, -1)
                        for gi, p in zip(g, params)], dim=1)
                except (RuntimeError, NotImplementedError) as e:
                    use_batched = False
                    bar.write(f"  note: vectorised per-event score unavailable "
                              f"({type(e).__name__}); using the per-event loop")
            if S is None:
                rows = []
                for j in range(s0, s1):
                    gj = torch.autograd.grad(per_d[j], params, retain_graph=True,
                                             allow_unused=True)
                    rows.append(torch.cat([
                        (x if x is not None else torch.zeros_like(p)).reshape(-1)
                        for x, p in zip(gj, params)]))
                S = torch.stack(rows)
            wc = w_d[s0:s1]
            J += (S * wc.unsqueeze(1)).t() @ S
        seen += nd; sw += float(w_d[:nd].sum())
        bar.set_postfix_str(f"events={seen:,}")
    bar.close()
    if seen == 0:
        raise RuntimeError("empirical Fisher (mlp): zero data-branch events seen.")
    J = 0.5 * (J + J.T)
    layout = {"n_net": n_net, "n_bg": n_bg, "sw": sw, "seen": seen,
              "hit_cap": hit_cap}
    return J.detach().cpu(), layout


def _propagate_net_cov_to_outputs(model, cov_w, eta_edges, *, n_phi_avg=16,
                                  device="cpu"):
    """Propagate the θ-net weight covariance ``cov_w`` [n_net, n_net] to the
    per-η (A,e,M) and effective (a,c) output covariances via the Jacobian
    ``G = ∂o/∂w`` of the φ-AVERAGED outputs (the closure plot's central value),
    delta-method ``C_o = G cov_w Gᵀ``. Returns ``(cov_scale_24_3_24_3,
    sigma_scale_24_3, sigma_smear_eff_24_2, cov_smear_24_2_24_2)`` (physical
    units), matching the binned Fisher keys consumed by the diagnostics."""
    net_params = list(model.theta_net.parameters())
    cov_w = cov_w.to(torch.float64)
    centers = 0.5 * (np.asarray(eta_edges[:-1]) + np.asarray(eta_edges[1:]))
    n_eta = int(centers.shape[0])
    centers_t = torch.as_tensor(centers, dtype=torch.float32, device=device)
    # φ-average grid (uniform on the circle → exact mean), both muons share (η,φ);
    # muon-0 output is the plotted one (per-muon net, symmetry implicit).
    phi_avg = torch.linspace(0.0, 2 * np.pi * (1.0 - 1.0 / n_phi_avg), n_phi_avg,
                             dtype=torch.float32, device=device)
    eta_grid = centers_t[:, None, None].expand(n_eta, n_phi_avg, 2).reshape(-1, 2)
    phi_grid = phi_avg[None, :, None].expand(n_eta, n_phi_avg, 2).reshape(-1, 2)
    smear_scale = torch.tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C],
                               dtype=torch.float32, device=device)
    AeM_g, ac_g = model.theta_net(eta_grid, phi_grid)            # physical A,e,M; O(1) a,c
    AeM = AeM_g[:, 0, :].view(n_eta, n_phi_avg, 3).mean(dim=1)   # [n_eta,3] φ-mean
    ac_eff = model._smear_raw_to_effective(ac_g[:, 0, :])        # softplus·mask
    ac_phys = (ac_eff * smear_scale).view(n_eta, n_phi_avg, 2).mean(dim=1)  # [n_eta,2]
    o_scale = AeM.reshape(-1)                                    # [n_eta*3]
    o_smear = ac_phys.reshape(-1)                                # [n_eta*2]

    def _jac(o_flat):
        # Autograd runs on the net's device; move each row to CPU float64 so the
        # delta-method matmul below matches cov_w (returned on CPU from the
        # CPU-side Fisher) — avoids a cuda/cpu device mismatch on GPU runs.
        rows = []
        for k in range(int(o_flat.numel())):
            g = torch.autograd.grad(o_flat[k], net_params, retain_graph=True,
                                    allow_unused=True)
            rows.append(torch.cat([
                (gi if gi is not None else torch.zeros_like(p)).reshape(-1)
                for gi, p in zip(g, net_params)]).to(torch.float64).cpu())
        return torch.stack(rows)                                # [n_out, n_net] cpu

    G_s = _jac(o_scale)
    G_c = _jac(o_smear)
    cov_scale = (G_s @ cov_w @ G_s.t())                         # [n_eta*3, n_eta*3]
    cov_smear = (G_c @ cov_w @ G_c.t())                         # [n_eta*2, n_eta*2]
    sigma_scale = torch.sqrt(torch.clamp(torch.diag(cov_scale), min=0.0)
                             ).view(n_eta, 3).float()
    sigma_smear = torch.sqrt(torch.clamp(torch.diag(cov_smear), min=0.0)
                             ).view(n_eta, 2).float()
    return (cov_scale.view(n_eta, 3, n_eta, 3).float(), sigma_scale, sigma_smear,
            cov_smear.view(n_eta, 2, n_eta, 2).float())


def _run_empirical_fisher_mlp(args, model, shard_files, stats, device) -> None:
    """``--theta-mlp`` empirical Fisher: weight-space Fisher over θ_net (+ bkg
    MLP), ridge-inverted, propagated through the net's output Jacobian to the
    per-η (A,e,M)/(a,c) → ``<output>/empirical_fisher.pt`` with the diagnostics
    σ-band keys (covariance_24_3_24_3, sigma_smear_eff_24_2, theta_mode=mlp)."""
    half = _validation_half(args, "fit")
    inj = _inject_theta_np(args, len(stats.eta_edges) - 1) if args.validation else None
    inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1) if args.validation else None
    fisher_bs = args.batch_size
    if args.empirical_fisher_max_events > 0:
        fisher_bs = min(fisher_bs, max(1, args.empirical_fisher_max_events))
    loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=fisher_bs, split=args.fisher_split,
        val_fraction=0.0, holdout_fraction=0.0,   # all events, matching the all-events fit
        drop_last=False, half=half, inject_theta_scale=inj,
        inject_theta_smear=inj_sm, inject_seed=int(args.inject_smear_seed),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=int(getattr(args, "max_events", 0) or 0),
        event_fraction=float(getattr(args, "event_fraction", 1.0) or 1.0))
    ridge = float(args.empirical_fisher_ridge)
    print(f"\ncomputing θ-NET-weight empirical Fisher (per-event scores → ridge="
          f"{ridge:g} inverse → output Jacobian) on split={args.fisher_split}"
          + ("  half=%s (MC pseudo-data)" % ('all' if half is None else half)
             if args.validation else "")
          + (f"; ≤{args.empirical_fisher_max_events:,} events"
             if args.empirical_fisher_max_events > 0 else ""))
    t0 = time.time()
    J, layout = compute_empirical_fisher_net(
        model, loader, device, mc_as_data=args.validation,
        n_iter=args.continuity_n_iter, chunk_events=args.empirical_fisher_chunk,
        max_events=args.empirical_fisher_max_events,
        progress=args.progress, vectorized=args.fisher_vectorized)
    # Scale J to full statistics if the event budget capped it (J ∝ Σw).
    if layout["hit_cap"] and layout["sw"] > 0:
        sw_total = 0.0
        for batch in loader:
            dm = (~batch["is_data_mask"] if args.validation else batch["is_data_mask"])
            sw_total += float((batch["w"] * dm.to(batch["w"].dtype)).sum())
        if sw_total > layout["sw"]:
            sc = sw_total / layout["sw"]
            J = J * sc
            print(f"  scaled J by Σw_total/Σw_seen = {sc:.2f} "
                  f"(subsampled {layout['seen']:,} events)")
            layout["sw"] = sw_total
    n_net = layout["n_net"]
    rank = int(torch.linalg.matrix_rank(J).item())
    # θ_net block of the ridge-regularised inverse = the weight covariance.
    cov_w = _empirical_cov_theta_block(J, n_net, ridge)          # [n_net, n_net]
    (cov_scale_24_3, sigma_scale_24_3, sigma_smear_eff,
     cov_smear_24_2) = _propagate_net_cov_to_outputs(
        model, cov_w, stats.eta_edges, device=device)
    out = {
        "method": (f"theta-net empirical Fisher (ridge={ridge:g}) propagated to "
                   f"per-η outputs via the net Jacobian"),
        "ridge": ridge, "theta_mode": "mlp",
        "covariance_24_3_24_3": cov_scale_24_3,
        "sigma_scale_24_3": sigma_scale_24_3,
        "sigma_smear_eff_24_2": sigma_smear_eff,
        "covariance_smear_24_2_24_2": cov_smear_24_2,
        "smear_fit_params": model.smear_fit_params,
        "n_events": layout["seen"], "sum_weight": layout["sw"],
        "n_net": n_net, "n_bg": layout["n_bg"],
        "joint_dim": n_net + layout["n_bg"], "joint_rank": rank,
        "param_space": "theta_net weights → physical (A,e,M) and effective (a,c)",
    }
    path = os.path.join(args.output, "empirical_fisher.pt")
    torch.save(out, path)
    print(f"  wrote {path}: θ_net-weight {n_net}×{n_net} (+bkg {layout['n_bg']}), "
          f"joint rank {rank}/{n_net + layout['n_bg']} via ridge={ridge:g}, "
          f"propagated to per-η σ ({layout['seen']:,} events, "
          f"Σw={layout['sw']:.2e}) in {time.time()-t0:.1f}s")
    ss = sigma_scale_24_3
    print(f"  propagated σ(A,e,M) median over bins = "
          f"({float(ss[:,0].median()):.2e}, {float(ss[:,1].median()):.2e}, "
          f"{float(ss[:,2].median()):.2e})")


def compute_output_fisher_2d(model, loader, device, *, method="empirical",
                             n_phi=4, eta_edges=None, mc_as_data=False, n_iter=2,
                             chunk_events=64, max_events=0, progress=True,
                             marginalize_bkg=False, vectorized=True):
    """OUTPUT-space Fisher for --theta-mlp: reparametrise θ as a 2-D (η-bin ×
    φ-bin) TABLE seeded from the net at the bin centres, with the per-muon lookup
    reading that table, so the data's information about the per-(η,φ) OUTPUTS is
    measured directly — no 1349-weight over-parameterisation, no arbitrary
    weight-space ridge. Accumulates, over the table's ACTIVE entries (the fitted
    A/e/M and a/c columns):
      • observed  H = Σ w ∂²(−ln L)/∂θ²  (double-backward through data_nll);
      • empirical J = Σ w (∂ln L/∂θ)(…)ᵀ (per-event score outer products);
      • sandwich → both.
    The per-η φ-mean covariance is recovered downstream by averaging the table
    covariance over its φ axes. θ_scale is physical (A,e,M); θ_smear is the O(1)
    (a,c).

    ``marginalize_bkg``: ALSO float the background MLP weights and accumulate the
    JOINT [θ-table ⊕ bkg-weight] information, so the caller can Schur-complement
    the background block out → the background-marginalized θ covariance (the
    fixed-background default treats f_data(c) as known). The returned matrices are
    then ``[n_joint, n_joint]`` (θ-active block first, bkg block after), with the
    split recorded in ``layout`` (``n_theta_act``, ``n_bg``).

    Returns ``(H|None, J|None, layout)`` (cpu)."""
    need_H = method in ("observed", "sandwich")
    need_J = method in ("empirical", "sandwich")
    model.eval()
    for p in model.flow.parameters():
        p.requires_grad_(False)
    bg_params = []
    if marginalize_bkg and getattr(model, "background_enabled", True):
        bg_params = list(model.mlp.parameters())
        for p in bg_params:
            p.requires_grad_(True)
    else:
        for p in model.mlp.parameters():
            p.requires_grad_(False)
    model.theta_scale.requires_grad_(False)
    model.theta_smear.requires_grad_(False)
    centres = 0.5 * (np.asarray(eta_edges[:-1]) + np.asarray(eta_edges[1:]))
    n_eta = int(len(centres))
    eta_c = torch.as_tensor(centres, dtype=torch.float32, device=device)
    phi_edges = torch.linspace(-float(np.pi), float(np.pi), n_phi + 1, device=device)
    phi_c = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    eg = eta_c[:, None, None].expand(n_eta, n_phi, 2).reshape(-1, 2)
    pg = phi_c[None, :, None].expand(n_eta, n_phi, 2).reshape(-1, 2)
    with torch.no_grad():
        AeM, ac = model.theta_net(eg, pg)
        t2s = AeM[:, 0, :].view(n_eta, n_phi, 3).contiguous()          # physical A,e,M
        t2c = model._smear_raw_to_effective(ac[:, 0, :]).view(
            n_eta, n_phi, 2).contiguous()                             # O(1) a,c
    t2s = t2s.detach().clone().requires_grad_(model.scale_enabled)
    t2c = t2c.detach().clone().requires_grad_(model.smearing_enabled)
    bnd = phi_edges[1:-1].contiguous()

    def _phibin(phi):
        return torch.clamp(torch.bucketize(phi, bnd), 0, n_phi - 1)

    orig_s, orig_c = model._scale_AeM_pm, model._smear_ac_pm
    model._scale_AeM_pm = lambda e, p, b: t2s[b, _phibin(p)] * model.scale_param_mask
    model._smear_ac_pm = lambda e, p, b: t2c[b, _phibin(p)] * model.smear_param_mask

    scale_cols = ([c for c in range(3) if float(model.scale_param_mask[c]) != 0.0]
                  if model.scale_enabled else [])
    smear_cols = _smear_active_cols(model) if model.smearing_enabled else []
    params, active, off = [], [], 0
    if model.scale_enabled:
        params.append(t2s)
        for ie in range(n_eta):
            for ip in range(n_phi):
                for c in scale_cols:
                    active.append(off + (ie * n_phi + ip) * 3 + c)
        off += t2s.numel()
    if smear_cols:
        params.append(t2c)
        for ie in range(n_eta):
            for ip in range(n_phi):
                for c in smear_cols:
                    active.append(off + (ie * n_phi + ip) * 2 + c)
    if not active:
        for p in (t2s, t2c):
            p.requires_grad_(False)
        model._scale_AeM_pm, model._smear_ac_pm = orig_s, orig_c
        raise RuntimeError("output Fisher: no active θ columns.")
    n_theta_act = len(active)
    # Background marginalisation: append ALL bkg-MLP weights to the differentiated
    # set, with EVERY entry active (the joint block to be Schur-complemented out).
    off_after_theta = (t2s.numel() if model.scale_enabled else 0) \
        + (t2c.numel() if smear_cols else 0)
    n_bg = 0
    if bg_params:
        params = params + bg_params
        n_bg = int(sum(p.numel() for p in bg_params))
        active += [off_after_theta + i for i in range(n_bg)]
    active_idx = torch.tensor(active, dtype=torch.long, device=device)
    n_act = int(active_idx.numel())

    H = torch.zeros((n_act, n_act), device=device) if need_H else None
    J = torch.zeros((n_act, n_act), device=device) if need_J else None
    sw = 0.0; seen = 0; hit = False
    # Per-(η,φ)-cell Σw occupancy (both muons of each data event), for the
    # SAMPLE-weighted φ-collapse + η-average of the covariance downstream.
    cell_w = torch.zeros((n_eta, n_phi), device=device)
    use_batched = bool(vectorized) and need_H  # vectorise the observed-H rows
    try:
        bar = tqdm(loader, desc=f"out-fisher({method})", disable=not progress, unit="batch")
        for batch in bar:
            if max_events > 0 and seen >= max_events:
                hit = True; break
            batch = _move_batch(batch, device)
            dm = ~batch["is_data_mask"] if mc_as_data else batch["is_data_mask"]
            di = dm.nonzero(as_tuple=True)[0]
            if di.numel() == 0:
                continue
            per = model.data_nll_continuity(
                batch["mll"], batch["pt_pm"], batch["eta_pm"], batch["phi_pm"],
                batch["q_pm"], batch["b_pm"], batch["cond_std"], dm, n_iter=n_iter)
            per_d = per[di]; w_d = batch["w"][di].detach()
            nd = int(per_d.shape[0])
            if max_events > 0 and seen + nd > max_events:
                nd = max(0, max_events - seen); hit = True
            if nd == 0:
                break
            per_d = per_d[:nd]; w_d = w_d[:nd]
            # Accumulate per-cell Σw from both muons of the kept data events.
            with torch.no_grad():
                eb = batch["b_pm"][di][:nd].long()                 # [nd, 2] η-bins
                pb = _phibin(batch["phi_pm"][di][:nd])             # [nd, 2] φ-bins
                flat = (eb * n_phi + pb).reshape(-1)
                cell_w.view(-1).index_add_(
                    0, flat, w_d.unsqueeze(1).expand(-1, 2).reshape(-1))
            if need_J:                                  # per-event scores first
                for s0 in range(0, nd, max(1, chunk_events)):
                    s1 = min(s0 + chunk_events, nd)
                    rows = []
                    for j in range(s0, s1):
                        gj = torch.autograd.grad(per_d[j], params, retain_graph=True,
                                                 allow_unused=True)
                        rows.append(torch.cat([
                            (x if x is not None else torch.zeros_like(p)).reshape(-1)
                            for x, p in zip(gj, params)]))
                    S = torch.stack(rows)[:, active_idx]
                    J += (S * w_d[s0:s1].unsqueeze(1)).t() @ S
            if need_H:
                nll = (w_d * per_d).sum()
                if torch.isfinite(nll):
                    g = torch.autograd.grad(nll, params, create_graph=True,
                                            retain_graph=True)
                    g_full = torch.cat([gi.reshape(-1) for gi in g])
                    # Vectorised second backward (one vmapped vjp over the n_act
                    # identity rows) — collapses n_act serial row-passes to one,
                    # at ~n_act× backward-graph memory. Fall back to the per-row
                    # loop once on any engine/vmap failure or OOM (the inner
                    # autograd.grad of the gh_qop change-of-variables Jacobian +
                    # fixed-point clamps may be unsupported by vmap).
                    if use_batched:
                        try:
                            g_active = g_full[active_idx]
                            Hb = _hessian_block_batched(
                                g_active, params, active_idx, n_act,
                                chunk=max(1, chunk_events))
                        except (RuntimeError, NotImplementedError) as e:
                            use_batched = False
                            if str(device).startswith("cuda"):
                                torch.cuda.empty_cache()
                            bar.write(
                                f"  note: vectorised Hessian unavailable "
                                f"({type(e).__name__}: "
                                f"{str(e).splitlines()[0][:80]}); "
                                f"using the per-row loop")
                            Hb = _hessian_block_loop(g_full, params, active_idx, n_act)
                    else:
                        Hb = _hessian_block_loop(g_full, params, active_idx, n_act)
                    H += Hb.detach()
            sw += float(w_d.sum()); seen += nd
            bar.set_postfix_str(f"events={seen:,}")
        bar.close()
    finally:
        model._scale_AeM_pm, model._smear_ac_pm = orig_s, orig_c
        for p in bg_params:
            p.requires_grad_(False)
    if seen == 0:
        raise RuntimeError("output Fisher: zero data-branch events seen.")
    if H is not None:
        H = (0.5 * (H + H.T)).detach().cpu()
    if J is not None:
        J = (0.5 * (J + J.T)).detach().cpu()
    layout = {"n_eta": n_eta, "n_phi": n_phi, "scale_cols": scale_cols,
              "smear_cols": smear_cols, "n_sa": n_eta * n_phi * len(scale_cols),
              "n_ca": n_eta * n_phi * len(smear_cols), "sw": sw, "seen": seen,
              "hit_cap": hit, "n_theta_act": n_theta_act, "n_bg": n_bg,
              "cell_w": cell_w.detach().cpu().numpy()}   # [n_eta, n_phi] Σw occupancy
    return H, J, layout


def _table_output_jacobian(model, eta_edges, n_phi, scale_cols, smear_cols,
                           device):
    """Jacobian ``G = ∂o/∂w`` [n_act, n_w] of the ACTIVE η×φ table OUTPUTS w.r.t.
    the θ-net weights, in the SAME active ordering ``compute_output_fisher_2d``
    uses (scale block then smear block; each (η outer, φ middle, col inner)).

    The table outputs are exactly the net evaluated at the bin centres: scale =
    physical (A,e,M); smear = the O(1) EFFECTIVE (a,c) = ``_smear_raw_to_effective``
    — matching the units the output-space information ``J_table`` is built in, so
    ``G`` and ``J_table`` are consistent for the projection ``Uᵀ J_table U``.
    Returned on CPU float64 (to match the CPU-side information matrices)."""
    net_params = list(model.theta_net.parameters())
    centres = 0.5 * (np.asarray(eta_edges[:-1]) + np.asarray(eta_edges[1:]))
    n_eta = int(len(centres))
    eta_c = torch.as_tensor(centres, dtype=torch.float32, device=device)
    phi_edges = torch.linspace(-float(np.pi), float(np.pi), n_phi + 1, device=device)
    phi_c = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    eg = eta_c[:, None, None].expand(n_eta, n_phi, 2).reshape(-1, 2)
    pg = phi_c[None, :, None].expand(n_eta, n_phi, 2).reshape(-1, 2)
    AeM, ac = model.theta_net(eg, pg)                    # physical A,e,M; raw a,c
    AeM = AeM[:, 0, :]                                   # [n_eta*n_phi, 3]
    ac_eff = model._smear_raw_to_effective(ac[:, 0, :])  # [n_eta*n_phi, 2] effective
    # Assemble the active output vector in the exact active_idx order.
    o_rows = []
    for ie in range(n_eta):
        for ip in range(n_phi):
            row = ie * n_phi + ip
            for c in scale_cols:
                o_rows.append(AeM[row, c])
    for ie in range(n_eta):
        for ip in range(n_phi):
            row = ie * n_phi + ip
            for c in smear_cols:
                o_rows.append(ac_eff[row, c])
    o_active = torch.stack(o_rows)                       # [n_act]
    rows = []
    for k in range(int(o_active.numel())):
        g = torch.autograd.grad(o_active[k], net_params, retain_graph=True,
                                allow_unused=True)
        rows.append(torch.cat([
            (gi if gi is not None else torch.zeros_like(p)).reshape(-1)
            for gi, p in zip(g, net_params)]).to(torch.float64).cpu())
    return torch.stack(rows)                             # [n_act, n_w] cpu float64


def _marginalize_bkg_block(M, n_theta, ridge):
    """Background-marginalise a JOINT information matrix ``M`` (θ-active block
    [:n_theta] then bkg-weight block [n_theta:]) → the θ-only information via the
    Schur complement ``S = M_θθ − M_θb (M_bb)⁻¹ M_bθ``.

    ``S`` is the information about θ AFTER profiling out the background — its
    inverse is the background-marginalised θ covariance (the off-diagonal θ↔bkg
    coupling inflates θ exactly as a joint-then-invert would). The bkg block is
    rank-deficient (≫ events worth of MLP weights), so ``M_bb⁻¹ M_bθ`` is a
    RIDGE-regularised solve in bkg-weight space (``--empirical-fisher-ridge``,
    scale-aware): that ridge only damps the PSD correction term ``M_θb M_bb⁻¹
    M_bθ`` (which can only inflate θ), so over-ridging → falls back to the fixed-
    background θ block (conservative), never destabilising the primary θ inverse.
    Returns the θ-only matrix ``S`` [n_theta, n_theta]."""
    M = 0.5 * (M + M.T)
    Mtt = M[:n_theta, :n_theta]
    if M.shape[0] == n_theta:                     # no bkg block → nothing to do
        return Mtt
    Mtb = M[:n_theta, n_theta:]
    Mbb = M[n_theta:, n_theta:]
    nb = Mbb.shape[0]
    # Scale-aware ridge solve  X = (Mbb + ridge·diag(Mbb))⁻¹ Mbθ  in standardised
    # bkg space, mirroring _empirical_cov_theta_block's conditioning.
    Mbb = 0.5 * (Mbb + Mbb.T)
    d = torch.sqrt(torch.clamp(torch.diag(Mbb), min=0.0))
    dmax = float(d.max()) if d.numel() else 1.0
    dinv = 1.0 / torch.clamp(d, min=1e-12 * (dmax if dmax > 0 else 1.0))
    Mbb_s = Mbb * dinv.unsqueeze(0) * dinv.unsqueeze(1)
    eye = torch.eye(nb, dtype=M.dtype, device=M.device)
    rhs = (dinv.unsqueeze(1) * Mtb.t())           # D⁻¹ Mbθ  [nb, n_theta]
    x = torch.linalg.solve(Mbb_s + max(ridge, 1e-12) * eye, rhs)
    correction = (dinv.unsqueeze(1) * Mtb.t()).t() @ x   # Mθb D⁻¹ (…)⁻¹ D⁻¹ Mbθ
    return Mtt - 0.5 * (correction + correction.T)


def _project_to_net_subspace(M_dict, G, ridge, svd_rtol):
    """Option (b): restrict the output-space covariance to the network-REACHABLE
    subspace, so the band carries the cross-bin (smoothness) correlations the MLP
    induces while keeping the ridge in well-conditioned output space.

    ``G = U Σ Vᵀ``; ``U_r`` = left singular vectors with σ_k/σ_max > svd_rtol
    (the SMOOTHNESS CUTOFF — directions the net can only produce with large weight
    excursions are dropped). Project each information matrix ``M_a = U_rᵀ M U_r``,
    invert in that r-dim space with the scale-aware output ridge, map back
    ``C = U_r C_a U_rᵀ``. With no ridge + no truncation this equals the
    pseudoinverse limit of the weight-propagated covariance
    ``G (Gᵀ J U_r-style)⁺ Gᵀ`` (network-consistent), but computed in a clean basis.

    ``M_dict`` carries the matrices the method needs: {'J':…} (empirical),
    {'H':…} (observed), or both (sandwich). Returns ``(C [n_act,n_act], r_kept)``.
    """
    U_full, S, _ = torch.linalg.svd(G, full_matrices=False)   # U[n_act,k] S[k]
    smax = float(S[0]) if S.numel() else 0.0
    keep = (S > svd_rtol * smax) if smax > 0 else (S > 0)
    r = int(keep.sum())
    U = U_full[:, :r]                                          # [n_act, r]

    def _proj(M):
        return U.t() @ M.to(torch.float64) @ U                # [r, r]

    if "J" in M_dict and "H" in M_dict:                        # sandwich
        Ja, Ha = _proj(M_dict["J"]), _proj(M_dict["H"])
        Ha_inv = _empirical_cov_theta_block(Ha, r, ridge)
        C_a = Ha_inv @ Ja @ Ha_inv
    else:                                                      # empirical / observed
        Ma = _proj(next(iter(M_dict.values())))
        C_a = _empirical_cov_theta_block(Ma, r, ridge)
    C = (U @ C_a.to(torch.float64) @ U.t())                    # [n_act, n_act]
    return C.float(), r


def _run_output_fisher_mlp(args, model, shard_files, stats, device) -> None:
    """--theta-mlp OUTPUT-space Fisher (observed / empirical / sandwich over a
    2-D η×φ θ-table) → ``<output>/empirical_fisher.pt`` with the diagnostics
    σ-band keys, the φ-mean per-η covariance obtained by averaging the table
    covariance over φ. The smarter alternative to the weight-space empirical
    Fisher + ridge: the degeneracy structure is physical (output space), the
    regularisation scale is physical, and observed/sandwich are available.

    ``--output-fisher-project``: 'free' (default) treats every (η,φ) cell as an
    independent parameter — model-agnostic, well-conditioned, but BLIND to the
    cross-bin correlations the smooth net induces (so the φ-mean band averages
    down like 1/n_φ, often too small). 'net' restricts the covariance to the
    network-reachable output subspace (SVD of the table-output Jacobian G), so
    the band carries the MLP's smoothness correlations — the network-consistent
    covariance — while the ridge stays in clean output space."""
    method = getattr(args, "output_fisher_method", "empirical")
    project = getattr(args, "output_fisher_project", "free")
    svd_rtol = float(getattr(args, "output_fisher_svd_rtol", 1e-2))
    n_phi = max(1, int(getattr(args, "output_fisher_nphi", 4)))
    ridge = float(args.empirical_fisher_ridge)
    # Background marginalisation: tri-state CLI 'auto' (default) → ON for BOTH
    # projections (float the bkg MLP into the joint output-space info and Schur-
    # complement it out, so f_data(c) is profiled not held fixed — making the
    # band comparable to the net-weight / bootstrap bands). 'on'/'off' force it.
    mb_opt = getattr(args, "output_fisher_marginalize_bkg", "auto")
    marg_bkg = True if mb_opt == "auto" else (mb_opt == "on")
    half = _validation_half(args, "fit")
    inj = _inject_theta_np(args, len(stats.eta_edges) - 1) if args.validation else None
    inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1) if args.validation else None
    fbs = args.batch_size
    if args.empirical_fisher_max_events > 0:
        fbs = min(fbs, max(1, args.empirical_fisher_max_events))
    loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=fbs, split=args.fisher_split,
        val_fraction=0.0, holdout_fraction=0.0, drop_last=False, half=half,
        inject_theta_scale=inj, inject_theta_smear=inj_sm,
        inject_seed=int(args.inject_smear_seed),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=int(getattr(args, "max_events", 0) or 0),
        event_fraction=float(getattr(args, "event_fraction", 1.0) or 1.0))
    print(f"\ncomputing OUTPUT-space Fisher (method={method}, project={project}"
          + (f", svd_rtol={svd_rtol:g}" if project == "net" else "")
          + (", marginalize_bkg" if marg_bkg else "")
          + f", 2-D η×φ table n_phi={n_phi}, ridge={ridge:g}) "
          f"on split={args.fisher_split}"
          + ("  half=%s (MC pseudo-data)" % ('all' if half is None else half)
             if args.validation else "")
          + (f"; ≤{args.empirical_fisher_max_events:,} events"
             if args.empirical_fisher_max_events > 0 else ""))
    t0 = time.time()
    H, J, layout = compute_output_fisher_2d(
        model, loader, device, method=method, n_phi=n_phi,
        eta_edges=stats.eta_edges, mc_as_data=args.validation,
        n_iter=args.continuity_n_iter, chunk_events=args.empirical_fisher_chunk,
        max_events=args.empirical_fisher_max_events, progress=args.progress,
        marginalize_bkg=marg_bkg, vectorized=args.fisher_vectorized)
    n_act = layout["n_sa"] + layout["n_ca"]
    n_theta_act = layout.get("n_theta_act", n_act)
    # Scale H/J to full subset Σw if the event budget capped them (both ∝ Σw).
    if layout["hit_cap"] and layout["sw"] > 0:
        sw_total = 0.0
        for batch in loader:
            dm = (~batch["is_data_mask"] if args.validation else batch["is_data_mask"])
            sw_total += float((batch["w"] * dm.to(batch["w"].dtype)).sum())
        if sw_total > layout["sw"]:
            sc = sw_total / layout["sw"]
            if H is not None:
                H = H * sc
            if J is not None:
                J = J * sc
            print(f"  scaled Fisher by Σw_total/Σw_seen = {sc:.2f}")
            layout["sw"] = sw_total
    # Background marginalisation: Schur-complement the joint [θ ⊕ bkg] info down
    # to the θ-only block (commutes with the global Σw scalar above, and — for the
    # net case — with the U-projection, since the bkg block is disjoint from the
    # θ-table indices U spans). Yields the background-marginalised θ information.
    if marg_bkg and layout.get("n_bg", 0) > 0:
        if H is not None:
            H = _marginalize_bkg_block(H.to(torch.float64), n_theta_act, ridge).float()
        if J is not None:
            J = _marginalize_bkg_block(J.to(torch.float64), n_theta_act, ridge).float()
        print(f"  marginalised background ({layout['n_bg']} MLP weights) via "
              f"Schur complement → {n_theta_act}-param θ information")
    n_act = n_theta_act
    # Invert per method → active θ-table covariance C [n_act, n_act].
    r_kept = None
    if project == "net":
        # Option (b): project the output information into the network-reachable
        # subspace (smoothness correlations retained), then invert there.
        G = _table_output_jacobian(
            model, stats.eta_edges, n_phi, layout["scale_cols"],
            layout["smear_cols"], device)
        M_dict = {}
        if method in ("empirical", "sandwich"):
            M_dict["J"] = J
        if method in ("observed", "sandwich"):
            M_dict["H"] = H
        C, r_kept = _project_to_net_subspace(M_dict, G, ridge, svd_rtol)
        print(f"  projected to net subspace: kept {r_kept}/{n_act} singular "
              f"directions (svd_rtol={svd_rtol:g})")
    elif method == "empirical":
        C = _empirical_cov_theta_block(J, n_act, ridge)
    elif method == "observed":
        C = _empirical_cov_theta_block(H, n_act, ridge)      # robust ridge inverse
    else:  # sandwich  C = H⁻¹ J H⁻¹
        Hinv = _empirical_cov_theta_block(H, n_act, ridge)
        C = Hinv @ J @ Hinv
    C = C.numpy()
    n_sa, n_eta, nphi = layout["n_sa"], layout["n_eta"], layout["n_phi"]
    sc_cols, cc_cols = layout["scale_cols"], layout["smear_cols"]
    proj_tag = (f", project=net[{r_kept}/{n_act}], svd_rtol={svd_rtol:g}"
                if project == "net" else ", project=free")
    bkg_tag = (f", bkg-marginalised[{layout['n_bg']}w]"
               if (marg_bkg and layout.get("n_bg", 0) > 0) else "")
    out = {
        "method": (f"output-space Fisher ({method}, 2-D η×φ n_phi={nphi}, "
                   f"ridge={ridge:g}{proj_tag}{bkg_tag})"),
        "ridge": ridge, "theta_mode": "mlp", "output_fisher_method": method,
        "output_fisher_project": project,
        "output_fisher_svd_rtol": (svd_rtol if project == "net" else None),
        "net_subspace_rank": (r_kept if project == "net" else None),
        "background_marginalized": bool(marg_bkg and layout.get("n_bg", 0) > 0),
        "n_bkg_weights_marginalized": (layout.get("n_bg", 0) if marg_bkg else 0),
        "smear_fit_params": model.smear_fit_params,
        "n_events": layout["seen"], "sum_weight": layout["sw"],
        "param_space": "2-D (η,φ) θ-table outputs (physical A,e,M; O(1) a,c)",
        "n_phi": nphi,
    }
    # SAMPLE-weighted φ-collapse: the per-η φ-mean is the OCCUPANCY-weighted mean
    # θ̄_η = Σ_φ a_{η,φ} θ_{η,φ}, a_{η,φ} = w_{η,φ}/Σ_φ w_{η,φ} (cell Σw from the
    # data, not uniform 1/n_φ). Cov of that weighted mean over φ:
    #   cov_φmean[η,η'] = Σ_{φ,φ'} a_{η,φ} a_{η',φ'} C[(η,φ),(η',φ')].
    # Also save the FULL 2-D (η,φ) covariance + the φ-weights so the diagnostics
    # can draw a φ-RESOLVED band; fall back to uniform weights for empty cells.
    cell_w = np.asarray(layout.get("cell_w"))                # [n_eta, n_phi] Σw
    if cell_w is None or cell_w.shape != (n_eta, nphi):
        cell_w = np.ones((n_eta, nphi), dtype=np.float64)
    aw = cell_w.astype(np.float64).copy()
    rs = aw.sum(axis=1, keepdims=True)
    unif = (rs <= 0).reshape(-1)
    aw[unif] = 1.0; rs[unif.reshape(-1, 1)] = nphi          # empty η → uniform φ
    aw = aw / rs                                            # [n_eta, n_phi] weights

    def _phi_collapse(block, cols):
        # block: C sub-block [n_eta*nphi*ncol]² → weighted φ-mean [n_eta,ncol,n_eta,ncol]
        nc = len(cols)
        B = block.reshape(n_eta, nphi, nc, n_eta, nphi, nc)
        # Σ_{φ,φ'} a[η,φ] a[η',φ'] B[η,φ,·,η',φ',·]
        return np.einsum('ep,EP,epcEPC->ecEC', aw, aw, B)

    out["phi_weights"] = torch.tensor(aw, dtype=torch.float32)
    out["cell_w"] = torch.tensor(cell_w, dtype=torch.float32)   # [n_eta,n_phi] Σw
    if sc_cols:
        Cs = _phi_collapse(C[:n_sa, :n_sa], sc_cols)
        cov_s = np.zeros((n_eta, 3, n_eta, 3), dtype=np.float64)
        for i, ci in enumerate(sc_cols):
            for j, cj in enumerate(sc_cols):
                cov_s[:, ci, :, cj] = Cs[:, i, :, j]
        cs2 = cov_s.reshape(n_eta * 3, n_eta * 3)         # symmetrise ULP asymmetry
        cov_s = (0.5 * (cs2 + cs2.T)).reshape(n_eta, 3, n_eta, 3)
        out["covariance_24_3_24_3"] = torch.tensor(cov_s, dtype=torch.float32).view(n_eta, 3, n_eta, 3)
        out["sigma_scale_24_3"] = torch.sqrt(torch.clamp(
            torch.tensor(np.einsum('icic->ic', cov_s)), min=0.0)).float()
        # FULL 2-D (η,φ) covariance [n_eta,nphi,ncol-in-3, …] embedded in 3-col.
        Cs2d = C[:n_sa, :n_sa].reshape(n_eta, nphi, len(sc_cols),
                                       n_eta, nphi, len(sc_cols))
        full_s = np.zeros((n_eta, nphi, 3, n_eta, nphi, 3), dtype=np.float64)
        for i, ci in enumerate(sc_cols):
            for j, cj in enumerate(sc_cols):
                full_s[:, :, ci, :, :, cj] = Cs2d[:, :, i, :, :, j]
        out["covariance_scale_2d"] = torch.tensor(full_s, dtype=torch.float32)
    if cc_cols:
        sv = [SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C]
        Cc = _phi_collapse(C[n_sa:, n_sa:], cc_cols)
        cov_c = np.zeros((n_eta, 2, n_eta, 2), dtype=np.float64)
        for i, ci in enumerate(cc_cols):
            for j, cj in enumerate(cc_cols):
                cov_c[:, ci, :, cj] = Cc[:, i, :, j] * sv[ci] * sv[cj]   # → physical
        cc2 = cov_c.reshape(n_eta * 2, n_eta * 2)         # symmetrise ULP asymmetry
        cov_c = (0.5 * (cc2 + cc2.T)).reshape(n_eta, 2, n_eta, 2)
        out["covariance_smear_24_2_24_2"] = torch.tensor(cov_c, dtype=torch.float32).view(n_eta, 2, n_eta, 2)
        sig_sm = torch.sqrt(torch.clamp(
            torch.tensor(np.einsum('icic->ic', cov_c)), min=0.0)).float()
        out["sigma_smear_eff_24_2"] = sig_sm
        Cc2d = C[n_sa:, n_sa:].reshape(n_eta, nphi, len(cc_cols),
                                       n_eta, nphi, len(cc_cols))
        full_c = np.zeros((n_eta, nphi, 2, n_eta, nphi, 2), dtype=np.float64)
        for i, ci in enumerate(cc_cols):
            for j, cj in enumerate(cc_cols):
                full_c[:, :, ci, :, :, cj] = Cc2d[:, :, i, :, :, j] * sv[ci] * sv[cj]
        out["covariance_smear_2d"] = torch.tensor(full_c, dtype=torch.float32)
    path = os.path.join(args.output, "empirical_fisher.pt")
    torch.save(out, path)
    print(f"  wrote {path}: output-space {method} Fisher ({n_act}-param η×φ table, "
          f"n_phi={nphi}) → per-η φ-mean σ ({layout['seen']:,} events, "
          f"Σw={layout['sw']:.2e}) in {time.time()-t0:.1f}s")
    if "sigma_scale_24_3" in out:
        ss = out["sigma_scale_24_3"]
        print(f"  σ(A,e,M) median over bins = "
              f"({float(ss[:,0].median()):.2e}, {float(ss[:,1].median()):.2e}, "
              f"{float(ss[:,2].median()):.2e})")


def _run_empirical_fisher(args, model, shard_files, stats, device) -> None:
    """Joint (θ,φ) empirical Fisher → ``<output>/empirical_fisher.pt`` (PSD,
    background-included θ covariance via pinv; diagnostics-compatible)."""
    if not (model.scale_enabled or model.smearing_enabled):
        print("skipping empirical Fisher: both --disable-scale and --disable-smearing.")
        return
    if model.theta_mode == "mlp":
        if getattr(args, "output_fisher", False):
            _run_output_fisher_mlp(args, model, shard_files, stats, device)
        else:
            _run_empirical_fisher_mlp(args, model, shard_files, stats, device)
        return
    half = _validation_half(args, "fit")
    inj = _inject_theta_np(args, len(stats.eta_edges) - 1) if args.validation else None
    inj_sm = _inject_smear_np(args, len(stats.eta_edges) - 1) if args.validation else None
    # If max_events is below one batch, shrink the loader's batch_size so we
    # neither (a) compute data_nll_continuity over a 65k-event batch only to
    # use 20k of it, nor (b) loop over more per-event score chunks than the
    # cap allows. The chunked score loop dominates the cost (O(n_active)×
    # per-event), so this matters; the inner per-batch truncation below is a
    # belt-and-braces guard for the case where max_events falls mid-batch.
    fisher_bs = args.batch_size
    if args.empirical_fisher_max_events > 0:
        fisher_bs = min(fisher_bs, max(1, args.empirical_fisher_max_events))
    loader = JpsiMassArrowLoader(
        shard_files, stats, batch_size=fisher_bs, split=args.fisher_split,
        val_fraction=0.0, holdout_fraction=0.0,   # all events, matching the all-events fit
        drop_last=False, half=half, inject_theta_scale=inj,
        inject_theta_smear=inj_sm, inject_seed=int(args.inject_smear_seed),
        cond_basis=getattr(args, "cond_basis", "muon_kin"),
        max_events=int(getattr(args, "max_events", 0) or 0),
        event_fraction=float(getattr(args, "event_fraction", 1.0) or 1.0))
    print(f"\ncomputing joint (θ,φ) empirical Fisher (per-event scores → pinv) on "
          f"split={args.fisher_split}"
          + ("  half=%s (MC pseudo-data)" % ('all' if half is None else half)
             if args.validation else "")
          + f"  [smear_fit={model.smear_fit_params}]"
          + (f"; ≤{args.empirical_fisher_max_events:,} events"
             if args.empirical_fisher_max_events > 0 else ""))
    t0 = time.time()
    J, layout = compute_empirical_fisher_joint(
        model, loader, device, mc_as_data=args.validation,
        n_iter=args.continuity_n_iter,
        chunk_events=args.empirical_fisher_chunk,
        max_events=args.empirical_fisher_max_events,
        progress=args.progress, vectorized=args.fisher_vectorized)
    # If capped, scale J to full statistics (J ∝ Σw): a cheap weight-only pass
    # for Σw_total, then J ← J · Σw_total/Σw_seen so the covariance has the
    # correct 1/N_fit scale.
    if layout["hit_cap"] and layout["sw"] > 0:
        sw_total = 0.0
        for batch in loader:
            dm = (~batch["is_data_mask"] if args.validation else batch["is_data_mask"])
            sw_total += float((batch["w"] * dm.to(batch["w"].dtype)).sum())
        if sw_total > layout["sw"]:
            scale = sw_total / layout["sw"]
            J = J * scale
            print(f"  scaled J by Σw_total/Σw_seen = {scale:.2f} "
                  f"(subsampled {layout['seen']:,} events)")
            layout["sw"] = sw_total
    jd = layout["n_theta_active"] + layout["n_mlp"]
    rank = int(torch.linalg.matrix_rank(J).item())
    n_t = layout["n_theta_active"]
    ridge = float(args.empirical_fisher_ridge)
    cov_theta = _empirical_cov_theta_block(J, n_t, ridge)
    out = {
        "method": (f"joint (theta,phi) empirical Fisher, ridge={ridge:g}"
                   if ridge > 0 else "joint (theta,phi) empirical Fisher, pinv"),
        "ridge": ridge,
        "covariance": cov_theta,
        "labels": _active_param_labels(model, layout["smear_cols"]),
        "n_scale": layout["n_scale"], "smear_cols": layout["smear_cols"],
        "smear_fit_params": model.smear_fit_params,
        "n_events": layout["seen"], "sum_weight": layout["sw"],
        "n_theta_active": n_t, "n_mlp": layout["n_mlp"],
        "joint_dim": jd, "joint_rank": rank,
        "param_space": "scale: linear (A,e,M); smear: raw pre-softplus theta_smear",
    }
    out.update(_theta_cov_extras(cov_theta, model, layout["smear_cols"], layout["n_scale"]))
    path = os.path.join(args.output, "empirical_fisher.pt")
    torch.save(out, path)
    print(f"  wrote {path}: joint {jd}×{jd} (θ:{n_t}, φ:{layout['n_mlp']}), "
          f"rank {rank}/{jd}; θ-block {n_t}×{n_t} via "
          f"{'ridge=' + format(ridge, 'g') if ridge > 0 else 'pinv'} "
          f"({layout['seen']:,} events, Σw={layout['sw']:.2e}) in {time.time()-t0:.1f}s")
    if "sigma_scale_24_3" in out:
        ss = out["sigma_scale_24_3"]
        print(f"  empirical σ(A,e,M) median over bins = "
              f"({float(ss[:,0].median()):.2e}, {float(ss[:,1].median()):.2e}, "
              f"{float(ss[:,2].median()):.2e})")


def _load_full_fit(args, device):
    """Load a FULL stage-2 fit (flow + MLP + θ) from ``--checkpoint`` for
    ``--stage uncertainties``. Adopts the model-defining settings (flow arch,
    MLP size, smear-fit choice, scale/smear enables, validation) and the
    standardisation stats from the checkpoint so the rebuilt model matches the
    saved weights exactly. Returns ``(path, stats, model)`` or ``None``."""
    ck_path = args.checkpoint or os.path.join(args.output, "fit_best.pt")
    if not os.path.exists(ck_path):
        print(f"error: --stage uncertainties needs a fitted checkpoint; {ck_path!r} "
              f"not found (point --checkpoint at a fit_best.pt / fit_last.pt).",
              file=sys.stderr)
        return None
    print(f"loading fit checkpoint: {ck_path}")
    ck = torch.load(ck_path, map_location=device, weights_only=False)
    ck_args = ck.get("args", {}) or {}
    if ck_args.get("stage") == "flow":
        print("  warning: --checkpoint is a stage-1 FLOW checkpoint (θ not fit); "
              "the uncertainty will be evaluated at the un-fit θ.", file=sys.stderr)
    # Adopt the model-defining settings from the checkpoint.
    _apply_flow_arch_from_ckpt(args, ck_args)
    for k in ("mlp_hidden", "mlp_n_layers", "smear_fit_params", "scale_fit_params",
              "smear_flow_steps", "smear_operator", "n_gh_nodes",
              "jacobian_form", "smear_param_form",
              "norm_correction", "no_background",
              "qop_floor_frac", "theta_mlp", "theta_mlp_hidden", "theta_mlp_layers",
              "cond_basis",
              "theta_whiten", "theta_whiten_max_rho",
              "disable_scale", "disable_smearing", "validation",
              "no_validation_split",
              "inject_A", "inject_e", "inject_M",
              "inject_a", "inject_c", "inject_smear_seed", "inject_nonuniform"):
        if k in ck_args:
            setattr(args, k, ck_args[k])
    # Stats: --stats-in overrides; else the fit's own stats.
    if args.stats_in is not None and os.path.exists(args.stats_in):
        with open(args.stats_in) as f:
            stats = _stats_from_dict(json.load(f))
        print(f"  preproc stats from {args.stats_in}")
    elif "stats" in ck:
        stats = _stats_from_dict(ck["stats"])
        print("  preproc stats from the checkpoint")
    else:
        print("error: no stats in checkpoint and no --stats-in given.", file=sys.stderr)
        return None
    model = _build_model(args, stats, device)
    model.load_state_dict(ck["state_dict"])
    model.eval()
    return ck_path, stats, model


def _run_uncertainties_stage(args, device) -> int:
    """--stage uncertainties: load an existing fit and run only the Fisher info
    and/or the warm-start bootstrap on it (no training)."""
    loaded = _load_full_fit(args, device)
    if loaded is None:
        return 1
    ck_path, stats, model = loaded
    print(f"=== stage uncertainties: full fit loaded (flow + MLP + θ) — "
          f"scale={'on' if model.scale_enabled else 'off'}, "
          f"smear={'on' if model.smearing_enabled else 'off'}, "
          f"smear_fit={model.smear_fit_params}"
          + ("  [validation: MC pseudo-data]" if args.validation else "") + " ===")
    shard_files = discover_shards(args.inputs)
    if not shard_files:
        print("error: no Arrow shards found under inputs", file=sys.stderr)
        return 1
    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, "preproc_stats.json"), "w") as f:
        json.dump(_stats_to_dict(stats), f, indent=2)
    if (not args.fisher_info and not args.empirical_fisher
            and not getattr(args, "output_fisher", False) and args.bootstrap <= 0):
        print("warning: --stage uncertainties but none of --fisher-info / "
              "--empirical-fisher / --output-fisher / --bootstrap (>0) requested "
              "— nothing to compute.", file=sys.stderr)
        return 0
    if args.fisher_info:
        _run_fisher_continuity(args, model, shard_files, stats, device)
    if args.empirical_fisher or getattr(args, "output_fisher", False):
        _run_empirical_fisher(args, model, shard_files, stats, device)
    if args.bootstrap > 0:
        run_bootstrap_continuity(args, model, shard_files, stats, device,
                                 mc_as_data=args.validation)
    return 0

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="train_jpsi_mass_fit",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "inputs",
        nargs="+",
        help="Arrow file(s) or director(y/ies) of .arrow files produced "
        "by jpsi_mass_fit_snapshot.py (MC and data combined).",
    )
    p.add_argument("--output", required=True, help="Output directory.")
    p.add_argument(
        "--stats-in", default=None,
        help="Path to a precomputed preproc-stats JSON; if unset, stats are "
        "computed from the input shards.",
    )
    p.add_argument(
        "--device", default=("cuda:0" if torch.cuda.is_available() else "cpu"),
        help="Torch device.",
    )
    # Two-stage continuity pipeline.
    p.add_argument(
        "--stage", choices=["both", "flow", "fit", "uncertainties"],
        default="both",
        help="Two-stage continuity training: 'flow' = stage 1 (nominal flow on "
        "simulation, no θ conditioning); 'fit' = stage 2 (freeze flow, fit θ + "
        "background on data via the analytic continuity tilt); 'both' = run 1 "
        "then 2 in-process (default); 'uncertainties' = load an existing FULL fit "
        "from --checkpoint and run only the Fisher info (--fisher-info) and/or "
        "warm-start bootstrap (--bootstrap), no training.",
    )
    p.add_argument(
        "--checkpoint", type=str, default=None,
        help="Checkpoint to load. '--stage fit': a stage-1 FLOW checkpoint "
        "(default <output>/flow_best.pt); only its flow weights are loaded (MLP + θ "
        "start fresh). '--stage uncertainties': a FULL fit checkpoint "
        "(default <output>/fit_best.pt; point it at fit_last.pt to use the latest) "
        "— flow + MLP + θ are all loaded. Either way the flow architecture and the "
        "preproc stats are read from the checkpoint automatically (so they need not "
        "be re-specified); pass --stats-in to override the stats.",
    )
    p.add_argument(
        "--validation", action="store_true",
        help="MC-closure validation mode: use simulation for BOTH stages "
        "instead of data for stage 2. A deterministic disjoint half of the "
        "simulation events trains the stage-1 flow (half 0); the other half "
        "(half 1) is treated as pseudo-data for the stage-2 θ fit. The disjoint "
        "split prevents stage 2 from fitting θ against the events stage 1 "
        "trained on; the closure target is θ → 0. Any real data in the inputs "
        "is unused.",
    )
    p.add_argument(
        "--no-validation-split", action="store_true",
        help="In --validation mode, SKIP the half-split: stage 1 and stage 2 "
        "(and any --fisher-info / --empirical-fisher / --bootstrap runs) all "
        "use ALL simulation events instead of disjoint halves. Doubles the "
        "stats available to each stage at the cost of stage 2 fitting θ on "
        "the same events stage 1 trained on — the closure-target uncertainty "
        "from the flow's stat error effectively collapses to zero, so this "
        "is for closure-machinery testing (recovering an injection) rather "
        "than for representative real-data uncertainty estimates. Default OFF "
        "(disjoint-half split, the standard MC-closure setup).",
    )
    # Inject a known θ_scale shift into the validation pseudo-data (closure with
    # a non-zero target): the stage-2 pseudo-data m_ll is advected by this scale,
    # so the fit should recover it. A constant shift per component over all η
    # bins. Only active with --validation.
    p.add_argument("--inject-A", type=float, default=0.0,
                   help="(--validation) Inject this constant A scale shift into "
                   "the pseudo-data; the fit should recover it (closure target).")
    p.add_argument("--inject-e", type=float, default=0.0,
                   help="(--validation) Inject this constant e [GeV] scale shift.")
    p.add_argument("--inject-M", type=float, default=0.0,
                   help="(--validation) Inject this constant M scale shift.")
    p.add_argument("--inject-a", type=float, default=0.0,
                   help="(--validation) Inject this constant PHYSICAL qop-variance "
                   "coefficient 'a' (σ²_qop = a + c·k², the constant hit-resolution "
                   "term; physical scale ~1e-7) into the pseudo-data via the same "
                   "per-muon qop fold as the validation plots — a Gaussian "
                   "√(a+c·k²) kick, m_ll recomputed. Physical units throughout "
                   "(the O(1) optimizer rescaling is internal); the fit recovers "
                   "the injected value (shown on the θ_smear plot).")
    p.add_argument("--inject-c", type=float, default=0.0,
                   help="(--validation) Inject this constant PHYSICAL qop-variance "
                   "coefficient 'c' (the ∝k²=1/pt² multiple-scattering term; "
                   "physical scale ~1e-6; see --inject-a).")
    p.add_argument("--inject-smear-seed", type=int, default=12345,
                   help="Seed for the injected-smear Gaussian qop kick, so the "
                   "pseudo-data realisation is reproducible across epochs/runs.")
    p.add_argument("--inject-nonuniform", action="store_true",
                   help="(--validation) Make the injected θ NON-UNIFORM across the "
                   "detector: multiply each injected term by a quadratic-in-η × "
                   "sinusoidal-in-φ factor (the φ sinusoid does ~2 oscillations "
                   "over 2π). Each factor varies ~±50% and is ~zero-mean (φ "
                   "exactly; η centred over uniform η), so the detector-average "
                   "injection stays ≈ the nominal --inject-* values. Tests the "
                   "η/φ-dependent (MLP) θ recovery; the constant --inject-* are "
                   "the base values that get modulated.")
    p.add_argument("--flow-epochs", type=int, default=0,
                   help="Max epochs for stage 1 (0 → use --epochs).")
    p.add_argument("--fit-epochs", type=int, default=0,
                   help="Max epochs for stage 2 (0 → use --epochs).")
    p.add_argument("--fit-scale-lr", type=float, default=0.1,
                   help="Stage-2 Adam lr for θ_scale. θ_scale is O(1) (physical "
                   "A,e,M = θ·THETA_SCALE_REF=(1e-4,1e-3,1e-5)), so all three "
                   "components share a well-conditioned step at this O(1) lr. "
                   "(Adam's step is ~lr in the param's units, so this matches the "
                   "physical convergence of the old physical-θ lr 1e-5 = 0.1·REF_A; "
                   "1e-3 was ~100× too slow and left θ crawling short of the optimum.)")
    p.add_argument("--fit-smear-lr", type=float, default=0.1,
                   help="Stage-2 Adam lr for θ_smear (O(1); σ²_qop = θ·SMEAR_VAR_SCALE, "
                   "with SCALE_C calibrated so a realistic c → θ_c≈1). Matches "
                   "--fit-scale-lr: the smear signal is weak (small per-event "
                   "gradient S/N), so an O(1) lr is needed to reach the optimum in "
                   "a reasonable number of steps.")
    p.add_argument("--init-theta-a", type=float, default=0.0,
                   help="Initial value for the RAW θ_smear[:, 0] ('a' column), "
                   "broadcast to all η-bins. Linear form: physical a_init = θ·SMEAR_VAR_"
                   "SCALE_A. Softplus form: physical a_init = softplus(θ)·SMEAR_VAR_SCALE_A — "
                   "so the default θ=0 starts at softplus(0)=0.69, only 0.7σ from the "
                   "negative-tail saturation knee; raise to e.g. 1.0 to keep the "
                   "parameter well inside softplus's active region (where its gradient "
                   "is not vanishing). Square form: physical a_init = θ²·SMEAR_VAR_"
                   "SCALE_A, so θ=0 → 0 exactly (identity init, no offset needed). The "
                   "frozen column per --smear-fit-params is masked to 0 in the forward "
                   "pass regardless of init.")
    p.add_argument("--init-theta-c", type=float, default=0.0,
                   help="Initial value for the RAW θ_smear[:, 1] ('c' column), "
                   "broadcast to all η-bins. See --init-theta-a for the linear / "
                   "softplus / square interpretation; in 'softplus' mode start at e.g. "
                   "1.0–2.0 to avoid edge-bin softplus saturation traps. In 'square' "
                   "mode the default 0 already maps to physical 0.")
    p.add_argument("--fit-mlp-lr", type=float, default=1e-3,
                   help="Stage-2 Adam lr for the background-fraction MLP.")
    p.add_argument("--fit-theta-mlp-lr", type=float, default=1e-3,
                   help="Stage-2 Adam lr for the θ ThetaNet (--theta-mlp). One lr "
                   "for all of (A,e,M,a,c); the net's output reference scaling "
                   "sets the relative A,e,M vs a,c magnitudes.")
    p.add_argument("--fit-optimizer",
                   choices=("adam", "soap", "lbfgs", "adam+lbfgs", "soap+lbfgs",
                            "trust-krylov", "trust-ncg", "trust-exact",
                            "adam+trust-krylov", "soap+trust-krylov",
                            "adam+trust-ncg", "soap+trust-ncg",
                            "adam+trust-exact", "soap+trust-exact"),
                   default="adam",
                   help="Stage-2 (θ + background) optimizer. 'adam' (default): "
                   "torch.optim.Adam, the historical choice. 'soap': SOAP "
                   "(Shampoo in the Adam eigenbasis, from pytorch_optimizer) — a "
                   "Kronecker-factored curvature preconditioner that whitens each "
                   "parameter tensor's gradient, so ill-scaled/correlated "
                   "directions (the θ_net and background MLP weights especially) "
                   "converge more completely toward the local minimum with less "
                   "per-group lr tuning. weight_decay is forced to 0 (a "
                   "calibration fit must not be pulled toward θ=0) and "
                   "precondition_1d is forced True (SOAP defaults it False, which "
                   "would leave every 1-D bias — incl. the θ_net final-layer "
                   "(A,e,M,a,c) bias — un-preconditioned). Preconditions WITHIN "
                   "each tensor, so it complements (not replaces) --theta-whiten, "
                   "which couples the small cross-tensor degeneracies (A/e, a/c, "
                   "scale/smear); the two can be combined. Per-group lrs / plateau "
                   "schedule / bootstrap reset are unchanged. Requires the "
                   "pytorch_optimizer package. 'lbfgs': torch.optim.LBFGS with a "
                   "strong-Wolfe line search — a quasi-Newton FULL-BATCH minimiser "
                   "for exhaustive descent on the (deterministic) stage-2 NLL: "
                   "each epoch runs one optim.step(closure) of up to "
                   "--lbfgs-max-iter inner iterations, re-evaluating the "
                   "weighted-mean NLL+grad over the WHOLE loader, converging to a "
                   "tiny gradient norm (reported as |g|=… in place of lr). Best "
                   "for the small, degenerate BINNED θ fit; the line search "
                   "auto-scales the step so the per-group lrs / plateau schedule "
                   "don't apply. The bootstrap refit (minibatch) falls back to "
                   "Adam. Pairs well with --theta-whiten. 'trust-krylov' / "
                   "'trust-ncg' / 'trust-exact': scipy second-order trust region "
                   "(full-sample exact gradient, exact Hessian/HVP over "
                   "--hess-subsample-events). trust-krylov/ncg are Hessian-FREE "
                   "(HVP); trust-exact builds the full Hessian (binned θ only, "
                   "see --trust-exact-max-par). HYBRIDS "
                   "'{adam,soap}+{lbfgs,trust-krylov,trust-ncg,trust-exact}': "
                   "two-phase — Adam/SOAP to its normal stopping (robust bulk "
                   "descent; escapes the near-flat θ≈0 start that stalls L-BFGS's "
                   "line search), THEN the second-order polish warm-started from "
                   "it (L-BFGS for --lbfgs-final-epochs, or the trust-region "
                   "driver for --fit-epochs iters) to reach a tight gradient norm "
                   "Adam/SOAP — normalising by the gradient RMS — structurally "
                   "cannot. If the polish regresses, the phase-1 result is "
                   "restored.")
    p.add_argument("--lbfgs-lr", type=float, default=1.0,
                   help="(--fit-optimizer lbfgs) Initial step scale; with the "
                   "strong-Wolfe line search 1.0 is standard (the search rescales "
                   "it). The O(1) THETA_SCALE_REF/SMEAR_VAR_SCALE rescaling keeps "
                   "θ well-conditioned for this.")
    p.add_argument("--lbfgs-max-iter", type=int, default=20,
                   help="(--fit-optimizer lbfgs) Quasi-Newton inner iterations per "
                   "epoch (per optim.step(closure)); each is one full-batch "
                   "NLL+grad pass. Total work ≈ --fit-epochs × this.")
    p.add_argument("--lbfgs-history-size", type=int, default=20,
                   help="(--fit-optimizer lbfgs) Number of (s,y) curvature pairs "
                   "kept for the inverse-Hessian approximation.")
    p.add_argument("--lbfgs-tolerance-grad", type=float, default=1e-9,
                   help="(--fit-optimizer lbfgs) Gradient-norm convergence "
                   "tolerance for the inner iterations (exhaustive → small).")
    p.add_argument("--lbfgs-tolerance-change", type=float, default=1e-12,
                   help="(--fit-optimizer lbfgs) Min objective/param change per "
                   "inner step before the line search declares convergence.")
    p.add_argument("--lbfgs-final-epochs", type=int, default=1,
                   help="(--fit-optimizer adam+lbfgs / soap+lbfgs) Number of "
                   "L-BFGS polish epochs in phase 2 (each = one optim.step(closure) "
                   "of up to --lbfgs-max-iter inner iterations), warm-started from "
                   "the phase-1 Adam/SOAP result.")
    # SciPy second-order trust region (trust-krylov / trust-ncg / trust-exact).
    # trust-krylov/ncg are Hessian-FREE (HVP, --trust-hvp); trust-exact builds
    # the FULL n_par×n_par Hessian per iteration (vectorised second-backward,
    # float64) and solves the subproblem exactly — best for the small binned θ
    # (few params), blocked above --trust-exact-max-par (e.g. --theta-mlp).
    p.add_argument("--trust-gtol", type=float, default=1e-6,
                   help="(--fit-optimizer trust-*) Gradient-norm convergence "
                   "tolerance (the exhaustive-minimisation stopping criterion; "
                   "the gradient is exact on the FULL sample). NOTE: the float32 "
                   "model evals floor ‖g‖ at ~1e-4..1e-5, so a much smaller gtol "
                   "just runs to the iter cap.")
    p.add_argument("--trust-exact-max-par", type=int, default=400,
                   help="(--fit-optimizer trust-exact) Max number of fitted "
                   "parameters for which the full Hessian is built; above this "
                   "trust-exact errors out (use trust-krylov instead). Binned θ "
                   "is well under; --theta-mlp (~1000s of net weights) is not.")
    p.add_argument("--trust-hvp", choices=("reuse", "recompute"), default="reuse",
                   help="(trust-krylov/ncg) Exact (true-Hessian) HVP scheme — H is "
                   "never formed. 'reuse' (scheme A, default): build the "
                   "differentiable gradient g_S once per step (create_graph) and "
                   "RETAIN its graph; each Krylov HVP is one second backward, NO "
                   "data pass/forward per HVP — cheapest, but the whole subset "
                   "graph stays resident (memory ∝ subset, not chunkable). "
                   "'recompute' (scheme B): re-do forward+double-backward per HVP, "
                   "chunked over events (--trust-hvp-chunk) so peak memory is "
                   "bounded. (Ignored by trust-exact, which builds the full H.)")
    p.add_argument("--hess-subsample-events", type=int, default=0,
                   help="(trust-*) Cap on events used for the HVP/Hessian subset "
                   "(0 = full sample). The GRADIENT is always exact on the full "
                   "sample; subsampling only the Hessian is safe because the "
                   "full-objective trust ratio test corrects a noisy H_S (it only "
                   "costs iterations, never convergence). The subset is FIXED "
                   "across all HVPs of a step (operator consistency for Krylov).")
    p.add_argument("--trust-hvp-chunk", type=int, default=0,
                   help="(--trust-hvp recompute / trust-ncg,krylov) Events per "
                   "forward+double-backward chunk. 0 (default) = WHOLE BATCH (no "
                   "sub-batch chunking): a double-backward HVP peaks at only ~2× a "
                   "plain gradient, so if a full-batch gradient fits the full-batch "
                   "HVP fits too — sub-chunking is pure overhead for the same "
                   "compute. Set >0 only if genuinely memory-bound.")
    p.add_argument("--soap-precondition-frequency", type=int, default=10,
                   help="(--fit-optimizer soap) Optimizer steps between SOAP's "
                   "preconditioner eigendecompositions. Cheap here (tiny θ "
                   "tensors); cost is dominated by data_nll_continuity.")
    p.add_argument("--soap-shampoo-beta", type=float, default=-1.0,
                   help="(--fit-optimizer soap) EMA decay for the Shampoo "
                   "preconditioner. <0 (default) → use SOAP's own default (the "
                   "second β of --betas).")
    p.add_argument("--soap-eps", type=float, default=1e-8,
                   help="(--fit-optimizer soap) Denominator floor on the "
                   "second-moment √(exp_avg_sq) in SOAP's ROTATED (preconditioned) "
                   "space — effectively a ridge on the preconditioner: LARGER eps "
                   "→ less aggressive whitening of small-curvature (sloppy / "
                   "near-degenerate, e.g. A/e) directions → more Adam-like and "
                   "more conservative steps along flat directions; smaller → "
                   "stronger whitening. Default 1e-8 (SOAP's own default).")
    p.add_argument("--continuity-n-iter", type=int, default=2,
                   help="Fixed-point iterations for the #2 source solve "
                   "(advection+smear pre-image).")
    p.add_argument("--epochs", type=int, default=50, help="Maximum training epochs.")
    p.add_argument("--batch-size", type=int, default=65536, help="Events per batch.")
    p.add_argument("--lr", type=float, default=1e-3,
                   help="Adam lr for flow + MLP.")
    p.add_argument("--weight-decay", type=float, default=0.0,
                   help="Adam weight decay (L2) on all optimized parameters.")
    p.add_argument("--patience", type=int, default=8,
                   help="Early-stop after this many epochs without val improvement.")
    p.add_argument("--patience-threshold", type=float, default=1e-4,
                   help="Minimum val-NLL decrease that counts as an improvement.")
    p.add_argument("--no-early-stop", action="store_true",
                   help="Disable early stopping (train the full --epochs).")
    p.add_argument("--lr-schedule", choices=["plateau", "cosine", "none"],
                   default="plateau",
                   help="LR schedule (both stages). 'plateau': reduce "
                   "lr on val plateau (ReduceLROnPlateau); 'cosine': decay to "
                   "--min-lr over --epochs; 'none': fixed lr.")
    p.add_argument("--lr-reduce-factor", type=float, default=0.3,
                   help="Plateau schedule: multiply lr by this on a plateau.")
    p.add_argument("--lr-reduce-patience", type=int, default=4,
                   help="Plateau schedule: epochs without val improvement before "
                   "reducing lr (keep < --patience so lr drops before early-stop).")
    p.add_argument("--min-lr", type=float, default=1e-7,
                   help="Floor on the scheduled lr; early-stop fires once the lr "
                   "has reached this (plateau) reductions are exhausted.")
    p.add_argument("--val-fraction", type=float, default=0.10,
                   help="Fraction of events held out for validation.")
    p.add_argument("--holdout-fraction", type=float, default=0.05,
                   help="Fraction held out from train+val (e.g. for Fisher info).")
    p.add_argument("--max-events", type=int, default=0,
                   help="Subsample to ~this many events for the flow + fit stages "
                   "(0 = use all). Applied per shard AFTER the --validation "
                   "half-split and BEFORE the train/val/holdout split, so each "
                   "half is independently capped and the splits stay proportional. "
                   "Rows are pre-shuffled, so this is an unbiased subset.")
    p.add_argument("--event-fraction", type=float, default=1.0,
                   help="Keep this fraction of events per shard for the flow + fit "
                   "stages (1.0 = all), composed the same way as --max-events "
                   "(after the validation half-split, before the split; the "
                   "tighter of the two wins per shard).")
    p.add_argument("--m-lo", type=float, default=2.92, dest="m_lo",
                   help="Lower edge of the m_ll fit window [GeV].")
    p.add_argument("--m-hi", type=float, default=3.28, dest="m_hi",
                   help="Upper edge of the m_ll fit window [GeV].")
    p.add_argument("--n-eta-bins", type=int, default=None,
                   help="Number of η bins for the BINNED θ parameters (and the "
                   "diagnostics' per-η tables), spanning ±--eta-range. Default 24 "
                   "(uniform over ±2.4). The model's binned θ tables, the per-event "
                   "η-bin index, the Fisher/bootstrap per-η covariance, and the "
                   "diagnostics all derive their bin count from the resulting "
                   "stats.eta_edges. The binned-θ η binning is independent of the "
                   "flow (which conditions on continuous η), so passing this in the "
                   "fit stage rebuilds stats.eta_edges even when stats come from a "
                   "flow --checkpoint or --stats-in. Has no effect on --theta-mlp "
                   "(continuous in η).")
    p.add_argument("--eta-range", type=float, default=None,
                   help="Half-range of the η binning: bins span [−R, +R] with "
                   "--n-eta-bins uniform bins (default R=2.4, the muon acceptance).")
    # Scale transform + smear init + noise sampling
    p.add_argument(
        "--qop-floor-frac", type=float, default=0.0,
        help="DEPRECATED / inert. The qop→pt inversion is now pt = |sinθ/qop| "
        "with only the qop=0 pole guarded (QOP_EPS): pt is a magnitude, so a "
        "kick large enough to flip the sign of qop is kept as the physical "
        "charge mis-reconstruction it is, not floored away. Accepted for "
        "checkpoint/back-compat but no longer affects the fold.",
    )
    p.add_argument(
        "--scale-fit-params", default="AM",
        help="Which per-η-bin SCALE terms to FIT — a subset of 'AeM' (e.g. 'AM', "
        "'A', 'AeM'). A (constant δqop) and e (∝1/pt) are NEARLY DEGENERATE over "
        "the narrow J/ψ pt range, so fitting both from J/ψ alone is ill-posed — "
        "the fit slides into large opposite-sign (A,e). Default 'AM' drops e (the "
        "J/ψ-identifiable subset: constant scale A + charge-odd sagitta M). For a "
        "single-parameter closure, fit only the injected term (e.g. inject-e → "
        "--scale-fit-params e). Dropped terms are held inert at 0.",
    )
    p.add_argument(
        "--smear-fit-params", choices=["both", "a", "c"], default="both",
        help="Which per-η-bin smear terms to FIT: 'both', 'a' (constant term "
        "only), or 'c' (∝1/pt² term only). The constant a and the c·k² term are "
        "nearly degenerate over the narrow J/ψ pt range, so fitting both per "
        "bin is ill-posed and yields the unphysical bin-to-bin zig-zag (use 'a' "
        "or 'c' to break it). The non-fitted term is zeroed — removed from the "
        "σ_qop variance entirely, not floated.",
    )
    p.add_argument(
        "--smear-flow-steps", type=int, default=1,
        help="Euler steps integrating the smear's probability-flow ODE in the "
        "density (the score-driven deterministic-diffusion change of variables). "
        "1 = first-order (single score displacement); more steps integrate the "
        "score flow more finely (more robust/accurate) at a higher nested-"
        "autograd cost. The per-muon qop fold (closure plots) and the injection "
        "are exact convolutions regardless. Ignored when --smear-operator="
        "gh_convolution.",
    )
    p.add_argument(
        "--smear-operator",
        choices=("pf_ode", "gh_convolution", "gh_convolution_qop"), default="pf_ode",
        help="Smear operator for the continuity density. 'pf_ode' (default): "
        "deterministic probability-flow ODE — fast but over-broadens at large "
        "V/σ² (the model density is wider than the actual Gaussian convolution "
        "the pseudo-data uses, biasing the fit to under-recover θ_smear at "
        "forward |η|). 'gh_convolution': stochastic Gaussian convolution in MASS "
        "space via Gauss-Hermite quadrature (--n-gh-nodes nodes) — the small-kick "
        "linearisation; accurate at small σ_qop/|qop| but under-represents the "
        "skewed/heavy mass tails at forward |η| + large smearing (so it still "
        "under-recovers c there). 'gh_convolution_qop': EXACT per-muon qop smear "
        "as a 2-D GH over the two independent muon kicks — kicks each qop and "
        "recomputes m_ll, reproducing the true non-Gaussian mass smearing the "
        "pseudo-data fold uses at ALL V/|η| (incl. the qop floor). Cost: n_gh² "
        "per-node evals (each muon kick is one 1-D Gaussian, so n_gh≈5–6 is "
        "enough → 25–36 nodes). Propagates the smear to the conditioning ρ "
        "per node. Requires V ≥ 0 (clamps internally; pair with "
        "--smear-param-form softplus for the cleanest guarantee).",
    )
    p.add_argument(
        "--n-gh-nodes", type=int, default=8,
        help="Number of Gauss-Hermite quadrature nodes PER smearing dimension "
        "(ignored for pf_ode). For 'gh_convolution' (1-D, mass space) 8 is "
        "plenty (truncation error is exponential in n_gh for smooth p_0). For "
        "'gh_convolution_qop' the cost is n_gh² (a 2-D grid over the two muon "
        "kicks) — each kick is a single 1-D Gaussian, so 4–5 (16–25 nodes) is "
        "usually enough; large n_gh here (8 → 64 nodes) multiplies the flow-eval "
        "rows and memory accordingly (flow evals are auto-chunked to bound the "
        "single-allocation size, but total memory still scales with n_gh²).",
    )
    p.add_argument(
        "--norm-correction", choices=("none", "linear", "flow_cdf"), default="none",
        help="Per-event normalisation correction for the transformed-flow "
        "density. The forward map T_θ is NOT boundary-preserving on [m_lo, "
        "m_hi]: broadening pushes some probability mass outside the window, so "
        "Z(θ;c) = ∫_window p_θ(x|c) dx < 1 and `log p_θ(x)` carries a -log Z "
        "bias that pulls the fit toward smaller broadening. "
        "'none' (default): no correction, preserves current behaviour. "
        "'linear': leading-order boundary expansion (PF-ODE forward map) — 2 "
        "forward + 2 flow evals per event, valid for small V. "
        "'flow_cdf': EXACT via Z = F_0(T⁻¹(m_hi)) - F_0(T⁻¹(m_lo)), the GF "
        "flow's CDF at boundary preimages — ~2× the per-event work of the bare "
        "density. NOTE: when --smear-operator=gh_convolution, BOTH 'linear' and "
        "'flow_cdf' route to the GH-specific exact formula Z = Σ_i W_i [F_0(m'_"
        "hi(ξ_i)) - F_0(m'_lo(ξ_i))] (per-GH-node boundary inversion); the "
        "PF-ODE-based forms above use the wrong operator at the boundary.")
    p.add_argument(
        "--smear-param-form", choices=("linear", "softplus", "square"),
        default="linear",
        help="Positivity reparameterisation for θ_smear. 'linear' (default): "
        "the O(1) θ_smear is the coefficient directly — signed, supports both "
        "broadening (V>0) and unsmear (V<0). 'softplus': each of (a, c) "
        "INDIVIDUALLY constrained to ≥ 0 via physical = softplus(θ)·SMEAR_VAR_"
        "SCALE; the per-η fit mask is applied AFTER softplus so frozen "
        "params (or 'a'-/'c'-only modes) remain EXACTLY zero in the per-muon "
        "σ_qop and downstream transformations. Use 'softplus' to defend "
        "against the negative-c drift (issue #2) by construction, at the cost "
        "of losing the two-sided fit (the model can no longer represent MC "
        "that's too broad vs data). The Fisher σ accounts for the sigmoid(θ̂) "
        "delta-method Jacobian in 'softplus'. 'square': physical = θ²·SMEAR_"
        "VAR_SCALE — same individual positivity as softplus but better "
        "convergence (θ=0 → 0 exactly, like 'linear', so identity init with no "
        "softplus ln2 offset / saturation knee, and a non-saturating Jacobian "
        "2·SMEAR_VAR_SCALE·θ). Caveat: θ↔−θ degenerate and ∂physical/∂θ→0 at "
        "θ=0, so the raw-θ covariance (--fisher-info / --empirical-fisher / "
        "--bootstrap) is singular for a smear bin pinned at zero — take the "
        "smear σ from --output-fisher there instead. --init-theta-{a,c} sets a "
        "nonzero raw θ init (interpreted in the chosen form).",
    )
    p.add_argument(
        "--jacobian-form", choices=("softlog", "exp"), default="softlog",
        help="Smear-Jacobian formula for the continuity density. 'softlog' "
        "(default): the exact autograd Jacobian G' = dx/dm' of the N-step Euler "
        "forward map, with a C¹ linear-tangent extension below SMEAR_GP_FLOOR "
        "(keeps the gradient alive past the floor → the optimiser is pulled BACK "
        "from G' < 0 instead of drifting into the previously-flat floored basin "
        "where unphysical θ_smear<0 was rewarded). 'exp': frozen-score "
        "continuous-flow approximation log G' = log(1+s_adv') − V·∂²log p₀/2 "
        "(always finite, no floor) — but a DIFFERENT operator approximation "
        "(assumes ∂²log p₀ constant along the smear trajectory), and rewards "
        "unphysical sharpening UNBOUNDEDLY by −log G' = +V·∂²/2; prefer "
        "'softlog' unless you specifically want the closed-form continuous-flow "
        "Jacobian. See _continuity_logp for the full caveats.",
    )
    p.add_argument(
        "--theta-mlp", action="store_true",
        help="Replace the per-η-bin (A,e,M,a,c) tables with a small MLP mapping "
        "each muon's (η, φ) → (A,e,M,a,c) CONTINUOUSLY (trained in stage 2 like "
        "the background MLP). Note: the observed/empirical Fisher and bootstrap "
        "uncertainties are binned-θ-only and are skipped in this mode.",
    )
    p.add_argument(
        "--cond-basis", choices=["muon_kin", "event_level"], default="muon_kin",
        help="Flow + background-MLP conditioning basis. 'muon_kin' (default): "
        "the leak-free per-muon (η±, cosφ±, sinφ±, ρ). 'event_level': the "
        "dilepton vars (yll, ln ptll, cosPhill, sinPhill, cosθ*, sinφ*, cosφ*) — "
        "all pt-dependent, so the qop scale/smear propagates to the whole "
        "conditioning (and the background MLP then sees ln ptll). event_level "
        "with smearing requires --smear-operator gh_convolution_qop; the basis "
        "must be the SAME for stage 1 and stage 2.")
    p.add_argument("--theta-mlp-hidden", type=int, default=32,
                   help="(--theta-mlp) Hidden width of the θ ThetaNet.")
    p.add_argument("--theta-mlp-layers", type=int, default=2,
                   help="(--theta-mlp) Number of hidden layers of the θ ThetaNet.")
    p.add_argument(
        "--theta-whiten", action="store_true",
        help="Decorrelate the degenerate (A,e) and (a,c) parameter pairs with a "
        "fixed gradient-whitening preconditioner built from the per-η-bin "
        "curvature moments (⟨k⟩, ⟨k²⟩, ⟨k⁴⟩, k=1/pt). The (A,e) scale response "
        "(1,−k) and (a,c) variance response (1,k²) are near-collinear over the "
        "narrow J/ψ pt range, giving a tilted loss valley that Adam (a diagonal "
        "preconditioner) cannot navigate — forcing one term of each pair to be "
        "dropped. This rotates the gradient of each pair into a unit-conditioned "
        "basis so BOTH terms can be fit jointly. Backward-only: the forward map, "
        "likelihood, and observed Fisher are unchanged (it is inert during "
        "eval/Fisher/diagnostics). Only acts where both members of a pair are "
        "fit — pair it with --scale-fit-params AeM and/or --smear-fit-params both. "
        "Works in both binned (per-η-bin rotation) and --theta-mlp (single global "
        "rotation, no η encoding) modes. Requires stats with k_moments (recomputed "
        "automatically; older stats.json disable it with a warning).",
    )
    p.add_argument(
        "--theta-whiten-max-rho", type=float, default=0.99,
        help="(--theta-whiten) Cap on the |correlation| used to build each "
        "whitening matrix (bounds how hard the near-degenerate direction is "
        "amplified — cond(C⁻¹)=(1+ρ)/(1−ρ)). Default 0.99.",
    )
    # θ_scale sampling widths, split per component (A, e, M) since they live
    # in different physical units. Each is the σ of the Gaussian added to that
    # component (additive, physical units) — the fixed width with
    # --fixed-theta-sampling / --no-adaptive-sigma, and the adaptive-σ fallback
    # during warmup.
    # Flow / MLP hyperparams
    p.add_argument(
        "--flow-arch", choices=("gf", "nsf"), default="gf",
        help="Signal flow architecture: 'gf' = Gaussianization flow (default) — "
        "C∞-smooth density, so the continuity score/Hessian have no knot kinks. "
        "'nsf' = neural rational-quadratic spline flow — bounded (linear tails "
        "outside ±5), avoids the erf/exp saturation the GF needs guards for, but "
        "only C¹ (the score kinks at the spline knots).",
    )
    p.add_argument(
        "--nsf-bins", type=int, default=8,
        help="(--flow-arch nsf only) Number of rational-quadratic spline "
        "knots per transform.",
    )
    p.add_argument("--flow-n-transforms", type=int, default=5,
                   help="Number of stacked flow transforms.")
    p.add_argument("--flow-hidden", type=int, default=128,
                   help="Hidden width of each flow conditioner MLP.")
    p.add_argument("--flow-n-hidden", type=int, default=3,
                   help="Number of hidden layers in each flow conditioner MLP.")
    p.add_argument("--gf-components", type=int, default=8,
                   help="(--flow-arch gf only) Gaussian-mixture components per layer.")
    p.add_argument("--mlp-hidden", type=int, default=32,
                   help="Hidden width of the data-branch mixture MLP.")
    p.add_argument("--mlp-n-layers", type=int, default=2,
                   help="Number of hidden layers in the mixture MLP.")
    # Fisher
    p.add_argument("--fisher-info", action="store_true",
                   help="After training, compute + save the observed Fisher "
                   "information (Hessian + covariance) → <output>/fisher_info.pt. "
                   "Two-stage pipeline: over θ_scale + the ACTIVE θ_smear column(s) "
                   "jointly, with the flow and background MLP held fixed (fixed-φ / "
                   "conditional). Legacy pipeline: θ_scale only.")
    p.add_argument("--fisher-split", default="train",
                   choices=("train", "val", "holdout", "all"),
                   help="(two-stage) Loader split the Fisher info is summed over. "
                   "Default 'train' = the data the fit used, so the covariance has "
                   "the correct statistical scale (∝ 1/N_fit). In --validation mode "
                   "the half-1 MC pseudo-data of this split is used. NOTE: the fit "
                   "and the Fisher/bootstrap use ALL events (no val/holdout split), "
                   "so 'train' == 'all' here and 'val'/'holdout' are EMPTY — keep "
                   "the default 'train'.")
    p.add_argument("--fisher-vectorized", default=True,
                   action=argparse.BooleanOptionalAction,
                   help="(two-stage) Compute the Hessian with one vmapped "
                   "(is_grads_batched) second backward instead of a per-row loop. "
                   "Faster but holds ~n_param copies of the backward graph; "
                   "automatically falls back to the loop on engine failure / OOM. "
                   "Use --no-fisher-vectorized to force the loop. Also governs the "
                   "per-event score in --empirical-fisher and the observed-Hessian "
                   "rows in --output-fisher-method observed/sandwich.")
    # Joint (theta, phi) empirical Fisher: J = Σ w_i s_i s_iᵀ over θ + the MLP,
    # PSD → pinv covariance (background-included, Hessian-free, no negative σ).
    p.add_argument("--empirical-fisher", action="store_true",
                   help="(two-stage) After the fit, compute the JOINT (θ, MLP) "
                   "empirical Fisher J = Σ w_i s_i s_iᵀ from per-event scores and "
                   "write the θ-block of pinv(J) → <output>/empirical_fisher.pt. "
                   "PSD by construction (no negative/clamped σ), background "
                   "uncertainty included (the MLP is in the joint), and Hessian-"
                   "free (first derivatives only). Cost is O(N_events) per-event "
                   "gradients — see --empirical-fisher-max-events.")
    p.add_argument("--empirical-fisher-max-events", type=int, default=0,
                   help="Cap on events used for the empirical Fisher (default 0 = "
                   "all events). If set >0, J (an average per-event quantity) is "
                   "computed on that representative subset and rescaled by "
                   "Σw_total/Σw_seen to the full-statistics covariance scale.")
    p.add_argument("--empirical-fisher-chunk", type=int, default=64,
                   help="Events per is_grads_batched call when extracting "
                   "per-event scores (memory ≈ chunk × backward graph).")
    p.add_argument("--empirical-fisher-ridge", type=float, default=0.0,
                   help="Ridge for the empirical-Fisher covariance. 0 (default) → "
                   "Moore–Penrose pinv: unconstrained / degenerate directions get "
                   "σ=0 (e.g. the A/e near-degeneracy over the J/ψ pt range reads "
                   "as σ(e)≈0). >0 → scale-aware ridge cov=(J + ridge·diag(J))⁻¹, "
                   "so a degenerate direction reads as a LARGE finite σ instead. "
                   "Dimensionless (relative to each parameter's own information); "
                   "try ~1e-3.")
    # OUTPUT-space Fisher for --theta-mlp: measure the information on the per-(η,φ)
    # OUTPUTS directly via a 2-D η×φ θ-table shadow, instead of the 1349-weight
    # net + weight-space ridge. Physical degeneracy structure + physical ridge
    # scale; observed / sandwich available (the net empirical Fisher is empirical-
    # only). The φ-mean per-η σ band is recovered by averaging the table cov over φ.
    p.add_argument("--output-fisher", action="store_true",
                   help="(--theta-mlp) After the fit, compute the OUTPUT-space "
                   "Fisher over a 2-D η×φ θ-table seeded from the net at the bin "
                   "centres (physical A/e/M; O(1) a/c) → <output>/empirical_fisher.pt "
                   "(same keys the diagnostics σ-band consumes). Measures the data's "
                   "information about the per-(η,φ) OUTPUTS directly, so the ridge "
                   "scale and degeneracy structure are physical (output space), not "
                   "the 1349 net weights. The per-η φ-mean covariance is the table "
                   "covariance averaged over its φ axes. Falls back to the net "
                   "empirical Fisher when off. Ignored for binned θ.")
    p.add_argument("--output-fisher-method", default="empirical",
                   choices=("observed", "empirical", "sandwich"),
                   help="(--output-fisher) Information estimator over the η×φ table: "
                   "'empirical' J=Σw·ssᵀ (per-event scores, default); 'observed' "
                   "H=Σw·∂²(−lnL)/∂θ² (Hessian, double-backward); 'sandwich' "
                   "H⁻¹JH⁻¹ (robust). Empirical is the cheapest; observed/sandwich "
                   "cost a per-batch Hessian block over the active table entries.")
    p.add_argument("--output-fisher-nphi", type=int, default=16,
                   help="(--output-fisher) Number of φ bins in the η×φ table "
                   "(default 16). More bins resolve φ structure (and give a finer "
                   "φ-resolved Fisher band) but split the statistics per cell; the "
                   "per-η φ-collapse + the φ-band average over them, occupancy-"
                   "weighted by the per-cell Σw.")
    p.add_argument("--output-fisher-project", default="free",
                   choices=("free", "net"),
                   help="(--output-fisher) Output-covariance subspace. 'free' "
                   "(default): every (η,φ) cell is an independent parameter — "
                   "model-agnostic and well-conditioned, but BLIND to the cross-"
                   "bin (smoothness) correlations the MLP induces, so the φ-mean "
                   "band averages down like 1/n_φ (often too small). 'net': "
                   "restrict the covariance to the network-REACHABLE output "
                   "subspace (SVD of the table-output Jacobian G = ∂o/∂w), so the "
                   "band carries the MLP's smoothness correlations — the network-"
                   "consistent covariance. With no ridge/truncation this equals "
                   "the pseudoinverse limit of the weight-propagated covariance, "
                   "but computed in a clean low-dim basis (no rank-deficient "
                   "weight-space ridge lottery).")
    p.add_argument("--output-fisher-svd-rtol", type=float, default=1e-2,
                   help="(--output-fisher --output-fisher-project net) Relative "
                   "singular-value cutoff σ_k/σ_max for the net-subspace "
                   "projection — the SMOOTHNESS CUTOFF: output directions the net "
                   "can only reach with large weight excursions (σ_k below this) "
                   "are dropped. Default 1e-2; smaller keeps more (rougher) "
                   "network modes.")
    p.add_argument("--output-fisher-marginalize-bkg", default="auto",
                   choices=("auto", "on", "off"),
                   help="(--output-fisher) Marginalise the background-fraction "
                   "MLP into the θ uncertainty by floating its weights in the "
                   "JOINT output-space information and Schur-complementing them "
                   "out (so f_data(c) is profiled, not held fixed). 'auto' "
                   "(default): ON for BOTH projections — making the band "
                   "comparable to the net-weight / bootstrap bands, which also "
                   "marginalise the background. 'off' gives the fixed-background "
                   "conditional band; 'on' is the same as auto. The bkg Schur "
                   "block uses the same --empirical-fisher-ridge (it only damps "
                   "the PSD inflation term, so over-ridging falls back to the "
                   "fixed-bkg band). NOTE: with project=free the per-cell inverse "
                   "is already ill-conditioned, so the marginalised free band is "
                   "the noisiest combination — project=net is recommended.")
    # Warm-start Poisson bootstrap (Hessian-free covariance incl. the background)
    p.add_argument("--bootstrap", type=int, default=0,
                   help="(two-stage) After the stage-2 fit, run this many warm-start "
                   "Poisson-bootstrap replicas → <output>/bootstrap_cov.pt. Each "
                   "replica restarts from the nominal (θ̂, φ̂), re-fits θ + the "
                   "background MLP jointly on the data Poisson(1)-reweighted, and the "
                   "covariance of {θ̂_b} folds in the background uncertainty with no "
                   "Hessian. 0 = off.")
    p.add_argument("--bootstrap-epochs", type=int, default=None,
                   help="Max epochs per bootstrap replica. Default (None) inherits "
                   "the nominal stage-2 cap (--fit-epochs or --epochs). Replicas "
                   "hitting the cap without plateauing are flagged as possibly "
                   "under-converged.")
    p.add_argument("--bootstrap-patience", type=int, default=None,
                   help="Per-replica early-stop patience (epochs without "
                   "reweighted-NLL improvement > --patience-threshold). Default "
                   "(None) inherits the nominal fit's --patience, so each replica "
                   "converges with the SAME early-stop / threshold / LR-schedule "
                   "(--lr-schedule + --lr-reduce-*) as the nominal fit — monitored "
                   "on the replica's reweighted NLL (a replica has no separate val).")
    p.add_argument("--bootstrap-seed", type=int, default=12345,
                   help="Base seed for the per-replica Poisson(1) event weights "
                   "(replica b uses seed + b; reset each epoch so an event keeps its "
                   "count across the replica's epochs).")
    # Mixed-precision + torch.compile (same convention as
    # train_muon_response_flow.py).
    p.add_argument(
        "--compile", action="store_true",
        help="Wrap the model in torch.compile. Currently SKIPPED — the "
        "GF base flow's Jacobian uses an inner torch.autograd.grad call "
        "(zuko's MonotonicTransform.call_and_ladj) that dynamo cannot "
        "trace. Flag is accepted for parity with the response-flow "
        "trainer; a one-line note is printed at startup.",
    )
    p.add_argument(
        "--precision", choices=("fp32", "bf16", "fp16"), default="fp32",
        help="Training + validation forward-pass precision. fp32 = no "
        "autocast. bf16 = bfloat16 autocast (Ampere+, no GradScaler). "
        "fp16 = float16 autocast + enabled GradScaler for loss scaling.",
    )
    p.add_argument(
        "--progress", default=True, action=argparse.BooleanOptionalAction,
        help="Show a per-epoch tqdm progress bar with running NLL "
        "(--no-progress to disable).",
    )
    p.add_argument(
        "--profile-steps", type=int, default=0,
        help="If >0, time the first N training steps of each epoch and print a "
        "per-phase breakdown (data+H2D / forward / backward+step), to diagnose "
        "whether training is CPU/dataloader-bound vs compute-bound. On CUDA each "
        "phase is bracketed with torch.cuda.synchronize() (so async kernels are "
        "attributed correctly) and a GPU util/peak-mem readout is appended; "
        "reading — data≫compute → data/host-bound, compute high + util high → "
        "GPU-bound, compute high + util low → launch/sync-bound. 0 = off.",
    )
    # Adaptive σ from Adam's second moment.
    # Adaptive-σ clamps are split scale/smear because the two parameters live
    # in different units: θ_scale is the linear (A,e,M) ~1e-4…1e-2, while
    # θ_smear is the *raw* (pre-softplus) param whose natural sampling scale is
    # O(1). The floor prevents σ→0 collapse near sharp optima; the ceiling
    # prevents σ→∞ on parameters that haven't received gradient yet (v_i = 0).
    p.add_argument(
        "--val-seed", type=int, default=42,
        help="Fixed RNG seed used at the start of every validation pass. "
        "Keeps the σ̃_smear noise (and the per-event smearing-kernel ε) "
        "deterministic across epochs so val_nll has no Monte Carlo "
        "fluctuation and can be compared apples-to-apples with train_nll.",
    )
    # Debug flags — bisect the model down to the smallest piece that
    # reproduces the failure (NaN, instability, etc.).
    p.add_argument(
        "--disable-smearing", action="store_true",
        help="Drop the residual-smearing kernel and the θ_smear "
        "conditioning of the flow. θ_smear is excluded from the optimizer "
        "and never sampled (it stays an inert, fixed Parameter).",
    )
    p.add_argument(
        "--disable-scale", action="store_true",
        help="Drop the T_scale forward-fold and the θ_scale conditioning of "
        "the flow. θ_scale is excluded from the optimizer and never sampled "
        "(inert, fixed at 0); Fisher info is skipped. With BOTH "
        "--disable-scale and --disable-smearing, only the flow and the MLP "
        "(background normalisation) are trained.",
    )
    p.add_argument(
        "--no-background", action="store_true",
        help="Disable the data-branch background mixture: the data NLL "
        "reduces to pure signal (−log p_signal), the f_data MLP is bypassed "
        "and its parameters are excluded from the stage-2 optimiser. Intended "
        "for VALIDATION CLOSURES where the truth f_bkg = 0 by construction "
        "(MC pseudo-data is signal-only) — removes the degeneracy where the "
        "MLP grows f_bkg in forward |η| bins to absorb tail events the signal "
        "model can't broaden into (and to absorb injection-induced "
        "out-of-window pollution that the Bernstein basis extrapolates to). "
        "Leave OFF for real-data fits, which need the MLP for genuine bkg.",
    )
    p.add_argument(
        "--detect-anomaly", action="store_true",
        help="Wrap training in torch.autograd.set_detect_anomaly so the "
        "first NaN in the backward graph throws with a stack trace "
        "pointing to the offending op. Slow — debug only.",
    )
    p.add_argument(
        "--nan-on-step", default="raise", choices=("raise", "skip", "ignore"),
        help="What to do when the per-batch loss is NaN/Inf. 'raise' "
        "stops training immediately with the batch index. 'skip' drops "
        "the batch and continues. 'ignore' lets the NaN propagate (the "
        "original behaviour).",
    )
    return p.parse_args(argv)


def main() -> int:
    return train_loop(parse_args())


if __name__ == "__main__":
    sys.exit(main())
