#!/usr/bin/env python
"""Decompose the baseline ∂NLL/∂A gradient of the single-η-bin J/ψ scale fit
into its three model contributions, to localise the A-closure bias:

    NLL = -logp0(m_t) - logJ(m_t) + logZ        (per event, weighted)
    ∂NLL/∂A = g_p0 + g_jac + g_norm
      g_p0   = ∂/∂A Σw·(-logp0)   flow DENSITY at the un-kicked source mass
      g_jac  = ∂/∂A Σw·(-logJ)    change-of-variables JACOBIAN |∂m_t/∂x|
      g_norm = ∂/∂A Σw·(+logZ)    truncated-window NORMALISATION (flow_cdf)

At A=0 on the NOMINAL (un-injected) data — the flow's own training events — a
perfect (truncated-MLE) flow gives total ≈ 0. A non-zero total = E_model[score]
− E_data[score] is a direct measure of the frozen flow failing to reproduce the
nominal data in the scale direction; this script says WHICH term carries it
(flow template vs Jacobian vs normalisation). In practice the flow density
(g_p0) dominates, and it is reproduced even by a pure additive mass shift
(--shift-mode additive, g_jac=0), pinning the bias on the flow template — see the
stage-1 --no-flow-window-norm note in train_jpsi_mass_fit.py (a full-support
vs truncated normalisation mismatch was the identified cause).

Assumes --disable-smearing (single GH node, so logp0 and logJ separate cleanly)
and --smear-operator gh_convolution_qop, matching the testvalidation3x closures.
Runs the model in float64 (the ~6e-5 signal is below the fp32 gradient floor).

Usage (inside the wmassdevrolling container + `source setup.sh`):
    python scripts/corrections/muon_calibration/decompose_A_bias.py \
        --checkpoint /…/fit_last.pt \
        --shards     /…/jpsi_v1/shards
Options: --shift-mode {qop,additive}, --A-raw (scan point in raw θ units,
default 0), --inject (inject A into the pseudo-data, default none = nominal),
--device, --max-batches (quick subset; NB the shards are inhomogeneous so a
prefix's absolute value is unreliable — use the full sample for the headline),
--batch-size, --repo (defaults to this script's directory).
"""
import argparse
import os
import sys
import time

import numpy as np
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--checkpoint", required=True)
ap.add_argument("--shards", required=True)
ap.add_argument("--repo", default=os.path.dirname(os.path.abspath(__file__)),
                help="path to the muon_calibration scripts dir (for imports); "
                     "defaults to this script's own directory.")
ap.add_argument("--A-raw", type=float, default=0.0,
                help="θ_scale[0,0] in RAW units to evaluate the gradient at "
                     "(physical A = A_raw·1e-4). 0 = the nominal-closure probe.")
ap.add_argument("--inject", type=float, default=None,
                help="inject this physical A into the pseudo-data (default none "
                     "= nominal data, the flow's own training events).")
ap.add_argument("--device", default=None)
ap.add_argument("--shift-mode", choices=("qop", "additive"), default="qop",
                help="How A acts on the mass. 'qop' = the real per-muon qop "
                     "scale continuity operator (logp0 + logJ + logZ). "
                     "'additive' = a pure constant mass shift m_src = m + A·M0 "
                     "(M0 = mean J/ψ mass): unit Jacobian, so the decomposition "
                     "is just flow density + window normalisation (g_jac = 0). "
                     "Isolates the flow-template bias with no operator/Jacobian.")
ap.add_argument("--batch-size", type=int, default=65536)
ap.add_argument("--eval-chunk", type=int, default=0,
                help="events per autograd slice within each loader batch (the "
                     "gradient sums are linear in events, so slicing is exact). "
                     "0 = auto: full batch for analytic-CDF flows; 8192 for "
                     "--flow-arch nce, whose quadrature window-Z "
                     "(--nce-quad-nodes MLP evals per event, fp64, activations "
                     "RETAINED for the autograd.grad) exhausts the CUDA "
                     "allocator at the default --batch-size (the NVML assert).")
ap.add_argument("--max-batches", type=int, default=0,
                help="0 = full sample (recommended). >0 truncates (fast but the "
                     "inhomogeneous shards make a prefix non-representative).")
ap.add_argument("--match", choices=("flow", "fit", "all"), default="flow",
                help="Which event set to evaluate on, reproduced from the "
                     "checkpoint's training config (the loader split/half/"
                     "event_fraction are deterministic, so this is the EXACT same "
                     "set). 'flow' (default) = the stage-1 flow TRAIN split "
                     "(split=train, val/holdout fractions and flow half from the "
                     "checkpoint) — so the baseline gradient isolates flow "
                     "in-window fidelity with no train/test statistical gap. "
                     "'fit' = the stage-2 fit events (all events of the fit half, "
                     "no val/holdout) — the actual injected-closure data. "
                     "'all' = every event (legacy behaviour).")
args = ap.parse_args()

sys.path.insert(0, args.repo)
import train_jpsi_mass_fit as T                       # noqa: E402
from jpsi_mass_fit_diagnostics import load_model_from_checkpoint, _move_batch  # noqa: E402
from jpsi_mass_arrow_loader import JpsiMassArrowLoader  # noqa: E402
from jpsi_mass_model import _event_mll, THETA_SCALE_REF  # noqa: E402

dev = args.device or ("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"device={dev}  shift-mode={args.shift_mode}")
model, stats, targs, ck = load_model_from_checkpoint(args.checkpoint, dev)
model.eval()
model.double()                                          # fp64: signal ~6e-5 < fp32 floor
if args.shift_mode == "qop":
    if model.smearing_enabled:
        sys.exit("ERROR: the qop decomposition assumes --disable-smearing "
                 "(single GH node so logp0/logJ separate). Rerun with smearing off.")
    if model.smear_operator != "gh_convolution_qop":
        sys.exit(f"ERROR: expected smear_operator=gh_convolution_qop, "
                 f"got {model.smear_operator}.")

n_eta = len(stats.eta_edges) - 1
inj = None
if args.inject is not None:
    inj = np.zeros((n_eta, 3), dtype=np.float64); inj[:, 0] = args.inject

# Reproduce the training event set from the checkpoint so the baseline gradient
# isolates flow fidelity (no train/test statistical gap). The loader's split /
# half / event_fraction are all deterministic, so matching these args yields the
# EXACT same events the flow (or fit) saw.
validation = bool(targs.get("validation", False))
no_split = bool(targs.get("no_validation_split", False))
def _half(which):  # mirrors train_jpsi_mass_fit._validation_half
    if not validation or no_split:
        return None
    return 0 if which == "flow" else 1
if args.match == "flow":
    half = _half("flow")
    vf = float(targs.get("val_fraction", 0.1))
    hf = float(targs.get("holdout_fraction", 0.05))
    sel = "stage-1 flow TRAIN split"
elif args.match == "fit":
    half = _half("fit"); vf = 0.0; hf = 0.0
    sel = "stage-2 fit events (all of the fit half)"
else:  # all
    half = None; vf = 0.0; hf = 0.0
    sel = "ALL events"
efrac = float(targs.get("event_fraction", 1.0) or 1.0)
print(f"event set: --match={args.match} → {sel}  "
      f"(split=train half={half} val_fraction={vf:g} holdout_fraction={hf:g} "
      f"event_fraction={efrac:g})")
loader = JpsiMassArrowLoader(
    T.discover_shards([args.shards]), stats, batch_size=args.batch_size,
    split="train", val_fraction=vf, holdout_fraction=hf, drop_last=False,
    half=half, event_fraction=efrac,
    inject_theta_scale=inj, inject_theta_smear=None,
    inject_seed=int(targs.get("inject_smear_seed", 12345)),
    cond_basis=targs.get("cond_basis", "muon_kin"))

with torch.no_grad():
    model.theta_scale[0, 0] = float(args.A_raw)
ts = model.theta_scale


def parts_logp(m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm):
    """Single-node (smearing-off) gh_qop continuity, returning (logp0, logJ) per
    event [B] — a faithful copy of _continuity_logp_gh_qop's ng=1 path."""
    B = m_obs.shape[0]
    eps = torch.zeros(1, 1, 2, device=m_obs.device, dtype=m_obs.dtype)   # ξ=0
    pto, etao, phio = pt_obs.unsqueeze(1), eta_pm.unsqueeze(1), phi_pm.unsqueeze(1)
    qo, bpo = q_pm.unsqueeze(1), b_pm.unsqueeze(1)
    with torch.enable_grad():
        lam = torch.ones(B, 1, device=m_obs.device, dtype=m_obs.dtype,
                         requires_grad=True)
        pt_cfg = pto * lam.unsqueeze(-1)
        mo_lam = _event_mll(pt_cfg, etao, phio)
        m_t, pt_truth = model._gh_qop_unsmear(pt_cfg, etao, phio, qo, bpo, eps)
        g_mt = torch.autograd.grad(m_t.sum(), lam, create_graph=True,
                                   retain_graph=True)[0]
        g_mo = torch.autograd.grad(mo_lam.sum(), lam, create_graph=True,
                                   retain_graph=True)[0]
    logJ = (torch.log(g_mt.abs().clamp_min(1e-12))
            - torch.log(g_mo.abs().clamp_min(1e-12)))            # [B,1]
    mk_g = mk.unsqueeze(1).expand(B, 1, mk.shape[-1]).clone()
    mk_g = model._node_cond(mk_g, pt_truth, etao, phio, qo)
    logp0 = model._flow_eval_chunked(
        model.log_p_nominal, m_t.reshape(-1), mk_g.reshape(B, -1)).reshape(B, 1)
    return logp0.squeeze(1), logJ.squeeze(1)


M0 = float(stats.mll_mean)                 # constant shift scale (mean J/ψ mass)
REF0 = float(THETA_SCALE_REF[0])           # raw θ → physical A (1e-4)
M_LO, M_HI = float(model.m_lo), float(model.m_hi)


def additive_logp_logZ(m_obs, mk):
    """Pure additive mass shift m_src = m_obs + A_phys·M0 (A_phys = θ_raw·1e-4).
    Unit Jacobian (the shift is independent of the integration variable), so NO
    logJ term. Returns (logp0, logZ); the flow density and the EXACT window
    normalisation Z = F0(m_hi+s) − F0(m_lo+s) are both analytic in the flow.
    Conditioning is held at the observed value (an additive mass shift changes
    no kinematics). Sign matches the qop un-kick (source mass higher for A>0)."""
    s = ts[0, 0] * REF0 * M0                                  # scalar, differentiable
    logp0 = model.log_p_nominal(m_obs + s, mk)                # [B]
    mhi = m_obs.new_full(m_obs.shape, M_HI) + s
    mlo = m_obs.new_full(m_obs.shape, M_LO) + s
    Z = (model._flow_log_cdf(mhi, mk).exp()
         - model._flow_log_cdf(mlo, mk).exp()).clamp_min(1e-30)
    return logp0, Z.log()                                     # [B], [B]


ni = int(targs.get("continuity_n_iter", 2))
# Per-batch autograd slicing (memory bound, not a statistics knob): the three
# gradient sums are LINEAR in events, so per-slice accumulation is exact. The
# nce flow's quadrature window-Z retains --nce-quad-nodes MLP activations per
# event for the backward — at fp64 and the default --batch-size that exhausts
# the CUDA allocator (NVML assert); analytic-CDF flows are fine unsliced.
chunk = args.eval_chunk or (8192 if getattr(model, "flow_is_nce", False)
                            else args.batch_size)
if chunk != args.batch_size:
    print(f"eval-chunk: {chunk} events per autograd slice "
          f"({'auto: nce quadrature CDF' if not args.eval_chunk else 'user'})")
g_p0 = g_jac = g_norm = 0.0
sw = 0.0
t0 = time.time()
nb = 0
for i, b in enumerate(loader):
    if args.max_batches and i >= args.max_batches:
        break
    b = _move_batch(b, dev)
    dm = ~b["is_data_mask"]
    if not bool(dm.any()):
        continue
    cast = lambda x: x.double() if x.is_floating_point() else x
    m_b, pt_b = cast(b["mll"]), cast(b["pt_pm"])
    eta_b, phi_b = cast(b["eta_pm"]), cast(b["phi_pm"])
    q_b, bp_b, mk_b = b["q_pm"], b["b_pm"], cast(b["cond_std"])
    w_b = cast(b["w"]) * dm.double()

    for s0 in range(0, m_b.shape[0], chunk):
        sl = slice(s0, s0 + chunk)
        m, pt, eta, phi = m_b[sl], pt_b[sl], eta_b[sl], phi_b[sl]
        q, bp, mk, w = q_b[sl], bp_b[sl], mk_b[sl], w_b[sl]

        if args.shift_mode == "additive":
            logp0, lZ = additive_logp_logZ(m, mk)
            g_p0 += float(torch.autograd.grad((w * (-logp0)).sum(), ts,
                                              retain_graph=True)[0][0, 0])
            # g_jac stays 0 (unit Jacobian)
            g_norm += float(torch.autograd.grad((w * lZ).sum(), ts)[0][0, 0])
        else:
            logp0, logJ = parts_logp(m, mk, pt, eta, phi, q, bp)
            g_p0 += float(torch.autograd.grad((w * (-logp0)).sum(), ts,
                                              retain_graph=True)[0][0, 0])
            g_jac += float(torch.autograd.grad((w * (-logJ)).sum(), ts)[0][0, 0])
            lZ = model._norm_correction_log_Z(m, mk, pt, eta, phi, q, bp,
                                              n_iter=ni)
            g_norm += float(torch.autograd.grad((w * lZ).sum(), ts)[0][0, 0])

        sw += float(w.sum())
    nb += 1
    if i % 20 == 0:
        print(f"  batch {i:>4}  Σw={sw:.3e}  dt={time.time()-t0:.0f}s", flush=True)

tot = (g_p0 + g_jac + g_norm) / sw
print(f"\nshift-mode={args.shift_mode}  A_raw={args.A_raw}  inject={args.inject}  "
      f"batches={nb}  Σw={sw:.4e}")
print(f"  g_p0   [-logp0 (flow density)] /Σw = {g_p0/sw:+.4e}")
if args.shift_mode == "additive":
    print(f"  g_jac  [-logJ  (Jacobian)]     /Σw = {0.0:+.4e}  (unit J — none)")
else:
    print(f"  g_jac  [-logJ  (Jacobian)]     /Σw = {g_jac/sw:+.4e}")
print(f"  g_norm [+logZ  (normalisation)]/Σw = {g_norm/sw:+.4e}")
print(f"  TOTAL  ∂NLL/∂θ                 /Σw = {tot:+.4e}")
print("  (a non-zero TOTAL at A=0 on nominal data = the baseline A bias; with "
      "g_jac=0 in additive mode it is just flow-density vs window-normalisation.\n"
      "   pre-fix qop full-sample reference: TOTAL≈+6.2e-5 = +1.83e-4(flow) "
      "-1.0e-4(J) -2e-5(Z); retrain the flow without --no-flow-window-norm and "
      "this should shrink toward 0.)")
