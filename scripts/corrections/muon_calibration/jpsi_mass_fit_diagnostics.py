"""Diagnostic plots for the J/ψ unbinned mass-fit calibration.

Inputs:
  --checkpoint <run>/checkpoint_best.pt   model + stats + θ_scale + θ_smear
  --fisher     <run>/fisher_info.pt       Hessian + covariance (optional;
                                          omit → θ ±σ bands skipped)
  --shards     <run>/shards/              per-bucket Arrow files

Outputs (under --output, default ``<checkpoint_dir>/diagnostics/``):
  1. mll_closure_{eta,rho,cosalpha}.png  (+ mc_closure_{eta,rho,cosalpha}.png)
     m_ll histograms (data + MC, weighted) with overlaid model curves
     (signal + Bernstein backgrounds + total mixture). One panel figure per
     slice dimension — columns are the slices: inclusive + |η_+| (eta), and
     the conditional ρ (pt asymmetry) and cos α (opening angle) tertiles.
  2. theta_scale_vs_eta.png
     A, e, M per η-bin with ±1σ error bars (from --fisher: fisher_info.pt or
     bootstrap_cov.pt).
  3. theta_smear_vs_eta.png
     effective a, c per η-bin with ±1σ error bars (from --fisher).
  4. covariance_correlation.png
     full joint covariance (signed-log) + correlation heatmaps over
     θ_scale + the active θ_smear, with the scale/smear block separator.
     (Legacy θ_scale-only files fall back to fisher_correlation.png.)
  5. mll_pulls_inclusive.png  + mll_pulls_eta{0..3}.png
     Per-bin (data − model)/√model histograms; expect ~N(0,1).

The ±1σ bands and the matrix plot require a covariance file via --fisher
(fisher_info.pt from --fisher-info, or bootstrap_cov.pt from --bootstrap).
"""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
from typing import List

import matplotlib

matplotlib.use("Agg")  # noqa: E402 — must precede pyplot
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402
from tqdm import tqdm  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from jpsi_mass_arrow_loader import (  # noqa: E402
    JpsiMassArrowLoader,
    discover_shards,
    _event_cond_raw_np,
    _inject_modulation_eta_np,
    _INJECT_ETA_REF, _INJECT_PHI_NOSC, _INJECT_AMP,
)
from jpsi_mass_model import (  # noqa: E402
    JpsiMassMixtureModel, _event_mll, _event_cond_raw,
    N_THETA_SCALE, N_THETA_SCALE_PM, N_THETA_SMEAR_PM,
    SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C, THETA_SCALE_REF,
    bernstein_basis_n, exp_bkg_density,
)
from train_jpsi_mass_fit import _move_batch, _stats_from_dict  # noqa: E402
from compact_flow import infer_learn_weights  # noqa: E402


# ---------------------------------------------------------------------------
# Checkpoint loading
# ---------------------------------------------------------------------------


def _fit_select_from(targs):
    """(ptll_min, pt_lead_min, pt_both_min) loader selection from checkpoint
    args, or None — selects the analysis sample at diagnostics time so the
    plots/normalisation match the fit (the model also carries these as its
    per-event normalisation edge)."""
    sel = (targs.get("fit_ptll_min"), targs.get("fit_pt_lead_min"),
           targs.get("fit_pt_both_min"))
    if not any(sel):
        return None
    if targs.get("fit_cuts_as_ratio"):
        # Scale-invariant ratio mode: thresholds C/m_lo, cut on pt/m_ll (4th
        # element flags the loader). Matches the fit (no per-event m_min).
        m_lo = float(targs.get("fit_m_lo"))
        return tuple((c / m_lo if c else None) for c in sel) + (True,)
    return sel


def load_model_from_checkpoint(checkpoint_path: str, device: str):
    """Rebuild the trained model from a checkpoint dict."""
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    args = ckpt["args"]
    stats = _stats_from_dict(ckpt["stats"])
    # Per-stage mass windows (default = the shard window in stats): the model
    # window is the FIT window (window-Z, background, plots); the flow window
    # sets the compact/nce/dcb/ege [a, b].
    _fl_lo = float(args.get("flow_m_lo") or stats.m_lo)
    _fl_hi = float(args.get("flow_m_hi") or stats.m_hi)
    _ft_lo = float(args.get("fit_m_lo") or _fl_lo)
    _ft_hi = float(args.get("fit_m_hi") or _fl_hi)
    model = JpsiMassMixtureModel(
        m_lo=_ft_lo,
        m_hi=_ft_hi,
        flow_m_lo=_fl_lo,
        flow_m_hi=_fl_hi,
        mll_log_scale=stats.mll_log_scale,
        mll_mean=stats.mll_mean,
        mll_std=stats.mll_std,
        y_event_mean=torch.from_numpy(stats.y_event_mean),
        y_event_std_tensor=torch.from_numpy(stats.y_event_std),
        muon_kin_mean=torch.from_numpy(stats.muon_kin_mean),
        muon_kin_std_tensor=torch.from_numpy(stats.muon_kin_std),
        flow_arch=args.get("flow_arch", "gf"),
        flow_n_transforms=args["flow_n_transforms"],
        flow_hidden_features=args["flow_hidden"],
        flow_n_hidden_layers=args["flow_n_hidden"],
        flow_gf_components=args["gf_components"],
        flow_nsf_bins=args.get("nsf_bins", 8),
        mlp_hidden=args["mlp_hidden"],
        mlp_n_layers=args["mlp_n_layers"],
        smearing_enabled=not args.get("disable_smearing", False),
        scale_enabled=not args.get("disable_scale", False),
        qop_floor_frac=args.get("qop_floor_frac", 0.0),
        smear_fit_params=args.get("smear_fit_params", "both"),
        scale_fit_params=args.get("scale_fit_params", "AM"),
        smear_flow_steps=args.get("smear_flow_steps", 1),
        smear_operator=args.get("smear_operator", "pf_ode"),
        n_gh_nodes=args.get("n_gh_nodes", 8),
        jacobian_form=args.get("jacobian_form", "softlog"),
        smear_param_form=args.get("smear_param_form", "linear"),
        norm_correction=args.get("norm_correction", "none"),
        background_enabled=not bool(args.get("no_background", False)),
        # Ratio-mode cuts are conditioning-fixed → fixed-window normalisation, so
        # the model carries NO per-event m_min edge (matches the fit).
        fit_ptll_min=(None if args.get("fit_cuts_as_ratio")
                      else args.get("fit_ptll_min")),
        fit_pt_lead_min=(None if args.get("fit_cuts_as_ratio")
                         else args.get("fit_pt_lead_min")),
        fit_pt_both_min=(None if args.get("fit_cuts_as_ratio")
                         else args.get("fit_pt_both_min")),
        bkg_model=args.get("bkg_model", "bernstein"),
        bkg_degree=int(args.get("bkg_degree", 1)),
        theta_mode=("mlp" if args.get("theta_mlp", False)
                    else ("binned2d" if args.get("theta_2d_binned", False)
                          else "binned")),
        n_phi_bins=int(args.get("n_phi_bins", 16)),
        n_eta_bins=len(stats.eta_edges) - 1,
        # compact flow: adopt the mixture-weight mode from the checkpoint;
        # older compact checkpoints (always learnable, 3K heads) lack the key —
        # infer it from the saved conditioner head shape. The head-shape
        # heuristic only applies to LOGISTIC layers (bernstein heads have M
        # outputs, which could collide with 2K/3K — but every bernstein
        # checkpoint postdates the key, so the args.get always hits first).
        compact_learn_weights=(
            args.get("compact_learn_weights",
                     infer_learn_weights(ckpt["state_dict"],
                                         int(args.get("gf_components", 8))))
            if args.get("compact_layer", "logistic") == "logistic" else False),
        compact_layer=args.get("compact_layer", "logistic"),
        bernstein_degree=args.get("bernstein_degree", 16),
        nce_quad_nodes=args.get("nce_quad_nodes", 64),
        cond_basis=args.get("cond_basis", "muon_kin"),
        theta_mlp_hidden=args.get("theta_mlp_hidden", 32),
        theta_mlp_layers=args.get("theta_mlp_layers", 2),
    ).to(device)
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, stats, args, ckpt


# ---------------------------------------------------------------------------
# Model evaluation on a loader: per-event signal/bkg densities on a grid
# ---------------------------------------------------------------------------


@torch.no_grad()
@torch.no_grad()
def _tilt_density_on_grid(
    model, batch, idx, m_centers_dev, *, chunk_events: int = 2048,
    n_iter: int = 2,
) -> torch.Tensor:
    """``[len(idx), n_grid]`` log p_s(m_grid | c_e, θ_fit) — the #2 direct-eval
    signal density (``_continuity_logp``) evaluated at each grid mass, exactly
    what the stage-2 fit optimises. ``m_centers_dev`` is the physical bin grid.
    Per event, pt scales as ``pt·(m_grid/m_obs)`` (pt∝m at fixed conditioning)."""
    n = idx.shape[0]
    G = m_centers_dev.shape[0]
    eta = batch["eta_pm"][idx]; q = batch["q_pm"][idx]; b = batch["b_pm"][idx]
    phi = batch["phi_pm"][idx]
    # For event_level the operator recomputes the conditioning from the
    # grid-scaled pt_g (a pure dilation, under which the dimensionless basis is
    # invariant — the recompute is then a consistency no-op); the passed mk is
    # used only for muon_kin's pt-invariant η/φ.
    mk = batch["cond_std"][idx]; pt = batch["pt_pm"][idx]
    m_obs = batch["mll"][idx]
    out = torch.empty((n, G), device=m_centers_dev.device, dtype=mk.dtype)
    for start in range(0, n, max(1, chunk_events)):
        end = min(start + chunk_events, n); sub = end - start
        mg = m_centers_dev.view(1, G).expand(sub, G)                # [sub,G]
        scale = (mg / m_obs[start:end].view(sub, 1)).unsqueeze(-1)  # [sub,G,1]
        rep = lambda x: x[start:end].unsqueeze(1).expand(
            sub, G, *x.shape[1:]).reshape(sub * G, *x.shape[1:])
        pt_g = (pt[start:end].unsqueeze(1) * scale).reshape(sub * G, 2)
        # Refine the pt scale so the RECOMPUTED event mass hits the grid mass
        # EXACTLY: the operator rebuilds m from the per-muon kinematics, and
        # the muon-mass term breaks the m ∝ pt proportionality — the naive
        # scale misses the grid point by up to ~3 MeV at the far grid edge
        # (δm² = B·(1−(m_g/m_obs)²), B ~ m_μ²·(p₁/p₂+p₂/p₁+2)), i.e. an O(50%)
        # density error on the steep window edges (it made the θ=0 tilt curve
        # visibly disagree with the nominal overlay on flow-checkpoint
        # closures). Two fixed-point steps leave a sub-keV residual.
        eta_r, phi_r = rep(eta), rep(phi)
        mg_flat = mg.reshape(-1)
        for _ in range(2):
            m_cur = _event_mll(pt_g.unsqueeze(1), eta_r.unsqueeze(1),
                               phi_r.unsqueeze(1)).squeeze(1)
            pt_g = pt_g * (mg_flat / m_cur).unsqueeze(-1)
        lp = model._continuity_logp(
            mg_flat, rep(mk), pt_g, eta_r, phi_r, rep(q), rep(b),
            n_iter=n_iter)
        out[start:end] = lp.reshape(sub, G)
    # Window-normalize per event, exactly as the fit does in data_nll_continuity:
    # subtract logZ = log ∫_window p_θ (a per-event scalar, constant over the mass
    # grid), so the displayed signal density integrates to 1 over [m_lo,m_hi] and
    # the closure curve matches the data. Essential once the stage-1 flow is
    # window-normalized (the truncated training no longer pins ∫_window p₀ to 1,
    # so the bare density is off by Z); a no-op for a full-support flow (Z≈1).
    if getattr(model, "norm_correction", "none") != "none":
        logZ = model._norm_correction_log_Z(
            m_obs, mk, pt, eta, phi, q, b, n_iter=n_iter)            # [n]
        out = out - logZ.view(n, 1)
    return out


@torch.no_grad()
def _nominal_density_on_grid(model, batch, idx, m_centers_dev, *, chunk_events=4096):
    """``[len(idx), n_grid]`` log p₀(m_grid | c) — the *nominal* (θ=0) flow
    density, i.e. the stage-1 template with no scale/smear correction. Point
    evaluations of the frozen flow at the grid masses.

    For ``event_level`` the conditioning is pt-dependent, so it is recomputed
    per grid column from the grid-scaled pt (pt∝m) — matching how the tilt curve
    sweeps the conditioning, so the nominal/tilt overlay is consistent. For
    ``muon_kin`` the conditioning is invariant under the common pt scaling, so it
    is held fixed (the legacy fast path)."""
    n = idx.shape[0]; G = m_centers_dev.shape[0]
    mk = batch["cond_std"][idx]
    out = torch.empty((n, G), device=m_centers_dev.device, dtype=mk.dtype)
    event_level = getattr(model, "cond_basis", "muon_kin") == "event_level"
    if event_level:
        pt = batch["pt_pm"][idx]; eta = batch["eta_pm"][idx]
        phi = batch["phi_pm"][idx]; q = batch["q_pm"][idx]
        m_obs = batch["mll"][idx]
    for start in range(0, n, max(1, chunk_events)):
        end = min(start + chunk_events, n); sub = end - start
        mg = m_centers_dev.view(1, G).expand(sub, G)                # [sub,G]
        if event_level:
            scale = (mg / m_obs[start:end].view(sub, 1)).unsqueeze(-1)
            rep = lambda x: x[start:end].unsqueeze(1).expand(
                sub, G, *x.shape[1:]).reshape(sub * G, *x.shape[1:])
            pt_g = (pt[start:end].unsqueeze(1) * scale).reshape(sub * G, 2)
            # Same muon-mass refinement as the tilt grid (see
            # _tilt_density_on_grid) so the conditioning sweep is consistent.
            eta_r, phi_r = rep(eta), rep(phi)
            mg_flat = mg.reshape(-1)
            for _ in range(2):
                m_cur = _event_mll(pt_g.unsqueeze(1), eta_r.unsqueeze(1),
                                   phi_r.unsqueeze(1)).squeeze(1)
                pt_g = pt_g * (mg_flat / m_cur).unsqueeze(-1)
            mke = model._cond_from_muons(pt_g, eta_r, phi_r, rep(q))
        else:
            mke = mk[start:end].unsqueeze(1).expand(
                sub, G, mk.shape[-1]).reshape(sub * G, -1)
        out[start:end] = model.log_p_nominal(mg.reshape(-1), mke).reshape(sub, G)
    # Window-normalize the θ=0 nominal template the same way stage-1 trains it
    # (truncated): subtract logZ_window = log(F0(m_hi|c) − F0(m_lo|c)), the
    # per-event window mass of the frozen flow (constant over the grid), so the
    # nominal curve is a proper density on [m_lo,m_hi]. No-op for a full-support
    # flow (Z≈1); required once the flow is window-normalized (Z≠1). Uses the
    # observed (θ=0) conditioning.
    if getattr(model, "norm_correction", "none") != "none":
        m_hi = mk.new_full((n,), float(model.m_hi))
        m_lo = mk.new_full((n,), float(model.m_lo))
        # STABLE log window mass: the naive exp-difference of CDFs collapses
        # to the clamp once the (gauge-free) window mass drifts below the fp32
        # floor — the nominal curve then explodes by e⁶⁹·Z_true (the e21-scale
        # mc_closure pathology at forward η).
        logZ = model._flow_log_window_Z(m_lo, m_hi, mk)
        out = out - logZ.view(n, 1)
    return out


@torch.no_grad()
def _injected_raw_theta(model, inject_scale_phys, inject_smear_phys):
    """Convert physical injection values to RAW θ tensors (shape ``[n_eta, 3]``
    and ``[n_eta, 2]``) that reproduce them through the model's parameterisation.
    Components that were NOT injected are filled with zero physical, i.e. the
    truth in the closure is θ=0 for those (not the fitted value).

    * scale: physical (A, e, M) = θ · THETA_SCALE_REF (linear), so the raw θ
      that gives the injected physical value is ``inject_phys / REF``.
    * smear:
        - ``smear_param_form='linear'``: effective(θ) = θ; raw = inject_phys /
          SMEAR_VAR_SCALE.
        - ``smear_param_form='softplus'``: effective(θ) = softplus(θ); invert
          via raw = log(exp(eff) − 1) (= softplus⁻¹). For zero (un-injected)
          physical, eff is clipped to 1e-10 so the raw value is large-negative
          (~ −23) with softplus(.) ≈ 0 — effectively the no-smear identity.
        - ``smear_param_form='square'``: effective(θ) = θ²; invert via
          raw = √(eff) (positive root; θ↔−θ are equivalent). Zero injection →
          raw 0 exactly (no clipping needed)."""
    n_eta = model.theta_scale.shape[0]
    if inject_scale_phys is None:
        inject_scale_phys = np.zeros((n_eta, 3), dtype=np.float64)
    if inject_smear_phys is None:
        inject_smear_phys = np.zeros((n_eta, 2), dtype=np.float64)
    ref = np.asarray(THETA_SCALE_REF, dtype=np.float64)
    raw_scale = torch.tensor(inject_scale_phys / ref[None, :], dtype=torch.float32)
    scale = np.asarray([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C], dtype=np.float64)
    eff = inject_smear_phys / scale[None, :]   # target effective θ_smear
    form = getattr(model, "smear_param_form", "linear")
    if form == "softplus":
        # softplus⁻¹(y) = log(exp(y) − 1) = log(expm1(y)); clip eff at 1e-10 so
        # zero-injection columns give a large-negative raw value (softplus→0).
        eff_safe = np.clip(eff, 1e-10, None)
        raw_smear = torch.tensor(np.log(np.expm1(eff_safe)), dtype=torch.float32)
    elif form == "square":
        # square⁻¹(y) = √y (positive root); zero injection → raw 0 exactly.
        raw_smear = torch.tensor(np.sqrt(np.clip(eff, 0.0, None)), dtype=torch.float32)
    else:
        raw_smear = torch.tensor(eff, dtype=torch.float32)
    return raw_scale, raw_smear


def _degeneracy_eigbasis(k_moments):
    """Global (A,e) and (a,c) degeneracy eigenbases from the curvature moments.

    Built from the O(1)-space Gauss–Newton Fisher of each fitted pair, with the
    response weighted by the SAME reference scales the parameters carry (so the
    eigenvectors live in the O(1) θ coordinates the fit moves in):

      scale (θ_A, θ_e): response (REF_A·1, REF_e·k)  → M_s
      smear (θ_a, θ_c): response (SCALE_A·1, SCALE_C·k²) → M_c

    using the global moments ⟨k⟩,⟨k²⟩,⟨k⁴⟩. (Using the bare correlation
    [[1,ρ],[ρ,1]] instead — a unit-diagonal approximation — slightly rotates
    the eigenvectors and can spuriously put a (θ_A=θ_e) injection at stiff=0;
    the REF-weighted diagonal is the correct measured/degenerate split.) The
    eigenvector with the LARGE eigenvalue is the STIFF (well-measured)
    combination, the small one is the SLOPPY (degenerate) one.

    ``k_moments``: [n_eta, 4] sums (N, Σk, Σk², Σk⁴). Returns
    ``(E_scale, λ_scale, E_smear, λ_smear)`` with each ``E[:, 0]`` = stiff,
    ``E[:, 1]`` = sloppy and λ NORMALISED so the sloppy eigenvalue = 1 (so
    λ_stiff is the stiff/sloppy information ratio); or None if unavailable."""
    if k_moments is None:
        return None
    km = np.asarray(k_moments, dtype=np.float64).reshape(-1, 4)
    Ntot = max(float(km[:, 0].sum()), 1.0)
    k1, k2, k4 = km[:, 1].sum() / Ntot, km[:, 2].sum() / Ntot, km[:, 3].sum() / Ntot
    rA, re_, _ = THETA_SCALE_REF
    sA, sC = SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C
    M_s = np.array([[rA * rA,        -rA * re_ * k1],
                    [-rA * re_ * k1,  re_ * re_ * k2]], dtype=np.float64)
    M_c = np.array([[sA * sA,         sA * sC * k2],
                    [sA * sC * k2,    sC * sC * k4]], dtype=np.float64)

    def _eb(M):
        lam, E = np.linalg.eigh(M)        # ascending
        E, lam = E[:, ::-1], lam[::-1]    # → (stiff, sloppy)
        lam = lam / max(lam[1], 1e-300)   # normalise: sloppy = 1
        for j in range(2):                # sign: dominant component positive
            if E[np.argmax(np.abs(E[:, j])), j] < 0:
                E[:, j] = -E[:, j]
        return E, lam

    Es, ls = _eb(M_s)
    Ec, lc = _eb(M_c)
    return Es, ls, Ec, lc


def _project_whitened(phys, ref, E):
    """Physical (A,e)/(a,c) → O(1) (÷ref) → project onto eigenvectors E.

    ``phys`` last axis is size 2 (any leading shape); ``ref`` is the (2,) O(1)
    reference scale; ``E`` is [2, 2] with columns (stiff, sloppy). Returns the
    same shape with the last axis = (stiff, sloppy) dimensionless coordinates."""
    o1 = np.asarray(phys, dtype=np.float64) / np.asarray(ref, dtype=np.float64)
    return o1 @ E


class _override_theta:
    """Context manager that temporarily sets the model's θ_scale / θ_smear to
    the supplied raw values and restores on exit. Used to evaluate the model
    density at the INJECTED θ values for the validation-closure plots —
    showing where the fitted curve *should* converge to."""

    def __init__(self, model, raw_scale=None, raw_smear=None):
        self.model = model
        self.raw_scale = raw_scale
        self.raw_smear = raw_smear
        self._saved_scale = None
        self._saved_smear = None

    def __enter__(self):
        if self.raw_scale is not None:
            self._saved_scale = self.model.theta_scale.detach().clone()
            with torch.no_grad():
                self.model.theta_scale.copy_(
                    self.raw_scale.to(self.model.theta_scale.device,
                                      self.model.theta_scale.dtype))
        if self.raw_smear is not None:
            self._saved_smear = self.model.theta_smear.detach().clone()
            with torch.no_grad():
                self.model.theta_smear.copy_(
                    self.raw_smear.to(self.model.theta_smear.device,
                                      self.model.theta_smear.dtype))
        return self

    def __exit__(self, *_):
        if self._saved_scale is not None:
            with torch.no_grad():
                self.model.theta_scale.copy_(self._saved_scale)
        if self._saved_smear is not None:
            with torch.no_grad():
                self.model.theta_smear.copy_(self._saved_smear)


def _continuity_mc_fold(model, ptm, etam, phim, qm, bm):
    """Fold MC reco at the fitted θ — the empirical template the model signal
    curve should reproduce. The per-muon PHYSICAL fold, in the model's
    generative order: FIRST the forward Gaussian σ_qop smear at the truth pt
    (``model.fold_sigma_qop_pm``, the signed (a, c) clipped at 0), THEN the
    forward scale — the exact functional inverse of the model's defining
    BACKWARD (data → MC) map, via ``model._scale_apply_pt_forward`` (fixed
    point; the iteration cost lives on the fold side, where it is free). m_ll
    is RECOMPUTED from the folded 4-vectors. This is the exact operation the
    continuity density inverts, so the closure overlay exposes the residual
    approximations. (ρ is a flow condition, not histogrammed/binned here, so
    it is not recomputed — no effect on the m_ll closure plots.)"""
    pt_cur = ptm
    if model.smearing_enabled:
        sig = model.fold_sigma_qop_pm(pt_cur, etam, phim, bm)
        pt_cur = model.apply_smear_pt(pt_cur, etam, qm, sig, torch.randn_like(sig))
    if model.scale_enabled:
        pt_cur = model._scale_apply_pt_forward(pt_cur, etam, phim, qm, bm)
    return _event_mll(pt_cur, etam, phim).detach()


N_PHI_BKG_BINS = 16


def _accum_bkg_frac(acc, f, w, b_pm, phi_pm, label, slope=None):
    """Accumulate the background-fraction closure sums per η bin (via the
    loader's b_pm index) and per φ bin (``N_PHI_BKG_BINS`` uniform over
    [−π, π]): columns = (Σw, Σw·f_k for each background component,
    Σw·f_bkg_total, Σw·[label==1], Σw·[label==2]), each event filling BOTH
    muon legs (the MLP conditioning carries both). Component count follows
    the active background model (signal is the LAST mixture output)."""
    w_np = w.detach().cpu().numpy().astype(np.float64)
    f_np = f.detach().cpu().numpy().astype(np.float64)
    lab_np = label.detach().cpu().numpy()
    cols = ([np.ones_like(w_np)]
            + [f_np[:, k] for k in range(f_np.shape[1] - 1)]
            + [1.0 - f_np[:, -1],
               (lab_np == 1).astype(np.float64),
               (lab_np == 2).astype(np.float64)])
    if slope is not None:   # exp model: Σw·s(c) for the per-bin mean slope
        cols.append(slope.detach().cpu().numpy().astype(np.float64))
    vals = np.stack(cols, axis=1) * w_np[:, None]
    b_np = b_pm.detach().cpu().numpy()
    phi_np = phi_pm.detach().cpu().numpy()
    nb_eta = acc["eta"].shape[0]
    nb_phi = acc["phi"].shape[0]
    pbin = np.clip(((phi_np + np.pi) / (2.0 * np.pi) * nb_phi).astype(np.int64),
                   0, nb_phi - 1)
    for leg in (0, 1):
        bl = np.clip(b_np[:, leg], 0, nb_eta - 1)
        for k in range(vals.shape[1]):
            acc["eta"][:, k] += np.bincount(bl, weights=vals[:, k],
                                            minlength=nb_eta)[:nb_eta]
            acc["phi"][:, k] += np.bincount(pbin[:, leg], weights=vals[:, k],
                                            minlength=nb_phi)[:nb_phi]


@torch.no_grad()
def evaluate_predictions(
    model: JpsiMassMixtureModel,
    loader: JpsiMassArrowLoader,
    device: str,
    m_centers: torch.Tensor,
    m_grid_std: torch.Tensor,
    bin_width: float,
    *,
    chunk_events: int = 4096,
    max_events: int = 0,
    progress: bool = True,
    seed: int = 42,
    n_iter: int = 2,
    mc_as_data: bool = False,
):
    """Stream the loader once; collect per-event aggregates AND per-event
    × per-bin signal densities for the model curves.

    The signal grid density is the **#2 direct-eval** ``_tilt_density_on_grid``
    (the same forward-folded flow density the fit optimises) and the MC
    comparison is the **per-muon physical fold** (``_continuity_mc_fold``:
    scale δqop + Gaussian σ_qop smear, m_ll recomputed).

    Per data event: store m_ll, w, |η_+|, MLP outputs ``f = (f_0, f_1, f_s)``,
    and ``pred_signal_data[e, j] = p_flow(m_bin_j | y_e, fitted θ_scale+θ_smear)``.
    Per MC event: store m_fold (forward-fold of scale+smear at the fitted
    nuisances), w, |η_+|, and
    ``pred_signal_mc[e, j] = p_flow(m_bin_j | y_e, fitted θ_scale+θ_smear)``.

    Cost: O((N_data + N_mc) × n_grid) flow forwards — bounded by
    ``--grid-chunk-events`` per call and ``--max-events`` per pass.

    Returns dict (1-D unless noted):
      mll_data, w_data, eta_data           — per data event
      f_data                                — [N_data, 3]
      pred_signal_data                      — [N_data, n_grid]
      mll_mc_fold, w_mc, eta_mc            — per MC event (scale+smear fold)
      pred_signal_mc                        — [N_mc, n_grid]
      bin_width                              — scalar passed-through
    """
    # Save + seed RNG so the smearing-kernel ε is deterministic across
    # diagnostic runs.
    cpu_state = torch.random.get_rng_state()
    cuda_state = None
    cuda_avail = device.startswith("cuda") and torch.cuda.is_available()
    if cuda_avail:
        try:
            cuda_state = torch.cuda.get_rng_state(device)
        except Exception:
            cuda_state = None
    torch.manual_seed(seed)
    if cuda_avail:
        torch.cuda.manual_seed_all(seed)

    m_grid_std_dev = m_grid_std.to(device)
    m_centers_dev = m_centers.to(device)
    n_grid = m_centers.shape[0]
    print("  (continuity model: signal curve = #2 direct-eval tilt; "
          "MC = per-muon scale+smear fold with m_ll recomputed)")
    if mc_as_data:
        print("  (validation: simulation routed through the data branch as "
              "pseudo-data for the m_ll closure / pulls, in addition to the "
              "MC-branch closure)")

    def _sig_grid(idx):
        return _tilt_density_on_grid(
            model, batch, idx, m_centers_dev, chunk_events=chunk_events,
            n_iter=n_iter)

    # Conditional slice variables for the closure plots, chosen for the active
    # basis (muon_kin → ρ, cos α; event_level → p_T^ll, y_ll, cos θ*). The
    # tertile-mode ones are accumulated per event below; |η₊| (eta_edges mode)
    # uses eta_data.
    _slv = _diag_slice_vars(getattr(model, "cond_basis", "muon_kin"))
    _tkeys = [k for k, _, _, m in _slv if m == "tertile"]

    out = {
        "mll_data": [], "w_data": [], "eta_data": [], "f_data": [],
        "bkg_slope_data": [],
        "pred_signal_data": [],
        "mll_mc_fold": [], "w_mc": [], "eta_mc": [],
        "pred_signal_mc": [],
        # continuity only: the nominal (θ=0, unshifted/unsmeared) MC + flow,
        # to overlay the stage-1 closure alongside the folded stage-2 one.
        "mll_mc_nominal": [], "pred_nominal_mc": [],
        # validation-with-injection only: the INJECTED pseudo-data m_ll (the
        # closure target the fold should reproduce). Distinct from the nominal
        # whenever a θ/smear injection was replayed into the loader.
        "mll_mc_pseudodata": [],
        # validation-with-injection only: signal density on the grid evaluated
        # at the INJECTED θ values (the closure target curve — where the
        # fitted-θ density should converge to if the fit recovers the truth).
        "pred_signal_mc_at_inj": [],
    }
    for _k in _tkeys:
        out[f"sl_{_k}_data"] = []
        out[f"sl_{_k}_mc"] = []
    # Convert the loader's replayed injection (physical) to RAW θ tensors that
    # produce them through the model — used below to evaluate the model
    # density at the truth.
    inject_scale_phys = getattr(loader, "inject_theta_scale", None)
    inject_smear_phys = getattr(loader, "inject_theta_smear", None)
    has_inj = inject_scale_phys is not None or inject_smear_phys is not None
    if has_inj:
        # Full per-η-bin PHYSICAL injection (0 for non-injected terms — the truth
        # there is θ=0). Used to SET the model's effective per-muon θ for the
        # 'flow at injected θ' curve (works for binned AND MLP).
        n_eta = model.theta_scale.shape[0]
        inj_scale_full = (np.asarray(inject_scale_phys, dtype=np.float64)
                          if inject_scale_phys is not None else np.zeros((n_eta, 3)))
        inj_smear_full = (np.asarray(inject_smear_phys, dtype=np.float64)
                          if inject_smear_phys is not None else np.zeros((n_eta, 2)))

    # Background-fraction closure sums (η via b_pm; φ via N_PHI_BKG_BINS).
    # Columns: Σw, Σw·f_k (one per background component), Σw·f_bkg_total,
    # Σw·[label==1], Σw·[label==2].
    _nbk = int(getattr(model, "n_bkg_comp", 2))
    _nac = _nbk + 4 + (1 if getattr(model, "bkg_model", "bernstein") == "exp"
                       else 0)   # exp: + Σw·s(c) column (LAST)
    # η accumulator is indexed by b_pm (the per-muon η-bin), so it must be sized
    # by the true η-bin count — NOT theta_scale.shape[0], which is n_eta×n_phi
    # cells in the 2-D binned mode.
    bkg_acc = {"eta": np.zeros((getattr(model, "n_eta_bins",
                                        model.theta_scale.shape[0]), _nac)),
               "phi": np.zeros((N_PHI_BKG_BINS, _nac))}

    total_events = 0
    bar = tqdm(loader, desc="eval", disable=not progress, unit="batch")
    try:
        for batch in bar:
            if max_events > 0 and total_events >= max_events:
                break
            batch = _move_batch(batch, device)
            is_data = batch["is_data_mask"]
            is_mc = ~is_data
            # Validation (MC-closure) checkpoints fit θ on simulation as
            # pseudo-data, so route MC through the data branch too. The MC
            # branch always uses the simulation rows, so in that mode the same
            # events feed both closures.
            data_sel = is_mc if mc_as_data else is_data
            mc_sel = is_mc

            if bool(data_sel.any()):
                data_idx = data_sel.nonzero(as_tuple=True)[0]
                _slope_d = None
                if getattr(model, "background_enabled", True):
                    if getattr(model, "bkg_model", "bernstein") == "exp":
                        f, _slope_d = model.mlp.forward_with_slope(
                            batch["cond_std"][data_idx])
                    else:
                        f = model.f_data(batch["cond_std"][data_idx])
                else:
                    # --no-background: data branch is pure signal. The MLP is
                    # bypassed in data_nll_continuity; mirror that here so the
                    # mixture plot doesn't show whatever the random-init MLP
                    # happens to output. f ≡ [0, 0, 1] → no Bernstein bkg.
                    f = torch.zeros((data_idx.numel(),
                                     int(getattr(model, "n_bkg_comp", 2)) + 1),
                                    device=batch["cond_std"].device,
                                    dtype=batch["cond_std"].dtype)
                    f[:, -1] = 1.0
                if getattr(model, "background_enabled", True):
                    _lab = batch.get("bkg_label")
                    _lab = (_lab[data_idx] if _lab is not None else
                            torch.zeros(data_idx.numel(), dtype=torch.int8,
                                        device=data_idx.device))
                    _accum_bkg_frac(bkg_acc, f, batch["w"][data_idx],
                                    batch["b_pm"][data_idx],
                                    batch["phi_pm"][data_idx], _lab,
                                    slope=_slope_d)
                # Signal density at every bin centre for every data event at the
                # fitted θ: tilt (continuity) or θ-conditioned flow (legacy).
                log_p_grid = _sig_grid(data_idx)  # [n_data, n_grid] log-density (1/GeV)
                out["mll_data"].append(batch["mll"][data_idx].cpu().numpy())
                out["w_data"].append(batch["w"][data_idx].cpu().numpy())
                out["eta_data"].append(batch["eta_pm"][data_idx, 0].cpu().numpy())
                _sv_d = _slice_vals_np(
                    _tkeys, batch["pt_pm"][data_idx], batch["eta_pm"][data_idx],
                    batch["phi_pm"][data_idx])
                for _k in _tkeys:
                    out[f"sl_{_k}_data"].append(_sv_d[_k])
                out["f_data"].append(f.cpu().numpy())
                if _slope_d is not None:
                    out["bkg_slope_data"].append(_slope_d.detach().cpu().numpy())
                out["pred_signal_data"].append(log_p_grid.exp().cpu().numpy())

            if bool(mc_sel.any()):
                mc_idx = mc_sel.nonzero(as_tuple=True)[0]
                # Directly shift+smear the MC at the *fitted* θ — the empirical
                # template the model signal curve should reproduce. Use the
                # NOMINAL (pre-injection) pt: in validation the loader replaces
                # pt_pm with the injected/smeared pt (so mll = _event_mll(pt_pm)
                # is consistent), so folding pt_pm again would DOUBLE the
                # injection. pt_pm_nominal is the un-injected pt (== pt_pm when
                # no injection / older loaders without the field).
                ptm = batch.get("pt_pm_nominal", batch["pt_pm"])[mc_idx]
                etam = batch["eta_pm"][mc_idx]
                phim = batch["phi_pm"][mc_idx]
                qm = batch["q_pm"][mc_idx]
                bm = batch["b_pm"][mc_idx]
                mll_fold = _continuity_mc_fold(model, ptm, etam, phim, qm, bm)
                # Signal density on the grid for every MC event (#2 tilt).
                log_p_grid_mc = _sig_grid(mc_idx)  # [n_mc, n_grid]
                out["mll_mc_fold"].append(mll_fold.cpu().numpy())
                out["w_mc"].append(batch["w"][mc_idx].cpu().numpy())
                out["eta_mc"].append(batch["eta_pm"][mc_idx, 0].cpu().numpy())
                # Slice variables from the conditioning the model sees (batch
                # pt_pm = injected pt in validation), matching the data side so a
                # given bin compares the same physical region.
                _sv_m = _slice_vals_np(
                    _tkeys, batch["pt_pm"][mc_idx], etam, phim)
                for _k in _tkeys:
                    out[f"sl_{_k}_mc"].append(_sv_m[_k])
                out["pred_signal_mc"].append(log_p_grid_mc.exp().cpu().numpy())
                # Closure target: signal density at the INJECTED θ values.
                # Temporarily SET the model's effective per-muon θ to the
                # injection, evaluate the tilt density, restore. Uses
                # _set_theta_output so it works for the MLP too (the old
                # _override_theta set only the binned tensors → no-op for
                # --theta-mlp, leaving this curve meaningless). The fitted curve
                # should converge to this when the closure is good.
                if has_inj:
                    with _set_theta_output(
                            model, inj_scale_full, inj_smear_full,
                            nonuniform=getattr(loader, "inject_nonuniform", False)):
                        log_p_grid_mc_inj = _sig_grid(mc_idx)
                    out["pred_signal_mc_at_inj"].append(
                        log_p_grid_mc_inj.exp().cpu().numpy())
                # nominal (θ=0): the TRUE un-injected reco mass, recomputed from
                # the (un-injected) per-muon pt — NOT batch["mll"], which carries
                # the replayed validation injection. Plus the untilted flow p₀.
                out["mll_mc_nominal"].append(
                    _event_mll(ptm, etam, phim).cpu().numpy())
                # The injected pseudo-data m_ll (= nominal + replayed injection),
                # i.e. the closure target the fold should reproduce. Equals the
                # nominal when no injection was replayed.
                out["mll_mc_pseudodata"].append(batch["mll"][mc_idx].cpu().numpy())
                out["pred_nominal_mc"].append(
                    _nominal_density_on_grid(
                        model, batch, mc_idx, m_centers_dev,
                        chunk_events=chunk_events).exp().cpu().numpy())

            total_events += int(batch["mll"].shape[0])
            bar.set_postfix_str(f"n_events={total_events:,}")
    finally:
        bar.close()
        torch.random.set_rng_state(cpu_state)
        if cuda_state is not None:
            torch.cuda.set_rng_state(cuda_state, device)

    for k, lst in out.items():
        if lst:
            out[k] = np.concatenate(lst, axis=0)
        else:
            if k == "f_data":
                out[k] = np.zeros((0, 3))
            elif k in ("pred_signal_data", "pred_signal_mc", "pred_nominal_mc",
                       "pred_signal_mc_at_inj"):
                out[k] = np.zeros((0, n_grid))
            else:
                out[k] = np.zeros((0,))
    out["bin_width"] = bin_width
    out["bkg_frac"] = bkg_acc
    # Per-component ∫_bin p_k dm for the ACTIVE background model on the
    # plotting bins (consumed by _model_pred_histograms).
    _edges_np = np.concatenate([m_centers.cpu().numpy() - 0.5 * bin_width,
                                [float(m_centers[-1]) + 0.5 * bin_width]])
    out["bkg_bin_integrals"] = _bkg_bin_integrals_model(model, _edges_np)
    out["continuity"] = True
    out["mc_as_data"] = mc_as_data
    # Basis-aware slice variables (key, label, fmt, mode) consumed by the closure
    # plots; tertile-mode values are in out["sl_<key>_data/_mc"].
    out["slice_specs"] = _slv
    # Whether the loader replayed a validation injection into batch["mll"] — so
    # the closure plot draws the injected pseudo-data curve only when it is
    # genuinely distinct from the (un-injected) nominal.
    out["injected"] = bool(
        getattr(loader, "inject_theta_scale", None) is not None
        or getattr(loader, "inject_theta_smear", None) is not None)
    return out


def _bernstein_bin_integrals(m_lo: float, m_hi: float, m_edges: np.ndarray):
    """∫_bin p_0 dm and ∫_bin p_1 dm, closed form.

    p_0(m) = 2(1 − u)/width, p_1(m) = 2u/width, u = (m − m_lo)/width
    so ∫(1 − u) du = u − u²/2 and ∫u du = u²/2 ⇒ each per-bin integral
    is 2·ΔF(u). Returns ``(I_0_per_bin, I_1_per_bin)`` each shape
    ``[len(m_edges) − 1]``.
    """
    width = m_hi - m_lo
    u = (m_edges - m_lo) / width
    F0 = lambda x: x - 0.5 * x * x  # noqa: E731
    F1 = lambda x: 0.5 * x * x      # noqa: E731
    return 2.0 * (F0(u[1:]) - F0(u[:-1])), 2.0 * (F1(u[1:]) - F1(u[:-1]))


def _bkg_bin_integrals_model(model, m_edges: np.ndarray) -> np.ndarray:
    """Per-component ``∫_bin p_k(m) dm`` for the BERNSTEIN background by
    composite trapezoid on a 32-point subgrid per bin (plot-precision exact).
    Returns ``[n_bkg_comp, n_bins]``. The exp model's per-EVENT slope makes
    its bin contents event-dependent — handled by the closed-form CDF in
    ``_model_pred_histograms`` instead (returns zeros here)."""
    n_bins = len(m_edges) - 1
    n_bkg = int(getattr(model, "n_bkg_comp", 2))
    out = np.zeros((n_bkg, n_bins))
    if (not getattr(model, "background_enabled", True)
            or getattr(model, "bkg_model", "bernstein") == "exp"):
        return out
    n_sub = 32
    for j in range(n_bins):
        g = torch.linspace(float(m_edges[j]), float(m_edges[j + 1]), n_sub + 1,
                           dtype=torch.float64)
        with torch.no_grad():
            pk = bernstein_basis_n(
                g, model._m_lo_f, model._m_hi_f,
                int(getattr(model, "bkg_degree", 1))).T
        # np.trapz was removed in numpy 2 (renamed np.trapezoid); torch's
        # trapezoid exists in every torch we support — use it directly.
        out[:, j] = torch.trapezoid(pk, g, dim=1).numpy()
    return out


def _exp_bin_fractions_np(u_edges: np.ndarray, slopes: np.ndarray) -> np.ndarray:
    """Closed-form per-event exp-background bin fractions
    ``F(u_{j+1}; s_e) − F(u_j; s_e)`` with ``F(u; s) = (1−e^{−s·u})/(1−e^{−s})``
    (window CDF), stable through s → 0 via ``F ≈ u·(1 + s(1−u)/2)``.
    ``u_edges [n_bins+1]``, ``slopes [n_e]`` → ``[n_e, n_bins]``."""
    s = slopes[:, None]
    u = u_edges[None, :]
    small = np.abs(s) < 1e-4
    s_safe = np.where(small, 1.0, s)
    F_exact = -np.expm1(-s_safe * u) / (-np.expm1(-s_safe))
    F_small = u * (1.0 + s * (1.0 - u) / 2.0)
    F = np.where(small, F_small, F_exact)
    return np.diff(F, axis=1)


def _model_pred_histograms(
    evals, m_edges: np.ndarray, m_lo: float, m_hi: float,
    slice_mask_data,
):
    """``(signal, bkg_total)`` per-bin predicted counts for one slice.

    signal[j] = Δm · Σ_data w_e · f_s(y_e) · p_flow(m_bin_j | y_e, σ_e)
              (grid eval over data events — flow density at each bin centre)

    bkg[j] = Σ_k (Σ_data w_e · f_k(y_e)) · ∫_bin_j p_k(m) dm
              (per-component analytic/numeric integrals × MLP fractions;
              components follow the active background model — signal is the
              LAST mixture output)
    """
    n_bins = len(m_edges) - 1
    bin_width = float(m_edges[1] - m_edges[0])
    sig = np.zeros(n_bins)
    bkg = np.zeros(n_bins)

    if evals["mll_data"].size and slice_mask_data.any():
        w_d = evals["w_data"][slice_mask_data]
        f_d = evals["f_data"][slice_mask_data]
        pred = evals["pred_signal_data"][slice_mask_data]  # [n_d, n_bins]
        weights = (w_d * f_d[:, -1])[:, None]  # [n_d, 1]
        sig = bin_width * (pred * weights).sum(axis=0)
        slopes = np.asarray(evals.get("bkg_slope_data", np.zeros(0)))
        if slopes.size == evals["mll_data"].size and slopes.size > 0:
            # exp model (per-event slope array parallel to the data events):
            # closed-form per-event bin fractions × fraction weights.
            u_edges = (np.asarray(m_edges) - m_lo) / (m_hi - m_lo)
            dF = _exp_bin_fractions_np(u_edges, slopes[slice_mask_data])
            bkg = ((w_d * f_d[:, 0])[:, None] * dF).sum(axis=0)
        else:
            sum_fk_w = (w_d[:, None] * f_d[:, :-1]).sum(axis=0)   # [n_bkg]
            I = evals.get("bkg_bin_integrals")
            if I is not None and I.shape[0] == sum_fk_w.shape[0]:
                bkg = (sum_fk_w[:, None] * I).sum(axis=0)
            else:  # legacy fallback: degree-1 closed form
                I0, I1 = _bernstein_bin_integrals(m_lo, m_hi, m_edges)
                bkg = sum_fk_w[0] * I0 + (sum_fk_w[1] * I1 if len(sum_fk_w) > 1 else 0.0)

    return sig, bkg


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------


def _save_fig(fig, output_dir: str, stem: str, formats=("png", "pdf"), dpi: int = 110):
    """Write ``fig`` as ``<stem>.<ext>`` for each requested format.
    Returns the list of written paths (for the diagnostic stdout line).
    """
    out_paths = []
    for ext in formats:
        path = os.path.join(output_dir, f"{stem}.{ext}")
        fig.savefig(path, dpi=dpi)
        out_paths.append(path)
    plt.close(fig)
    return out_paths


def _chi2_compat_zero(theta: np.ndarray, cov: np.ndarray):
    """χ² for the compatibility of ``theta`` with zero given covariance ``cov``,
    via the eigendecomposition restricted to the POSITIVE-variance directions:
    ``χ² = Σ_{λ_k>0} (v_kᵀθ)² / λ_k``, ``dof = #{λ_k > 0}``. This is a proper
    (≥0) generalised χ² that drops both zero-variance (rank-deficient / pinv)
    and any negative-eigenvalue directions (a non-PD Hessian covariance), so it
    is well-defined for every estimator. Returns ``(chi2, dof, p_value)``."""
    theta = np.asarray(theta, dtype=np.float64).reshape(-1)
    cov = np.asarray(cov, dtype=np.float64)
    cov = 0.5 * (cov + cov.T)
    w, V = np.linalg.eigh(cov)
    wmax = float(w.max()) if w.size else 0.0
    tol = cov.shape[0] * np.finfo(np.float64).eps * max(wmax, 0.0)
    pos = w > tol
    dof = int(pos.sum())
    chi2 = 0.0
    if dof > 0:
        proj = V[:, pos].T @ theta            # θ along the positive-variance dirs
        chi2 = float(np.sum(proj * proj / w[pos]))
    p = float("nan")
    if dof > 0 and np.isfinite(chi2):
        try:
            from scipy.stats import chi2 as _chi2dist
            p = float(_chi2dist.sf(chi2, dof))
        except Exception:
            try:
                from scipy.special import gammaincc
                p = float(gammaincc(dof / 2.0, chi2 / 2.0))
            except Exception:
                p = float("nan")
    return chi2, dof, p


def _diag_slice_vars(cond_basis):
    """Slice variables for the closure / param-sensitivity plots, all DERIVABLE
    FROM THE ACTIVE CONDITIONING (so a slice never leaks the observable m_ll).
    Each entry: ``(key, label, fmt, mode)``; mode ``'eta_edges'`` uses the passed
    fixed η-bin edges, ``'tertile'`` uses adaptive tertiles.

    - ``muon_kin``: |η₊| (directions are in the basis), ρ (= muon_kin[-1]), and
      cos α (3-D opening angle, a pure function of the muon directions).
    - ``event_level``: the conditioning components themselves — p_T^ll/m_ll, y_ll, and
      cos θ* (CS polar decay angle). cos α / |η₊| / ρ are NOT used here: they
      depend on the lab→CS boost (hence m_ll), so they are not functions of the
      event-level conditioning alone. cos θ* is the event-level analog of ρ for
      M-sensitivity; p_T^ll/m_ll (the basis' dimensionless dilepton-pt ratio)
      discriminates the A/e and a/c pt-scale degeneracies."""
    if cond_basis == "event_level":
        return [("ptll_over_mll", "p_T^ll/m_ll", ".2f", "tertile"),
                ("yll", "y_ll", ".2f", "tertile"),
                ("costhetastar", "cos θ*", ".2f", "tertile")]
    return [("eta", "|η₊|", ".1f", "eta_edges"),
            ("rho", "ρ", ".2f", "tertile"),
            ("cosalpha", "cos α", ".3f", "tertile")]


def _slice_vals_np(keys, pt_pm, eta_pm, phi_pm):
    """``{key: [N] array}`` for the requested tertile slice keys, from per-muon
    pt/η/φ. ptll/mll, yll, cosθ* come from ``_event_cond_raw_np`` (cols 1/0/4
    — col 1 is u = ln(ptll/mll), so exp(u) is the dimensionless ratio); ρ and
    cos α are the pt-asymmetry and opening angle. cos α uses only directions
    (pt cancels): ``(cosΔφ + sinhη₊sinhη₋)/(coshη₊coshη₋)``."""
    keys = set(keys)
    pt = pt_pm.detach().cpu().numpy()
    e = eta_pm.detach().cpu().numpy()
    p = phi_pm.detach().cpu().numpy()
    out = {}
    if "rho" in keys:
        out["rho"] = (pt[:, 0] - pt[:, 1]) / (pt[:, 0] + pt[:, 1])
    if "cosalpha" in keys:
        num = np.cos(p[:, 0] - p[:, 1]) + np.sinh(e[:, 0]) * np.sinh(e[:, 1])
        out["cosalpha"] = num / (np.cosh(e[:, 0]) * np.cosh(e[:, 1]))
    if keys & {"ptll_over_mll", "yll", "costhetastar"}:
        ev = _event_cond_raw_np(pt.astype(np.float32), e.astype(np.float32),
                                p.astype(np.float32))
        if "ptll_over_mll" in keys:
            out["ptll_over_mll"] = np.exp(ev[:, 1].astype(np.float64))
        if "yll" in keys:
            out["yll"] = ev[:, 0]
        if "costhetastar" in keys:
            out["costhetastar"] = ev[:, 4]
    return out


def _select_slice(eta_abs: np.ndarray, slice_def):
    """Boolean mask for a |η| slice (lo, hi) or None for inclusive."""
    if slice_def is None:
        return np.ones_like(eta_abs, dtype=bool)
    lo, hi = slice_def
    return (eta_abs >= lo) & (eta_abs < hi)


def _closure_slice_dims(evals, eta_slice_edges):
    """Slice dimensions shared by the m_ll- and MC-closure PANEL plots, from the
    active basis's ``evals["slice_specs"]`` (all DERIVABLE FROM THE CONDITIONING):
      • muon_kin → |η₊| (fixed edges + leading inclusive panel), ρ, cos α tertiles.
      • event_level → p_T^ll, y_ll, cos θ* tertiles (the conditioning components;
        cos α / |η₊| / ρ are not functions of the event-level conditioning alone).
    The inclusive panel leads the FIRST dimension. Returns a list of
    ``(prefix, label, fmt, data_vals, mc_vals, columns)`` with ``columns`` the
    per-panel ``[(tag, slice_def), ...]``."""
    def _tertiles(d, mc):
        v = d if (d is not None and d.size) else mc
        if v is None or v.size == 0:
            return None
        v = v[np.isfinite(v)]
        if v.size == 0:
            return None
        e = np.unique(np.percentile(v, [0.0, 100.0 / 3.0, 200.0 / 3.0, 100.0]))
        return e if e.size >= 2 else None

    specs = evals.get("slice_specs") or _diag_slice_vars("muon_kin")
    dims = []
    for key, label, fmt, mode in specs:
        if mode == "eta_edges":
            dv = np.abs(evals["eta_data"]); mv = np.abs(evals["eta_mc"])
            cols = [(f"{key}{i}", (eta_slice_edges[i], eta_slice_edges[i + 1]))
                    for i in range(len(eta_slice_edges) - 1)]
        else:
            dv = evals.get(f"sl_{key}_data", np.zeros((0,)))
            mv = evals.get(f"sl_{key}_mc", np.zeros((0,)))
            edges = _tertiles(dv, mv)
            if edges is None:
                continue
            cols = [(f"{key}{i}", (edges[i], edges[i + 1]))
                    for i in range(len(edges) - 1)]
        dims.append((key, label, fmt, dv, mv, cols))
    # The shared inclusive panel leads the first dimension's figure.
    if dims:
        k, lb, fm, dv, mv, cols = dims[0]
        dims[0] = (k, lb, fm, dv, mv, [("inclusive", None)] + cols)
    return dims


def plot_mll_closure(
    evals, m_centers_np, eta_slice_edges, m_lo: float, m_hi: float,
    output_dir: str,
):
    """Plot data + forward-folded-MC histograms with overlaid model curves.
    One PANEL FIGURE per slice dimension (columns = slices): inclusive + |η_+|
    (fixed edges), ρ (pt asymmetry) tertiles, and cos α (3-D opening angle)
    tertiles — the conditional discriminants of plot_param_sensitivity (ρ → M;
    cos α → A vs e and a vs c). The model signal curve is the flow density
    evaluated on a per-event × per-bin grid at the fitted θ_scale + θ_smear; the
    green MC curve is the MC forward-folded (scale+smear) at the same fitted
    nuisances — independent estimates of the signal."""
    m_edges = np.concatenate([
        [m_centers_np[0] - evals["bin_width"] / 2],
        m_centers_np[:-1] + evals["bin_width"] / 2,
        [m_centers_np[-1] + evals["bin_width"] / 2],
    ])
    # Validation runs fit θ on simulation as pseudo-data — label accordingly.
    pseudo = bool(evals.get("mc_as_data", False))
    data_label = "MC (pseudo-data)" if pseudo else "data"
    cont = bool(evals.get("continuity", False))

    def _draw_panel(ax, axr, data_mask, mc_mask):
        """Draw one slice into (main, ratio) axes and set the ratio panel's
        own autoscaled y-range; return (ymax, has_curves)."""
        # Data hist.
        if data_mask.any():
            data_hist, _ = np.histogram(
                evals["mll_data"][data_mask], bins=m_edges,
                weights=evals["w_data"][data_mask])
            ax.errorbar(m_centers_np, data_hist, yerr=np.sqrt(np.abs(data_hist)),
                        fmt="o", color="k", markersize=3, label=data_label, zorder=3)
        else:
            data_hist = np.zeros(m_centers_np.shape[0])
        # Forward-folded MC, scaled to the data signal weight in the slice.
        if mc_mask.any():
            mc_hist_raw, _ = np.histogram(
                evals["mll_mc_fold"][mc_mask], bins=m_edges,
                weights=evals["w_mc"][mc_mask])
            w_mc_sum = float(evals["w_mc"][mc_mask].sum())
            scale = 0.0
            if data_mask.any() and w_mc_sum > 0:
                w_d = evals["w_data"][data_mask]
                f_d = evals["f_data"][data_mask]
                # signal fraction is the LAST mixture output (generalises
                # the historical f[:, 2] of the 3-way deg-1 Bernstein head)
                scale = float((w_d * f_d[:, -1]).sum()) / w_mc_sum
            mc_hist = mc_hist_raw * scale
            mc_label = ("MC (shifted+smeared, scaled to signal wt)" if cont
                        else "MC (scale+smear folded, scaled to signal wt)")
            ax.step(m_edges[:-1], mc_hist, where="post", color="C2", lw=1.2,
                    label=mc_label)
        else:
            mc_hist = np.zeros(m_centers_np.shape[0])
        # Model components — analytic bkg + grid-eval signal (tilt / flow).
        signal, bkg_tot = _model_pred_histograms(
            evals, m_edges, m_lo, m_hi, data_mask)
        total = signal + bkg_tot
        sig_label = ("model signal (tilt p₀·e^δ at fitted θ)" if cont
                     else "model signal (flow at fitted scale+smear)")
        if total.sum() > 0:
            ax.plot(m_centers_np, total, color="C0", lw=1.5, label="model total")
            ax.plot(m_centers_np, signal, color="C1", lw=1.0, ls="--",
                    label=sig_label)
            ax.fill_between(m_centers_np, 0, bkg_tot, alpha=0.3, color="C3",
                            label="model bkg (Bernstein)")
            with np.errstate(divide="ignore", invalid="ignore"):
                denom = np.where(total > 0, total, np.nan)
                ratio = data_hist / denom
                ratio_err = np.sqrt(np.abs(data_hist)) / denom
            axr.errorbar(m_centers_np, ratio, yerr=ratio_err, fmt="o",
                         color="k", markersize=3)
            # ratio-panel autoscale: a ROBUST symmetric half-window about 1.
            # Use only well-populated bins so the large-error tail bins don't
            # blow up the zoom; size the window from the 90th-percentile
            # deviation (or the typical per-bin error, whichever is larger).
            good = np.isfinite(ratio) & (
                data_hist > 0.10 * max(float(data_hist.max()), 1e-30))
            if good.any():
                dev = float(np.percentile(np.abs(ratio[good] - 1.0), 90))
                half = max(dev, float(np.median(ratio_err[good]))) * 1.3
            else:
                half = None
        else:
            half = None
        # Each slice gets its OWN symmetric ratio window (clamped so it is
        # neither absurdly tight nor wider than the legacy [0.6, 1.4]).
        if half is not None:
            h = min(max(half, 0.03), 0.4)
            axr.set_ylim(1.0 - h, 1.0 + h)
        else:
            axr.set_ylim(0.6, 1.4)
        axr.axhline(1.0, color="C0", lw=1)
        ymax = max(float((data_hist + np.sqrt(np.abs(data_hist))).max()),
                   float(mc_hist.max()), float(total.max()))
        return ymax, total.sum() > 0

    for prefix, label, fmt, dv, mv, cols in _closure_slice_dims(
            evals, eta_slice_edges):
        ncol = len(cols)
        fig, axes = plt.subplots(
            2, ncol, figsize=(max(6.0, 4.3 * ncol), 6.2), squeeze=False,
            sharex="col", gridspec_kw={"height_ratios": [3, 1]})
        leg_ax = None
        for ci, (tag, slice_def) in enumerate(cols):
            ax, axr = axes[0, ci], axes[1, ci]
            data_mask = _select_slice(dv, slice_def) if dv.size else np.zeros((0,), bool)
            mc_mask = _select_slice(mv, slice_def) if mv.size else np.zeros((0,), bool)
            if data_mask.sum() == 0 and mc_mask.sum() == 0:
                ax.set_visible(False); axr.set_visible(False); continue
            ymax, ok = _draw_panel(ax, axr, data_mask, mc_mask)
            if ok and leg_ax is None:
                leg_ax = ax
            ttl = ("inclusive" if slice_def is None
                   else f"{label} ∈ [{slice_def[0]:{fmt}}, {slice_def[1]:{fmt}}]")
            ax.set_title(ttl, fontsize=9)
            if ymax > 0:
                ax.set_ylim(0, ymax * 1.25)
            if ci == 0:
                ax.set_ylabel("events / bin (weighted)")
                axr.set_ylabel(f"{'pseudo-data' if pseudo else 'data'} / model")
            axr.set_xlabel("m_ll [GeV]")
        if leg_ax is not None:
            h, l = leg_ax.get_legend_handles_labels()
            fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.945),
                       ncol=len(l), fontsize=8, framealpha=0.9)
        fig.suptitle(f"m_ll closure{' (MC pseudo-data)' if pseudo else ''} "
                     f"— slices of {label}", y=0.998)
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        for p in _save_fig(fig, output_dir, f"mll_closure_{prefix}"):
            print(f"  wrote {p}")


def plot_theta_vs_eta(
    theta: np.ndarray,       # [n_eta, n_comp]
    sigma: "np.ndarray | None",  # [n_eta, n_comp] or None
    component_names: List[str],
    name: str,
    eta_edges: np.ndarray,    # [n_eta + 1]
    output_dir: str,
    edm: "float | None" = None,
    chi2_info=None,
    ref: "np.ndarray | None" = None,   # [n_eta, n_comp] reference (e.g. injected)
    band: "np.ndarray | None" = None,  # [n_eta, n_comp] φ-std for shaded band
    slices: "np.ndarray | None" = None,    # [n_eta, n_slices, n_comp] φ-slices
    slice_labels: "list | None" = None,    # length n_slices, e.g. ['φ=0', ...]
    slice_sigma: "np.ndarray | None" = None,  # [n_eta, n_slices, n_comp] per-cell σ
    sigma_band: bool = False,    # draw `sigma` as a continuous shaded band (MLP)
    sigma_total: "np.ndarray | None" = None,  # [n_eta, n_comp] data⊕flow ±1σ
    points: bool = False,   # 2-D binned: band → outer error bars, slices → markers
):
    n_eta, n_comp = theta.shape
    eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])
    dx = float(eta_centers[1] - eta_centers[0]) if n_eta > 1 else 1.0

    fig, axes = plt.subplots(n_comp, 1, sharex=True, figsize=(8, 2.5 * n_comp))
    if n_comp == 1:
        axes = [axes]
    slice_colors = ["C0", "C1", "C2", "C4", "C5"]   # skip C3 (reserved for ref)
    for i, ax in enumerate(axes):
        # Total (data⊕flow) ±1σ: a SHADED BAND for the continuous MLP curve,
        # but ERROR BARS for binned θ (the data-stat inner bar + marker is
        # drawn below; here we draw the wider OUTER total bar). Band drawn first
        # (behind) for MLP.
        if sigma_total is not None and sigma_band:
            ax.fill_between(
                eta_centers, theta[:, i] - sigma_total[:, i],
                theta[:, i] + sigma_total[:, i],
                color="C7", alpha=0.30, linewidth=0,
                label="±1σ (data⊕flow)")
        # φ-spread: MLP draws it as a shaded band; 2-D binned (`points`) as
        # wide light OUTER error bars (per-bin params — no continuity implied).
        if band is not None:
            if points:
                ax.errorbar(
                    eta_centers, theta[:, i], yerr=band[:, i], fmt="none",
                    ecolor="0.65", elinewidth=1.1, capsize=4,
                    label="±1σ over φ")
            else:
                ax.fill_between(
                    eta_centers, theta[:, i] - band[:, i],
                    theta[:, i] + band[:, i],
                    color="0.5", alpha=0.18, label="±1σ over φ", linewidth=0)
        # Individual φ slices: faint coloured lines (MLP) or x-offset markers
        # (2-D binned).
        if slices is not None:
            n_slices = slices.shape[1]
            for s_i in range(n_slices):
                lab = (slice_labels[s_i] if slice_labels and s_i < len(slice_labels)
                       else f"slice {s_i}")
                if points:
                    off = dx * 0.3 * ((s_i + 0.5) / n_slices - 0.5)
                    ax.errorbar(
                        eta_centers + off, slices[:, s_i, i],
                        yerr=(slice_sigma[:, s_i, i]
                              if slice_sigma is not None else None),
                        fmt="o", markersize=2.4, elinewidth=0.8, capsize=1.5,
                        color=slice_colors[s_i % len(slice_colors)],
                        alpha=0.65, label=lab)
                else:
                    ax.plot(eta_centers, slices[:, s_i, i],
                            color=slice_colors[s_i % len(slice_colors)],
                            lw=0.9, alpha=0.6, label=lab)
        # Main: per-bin error bars (binned θ) or a continuous φ-mean line (MLP).
        main_label = "fit" if band is None else "fit (φ-mean)"
        if sigma is not None and sigma_band:
            # MLP: the θ output is a continuous function of η, so show the
            # Fisher-propagated ±1σ as a SHADED BAND around the φ-mean curve
            # rather than discrete error bars (which would imply per-bin params).
            ax.fill_between(
                eta_centers, theta[:, i] - sigma[:, i], theta[:, i] + sigma[:, i],
                color="C0", alpha=0.30, linewidth=0, label="±1σ (Fisher)")
            ax.plot(eta_centers, theta[:, i], "-", color="k", lw=1.6,
                    label=main_label)
        elif sigma is not None:
            # binned: wider OUTER total (data⊕flow) error bar first, then the
            # data-stat error bar + marker on top — the gap is the flow.
            if sigma_total is not None:
                ax.errorbar(
                    eta_centers, theta[:, i], yerr=sigma_total[:, i], fmt="none",
                    ecolor="0.55", elinewidth=1.2, capsize=5,
                    label="±1σ (data⊕flow)")
            ax.errorbar(
                eta_centers, theta[:, i], yerr=sigma[:, i],
                fmt="o", color="k", markersize=4, capsize=2, label=main_label,
            )
        elif sigma_total is not None:
            # binned, flow-only (no --fisher): the total IS the error bar.
            ax.errorbar(
                eta_centers, theta[:, i], yerr=sigma_total[:, i],
                fmt="o", color="k", markersize=4, capsize=3,
                label="±1σ (data⊕flow)")
        else:
            ax.plot(eta_centers, theta[:, i], "o" if points else "o-",
                    color="k", markersize=4, label=main_label)
        ax.axhline(0, color="0.5", lw=0.8, ls=":")
        if ref is not None:
            ax.plot(eta_centers, ref[:, i], color="C3", ls="--", lw=1.3,
                    label="injected")
        if points:
            # Autoscale WITHOUT the slice error-bar extents (low-occupancy
            # cells have huge per-cell σ that would flatten everything else):
            # cover the mean ± its bars, the ref, and the slice point VALUES.
            cand = [theta[:, i]]
            for b in (sigma, band, sigma_total):
                if b is not None:
                    cand += [theta[:, i] - b[:, i], theta[:, i] + b[:, i]]
            if ref is not None:
                cand.append(ref[:, i])
            if slices is not None:
                cand.append(slices[:, :, i].ravel())
            v = np.concatenate([np.asarray(c, dtype=float).ravel() for c in cand])
            v = v[np.isfinite(v)]
            if v.size:
                lo, hi = float(v.min()), float(v.max())
                pad = 0.06 * ((hi - lo) or abs(hi) or 1.0)
                ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel(component_names[i])
        ax.grid(True, alpha=0.3)
    if (ref is not None or slices is not None or band is not None
            or sigma_total is not None):
        axes[0].legend(loc="best", fontsize=7, ncol=max(1,
            (1 + int(ref is not None) + int(band is not None)
             + int(sigma_total is not None)
             + (slices.shape[1] if slices is not None else 0)) // 4 + 1))
    axes[-1].set_xlabel("η-bin center")
    title = name
    if chi2_info is not None:
        chi2, dof, p = chi2_info
        cmp = "injected" if ref is not None else "0"
        title += (f"   (vs {cmp}: χ²/dof = {chi2:.1f}/{dof} = {chi2 / max(dof, 1):.2f}, "
                  f"p = {p:.3g})")
    axes[0].set_title(title)
    if edm is not None:
        fig.text(0.995, 0.005, f"EDM = {edm:.2e}", ha="right", va="bottom",
                 fontsize=8, color="0.4")
    fig.tight_layout()
    for p in _save_fig(fig, output_dir, name):
        print(f"  wrote {p}")


def plot_theta_grid_etaphi(grid, sigma, component_names, name, eta_edges,
                           n_phi, output_dir, chi2_info=None, ref=None,
                           edm=None):
    """θ as η×φ heatmaps for the 2-D binned (``--theta-2d-binned``) mode.

    ``grid`` ``[n_eta, n_phi, n_comp]`` physical values; ``sigma`` same shape
    (or None); ``ref`` ``[n_eta, n_phi, n_comp]`` / flat reference (injected
    truth or 0). One COLUMN per component: top row the fitted value, bottom row
    the pull ``(value − ref)/σ`` when both ``sigma`` and ``ref`` are available
    (else σ if only σ is present). φ runs over [−π, π); η over the bin edges."""
    import matplotlib.pyplot as plt
    n_eta, n_phi_g, n_comp = grid.shape
    eta_lo, eta_hi = float(eta_edges[0]), float(eta_edges[-1])
    extent = [eta_lo, eta_hi, -np.pi, np.pi]
    has_pull = sigma is not None and ref is not None
    nrow = 2 if (has_pull or sigma is not None) else 1
    fig, axes = plt.subplots(nrow, n_comp, figsize=(4.2 * n_comp, 3.4 * nrow),
                             squeeze=False)
    refg = None if ref is None else np.asarray(ref).reshape(n_eta, n_phi_g, n_comp)
    # heatmaps show φ on y, η on x: transpose [n_eta,n_phi] → [n_phi,n_eta].
    for c in range(n_comp):
        v = grid[:, :, c].T
        im = axes[0][c].imshow(v, origin="lower", aspect="auto", extent=extent,
                               cmap="RdBu_r",
                               vmax=np.abs(v).max() or 1.0,
                               vmin=-(np.abs(v).max() or 1.0))
        axes[0][c].set_title(component_names[c])
        axes[0][c].set_xlabel("η"); axes[0][c].set_ylabel("φ")
        fig.colorbar(im, ax=axes[0][c], fraction=0.046, pad=0.04)
        if nrow == 2:
            if has_pull:
                pull = ((grid[:, :, c] - refg[:, :, c])
                        / np.where(sigma[:, :, c] > 0, sigma[:, :, c], np.inf)).T
                lab, cmap, lim = "pull (val−ref)/σ", "RdBu_r", 5.0
                im2 = axes[1][c].imshow(pull, origin="lower", aspect="auto",
                                        extent=extent, cmap=cmap,
                                        vmin=-lim, vmax=lim)
            else:
                sg = sigma[:, :, c].T
                im2 = axes[1][c].imshow(sg, origin="lower", aspect="auto",
                                        extent=extent, cmap="viridis",
                                        vmin=0.0)
                lab = "σ"
            axes[1][c].set_title(lab)
            axes[1][c].set_xlabel("η"); axes[1][c].set_ylabel("φ")
            fig.colorbar(im2, ax=axes[1][c], fraction=0.046, pad=0.04)
    if chi2_info is not None:
        chi2, dof, pval = chi2_info
        fig.suptitle(f"{name}   χ²/dof = {chi2:.1f}/{dof} = "
                     f"{chi2 / max(dof, 1):.2f},  p = {pval:.3g}")
    if edm is not None:
        fig.text(0.995, 0.005, f"EDM = {edm:.2e}", ha="right", va="bottom",
                 fontsize=8, color="0.4")
    fig.tight_layout()
    for p in _save_fig(fig, output_dir, name):
        print(f"  wrote {p}")


def plot_theta_grid_projections(grid, sigma, component_names, prefix,
                                eta_edges, output_dir, ref=None, edm=None,
                                cov=None, n_phi_slices=4, n_eta_slices=5):
    """1-D projections of the 2-D binned θ grid (``--theta-2d-binned``):
    φ-averaged θ vs η and η-averaged θ vs φ, drawn as POINTS — black markers
    with inner error bars (σ of the weighted mean) and light outer error bars
    (the spread over the averaged coordinate) plus a few slices of it as
    offset markers — easier to read/compare than the η×φ heatmaps.

    ``grid``/``sigma``: ``[n_eta, n_phi, n_comp]`` physical values (σ may be
    None); ``ref``: injected-truth reference, per-cell ``[n_eta, n_phi,
    n_comp]`` (carries the f_φ sinusoid under --inject-nonuniform) or per-η
    ``[n_eta, n_comp]`` (φ-uniform). Averages are INVERSE-VARIANCE weighted
    when σ is available (so empty/unconstrained cells don't dominate), else
    uniform — the ref is projected with the SAME weights so the comparison is
    apples-to-apples. ``cov``: flat per-cell covariance
    ``[n_cells·n_comp, n_cells·n_comp]`` (η-major cell order) — when given,
    the means are treated as the linear map m = Wθ, so Cov(m) = W·C·Wᵀ with
    per-cell correlations fully included: the error bars on the means are
    √diag(W·C·Wᵀ) and each projection gets a PROPER χ²/dof vs the projected
    ref. (Only the WEIGHTS use the diagonal σ; without ``cov`` the σ on the
    mean falls back to the uncorrelated inverse-variance 1/√Σw.)"""
    n_eta, n_phi, n_comp = grid.shape
    if sigma is not None:
        w = 1.0 / np.clip(sigma, 1e-30, None) ** 2
        w = np.where(np.isfinite(w), w, 0.0)
    else:
        w = np.ones_like(grid)
    ref_g = None
    if ref is not None:
        ref = np.asarray(ref)
        ref_g = (np.broadcast_to(ref[:, None, :], grid.shape)
                 if ref.ndim == 2 else ref)

    def _wmean(a, axis):
        wsum = w.sum(axis=axis)
        return (w * a).sum(axis=axis) / np.where(wsum > 0, wsum, 1.0)

    def _proj(axis):
        mean = _wmean(grid, axis)
        spread = np.sqrt(np.clip(
            _wmean((grid - np.expand_dims(mean, axis)) ** 2, axis), 0.0, None))
        sig = (np.where(w.sum(axis=axis) > 0,
                        1.0 / np.sqrt(np.clip(w.sum(axis=axis), 1e-300, None)),
                        np.nan)
               if sigma is not None else None)
        return mean, spread, sig

    def _proj_cov(axis, out):
        # FULL covariance of the projected weighted means: m = Wθ is linear,
        # so Cov(m) = W·C·Wᵀ — per-cell correlations included. Returns the
        # flat [n_out, n_out] matrix (row-major in the `out` index order).
        if cov is None:
            return None
        wsum = w.sum(axis=axis, keepdims=True)
        wn = w / np.where(wsum > 0, wsum, 1.0)
        c6 = np.asarray(cov, dtype=np.float64).reshape(
            n_eta, n_phi, n_comp, n_eta, n_phi, n_comp)
        cm = np.einsum('epc,epcfqd,fqd->' + out, wn, c6, wn, optimize=True)
        n = int(np.sqrt(cm.size))
        return cm.reshape(n, n)

    def _proj_chi2(cm, mean, ref_m):
        # χ² of the projected means vs the (identically projected) ref — or 0
        # if no injection. Frozen / zero-variance directions are dropped by
        # the positive-variance eigendecomposition in _chi2_compat_zero.
        if cm is None:
            return None
        r = mean - (ref_m if ref_m is not None else 0.0)
        return _chi2_compat_zero(r.reshape(-1), cm)

    phi_ctr = (np.arange(n_phi) + 0.5) * (2.0 * np.pi / n_phi) - np.pi
    eta_ctr = 0.5 * (np.asarray(eta_edges[:-1]) + np.asarray(eta_edges[1:]))

    # φ-averaged θ(η): error bars = σ of the weighted φ-mean, grey band = the
    # φ STRUCTURE spread, faint curves = individual φ slices.
    mean_e, band_e, sig_e = _proj(1)
    ref_e = None if ref_g is None else _wmean(ref_g, 1)
    cm_e = _proj_cov(1, 'ecfd')
    if cm_e is not None:
        # error bar on the mean from the FULL projected covariance
        # (correlations included) — supersedes the uncorrelated 1/√Σw.
        sig_e = np.sqrt(np.clip(np.diag(cm_e), 0.0, None)).reshape(mean_e.shape)
    chi2_e = _proj_chi2(cm_e, mean_e, ref_e)
    sl_p = np.unique(np.linspace(0, n_phi - 1, n_phi_slices).round().astype(int))
    plot_theta_vs_eta(
        mean_e, sig_e, component_names, f"{prefix}_vs_eta", eta_edges,
        output_dir, edm=edm, ref=ref_e, band=band_e, chi2_info=chi2_e,
        slices=grid[:, sl_p, :],
        slice_sigma=(sigma[:, sl_p, :] if sigma is not None else None),
        slice_labels=[f"φ≈{phi_ctr[s]:+.2f}" for s in sl_p], points=True)

    # η-averaged θ(φ): same construction transposed.
    mean_p, band_p, sig_p = _proj(0)
    ref_p = None if ref_g is None else _wmean(ref_g, 0)
    cm_p = _proj_cov(0, 'pcqd')
    if cm_p is not None:
        sig_p = np.sqrt(np.clip(np.diag(cm_p), 0.0, None)).reshape(mean_p.shape)
    chi2_p = _proj_chi2(cm_p, mean_p, ref_p)
    sl_e = np.unique(np.linspace(0, n_eta - 1, n_eta_slices).round().astype(int))
    plot_theta_vs_phi(
        phi_ctr, grid[sl_e], component_names, f"{prefix}_vs_phi",
        eta_ctr[sl_e], output_dir, ref=ref_p,
        eta_mean=mean_p, eta_band=band_p, band_label="±1σ over η",
        fisher_sigma=sig_p, points=True, chi2_info=chi2_p,
        slice_sigma=(sigma[sl_e] if sigma is not None else None))
    for tag, ci in (("φ-average vs η", chi2_e), ("η-average vs φ", chi2_p)):
        if ci is not None:
            chi2, dof, pval = ci
            print(f"  {prefix} {tag}: χ²/dof = {chi2:.1f}/{dof} = "
                  f"{chi2 / max(dof, 1):.2f}, p = {pval:.3g}")


def plot_theta_vs_phi(
    phi_grid: np.ndarray,        # [n_phi] φ sample points
    theta: np.ndarray,           # [n_eta_slc, n_phi, n_comp] net output at η-slices
    component_names: List[str],
    name: str,
    eta_slice_vals: np.ndarray,  # [n_eta_slc] the |η| (centre) of each slice curve
    output_dir: str,
    ref: "np.ndarray | None" = None,   # [n_phi, n_comp] φ-dependent injected ref
    eta_mean: "np.ndarray | None" = None,   # [n_phi, n_comp] η-averaged curve
    eta_band: "np.ndarray | None" = None,   # [n_phi, n_comp] ±1σ band on the η-mean
    band_label: str = "±1σ over η",
    fisher_sigma: "np.ndarray | None" = None,  # [n_phi, n_comp] Fisher ±1σ(φ) on η-mean
    points: bool = False,   # 2-D binned: bands → error bars, slices → markers
    slice_sigma: "np.ndarray | None" = None,  # [n_eta_slc, n_phi, n_comp] per-cell σ
    chi2_info=None,         # (χ², dof, p) of the η-mean vs ref (title suffix)
):
    """θ output (A,e,M or a,c) as a function of φ, one (faint) curve per
    representative η slice — the φ-direction companion to plot_theta_vs_eta. Only
    meaningful for --theta-mlp (the binned θ has no φ dependence). Reveals the
    net's learned φ structure (and, in --inject-nonuniform validation, whether it
    tracks the injected sinusoidal-φ modulation).

    ``eta_mean``/``eta_band``: an η-AVERAGED φ curve (bold black) with a shaded
    ±1σ-over-η band (the η STRUCTURE spread). ``fisher_sigma`` [n_phi, n_comp]:
    the STATISTICAL ±1σ(φ) of the η-averaged output from the full 2-D (η,φ)
    Fisher covariance — φ-resolved (varies with φ)."""
    n_eta_slc, n_phi, n_comp = theta.shape
    fig, axes = plt.subplots(n_comp, 1, sharex=True, figsize=(8, 2.5 * n_comp))
    if n_comp == 1:
        axes = [axes]
    eta_colors = ["C0", "C1", "C2", "C4", "C5", "C6"]
    dphi = float(phi_grid[1] - phi_grid[0]) if n_phi > 1 else 1.0
    for i, ax in enumerate(axes):
        for s in range(n_eta_slc):
            if points:
                off = dphi * 0.3 * ((s + 0.5) / n_eta_slc - 0.5)
                ax.errorbar(
                    phi_grid + off, theta[s, :, i],
                    yerr=(slice_sigma[s, :, i]
                          if slice_sigma is not None else None),
                    fmt="o", markersize=2.4, elinewidth=0.8, capsize=1.5,
                    color=eta_colors[s % len(eta_colors)], alpha=0.6,
                    label=f"η≈{eta_slice_vals[s]:+.1f}")
            else:
                ax.plot(phi_grid, theta[s, :, i],
                        color=eta_colors[s % len(eta_colors)], lw=0.9,
                        alpha=0.55, label=f"η≈{eta_slice_vals[s]:+.1f}")
        # η-averaged curve: bold line + shaded bands (MLP) or black markers
        # with inner (statistical) + outer (η-spread) error bars (2-D binned).
        if eta_mean is not None and points:
            if eta_band is not None:
                ax.errorbar(
                    phi_grid, eta_mean[:, i], yerr=eta_band[:, i], fmt="none",
                    ecolor="0.65", elinewidth=1.1, capsize=4, label=band_label)
            fs = (np.nan_to_num(fisher_sigma[:, i], nan=0.0)
                  if fisher_sigma is not None and np.any(fisher_sigma[:, i] > 0)
                  else None)
            ax.errorbar(
                phi_grid, eta_mean[:, i], yerr=fs, fmt="o", color="k",
                markersize=4, capsize=2, label="η-mean")
        elif eta_mean is not None:
            # Fisher statistical band, φ-resolved — drawn first as the reference.
            if fisher_sigma is not None and np.any(fisher_sigma[:, i] > 0):
                fs = np.nan_to_num(fisher_sigma[:, i], nan=0.0)
                ax.fill_between(
                    phi_grid, eta_mean[:, i] - fs, eta_mean[:, i] + fs,
                    color="C0", alpha=0.20, linewidth=0, label="±1σ (Fisher)")
            if eta_band is not None:
                ax.fill_between(
                    phi_grid, eta_mean[:, i] - eta_band[:, i],
                    eta_mean[:, i] + eta_band[:, i],
                    color="0.4", alpha=0.25, linewidth=0, label=band_label)
            ax.plot(phi_grid, eta_mean[:, i], "-", color="k", lw=1.8,
                    label="η-mean")
        if ref is not None:
            ax.plot(phi_grid, ref[:, i], color="C3", ls="--", lw=1.3,
                    label="injected")
        ax.axhline(0, color="0.5", lw=0.8, ls=":")
        if points:
            # Autoscale WITHOUT the slice error-bar extents (see
            # plot_theta_vs_eta): mean ± its bars, ref, slice point values.
            cand = [theta[:, :, i].ravel()]
            if eta_mean is not None:
                cand.append(eta_mean[:, i])
                if eta_band is not None:
                    cand += [eta_mean[:, i] - eta_band[:, i],
                             eta_mean[:, i] + eta_band[:, i]]
                if fisher_sigma is not None:
                    fs0 = np.nan_to_num(fisher_sigma[:, i], nan=0.0)
                    cand += [eta_mean[:, i] - fs0, eta_mean[:, i] + fs0]
            if ref is not None:
                cand.append(ref[:, i])
            v = np.concatenate([np.asarray(c, dtype=float).ravel() for c in cand])
            v = v[np.isfinite(v)]
            if v.size:
                lo, hi = float(v.min()), float(v.max())
                pad = 0.06 * ((hi - lo) or abs(hi) or 1.0)
                ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel(component_names[i])
        ax.grid(True, alpha=0.3)
    nleg = (n_eta_slc + int(ref is not None) + 2 * int(eta_mean is not None)
            + int(fisher_sigma is not None))
    axes[0].legend(loc="best", fontsize=7, ncol=max(1, nleg // 4 + 1))
    axes[-1].set_xlabel("φ [rad]")
    title = name
    if chi2_info is not None:
        chi2, dof, p = chi2_info
        cmp = "injected" if ref is not None else "0"
        title += (f"   (vs {cmp}: χ²/dof = {chi2:.1f}/{dof} = "
                  f"{chi2 / max(dof, 1):.2f}, p = {p:.3g})")
    axes[0].set_title(title)
    fig.tight_layout()
    for p in _save_fig(fig, output_dir, name):
        print(f"  wrote {p}")


def plot_fisher_correlation(cov: np.ndarray, output_dir: str):
    """72×72 correlation heatmap, with η-bin grid lines + (A,e,M) tick labels."""
    n = cov.shape[0]
    d = np.sqrt(np.diag(cov))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(d, d)
    corr = np.where(np.isfinite(corr), corr, 0.0)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    # η-bin grid lines (3 cols per η-bin).
    for k in range(3, n, 3):
        ax.axhline(k - 0.5, color="k", lw=0.3, alpha=0.4)
        ax.axvline(k - 0.5, color="k", lw=0.3, alpha=0.4)
    ax.set_xticks(np.arange(1, n, 6))
    ax.set_xticklabels([f"η{j//3}" for j in range(1, n, 6)], fontsize=7, rotation=90)
    ax.set_yticks(np.arange(1, n, 6))
    ax.set_yticklabels([f"η{j//3}" for j in range(1, n, 6)], fontsize=7)
    ax.set_title("Fisher correlation matrix (A, e, M per η-bin)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    for p in _save_fig(fig, output_dir, "fisher_correlation"):
        print(f"  wrote {p}")


def plot_cov_corr(cov: np.ndarray, labels, n_scale: int, output_dir: str,
                  edm: "float | None" = None):
    """Side-by-side covariance (symmetric-log) + correlation ([-1,1]) heatmaps of
    the FULL joint parameter covariance (θ_scale + the active θ_smear), with the
    scale/smear block separator and sparse per-parameter tick labels.

    The covariance mixes units across blocks (A, e, M, smear), spanning many
    orders of magnitude, so it is shown on a signed-log (SymLog) colour scale;
    the correlation is the dimensionless, directly-readable companion.
    """
    import matplotlib.colors as mcolors

    n = cov.shape[0]
    d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(d, d)
    corr = np.where(np.isfinite(corr), corr, 0.0)

    fig, (axc, axr) = plt.subplots(1, 2, figsize=(15, 6.5))
    amax = float(np.abs(cov).max()) or 1.0
    norm = mcolors.SymLogNorm(linthresh=amax * 1e-6, vmin=-amax, vmax=amax, base=10)
    im0 = axc.imshow(cov, cmap="RdBu_r", norm=norm)
    axc.set_title("covariance (signed log)")
    fig.colorbar(im0, ax=axc, fraction=0.046, pad=0.04)
    im1 = axr.imshow(corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0)
    axr.set_title("correlation")
    fig.colorbar(im1, ax=axr, fraction=0.046, pad=0.04)

    # Sparse ticks straight from the parameter labels (~16 across the axis).
    if labels is not None and len(labels) == n:
        step = max(1, n // 16)
        ticks = list(range(0, n, step))
        tlabels = [labels[i] for i in ticks]
    else:
        ticks, tlabels = [], []
    for ax in (axc, axr):
        if 0 < n_scale < n:   # separate the θ_scale and θ_smear blocks
            ax.axhline(n_scale - 0.5, color="k", lw=1.0)
            ax.axvline(n_scale - 0.5, color="k", lw=1.0)
        ax.set_xticks(ticks)
        ax.set_xticklabels(tlabels, fontsize=6, rotation=90)
        ax.set_yticks(ticks)
        ax.set_yticklabels(tlabels, fontsize=6)
    suptitle = "parameter covariance / correlation (θ_scale + active θ_smear)"
    if edm is not None:
        suptitle += f"    EDM = {edm:.2e}"
    fig.suptitle(suptitle)
    fig.tight_layout()
    for p in _save_fig(fig, output_dir, "covariance_correlation"):
        print(f"  wrote {p}")


def plot_mc_closure(
    evals, m_centers_np, eta_slice_edges, output_dir: str,
):
    """MC closure: forward-folded-MC histogram (points) vs flow-density
    curve (line), both at the fitted scale + smearing. One PANEL FIGURE per
    slice dimension (columns = slices): inclusive + |η_+| (fixed edges), ρ (pt
    asymmetry) tertiles, cos α (3-D opening angle) tertiles.

    Empirical: histogram of m_fold_e — the MC reco forward-folded through
    scale then smearing at the fitted θ_scale / θ_smear — weighted by w_mc.

    Model curve: ``Δm · Σ_e w_e · p_flow(m_bin_j | y_e, fitted θ_scale+θ_smear)`` —
    per-event × per-bin flow density at the fitted nuisances, summed over
    MC events. The two should agree when the flow has correctly learned the
    forward-fold distribution; this is the analogue of the dashed orange
    "model signal" curve on the data plot (no Bernstein bkg, MC is
    signal-only).
    """
    bin_width = float(m_centers_np[1] - m_centers_np[0])
    m_edges = np.concatenate([
        [m_centers_np[0] - bin_width / 2],
        m_centers_np[:-1] + bin_width / 2,
        [m_centers_np[-1] + bin_width / 2],
    ])
    cont = bool(evals.get("continuity", False))
    injected = bool(evals.get("injected", False))

    def _draw_panel(ax, axr, mc_mask):
        """Draw one MC-closure slice into (main, ratio) axes and set the ratio
        panel's own autoscaled y-range; return ymax."""
        mc_hist, _ = np.histogram(
            evals["mll_mc_fold"][mc_mask], bins=m_edges,
            weights=evals["w_mc"][mc_mask])
        w = evals["w_mc"][mc_mask][:, None]
        model_curve = bin_width * (evals["pred_signal_mc"][mc_mask] * w).sum(axis=0)
        mc_label = ("MC (shifted+smeared, fitted θ)" if cont
                    else "MC (scale+smear folded)")
        model_label = ("flow (folded, tilt at fitted θ)" if cont
                       else "flow (at fitted scale+smear)")
        has_nom = cont and evals.get("mll_mc_nominal", np.zeros((0,))).size > 0
        if has_nom:
            nom_hist, _ = np.histogram(
                evals["mll_mc_nominal"][mc_mask], bins=m_edges,
                weights=evals["w_mc"][mc_mask])
            nom_curve = bin_width * (evals["pred_nominal_mc"][mc_mask] * w).sum(axis=0)
        show_pseudo = (cont and injected
                       and evals.get("mll_mc_pseudodata", np.zeros((0,))).size > 0)
        if show_pseudo:
            pseudo_hist, _ = np.histogram(
                evals["mll_mc_pseudodata"][mc_mask], bins=m_edges,
                weights=evals["w_mc"][mc_mask])
        show_inj_curve = (cont and injected
                          and evals.get("pred_signal_mc_at_inj", np.zeros((0,))).size > 0)
        if show_inj_curve:
            inj_curve = bin_width * (
                evals["pred_signal_mc_at_inj"][mc_mask] * w).sum(axis=0)

        if has_nom:
            ax.step(m_edges[:-1], nom_hist, where="post", color="0.6", lw=1.0,
                    label="MC (nominal, θ=0)", zorder=2)
            ax.plot(m_centers_np, nom_curve, color="C0", ls=":", lw=1.3,
                    label="flow (nominal p₀, θ=0)", zorder=2)
        if show_pseudo:
            ax.step(m_edges[:-1], pseudo_hist, where="post", color="C2", lw=1.3,
                    label="pseudo-data (injected θ)", zorder=2)
        ax.errorbar(m_centers_np, mc_hist, yerr=np.sqrt(np.abs(mc_hist)),
                    fmt="o", color="k", markersize=3, label=mc_label, zorder=3)
        ax.plot(m_centers_np, model_curve, color="C1", ls="--", lw=1.5,
                label=model_label)
        if show_inj_curve:
            ax.plot(m_centers_np, inj_curve, color="C3", ls="-.", lw=1.3,
                    label="flow (at INJECTED θ)", zorder=2)

        # Ratio panel: everything relative to the folded flow (the model).
        denom = np.where(model_curve > 0, model_curve, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = mc_hist / denom
            ratio_err = np.sqrt(np.abs(mc_hist)) / denom
        axr.errorbar(m_centers_np, ratio, yerr=ratio_err, fmt="o",
                     color="k", markersize=3, zorder=3)
        if show_pseudo:
            with np.errstate(divide="ignore", invalid="ignore"):
                axr.step(m_edges[:-1], pseudo_hist / denom, where="post",
                         color="C2", lw=1.3, zorder=2)
        if show_inj_curve:
            with np.errstate(divide="ignore", invalid="ignore"):
                axr.plot(m_centers_np, inj_curve / denom, color="C3", ls="-.",
                         lw=1.3, zorder=2)
        if has_nom:
            with np.errstate(divide="ignore", invalid="ignore"):
                axr.step(m_edges[:-1], nom_hist / denom, where="post",
                         color="0.6", lw=1.0, zorder=2)
                axr.plot(m_centers_np, nom_curve / denom, color="C0", ls=":",
                         lw=1.3, zorder=2)
        axr.axhline(1.0, color="C1", lw=1)   # folded flow (model) = reference
        # ratio-panel autoscale: robust symmetric half-window about 1, sized
        # from well-populated bins (see plot_mll_closure). Each slice gets its
        # OWN window.
        good = np.isfinite(ratio) & (
            mc_hist > 0.10 * max(float(mc_hist.max()), 1e-30))
        if good.any():
            dev = float(np.percentile(np.abs(ratio[good] - 1.0), 90))
            half = max(dev, float(np.median(ratio_err[good]))) * 1.3
            h = min(max(half, 0.03), 0.4)
            axr.set_ylim(1.0 - h, 1.0 + h)
        else:
            axr.set_ylim(0.6, 1.4)
        ymax = max(float((mc_hist + np.sqrt(np.abs(mc_hist))).max()),
                   float(model_curve.max()),
                   float(nom_hist.max()) if has_nom else 0.0,
                   float(nom_curve.max()) if has_nom else 0.0,
                   float(pseudo_hist.max()) if show_pseudo else 0.0,
                   float(inj_curve.max()) if show_inj_curve else 0.0)
        return ymax

    for prefix, label, fmt, _dv, mv, cols in _closure_slice_dims(
            evals, eta_slice_edges):
        ncol = len(cols)
        fig, axes = plt.subplots(
            2, ncol, figsize=(max(6.0, 4.3 * ncol), 6.2), squeeze=False,
            sharex="col", gridspec_kw={"height_ratios": [3, 1]})
        leg_ax = None
        for ci, (tag, slice_def) in enumerate(cols):
            ax, axr = axes[0, ci], axes[1, ci]
            mc_mask = _select_slice(mv, slice_def) if mv.size else np.zeros((0,), bool)
            if mc_mask.sum() == 0:
                ax.set_visible(False); axr.set_visible(False); continue
            ymax = _draw_panel(ax, axr, mc_mask)
            if leg_ax is None:
                leg_ax = ax
            ttl = ("inclusive" if slice_def is None
                   else f"{label} ∈ [{slice_def[0]:{fmt}}, {slice_def[1]:{fmt}}]")
            ax.set_title(ttl, fontsize=9)
            if ymax > 0:
                ax.set_ylim(0, ymax * 1.25)
            if ci == 0:
                ax.set_ylabel("events / bin (weighted)")
                axr.set_ylabel("ratio to folded flow")
            axr.set_xlabel("m_ll [GeV]")
        if leg_ax is not None:
            h, l = leg_ax.get_legend_handles_labels()
            fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.945),
                       ncol=min(len(l), 4), fontsize=8, framealpha=0.9)
        fig.suptitle(f"MC closure — slices of {label}", y=0.998)
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        for p in _save_fig(fig, output_dir, f"mc_closure_{prefix}"):
            print(f"  wrote {p}")


def plot_bkg_fractions(agg, eta_edges, output_dir, inject_bkg=None,
                       bkg_model="bernstein"):
    """Background-fraction closure vs η and vs φ: the Σw-weighted per-bin
    means of the fitted MLP background components (and their TOTAL), filled
    PER MUON over the data-branch events, overlaid with the EMPIRICAL
    injected fractions (per-bin Σw of the truth-labelled background events —
    present only with --inject-bkg-*) and the constant injected values
    (dashed). Per-component injected references are drawn only for the
    degree-1 Bernstein model, whose components match the injection's; for
    higher degrees / 'exp' the TOTAL closure is the meaningful comparison.
    Accumulator columns: (Σw, Σw·f_k…, Σw·f_bkg_tot, Σw·[lab==1],
    Σw·[lab==2][, Σw·s(c) for 'exp' — drawn as the per-bin mean slope on a
    twin axis]). Writes bkg_fraction_eta and bkg_fraction_phi."""
    import matplotlib.pyplot as plt
    eta_edges = np.asarray(eta_edges, dtype=float)
    has_slope = (bkg_model == "exp")
    n_bkg = agg["eta"].shape[1] - 4 - (1 if has_slope else 0)
    comp_refs = (bkg_model == "bernstein" and n_bkg == 2)
    for tag, acc, centers, xlabel in (
            ("eta", agg["eta"], 0.5 * (eta_edges[:-1] + eta_edges[1:]), r"$\eta_\mu$"),
            ("phi", agg["phi"],
             (-np.pi + (np.arange(agg["phi"].shape[0]) + 0.5)
              * (2.0 * np.pi / agg["phi"].shape[0])), r"$\phi_\mu$")):
        sw = np.clip(acc[:, 0], 1e-30, None)
        fk = acc[:, 1:1 + n_bkg] / sw[:, None]
        ftot = acc[:, 1 + n_bkg] / sw
        lab_off = 2 + n_bkg
        e0m, e1m = acc[:, lab_off] / sw, acc[:, lab_off + 1] / sw
        etot = e0m + e1m
        err = lambda p: np.sqrt(np.clip(p * (1 - p), 0, None)
                                / np.clip(sw, 1, None))
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        comp_colors = ["C0", "C3", "C4", "C5", "C6", "C8"]
        for k in range(n_bkg):
            lab = (r"fitted $\langle f_0\rangle$ (falling)" if (comp_refs and k == 0)
                   else r"fitted $\langle f_1\rangle$ (rising)" if (comp_refs and k == 1)
                   else f"fitted comp {k}")
            ax.plot(centers, fk[:, k], "o-", ms=3, lw=1.0,
                    color=comp_colors[k % len(comp_colors)], label=lab)
        ax.plot(centers, ftot, "k^-", ms=4, lw=1.4,
                label=r"fitted $\langle f_{bkg}\rangle$ (total)")
        has_inj = bool(acc[:, lab_off].sum() > 0 or acc[:, lab_off + 1].sum() > 0)
        if has_inj:
            if comp_refs:
                ax.errorbar(centers, e0m, yerr=err(e0m), fmt=".", color="C0",
                            alpha=0.55, capsize=2, label="injected (emp., comp 0)")
                ax.errorbar(centers, e1m, yerr=err(e1m), fmt=".", color="C3",
                            alpha=0.55, capsize=2, label="injected (emp., comp 1)")
            ax.errorbar(centers, etot, yerr=err(etot), fmt=".", color="k",
                        alpha=0.55, capsize=2, label="injected (emp., total)")
        if inject_bkg is not None:
            if comp_refs:
                ax.axhline(inject_bkg[0], color="C0", ls="--", lw=1, alpha=0.7)
                ax.axhline(inject_bkg[1], color="C3", ls="--", lw=1, alpha=0.7)
            ax.axhline(inject_bkg[0] + inject_bkg[1], color="k", ls="--",
                       lw=1, alpha=0.7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("background fraction")
        ax.set_ylim(bottom=0.0)
        ax.legend(fontsize=8, ncol=2)
        if has_slope:
            ax2 = ax.twinx()
            ax2.plot(centers, acc[:, -1] / sw, color="0.4", ls="-.", lw=1.2,
                     label=r"fitted $\langle s(c)\rangle$")
            ax2.set_ylabel("exp slope s = λ·width", color="0.35", fontsize=9)
            ax2.tick_params(axis="y", labelcolor="0.35")
            h2, l2 = ax2.get_legend_handles_labels()
            h1, l1 = ax.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, fontsize=8, ncol=2)
        ttl = "background-fraction closure (per-muon filled, data branch)"
        if bkg_model == "exp":
            ttl += "  [exp: fraction + conditioning-dependent slope]"
        elif bkg_model == "bernstein" and n_bkg > 2:
            ttl += f"  [bernstein deg {n_bkg - 1}]"
        ax.set_title(ttl, fontsize=10)
        fig.tight_layout()
        _save_fig(fig, output_dir, f"bkg_fraction_{tag}")
    print("wrote bkg_fraction_eta / bkg_fraction_phi")


def plot_pulls(
    evals, m_centers_np, eta_slice_edges, m_lo: float, m_hi: float,
    output_dir: str,
):
    """Per-bin (data − model)/√(Σw²) pulls; expect ~N(0,1) if model is OK. The
    denominator is the WEIGHTED per-bin error √(Σw²) (not the Poisson √model,
    which mis-scales pulls for weighted MC pseudo-data — inflating the std and
    biasing the mean)."""
    pseudo = bool(evals.get("mc_as_data", False))
    label = "pseudo-data" if pseudo else "data"
    m_edges = np.concatenate([
        [m_centers_np[0] - evals["bin_width"] / 2],
        m_centers_np[:-1] + evals["bin_width"] / 2,
        [m_centers_np[-1] + evals["bin_width"] / 2],
    ])

    slices = [("inclusive", None)] + [
        (f"eta{i}", (eta_slice_edges[i], eta_slice_edges[i + 1]))
        for i in range(len(eta_slice_edges) - 1)
    ]

    for tag, slice_def in slices:
        data_mask = (
            _select_slice(np.abs(evals["eta_data"]), slice_def)
            if evals["eta_data"].size else np.zeros((0,), bool)
        )
        mc_mask = (
            _select_slice(np.abs(evals["eta_mc"]), slice_def)
            if evals["eta_mc"].size else np.zeros((0,), bool)
        )
        if data_mask.sum() == 0:
            continue

        w_sel = evals["w_data"][data_mask]
        data_hist, _ = np.histogram(
            evals["mll_data"][data_mask], bins=m_edges, weights=w_sel,
        )
        # Per-bin error for WEIGHTED (pseudo-)data is √(Σw²), NOT √(Σw)=√model:
        # with event weights Var[bin] = Σw² ≠ mean, so the Poisson √model
        # normalisation mis-scales the pull (inflates the std away from 1 and
        # biases the mean). The model curve is the smooth expectation, so its own
        # variance is negligible vs the data's; floor by 1 effective count's
        # worth (the bin's mean weight²) so empty/√0 bins don't blow up.
        sumw2_hist, _ = np.histogram(
            evals["mll_data"][data_mask], bins=m_edges, weights=w_sel ** 2,
        )
        signal_curve, bkg_curve = _model_pred_histograms(
            evals, m_edges, m_lo, m_hi, data_mask,
        )
        total_curve = signal_curve + bkg_curve

        with np.errstate(divide="ignore", invalid="ignore"):
            err = np.sqrt(np.where(sumw2_hist > 0, sumw2_hist, np.nan))
            pulls = (data_hist - total_curve) / err
        pulls_finite = pulls[np.isfinite(pulls)]

        fig, (ax_b, ax_h) = plt.subplots(1, 2, figsize=(10, 4))
        ax_b.stem(m_centers_np, pulls, markerfmt="ko", basefmt="grey", linefmt="k-")
        ax_b.axhline(0, color="grey", lw=0.5)
        ax_b.set_xlabel("m_ll [GeV]")
        ax_b.set_ylabel(f"({label} − model) / √(Σw²)")
        ax_b.set_title(
            f"per-bin pulls{' (MC pseudo-data)' if pseudo else ''} — {tag}"
            + (f"  |η₊| ∈ [{slice_def[0]:.1f}, {slice_def[1]:.1f}]" if slice_def else "")
        )

        ax_h.hist(pulls_finite, bins=20, range=(-5, 5),
                  histtype="step", color="k", lw=1.2)
        # Overlay N(0,1) reference scaled to integral=n_bins.
        x = np.linspace(-5, 5, 200)
        ax_h.plot(
            x, len(pulls_finite) * 10 / 20 * np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi),
            color="C0", lw=1, label="N(0,1)",
        )
        ax_h.set_xlabel("pull")
        ax_h.set_ylabel("bins")
        ax_h.legend(fontsize=8)
        ax_h.text(
            0.05, 0.95,
            f"mean={np.nanmean(pulls_finite):+.2f}\n"
            f"std={np.nanstd(pulls_finite):.2f}\n"
            f"n={len(pulls_finite)}",
            transform=ax_h.transAxes, va="top", fontsize=8,
        )

        fig.tight_layout()
        for p in _save_fig(fig, output_dir, f"mll_pulls_{tag}"):
            print(f"  wrote {p}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="jpsi_mass_fit_diagnostics",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--checkpoint", required=True,
                   help="Path to checkpoint_best.pt from train_jpsi_mass_fit.")
    p.add_argument("--shards", required=True,
                   help="Path to the shard directory (or list).")
    p.add_argument("--fisher", default=None,
                   help="Path to fisher_info.pt (optional; enables θ_scale ±σ "
                   "bands and the correlation heatmap).")
    p.add_argument("--flow-uncertainty", default=None,
                   help="Path to flow_uncertainty.pt (optional). Adds the "
                   "FLOW (stage-1 template) uncertainty to the θ-vs-η plots: "
                   "the data-stat ±1σ (from --fisher) stays as the error bars / "
                   "inner band, and a wider ±1σ (data⊕flow) total band is "
                   "drawn behind it (combined in quadrature per component). The "
                   "per-component data/flow/total σ budget is also printed.")
    p.add_argument("--output", default=None,
                   help="Output directory (default: <checkpoint_dir>/diagnostics/).")
    p.add_argument("--device", default=("cuda:0" if torch.cuda.is_available() else "cpu"))
    p.add_argument("--batch-size", type=int, default=65536,
                   help="Loader batch size. Bigger is fine here — no "
                   "backward pass, only forward + collect to host.")
    p.add_argument("--n-mll-bins", type=int, default=50,
                   help="Number of m_ll bins. Cost of the per-event "
                   "grid eval scales linearly in this.")
    p.add_argument("--grid-chunk-events", type=int, default=4096,
                   help="Cap on events processed per flow call inside "
                   "the per-event × per-bin grid eval — bounds memory "
                   "(chunk_events × n_mll_bins expanded inputs per call).")
    p.add_argument("--max-events", type=int, default=0,
                   help="Stop after this many events (0 = run to end). "
                   "Bounds the grid-eval cost (O(N × n_mll_bins) flow "
                   "forwards) on full-statistics shards. Applies identically to "
                   "the main closure evaluation and the param-sensitivity pass.")
    p.add_argument("--split", default="holdout", choices=("train", "val", "holdout", "all"),
                   help="Which loader split to evaluate on. 'holdout' is the "
                   "untouched-by-training default and the canonical choice.")
    p.add_argument("--half", type=int, default=None, choices=(0, 1),
                   help="Validation deterministic event half (0=flow-train half, "
                   "1=fit half). When set, the loader uses ALL of that half "
                   "(split=train, val/holdout fraction=0), matching the fit's "
                   "stage-2 sample — for consistent in-sample (half 1) vs "
                   "out-of-sample (half 0) closure. Default: no half split.")
    p.add_argument("--eval-seed", type=int, default=42,
                   help="Fixed seed for the smearing-kernel ε on MC events. "
                   "Keeps the predicted-signal histogram deterministic across "
                   "diagnostic runs.")
    p.add_argument("--continuity-n-iter", type=int, default=2,
                   help="Fixed-point iterations for the #2 source solve.")
    p.add_argument("--param-shift", type=float, default=3.0,
                   help="Representative shift (in units of each parameter's "
                   "reference scale: THETA_SCALE_REF for A/e/M, SMEAR_VAR_SCALE "
                   "for a/c) used for the parameter-sensitivity overlay curves. "
                   "Default 3.0 → A±3e-4, e±3e-3, M±3e-5, a±3e-7, c±6e-5 — large "
                   "enough that the (otherwise sub-MeV) A/e/M peak shifts are "
                   "visible. The reference scales already roughly equalise the "
                   "inclusive m_ll effect across A/e/M (REF ratio ≈ k̄), so a "
                   "single global factor keeps them comparable.")
    p.add_argument("--no-param-sensitivity", action="store_true",
                   help="Skip the parameter-sensitivity slice plots "
                   "(param_sensitivity_*). They re-iterate the loader and do "
                   "~2·n_fitted extra grid evaluations per event, which is the "
                   "slowest diagnostic for the qop smear operator — pair with "
                   "--max-events for a quick look.")
    p.add_argument("--theta-scan-points", type=int, default=25,
                   help="Number of scan points for the single-η-bin θ_scale "
                   "likelihood scan (theta_scale_likelihood_scan.png). The scan "
                   "is auto-enabled only for binned θ with a single η bin "
                   "(--n-eta-bins 1); each point re-evaluates the data NLL + its "
                   "gradient over the FULL fitted sample (a dedicated loader "
                   "matching the fit — independent of --split and --max-events — "
                   "so the scan diagnoses the actual fit minimum).")
    p.add_argument("--theta-scan-nsigma", type=float, default=4.0,
                   help="Half-width of the θ_scale likelihood scan, in Fisher σ "
                   "(falls back to ±10 raw units when no covariance is available).")
    p.add_argument("--no-theta-scan", action="store_true",
                   help="Skip the single-η-bin θ_scale likelihood scan entirely "
                   "(the loader build + all passes). This is the slowest "
                   "diagnostic — it runs n_points forward+backward passes over the "
                   "FULL fitted sample — so use this to skip it when iterating.")
    return p.parse_args(argv)


def _inject_modulation_torch(eta, phi):
    """Per-muon non-uniform injection factor f(η,φ) (torch twin of the loader's
    ``_inject_modulation_np``) — quadratic-in-η × sinusoidal-in-φ, ~±50%."""
    u = eta / _INJECT_ETA_REF
    f_eta = 1.0 + _INJECT_AMP * (2.0 * u * u - 2.0 / 3.0)
    f_phi = 1.0 + _INJECT_AMP * torch.sin(_INJECT_PHI_NOSC * phi)
    return f_eta * f_phi


@contextlib.contextmanager
def _set_theta_output(model, scale_phys=None, smear_phys=None, nonuniform=False):
    """Temporarily SET the model's per-muon parameter output to fixed per-η-bin
    PHYSICAL values (masked to the fitted terms), for BOTH binned and MLP θ —
    used for the 'flow at injected θ' closure-target curve. Replaces
    ``_override_theta``, which sets the binned tensors and is a NO-OP under
    ``--theta-mlp`` (so that overlay was previously meaningless for MLP fits).
    ``scale_phys`` is [n_eta, 3] = (A, e, M); ``smear_phys`` is [n_eta, 2] =
    (a, c). The smear is converted to the O(1) coefficient (÷ SMEAR_VAR_SCALE)
    since ``_smear_ac_pm`` returns O(1) ac. The model's masks are applied so
    non-fitted terms are exactly 0 (matching the fit's parameterisation).

    ``nonuniform``: multiply the per-muon output by the (η,φ) modulation factor
    so the closure target matches a non-uniform injection (--inject-nonuniform).
    """
    dev = model.m_lo.device

    def _fmod(e, p):
        return _inject_modulation_torch(e, p).unsqueeze(-1) if nonuniform else 1.0

    orig_s, orig_c = model._scale_AeM_pm, model._smear_ac_pm
    if scale_phys is not None:
        S = torch.as_tensor(scale_phys, dtype=torch.float32, device=dev)  # [n_eta,3]
        model._scale_AeM_pm = lambda e, p, b: S[b] * model.scale_param_mask * _fmod(e, p)
    if smear_phys is not None:
        sc = torch.tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C],
                          dtype=torch.float32, device=dev)
        C = torch.as_tensor(smear_phys, dtype=torch.float32, device=dev) / sc  # O(1)
        model._smear_ac_pm = lambda e, p, b: C[b] * model.smear_param_mask * _fmod(e, p)
    try:
        yield
    finally:
        if scale_phys is not None:
            del model.__dict__["_scale_AeM_pm"]
        if smear_phys is not None:
            del model.__dict__["_smear_ac_pm"]


@contextlib.contextmanager
def _shift_theta_output(model, dscale=None, dsmear=None):
    """Temporarily ADD a physical (A,e,M)/(a,c) shift to the model's per-muon
    parameter output. Works for BOTH binned and MLP θ because it wraps the
    output methods (``_scale_AeM_pm`` / ``_smear_ac_pm``), shifting the effective
    θ however it is produced — unlike ``_override_theta``, which only sets the
    binned tensors and is a no-op under ``--theta-mlp``. ``dscale`` is a physical
    [A, e, M] shift; ``dsmear`` a physical [a, c] shift (converted to the O(1)
    coefficient via SMEAR_VAR_SCALE before being added to the O(1) ac output)."""
    dev = model.m_lo.device
    set_s, set_c = dscale is not None, dsmear is not None
    orig_s, orig_c = model._scale_AeM_pm, model._smear_ac_pm
    if set_s:
        ds = torch.as_tensor(dscale, dtype=torch.float32, device=dev)
        model._scale_AeM_pm = lambda e, p, b: orig_s(e, p, b) + ds
    if set_c:
        sc = torch.tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C],
                          dtype=torch.float32, device=dev)
        dc = torch.as_tensor(dsmear, dtype=torch.float32, device=dev) / sc
        model._smear_ac_pm = lambda e, p, b: orig_c(e, p, b) + dc
    try:
        yield
    finally:
        if set_s:
            del model.__dict__["_scale_AeM_pm"]
        if set_c:
            del model.__dict__["_smear_ac_pm"]


@torch.no_grad()
def plot_param_sensitivity(model, loader, stats, m_centers, out_dir, *,
                           shift_scale=1.0, max_events=0, chunk_events=2048,
                           n_iter=2, device="cpu", mc_as_data=False, progress=True):
    """Per-slice m_ll closure with overlays of the model signal density at ±a
    representative shift of each FITTED parameter — to expose parameter
    sensitivities and degeneracies. The data (pseudo-data in validation) is
    histogrammed per slice; the model density at the fitted θ and at θ_fit ± Δ
    (one parameter at a time) are slice-weighted-averaged and overlaid.

    Slice variables are CONDITIONAL quantities (or functions of them) — the flow
    conditions on muon_kin = (η_±, φ_±, ρ), so slicing in those keeps the
    conditioning fixed within a slice and the closure cleanly isolates the model
    behaviour. (Slicing in a NON-conditional variable such as pt_avg — the
    absolute pt scale, which the leak-free conditioning deliberately omits —
    would conflate the flow's pt-marginalisation with the θ-effects.) The
    discriminating slices:
      |η|max — reference (the existing closure axis), from η_±;
      ρ = (pt₊−pt₋)/(pt₊+pt₋) — charge-odd → isolates M (mass effect ∝ pt₊−pt₋);
      cos α = (cosΔφ + sinhη₊ sinhη₋)/(coshη₊ coshη₋) — the 3-D opening angle.
              The pt CANCELS, so it is a pure function of (η_±, Δφ) ⊂ muon_kin;
              and at fixed m, m² = 2p₊p₋(1−cosα) ⇒ α fixes p₊p₋, i.e. the pt
              scale (k̄) — so it separates A vs e (peak shift ∝ A−e·k̄) and
              a vs c (mass-variance ∝ a·pt²+c) using ONLY conditional info.

    Slice edges are adaptive (tertiles) per variable. The shift Δ per parameter
    is ``shift_scale`` × its reference scale (THETA_SCALE_REF for A/e/M,
    SMEAR_VAR_SCALE for a/c)."""
    G = m_centers.shape[0]
    mdev = m_centers.to(device)
    mc = m_centers.cpu().numpy()
    dmb = float(m_centers[1] - m_centers[0])
    m_edges = np.concatenate([[mc[0] - dmb / 2], mc[:-1] + dmb / 2, [mc[-1] + dmb / 2]])

    # ± representative shift per FITTED parameter.
    ref_s = THETA_SCALE_REF
    ref_c = (SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C)
    pal = {"A": "C0", "e": "C1", "M": "C2", "a": "C5", "c": "C6"}
    shifts = []  # (label, color, '+'/'-', dscale|None, dsmear|None)
    for j, nm in enumerate("AeM"):
        if model.scale_enabled and nm in model.scale_fit_params:
            for s in (1, -1):
                d = [0.0, 0.0, 0.0]; d[j] = s * shift_scale * ref_s[j]
                shifts.append((nm, pal[nm], "+" if s > 0 else "-", d, None))
    for j, nm in enumerate("ac"):
        if model.smearing_enabled and model.smear_fit_params in ("both", nm):
            for s in (1, -1):
                d = [0.0, 0.0]; d[j] = s * shift_scale * ref_c[j]
                shifts.append((nm, pal[nm], "+" if s > 0 else "-", None, d))
    if not shifts:
        print("  no fitted parameters → skipping param-sensitivity plots")
        return

    def _cos_alpha(b):
        e, p = b["eta_pm"], b["phi_pm"]
        num = torch.cos(p[:, 0] - p[:, 1]) + torch.sinh(e[:, 0]) * torch.sinh(e[:, 1])
        return num / (torch.cosh(e[:, 0]) * torch.cosh(e[:, 1]))

    # Slice variables, chosen for the active basis so each is a function of the
    # CONDITIONING (never the observable m_ll). For event_level, cos α / |η₊| / ρ
    # depend on the lab→CS boost (hence m_ll) and are NOT conditional, so use the
    # conditioning components themselves: p_T^ll/m_ll (the dimensionless ratio in
    # the basis) discriminates the A/e and a/c pt-scale degeneracies; cos θ*
    # (CS polar angle) is the M-sensitive analog of ρ (it controls the
    # charge-odd pt sharing).
    if getattr(model, "cond_basis", "muon_kin") == "event_level":
        def _ev(b):
            return _event_cond_raw(b["pt_pm"], b["eta_pm"], b["phi_pm"])
        slice_vars = {
            "ptll_over_mll": (lambda b: torch.exp(_ev(b)[:, 1]),
                              "p_T^ll/m_ll"),
            "yll": (lambda b: _ev(b)[:, 0], "y_ll"),
            "costhetastar": (lambda b: _ev(b)[:, 4], "cos θ*"),
        }
    else:
        slice_vars = {
            "abseta": (lambda b: b["eta_pm"].abs().max(1).values, "|η|max"),
            "rho": (lambda b: (b["pt_pm"][:, 0] - b["pt_pm"][:, 1])
                    / (b["pt_pm"][:, 0] + b["pt_pm"][:, 1]), "ρ"),
            "cosalpha": (_cos_alpha, "cos α (opening angle)"),
        }

    def _select(b):
        s = (~b["is_data_mask"]) if mc_as_data else b["is_data_mask"]
        return s.nonzero(as_tuple=True)[0]

    # Pre-pass: collect the slice-variable values to set adaptive (tertile) edges
    # — cheap (no flow evals), and robust to each variable's unknown range.
    vals = {v: [] for v in slice_vars}
    seen = 0
    for batch in loader:
        if max_events > 0 and seen >= max_events:
            break
        batch = _move_batch(batch, device)
        # Count ALL events per batch (as the main closure evaluate_predictions
        # does), so the same --max-events selects the same events in both passes.
        seen += int(batch["mll"].shape[0])
        idx = _select(batch)
        if idx.numel() == 0:
            continue
        for v, (fn, _) in slice_vars.items():
            vals[v].append(fn(batch)[idx].cpu().numpy())
    # `seen` now counts all events (for the cap); guard on whether any SELECTED
    # events were collected for the tertile edges.
    if all(len(lst) == 0 for lst in vals.values()):
        print("  no events selected → skipping param-sensitivity plots")
        return
    edges = {}
    for v in slice_vars:
        a = np.concatenate(vals[v])
        e = np.quantile(a, [0.0, 1.0 / 3, 2.0 / 3, 1.0])
        e[0] -= 1e-6; e[-1] += 1e-6   # include the extremes
        edges[v] = e
    del vals

    slice_specs = {v: (slice_vars[v][0], edges[v], slice_vars[v][1])
                   for v in slice_vars}
    acc = {v: [{"data": np.zeros(G), "fit": np.zeros(G), "w": 0.0,
                "sh": {i: np.zeros(G) for i in range(len(shifts))}}
               for _ in range(len(edges[v]) - 1)]
           for v in slice_vars}

    seen = 0
    bar = tqdm(loader, desc="param-sens", disable=not progress, unit="batch")
    for batch in bar:
        if max_events > 0 and seen >= max_events:
            break
        batch = _move_batch(batch, device)
        # Count ALL events per batch (matching the main closure), so the same
        # --max-events selects the same events here as in evaluate_predictions.
        seen += int(batch["mll"].shape[0])
        sel = (~batch["is_data_mask"]) if mc_as_data else batch["is_data_mask"]
        idx = sel.nonzero(as_tuple=True)[0]
        if idx.numel() == 0:
            continue
        mll = batch["mll"][idx].cpu().numpy()
        w = batch["w"][idx].cpu().numpy()
        p_fit = torch.exp(_tilt_density_on_grid(
            model, batch, idx, mdev, chunk_events=chunk_events, n_iter=n_iter)).cpu().numpy()
        p_sh = []
        for _, _, _, ds, dc in shifts:
            with _shift_theta_output(model, ds, dc):
                p_sh.append(torch.exp(_tilt_density_on_grid(
                    model, batch, idx, mdev, chunk_events=chunk_events,
                    n_iter=n_iter)).cpu().numpy())
        for v, (fn, edges, _) in slice_specs.items():
            val = fn(batch)[idx].cpu().numpy()
            for si in range(len(edges) - 1):
                sm = (val >= edges[si]) & (val < edges[si + 1])
                if not sm.any():
                    continue
                a = acc[v][si]
                a["data"] += np.histogram(mll[sm], bins=m_edges, weights=w[sm])[0]
                a["fit"] += (p_fit[sm] * w[sm, None]).sum(0)
                a["w"] += float(w[sm].sum())
                for i in range(len(shifts)):
                    a["sh"][i] += (p_sh[i][sm] * w[sm, None]).sum(0)
    bar.close()

    for v, (fn, edges, xlabel) in slice_specs.items():
        ns = len(edges) - 1
        # Two rows per slice: density (top) + ratio-to-fit (bottom). The ratio
        # panel makes the (small) shift sensitivity and the data/fit closure
        # legible in fractional terms — shift/fit isolates each parameter's
        # effect and data/fit shows the residual mis-closure.
        fig, axes = plt.subplots(
            2, ns, figsize=(4.8 * ns, 4.6), squeeze=False, sharex="col",
            gridspec_kw={"height_ratios": [3, 1]})
        for si in range(ns):
            ax, axr, a = axes[0, si], axes[1, si], acc[v][si]
            if a["w"] <= 0:
                ax.set_visible(False); axr.set_visible(False); continue
            fit_cnt = a["fit"] * dmb                  # model expected counts
            ax.step(mc, a["data"], where="mid", color="k", lw=1.2,
                    label="pseudo-data" if mc_as_data else "data")
            ax.plot(mc, fit_cnt, color="0.3", lw=2.0, label="model (fit)")
            for i, (nm, col, sgn, _, _) in enumerate(shifts):
                ax.plot(mc, a["sh"][i] * dmb, color=col, lw=1.0,
                        ls="--" if sgn == "+" else ":", alpha=0.85,
                        label=f"{nm}{sgn}Δ")
            hi = "∞" if edges[si + 1] > 1e8 else f"{edges[si + 1]:.2f}"
            ax.set_title(f"{xlabel} ∈ [{edges[si]:.2f}, {hi})", fontsize=9)
            ax.grid(alpha=0.3)
            if si == 0:
                ax.legend(fontsize=6, ncol=2)
                axr.set_ylabel("ratio / fit", fontsize=8)
            # --- ratio panel: everything ÷ the fitted model ---
            fok = fit_cnt > 0
            r_data = np.divide(a["data"], fit_cnt,
                               out=np.full_like(fit_cnt, np.nan), where=fok)
            axr.axhline(1.0, color="0.3", lw=1.2)
            axr.step(mc, r_data, where="mid", color="k", lw=1.0)
            r_all = [r_data]
            for i, (nm, col, sgn, _, _) in enumerate(shifts):
                r_sh = np.divide(a["sh"][i], a["fit"],
                                 out=np.full_like(a["fit"], np.nan), where=fok)
                axr.plot(mc, r_sh, color=col, lw=1.0,
                         ls="--" if sgn == "+" else ":", alpha=0.85)
                r_all.append(r_sh)
            # Robust y-range over the data-supported region (ignore empty tails).
            sup = a["data"] > 0.02 * (a["data"].max() if a["data"].max() > 0 else 1.0)
            rr = np.concatenate([r[sup] for r in r_all]) if sup.any() else np.array([])
            rr = rr[np.isfinite(rr)]
            if rr.size:
                lo, hi_r = np.percentile(rr, [1, 99])
                pad = 0.1 * max(hi_r - lo, 1e-3)
                axr.set_ylim(max(0.0, lo - pad), hi_r + pad)
            axr.set_xlabel("m_ll [GeV]"); axr.grid(alpha=0.3)
        fig.suptitle(f"parameter sensitivity — slices of {xlabel} "
                     f"(Δ = {shift_scale:g}× ref scale)")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(out_dir, f"param_sensitivity_{v}.{ext}"), dpi=110)
        plt.close(fig)
        print(f"  wrote param_sensitivity_{v} ({ns} slices)")


def plot_theta_scale_likelihood_scan(
    model, loader, device, output_dir, *,
    scale_fit_params: str,
    mc_as_data: bool,
    n_iter: int = 2,
    max_events: int = 0,
    sigma_scale: "np.ndarray | None" = None,
    inject_ref: "np.ndarray | None" = None,
    n_points: int = 25,
    n_sigma: float = 4.0,
    progress: bool = True,
) -> None:
    """1-D NLL likelihood scan (+ gradient) over each active θ_scale component.

    Intended for the SINGLE-η-bin binned-θ fit (``--n-eta-bins 1``), where the
    three scale nuisances (A, e, M) are global numbers and a direct scan of the
    objective is both cheap and the most transparent uncertainty diagnostic.

    Caches the (pseudo-)data once, then for every fit parameter listed in
    ``scale_fit_params`` sweeps that component over a window — ±``n_sigma``·σ from
    the Fisher covariance when available, otherwise a default ±10 raw units —
    holding the others at the fit value, and re-accumulates BOTH the weighted
    data NLL and its gradient ∂NLL/∂θ (``data_nll_continuity`` + autograd, exactly
    the stage-2 objective and its gradient). The NLL scan, the Fisher parabola,
    and the gradient (twin axis) are overlaid with the fit value and injected
    truth, directly exposing curvature, non-parabolicity, bias, AND convergence:
    a non-zero ∂NLL/∂θ at the fit marker (the gradient zero-crossing displaced
    from the fit line) is the under-convergence signature.

    IMPORTANT: ``loader`` must yield the SAME events the fit used (the FULL
    fitted sample) and ``max_events`` should be 0. The objective minimum is a
    property of the fitted sample; the shards are inhomogeneous, so a subset
    (a --split slice, or a --max-events truncation) has its OWN minimum that can
    sit anywhere — diagnosing it would mislocate the fit minimum and fake or hide
    a bias. The caller builds a dedicated full-sample loader for exactly this.

    The scan runs the model in **float64**: the per-event NLL is otherwise
    computed in fp32 (the training --precision), and summed over the full
    (often millions of events) dataset the fp32 per-term round-off (~N·ε) shifts
    unpredictably between adjacent scan points, producing a spuriously spiky
    curve. fp64 removes that round-off so the curve reflects the true objective.
    """
    comp_index = {"A": 0, "e": 1, "M": 2}
    comp_label = {"A": "A", "e": "e [GeV]", "M": "M"}
    active = [c for c in ("A", "e", "M") if c in (scale_fit_params or "")]
    if not active:
        print("  θ_scale scan: no active scale fit params; skipping.")
        return
    ref_phys = np.asarray(THETA_SCALE_REF, dtype=np.float64)   # (1e-4, 1e-3, 1e-5)

    # Cache the (pseudo-)data the fit used (one host→device pass), keeping only
    # the tensors the continuity NLL needs. Float tensors are upcast to float64
    # to match the float64 model forward (see the docstring on spikiness).
    cache = []
    seen = 0
    keys = ("mll", "pt_pm", "eta_pm", "phi_pm", "q_pm", "b_pm", "cond_std", "w")
    for batch in tqdm(loader, desc="θ_scale scan: caching events", unit="batch",
                      leave=False, disable=not progress):
        if max_events > 0 and seen >= max_events:
            break
        batch = _move_batch(batch, device)
        dm = (~batch["is_data_mask"]) if mc_as_data else batch["is_data_mask"]
        seen += int(batch["mll"].shape[0])
        if not bool(dm.any()):
            continue
        entry = {k: (batch[k].double() if batch[k].is_floating_point() else batch[k])
                 for k in keys}
        entry["dm"] = dm
        cache.append(entry)
    if not cache:
        print("  θ_scale scan: no data rows found; skipping.")
        return
    print(f"  θ_scale scan: cached {len(cache)} batches ({seen} events); "
          f"{int(n_points)} scan pts × {len(active)} param(s) "
          f"(each pt = 1 fwd+bwd over the full sample at the OBSERVED mass; no grid)")

    # Run the model forward in float64 for the scan to kill fp32 round-off noise
    # (restored to the original dtype at the end).
    orig_dtype = next(model.parameters()).dtype
    model.double()

    def nll_and_grad(j: int, bar=None):
        """Total weighted NLL and ∂NLL/∂θ_scale[0, j] (raw units) summed over the
        cached (pseudo-)data — the exact stage-2 objective + gradient (autograd).
        Both come from a single forward+backward per call. Each event is evaluated
        ONLY at its observed mass (data_nll_continuity) — no mass grid."""
        nll_acc = 0.0
        g_acc = 0.0
        for b in cache:
            per = model.data_nll_continuity(
                b["mll"], b["pt_pm"], b["eta_pm"], b["phi_pm"], b["q_pm"],
                b["b_pm"], b["cond_std"], b["dm"], n_iter=n_iter)
            w = b["w"] * b["dm"].to(b["w"].dtype)
            loss = (w * per).sum()
            g, = torch.autograd.grad(loss, model.theta_scale)
            nll_acc += float(loss.item())
            g_acc += float(g[0, j].item())
            if bar is not None:
                bar.update(1)
        return nll_acc, g_acc

    fit_raw = model.theta_scale.detach().clone()    # [1, 3] raw O(1) params

    results = []
    for c in active:
        j = comp_index[c]
        center_raw = float(fit_raw[0, j].item())
        sig_phys = None
        if sigma_scale is not None:
            s = float(sigma_scale[0, j])
            if np.isfinite(s) and s > 0:
                sig_phys = s
        if sig_phys is not None:
            sig_raw = sig_phys / ref_phys[j]
            hw_raw = n_sigma * sig_raw
        else:
            sig_raw = None
            hw_raw = 10.0   # default raw half-window when no Fisher σ available
        # Always keep the injected truth in frame (with 30% margin): the Fisher σ
        # can be far too tight when the likelihood is flat, which would otherwise
        # push the injected value — and the closure penalty ΔNLL(injected) — off
        # the edge of the plot.
        inj_raw = (None if inject_ref is None else float(inject_ref[0, j]) / ref_phys[j])
        if inj_raw is not None:
            hw_raw = max(hw_raw, 1.3 * abs(inj_raw - center_raw))
        xs_raw = np.linspace(center_raw - hw_raw, center_raw + hw_raw, int(n_points))
        nlls = np.empty(xs_raw.shape[0], dtype=np.float64)
        grads = np.empty(xs_raw.shape[0], dtype=np.float64)   # ∂NLL/∂θ_raw
        bar = tqdm(total=int(n_points) * len(cache), unit="batch",
                   desc=f"θ_scale scan [{c}] ({len(cache)} batches × {int(n_points)} pts)",
                   leave=False, disable=not progress)
        for k, xr in enumerate(xs_raw):
            with torch.no_grad():
                model.theta_scale[0, j] = float(xr)
            nlls[k], grads[k] = nll_and_grad(j, bar)
            bar.set_postfix_str(f"A={xr * ref_phys[j]:.3e} ΔNLL={nlls[k]-nlls[:k+1].min():.2f}")
        bar.close()
        with torch.no_grad():
            model.theta_scale[0, j] = center_raw   # restore the fit value
        g_fit = float(np.interp(center_raw, xs_raw, grads))   # ∂NLL/∂θ at the fit
        dnll = nlls - float(np.min(nlls))
        parab = None if sig_raw is None else 0.5 * ((xs_raw - center_raw) / sig_raw) ** 2
        results.append(dict(
            comp=c, x_phys=xs_raw * ref_phys[j], dnll=dnll, grad=grads,
            fit_phys=center_raw * ref_phys[j], sig_phys=sig_phys, g_fit=g_fit,
            inj_phys=(float(inject_ref[0, j]) if inject_ref is not None else None),
            parab=parab))

    # Restore the model to its original (training) precision for the rest of the
    # diagnostics.
    if orig_dtype == torch.float32:
        model.float()

    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 3.2 * n), squeeze=False)
    axes = axes[:, 0]
    for ax, r in zip(axes, results):
        h_nll, = ax.plot(r["x_phys"], r["dnll"], "o-", color="k", ms=3, lw=1.2,
                         label="NLL scan")
        handles = [h_nll]
        if r["parab"] is not None:
            h_p, = ax.plot(r["x_phys"], r["parab"], color="C0", ls="-", lw=1.2,
                           label="Gaussian (Fisher)")
            handles.append(h_p)
        ax.axhline(0.5, color="0.6", lw=0.8, ls=":")   # 1σ crossing
        ax.set_ylabel("ΔNLL")
        ax.set_ylim(bottom=0.0)
        ax.grid(True, alpha=0.3)
        # Gradient ∂NLL/∂θ (raw fit-param units) on a twin axis — the quantity
        # the optimiser drives to zero; its zero-crossing is the true minimum, so
        # a non-zero value at the fit line = under-convergence.
        ax2 = ax.twinx()
        h_g, = ax2.plot(r["x_phys"], r["grad"], "-", color="C4", lw=1.3,
                        label="∂NLL/∂θ (raw)")
        ax2.axhline(0.0, color="C4", lw=0.7, ls=":")
        ax2.set_ylabel("∂NLL/∂θ  (raw θ units)", color="C4")
        ax2.tick_params(axis="y", labelcolor="C4")
        gmax = float(np.max(np.abs(r["grad"]))) or 1.0
        ax2.set_ylim(-1.15 * gmax, 1.15 * gmax)   # symmetric so the zero line is centred
        handles.append(h_g)
        # fit + injected markers (drawn on the NLL axis).
        handles.append(ax.axvline(r["fit_phys"], color="C2", lw=1.0,
                                   label=f"fit (∂NLL/∂θ|fit={r['g_fit']:+.2g})"))
        if r["inj_phys"] is not None:
            handles.append(ax.axvline(r["inj_phys"], color="C3", lw=1.2, ls="--",
                                      label="injected"))
        xlab = comp_label[r["comp"]]
        if r["sig_phys"] is not None:
            xlab += f"   (σ_Fisher = {r['sig_phys']:.2e})"
        ax.set_xlabel(xlab)
        ax.legend(handles=handles, loc="upper center", fontsize=7, ncol=2)
    axes[0].set_title("θ_scale likelihood scan (single η-bin)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(output_dir, f"theta_scale_likelihood_scan.{ext}"),
                    dpi=120)
    plt.close(fig)
    print(f"  wrote {os.path.join(output_dir, 'theta_scale_likelihood_scan.png')} (+ .pdf)")


def main() -> int:
    args = parse_args()
    out_dir = args.output or os.path.join(
        os.path.dirname(args.checkpoint), "diagnostics",
    )
    os.makedirs(out_dir, exist_ok=True)

    device = args.device
    if device.startswith("cuda") and not torch.cuda.is_available():
        print("CUDA requested but unavailable; using CPU.")
        device = "cpu"

    print(f"loading checkpoint: {args.checkpoint}")
    model, stats, train_args, ckpt = load_model_from_checkpoint(args.checkpoint, device)

    # Stage-1 FLOW checkpoint (no fit ran): the meaningful closure is the RAW
    # template vs NOMINAL MC at θ = 0. The checkpoint's args still carry the
    # full command line (--inject-A/... apply only to the STAGE-2 pseudo-data;
    # stage 1 trains on un-injected events), and the saved θ tables are at
    # their INIT values — which for --smear-param-form softplus is a NONZERO
    # effective smearing (softplus(0)·scale ≈ 0.69·scale). Force exact θ = 0
    # (scale AND smear, binned tables and θ-net alike) and skip the injection
    # replay below, so the closure plots show pure template fidelity.
    is_flow_ckpt = str(ckpt.get("stage", "")) == "flow"
    if is_flow_ckpt:
        smear_zero_raw = (-30.0   # softplus(−30) ≈ 1e−13 → effective smear ≈ 0
                          if getattr(model, "smear_param_form", "linear") == "softplus"
                          else 0.0)
        with torch.no_grad():
            model.theta_scale.zero_()
            model.theta_smear.fill_(smear_zero_raw)
            if getattr(model, "theta_net", None) is not None:
                last = model.theta_net.net[-1]
                last.weight.zero_()
                last.bias.zero_()
                last.bias[N_THETA_SCALE:].fill_(smear_zero_raw)
        print("stage-1 FLOW checkpoint detected: forcing θ = 0 (no scale "
              "shift, no smearing) and skipping the injection replay — the "
              "closure plots compare the raw template against NOMINAL MC")

    # Loader.
    shard_files = discover_shards([args.shards])
    if not shard_files:
        print(f"error: no .arrow shards found under {args.shards!r}", file=sys.stderr)
        return 1
    # Validation (MC-closure) checkpoints fit θ on simulation as pseudo-data;
    # auto-detect so the m_ll closure / pulls treat simulation as the data
    # branch (otherwise those plots are empty — there are no is_data rows).
    mc_as_data = bool(train_args.get("validation", False))
    if mc_as_data:
        print("checkpoint was trained with --validation: routing simulation "
              "through the data branch as pseudo-data for the m_ll closure / pulls")

    # Injected θ_scale shift (validation closure with a non-zero target): read
    # the values the fit was trained with, replay the same m_ll injection in the
    # pseudo-data, and use them as the χ² reference + dashed line on the θ plot.
    inject_np = None
    inject_smear_np = None
    if mc_as_data and not is_flow_ckpt:
        n_eta = len(stats.eta_edges) - 1
        ia = float(train_args.get("inject_A", 0.0) or 0.0)
        ie = float(train_args.get("inject_e", 0.0) or 0.0)
        im = float(train_args.get("inject_M", 0.0) or 0.0)
        if ia or ie or im:
            inject_np = np.zeros((n_eta, 3), dtype=np.float64)
            inject_np[:, 0] = ia; inject_np[:, 1] = ie; inject_np[:, 2] = im
            print(f"checkpoint injected θ_scale (A,e,M)=({ia:g},{ie:g},{im:g}) — "
                  f"replaying the pseudo-data injection; closure target = injected")
        isa = float(train_args.get("inject_a", 0.0) or 0.0)
        isc = float(train_args.get("inject_c", 0.0) or 0.0)
        if isa or isc:
            # --inject-a/-c are PHYSICAL σ²_qop coefficients — used directly.
            inject_smear_np = np.zeros((n_eta, 2), dtype=np.float64)
            inject_smear_np[:, 0] = isa
            inject_smear_np[:, 1] = isc
            print(f"checkpoint injected smear (a,c)=({isa:g},{isc:g}) — replaying "
                  f"the per-muon qop fold into the pseudo-data")

    inject_bkg_np = None
    if mc_as_data and not is_flow_ckpt:
        bf0 = float(train_args.get("inject_bkg_f0", 0.0) or 0.0)
        bf1 = float(train_args.get("inject_bkg_f1", 0.0) or 0.0)
        if bf0 > 0.0 or bf1 > 0.0:
            inject_bkg_np = (bf0, bf1)
            print(f"checkpoint injected Bernstein background (f0,f1)=({bf0:g},"
                  f"{bf1:g}) — replaying into the pseudo-data; the fitted MLP "
                  f"f(c) should recover these (see bkg_fraction closure plots)")

    inject_prod_np = None
    if mc_as_data and not is_flow_ckpt:
        sp = (float(train_args.get("inject_prod_ptll_slope", 0.0) or 0.0),
              float(train_args.get("inject_prod_yll_slope", 0.0) or 0.0),
              float(train_args.get("inject_prod_costheta_slope", 0.0) or 0.0))
        if any(v != 0.0 for v in sp):
            inject_prod_np = sp
            print(f"checkpoint injected production/decay BIAS slopes "
                  f"(ln ptll, yll, cosθ*)=({sp[0]:g},{sp[1]:g},{sp[2]:g}) — "
                  f"replaying the signal-MC reweight into the pseudo-data so "
                  f"the closure plots show the biased data")

    nonuniform = bool(train_args.get("inject_nonuniform", False))
    if nonuniform and (inject_np is not None or inject_smear_np is not None):
        print("  injection is NON-UNIFORM (quadratic-η × sinusoidal-φ); the "
              "θ-vs-η reference uses the φ-averaged truth base·f_η(η)")
    # --half: match the validation fit's stage-2 sample (ALL of the half, no
    # train/val/holdout sub-split) so the closure is consistent in-sample
    # (half 1) / out-of-sample (half 0). Otherwise use the requested split.
    _half = getattr(args, "half", None)
    _split = "train" if _half is not None else args.split
    _vf = 0.0 if _half is not None else float(train_args.get("val_fraction", 0.10))
    _hf = 0.0 if _half is not None else float(train_args.get("holdout_fraction", 0.05))
    print(f"found {len(shard_files)} shard(s); split={_split}"
          + (f"; half={_half} (all of the half, no sub-split)" if _half is not None else ""))
    loader = JpsiMassArrowLoader(
        shard_files, stats,
        batch_size=args.batch_size,
        split=_split, half=_half,
        val_fraction=_vf,
        holdout_fraction=_hf,
        drop_last=False,
        inject_theta_scale=inject_np,
        inject_theta_smear=inject_smear_np,
        inject_seed=int(train_args.get("inject_smear_seed", 12345)),
        cond_basis=train_args.get("cond_basis", "muon_kin"),
        inject_nonuniform=nonuniform,
        inject_bkg=inject_bkg_np,
        m_window=(model._m_lo_f, model._m_hi_f),
        inject_prod=inject_prod_np,
        reco_ptll_min=train_args.get("reco_ptll_min"),
        reco_ptll_max=train_args.get("reco_ptll_max"),
        fit_select=_fit_select_from(train_args),
    )
    # φ-AVERAGED injected reference for the θ-vs-η plots + χ²: the φ sinusoid
    # averages to 1 over the plotted φ-mean, leaving base·f_η(η) per η-bin
    # centre. (The loader keeps the un-modulated base; it applies f(η,φ) per muon
    # to the pseudo-data itself.)
    _eta_c = 0.5 * (np.asarray(stats.eta_edges[:-1]) + np.asarray(stats.eta_edges[1:]))
    _feta = _inject_modulation_eta_np(_eta_c) if nonuniform else np.ones_like(_eta_c)
    inject_ref_np = None if inject_np is None else inject_np * _feta[:, None]
    inject_smear_ref_np = (None if inject_smear_np is None
                           else inject_smear_np * _feta[:, None])

    def _ref_cells(ref_eta, n_phi):
        # Per-CELL injected reference for the 2-D binned θ grid
        # [n_eta, n_phi, ncol]: the per-η base·f_η times the φ-BIN-AVERAGED
        # sinusoid f_φ(φ) = 1 + amp·sin(Nφ) when --inject-nonuniform (f_φ ≡ 1
        # for a uniform injection). A per-cell θ absorbs the cell AVERAGE of
        # the modulation, so the bin-centre sinusoid is damped by
        # sinc(N/n_phi) = sin(Nw/2)/(Nw/2), w = 2π/n_phi (11% for N=2, 8 bins).
        phi_ctr = (np.arange(n_phi) + 0.5) * (2.0 * np.pi / n_phi) - np.pi
        fphi = (1.0 + _INJECT_AMP * np.sinc(_INJECT_PHI_NOSC / n_phi)
                * np.sin(_INJECT_PHI_NOSC * phi_ctr)
                if nonuniform else np.ones_like(phi_ctr))
        return np.asarray(ref_eta)[:, None, :] * fphi[None, :, None]

    # m_ll grid.
    m_edges = torch.linspace(model._m_lo_f, model._m_hi_f, args.n_mll_bins + 1)
    bin_width = float((m_edges[1] - m_edges[0]).item())
    m_centers = 0.5 * (m_edges[:-1] + m_edges[1:])
    m_centers_np = m_centers.cpu().numpy()
    m_grid_std = (m_centers - stats.mll_mean) / stats.mll_std

    # η slices for the per-slice closure plots.
    eta_slice_edges = np.array([0.0, 0.6, 1.2, 1.8, 2.4])

    print(
        f"evaluating model on the loader (batch_size={args.batch_size}, "
        f"n_mll_bins={args.n_mll_bins}, "
        f"grid_chunk_events={args.grid_chunk_events}"
        f"{', max_events=' + str(args.max_events) if args.max_events else ''})..."
    )
    evals = evaluate_predictions(
        model, loader, device, m_centers, m_grid_std, bin_width,
        chunk_events=args.grid_chunk_events,
        max_events=args.max_events,
        progress=True,
        seed=args.eval_seed,
        n_iter=args.continuity_n_iter,
        mc_as_data=mc_as_data,
    )
    if mc_as_data:
        print(
            f"  collected {evals['mll_data'].shape[0]} MC pseudo-data events "
            f"(same simulation also drives the MC-branch closure)"
        )
    else:
        print(
            f"  collected {evals['mll_data'].shape[0]} data events, "
            f"{evals['mll_mc_fold'].shape[0]} MC events"
        )

    # Plot 1: m_ll closure.
    print("plotting m_ll closure...")
    plot_mll_closure(
        evals, m_centers_np, eta_slice_edges,
        model._m_lo_f, model._m_hi_f, out_dir,
    )

    # Fisher info → ±1σ for θ_scale (and θ_smear, when present).
    sigma_scale = None
    sigma_smear = None
    cov_scale_flat = None   # 72×72 θ_scale covariance block (for the χ² test)
    cov_smear_flat = None    # 48×48 θ_smear (a,c) covariance (for the whitened band)
    cov_scale_2d = cov_smear_2d = None   # full (η,φ) output cov (φ-resolved band)
    fisher_cell_w = None     # [n_eta,n_phi] Σw occupancy (sample-weighted avgs)
    edm = None
    mlp_fisher = False   # the file carries θ-net-weight Fisher propagated to per-η σ
    if args.fisher and os.path.exists(args.fisher):
        f = torch.load(args.fisher, weights_only=False)
        mlp_fisher = (f.get("theta_mode") == "mlp")
        if f.get("n_negative_eig", 0) not in (0, None):
            print(f"  warning: fisher_info.pt has {f['n_negative_eig']} "
                  f"non-positive eigenvalue(s) (min={f.get('min_eig', float('nan')):.2e}) "
                  f"— σ bands for affected params are unreliable.")
        edm = f.get("edm")
        if edm is not None:
            print(f"  fit EDM (½ gᵀV g) = {edm:.3e}")
        cov_pt = f.get("covariance_24_3_24_3")
        if cov_pt is not None:
            # generic η-bin count from the saved [n_eta,3,n_eta,3] tensor (the
            # "24_3" key name is historical; n_eta is set by --n-eta-bins).
            n_eta_s = int(cov_pt.shape[0])
            cov_scale_flat = cov_pt.reshape(3 * n_eta_s, 3 * n_eta_s).cpu().numpy()
            sigma_scale = np.sqrt(np.maximum(np.diag(cov_scale_flat), 0.0)).reshape(n_eta_s, 3)
        ss = f.get("sigma_smear_eff_24_2")
        if ss is not None:
            sigma_smear = ss.cpu().numpy()
        cs_pt = f.get("covariance_smear_24_2_24_2")
        if cs_pt is not None:
            n_eta_c = int(cs_pt.shape[0])
            cov_smear_flat = cs_pt.reshape(2 * n_eta_c, 2 * n_eta_c).cpu().numpy()
        # Full 2-D (η,φ) output covariance + occupancy weights (output-fisher),
        # for the φ-RESOLVED Fisher band and the SAMPLE-weighted η-average on the
        # θ-vs-φ plots. Present only for the output-space Fisher.
        cov_scale_2d = (f["covariance_scale_2d"].cpu().numpy()
                        if f.get("covariance_scale_2d") is not None else None)
        cov_smear_2d = (f["covariance_smear_2d"].cpu().numpy()
                        if f.get("covariance_smear_2d") is not None else None)
        fisher_cell_w = (f["cell_w"].cpu().numpy()
                         if f.get("cell_w") is not None else None)
        # Full joint covariance (θ_scale + active θ_smear); the matrix plot is
        # drawn later from the TOTAL (data⊕flow) when --flow-uncertainty is set.
        full_cov = f.get("covariance")
        full_cov_labels = f.get("labels")
        full_cov_nscale = int(f.get("n_scale", 72))
    else:
        full_cov = full_cov_labels = None
        full_cov_nscale = 72
        print("no fisher_info.pt → skipping θ_scale ±σ bands + correlation plot.")

    # Flow (stage-1 template) uncertainty: ADD it to the data-stat covariance to
    # get the TOTAL (data⊕flow), used everywhere downstream — the θ-vs-η error
    # bars (binned) / band (mlp), the χ² compatibility test, the covariance /
    # correlation matrices, and the whitened bands. If --fisher was not given,
    # the "total" is the flow contribution alone. (Per-block keys add exactly:
    # cov_total = cov_data + cov_flow.) sigma_*_total = √diag(cov_total).
    sigma_scale_data = sigma_scale          # data-stat σ (inner error bar)
    sigma_smear_data = sigma_smear
    sigma_scale_total = sigma_smear_total = None
    cov_scale_total = cov_scale_flat        # default = data-stat (no flow)
    cov_smear_total = cov_smear_flat
    full_cov_total = full_cov
    fu_path = getattr(args, "flow_uncertainty", None)
    if fu_path and os.path.exists(fu_path):
        fu = torch.load(fu_path, weights_only=False)
        print(f"flow uncertainty from {os.path.basename(fu_path)} "
              f"(w-solver={fu.get('w_solver')}, flow-solver={fu.get('flow_solver')}):")

        def _block(cov_data_ext, key, ncol):
            """Resolve (cov_data, cov_total) [flat n,n] for a block. Prefer the
            file's OWN <key>_data / <key>_total (the --w-solver fisher run stores
            them); else combine the file's <key>_flow with the --fisher data."""
            def _flat(t):
                if t is None:
                    return None
                a = t.cpu().numpy()
                ne = a.shape[0]
                return a.reshape(ne * ncol, ne * ncol)
            f_tot, f_dat = _flat(fu.get(f"{key}_total")), _flat(fu.get(f"{key}_data"))
            if f_tot is not None:                     # self-sufficient file
                return (f_dat if f_dat is not None else cov_data_ext), f_tot
            f_flow = _flat(fu.get(f"{key}_flow"))
            if f_flow is None:
                return cov_data_ext, None
            tot = (cov_data_ext + f_flow) if (cov_data_ext is not None
                   and cov_data_ext.shape == f_flow.shape) else f_flow
            return cov_data_ext, tot

        d_s, cov_scale_total = _block(cov_scale_flat, "covariance_24_3_24_3", 3)
        if cov_scale_total is not None:
            ne = cov_scale_total.shape[0] // 3
            sigma_scale_total = np.sqrt(np.maximum(np.diag(cov_scale_total), 0.0)).reshape(ne, 3)
            if d_s is not None:                       # data-stat from the file
                sigma_scale_data = np.sqrt(np.maximum(np.diag(d_s), 0.0)).reshape(ne, 3)
                cov_scale_flat = d_s
        d_c, cov_smear_total = _block(cov_smear_flat, "covariance_smear_24_2_24_2", 2)
        if cov_smear_total is not None:
            ne = cov_smear_total.shape[0] // 2
            sigma_smear_total = np.sqrt(np.maximum(np.diag(cov_smear_total), 0.0)).reshape(ne, 2)
            if d_c is not None:
                sigma_smear_data = np.sqrt(np.maximum(np.diag(d_c), 0.0)).reshape(ne, 2)
                cov_smear_flat = d_c
        # Full joint covariance (raw active-θ layout): file's own total/data, or
        # data(--fisher) + flow.
        ftot = fu.get("covariance_total")
        if ftot is not None:
            full_cov_total = ftot.double()
            full_cov_labels = full_cov_labels or fu.get("labels")
            full_cov_nscale = int(fu.get("n_scale", full_cov_nscale))
        else:
            fcf = fu.get("covariance_flow")
            if fcf is not None:
                if full_cov is not None and tuple(full_cov.shape) == tuple(fcf.shape):
                    full_cov_total = (full_cov.double() + fcf.double())
                else:
                    full_cov_total = fcf
                    full_cov_labels = full_cov_labels or fu.get("labels")
                    full_cov_nscale = int(fu.get("n_scale", full_cov_nscale))

        def _budget(tag, sig_data, total):
            d = (float(np.median(sig_data)) if sig_data is not None else float("nan"))
            t = (float(np.median(total)) if total is not None else float("nan"))
            fl = (float(np.median(np.sqrt(np.maximum(total ** 2 - (sig_data ** 2
                  if sig_data is not None else 0.0), 0.0))))
                  if total is not None else float("nan"))
            print(f"  {tag}: median σ — data {d:.3e} | flow {fl:.3e} | total {t:.3e}")
        _budget("scale(A,e,M)", sigma_scale_data, sigma_scale_total)
        _budget("smear(a,c)", sigma_smear_data, sigma_smear_total)
    elif fu_path:
        print(f"  warning: --flow-uncertainty {fu_path} not found; skipping the "
              f"total band.")

    # Covariance / correlation matrices (from the TOTAL when --flow-uncertainty
    # is set, else the data-stat Fisher).
    _cov_tag = "total (data⊕flow)" if (fu_path and sigma_scale_total is not None
                                       ) else "data-stat"
    if full_cov_total is not None:
        print(f"plotting covariance / correlation matrix ({_cov_tag})...")
        plot_cov_corr(full_cov_total.detach().cpu().numpy(), full_cov_labels,
                      full_cov_nscale, out_dir, edm=edm)
    elif cov_scale_total is not None:
        print(f"plotting correlation matrix (θ_scale block, {_cov_tag})...")
        plot_fisher_correlation(cov_scale_total, out_dir)
    else:
        print("  no covariance available; skipping matrix plot.")

    # Plots 2, 3: θ vs η — only for the *enabled* nuisances (a disabled one
    # is an inert, fixed parameter; plotting it would be misleading).
    print("plotting θ vs η...")
    # 'mlp' θ: sample the continuous ThetaNet at the η-bin centres (φ-mean + a
    # φ-spread band) so the per-bin plot shows the learned function. The
    # statistical ±1σ error bars come from --fisher empirical_fisher.pt when it
    # carries the θ-net-weight Fisher propagated to the per-η outputs (the φ-std
    # band is a separate, systematic-like φ-variation overlay).
    mlp_scale_grid = mlp_smear_grid = None
    mlp_scale_band = mlp_smear_band = None
    mlp_scale_slices = mlp_smear_slices = None
    mlp_scale_avg_samples = mlp_smear_avg_samples = None
    mlp_slice_labels = None
    mlp_phi = None
    if model.theta_mode == "mlp":
        # Sample the ThetaNet on a 2D (η-centre, φ) grid: 16 points uniformly
        # spaced over [0, 2π) for the φ-average and ±std band, plus 4 cardinal
        # slices {0, π/2, π, −π/2} for the overlaid curves. Both muons are
        # given the same (η, φ) since the ThetaNet is per-muon (the symmetry
        # is implicit). The conditioner sees (cos φ, sin φ), so uniform φ on
        # the circle gives an exact mean + std with no statistical noise.
        centers = 0.5 * (np.asarray(stats.eta_edges[:-1]) + np.asarray(stats.eta_edges[1:]))
        centers_t = torch.as_tensor(centers, dtype=torch.float32, device=device)
        n_eta_c = int(centers_t.shape[0])
        # φ-average sampling (16 uniform points → exact integral over the circle).
        n_phi_avg = 16
        phi_avg = torch.linspace(
            0.0, 2 * np.pi * (1.0 - 1.0 / n_phi_avg), n_phi_avg,
            dtype=torch.float32, device=device)
        # Overlaid φ-slices for the reader.
        phi_slice_vals = torch.tensor(
            [0.0, np.pi / 2, np.pi, -np.pi / 2],
            dtype=torch.float32, device=device)
        mlp_slice_labels = ["φ=0", "φ=π/2", "φ=π", "φ=−π/2"]
        all_phi = torch.cat([phi_avg, phi_slice_vals])
        n_phi_all = int(all_phi.shape[0])
        # Build the (n_eta * n_phi, 2) per-muon η/φ grids (both muons share
        # (η, φ) — the ThetaNet's output for muon 0 is what gets plotted).
        eta_grid = centers_t[:, None, None].expand(
            n_eta_c, n_phi_all, 2).reshape(-1, 2)
        phi_grid = all_phi[None, :, None].expand(
            n_eta_c, n_phi_all, 2).reshape(-1, 2)
        with torch.no_grad():
            AeM_g, ac_g = model.theta_net(eta_grid, phi_grid)
        # Keep muon-0 output and reshape back to (n_eta, n_phi_all, n_comp).
        AeM_g = AeM_g[:, 0, :].view(n_eta_c, n_phi_all, 3)              # physical scale (×scale_ref inside the net)
        # Effective smear: route the raw MLP (a, c) through
        # _smear_raw_to_effective so the positivity reparam (softplus / square,
        # per --smear-param-form) is applied here too — same operator the model
        # uses internally for the MLP smear branch in _smear_ac_pm. Then
        # multiply by SMEAR_VAR_SCALE for the physical units.
        smear_scale = ac_g.new_tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C])
        ac_eff = model._smear_raw_to_effective(ac_g[:, 0, :])           # reparam(ac)·mask
        ac_phys = (ac_eff * smear_scale).view(n_eta_c, n_phi_all, 2)
        # Split: first n_phi_avg are the φ-average grid; rest are the slices.
        AeM_avg, AeM_slc = AeM_g[:, :n_phi_avg, :], AeM_g[:, n_phi_avg:, :]
        ac_avg,  ac_slc  = ac_phys[:, :n_phi_avg, :], ac_phys[:, n_phi_avg:, :]
        mlp_scale_grid = AeM_avg.mean(dim=1).cpu().numpy()              # [n_eta, 3]  φ-mean
        mlp_scale_band = AeM_avg.std(dim=1).cpu().numpy()               # [n_eta, 3]  φ-std
        mlp_smear_grid = ac_avg.mean(dim=1).cpu().numpy()               # [n_eta, 2]
        mlp_smear_band = ac_avg.std(dim=1).cpu().numpy()                # [n_eta, 2]
        mlp_scale_slices = AeM_slc.cpu().numpy()                        # [n_eta, 4, 3]
        mlp_smear_slices = ac_slc.cpu().numpy()                         # [n_eta, 4, 2]
        # Raw φ-samples (physical) kept for the whitened-basis band/slices —
        # the φ-std must be computed AFTER projecting onto the eigenbasis (a
        # linear combination's spread ≠ the combination of per-component stds).
        mlp_scale_avg_samples = AeM_avg.cpu().numpy()                   # [n_eta, n_phi, 3]
        mlp_smear_avg_samples = ac_avg.cpu().numpy()                    # [n_eta, n_phi, 2]
        # Keep the Fisher-propagated ±1σ bands ONLY when the file was produced
        # for the MLP (θ-net-weight covariance propagated to the per-η outputs);
        # a binned fisher file's per-bin σ does not apply to the net outputs.
        if not mlp_fisher:
            sigma_scale = sigma_smear = None
        # θ-vs-φ plots: sample the net on a FINE φ grid over ALL η-bin centres
        # (the φ-direction companion to the θ-vs-η plots; meaningful only for the
        # continuous MLP θ). A few |η| slices are shown faint; the η-AVERAGE over
        # all bins (+ ±1σ-over-η band) is the headline curve.
        n_phi_fine = 49
        phi_fine = torch.linspace(-np.pi, np.pi, n_phi_fine,
                                  dtype=torch.float32, device=device)
        eg = centers_t[:, None, None].expand(n_eta_c, n_phi_fine, 2).reshape(-1, 2)
        pg = phi_fine[None, :, None].expand(n_eta_c, n_phi_fine, 2).reshape(-1, 2)
        with torch.no_grad():
            AeM_p, ac_p = model.theta_net(eg, pg)
        # [n_eta, n_phi, n_comp] physical, over all η bins
        AeM_all = AeM_p[:, 0, :].view(n_eta_c, n_phi_fine, 3).cpu().numpy()
        ac_all = (model._smear_raw_to_effective(ac_p[:, 0, :]) * smear_scale
                  ).view(n_eta_c, n_phi_fine, 2).cpu().numpy()
        # representative slices (faint curves)
        sl_idx = np.unique(np.linspace(0, n_eta_c - 1, 5).round().astype(int))
        # SAMPLE-weighted η-average: weight each η bin by its Σw occupancy (from
        # the Fisher file's cell_w, summed over φ) so the mean/band reflect where
        # the data actually is; uniform if no occupancy info.
        if fisher_cell_w is not None and fisher_cell_w.shape[0] == n_eta_c:
            wgt = fisher_cell_w.sum(axis=1).astype(np.float64)     # [n_eta]
        else:
            wgt = np.ones(n_eta_c, dtype=np.float64)
        if wgt.sum() <= 0:
            wgt = np.ones(n_eta_c, dtype=np.float64)
        wn = wgt / wgt.sum()

        def _wavg(a):   # a: [n_eta, n_phi, n_comp] → weighted mean/std over η
            mean = np.einsum('e,epc->pc', wn, a)
            var = np.einsum('e,epc->pc', wn, (a - mean[None]) ** 2)
            return mean, np.sqrt(np.clip(var, 0.0, None))
        AeM_m, AeM_b = _wavg(AeM_all)
        ac_m, ac_b = _wavg(ac_all)
        mlp_phi = dict(
            phi=phi_fine.cpu().numpy(),
            AeM=AeM_all[sl_idx], ac=ac_all[sl_idx],     # [n_slc, n_phi, n_comp]
            eta_sl=centers[sl_idx],
            # η-averaged (sample-weighted) curve + spread-over-η band [n_phi,n_comp]
            AeM_mean=AeM_m, AeM_band=AeM_b, ac_mean=ac_m, ac_band=ac_b,
            eta_w=wn)
    if model.scale_enabled:
        # binned θ_scale is the O(1) fit param → ×THETA_SCALE_REF for physical
        # (A,e,M); the MLP grid is already physical (scale_ref inside the net).
        theta_scale = (mlp_scale_grid if mlp_scale_grid is not None
                       else ckpt.get("theta_scale", model.theta_scale.detach()).cpu().numpy()
                       * np.asarray(THETA_SCALE_REF))
        # χ² for compatibility of all θ_scale (A,e,M over the η bins) with the
        # reference (the injected values if a shift was injected, else 0), using
        # the full θ_scale covariance block (correlations included).
        ref = inject_ref_np if inject_ref_np is not None else None
        # 2-D binned: expand the [n_eta,3] φ-averaged reference to per-cell
        # [n_cells,3] (× the injected f_φ sinusoid when --inject-nonuniform) so
        # it aligns with the per-cell θ table / covariance.
        if ref is not None and model.theta_grid:
            ref = _ref_cells(ref, model.n_phi_bins).reshape(-1, 3)
        chi2_info = None
        # Use the TOTAL (data⊕flow) covariance when --flow-uncertainty is set,
        # so the compatibility test accounts for the flow uncertainty too.
        if cov_scale_total is not None:
            resid = theta_scale.reshape(-1)
            if ref is not None:
                resid = resid - ref.reshape(-1)
            chi2, dof, pval = _chi2_compat_zero(resid, cov_scale_total)
            chi2_info = (chi2, dof, pval)
            tgt = "injected" if ref is not None else "0"
            ctag = ("total" if sigma_scale_total is not None else "data-stat")
            print(f"  θ_scale compatibility with {tgt} ({ctag} cov): "
                  f"χ²/dof = {chi2:.1f}/{dof} = "
                  f"{chi2 / max(dof, 1):.2f}, p = {pval:.3g}")
        if model.theta_grid:
            # 2-D binned: per-cell table → η×φ heatmaps (value + pull-vs-ref).
            n_phi_g = model.n_phi_bins
            n_eta_g = theta_scale.shape[0] // n_phi_g
            g = theta_scale.reshape(n_eta_g, n_phi_g, 3)
            # Per-cell σ: the TOTAL (data⊕flow) when --flow-uncertainty is set
            # (consistent with the χ², which uses the total covariance), else
            # the data-stat --fisher σ. sigma_scale alone misses the
            # flow-uncertainty-only case (no --fisher): the per-cell σ then
            # lives in sigma_scale_total/_data, not sigma_scale.
            sig_cell = (sigma_scale_total if sigma_scale_total is not None
                        else sigma_scale)
            sg = (sig_cell.reshape(n_eta_g, n_phi_g, 3)
                  if sig_cell is not None else None)
            plot_theta_grid_etaphi(
                g, sg, ["A", "e [GeV]", "M"], "theta_scale_etaphi",
                stats.eta_edges, n_phi_g, out_dir, chi2_info=chi2_info,
                ref=ref, edm=edm)
            # 1-D projections (φ-mean vs η, η-mean vs φ, + slices) — easier to
            # read/compare than the heatmaps. Pass the full per-cell ref so the
            # vs-φ curve shows the injected f_φ sinusoid under
            # --inject-nonuniform.
            plot_theta_grid_projections(
                g, sg, ["A", "e [GeV]", "M"], "theta_scale",
                stats.eta_edges, out_dir,
                ref=(None if ref is None else ref.reshape(n_eta_g, n_phi_g, 3)),
                cov=cov_scale_total, edm=edm)
        else:
            plot_theta_vs_eta(
                theta_scale, sigma_scale, ["A", "e [GeV]", "M"],
                "theta_scale_vs_eta", stats.eta_edges, out_dir, edm=edm,
                chi2_info=chi2_info, ref=ref,
                band=mlp_scale_band, slices=mlp_scale_slices,
                slice_labels=mlp_slice_labels,
                sigma_band=(model.theta_mode == "mlp"),
                sigma_total=sigma_scale_total,
            )
        # Single-η-bin binned θ: the global (A, e, M) scale params are cheap to
        # scan directly, so add a 1-D NLL likelihood scan (the most transparent
        # uncertainty diagnostic — exposes curvature + any non-parabolicity).
        if (not args.no_theta_scan and model.theta_mode == "binned"
                and model.theta_scale.shape[0] == 1):
            print("plotting θ_scale likelihood scan (single η-bin)...")
            # The scan/gradient must diagnose the fit MINIMUM, so it has to run on
            # the SAME events the fit used — the FULL fitted sample. Build a
            # dedicated loader matching the stage-2 fit (split='train', NO
            # val/holdout carve-out, the fit's validation half) and run on ALL of
            # it (max_events=0), NOT the diagnostics --split (default 'holdout', a
            # non-representative ~5%) and NOT a --max-events subsample. The shards
            # are inhomogeneous, so a partial/wrong sample mislocates the minimum
            # (and the Fisher parabola is also full-sample, for comparability).
            fit_half = (None if (not train_args.get("validation", False)
                                 or train_args.get("no_validation_split", False))
                        else 1)
            scan_loader = JpsiMassArrowLoader(
                shard_files, stats, batch_size=args.batch_size, split="train",
                val_fraction=0.0, holdout_fraction=0.0, drop_last=False,
                half=fit_half,
                inject_theta_scale=inject_np, inject_theta_smear=inject_smear_np,
                inject_seed=int(train_args.get("inject_smear_seed", 12345)),
                cond_basis=train_args.get("cond_basis", "muon_kin"),
                inject_nonuniform=nonuniform,
                inject_prod=inject_prod_np,
                m_window=(model._m_lo_f, model._m_hi_f),
                reco_ptll_min=train_args.get("reco_ptll_min"),
                reco_ptll_max=train_args.get("reco_ptll_max"),
                fit_select=_fit_select_from(train_args))
            plot_theta_scale_likelihood_scan(
                model, scan_loader, device, out_dir,
                scale_fit_params=model.scale_fit_params,
                mc_as_data=mc_as_data, n_iter=args.continuity_n_iter,
                max_events=0, sigma_scale=sigma_scale,
                inject_ref=inject_ref_np, n_points=args.theta_scan_points,
                n_sigma=args.theta_scan_nsigma, progress=True)
    else:
        print("  --disable-scale: skipping theta_scale_vs_eta")
    if model.smearing_enabled:
        # PHYSICAL qop-variance coefficients (a, c) (σ²_qop = a + c·k²); MLP mode
        # samples the net. effective_theta_smear() already applies SMEAR_VAR_SCALE.
        theta_smear_eff = (mlp_smear_grid if mlp_smear_grid is not None
                           else model.effective_theta_smear().detach().cpu().numpy())
        if model.theta_grid:
            n_phi_g = model.n_phi_bins
            n_eta_g = theta_smear_eff.shape[0] // n_phi_g
            g = theta_smear_eff.reshape(n_eta_g, n_phi_g, 2)
            sig_cell = (sigma_smear_total if sigma_smear_total is not None
                        else sigma_smear)
            sg = (sig_cell.reshape(n_eta_g, n_phi_g, 2)
                  if sig_cell is not None else None)
            sref = (_ref_cells(inject_smear_ref_np, n_phi_g)
                    if inject_smear_ref_np is not None else None)
            plot_theta_grid_etaphi(
                g, sg, ["a [qop²]", "c [qop²·GeV²]"], "theta_smear_etaphi",
                stats.eta_edges, n_phi_g, out_dir, edm=edm, ref=sref)
            plot_theta_grid_projections(
                g, sg, ["a [qop²]", "c [qop²·GeV²]"], "theta_smear",
                stats.eta_edges, out_dir, ref=sref, cov=cov_smear_total,
                edm=edm)
        else:
            plot_theta_vs_eta(
                theta_smear_eff, sigma_smear, ["a [qop²]", "c [qop²·GeV²]"],
                "theta_smear_vs_eta", stats.eta_edges, out_dir, edm=edm,
                ref=(inject_smear_ref_np if inject_smear_ref_np is not None else None),
                band=mlp_smear_band, slices=mlp_smear_slices,
                slice_labels=mlp_slice_labels,
                sigma_band=(model.theta_mode == "mlp"),
                sigma_total=sigma_smear_total,
            )
    else:
        print("  --disable-smearing: skipping theta_smear_vs_eta")

    # Plots 2b/3b: θ vs φ at representative |η| slices (MLP only — binned θ has
    # no φ dependence). Companion to the θ-vs-η plots; shows the net's learned
    # φ structure and, under --inject-nonuniform, whether it tracks the injected
    # f_phi = 1 + amp·sin(2φ) sinusoid (the injected ref scales the per-η base by
    # f_η(η_slice)·f_φ(φ)).
    if mlp_phi is not None:
        phi = mlp_phi["phi"]
        # η-averaged f_η so the injected ref matches the η-MEAN closure curve
        # (the η-mean of base·f_η(η)·f_φ(φ) = base·⟨f_η⟩·f_φ(φ)).
        _eta_c2 = 0.5 * (np.asarray(stats.eta_edges[:-1])
                         + np.asarray(stats.eta_edges[1:]))
        feta_mean = (float(_inject_modulation_eta_np(_eta_c2).mean())
                     if nonuniform else 1.0)

        def _phi_ref(base_row):
            # [n_phi, n_comp] = base · ⟨f_η⟩ · f_φ(φ), the η-averaged injected ref.
            if not nonuniform:
                return np.broadcast_to(base_row, (len(phi), len(base_row)))
            fphi = 1.0 + _INJECT_AMP * np.sin(_INJECT_PHI_NOSC * phi)   # [n_phi]
            return base_row[None, :] * (feta_mean * fphi)[:, None]

        eta_w = mlp_phi["eta_w"]                            # [n_eta] sample weights

        def _fisher_sigma_phi(cov2d):
            # φ-RESOLVED statistical ±1σ of the SAMPLE-weighted η-averaged output,
            # from the full 2-D (η,φ_save) Fisher covariance. At each saved φ-bin:
            #   var_c(φ) = Σ_{η,η'} wn_η wn_η' Cov[(η,φ,c),(η',φ,c)],
            # then broadcast/interpolate from the saved φ-bin centres onto the
            # fine plotted φ grid. cov2d: [n_eta,n_phi_s,ncol,n_eta,n_phi_s,ncol].
            if cov2d is None:
                return None
            n_eta_s, n_phi_s, ncol = cov2d.shape[:3]
            if n_eta_s != len(eta_w):
                return None
            sig_s = np.zeros((n_phi_s, ncol))
            for pi in range(n_phi_s):
                for c in range(ncol):
                    blk = cov2d[:, pi, c, :, pi, c]          # [n_eta, n_eta]
                    var = float(np.einsum('e,ef,f->', eta_w, blk, eta_w))
                    sig_s[pi, c] = np.sqrt(max(var, 0.0))
            # map saved φ-bin centres → fine grid by nearest bin (σ is smooth in φ)
            phi_edges_s = np.linspace(-np.pi, np.pi, n_phi_s + 1)
            phi_ctr_s = 0.5 * (phi_edges_s[:-1] + phi_edges_s[1:])
            idx = np.clip(np.searchsorted(phi_edges_s, phi) - 1, 0, n_phi_s - 1)
            return sig_s[idx]                                # [n_phi_fine, ncol]
        if model.scale_enabled:
            sref = (None if inject_np is None
                    else _phi_ref(np.asarray(inject_np)[0]))
            plot_theta_vs_phi(
                phi, mlp_phi["AeM"], ["A", "e [GeV]", "M"],
                "theta_scale_vs_phi", mlp_phi["eta_sl"], out_dir, ref=sref,
                eta_mean=mlp_phi["AeM_mean"], eta_band=mlp_phi["AeM_band"],
                fisher_sigma=_fisher_sigma_phi(cov_scale_2d))
        if model.smearing_enabled:
            cref = (None if inject_smear_np is None
                    else _phi_ref(np.asarray(inject_smear_np)[0]))
            plot_theta_vs_phi(
                phi, mlp_phi["ac"], ["a [qop²]", "c [qop²·GeV²]"],
                "theta_smear_vs_phi", mlp_phi["eta_sl"], out_dir, ref=cref,
                eta_mean=mlp_phi["ac_mean"], eta_band=mlp_phi["ac_band"],
                fisher_sigma=_fisher_sigma_phi(cov_smear_2d))
    elif model.theta_grid:
        print("  --theta-2d-binned: theta_*_vs_eta/_vs_phi show the grid "
              "projections (weighted φ-/η-means + slices); η×φ detail in the "
              "heatmaps")
    elif model.theta_mode != "mlp":
        print("  binned θ (no φ dependence): skipping theta_*_vs_phi")

    # Plots 3b/3c: closure in the DEGENERACY-WHITENED basis. (A,e) and (a,c)
    # are each near-degenerate over the J/ψ pt range, so the m_ll likelihood
    # constrains only the STIFF combination; the orthogonal SLOPPY combination
    # drifts. Projecting the fitted + injected θ onto the (stiff, sloppy)
    # eigenvectors separates "what J/ψ can measure" (stiff — should close) from
    # "what it cannot" (sloppy — large spread, may not close), which is exactly
    # the right way to read the closure when fitting both members of a pair.
    eigb = _degeneracy_eigbasis(getattr(stats, "k_moments", None))
    if model.theta_grid:
        # per-cell (η×φ) θ: the whitened η-curve plots assume a per-η table;
        # skip them (the η×φ heatmaps above carry the 2-D closure instead).
        eigb = None
        print("  --theta-2d-binned: skipping whitened-basis η-curve plots "
              "(see the η×φ heatmaps)")
    if eigb is None:
        if not model.theta_grid:
            print("  no k_moments in stats → skipping whitened-basis closure plots")
    else:
        E_s, l_s, E_c, l_c = eigb
        print("  degeneracy eigenbasis (global, O(1) θ coords; info ratio stiff/sloppy):")
        print("    scale (θ_A,θ_e): stiff×%.0f vec=[%+.3f,%+.3f]  sloppy vec=[%+.3f,%+.3f]"
              % (l_s[0], E_s[0, 0], E_s[1, 0], E_s[0, 1], E_s[1, 1]))
        print("    smear (θ_a,θ_c): stiff×%.0f vec=[%+.3f,%+.3f]  sloppy vec=[%+.3f,%+.3f]"
              % (l_c[0], E_c[0, 0], E_c[1, 0], E_c[0, 1], E_c[1, 1]))

        def _whitened_plot(phys_grid, phys_samples, phys_slices, ref_phys,
                           inj_phys, refE, lam, name, unit_pair, cov_phys=None):
            """Project a 2-component (A,e)/(a,c) set onto the eigenbasis and plot
            stiff/sloppy closure. ``phys_grid`` [n_eta,2]; ``phys_samples``
            [n_eta,n_phi,2] or None (→ φ-spread band); ``phys_slices``
            [n_eta,n_s,2] or None; ``inj_phys`` [n_eta,2] or None. ``cov_phys``
            [n_eta,2,2] (per-bin PHYSICAL covariance of the pair) → the
            Fisher ±1σ band IN the stiff/sloppy basis: standardise (÷ref⊗ref)
            then rotate, σ_w = √diag(Eᵀ C_o E)."""
            grid_w = _project_whitened(phys_grid, ref_phys, refE)        # [n_eta,2]
            band_w = None
            if phys_samples is not None:
                samp_w = _project_whitened(phys_samples, ref_phys, refE)  # [n_eta,n_phi,2]
                grid_w = samp_w.mean(axis=1)
                band_w = samp_w.std(axis=1)
            slices_w = (None if phys_slices is None
                        else _project_whitened(phys_slices, ref_phys, refE))
            ref_w = None if inj_phys is None else _project_whitened(inj_phys, ref_phys, refE)
            sigma_w = None
            if cov_phys is not None:
                rf = np.asarray(ref_phys, dtype=np.float64)
                co1 = (np.asarray(cov_phys, dtype=np.float64)
                       / (rf[None, :, None] * rf[None, None, :]))      # ÷ref⊗ref
                cw = np.einsum('ij,nik,kl->njl', refE, co1, refE)      # Eᵀ C_o E
                sigma_w = np.sqrt(np.clip(
                    np.diagonal(cw, axis1=1, axis2=2), 0.0, None))     # [n_eta,2]
            names = [f"STIFF (measured, info×{lam[0]:.0f})",
                     f"SLOPPY (degenerate, info×1)"]
            plot_theta_vs_eta(
                grid_w, sigma_w, names, name, stats.eta_edges, out_dir, edm=edm,
                ref=ref_w, band=band_w, slices=slices_w,
                slice_labels=mlp_slice_labels,
                sigma_band=(model.theta_mode == "mlp"))

        scale_pair_fit = ("A" in model.scale_fit_params and "e" in model.scale_fit_params)
        if model.scale_enabled and scale_pair_fit:
            sg = (mlp_scale_grid[:, :2] if mlp_scale_grid is not None
                  else (ckpt.get("theta_scale", model.theta_scale.detach()).cpu().numpy()
                        * np.asarray(THETA_SCALE_REF))[:, :2])
            ss = None if mlp_scale_avg_samples is None else mlp_scale_avg_samples[:, :, :2]
            sl = None if mlp_scale_slices is None else mlp_scale_slices[:, :, :2]
            ij = None if inject_ref_np is None else inject_ref_np[:, :2]
            # Per-bin (A,e) physical covariance for the band (diagonal blocks of
            # the θ_scale covariance — TOTAL data⊕flow when --flow-uncertainty).
            cov_ae = None
            if cov_scale_total is not None:
                ne = cov_scale_total.shape[0] // 3
                c4 = cov_scale_total.reshape(ne, 3, ne, 3)
                cov_ae = np.stack([c4[b, :2, b, :2] for b in range(c4.shape[0])])
            _whitened_plot(sg, ss, sl, THETA_SCALE_REF[:2], ij, E_s, l_s,
                           "theta_scale_whitened_vs_eta", ("A", "e"), cov_phys=cov_ae)
        if model.smearing_enabled and model.smear_fit_params == "both":
            cg = (mlp_smear_grid if mlp_smear_grid is not None
                  else model.effective_theta_smear().detach().cpu().numpy())
            cs = mlp_smear_avg_samples
            csl = mlp_smear_slices
            ij = inject_smear_ref_np
            cov_ac = None
            if cov_smear_total is not None:
                ne = cov_smear_total.shape[0] // 2
                c4 = cov_smear_total.reshape(ne, 2, ne, 2)
                cov_ac = np.stack([c4[b, :, b, :] for b in range(c4.shape[0])])
            _whitened_plot(cg, cs, csl, [SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C],
                           ij, E_c, l_c, "theta_smear_whitened_vs_eta", ("a", "c"),
                           cov_phys=cov_ac)

    # Plot 5: pulls.
    print("plotting per-bin pulls...")
    plot_pulls(
        evals, m_centers_np, eta_slice_edges,
        model._m_lo_f, model._m_hi_f, out_dir,
    )

    # Plot 6: MC closure (forward-folded MC vs flow density curve, both at
    # the fitted scale + smearing). Re-uses pred_signal_mc — no second pass.
    print("plotting MC closure...")
    plot_mc_closure(evals, m_centers_np, eta_slice_edges, out_dir)
    if getattr(model, "background_enabled", True) and "bkg_frac" in evals:
        plot_bkg_fractions(evals["bkg_frac"], stats.eta_edges, out_dir,
                           inject_bkg=inject_bkg_np,
                           bkg_model=getattr(model, "bkg_model", "bernstein"))

    # Plot 7: parameter-sensitivity slices — model density at ±shifts of each
    # fitted param, in conditional slices chosen to break degeneracies.
    # Re-iterates the loader (own pass; extra grid evals per shift). Uses the
    # same --max-events as the main closure above, counted the same way (all
    # events per batch), so both passes select the same events.
    if not bool(getattr(args, "no_param_sensitivity", False)):
        print("plotting parameter-sensitivity slices...")
        plot_param_sensitivity(
            model, loader, stats, m_centers, out_dir,
            shift_scale=args.param_shift, max_events=args.max_events,
            chunk_events=args.grid_chunk_events, n_iter=args.continuity_n_iter,
            device=device, mc_as_data=mc_as_data,
            progress=getattr(args, "progress", True))

    print(f"\nall diagnostics written under {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
