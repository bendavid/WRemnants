"""Mixture model + flow for the unbinned J/ψ mass calibration.

Two-stage continuity design. The flow models only the NOMINAL (θ=0) mass shape
``p₀(m | muon_kin)`` and never conditions on θ; the θ-dependence is supplied
analytically in stage 2 (the continuity tilt, ``data_nll_continuity``).

* ``theta_scale`` ∈ R^{24×3}  — per-η-bin (A, e, M) muon scale nuisances.
* ``theta_smear`` ∈ R^{24×2}  — per-η-bin (a, c) signed width-smear coefficients.

Stage 1 (``log_p_nominal``): train the conditional flow ``p₀(m | muon_kin)`` on
simulation at θ=0 (leak-free kinematic conditioning only).

Stage 2 (``data_nll_continuity``, frozen flow): the signal density is the
nominal flow forward-folded analytically by a deterministic, invertible map
``x = μ + (1+s)(m' + s_adv(m') − μ)`` — a scale advection ``s_adv`` plus a signed
mass-density stretch ``s`` (μ = mean m_ll). ``s`` is the variance-equivalent of a
per-muon qop smear with VARIANCE ``σ²_qop = a + c·k²`` (signed/two-sided): the
smear is applied as a score-driven probability-flow displacement that reproduces
a Gaussian qop smear of mass-variance ``σ²_qop = a + c·k²`` (``_continuity_logp``),
and the SAME (a, c) drive the per-muon qop fold in the validation plots. Evaluated at
the source pre-image with the change-of-variables Jacobian (``_continuity_logp``)
— no flow derivatives, normalised by construction; the smear is a pure mass-space
stretch so it leaves the ρ conditioning untouched. Mixed with a degree-1
Bernstein background via the MLP ``f(c)``:

  data event:  p(m | c, θ) = f_0(c) p_0 + f_1(c) p_1 + (1 − f_0 − f_1)·p_s(m | c, θ)
  MC event:    p(m | c, θ) = p_s(m | c, θ)
"""

from __future__ import annotations

import contextlib
import math
import os
import sys
import warnings

import torch
import torch.nn as nn
import torch.nn.functional as F

# Re-use ``build_flow`` from the existing trainer.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from train_muon_response_flow import FlowWithLogProb, build_flow  # noqa: E402
from compact_flow import CompactMatchedFlow  # noqa: E402
from mixture_flow import MixtureFlow  # noqa: E402
from dcb_density import DCBDensity, EGEDensity  # noqa: E402
from nce_density import NCEDensity  # noqa: E402


# ---------------------------------------------------------------------------
# zuko monkey-patch: floor the jacobian before .log() in MonotonicTransform.
# ---------------------------------------------------------------------------
#
# zuko's ``MonotonicTransform.call_and_ladj`` returns ``y, jacobian.log()``.
# For the GF flow's GMM-CDF base, the mathematical jacobian is always
# strictly positive, but it underflows to *exactly* 0 in float32 when the
# input is far from every mixture component — typical at random init.
#  → ``log(0) = -inf`` forward → ``log_p = -inf`` → ``-log p = +inf``
#    poisons the per-batch NLL on the very first batch.
# The backward of ``log(0)`` is also degenerate (``grad_out / 0`` → ±inf
# or NaN), which is the path the anomaly tracer initially reported.
#
# Floor the jacobian at ``1e-30`` before the log: forward becomes finite
# (worst-case ``log(1e-30) ≈ -69``) and the backward is bounded
# (``clamp_min`` gradient is 0 in the clamped region, so the second-order
# autograd through the GF's f-derivative stays clean). Events at the
# clamp boundary contribute zero gradient through this layer — they're
# already in a "no signal" regime where the flow can't tell which way
# to move its mixture, so dropping their gradient is the right behaviour
# anyway; the conditioner network still receives gradient through the
# other (un-clamped) layers.
#
# We additionally clamp the transform INPUT ``x`` before ``self.f(x)``.
# The GF base maps via ``erf(x / √2)``, which in float32 saturates to
# exactly ±1 for ``|x| ≳ 5.7`` — sending the subsequent inverse-CDF / GMM
# map to ±∞ and NaN-ing the backward (the ErfBackward0 → Mul/Exp NaN the
# anomaly tracer reports). ``x`` here is *post-conditioner-affine*, so an
# external mass clamp is not enough: a sharp conditioner (or a pathological
# forward-fold mass) can push ``x`` past the saturation point even for an
# in-window event. Clamping at ±5 (erf(5/√2)=erf(3.54)=0.9999994, finite)
# is well beyond the ~N(0,1) regime the flow maps real events to, so it
# only ever bounds the pathological tail; the jacobian is taken w.r.t. the
# clamped value so the density stays self-consistent there.
def _install_zuko_jacobian_floor(
    min_jac: float = 1e-30, x_clamp: float = 5.0
) -> None:
    import zuko.transforms as _zt

    if getattr(_zt.MonotonicTransform, "_jpsi_jacobian_floored", False):
        return

    _orig_call_and_ladj = _zt.MonotonicTransform.call_and_ladj

    def _safe_call_and_ladj(self, x):
        # Reproduce the zuko logic but (a) clamp the erf input and (b) clamp
        # the jacobian before log.
        create_graph = torch.is_grad_enabled() and (
            x.requires_grad or bool(self.phi)
        )
        with torch.enable_grad():
            x = x.clone().requires_grad_()
            x_safe = x.clamp(-x_clamp, x_clamp)
            y = self.f(x_safe)
        jacobian = torch.autograd.grad(
            y, x_safe, torch.ones_like(y), create_graph=create_graph,
        )[0]
        return y, jacobian.clamp_min(min_jac).log()

    _zt.MonotonicTransform.call_and_ladj = _safe_call_and_ladj
    _zt.MonotonicTransform._jpsi_jacobian_floored = True


_install_zuko_jacobian_floor()


# ---------------------------------------------------------------------------
# zuko monkey-patch: guard the Gaussianization transform's internal overflow.
# ---------------------------------------------------------------------------
#
# ``GaussianizationTransform.f`` computes
#     erf((x · exp(scale_i) + shift_i) / √2)   →  mean over i  →  erfinv
# where (shift_i, scale_i) are the per-component parameters the *conditioner*
# predicts. Two float32 failure modes, both reached when the conditioner
# outputs large values (e.g. when the standardised θ-conditioning grows):
#   1. ``self.scale = exp(scale_i)`` overflows to +inf for scale_i ≳ 88 →
#      ``x · inf`` is inf (or 0·inf = NaN when x = 0) → poisons forward AND
#      backward (0·inf in the Mul/Exp backward — the reported MulBackward0).
#   2. the erf argument saturates erf→±1 for |arg| ≳ 5.7 → its backward is a
#      0·inf NaN even though the forward (capped by the ·(1−1e-6) term) is OK.
# Clamp scale_i before the exp (keeps self.scale finite, huge margin to inf)
# and clamp the erf argument (keeps erf away from saturation). Both clamps
# have zero gradient in the clamped region, so pathological events stop
# contributing gradient cleanly instead of NaN-ing the batch. Healthy flows
# map real events to |arg| ~ O(1) with scale_i ~ O(1), so neither clamp ever
# fires there.
def _install_gaussianization_guards(
    scale_param_clamp: float = 30.0, erf_arg_clamp: float = 5.0
) -> None:
    import zuko.transforms as _zt

    G = _zt.GaussianizationTransform
    if getattr(G, "_jpsi_gf_guarded", False):
        return

    _orig_init = G.__init__

    def _safe_init(self, shift, scale, **kwargs):
        scale = scale.clamp(-scale_param_clamp, scale_param_clamp)
        _orig_init(self, shift, scale, **kwargs)

    def _safe_f(self, x):
        arg = x[..., None] * self.scale + self.shift
        arg = arg.clamp(-erf_arg_clamp, erf_arg_clamp)
        y = torch.erf(arg / math.sqrt(2))
        y = torch.mean(y, dim=-1) * (1 - 1e-6)
        y = torch.erfinv(y) * math.sqrt(2)
        return y

    G.__init__ = _safe_init
    G.f = _safe_f
    G._jpsi_gf_guarded = True


_install_gaussianization_guards()


def _noise_active(sigma) -> bool:
    """True if the given σ should trigger noise sampling.

    Scalar 0 / None → no noise (use parameter as-is). Any Tensor → noise
    (caller is responsible for non-negative values). Positive scalar →
    noise. The trainer's adaptive-σ helper returns Tensors when adaptive
    sampling is active and the fixed scalar during warmup.
    """
    if sigma is None:
        return False
    if isinstance(sigma, torch.Tensor):
        return True
    return float(sigma) > 0.0


# Hard clamp on the standardised mass fed to the GF flow. In float32,
# erf(x/√2) saturates to exactly ±1 for |x| ≳ 5.7, which sends the flow's
# inverse-CDF map to ±∞ and NaNs the exp(−z²/2) in the backward pass. The
# J/ψ window cut (2.92–3.28 GeV) keeps real events at |mll_std| ≲ 3.9, so
# this only ever bounds pathological MC forward-fold tails (a near-zero qop
# from a large sampled smear/scale → huge mass); those tail events get the
# clamp-edge density and contribute no gradient, instead of crashing.
MLL_STD_FLOW_CLAMP = 5.0

# Match the convention of make_jpsi_crctn_helper: 24 η-bins, (A, e, M).
N_ETA_BINS = 24
N_THETA_SCALE = 3      # (A, e, M)
N_THETA_SMEAR = 2      # (a, c)
N_THETA_SCALE_PM = 2 * N_THETA_SCALE   # 6 — flat per-event scale vector
N_THETA_SMEAR_PM = 2 * N_THETA_SMEAR   # 4 — flat per-event smear vector

# Conditioning sizes.
#
#   muon_kin_std   7  = (η_+, η_-, cos φ_+, sin φ_+, cos φ_-, sin φ_-, ρ)
#                       with ρ = (pt_+ − pt_-)/(pt_+ + pt_-)
#   theta_scale_pm 6  = (A_+, e_+, M_+, A_-, e_-, M_-)         ← scale conditioning
#   theta_smear_pm 4  = (a_+, c_+, a_-, c_-)                   ← smear conditioning
#
# ``muon_kin`` is the *leak-free* kinematic conditioning: η_±, φ_± (as cos/sin
# to be wrap-free and keep φ_± recoverable for the φ-dependent detector
# response) and the pt asymmetry ρ. These span 5 of the 6 dimuon DOF, leaving
# the pt *scale* ↔ m_ll free — so the flow's target is not determined by its
# conditioning. The SIGNAL FLOW conditions on (muon_kin, θ_scale_pm,
# θ_smear_pm); the background-fraction MLP conditions on muon_kin alone (same
# kinematics, no nuisances). ``y_event`` (dilepton-level vars) is the OPTIONAL
# alternative basis selected by ``cond_basis='event_level'`` (--cond-basis):
# (yll, ln(ptll/mll), cosPhill, sinPhill, cosθ*, sinφ*, cosφ*). Every component
# is DIMENSIONLESS → invariant under the common-pt dilation that is the mass
# direction of the 1-D conditional (plain ln ptll moves 1:1 with the mass there
# and is forbidden — it makes the fixed-conditioning density model inconsistent
# at first order in θ). Residual pt-dependence (e/M/smear, ρ-class) propagates
# through the per-node recompute (see _cond_from_muons). Default stays the
# leak-free muon_kin.
N_Y_EVENT = 7    # event-level basis dim (== N_MUON_KIN, so no flow-shape change)
N_MUON_KIN = 7
N_FLOW_COND = N_MUON_KIN + N_THETA_SCALE_PM + N_THETA_SMEAR_PM  # 17

# Muon rest mass in GeV (J/ψ analyses use this everywhere).
MUON_MASS_GEV = 0.1056583755

# Fixed reference scales for the per-muon qop-resolution VARIANCE parameters
# (σ²_qop = a·SCALE_A + c·SCALE_C·k², k=1/pt). The physical qop variance is
# O(1e-7) (σ_qop ~ 3e-4), a terrible optimizer scale; these put the fitted (and
# injected) (a, c) at O(1) so the standard smear LR works and a runaway is
# bounded to a physical σ_qop. Calibrated so θ≈1 ≈ a +20% m_ll-variance smear on
# the J/ψ sample: a≈8e-8 (c=0) and c≈2.3e-5 (a=0) each broaden Var(m_ll) by ~20%,
# i.e. θ_a≈0.8 at SCALE_A and θ_c≈1.1 at SCALE_C. (SCALE_C was 1e-6, which left
# θ_c≈20 — far from O(1) — so θ_c chronically under-converged at the smear LR.)
SMEAR_VAR_SCALE_A = 1e-7
SMEAR_VAR_SCALE_C = 2e-5

# Default RAW θ_smear init for the 'square' reparam (physical = raw²·SCALE).
# At raw=0 the square map is degenerate: effective=0 AND ∂effective/∂raw=2·raw=0,
# so the smear branch gets ZERO gradient and is frozen at the zero init (observed:
# θ_smear‖∞ stuck at 0 across a real fit). Initialising raw at this nonzero value
# starts the fit off that saddle; 0.5 is chosen so the square Jacobian 2·raw = 1
# at init — the same gradient scale the 'linear'/'softplus' forms have at their
# init — so the smear lr is calibrated identically across forms. (effective at
# init = 0.25, i.e. physical a,c = 0.25·SCALE, a small nonzero start.)
SMEAR_SQUARE_INIT_RAW = 0.5

# Singularity guard for the qop→pt inversion pt = |sinθ / qop|. A scale/smear
# shift adds to qop; if it drives qop EXACTLY through zero, pt → ∞. We take the
# magnitude (pt is positive; a sign flip of qop is a PHYSICAL charge mis-
# reconstruction, kept) and only clamp |qop| at this floor to keep pt finite —
# NOT a resolution-suppressing floor on the shift itself. Chosen tiny so it
# bites only on the measure-zero qop=0 case (|qop| is ~1e-3..1e-1 physically).
QOP_EPS = 1e-9

# Invertibility floor on the probability-flow smear Jacobian G' = dx/dm'. A valid
# forward (broadening) transport has G' > 0; the floor catches the fold/over-
# sharpen region (V·∂²_m log p₀ large) gracefully, mirroring the old (1+s) ≥ 0.05.
# Used by jacobian_form="softlog" as the seam where the analytic log is replaced
# by a C¹ linear tangent extension (see _softlog_below_floor). Unused by "exp".
SMEAR_GP_FLOOR = 0.05


# Gauss-Hermite nodes for the analytic Gaussian-convolution smear operator
# (smear_operator="gh_convolution"). Cached per (n, device, dtype).
_GH_CACHE: dict = {}

# Max rows per flow evaluation. The qop operator flattens B·n_gh² rows into a
# single flow call (8× the mass-space operator's B·n_gh at n_gh=8), big enough
# that the CUDA allocator can fail on one matmul (NVML_SUCCESS assert) at large
# batch. Flow evals on more than this many rows are split into chunks (the flow
# is per-row independent, so this is exact). ~5·10⁵ matches the mass-space
# operator's working size, which runs fine.
_FLOW_EVAL_CHUNK = 1 << 19  # 524288


def _gh_nodes(n: int, device, dtype):
    """Gauss–Hermite nodes/log-weights for ``E_{ε~N(0,1)}[f] ≈ Σ_i W_i f(ξ_i)``.

    ``∫ f(t) e^{-t²} dt ≈ Σ w_i^H f(t_i^H)`` ⇒ with ``ξ = √2 t^H`` and
    ``W = w^H/√π`` (so Σ W = 1) we get the standard-normal expectation.
    Cached. Returns ``(ξ [n], logW [n])``."""
    key = (int(n), str(device), str(dtype))
    if key not in _GH_CACHE:
        import numpy as _np
        t, w = _np.polynomial.hermite.hermgauss(int(n))
        xi = _np.sqrt(2.0) * t
        logW = _np.log(w) - 0.5 * _np.log(_np.pi)
        _GH_CACHE[key] = (
            torch.as_tensor(xi, device=device, dtype=dtype),
            torch.as_tensor(logW, device=device, dtype=dtype),
        )
    return _GH_CACHE[key]


def _softlog_below_floor(g: torch.Tensor, floor: float = SMEAR_GP_FLOOR) -> torch.Tensor:
    """``log(g.clamp_min(floor))`` PLUS a quadratic barrier ``(floor − g)⁺² /
    (2·floor²)`` for ``g < floor``. Drop-in replacement for the hard
    ``torch.log(g.clamp_min(floor))`` that:

    - matches it EXACTLY in the physical region (``g ≥ floor`` → barrier = 0),
    - is BOUNDED from below by ``log(floor)`` for any g (so the density
      cannot diverge to ±∞ from a degenerate or extreme-negative Jacobian
      — the linear-tangent variant ran away to ``−∞`` for very negative g,
      driving NLL ``→ −∞`` — a *reward* for unphysical G' < 0, the opposite
      of what we want),
    - is C¹ at the seam (barrier value and slope both zero at ``g = floor``)
      and grows quadratically past it, so the gradient w.r.t. g is
      ``−(floor − g)/floor²`` for ``g < floor`` — pulling g BACK toward the
      physical region, replacing the hard clamp's flat-NLL basin with a
      gentle restoring force whose strength scales with how far past the
      floor we've drifted.

    Returns ``log_Gp_effective`` to plug into ``−log p_θ = log_Gp_effective
    − log p₀(m')``."""
    log_safe = torch.log(g.clamp_min(floor))
    deficit = (floor - g).clamp_min(0.0)
    barrier = deficit * deficit / (2.0 * floor * floor)
    return log_safe + barrier

# Fixed reference scales for the per-muon SCALE parameters (A, e, M). The
# physical values are O(1e-4 / 1e-3 / 1e-5); these put the fitted (and injected)
# θ_scale at O(1) so the standard scale LR works and all three components share a
# well-conditioned step (physical A,e,M = θ_scale · THETA_SCALE_REF).
THETA_SCALE_REF = (1e-4, 1e-3, 1e-5)


# ---------------------------------------------------------------------------
# Degeneracy-whitening preconditioner
# ---------------------------------------------------------------------------
# (A, e) and (a, c) are each strongly degenerate over the narrow J/ψ pt range:
# the scale response basis (1, −k) for (A, e) and the smear-variance basis
# (1, k²) for (a, c) are nearly collinear, so the loss has a long tilted valley
# that Adam — a purely DIAGONAL preconditioner — cannot navigate. We rotate the
# gradient of each pair into a decorrelated basis with a fixed matrix L built
# from the per-η-bin curvature moments (⟨k⟩, ⟨k²⟩, ⟨k⁴⟩).
#
# The rotation is applied as a BACKWARD-ONLY transform (``_WhitenGradFn``):
# the forward pass is the identity, so the forward map, the softplus positivity
# reparam, the likelihood value, and the observed Fisher are ALL unchanged — only
# the optimizer's step direction is preconditioned. It is gated on ``self.training``
# so it is fully inert (identity backward too) during Fisher / diagnostics / eval,
# leaving the covariance in the physical θ basis exact. Adam's subsequent diagonal
# rescaling acts on the already-decorrelated axes and cannot reintroduce the
# cross-correlation, so Adam is kept for all parameters.


class _WhitenGradFn(torch.autograd.Function):
    """Identity forward; in backward left-multiplies the gradient w.r.t. the
    last (size-2) axis by a fixed lower-triangular whitening matrix ``L`` (the
    same gradient a forward reparam ``p = L·θ'`` would produce, but without
    touching the forward value). ``L`` is either ``[2, 2]`` (global, MLP mode)
    or batched ``[..., 2, 2]`` (per-η-bin, binned mode)."""

    @staticmethod
    def forward(ctx, x, L):
        ctx.save_for_backward(L)
        return x

    @staticmethod
    def backward(ctx, g):
        (L,) = ctx.saved_tensors
        if L.dim() == 2:
            gx = torch.einsum("...i,ij->...j", g, L)
        else:
            gx = torch.einsum("...i,...ij->...j", g, L)
        return gx, None


def _whitening_L(rho: float, max_rho: float = 0.99) -> torch.Tensor:
    """Lower-triangular ``L`` with ``L Lᵀ = C⁻¹`` for the 2×2 correlation
    ``C = [[1, ρ],[ρ, 1]]`` (ρ clamped to ±``max_rho`` to bound the whitening
    of the near-degenerate direction). Satisfies ``Lᵀ C L = I``, so a parameter
    pair whose loss curvature has correlation ρ becomes unit-conditioned in the
    rotated coordinates. Returns ``[2, 2]`` float32."""
    r = float(max(-max_rho, min(max_rho, rho)))
    C = torch.tensor([[1.0, r], [r, 1.0]], dtype=torch.float64)
    Cinv = torch.linalg.inv(C)
    L = torch.linalg.cholesky(Cinv)
    return L.to(torch.float32)


# ---------------------------------------------------------------------------
# Bernstein degree-1 background basis (unchanged)
# ---------------------------------------------------------------------------


def bernstein_d1(
    mll: torch.Tensor, m_lo: float, m_hi: float
) -> tuple[torch.Tensor, torch.Tensor]:
    """Degree-1 Bernstein basis on ``[m_lo, m_hi]``, normalised to ∫=1.

    Densities in raw m_ll units (1/GeV) so the mixture sum is dimensionally
    consistent with ``exp(log p_signal)``.
    """
    width = m_hi - m_lo
    u = (mll - m_lo) / width
    p0 = (1.0 - u) * 2.0 / width
    p1 = u * 2.0 / width
    return p0, p1


_BKG_BINOM_CACHE: dict = {}


def bernstein_basis_n(mll: torch.Tensor, m_lo: float, m_hi: float,
                      degree: int) -> torch.Tensor:
    """Degree-n Bernstein basis DENSITIES on ``[m_lo, m_hi]``, each ∫ = 1:

        p_i(m) = (n+1)·C(n,i)·u^i·(1−u)^{n−i} / width,  i = 0..n.

    Returns ``[..., n+1]`` in 1/GeV (mixture-compatible with the signal
    density). degree=1 reproduces ``bernstein_d1``. The binomial row is
    cached per (degree, device, dtype) — pre-warm under torch.compile (the
    miss branch is an untraceable constant build, same pattern as
    ``_gh_nodes``)."""
    key = (int(degree), str(mll.device), str(mll.dtype))
    if key not in _BKG_BINOM_CACHE:
        import math as _math
        _BKG_BINOM_CACHE[key] = (
            torch.tensor([_math.comb(degree, k) for k in range(degree + 1)],
                         device=mll.device, dtype=mll.dtype),
            torch.arange(degree + 1, device=mll.device, dtype=mll.dtype))
    binom, k = _BKG_BINOM_CACHE[key]
    width = m_hi - m_lo
    u = ((mll - m_lo) / width).unsqueeze(-1).clamp(0.0, 1.0)
    # m_lo/m_hi may be per-event [B] tensors (fit-time cut → per-event window);
    # the (degree+1)/width prefactor then needs a trailing axis to broadcast
    # against the [..., n+1] basis. Scalar edges keep the float fast path.
    pref = (degree + 1.0) / width
    if torch.is_tensor(pref):
        pref = pref.unsqueeze(-1)
    return pref * binom * u.pow(k) * (1.0 - u).pow(degree - k)


def exp_bkg_density(mll: torch.Tensor, m_lo: float, m_hi: float,
                    slope: torch.Tensor) -> torch.Tensor:
    """Window-normalised EXPONENTIAL background density, parameterised by the
    DIMENSIONLESS slope ``s`` (= λ·width; scalar or per-event [..] tensor —
    the conditioning-dependent MLP head):

        p(m) = s·e^{−s·u} / (width·(1 − e^{−s})),  u = (m − m_lo)/width,

    signed (s > 0 falls with mass, s < 0 rises), ∫_window = 1 exactly for any
    s. The s → 0 limit is the uniform density 1/width — taken through a
    2nd-order expansion below |s| < 1e-4 (the exact form is 0/0 there; the
    expansion keeps the s-gradient alive through the limit, unlike a
    where-switch to a constant)."""
    width = m_hi - m_lo
    u = ((mll - m_lo) / width).clamp(0.0, 1.0)
    small = slope.abs() < 1e-4
    s_safe = torch.where(small, torch.ones_like(slope), slope)
    p_exact = s_safe * torch.exp(-s_safe * u) / (-torch.expm1(-s_safe) * width)
    # p(u; s) = (1/width)·(1 + s·(½ − u) + s²·(u²/2 − u/2 + 1/12) + O(s³))
    p_small = (1.0 + slope * (0.5 - u)
               + slope * slope * (0.5 * u * u - 0.5 * u + 1.0 / 12.0)) / width
    return torch.where(small, p_small, p_exact)


# ---------------------------------------------------------------------------
# Event-level kinematic helpers (autograd-friendly)
# ---------------------------------------------------------------------------


def _sintheta_from_eta(eta_pm: torch.Tensor) -> torch.Tensor:
    return 1.0 / torch.cosh(eta_pm)


def _event_mll(
    pt_pm: torch.Tensor, eta_pm: torch.Tensor, phi_pm: torch.Tensor
) -> torch.Tensor:
    """Two-body invariant mass for muons of mass ``MUON_MASS_GEV``.

    Inputs ``[B, 2]`` per (+, −). Output ``[B]``. Autograd-friendly.
    """
    px = pt_pm * torch.cos(phi_pm)
    py = pt_pm * torch.sin(phi_pm)
    pz = pt_pm * torch.sinh(eta_pm)
    p2 = px * px + py * py + pz * pz
    E = torch.sqrt(p2 + MUON_MASS_GEV * MUON_MASS_GEV)
    Etot = E.sum(-1)
    Px = px.sum(-1)
    Py = py.sum(-1)
    Pz = pz.sum(-1)
    m2 = Etot * Etot - (Px * Px + Py * Py + Pz * Pz)
    return torch.sqrt(m2.clamp_min(1e-12))


def _dm_dpt_analytic(
    pt_pm: torch.Tensor, eta_pm: torch.Tensor, phi_pm: torch.Tensor
) -> torch.Tensor:
    """Closed-form ``∂m_ll/∂pt_μ`` at fixed (η, φ) — the exact derivative of
    ``_event_mll``:

        ∂m/∂pt_μ = [E_tot·pt_μ·cosh²η_μ/E_μ
                    − (P_x cosφ_μ + P_y sinφ_μ + P_z sinhη_μ)] / m

    (from m² = E_tot² − |ΣP|², with dE/dpt = pt·cosh²η/E and dp/dpt =
    (cosφ, sinφ, sinhη)). Inputs ``[..., 2]`` per muon; output ``[..., 2]``.
    Used by the closed-form operator Jacobians (replacing autograd through
    the kinematic algebra)."""
    cphi, sphi, sheta = torch.cos(phi_pm), torch.sin(phi_pm), torch.sinh(eta_pm)
    px = pt_pm * cphi
    py = pt_pm * sphi
    pz = pt_pm * sheta
    E = torch.sqrt(px * px + py * py + pz * pz
                   + MUON_MASS_GEV * MUON_MASS_GEV)
    Etot = E.sum(-1, keepdim=True)
    Px = px.sum(-1, keepdim=True)
    Py = py.sum(-1, keepdim=True)
    Pz = pz.sum(-1, keepdim=True)
    m = torch.sqrt((Etot * Etot - (Px * Px + Py * Py + Pz * Pz)).clamp_min(1e-12))
    cheta2 = torch.cosh(eta_pm) ** 2
    return (Etot * pt_pm * cheta2 / E - (Px * cphi + Py * sphi + Pz * sheta)) / m


def _event_cond_raw(
    pt_pm: torch.Tensor, eta_pm: torch.Tensor, phi_pm: torch.Tensor,
    eps: float = 1e-6,
) -> torch.Tensor:
    """Event-level conditioning ``(yll, ln(ptll/mll), cosPhill, sinPhill, cosθ*,
    sinφ*, cosφ*)`` from the two muon momenta — a differentiable torch port of
    the ROOT snapshot's dilepton + Collins–Soper computation
    (wremnants/production/include/csVariables.hpp ``csSineCosThetaPhi``).

    The dilepton-pt component is the DIMENSIONLESS ``ln(ptll/mll)``: along the
    mass direction (a common pt dilation, the slice the 1-D conditional + its
    window norm are taken on) ptll scales 1:1 with mll, so the ratio — like
    yll, φll and the CS angles — is dilation-INVARIANT. This is the leak-free
    requirement: a conditioner that shifts coherently with the mass coordinate
    (e.g. plain ln ptll) makes the fixed-conditioning 1-D density model
    inconsistent at first order in θ (missing conditioning-Jacobian/marginal
    terms on the steep ptll spectrum → the A/e sloppy-direction blow-up seen
    in validation). Residual e/M/smear shifts of the ratio are the same class
    as muon_kin's ρ and are carried by the per-node recompute.

    Inputs ``[..., 2]`` per (+, −); ``eta_pm``/``phi_pm`` broadcast against
    ``pt_pm`` (so a per-GH-node ``pt`` [.,G²,2] works with [.,1,2] η/φ). Output
    ``[..., 7]`` in the order of the loader's ``_Y_EVENT_FEATURES``. The μ+ is the
    antilepton (index 0), the μ− the lepton (index 1) — matching the snapshot."""
    cphi = torch.cos(phi_pm)
    sphi = torch.sin(phi_pm)
    px = pt_pm * cphi
    py = pt_pm * sphi
    pz = pt_pm * torch.sinh(eta_pm)
    E = torch.sqrt(px * px + py * py + pz * pz + MUON_MASS_GEV * MUON_MASS_GEV)
    Px, Py, Pz, Etot = px.sum(-1), py.sum(-1), pz.sum(-1), E.sum(-1)
    ptll = torch.sqrt((Px * Px + Py * Py).clamp_min(eps * eps))
    mll = torch.sqrt(
        (Etot * Etot - (Px * Px + Py * Py + Pz * Pz)).clamp_min(eps * eps))
    yll = 0.5 * torch.log(((Etot + Pz) / (Etot - Pz)).clamp_min(eps))
    cosPhill = Px / ptll
    sinPhill = Py / ptll
    # Boost everything into the dilepton rest frame: velocity b = −P/Etot
    # (ROOT BoostToCM), γ = 1/√(1−b²).
    bx, by, bz = -Px / Etot, -Py / Etot, -Pz / Etot
    b2 = (bx * bx + by * by + bz * bz).clamp(max=1.0 - 1e-9)
    gamma = torch.rsqrt(1.0 - b2)
    b2s = b2.clamp_min(1e-30)

    def _boost_unit(qx, qy, qz, qE):
        bdotp = bx * qx + by * qy + bz * qz
        fac = (gamma - 1.0) * bdotp / b2s + gamma * qE
        rx, ry, rz = qx + fac * bx, qy + fac * by, qz + fac * bz
        r = torch.sqrt((rx * rx + ry * ry + rz * rz).clamp_min(eps * eps))
        return rx / r, ry / r, rz / r

    # Proton beams: E=6500, m_p=0.938; beam ±z sign tracks the dilepton z
    # (copysign(1, Pz): +1 at Pz==0). lepton = μ− (index 1).
    zsign = torch.where(Pz >= 0, torch.ones_like(Pz), -torch.ones_like(Pz))
    pbeam = math.sqrt(6500.0 * 6500.0 - 0.93827208816 * 0.93827208816)
    zero = torch.zeros_like(Pz)
    p1x, p1y, p1z = _boost_unit(zero, zero, zsign * pbeam, 6500.0)
    p2x, p2y, p2z = _boost_unit(zero, zero, -zsign * pbeam, 6500.0)
    lx, ly, lz = _boost_unit(px[..., 1], py[..., 1], pz[..., 1], E[..., 1])

    def _unit(ax, ay, az):
        n = torch.sqrt((ax * ax + ay * ay + az * az).clamp_min(eps * eps))
        return ax / n, ay / n, az / n

    def _cross(ax, ay, az, bx_, by_, bz_):
        return (ay * bz_ - az * by_, az * bx_ - ax * bz_, ax * by_ - ay * bx_)

    fx, fy, fz = _unit(p1x - p2x, p1y - p2y, p1z - p2z)             # csFrame
    yx, yy, yz = _unit(*_cross(p1x, p1y, p1z, -p2x, -p2y, -p2z))    # csYaxis
    xx, xy, xz = _unit(*_cross(yx, yy, yz, fx, fy, fz))            # csXaxis
    costheta = fx * lx + fy * ly + fz * lz
    crx, cry, crz = _cross(fx, fy, fz, lx, ly, lz)
    sintheta = torch.sqrt((crx * crx + cry * cry + crz * crz).clamp_min(eps * eps))
    sinphi = (yx * lx + yy * ly + yz * lz) / sintheta
    cosphi = (xx * lx + xy * ly + xz * lz) / sintheta
    return torch.stack(
        [yll, torch.log(ptll) - torch.log(mll), cosPhill, sinPhill, costheta,
         sinphi, cosphi],
        dim=-1)


# ---------------------------------------------------------------------------
# Mixture MLP — unchanged
# ---------------------------------------------------------------------------


class MixtureMLP(nn.Module):
    """Softmax mixture fractions (f_bkg_0, …, f_bkg_{n−1}, f_signal) as a
    function of the conditioning — SIGNAL LAST (``f[:, -1]``; for the
    historical degree-1 Bernstein this is the (f_0, f_1, f_s) 3-way softmax).
    ``n_frac`` = number of background components + 1. With ``exp_slope`` the
    final layer carries one extra UNCONSTRAINED output: the
    conditioning-dependent dimensionless slope s(c) = λ(c)·width of the
    exponential background (sharing the MLP body with the fractions, so every
    "background parameters" collection — optimiser group, Fisher, freezing —
    covers it automatically)."""

    def __init__(self, n_input: int = N_MUON_KIN, hidden: int = 32,
                 n_layers: int = 2, n_frac: int = 3, exp_slope: bool = False):
        super().__init__()
        self.n_frac = int(n_frac)
        self.exp_slope = bool(exp_slope)
        layers: list[nn.Module] = []
        d_in = n_input
        for _ in range(n_layers):
            layers += [nn.Linear(d_in, hidden), nn.GELU()]
            d_in = hidden
        layers.append(nn.Linear(d_in, self.n_frac + (1 if exp_slope else 0)))
        self.net = nn.Sequential(*layers)
        if exp_slope:
            # The slope head s(c) = λ(c)·width (dimensionless, signed,
            # UNconstrained — no softmax) shares the MLP body with the
            # fractions. Zero-init its final-layer row so s(c) ≡ 0 at init
            # (uniform background), matching the fractions' near-uniform
            # softmax start.
            with torch.no_grad():
                self.net[-1].weight[-1].zero_()
                self.net[-1].bias[-1].zero_()

    def forward(self, y_std: torch.Tensor) -> torch.Tensor:
        out = self.net(y_std)
        return F.softmax(out[..., :self.n_frac], dim=-1)

    def forward_with_slope(self, y_std: torch.Tensor):
        """(fractions [.., n_frac], slope s(c) [..]) in ONE body evaluation —
        the exp-background hot path."""
        out = self.net(y_std)
        return F.softmax(out[..., :self.n_frac], dim=-1), out[..., -1]


class ThetaNet(nn.Module):
    """Per-muon calibration parameters (A, e, M, a, c) as a CONTINUOUS function
    of the muon (η, φ), replacing the η-binned θ. Input features (η, cosφ, sinφ),
    5 outputs. The (A, e, M) outputs are scaled by fixed references so the net
    outputs sit at O(1) (A,e ~ 1e-3, M ~ 1e-5); (a, c) are the qop-resolution
    variance coefficients in their O(1) units (SMEAR_VAR_SCALE_* applied
    downstream). The final layer is ZERO-INITIALISED, so the net outputs 0 at
    init — the binned θ=0 init (no scale/smear correction) — EXCEPT the smear
    (a, c) outputs, whose final-layer BIAS is set to ``smear_bias_init`` so the
    raw (a, c) start at that constant for every (η, φ). This is required for the
    'square' reparam, where a zero smear init is a dead saddle (∂effective/∂raw
    = 0); see ``SMEAR_SQUARE_INIT_RAW``. Default 0.0 reproduces the original
    zero-output init (correct for 'linear'/'softplus')."""

    def __init__(self, hidden: int = 32, n_layers: int = 2,
                 scale_ref=THETA_SCALE_REF, smear_bias_init: float = 0.0):
        super().__init__()
        layers: list[nn.Module] = []
        d_in = 3  # (η, cosφ, sinφ)
        for _ in range(max(1, n_layers)):
            layers += [nn.Linear(d_in, hidden), nn.GELU()]
            d_in = hidden
        last = nn.Linear(d_in, N_THETA_SCALE + N_THETA_SMEAR)  # 5 = (A,e,M,a,c)
        nn.init.zeros_(last.weight)
        nn.init.zeros_(last.bias)
        if smear_bias_init != 0.0:
            # Smear (a, c) outputs are the trailing N_THETA_SMEAR. With the
            # weight zeroed, the raw (a, c) output equals this bias for ALL
            # inputs at init → uniform nonzero raw smear, off the square saddle.
            with torch.no_grad():
                last.bias[N_THETA_SCALE:].fill_(float(smear_bias_init))
        layers.append(last)
        self.net = nn.Sequential(*layers)
        self.register_buffer(
            "scale_ref", torch.tensor(list(scale_ref), dtype=torch.float32))

    def forward(self, eta_pm: torch.Tensor, phi_pm: torch.Tensor):
        """``eta_pm``, ``phi_pm``: ``[B, 2]``. Returns ``(AeM [B,2,3], ac [B,2,2])``
        — per-muon physical (A, e, M) and O(1) (a, c)."""
        feat = torch.stack(
            [eta_pm, torch.cos(phi_pm), torch.sin(phi_pm)], dim=-1)  # [B,2,3]
        out = self.net(feat)                                        # [B,2,5]
        return out[..., :N_THETA_SCALE] * self.scale_ref, out[..., N_THETA_SCALE:]


@contextlib.contextmanager
def _freeze_param_grads(params):
    """Temporarily set ``requires_grad=False`` on ``params`` for the duration
    of the ``with`` block, restoring the previous state afterwards.

    Used to exclude a sub-module's *parameters* from a backward while still
    propagating gradient to that sub-module's *inputs* — a frozen layer still
    passes gradient to earlier tensors. We only flip (and later restore)
    params that were ``requires_grad=True`` on entry, so this is a no-op for
    already-frozen params (e.g. during Fisher-info, where the flow is frozen).
    """
    changed = [p for p in params if p.requires_grad]
    for p in changed:
        p.requires_grad_(False)
    try:
        yield
    finally:
        for p in changed:
            p.requires_grad_(True)


# ---------------------------------------------------------------------------
# Main model
# ---------------------------------------------------------------------------


class JpsiMassMixtureModel(nn.Module):
    """Unbinned J/ψ mass-fit model — two-stage continuity design (the flow models
    the nominal shape p₀(m|muon_kin); θ enters analytically in stage 2)."""

    def __init__(
        self,
        m_lo: float,
        m_hi: float,
        mll_log_scale: float,
        # Stats for standardising the flow's mass input + conditioning.
        # Passed in as plain floats / tensors; stored as buffers.
        mll_mean: float = 0.0,
        mll_std: float = 1.0,
        y_event_mean: torch.Tensor | None = None,
        y_event_std_tensor: torch.Tensor | None = None,
        muon_kin_mean: torch.Tensor | None = None,
        muon_kin_std_tensor: torch.Tensor | None = None,
        # Flow / MLP hyperparameters.
        flow_arch: str = "gf",
        flow_n_transforms: int = 5,
        flow_hidden_features: int = 128,
        flow_n_hidden_layers: int = 3,
        flow_gf_components: int = 8,
        flow_nsf_bins: int = 8,
        compact_learn_weights: bool = False,
        compact_layer: str = "logistic",
        bernstein_degree: int = 16,
        mixture_layer: str = "gaussian",
        nce_quad_nodes: int = 64,
        mlp_hidden: int = 32,
        mlp_n_layers: int = 2,
        n_eta_bins: int = N_ETA_BINS,
        # Debug toggle: drop the residual-smearing kernel + σ_qop_pm
        # conditioning entirely. ``theta_smear`` stays as a Parameter
        # (for state_dict shape consistency) but is unused; trainer is
        # expected to exclude it from the optimizer.
        smearing_enabled: bool = True,
        # Symmetric toggle for the scale: drop the T_scale forward-fold + the
        # θ_scale_pm conditioning entirely. ``theta_scale`` stays as a (zero,
        # inert) Parameter; trainer excludes it from the optimizer. With both
        # scale and smearing disabled only the flow + MLP (background) train.
        scale_enabled: bool = True,
        # DEPRECATED / inert. The qop→pt inversion is now ``pt = |sinθ/qop|``
        # with only the ``qop=0`` pole guarded (``QOP_EPS``): pt is a magnitude,
        # so a qop sign flip from a large kick is treated as the PHYSICAL charge
        # mis-reconstruction it is, not floored away. This argument is accepted
        # for checkpoint/back-compat but no longer affects the inversion.
        qop_floor_frac: float = 0.0,
        # Which per-bin smear terms to *fit*: "both" (a and c), "a" (constant
        # term only), or "c" (∝1/pt term only). The constant a and the c·k
        # term are nearly degenerate over the narrow J/ψ pt range, so fitting
        # both per η-bin is ill-posed and yields the unphysical bin-to-bin
        # zig-zag. Fitting one removes the degeneracy; the non-fitted term is
        # zeroed (``smear_param_mask``) so it contributes exactly 0 to the width
        # factor s and receives no gradient (inert).
        smear_fit_params: str = "both",
        # Which per-η-bin SCALE terms to fit — a subset of "AeM". A (constant
        # δqop) and e (∝1/pt) are NEARLY DEGENERATE over the narrow J/ψ pt range
        # (A's uniform m_ll scaling vs e's ∝(k₊+k₋) shift are ~collinear), so
        # fitting both from J/ψ alone is ill-posed — the fit slides into large
        # opposite-sign (A, e). Default "AM" drops the degenerate e (the J/ψ-
        # identifiable subset: constant scale A + charge-odd sagitta M). The
        # dropped term is zeroed (scale_param_mask) → 0 advection, no gradient.
        scale_fit_params: str = "AM",
        # Number of Euler steps integrating the smear's probability-flow ODE in
        # the density (``_continuity_logp``). 1 = first-order (single score
        # displacement); more steps integrate the score-driven flow more finely
        # (smaller per-step Jacobian → more robust + accurate) at a higher
        # nested-autograd cost. Frozen p₀ score per step (exact diffusion to
        # first order in V).
        smear_flow_steps: int = 1,
        # Smear operator for the continuity density. 'pf_ode' (default): the
        # deterministic probability-flow ODE ``y ← y − (V/2N)·∂_m log p_0(y)``
        # (N = smear_flow_steps), with Jacobian via jacobian_form. Cheap but
        # diverges from the true Gaussian convolution at large V/σ² — for
        # Gaussian p_0 the PF-ODE gives σ_out = σ·exp(V/2σ²), while the true
        # convolution gives σ_out = √(σ²+V); they agree to leading order but
        # the PF-ODE over-broadens exponentially at large V (e.g. +50% σ at
        # V/σ² ≈ 2). 'gh_convolution': EXACT stochastic Gaussian convolution
        # ``p_θ(x|c) = E_ε[p_0(m'(ε)|c)/|G'|]`` via Gauss-Hermite quadrature
        # over ε ~ N(0,1) (``n_gh_nodes`` nodes), with x = m' + s_adv(m') +
        # √V(m')·ε and G'(m') = 1 + s_adv'(m') + (V'/(2√V))·ε computed via
        # autograd. Matches the per-muon qop fold operator (which IS Gaussian
        # convolution); the closure-target curve (flow at injected θ) overlaps
        # the pseudo-data by construction at large V. Requires V ≥ 0 (clamps
        # internally); use with smear_param_form='softplus' to guarantee it
        # without clamping.
        smear_operator: str = "pf_ode",
        n_gh_nodes: int = 8,
        # Positivity reparameterisation for θ_smear. 'linear' (default): the
        # raw θ_smear is the O(1) coefficient directly (physical (a, c) =
        # θ·SMEAR_VAR_SCALE, signed — supports both broadening V>0 and the
        # 'unsmear' V<0 region). 'softplus': constrains each of (a, c)
        # INDIVIDUALLY to ≥ 0 via physical = softplus(θ)·SMEAR_VAR_SCALE; the
        # per-η mask is applied AFTER softplus so frozen params (or 'a'-/'c'-
        # only modes) are EXACTLY zero in the per-muon σ_qop and the
        # transformations. Useful when you want to defend against the
        # negative-c drift (issue #2) by construction, at the cost of losing
        # the two-sided fit (the model can no longer represent MC that is too
        # broad vs data).
        # 'square': physical = θ²·SMEAR_VAR_SCALE ≥ 0 — same individual-
        # positivity guarantee as softplus, but with better convergence: θ=0
        # maps to physical=0 EXACTLY (identity init, like 'linear'; no softplus
        # ln2 offset / saturation knee), and ∂physical/∂θ = 2·SMEAR_VAR_SCALE·θ
        # grows away from zero rather than saturating, so a fixed lr isn't
        # throttled in the small-variance regime. Caveats: θ↔−θ are degenerate
        # (the map is even) and ∂physical/∂θ → 0 at θ=0, so the RAW-θ covariance
        # paths (binned observed Fisher / net-weight empirical Fisher / boot-
        # strap) become singular for any smear bin pinned at zero — take the
        # smear σ from the reparam-robust --output-fisher there instead.
        smear_param_form: str = "linear",
        # Per-event normalisation correction for the transformed-flow density.
        # The forward map T_θ is NOT boundary-preserving on [m_lo, m_hi]:
        # broadening T pushes some probability mass outside the window, so
        # Z(θ;c) = ∫_window p_θ(x|c) dx < 1 and the bare `log p_θ(x)` carries a
        # `-log Z` bias per event that always pulls the fit toward smaller
        # broadening (the bias was measured at +0.34 per event at the truth in
        # forward |η| bins for a c=5e-5 injection). Three modes:
        #
        # "none" (default): no correction (current behaviour preserved).
        # "linear": leading-order boundary expansion `1 - Z ≈ p_0(m_lo)·(m_lo
        #   − T(m_lo))_+ + p_0(m_hi)·(T(m_hi) − m_hi)_+` — 2 boundary forward
        #   evals + 2 flow evals per event. Valid for small V.
        # "flow_cdf": exact via Z = F_0(T⁻¹(m_hi)|c) − F_0(T⁻¹(m_lo)|c), where
        #   F_0 is the flow's CDF (the GF's monotonic transform composed with
        #   the standard normal CDF Φ). Inverts T at the boundaries (~2 extra
        #   fixed-point inversions per event) and evaluates the CDF at the
        #   preimages. Exact up to inversion discretisation.
        norm_correction: str = "none",
        # Background mixture. True (default): the data branch's per-event NLL
        # is the full f_data(c)-weighted signal/Bernstein mixture (the model
        # for real data, which has genuine non-resonant background). False:
        # the data branch reduces to pure signal (NLL = −log p_signal); the
        # MLP `f_data` is bypassed entirely and its parameters are excluded
        # from the optimiser in stage 2. Use for validation closures (truth
        # f_bkg = 0 by construction) to remove the bkg ↔ smear degeneracy
        # where the MLP grows f_bkg in forward |η| bins to absorb tail
        # events the signal model can't broaden into.
        background_enabled: bool = True,
        # Background fractions as GLOBAL learnable scalars (constant across all
        # events / conditioning) instead of the per-event MLP f(c). The fractions
        # are softmax(global_logits) → [n_bkg_comp+1] (signal LAST), broadcast to
        # every event. Diagnostic mode to remove the MLP's per-event flexibility
        # (and its overfitting / smear↔bkg degeneracy sink): with the background a
        # pure pair of global numbers, an injected GLOBAL (f0, f1) must be
        # recovered exactly — isolating the signal-density / window-normalisation.
        bkg_global: bool = False,
        # FIT-TIME reco pt selection cuts (GeV). At fixed conditioning pt ∝ m,
        # so a fixed pt cut C on a pt quantity q (observed ratio R=q_obs/m_obs)
        # forbids m < C/R = C·m_obs/q_obs — a per-event LOWER mass edge
        # m_min(c). When set, the signal+background are normalised over the
        # per-event window [max(m_lo, m_min(c)), m_hi] instead of the fixed
        # [m_lo, m_hi]. None (default) → no extra edge → exact legacy behaviour.
        # Used when the cuts are deferred from shard production to the fit (so
        # the flow is trained data-constrained across the forbidden region).
        # The |η| cut is conditioning-fixed (no mass edge) and is applied as an
        # event mask in the loader, not here.
        fit_ptll_min: "float | None" = None,
        fit_pt_lead_min: "float | None" = None,
        fit_pt_both_min: "float | None" = None,
        # FLOW-stage mass window override (defaults to the model/fit window
        # m_lo/m_hi). The compact/nce/dcb/ege flows are defined and
        # normalised on THIS window (their standardised [a, b]); the fit's
        # truncated likelihood then renormalises the frozen flow over the
        # (possibly tighter) FIT window m_lo/m_hi via the window-Z machinery
        # exactly as for any transformed window — so the per-stage windows
        # compose without special cases as long as fit ⊆ flow.
        flow_m_lo: "float | None" = None,
        flow_m_hi: "float | None" = None,
        # Background model on the observed window: 'bernstein' = positive
        # degree-`bkg_degree` Bernstein mixture (bkg_degree+1 fractions from
        # the MLP; degree 1 is the historical default); 'exp' = window-
        # normalised exponential with ONE fraction plus a CONDITIONING-
        # DEPENDENT dimensionless slope s(c) = lambda(c)*width — an extra
        # unconstrained MLP head sharing the body with the fractions.
        bkg_model: str = "bernstein",
        bkg_degree: int = 1,
        # Smear-Jacobian formula for the continuity density. 'softlog' (default):
        # autograd-derived G' = dx/dm' of the actual forward map, with a C¹
        # tangent extension below SMEAR_GP_FLOOR so the optimiser is pulled BACK
        # from G' < floor instead of sliding into a flat-NLL basin. 'exp':
        # frozen-score continuous-flow approximation log G' = log(1+s_adv'(m'))
        # − V·∂²_m log p₀(m')/2 — always finite, no floor, but approximates a
        # *different* operator (extra "score constant along the trajectory"
        # assumption) — see _continuity_logp docstring caveats.
        jacobian_form: str = "softlog",
        # θ parameterisation. 'binned' (default): per-η-bin (A,e,M,a,c) tables
        # indexed by the muon's η-bin. 'mlp': a small ThetaNet maps each muon's
        # (η, φ) → (A,e,M,a,c) CONTINUOUSLY (trained in stage 2 like the
        # background MLP). The binned tables stay registered but inert in 'mlp'.
        theta_mode: str = "binned",
        theta_mlp_hidden: int = 32,
        theta_mlp_layers: int = 2,
        # 'binned2d' only: number of uniform φ bins over [-π, π). The θ tables
        # then have n_eta_bins × n_phi_bins CELLS (η outer, φ inner), indexed
        # per-muon by (η-bin, φ-bin). 1 ⇒ identical to plain 'binned'.
        n_phi_bins: int = 16,
        # Flow + background-MLP conditioning basis. 'muon_kin' (default): the
        # leak-free per-muon (η±, cosφ±, sinφ±, ρ). 'event_level': the dilepton
        # vars (yll, ln(ptll/mll), cosPhill, sinPhill, cosθ*, sinφ*, cosφ*) —
        # all DIMENSIONLESS, hence invariant under the common-pt-dilation mass
        # direction (the leak-free criterion); residual e/M/smear shifts
        # propagate through the per-node recompute (see _cond_from_muons).
        # Both are 7-dim (no flow-shape change); the basis must match between
        # stage-1 and stage-2.
        cond_basis: str = "muon_kin",
        # Degeneracy-whitening preconditioner (see _WhitenGradFn). When True,
        # the gradient of the (A, e) and (a, c) pairs is rotated into a
        # decorrelated basis built from the per-η-bin curvature moments, so the
        # near-collinear pairs can be fit JOINTLY without the tilted-valley
        # stalling that forces dropping a degenerate term. Backward-only: the
        # forward map / likelihood / Fisher are unchanged. Only acts when BOTH
        # members of a pair are fit (scale_fit_params ⊇ {A,e} / smear_fit_params
        # == 'both'); otherwise there is no degeneracy to whiten and it is inert.
        theta_whiten: bool = False,
        theta_whiten_max_rho: float = 0.99,
        # Per-η-bin curvature moment SUMS [n_eta, 4] = (N, Σk, Σk², Σk⁴) from
        # JpsiMassPreprocStats.k_moments. Required to build the whitening; if
        # None while theta_whiten=True, whitening is disabled with a warning.
        k_moments=None,
    ):
        super().__init__()
        self.smearing_enabled = bool(smearing_enabled)
        self.scale_enabled = bool(scale_enabled)
        self.qop_floor_frac = float(qop_floor_frac)
        self.smear_flow_steps = max(1, int(smear_flow_steps))
        if smear_operator not in ("pf_ode", "gh_convolution", "gh_convolution_qop"):
            raise ValueError(
                f"smear_operator must be 'pf_ode', 'gh_convolution', or "
                f"'gh_convolution_qop'; got {smear_operator!r}")
        self.smear_operator = str(smear_operator)
        self.n_gh_nodes = max(2, int(n_gh_nodes))
        if jacobian_form not in ("softlog", "exp"):
            raise ValueError(
                f"jacobian_form must be 'softlog' or 'exp'; got {jacobian_form!r}")
        self.jacobian_form = str(jacobian_form)
        if smear_param_form not in ("linear", "softplus", "square"):
            raise ValueError(
                f"smear_param_form must be 'linear', 'softplus', or 'square'; "
                f"got {smear_param_form!r}")
        self.smear_param_form = str(smear_param_form)
        if norm_correction not in ("none", "linear", "flow_cdf"):
            raise ValueError(
                f"norm_correction must be 'none', 'linear', or 'flow_cdf'; "
                f"got {norm_correction!r}")
        self.norm_correction = str(norm_correction)
        self.background_enabled = bool(background_enabled)
        # Fit-time reco pt cuts → per-event lower mass edge (None = inactive).
        self._fit_ptll_min = (float(fit_ptll_min) if fit_ptll_min else None)
        self._fit_pt_lead_min = (float(fit_pt_lead_min) if fit_pt_lead_min else None)
        self._fit_pt_both_min = (float(fit_pt_both_min) if fit_pt_both_min else None)
        if theta_mode not in ("binned", "mlp", "binned2d"):
            raise ValueError(
                f"theta_mode must be 'binned', 'binned2d', or 'mlp', "
                f"got {theta_mode!r}")
        self.theta_mode = str(theta_mode)
        # θ-table cell layout: plain 'binned'/'mlp' use the per-η table (n_phi=1
        # conceptually); 'binned2d' tiles each η-bin into n_phi_bins φ-cells,
        # ordered η-outer/φ-inner (matching the MLP grid-output convention).
        self.theta_grid = (self.theta_mode == "binned2d")
        self.n_eta_bins = int(n_eta_bins)
        self.n_phi_bins = int(n_phi_bins) if self.theta_grid else 1
        n_theta_cells = self.n_eta_bins * self.n_phi_bins
        if cond_basis not in ("muon_kin", "event_level"):
            raise ValueError(
                f"cond_basis must be 'muon_kin' or 'event_level'; got {cond_basis!r}")
        self.cond_basis = str(cond_basis)
        # event_level propagates the per-muon qop kick to ALL conditioning vars,
        # which requires reconstructing per-muon pt per smear node — only the
        # qop operator does that. The mass-space operators transport smear in
        # mass space (no per-muon smeared pt), so they cannot carry the smear
        # into the (pt-dependent) event vars; require the qop operator there.
        if (self.cond_basis == "event_level" and self.smearing_enabled
                and self.smear_operator != "gh_convolution_qop"):
            raise ValueError(
                "cond_basis='event_level' with smearing enabled requires "
                "--smear-operator gh_convolution_qop (the mass-space operators "
                "cannot propagate the smear to the event-level conditioning); "
                "use gh_convolution_qop, or --disable-smearing for scale-only.")
        self.flow_arch = str(flow_arch)

        # Per-bin smear fit mask: which of (a, c) float. The frozen column gets
        # no gradient (held at init) and is not perturbed in the MC sampling.
        if smear_fit_params not in ("both", "a", "c"):
            raise ValueError(
                f"smear_fit_params must be 'both', 'a', or 'c'; got {smear_fit_params!r}"
            )
        self.smear_fit_params = str(smear_fit_params)
        _mask = {"both": [1.0, 1.0], "a": [1.0, 0.0], "c": [0.0, 1.0]}[smear_fit_params]
        # Non-persistent: reconstructed from smear_fit_params at __init__, so
        # checkpoints without this buffer still load.
        self.register_buffer(
            "smear_param_mask", torch.tensor(_mask, dtype=torch.float32),
            persistent=False,
        )
        # Per-bin scale fit mask: which of (A, e, M) float (breaks the A/e
        # degeneracy). The dropped term gets no gradient (multiply by 0) → held
        # inert at 0. Non-persistent (reconstructed at __init__).
        if not scale_fit_params or any(ch not in "AeM" for ch in scale_fit_params):
            raise ValueError(
                f"scale_fit_params must be a non-empty subset of 'AeM'; got "
                f"{scale_fit_params!r}")
        self.scale_fit_params = str(scale_fit_params)
        _csmask = [1.0 if ch in scale_fit_params else 0.0 for ch in ("A", "e", "M")]
        self.register_buffer(
            "scale_param_mask", torch.tensor(_csmask, dtype=torch.float32),
            persistent=False,
        )

        # Two-stage continuity design: the flow models only the nominal shape
        # p₀(m|muon_kin) at θ=0 — it never conditions on θ (the θ-dependence is
        # supplied analytically in stage 2, see ``data_nll_continuity``).
        # ``compact``: uniform-base compact 1-D flow on the standardised mass
        # window with C²-matched analytic tails (see compact_flow.py). Exactly
        # normalised over the window (Z=1, no out-of-window mass gauge freedom →
        # no spurious far-tail structure), C∞ interior, and smoothly evaluable
        # just outside the window (the scale un-kick / smear / norm-Z all reach
        # there). The window edges are the standardised [m_lo, m_hi].
        self.flow_is_compact = (self.flow_arch == "compact")
        # ``nce``: classifier-vs-uniform (noise-contrastive) density — an
        # unconstrained MLP logit trained with BCE against window-uniform
        # twins in stage 1 (see nce_density.py). Absolute window-normalised
        # density in expectation (no Z gauge drift, truncation exact by
        # construction); the CDF (window Z) is per-event Gauss-Legendre
        # quadrature instead of analytic.
        self.flow_is_nce = (self.flow_arch == "nce")
        # ``dcb``: double-sided Crystal Ball with all six parameters
        # (μ, σ, α_L, n_L, α_R, n_R) conditional on c through an MLP
        # (see dcb_density.py). Fully analytic density AND CDF,
        # window-normalised by construction (Z ≡ 1 — no out-of-window gauge),
        # power-law tails evaluable on all of ℝ for the operator probes.
        # C¹ at the two core/tail junctions (score continuous; curvature
        # jumps at two isolated points).
        self.flow_is_dcb = (self.flow_arch == "dcb")
        # ``ege``: Gaussian core + C¹-MATCHED exponential tails (ExpGaussExp;
        # see dcb_density.EGEDensity) — like dcb but with only (μ, σ, α_L,
        # α_R) conditional parameters: the tail slopes are FORCED to α by the
        # C⁰+C¹ matching (the maximum possible for a log-linear tail). Same
        # window-normalised, analytic-CDF, full-ℝ-tails design.
        self.flow_is_ege = (self.flow_arch == "ege")
        # ``mixture``: normal-base, FULL-support 1-D flow — a stack of Gaussian-
        # or logistic-mixture-CDF Gaussianisation layers (see mixture_flow.py).
        # Analytic density + CDF (→ exact cheap window-Z) and analytic Jacobian
        # (torch.compile-traceable, unlike gf). Support = ℝ with LEARNED tails
        # (no compact matched tails); normal base, so it uses the same fp64-stable
        # window-Z path as gf/nsf (the window mass can be < 1 / drift small).
        self.flow_is_mixture = (self.flow_arch == "mixture")
        self._flow_m_lo_f = float(flow_m_lo if flow_m_lo is not None else m_lo)
        self._flow_m_hi_f = float(flow_m_hi if flow_m_hi is not None else m_hi)
        if self.flow_is_compact:
            a_std = (self._flow_m_lo_f - float(mll_mean)) / float(mll_std)
            b_std = (self._flow_m_hi_f - float(mll_mean)) / float(mll_std)
            self.flow = CompactMatchedFlow(
                n_cond=N_MUON_KIN, a=a_std, b=b_std,
                hidden_features=flow_hidden_features,
                n_layers=flow_n_hidden_layers,
                n_components=flow_gf_components,
                n_transforms=flow_n_transforms,   # composed depth, as for gf
                # learn_weights is a logistic-layer option; bernstein layers
                # have no mixture weights (sanitise rather than raise so old
                # configs combine freely with --compact-layer bernstein).
                learn_weights=(compact_learn_weights
                               and compact_layer == "logistic"),
                layer_type=compact_layer,
                bernstein_degree=bernstein_degree,
            )
        elif self.flow_is_nce:
            a_std = (self._flow_m_lo_f - float(mll_mean)) / float(mll_std)
            b_std = (self._flow_m_hi_f - float(mll_mean)) / float(mll_std)
            self.flow = NCEDensity(
                n_cond=N_MUON_KIN, a=a_std, b=b_std,
                hidden_features=flow_hidden_features,
                n_layers=flow_n_hidden_layers,
                quad_nodes=nce_quad_nodes,
            )
        elif self.flow_is_dcb:
            a_std = (self._flow_m_lo_f - float(mll_mean)) / float(mll_std)
            b_std = (self._flow_m_hi_f - float(mll_mean)) / float(mll_std)
            self.flow = DCBDensity(
                n_cond=N_MUON_KIN, a=a_std, b=b_std,
                hidden_features=flow_hidden_features,
                n_layers=flow_n_hidden_layers,
            )
        elif self.flow_is_ege:
            a_std = (self._flow_m_lo_f - float(mll_mean)) / float(mll_std)
            b_std = (self._flow_m_hi_f - float(mll_mean)) / float(mll_std)
            self.flow = EGEDensity(
                n_cond=N_MUON_KIN, a=a_std, b=b_std,
                hidden_features=flow_hidden_features,
                n_layers=flow_n_hidden_layers,
            )
        elif self.flow_is_mixture:
            # Full-support (ℝ), no window edges — the normal base supplies the
            # tails. Depth = flow_n_transforms, components/layer = flow_gf_components
            # (shared with gf for parameter parity).
            self.flow = MixtureFlow(
                n_cond=N_MUON_KIN,
                hidden_features=flow_hidden_features,
                n_layers=flow_n_hidden_layers,
                n_components=flow_gf_components,
                n_transforms=flow_n_transforms,
                layer_type=mixture_layer,
            )
        else:
            flow_inner = build_flow(
                n_features=1,
                n_cond=N_MUON_KIN,
                n_transforms=flow_n_transforms,
                hidden_features=flow_hidden_features,
                n_hidden_layers=flow_n_hidden_layers,
                architecture=self.flow_arch,
                gf_components=flow_gf_components,
                nsf_bins=flow_nsf_bins,
            )
            self.flow = FlowWithLogProb(flow_inner)

        # Background model: 'bernstein' (degree-n positive mixture, n+1
        # components) or 'exp' (one fraction + a fitted dimensionless slope).
        # Signal fraction is always the LAST softmax output.
        if bkg_model not in ("bernstein", "exp"):
            raise ValueError(f"bkg_model must be 'bernstein' or 'exp', got {bkg_model.__repr__()}")
        if bkg_model == "bernstein" and int(bkg_degree) < 1:
            raise ValueError("bkg_degree must be >= 1")
        self.bkg_model = str(bkg_model)
        self.bkg_degree = int(bkg_degree)
        self.n_bkg_comp = 1 if bkg_model == "exp" else self.bkg_degree + 1
        self.bkg_global = bool(bkg_global)
        # Background-fraction MLP conditions on the same kinematics as the
        # flow (muon_kin), minus the nuisances. ``bkg_global`` replaces it with a
        # single global logit vector (constant fractions for all events).
        if self.bkg_global:
            if bkg_model == "exp":
                raise ValueError("bkg_global is only supported for the bernstein "
                                 "background (global f0..f_{n} fractions)")
            # softmax over [n_bkg_comp+1] (signal LAST) → constant fractions;
            # init at uniform (logits 0) so f_sig and each f_bkg start at 1/(n+1).
            self.bkg_global_logits = nn.Parameter(
                torch.zeros(self.n_bkg_comp + 1, dtype=torch.float32))
            self.mlp = None
        else:
            self.mlp = MixtureMLP(
                n_input=N_MUON_KIN, hidden=mlp_hidden, n_layers=mlp_n_layers,
                n_frac=self.n_bkg_comp + 1, exp_slope=(bkg_model == "exp")
            )

        # Learnable nuisances.
        self.theta_scale = nn.Parameter(
            torch.zeros(n_theta_cells, N_THETA_SCALE, dtype=torch.float32)
        )
        # Default RAW θ_smear init. 'linear'/'softplus' start at 0 (the historical
        # identity init); 'square' starts at SMEAR_SQUARE_INIT_RAW because raw=0
        # is a dead saddle for that form (∂effective/∂raw = 0 → no smear gradient,
        # observed as θ_smear frozen at 0 in a real fit). Applied to BOTH the
        # binned table (here) and the ThetaNet smear bias (below) so MLP θ — which
        # has no --init-theta-{a,c} override path — is also seeded off the saddle.
        smear_init_raw = (SMEAR_SQUARE_INIT_RAW
                          if self.smear_param_form == "square" else 0.0)
        # θ_smear are signed per-η-bin qop-resolution VARIANCE coefficients
        # (a, c): σ²_qop = a + c·k² (two-sided). They drive BOTH the per-muon qop
        # fold (validation) and the mass-density stretch (density), consistently.
        # Init: 0 (linear/softplus) → σ²_qop=0 identity; SMEAR_SQUARE_INIT_RAW
        # (square) → small nonzero, off the saddle. The binned --init-theta-{a,c}
        # CLI overrides this downstream when explicitly set.
        self.theta_smear = nn.Parameter(
            torch.full((n_theta_cells, N_THETA_SMEAR), smear_init_raw,
                       dtype=torch.float32)
        )
        # 'binned2d': uniform φ-bin edges over [-π, π) for the per-muon cell index.
        self.register_buffer(
            "_theta_phi_edges",
            torch.linspace(-math.pi, math.pi, self.n_phi_bins + 1,
                           dtype=torch.float32),
            persistent=False)
        # 'mlp' θ: a small net maps each muon's (η, φ) → (A,e,M,a,c) continuously
        # (scale zero-init → 0; smear bias = smear_init_raw). Replaces the binned
        # tables above (which stay registered but inert). Trained in stage 2 like
        # the background MLP.
        self.theta_net = (
            ThetaNet(hidden=theta_mlp_hidden, n_layers=theta_mlp_layers,
                     smear_bias_init=smear_init_raw)
            if self.theta_mode == "mlp" else None
        )

        # Degeneracy-whitening preconditioner (backward-only; see _WhitenGradFn).
        self._build_whitening(
            bool(theta_whiten), float(theta_whiten_max_rho), k_moments, n_eta_bins)

        # Buffers — Bernstein window, density-rescale, standardisation stats.
        self.register_buffer("m_lo", torch.tensor(float(m_lo)))
        self.register_buffer("m_hi", torch.tensor(float(m_hi)))
        # Python-float twins of the window buffers for the hot paths:
        # float(<buffer>) is a hidden .item() — a per-call GPU sync in eager
        # mode AND a dynamo graph break under --fit-compile.
        self._m_lo_f = float(m_lo)
        self._m_hi_f = float(m_hi)
        self.register_buffer("mll_log_scale", torch.tensor(float(mll_log_scale)))
        self.register_buffer("mll_mean_buf", torch.tensor(float(mll_mean)))
        self.register_buffer("mll_std_buf", torch.tensor(float(mll_std)))

        def _buf(t, n):
            if t is None:
                return torch.zeros(n, dtype=torch.float32), torch.ones(n, dtype=torch.float32)
            return torch.as_tensor(t, dtype=torch.float32), None

        if y_event_mean is None:
            y_event_mean = torch.zeros(N_Y_EVENT)
        if y_event_std_tensor is None:
            y_event_std_tensor = torch.ones(N_Y_EVENT)
        if muon_kin_mean is None:
            muon_kin_mean = torch.zeros(N_MUON_KIN)
        if muon_kin_std_tensor is None:
            muon_kin_std_tensor = torch.ones(N_MUON_KIN)
        self.register_buffer(
            "y_event_mean", torch.as_tensor(y_event_mean, dtype=torch.float32)
        )
        self.register_buffer(
            "y_event_std", torch.as_tensor(y_event_std_tensor, dtype=torch.float32)
        )
        self.register_buffer(
            "muon_kin_mean", torch.as_tensor(muon_kin_mean, dtype=torch.float32)
        )
        self.register_buffer(
            "muon_kin_std", torch.as_tensor(muon_kin_std_tensor, dtype=torch.float32)
        )

    # ------------------------------------------------------------------
    # Per-event helpers
    # ------------------------------------------------------------------

    def _build_whitening(self, enabled, max_rho, k_moments, n_eta) -> None:
        """Build the (A,e) and (a,c) gradient-whitening matrices from the per-η-bin
        curvature moments ``k_moments`` [n_eta, 4] = (N, Σk, Σk², Σk⁴). Registers
        per-bin ``[n_eta, 2, 2]`` and global ``[2, 2]`` buffers for each pair and
        sets the per-pair active flags. The buffers are identity (inert) when
        whitening is disabled, when the pair is not fully fit (no degeneracy to
        break), or when moments are unavailable (older stats.json)."""
        cm = self.scale_param_mask
        sm = self.smear_param_mask
        scale_pair_fit = bool(cm[0] > 0 and cm[1] > 0)   # both A and e fit
        smear_pair_fit = bool(sm[0] > 0 and sm[1] > 0)   # both a and c fit
        self.theta_whiten = bool(enabled)
        self.theta_whiten_max_rho = float(max_rho)

        eye = torch.eye(2, dtype=torch.float32)
        sW_b = eye.expand(n_eta, 2, 2).clone()
        sW_g = eye.clone()
        cW_b = eye.expand(n_eta, 2, 2).clone()
        cW_g = eye.clone()

        have_moments = bool(enabled and k_moments is not None)
        if enabled and k_moments is None:
            warnings.warn(
                "theta_whiten=True but k_moments is None (older stats.json); "
                "degeneracy whitening DISABLED.", RuntimeWarning)
        if have_moments:
            km = torch.as_tensor(k_moments, dtype=torch.float64).reshape(-1, 4)
            N = km[:, 0].clamp_min(1.0)
            k1, k2, k4 = km[:, 1] / N, km[:, 2] / N, km[:, 3] / N
            Ntot = km[:, 0].sum().clamp_min(1.0)
            g1 = km[:, 1].sum() / Ntot
            g2 = km[:, 2].sum() / Ntot
            g4 = km[:, 3].sum() / Ntot
            # ρ_scale = −⟨k⟩/√⟨k²⟩  (A,e response basis (1, −k));
            # ρ_smear = +⟨k²⟩/√⟨k⁴⟩ (a,c variance-response basis (1, k²)).
            def _rs(m1, m2):
                return float(-m1 / torch.sqrt(m2.clamp_min(1e-30)))

            def _rc(m2, m4):
                return float(m2 / torch.sqrt(m4.clamp_min(1e-30)))

            for b in range(n_eta):
                sW_b[b] = _whitening_L(_rs(k1[b], k2[b]), max_rho)
                cW_b[b] = _whitening_L(_rc(k2[b], k4[b]), max_rho)
            sW_g = _whitening_L(_rs(g1, g2), max_rho)
            cW_g = _whitening_L(_rc(g2, g4), max_rho)

        self.register_buffer("_scale_W_binned", sW_b, persistent=False)
        self.register_buffer("_scale_W_global", sW_g, persistent=False)
        self.register_buffer("_smear_W_binned", cW_b, persistent=False)
        self.register_buffer("_smear_W_global", cW_g, persistent=False)
        self._whiten_scale_active = bool(have_moments and scale_pair_fit)
        self._whiten_smear_active = bool(have_moments and smear_pair_fit)

    def _theta_cell_idx(self, b_pm, phi_pm):
        """Flat θ-table index per muon: the η-bin ``b_pm`` ('binned'), or
        ``b_pm·n_phi + φ-bin`` ('binned2d', η-outer/φ-inner). φ is binned
        uniformly over [-π, π); out-of-range φ clamps to the edge cells."""
        if not self.theta_grid:
            return b_pm
        pf = (phi_pm + math.pi) * (self.n_phi_bins / (2.0 * math.pi))
        p = torch.clamp(pf.floor().to(torch.long), 0, self.n_phi_bins - 1)
        return b_pm * self.n_phi_bins + p

    def _scale_AeM_pm(self, eta_pm, phi_pm, b_pm) -> torch.Tensor:
        """Per-muon PHYSICAL scale params ``[B, 2, 3] = (A, e, M)``, masked to the
        fitted terms (scale_param_mask, default A,M only — drops the A/e-
        degenerate e). 'binned': the η-bin table θ_scale[b] (O(1) fit param) ×
        THETA_SCALE_REF; 'binned2d': the (η,φ)-cell table; 'mlp': the
        ThetaNet(η, φ), which already applies REF."""
        if self.theta_mode == "mlp":
            aem = self.theta_net(eta_pm, phi_pm)[0]
        else:
            idx = self._theta_cell_idx(b_pm, phi_pm)
            aem = self.theta_scale[idx] * self.theta_scale.new_tensor(THETA_SCALE_REF)
        # Whiten the (A, e) gradient (backward-only; identity in eval/forward).
        if self.training and self._whiten_scale_active:
            L = (self._scale_W_global if self.theta_mode == "mlp"
                 else self._scale_W_binned[b_pm])
            # Whiten in the O(1) θ space — where THETA_SCALE_REF calibrates the
            # curvature diagonal to ≈unit, which is the only space the
            # unit-diagonal L = chol(C⁻¹) is valid in. In PHYSICAL (A,e) space
            # the A/e diagonal ratio is ⟨(∂m/∂A)²⟩/⟨(∂m/∂e)²⟩ = 1/⟨k²⟩ ≈ 140,
            # and a unit-diagonal whitening there grossly mis-scales the step
            # (→ non-finite loss at large lr). Divide out REF, whiten, multiply
            # back: the forward is the identity, only the O(1)-space gradient is
            # rotated. (The smear (a,c) is already whitened in its O(1) space.)
            ref_ae = aem.new_tensor(THETA_SCALE_REF[:2])
            ae = _WhitenGradFn.apply(aem[..., :2] / ref_ae, L) * ref_ae
            aem = torch.cat([ae, aem[..., 2:]], dim=-1)
        return aem * self.scale_param_mask

    def _smear_raw_to_effective(self, raw: torch.Tensor) -> torch.Tensor:
        """Apply the positivity reparameterisation (if any) and the per-bin
        fit mask to the raw O(1) ``θ_smear`` tensor (shape ``[..., 2]``).
        'linear' (default): identity (signed). 'softplus': ``softplus(raw)``
        so each of (a, c) ≥ 0 INDIVIDUALLY. 'square': ``raw²`` — same individual
        positivity, but θ=0 → 0 exactly (identity init) and a non-saturating
        Jacobian (better convergence; see the constructor note on the raw-θ
        covariance caveat at zero). The mask is applied AFTER the transform so
        frozen params (or 'a'-/'c'-only modes) are EXACTLY zero regardless of
        the raw value — keeping the per-muon σ_qop and all downstream
        transformations evaluated to zero for the inactive term."""
        if self.smear_param_form == "softplus":
            raw = F.softplus(raw)
        elif self.smear_param_form == "square":
            raw = raw * raw
        return raw * self.smear_param_mask

    def _smear_ac_pm(self, eta_pm, phi_pm, b_pm) -> torch.Tensor:
        """Per-muon EFFECTIVE qop-resolution variance coefficients ``[B, 2, 2]
        = (a, c)`` (σ²_qop = a + c·k²), with the positivity reparam (if any)
        applied and masked to the fitted term(s). 'binned': ``θ_smear[b]``;
        'mlp': the continuous ``ThetaNet(η, φ)``."""
        if self.theta_mode == "mlp":
            ac = self.theta_net(eta_pm, phi_pm)[1]
        else:
            ac = self.theta_smear[self._theta_cell_idx(b_pm, phi_pm)]
        ac = self._smear_raw_to_effective(ac)
        # Whiten the (a, c) gradient (backward-only; identity in eval/forward).
        # Applied AFTER softplus — the forward pass is the identity so positivity
        # is untouched; only the gradient is rotated.
        if self.training and self._whiten_smear_active:
            L = (self._smear_W_global if self.theta_mode == "mlp"
                 else self._smear_W_binned[b_pm])
            ac = _WhitenGradFn.apply(ac, L)
        return ac

    def effective_theta_smear(self) -> torch.Tensor:
        """Per-η-bin PHYSICAL smear coefficients (a, c), masked (BINNED mode only
        — used for the diagnostics curve / bootstrap σ). The fit parameter
        ``theta_smear`` is O(1) for the optimizer; the physical qop-variance
        coefficients (σ²_qop = a + c·k²) are ``effective(θ_smear) · SMEAR_VAR_SCALE``,
        where ``effective`` applies the positivity reparam (linear ↔ identity,
        softplus ↔ softplus). In 'mlp' mode evaluate ``theta_net`` on an η grid."""
        scale = self.theta_smear.new_tensor([SMEAR_VAR_SCALE_A, SMEAR_VAR_SCALE_C])
        return self._smear_raw_to_effective(self.theta_smear) * scale

    def _scale_per_event(self, eta_pm, phi_pm, b_pm) -> torch.Tensor:
        """Per-event ``[B, 6] = (A_+, e_+, M_+, A_-, e_-, M_-)``."""
        return self._scale_AeM_pm(eta_pm, phi_pm, b_pm).reshape(b_pm.shape[0], -1)

    def _qop_var_pm(self, eta_pm, phi_pm, b_pm, pt_pm) -> torch.Tensor:
        """Per-muon SIGNED qop-resolution variance
        ``σ²_qop,μ = a·SCALE_A + c·SCALE_C·k_μ²`` (k = 1/pt), from the per-muon
        (a, c). ``a``, ``c`` are COMBINED here, before any clipping. ``[B, 2]``."""
        ac = self._smear_ac_pm(eta_pm, phi_pm, b_pm)        # [B,2,2] masked
        a_pm = ac[..., 0] * SMEAR_VAR_SCALE_A
        c_pm = ac[..., 1] * SMEAR_VAR_SCALE_C
        k2 = (1.0 / pt_pm) ** 2
        return a_pm + c_pm * k2                             # [B, 2], signed

    def fold_sigma_qop_pm(
        self, pt_pm: torch.Tensor, eta_pm: torch.Tensor, phi_pm: torch.Tensor,
        b_pm: torch.Tensor,
    ) -> torch.Tensor:
        """Per-muon σ_qop for the validation FOLD, from the fitted (a, c): the
        combined qop variance ``σ²_qop = a + c·k²`` CLIPPED AT 0 *after*
        combining (a stochastic Gaussian qop kick can only broaden, so the
        unsmearing region σ² < 0 → no kick). Returns ``[B, 2]``."""
        return torch.sqrt(self._qop_var_pm(eta_pm, phi_pm, b_pm, pt_pm).clamp_min(0.0))

    # ------------------------------------------------------------------
    # T_scale (analytic + linearized)
    # ------------------------------------------------------------------

    def _delta_qop_analytic(
        self,
        AeM_pm: torch.Tensor,
        pt_pm: torch.Tensor,
        eta_pm: torch.Tensor,
        q_pm: torch.Tensor,
    ) -> torch.Tensor:
        """Analytic δqop per muon (matches ``calculateQopUnc`` in
        ``muon_calibration.hpp``)::

            δqop_i = q_i · sinθ_i · [(A_i − e_i k_i) k_i + q_i M_i]

        ``AeM_pm`` is the per-muon ``[B, 2, 3] = (A, e, M)`` (from
        ``_scale_AeM_pm`` — binned table or ThetaNet).
        """
        sintheta = _sintheta_from_eta(eta_pm)
        k_pm = 1.0 / pt_pm
        A_pm = AeM_pm[..., 0]
        e_pm = AeM_pm[..., 1]
        M_pm = AeM_pm[..., 2]
        k_unc = (A_pm - e_pm * k_pm) * k_pm + q_pm * M_pm
        return q_pm * sintheta * k_unc

    def _qop_new_to_pt(
        self,
        qop: torch.Tensor,
        qop_new: torch.Tensor,
        q_pm: torch.Tensor,
        sintheta: torch.Tensor,
    ) -> torch.Tensor:
        """Invert a shifted qop back to pt as ``pt = |sinθ / qop_new|``.

        pt is a positive magnitude, so we take ``abs`` rather than carrying the
        qop sign: if a large scale/smear shift flips the sign of ``qop_new``,
        that is a PHYSICAL charge mis-reconstruction (real at large measurement
        uncertainty), not something to forbid — the magnitude |1/qop| is the
        right pt either way. The ONLY pathology is ``qop_new == 0`` → pt = ∞, so
        we clamp |qop_new| at ``QOP_EPS`` and nothing else (no resolution-
        suppressing floor on the shift). ``q_pm`` is unused (its sign is absorbed
        into the magnitude); kept in the signature for call-site symmetry.
        """
        return sintheta / qop_new.abs().clamp_min(QOP_EPS)   # sinθ > 0 → pt > 0

    def _apply_scale_pt(
        self,
        pt_pm: torch.Tensor,
        eta_pm: torch.Tensor,
        q_pm: torch.Tensor,
        delta_qop: torch.Tensor,
        sign: float,
    ) -> torch.Tensor:
        """Shift qop by ``sign·δqop`` and convert back to pt.

        Conventions: ``qop = q · sinθ / pt``. The scale transform is DEFINED in
        the BACKWARD (data → MC) direction: ``qop_mc = qop_obs − δqop(pt_obs)``
        with δqop evaluated at the OBSERVED pt — so ``sign=−1`` with the input
        pt is the exact, closed-form defining map (no fixed point anywhere in
        the fit). ``sign=+1`` only realises a single step of the implicit
        forward (MC → data) inverse — see ``_scale_apply_pt_forward``.
        """
        sintheta = _sintheta_from_eta(eta_pm)
        qop = q_pm * sintheta / pt_pm
        qop_new = qop + sign * delta_qop
        return self._qop_new_to_pt(qop, qop_new, q_pm, sintheta)

    def _scale_apply_pt_forward(
        self,
        pt_in: torch.Tensor,
        eta_pm: torch.Tensor,
        phi_pm: torch.Tensor,
        q_pm: torch.Tensor,
        b_pm: torch.Tensor,
        n_iter: int = 8,
    ) -> torch.Tensor:
        """FORWARD (MC → data) scale: the exact functional inverse of the
        defining backward map ``qop_mc = qop_obs − δqop(pt_obs)``, solved by
        fixed point ``qop_obs ← qop_in + δqop(pt_obs)`` (contraction rate
        ~|∂δqop/∂qop| ≲ 2e·k ~ 2e-3; 8 iterations → ~1e-12 relative). The FIT never needs
        this — only the forward consumers do (the diagnostics' physical fold;
        the loader's numpy injection twin), where the iteration cost is free."""
        if not self.scale_enabled:
            return pt_in
        sintheta = _sintheta_from_eta(eta_pm)
        qop_in = q_pm * sintheta / pt_in
        AeM_pm = self._scale_AeM_pm(eta_pm, phi_pm, b_pm)
        pt_obs = pt_in
        for _ in range(n_iter):
            dq = self._delta_qop_analytic(AeM_pm, pt_obs, eta_pm, q_pm)
            pt_obs = self._qop_new_to_pt(qop_in, qop_in + dq, q_pm, sintheta)
        return pt_obs

    def jacobian_mll_linearized(
        self,
        mll: torch.Tensor,
        pt_pm: torch.Tensor,
        q_pm: torch.Tensor,
        b_pm: torch.Tensor,
    ) -> torch.Tensor:
        """Per-event closed-form Jacobian J = ∂m_ll/∂θ_scale_pm.

        From ``m_ll ≈ m_ll · (1 − ½ Σ_i δqop_i / qop_i)`` and the analytic
        ``δqop_i = q_i sinθ_i [(A − e k) k + q_i M]`` (so ``δqop_i/qop_i =
        A − e k_i + q_i M pt_i``, since ``qop_i = q_i sinθ_i k_i``)::

            ∂m_ll/∂A_i = −½ m_ll               (charge-even scale)
            ∂m_ll/∂e_i = +½ m_ll · k_i         (charge-even)
            ∂m_ll/∂M_i = −½ m_ll · q_i · pt_i  (charge-odd / sagitta)

        Returns ``[B, 6]`` ordered ``(A_+, e_+, M_+, A_-, e_-, M_-)``; independent
        of ``b_pm`` (those route the gradient to ``theta_scale`` via scatter).
        """
        k_pm = 1.0 / pt_pm
        m_half = 0.5 * mll
        dA = -m_half.unsqueeze(-1).expand_as(q_pm)
        de = +m_half.unsqueeze(-1) * k_pm
        dM = -m_half.unsqueeze(-1) * q_pm * pt_pm
        return torch.stack(
            [dA[:, 0], de[:, 0], dM[:, 0], dA[:, 1], de[:, 1], dM[:, 1]], dim=-1
        )

    # ------------------------------------------------------------------
    # T_smear
    # ------------------------------------------------------------------

    def apply_smear_pt(
        self,
        pt_pm: torch.Tensor,
        eta_pm: torch.Tensor,
        q_pm: torch.Tensor,
        sigma_qop_pm: torch.Tensor,
        eps_pm: torch.Tensor,
    ) -> torch.Tensor:
        """Gaussian smear of qop_i → return smeared pt_i. ``eps_pm`` is
        the unit-normal noise the caller supplies (lets the trainer
        couple the same draw between an MLE branch and a flow-aux branch
        if it wants).
        """
        sintheta = _sintheta_from_eta(eta_pm)
        qop = q_pm * sintheta / pt_pm
        delta = sigma_qop_pm * eps_pm
        qop_new = qop + delta
        return self._qop_new_to_pt(qop, qop_new, q_pm, sintheta)

    # ------------------------------------------------------------------
    # Conditioning vector for the flow
    # ------------------------------------------------------------------

    def _standardise_mll(self, mll: torch.Tensor) -> torch.Tensor:
        return (mll - self.mll_mean_buf) / self.mll_std_buf

    # ------------------------------------------------------------------
    # Nominal (θ=0) flow density — the stage-1 template
    # ------------------------------------------------------------------

    def log_p_nominal(
        self, mll_obs: torch.Tensor, muon_kin_std_obs: torch.Tensor
    ) -> torch.Tensor:
        """``log p₀(m | muon_kin)`` — the θ=0 nominal flow density (1/GeV).

        The stage-1 target: a plain conditional density on the (uncorrected)
        reco mass, conditioned on the leak-free kinematics only — the flow
        never sees θ (the two-stage continuity design).
        """
        mll_std = self._standardise_mll(mll_obs).clamp(
            -MLL_STD_FLOW_CLAMP, MLL_STD_FLOW_CLAMP
        )
        return self.flow(mll_std.unsqueeze(-1), muon_kin_std_obs) - self.mll_log_scale

    # ------------------------------------------------------------------
    # MLP coefficients (data branch only)
    # ------------------------------------------------------------------

    def f_data(self, muon_kin_std: torch.Tensor) -> torch.Tensor:
        if self.bkg_global:
            # Constant fractions for every event: softmax(global_logits) → [n+1]
            # (signal LAST), broadcast to [B, n+1].
            f = torch.softmax(self.bkg_global_logits, dim=-1)
            return f.to(muon_kin_std.dtype).unsqueeze(0).expand(
                muon_kin_std.shape[0], -1)
        return self.mlp(muon_kin_std)

    # ------------------------------------------------------------------
    # Conditioning basis
    # ------------------------------------------------------------------

    def _cond_from_muons(self, pt_pm, eta_pm, phi_pm, q_pm) -> torch.Tensor:
        """STANDARDISED flow/MLP conditioning recomputed from per-muon momenta,
        for ``self.cond_basis``. Used wherever the qop scale/smear changes pt and
        the conditioning must follow (operator un-kick nodes, pseudo-data).

        - ``muon_kin``: ``(η±, cosφ±, sinφ±, ρ)`` — only ρ depends on pt; η/φ are
          pt-invariant, so recomputing the full vector reproduces the ρ-only
          update (kept as a fast path in ``_node_cond``).
        - ``event_level``: the dilepton vars via ``_event_cond_raw`` (all
          dilation-invariant, residually pt-dependent through the
          charge-differential terms). ``q_pm`` is unused (μ± identified by
          index)."""
        if self.cond_basis == "event_level":
            raw = _event_cond_raw(pt_pm, eta_pm, phi_pm)
            return (raw - self.y_event_mean) / self.y_event_std
        cphi = torch.cos(phi_pm)
        sphi = torch.sin(phi_pm)
        rho = ((pt_pm[..., 0] - pt_pm[..., 1])
               / (pt_pm[..., 0] + pt_pm[..., 1]))
        raw = torch.stack(
            [eta_pm[..., 0], eta_pm[..., 1], cphi[..., 0], sphi[..., 0],
             cphi[..., 1], sphi[..., 1], rho], dim=-1)
        return (raw - self.muon_kin_mean) / self.muon_kin_std

    def _node_cond(self, mk_g, pt_truth, eta_pm, phi_pm, q_pm) -> torch.Tensor:
        """Per-node conditioning from un-kicked ``pt_truth``. For ``muon_kin``,
        patch only ρ into the (already-expanded, observed) ``mk_g`` — η/φ are
        pt-invariant, so this is bit-for-bit the legacy behaviour. For
        ``event_level``, recompute the whole vector from the un-kicked momenta."""
        if self.cond_basis == "event_level":
            return self._cond_from_muons(pt_truth, eta_pm, phi_pm, q_pm)
        rho = ((pt_truth[..., 0] - pt_truth[..., 1])
               / (pt_truth[..., 0] + pt_truth[..., 1]))
        idx = N_MUON_KIN - 1
        mk_g[..., idx] = (rho - self.muon_kin_mean[idx]) / self.muon_kin_std[idx]
        return mk_g

    # ------------------------------------------------------------------
    # Stage-2 continuity-equation data fit (frozen flow + analytic v, κ)
    #
    #   log p_s(m|c,θ) = log p₀(m|c) + δ(m,c;θ)        [+ O(θ²) renorm]
    #   δ = Σ_k θ_k g_k,   g_scale = −v′ − v·s,   g_smear = ½ κ (s′ + s²)
    #
    # p₀ is the frozen θ=0 flow (trained in stage 1, conditioned on muon_kin
    # only). v = ∂m/∂θ_scale (analytic linearised Jacobian), v′ = ∂_m v (through
    # the (m,c)→pt reconstruction), κ = ∂Var[m]/∂θ_smear (analytic). The data
    # θ-gradient flows by autograd through this δ; the flow gets no gradient.
    # ------------------------------------------------------------------

    def _flow_logp_score(self, mll: torch.Tensor, muon_kin_std: torch.Tensor):
        """Nominal density + score from the frozen θ=0 flow.

        Returns ``(log p₀, s, s′)`` (each ``[B]``): p₀ in 1/GeV,
        ``s = ∂_m log p₀``, ``s′ = ∂²_m log p₀``, via autograd in ``m``.
        """
        m = mll.detach().requires_grad_(True)
        mll_std = self._standardise_mll(m).clamp(-MLL_STD_FLOW_CLAMP, MLL_STD_FLOW_CLAMP)
        cond = self._build_flow_cond(muon_kin_std, None, None)
        logp = self.flow(mll_std.unsqueeze(-1), cond) - self.mll_log_scale  # [B]
        s = torch.autograd.grad(logp.sum(), m, create_graph=True)[0]
        s_prime = torch.autograd.grad(s.sum(), m, create_graph=True)[0]
        return logp, s, s_prime

    @staticmethod
    def _reconstruct_pt(mll, eta_pm, phi_pm, rho):
        """``(m, η_±, φ_±, ρ) → (pt₊, pt₋)`` for the (massless) dimuon.

        ``ρ = (pt₊−pt₋)/(pt₊+pt₋)``; with ``S = pt₊+pt₋`` and
        ``m² = ½ S² (1−ρ²)(cosh Δη − cos Δφ)`` this pins the pt scale. Used to
        propagate ``v`` along the `m`-direction at fixed conditioning.
        """
        d_eta = eta_pm[:, 0] - eta_pm[:, 1]
        d_phi = phi_pm[:, 0] - phi_pm[:, 1]
        ang = (torch.cosh(d_eta) - torch.cos(d_phi)).clamp_min(1e-6)
        S2 = 2.0 * mll * mll / ((1.0 - rho * rho).clamp_min(1e-6) * ang)
        S = torch.sqrt(S2.clamp_min(1e-12))
        return torch.stack([S * (1.0 + rho) * 0.5, S * (1.0 - rho) * 0.5], dim=-1)

    def _v_and_vprime(self, mll, eta_pm, phi_pm, q_pm, b_pm, rho):
        """Advective velocity ``v = ∂m/∂θ_scale_pm`` ``[B,6]`` and its `m`-
        derivative ``v′ = ∂_m v`` ``[B,6]`` at fixed conditioning.

        ``v`` is the analytic linearised mass-Jacobian; ``v′`` is its *total*
        derivative in ``m`` (the continuity "missing-Jacobian" term ``∇·v``),
        taken through the ``(m,c)→pt`` reconstruction so the kinematics track
        ``m`` at fixed ``c``.
        """
        m = mll.detach().requires_grad_(True)
        pt = self._reconstruct_pt(m, eta_pm, phi_pm, rho)
        v = self.jacobian_mll_linearized(m, pt, q_pm, b_pm)  # [B,6]
        vprime = torch.stack(
            [torch.autograd.grad(v[:, k].sum(), m, create_graph=True, retain_graph=True)[0]
             for k in range(v.shape[-1])],
            dim=-1,
        )
        return v, vprime

    def _kappa_smear(self, mll, pt_pm, eta_pm, phi_pm, q_pm):
        """Diffusion coefficients ``κ = ∂Var[m]/∂θ_smear`` per per-muon term,
        ``[B,4] = (a₊,c₊,a₋,c₋)``.

        With the variance basis ``σ_qop² = a²·1 + c²·k²`` (so ``∂σ²/∂a²=1``,
        ``∂σ²/∂c²=k²``): ``κ_{a,i}=(∂m/∂qop_i)²`` and ``κ_{c,i}=(∂m/∂qop_i)² k_i²``.
        ``∂m/∂qop_i`` is taken by autograd through ``pt_i = q_i sinθ_i / qop_i``.
        """
        sintheta = _sintheta_from_eta(eta_pm)
        qop = (q_pm * sintheta / pt_pm).detach().requires_grad_(True)
        pt_from_qop = q_pm * sintheta / qop
        m = _event_mll(pt_from_qop, eta_pm, phi_pm)
        dm_dqop = torch.autograd.grad(m.sum(), qop, create_graph=True)[0]  # [B,2]
        d2 = dm_dqop * dm_dqop  # [B,2]
        k2 = (1.0 / pt_pm) ** 2  # [B,2]
        # order (a₊, c₊, a₋, c₋)
        return torch.stack([d2[:, 0], d2[:, 0] * k2[:, 0],
                            d2[:, 1], d2[:, 1] * k2[:, 1]], dim=-1)

    def _smear_per_event_linear(self, b_pm: torch.Tensor) -> torch.Tensor:
        """Per-event smear *increments* ``[B,4] = (a₊,c₊,a₋,c₋)`` for the
        continuity tilt — the EFFECTIVE (a, c) after the positivity reparam
        (if any) and masked to the fitted term(s)."""
        th = self._smear_raw_to_effective(self.theta_smear)  # [n_eta, 2]
        return th[b_pm].reshape(b_pm.shape[0], -1)

    def _continuity_g(self, m, mk, eta_pm, phi_pm, q_pm, b_pm, rho, pt_pm):
        """First-order continuity sensitivities at ``(m, c)``:
        ``g_scale = −v′ − v·s`` ``[P,6]`` and ``g_smear = ½ κ (s′+s²)`` ``[P,4]``.
        ``pt_pm`` is the per-point pt (observed at the data mass; reconstructed
        from ``(m,c)`` on the normalisation grid)."""
        _, s, s_prime = self._flow_logp_score(m, mk)
        if self.scale_enabled:
            v, vprime = self._v_and_vprime(m, eta_pm, phi_pm, q_pm, b_pm, rho)
            g_scale = -vprime - v * s.unsqueeze(-1)
        else:
            g_scale = m.new_zeros((m.shape[0], N_THETA_SCALE_PM))
        if self.smearing_enabled:
            kappa = self._kappa_smear(m, pt_pm, eta_pm, phi_pm, q_pm)
            g_smear = 0.5 * kappa * (s_prime + s * s).unsqueeze(-1)
        else:
            g_smear = m.new_zeros((m.shape[0], N_THETA_SMEAR_PM))
        return g_scale, g_smear

    def _continuity_logZ(self, mk, eta_pm, phi_pm, q_pm, b_pm, rho,
                         theta_pm, n_grid: int = 32):
        """2nd-order cumulant log-normalisation per event:
        ``logZ ≈ E_{p₀}[δ] + ½ Var_{p₀}[δ|c]`` with ``δ = θ_pm · g(m,c)``.

        This is ``log E_{p₀}[e^δ]`` to ``O(δ³)``; the mean term corrects for the
        grid truncation/discretisation (analytically ``E_{p₀}[g]=0`` over full
        support, but not on a finite grid). ``g`` is detached (θ-independent,
        flow-frozen), so only ``θ_pm`` carries the gradient, giving the proper
        score centering ``∂_θ logZ = E[g] + Cov[g]·θ``. Moments are estimated by
        ``p₀``-weighted quadrature on an ``n_grid`` mass grid.
        ``theta_pm = [θ_scale_pm (6), θ_smear_pm (4)]`` ``[B,10]``.
        """
        B = mk.shape[0]
        dev, dt = mk.device, mk.dtype
        mg = torch.linspace(self._m_lo_f + 1e-3, self._m_hi_f - 1e-3,
                            n_grid, device=dev, dtype=dt)            # [G]
        # expand per-event conditioning across the grid, flatten to [B*G, ...]
        def rep(x):
            return x.unsqueeze(1).expand(B, n_grid, *x.shape[1:]).reshape(
                B * n_grid, *x.shape[1:])
        m_f = mg.unsqueeze(0).expand(B, n_grid).reshape(-1)          # [B*G]
        eta_f, phi_f, q_f, mk_f = rep(eta_pm), rep(phi_pm), rep(q_pm), rep(mk)
        b_f, rho_f = rep(b_pm), rep(rho)
        pt_f = self._reconstruct_pt(m_f, eta_f, phi_f, rho_f)
        with torch.no_grad():
            logp0_f = self.log_p_nominal(m_f, mk_f)                  # [B*G]
        gs, gsm = self._continuity_g(m_f, mk_f, eta_f, phi_f, q_f, b_f, rho_f, pt_f)
        g = torch.cat([gs, gsm], dim=-1).detach().reshape(B, n_grid, -1)  # [B,G,10]
        w = torch.softmax(logp0_f.reshape(B, n_grid), dim=-1)        # p₀ weights, Σ=1
        d = (g * theta_pm.unsqueeze(1)).sum(-1)                      # δ on grid [B,G]
        mean = (w * d).sum(-1)
        var = (w * d * d).sum(-1) - mean * mean
        return mean + 0.5 * var.clamp_min(0.0)                       # 2nd-order cumulant

    # ------------------------------------------------------------------
    # Stage-2 continuity data fit — #2 "forward-fold the flow's eval point"
    #
    #   p_θ(x) = E_{ε~N(0,1)}[ p₀(m'(ε)|c) / |G'(m'(ε))| ],
    #   x = m' + s_adv(m') + √(V(m'))·ε ,   G'(m') = ∂x/∂m'
    #     = 1 + s_adv'(m') + (V'/2√V)·ε                       (kernel Jacobian)
    #
    # s_adv(m') = Σ_k v_k(m')·θ_scale_k  (advective mass shift; v = analytic J),
    # V(m')     = Σ_k κ_k(m')·softplus(θ_smear)_k² ≥ 0  (smear variance; the
    #             effective softplus(θ) are qop STDs, so σ_qop² = a²+c²k²).
    # This evaluates the FROZEN flow only as point values at the source m'(ε)
    # (no flow derivatives), captures advection+smear to all orders in θ and the
    # x-variation of v, V (source-evaluation + Jacobian), and is normalised by
    # construction. v, κ are evaluated at the source via the pt∝m scaling at
    # fixed conditioning — swappable to a learned v/κ MLP without touching this.
    # softplus/square on θ_smear keeps V ≥ 0 (no ill-posed de-convolution /
    # sharpening); 'linear' is two-sided.
    # ------------------------------------------------------------------

    def _continuity_response(self, m_eval, m_obs, pt_obs, eta_pm, q_pm,
                             theta_scale_pm):
        """Advective mass shift ``s_adv`` at evaluation mass ``m_eval``
        (broadcasting against the per-event observables), for the analytic scale
        transform with ``pt(m_eval) = pt_obs · m_eval/m_obs``.

        ``v = (−½m, ½m k, −½m q pt)`` per muon (the corrected scale Jacobian).
        Returns ``s_adv`` shaped like ``m_eval``. (Replaceable by a learned MLP.)
        """
        scale = (m_eval / m_obs).unsqueeze(-1)             # [...,1]
        pt = pt_obs * scale                                # [...,2]
        k = 1.0 / pt
        mh = (0.5 * m_eval).unsqueeze(-1)                  # [...,1]
        dA = -mh                                           # charge-even (per muon)
        de = mh * k                                        # [...,2]
        dM = -mh * q_pm * pt                               # [...,2]
        v = torch.stack([dA[..., 0], de[..., 0], dM[..., 0],
                         dA[..., 0], de[..., 1], dM[..., 1]], dim=-1)
        return (v * theta_scale_pm).sum(-1)

    def _smear_mass_var(self, eta_pm, phi_pm, b_pm, pt_obs, m_eval, m_obs):
        """Per-event SIGNED m_ll variance added by the per-muon qop smear,
        evaluated at the mass ``m_eval`` (NOT a fixed reference): with
        ``pt(m_eval) = pt_obs·m_eval/m_obs`` (pt ∝ m at fixed angles, as in
        ``_continuity_response``) and ``σ²_qop,μ = a + c·k²``,
        ``V = (m_eval/2)² Σ_μ σ²_qop,μ / qop_μ²`` (``∂m/∂qop = −m/2qop``,
        ``1/qop² = pt²/sin²θ``) ``∝ a·m⁴ + c·m²``. Two-sided (V<0 = unsmear).
        Evaluated at the source m' inside ``_continuity_logp`` so the smear's
        mass-dependence enters the change-of-variables Jacobian (consistent with
        the advection). This is the diffusion 'time' of the probability flow."""
        pt = pt_obs * (m_eval / m_obs).unsqueeze(-1)       # pt(m_eval) [B,2]
        vq = self._qop_var_pm(eta_pm, phi_pm, b_pm, pt)    # [B,2] signed σ²_qop(pt)
        sinth = _sintheta_from_eta(eta_pm)
        inv_qop2 = (pt * pt) / (sinth * sinth)             # 1/qop² = pt²/sin²θ
        return (0.5 * m_eval) ** 2 * (vq * inv_qop2).sum(-1)  # [B] signed

    def _flow_score(self, m, mk):
        """Flow score ``∂_m log p₀(m | mk)``, via autograd. Differentiable w.r.t.
        ``m`` when ``m`` carries grad — so the change-of-variables Jacobian of the
        probability-flow smear picks up the local curvature ``∂²_m log p₀``.
        Otherwise returns a detached value (the source fixed-point path)."""
        if m.requires_grad:
            lp = self.log_p_nominal(m, mk)
            return torch.autograd.grad(lp.sum(), m, create_graph=True)[0]
        with torch.enable_grad():
            ml = m.detach().requires_grad_(True)
            lp = self.log_p_nominal(ml, mk)
            return torch.autograd.grad(lp.sum(), ml, create_graph=True)[0].detach()

    def _scale_unapply_pt(self, pt_obs, eta_pm, phi_pm, q_pm, b_pm):
        """THE defining backward (data → MC) scale map on a per-muon config:
        ``qop_mc = qop_obs − δqop(pt_obs)`` with δqop evaluated at the OBSERVED
        pt — exact and iteration-free by construction. Identity if scale is
        disabled. Used by the operators (un-kick + source conditioning) and the
        boundary transport of the window normalisation."""
        if not self.scale_enabled:
            return pt_obs
        AeM_pm = self._scale_AeM_pm(eta_pm, phi_pm, b_pm)
        dqop_s = self._delta_qop_analytic(AeM_pm, pt_obs, eta_pm, q_pm)
        return self._apply_scale_pt(pt_obs, eta_pm, q_pm, dqop_s, sign=-1.0)

    def _scale_source_rho_std(self, pt_obs, eta_pm, phi_pm, q_pm, b_pm):
        """Standardised SOURCE ρ from un-applying the scale only (no smear) —
        the smear is a pure mass-space transport, so it leaves ρ untouched.
        Returns ``[B]``."""
        pt1 = self._scale_unapply_pt(pt_obs, eta_pm, phi_pm, q_pm, b_pm)
        rho = (pt1[:, 0] - pt1[:, 1]) / (pt1[:, 0] + pt1[:, 1])
        idx = N_MUON_KIN - 1
        return (rho - self.muon_kin_mean[idx]) / self.muon_kin_std[idx]

    def _scale_backward_mass_linear(self, m_eval, m_obs, pt_obs, eta_pm, q_pm,
                                    theta_scale_pm):
        """Explicit BACKWARD (data → MC) scale shift in mass space for the
        mass-space operators (pf_ode / gh_convolution), in their linearised
        representation: ``m_s(x) = x − s_adv(x)`` with the advective shift
        ``s_adv = v(x)·θ`` evaluated at the EVALUATION mass ``m_eval`` along the
        event's observed ray (``pt(x) = pt_obs·x/m_obs`` inside
        ``_continuity_response``) — exact and iteration-free in the defining
        direction. The exact ``log|dm_s/dx| = log|1 − ∂_x s_adv|`` is CLOSED
        FORM along the event ray (pt(x) = pt_obs·x/m_obs): the A term of
        s_adv is linear in x (slope −½ΣA_μ), the e term is CONSTANT
        ((x/2)·k(x) = m_obs/(2·pt_obs,μ)), and the M term is quadratic
        (−x²·q_μ·pt_obs,μ·M_μ/(2·m_obs)), so

            ∂s_adv/∂x = −½(A₊+A₋) − (x/m_obs)·Σ_μ q_μ·pt_obs,μ·M_μ.

        θ-gradients flow through the closed form by ordinary autograd.
        Returns ``(m_s, log_J)`` shaped like ``m_eval``."""
        if not self.scale_enabled:
            return m_eval, torch.zeros_like(m_eval)
        s_adv = self._continuity_response(
            m_eval, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)
        sp = (-0.5 * (theta_scale_pm[..., 0] + theta_scale_pm[..., 3])
              - (m_eval / m_obs)
              * (q_pm[..., 0] * pt_obs[..., 0] * theta_scale_pm[..., 2]
                 + q_pm[..., 1] * pt_obs[..., 1] * theta_scale_pm[..., 5]))
        m_s = m_eval - s_adv
        log_J = torch.log((1.0 - sp).abs().clamp_min(1e-6))
        return m_s, log_J

    def _continuity_logp(self, m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm,
                         n_iter: int = 2):
        """``log p_θ(x|c)`` — dispatches on ``self.smear_operator``.

        ``"gh_convolution"``: EXACT stochastic Gaussian convolution via
        Gauss-Hermite quadrature (matches the per-muon qop fold operator that
        generates the pseudo-data; converges to the convolution by construction
        at any V). See ``_continuity_logp_gh``.

        ``"pf_ode"`` (default, doc below): deterministic probability-flow
        ODE — cheaper but over-broadens at large V/σ²."""
        if self.smear_operator == "gh_convolution":
            return self._continuity_logp_gh(
                m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm,
                n_gh=self.n_gh_nodes, n_iter=n_iter)
        if self.smear_operator == "gh_convolution_qop":
            return self._continuity_logp_gh_qop(
                m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm,
                n_gh=self.n_gh_nodes, n_iter=n_iter)
        return self._continuity_logp_pf_ode(
            m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm, n_iter=n_iter)

    def _continuity_logp_pf_ode(self, m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm,
                                b_pm, n_iter: int = 2):
        """``log p_θ(x|c)``: the frozen nominal flow pushed through an exactly-
        normalized, INVERTIBLE TRANSPORT. The scale is the explicit BACKWARD
        (data → MC) shift ``m_s = x − s_adv(x)`` applied first at the OBSERVED
        mass (defining direction — closed form, exact Jacobian, no fixed
        point); the smear keeps its forward (MC → data) convention as a
        score-driven PROBABILITY-FLOW displacement
        ``y ← y − (V/2n)·∂_m log p₀(y)`` (``smear_flow_steps`` Euler steps; the
        deterministic equivalent of a Gaussian qop smear of mass-variance V,
        broaden V>0 / sharpen V<0), inverted from the scaled coordinate ``m_s``
        for the source m' by fixed point (smear shift only).

        Log-Jacobian — two forms (``self.jacobian_form``):

        * ``"softlog"`` (default): use the EXACT autograd Jacobian ``G' =
          dx/dm'`` of the N-step Euler forward map, fed through
          ``_softlog_below_floor`` — equal to ``log(G'.clamp_min(floor))`` in
          the physical region AND adds a C¹ quadratic barrier ``(floor−G')⁺²
          /(2·floor²)`` past the floor. The barrier value and slope are zero
          at the seam (so the physical region is unchanged) and grow
          quadratically past it, pulling G' back toward the physical region
          — replacing the hard clamp's flat-NLL basin (where the optimiser
          could drift into G' < 0 when V·∂²_m log p₀ is large at large |η|,
          with no gradient cost) with a gentle restoring force.
        * ``"exp"``: skip the autograd Jacobian and use the FROZEN-SCORE
          continuous-flow approximation ``log G' = log(1+s_adv'(m')) − V·∂²_m
          log p₀(m')/2`` — always finite, no floor needed. Approximates a
          DIFFERENT operator: assumes ∂²_m log p₀ is constant along the smear
          trajectory and is the analytic ``N→∞`` Euler-step limit of that
          frozen-score flow. Closer to the true Gaussian convolution at
          moderate V than 1-step Euler, but the assumption breaks near sharp
          features (the J/ψ peak crest) and the recovered θ may shift vs the
          'softlog' (different operator → different optimum). It also lacks
          the floor's natural cap on unphysical-sharpening rewards, so a
          large negative V is rewarded UNBOUNDEDLY by ``−log G' = +V·∂²/2``
          — only choose this when the V-too-large breakdown driving #1 is
          actually the dominant issue.

        The smear is a pure mass-space transport, so it leaves ρ untouched —
        only the scale's ρ shift is propagated to the conditioning."""
        B = m_obs.shape[0]
        theta_scale_pm = (self._scale_per_event(eta_pm, phi_pm, b_pm)
                          if self.scale_enabled
                          else m_obs.new_zeros((B, N_THETA_SCALE_PM)))
        n_step = max(1, int(self.smear_flow_steps))

        mk_src = mk
        if self.scale_enabled:
            if self.cond_basis == "event_level":
                pt1 = self._scale_unapply_pt(pt_obs, eta_pm, phi_pm, q_pm, b_pm)
                mk_src = self._cond_from_muons(pt1, eta_pm, phi_pm, q_pm)
            else:
                mk_src = mk.clone()
                mk_src[..., N_MUON_KIN - 1] = self._scale_source_rho_std(
                    pt_obs, eta_pm, phi_pm, q_pm, b_pm)

        # BACKWARD (data → MC) scale, explicit at the OBSERVED mass — the
        # defining direction: m_s = m_obs − s_adv(m_obs), with the exact scale
        # Jacobian log|dm_s/dm_obs|. The flow + smear inversion below act on
        # the scaled coordinate m_s (no scale term in the fixed point).
        m_s, log_J_s = self._scale_backward_mass_linear(
            m_obs, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)

        def forward_sm(mp):
            # smear-only forward map (MC → data convention unchanged): the
            # probability-flow score displacement. V (diffusion time) is
            # evaluated at the SOURCE mass mp — so its mass-dependence
            # (V ∝ a·m⁴ + c·m²) enters the autograd Jacobian G'.
            y = mp
            if self.smearing_enabled:
                V = self._smear_mass_var(eta_pm, phi_pm, b_pm, pt_obs, mp, m_obs)
                dt = V / (2.0 * n_step)
                for _ in range(n_step):
                    y = y - dt * self._flow_score(y, mk_src)
            return y

        if not self.smearing_enabled:
            # Scale-only: the source is the explicit backward image — no fixed
            # point, no smear Jacobian (one of the backward-direction
            # simplifications).
            log_p_theta = self.log_p_nominal(m_s, mk_src) + log_J_s
            return torch.nan_to_num(log_p_theta.clamp(max=50.0), nan=0.0)

        # invert the SMEAR-ONLY map for the source m': fixed point
        # m' = m_s − (forward_sm(m') − m'). (Graph is built in training so m'
        # carries the θ dependence, as for G'.)
        mp = m_s.clone()
        for _ in range(n_iter):
            mp = m_s - (forward_sm(mp) - mp)
        # change-of-variables log-Jacobian of the smear-only map: dispatch on
        # jacobian_form (the scale Jacobian is the explicit log_J_s above).
        if self.jacobian_form == "exp":
            # Frozen-score continuous-flow approximation, no floor needed.
            log_Gp = self._log_jacobian_exp(mp, mk_src, None, eta_pm, phi_pm,
                                            b_pm, pt_obs, m_obs)
        else:  # "softlog" (default): autograd Jacobian + tangent extension below floor.
            if mp.requires_grad:
                Gp = torch.autograd.grad(forward_sm(mp).sum(), mp, create_graph=True)[0]
            else:
                with torch.enable_grad():
                    mpj = mp.detach().requires_grad_(True)
                    Gp = torch.autograd.grad(forward_sm(mpj).sum(), mpj)[0].detach()
            log_Gp = _softlog_below_floor(Gp, SMEAR_GP_FLOOR)
        logp0 = self.log_p_nominal(mp, mk_src)
        log_p_theta = logp0 - log_Gp + log_J_s
        # Cap log_p_theta from ABOVE only (the overflow direction). The
        # downstream mixture in `data_nll_continuity` takes `.exp()` of this:
        # for jacobian_form='exp' the unbounded sharpening reward can drive
        # log_Gp → −∞ on a single event → log_p_theta → +∞ → `.exp()` → inf
        # → NaN loss. The physical log-density on this m-window lives in
        # ~[−30, +5], so cap = +50 leaves the operating regime untouched.
        # No LOWER cap: the softlog barrier deliberately drives log_p_theta
        # very negative to penalise unphysical G' < 0, and `.exp()` of a
        # large-negative number underflows cleanly to 0 (the mixture then
        # collapses to the Bernstein background — no NaN). nan_to_num catches
        # any residual NaN (rare, from autograd at boundary mp values).
        return torch.nan_to_num(log_p_theta.clamp(max=50.0), nan=0.0)

    def _continuity_logp_gh(self, m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm,
                             b_pm, n_gh: int = 8, n_iter: int = 2):
        """``log p_θ(x|c)`` via the EXACT Gaussian convolution operator,
        implemented as Gauss-Hermite quadrature of the per-event source map.
        The scale enters as the explicit BACKWARD (data → MC) shift
        ``m_s = x − s_adv(x)`` applied first at the OBSERVED mass (defining
        direction — closed form + exact Jacobian ``|dm_s/dx|``); the smear
        keeps its forward (MC → data) convention, inverted from the scaled
        coordinate:

            ``p_θ(x|c) = |dm_s/dx| · E_{ε~N(0,1)}[ p_0(m'(ε)|c_src) / |G'(m'(ε))| ]``
            ``       ≈ |dm_s/dx| · Σ_i W_i · p_0(m'_i | c_src,i) / |G'_i|``

        where, for each GH node ε = ξ_i:
            ``m_s = m'_i + √V(m'_i)·ξ_i``  (smear-only forward map)
            ``G'(m') = 1 + (V'(m')/(2√V(m')))·ε``  (autograd Jacobian)

        Equivalent — by construction — to the per-muon qop fold that generates
        the validation pseudo-data, so the closure-target curve (flow at
        injected θ) overlaps the pseudo-data at ALL V (no operator-level
        residual; only GH quadrature truncation, which is exponential in n_gh
        for smooth p_0). Cost: ``n_gh`` per-node `p_0` evaluations + the
        autograd through the source map.

        Requires V ≥ 0 (Gaussian variance is non-negative); we ``clamp_min(0)``
        the per-muon σ²_qop before the √V term. With smear_param_form='softplus'
        this is automatic; with 'linear' the clamp enforces it (the fit loses
        access to V < 0 sharpening, which is fine — that region is unphysical
        for a stochastic smear anyway). The ``jacobian_form`` switch is
        ignored here (the Jacobian is the source-map's exact autograd-derived
        ``G'`` — no PF-ODE-style score-flow Jacobian to soft-floor).

        The source ρ for the flow conditioning per GH node carries the smear's
        conditional-mean per-muon qop shift (``E[δqop_μ | δm] = (J^m_μ σ²_qop,μ
        /√V)·ε``, un-applied per node) — see ``_source_rho_std_gh``.
        """
        B = m_obs.shape[0]
        theta_scale_pm = (self._scale_per_event(eta_pm, phi_pm, b_pm)
                          if self.scale_enabled
                          else m_obs.new_zeros((B, N_THETA_SCALE_PM)))
        xi, logW = _gh_nodes(n_gh, m_obs.device, m_obs.dtype)            # [G], [G]
        G = xi.shape[0]

        # Broadcast all event conditioning with a singleton "GH node" dim.
        mo  = m_obs.unsqueeze(1)                                          # [B, 1]
        pto = pt_obs.unsqueeze(1)                                         # [B, 1, 2]
        etao = eta_pm.unsqueeze(1)                                        # [B, 1, 2]
        phio = phi_pm.unsqueeze(1)                                        # [B, 1, 2]
        qo  = q_pm.unsqueeze(1)                                           # [B, 1, 2]
        bpo = b_pm.unsqueeze(1).expand(B, G, b_pm.shape[-1])              # [B, G, 2]
        xig = xi.view(1, G)                                               # [1, G]

        # BACKWARD (data → MC) scale, explicit at the OBSERVED mass — the
        # defining direction (closed form + exact Jacobian); the GH smear
        # quadrature below acts on the scaled coordinate m_s.
        m_s, log_J_s = self._scale_backward_mass_linear(
            m_obs, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)
        ms_ = m_s.unsqueeze(1)                                            # [B, 1]

        # Smear mass-variance in POLYNOMIAL closed form — algebraically equal
        # to _smear_mass_var on the event ray (pt(me) = pt_obs·me/m_obs):
        #   V(me) = α·me⁴ + γ·me²,
        #   α = ¼ Σ_μ a_μ·pt_obs,μ² / (m_obs²·sin²θ_μ),  γ = ¼ Σ_μ c_μ/sin²θ_μ,
        # so V′(me) = 4α·me³ + 2γ·me is available WITHOUT autograd. θ-gradients
        # flow through (a, c) by ordinary autograd (incl. the whitening hook).
        if self.smearing_enabled:
            ac = self._smear_ac_pm(eta_pm, phi_pm, b_pm)            # [B, 2, 2]
            sin2 = _sintheta_from_eta(eta_pm) ** 2                  # [B, 2]
            al_ = (0.25 * (ac[..., 0] * SMEAR_VAR_SCALE_A
                           * pt_obs * pt_obs / sin2).sum(-1)
                   / (m_obs * m_obs)).unsqueeze(1)                  # [B, 1]
            ga_ = (0.25 * (ac[..., 1] * SMEAR_VAR_SCALE_C
                           / sin2).sum(-1)).unsqueeze(1)            # [B, 1]

        def _smear_disp(me):
            """The √V·ε displacement per GH node at evaluation mass me [B, G];
            vanishes exactly if disabled.

            The +EPS inside the sqrt keeps the V→0 behaviour finite (else
            V′/(2√V)→∞ in G′ and the θ-gradient hits 0·∞ = NaN, which
            nan_to_num turns into log p=0 → a spurious UNIFORM density at c≈0,
            which acts as a barrier blocking the fitted smear from reaching 0.
            EPS=1e-12 → a mass-displacement floor of 1e-6 GeV, negligible vs the
            physical √V ~ tens of MeV)."""
            if self.smearing_enabled:
                V = (al_ * me**4 + ga_ * me**2).clamp_min(0.0)
                return (V + 1e-12).sqrt() * xig
            return me_zero

        # Pre-allocate the zero displacement (used only when smearing_enabled=False).
        me_zero = m_obs.new_zeros((B, G))

        # Fixed-point source solve (smear shift ONLY — the scale is the explicit
        # backward image m_s above):  m'_i = m_s − √V(m'_i)·ξ_i. With smearing
        # disabled the source IS m_s — no fixed point at all (a backward-
        # direction simplification).
        mp = ms_.expand(B, G).clone()
        if self.smearing_enabled:
            for _ in range(n_iter):
                mp = ms_ - _smear_disp(mp)

        # Smear-map Jacobian G' = ∂m_s/∂m' in CLOSED FORM at the converged
        # source (≡ 1 with smearing disabled; the scale Jacobian is log_J_s):
        #   G′ = 1 + (V′(m′)/(2√(V(m′)+EPS)))·ξ,  V′ masked where the V ≥ 0
        # clamp is active (matching the old autograd through clamp_min).
        if not self.smearing_enabled:
            Gp = torch.ones_like(mp)
        else:
            V_raw = al_ * mp**4 + ga_ * mp**2
            V = V_raw.clamp_min(0.0)
            Vp = ((4.0 * al_ * mp**3 + 2.0 * ga_ * mp)
                  * (V_raw >= 0).to(mp.dtype))
            Gp = 1.0 + Vp / (2.0 * (V + 1e-12).sqrt()) * xig

        # Frozen-flow density at the per-node source points. ρ in the
        # conditioning is propagated per-node (scale un-applied + smear's
        # conditional-mean δqop un-applied per ε); η/φ are transform-invariant.
        mk_g = mk.unsqueeze(1).expand(B, G, mk.shape[-1]).clone()
        if self.scale_enabled or self.smearing_enabled:
            if self.cond_basis == "event_level":
                # smear is guarded to the qop operator, so this is scale-only:
                # the full event-level conditioning from the un-scaled pt.
                pt1 = self._scale_unapply_pt(pt_obs, eta_pm, phi_pm, q_pm, b_pm)
                cond1 = self._cond_from_muons(pt1, eta_pm, phi_pm, q_pm)
                mk_g = cond1.unsqueeze(1).expand(B, G, cond1.shape[-1]).contiguous()
            else:
                mk_g[..., N_MUON_KIN - 1] = self._source_rho_std_gh(
                    m_obs, pt_obs, eta_pm, phi_pm, q_pm, b_pm, xi)
        logp0 = self.log_p_nominal(
            mp.reshape(-1), mk_g.reshape(B * G, -1)).reshape(B, G)

        # log p_θ(x) = log|dm_s/dx| + logsumexp_i[ logW_i + log p_0(m'_i)
        # − log|G'_i| ] — the explicit backward-scale Jacobian times the
        # smear-only GH mixture at the scaled coordinate.
        log_terms = logW.view(1, G) + logp0 - torch.log(Gp.abs().clamp_min(1e-6))
        log_p_theta = torch.logsumexp(log_terms, dim=1) + log_J_s
        # Same upper-cap and NaN guard as the PF-ODE branch — see the
        # _continuity_logp_pf_ode tail for rationale.
        return torch.nan_to_num(log_p_theta.clamp(max=50.0), nan=0.0)

    def _source_rho_std_gh(self, m_obs, pt_obs, eta_pm, phi_pm, q_pm, b_pm, xig):
        """Standardised source ρ per GH node ``[B, G]`` for the
        ``gh_convolution`` operator's conditioning. Un-applies the scale (event-
        level, no GH dependence) and the smear's conditional-mean per-muon qop
        shift (per GH node ε = ξ_i, via ``E[δqop_μ | δm] = (J^m_μ σ²_qop,μ/√V)·ε``,
        with ``J^m_μ = ∂m/∂qop_μ = −m/(2 qop_μ)``). η, φ are exactly
        transform-invariant and don't need propagating."""
        B = m_obs.shape[0]
        G = xig.shape[-1]
        # 1) un-apply the scale (obs → truth), per event.
        pt1 = pt_obs
        if self.scale_enabled:
            AeM_pm = self._scale_AeM_pm(eta_pm, phi_pm, b_pm)
            dqop_s = self._delta_qop_analytic(AeM_pm, pt_obs, eta_pm, q_pm)
            pt1 = self._apply_scale_pt(pt_obs, eta_pm, q_pm, dqop_s, sign=-1.0)
        # 2) un-apply the smear's conditional-mean qop shift, per GH node.
        if self.smearing_enabled:
            sinth = _sintheta_from_eta(eta_pm)
            qop1 = q_pm * sinth / pt1                          # [B, 2]
            sig2 = self._qop_var_pm(
                eta_pm, phi_pm, b_pm, pt1).clamp_min(0.0)      # [B, 2] σ²_qop at source pt
            Jm = -0.5 * m_obs.unsqueeze(-1) / qop1             # [B, 2]  ∂m/∂qop
            V = (Jm * Jm * sig2).sum(-1, keepdim=True).clamp_min(1e-12)  # [B, 1]
            coef = Jm * sig2 / V.sqrt()                        # [B, 2]
            dqop_sm = coef.unsqueeze(1) * xig.view(1, G, 1)    # [B, G, 2]
            pt_src = self._apply_scale_pt(
                pt1.unsqueeze(1).expand(B, G, 2),
                eta_pm.unsqueeze(1).expand(B, G, 2),
                q_pm.unsqueeze(1).expand(B, G, 2),
                dqop_sm, sign=-1.0)                            # [B, G, 2]
            rho_src = (pt_src[..., 0] - pt_src[..., 1]) / (
                pt_src[..., 0] + pt_src[..., 1])
        else:
            r = (pt1[:, 0] - pt1[:, 1]) / (pt1[:, 0] + pt1[:, 1])
            rho_src = r.unsqueeze(1).expand(B, G)
        idx = N_MUON_KIN - 1
        return (rho_src - self.muon_kin_mean[idx]) / self.muon_kin_std[idx]

    # ------------------------------------------------------------------
    # Exact per-muon qop-space smear operator (2D Gauss–Hermite)
    # ------------------------------------------------------------------

    @staticmethod
    def _flow_eval_chunked(fn, m_flat, mk_flat):
        """Apply a per-row flow eval ``fn(m, mk)`` over the (flattened) first
        dim in chunks of ``_FLOW_EVAL_CHUNK`` rows and concatenate. The flow is
        row-independent so this is EXACT; it just caps the single-call
        allocation (the qop operator's B·n_gh² rows can otherwise trip the CUDA
        allocator on one matmul). Gradients flow normally through the chunks.

        BYPASSED under torch.compile (--fit-compile): the python chunk loop
        (a) UNROLLS — duplicating the flow subgraph once per chunk, so graph
        size and compile time grow ∝ batch size — and (b) SPECIALISES the
        graph on the exact row count, so every batch with a different
        data-mask count triggers a FULL recompile (observed as 'compiling on
        every batch' at large --fit-batch-size) until dynamo's recompile cap
        silently falls back to eager. The chunking is only an eager-mode
        allocation cap; under --fit-compile the single-call allocation is
        governed by --fit-batch-size directly."""
        n = m_flat.shape[0]
        if torch.compiler.is_compiling() or n <= _FLOW_EVAL_CHUNK:
            return fn(m_flat, mk_flat)
        return torch.cat([
            fn(m_flat[i:i + _FLOW_EVAL_CHUNK], mk_flat[i:i + _FLOW_EVAL_CHUNK])
            for i in range(0, n, _FLOW_EVAL_CHUNK)])

    def _gh_qop_unsmear(self, pt_cfg, etao, phio, qo, bpo, eps,
                        with_dlam: bool = False):
        """Un-kick an OBSERVED per-muon pt config to the nominal (truth) mass +
        ρ, per 2-D GH node. The model composes (generative, MC → data):

            1. forward smear at the truth pt:
               ``qop_sm,μ = qop_truth,μ + σ_qop,μ(pt_truth)·ξ_μ``
            2. forward scale — the implicit inverse of the DEFINING backward
               (data → MC) map ``qop_sm,μ = qop_obs,μ − δqop_μ(pt_obs)``
               (δqop at the OBSERVED pt; ``_inject_pt_np`` realises this exact
               forward by fixed point on the injection side).

        so the fit-side inversion here is, per node (ξ₊, ξ₋):

            1. BACKWARD scale, exact and closed-form (the defining direction —
               no iteration): ``qop_s,μ = qop_obs,μ − δqop_μ(pt_obs)``;
            2. un-smear with σ at the truth pt (the smear keeps its forward
               MC → data convention): the short fixed point
               ``qop_truth,μ = qop_s,μ − σ_qop,μ(pt_truth)·ξ_μ``,
               ``pt_truth,μ = |sinθ_μ / qop_truth,μ|`` — smear shift ONLY.

        No ordering ambiguity: the backward scale is applied first, and the
        flow/un-smear act on the scaled coordinates. pt = |sinθ/qop| is a
        magnitude, so a kick that flips the sign of qop is the physical charge
        mis-reco, kept (only the qop=0 pole guarded). The smear fixed point
        converges in a few iterations since σ ≪ |qop|; with smearing disabled
        the whole un-kick is a single closed-form expression.

        ``pt_cfg`` [B, G², 2] is the observed config (the mass-direction scaling
        of the old Jacobian leaf λ is evaluated at λ=1); ``eps`` [·, G², 2] the
        kicks. Returns ``(m_t [B, G²], pt_truth [B, G², 2])``; with
        ``with_dlam=True`` additionally returns the CLOSED-FORM per-muon
        ``dpt_truth/dλ`` [B, G², 2] for the change-of-variables Jacobian
        (λ scales the observed config along the mass direction, pt_cfg = pt·λ,
        evaluated at λ=1 — what the old autograd λ-leaf differentiated):

            dqop_obs/dλ = −qop_obs
            dδqop/dλ    = q·sinθ·(A − 2e·k)·(−k)            (∂δ/∂k · dk/dλ)
            dqop_t/dqop_s = 1/(1 + ξ·∂σ_eff/∂qop_t)         (implicit fn. thm.,
                ∂σ_eff/∂qop_t = [σ²>0]·c·k_t·sign(qop_t)/(sinθ·σ_eff))
            dpt/dqop    = −pt/qop  (zero where the qop=0 pole clamp engaged)

        The implicit-function smear factor is evaluated at the 3-step truncated
        fixed point rather than the exact root — a relative O(contraction³)
        ≈ 1e-8 mismatch vs the unrolled-autograd Jacobian, far below quadrature
        truncation.

        (+EPS inside the σ sqrt keeps the gradient finite as σ²→0, else
        d√v/dv→∞ times ∂v/∂λ→0 gives 0·∞ = NaN → a spurious uniform density at
        c≈0; EPS=1e-14 → σ floor 1e-7, negligible vs ~1e-3.)"""
        sinth = _sintheta_from_eta(etao)                          # [B,·,2]
        qop_obs = qo * sinth / pt_cfg                             # [B,G²,2] (bcast)
        dqop_dlam = -qop_obs if with_dlam else None               # dqop_obs/dλ|_{λ=1}

        # 1) BACKWARD (data → MC) scale — the defining direction: δqop at the
        #    OBSERVED config, explicit, exact, iteration-free.
        if self.scale_enabled:
            AeM = self._scale_AeM_pm(etao, phio, bpo)
            qop_s = qop_obs - self._delta_qop_analytic(AeM, pt_cfg, etao, qo)
            pt_s = self._qop_new_to_pt(qop_obs, qop_s, qo, sinth)
            if with_dlam:
                k = 1.0 / pt_cfg
                # dδqop/dλ = ∂δ/∂k · dk/dλ with ∂δ/∂k = q·sinθ·(A − 2e·k),
                # dk/dλ = −k at λ=1.
                dqop_dlam = dqop_dlam - (qo * sinth
                                         * (AeM[..., 0] - 2.0 * AeM[..., 1] * k)
                                         * (-k))
        else:
            qop_s, pt_s = qop_obs, pt_cfg
        qop_fin = qop_s if self.scale_enabled else qop_obs

        # 2) Un-smear from the scaled coordinates: σ at the truth pt (forward
        #    MC → data smear convention unchanged) → fixed point over the smear
        #    shift alone.
        pt_truth = pt_s
        if self.smearing_enabled:
            qop_truth = qop_s
            for _ in range(3):
                sig = (self._qop_var_pm(etao, phio, bpo, pt_truth).clamp_min(0.0)
                       + 1e-14).sqrt()
                qop_truth = qop_s - sig * eps
                pt_truth = self._qop_new_to_pt(qop_s, qop_truth, qo, sinth)
            qop_fin = qop_truth
            if with_dlam:
                # Implicit-function Jacobian of the smear stage at the truncated
                # fixed point: dqop_t/dqop_s = 1/(1 + ξ·∂σ_eff/∂qop_t).
                ac = self._smear_ac_pm(etao, phio, bpo)           # [B,·,2,2]
                c_phys = ac[..., 1] * SMEAR_VAR_SCALE_C
                var_raw = self._qop_var_pm(etao, phio, bpo, pt_truth)
                sig_eff = (var_raw.clamp_min(0.0) + 1e-14).sqrt()
                dsig = ((var_raw > 0).to(sig_eff.dtype) * c_phys
                        * torch.sign(qop_truth)
                        / (pt_truth * sinth * sig_eff))           # c·k_t·sign/(sinθ·σ)
                dqop_dlam = dqop_dlam / (1.0 + eps * dsig)
        m_t = _event_mll(pt_truth, etao, phio)                   # [B,G²]
        # Return the un-kicked nominal momenta; the caller builds the per-node
        # conditioning from them via _node_cond (ρ-only for muon_kin, the full
        # event-level vector for event_level).
        if not with_dlam:
            return m_t, pt_truth
        # dpt/dqop = −pt/qop (pt = sinθ/|qop|); zero where the qop=0 pole
        # clamp froze pt.
        pole = (qop_fin.abs() > QOP_EPS).to(pt_truth.dtype)
        dpt_dlam = -(pt_truth / qop_fin) * dqop_dlam * pole
        return m_t, pt_truth, dpt_dlam

    def _continuity_logp_gh_qop(self, m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm,
                                b_pm, n_gh: int = 8, n_iter: int = 2):
        """``log p_θ(x|c)`` via the EXACT per-muon qop smear, as a 2-D
        Gauss–Hermite quadrature over the two INDEPENDENT muon kicks (ξ₊, ξ₋),
        implemented as a DIRECT un-kick of the observed config:

            ``p_θ(x|c) = E_{ε₊,ε₋}[ p_0(m_t(ε₊,ε₋) | ρ_t) · |∂m_t/∂x| ]``
            ``       ≈ Σ_{i,j} W_i W_j · p_0(m_t,ij | ρ_t,ij) · |∂m_t/∂x|_{ij}``

        For each node, the observed per-muon config (the event's pt, scaled
        along the mass direction by a leaf ``λ``) is inverted through the
        COMBINED qop kick ``qop_obs = qop_truth + δqop + σ·ξ`` (scale shift and
        smear width both at the truth pt) to the nominal config → nominal mass
        ``m_t`` and ρ (see ``_gh_qop_unsmear``). ``p_0`` is the frozen flow at
        that nominal mass conditioned on the per-node nominal ρ; the change-of-
        variables Jacobian ``|∂m_t/∂x|`` is CLOSED FORM (analytic ∂m/∂pt and the
        un-kick chain derivative — see ``_gh_qop_unsmear`` with_dlam). This
        is the genuine inverse of the per-muon fold that generates the pseudo-
        data (the GH nodes live in the 2-D qop±-kick space), reproducing its full
        non-Gaussian shape — mean-shift and skew, not just the variance.

        Unlike ``gh_convolution`` (a single Gaussian convolution in MASS space,
        the small-kick linearisation), the kicks can't be factorised into two
        1-D integrals (``p_0`` is nonlinear in the combined config), so this is
        a genuine n_gh² quadrature. Each per-muon kick is a single 1-D Gaussian,
        so a modest ``n_gh`` (≈4–6) suffices. (``n_iter`` is unused — the
        un-kick runs its own short joint fixed point internally.)"""
        B = m_obs.shape[0]
        # --disable-smearing: no kick to integrate, so collapse the n_gh² GH
        # quadrature to a single node (ξ=0, W=1). _gh_qop_unsmear then does pure
        # scale un-application and the G²=1 logsumexp reduces to log p — the
        # right scale-only fallback (and avoids the n_gh²× redundant flow evals).
        ng = n_gh if self.smearing_enabled else 1
        xi, logW = _gh_nodes(ng, m_obs.device, m_obs.dtype)            # [G],[G]
        G = xi.shape[0]
        G2 = G * G
        xi_p = xi.view(G, 1).expand(G, G).reshape(G2)                   # ξ₊
        xi_m = xi.view(1, G).expand(G, G).reshape(G2)                   # ξ₋
        logW2 = (logW.view(G, 1) + logW.view(1, G)).reshape(G2)        # [G²]
        eps = torch.stack([xi_p, xi_m], dim=-1).view(1, G2, 2)         # [1,G²,2]

        pto = pt_obs.unsqueeze(1)                                       # [B,1,2]
        etao = eta_pm.unsqueeze(1)                                      # [B,1,2]
        phio = phi_pm.unsqueeze(1)                                      # [B,1,2]
        qo = q_pm.unsqueeze(1)                                          # [B,1,2]
        bpo = b_pm.unsqueeze(1)                                         # [B,1,2]
        # The observed config is self-consistent (mll = _event_mll(pt_obs)) —
        # the loader propagates the smeared pt + recomputed mll + ρ for the
        # validation pseudo-data, exactly as for real data — so the density is
        # evaluated at _event_mll(pt_obs) = m_obs directly, no rescaling needed.

        # d m_t / d m_obs Jacobian in CLOSED FORM (no autograd λ-leaf):
        # perturbing m_obs scales the observed config along pt∝m (fixed obs ρ),
        # so with pt_cfg = pt_obs·λ at λ=1, |∂m_t/∂x| = |∂m_t/∂λ|/|∂m_obs/∂λ|
        # with ∂m/∂λ = Σ_μ (∂m/∂pt_μ)·(dpt_μ/dλ) — ∂m/∂pt analytic
        # (_dm_dpt_analytic) and dpt_t/dλ from the un-kick chain
        # (_gh_qop_unsmear with_dlam; dpt_cfg/dλ = pt_cfg). θ-gradients flow
        # through the closed-form expressions by ordinary autograd.
        pt_cfg = pto.expand(B, G2, 2)                                   # [B,G²,2]
        m_t, pt_truth, dpt_dlam = self._gh_qop_unsmear(
            pt_cfg, etao, phio, qo, bpo, eps, with_dlam=True)
        g_mt = (_dm_dpt_analytic(pt_truth, etao, phio) * dpt_dlam).sum(-1)
        g_mo = (_dm_dpt_analytic(pt_cfg, etao, phio) * pt_cfg).sum(-1)
        # |∂m_t/∂x| = |∂m_t/∂λ| / |∂m_obs/∂λ|.
        logJ = (torch.log(g_mt.abs().clamp_min(1e-12))
                - torch.log(g_mo.abs().clamp_min(1e-12)))
        mk_g = mk.unsqueeze(1).expand(B, G2, mk.shape[-1]).clone()
        if self.scale_enabled or self.smearing_enabled:
            mk_g = self._node_cond(mk_g, pt_truth, etao, phio, qo)
        logp0 = self._flow_eval_chunked(
            self.log_p_nominal, m_t.reshape(-1), mk_g.reshape(B * G2, -1)).reshape(B, G2)
        log_terms = logW2.view(1, G2) + logp0 + logJ
        log_p_theta = torch.logsumexp(log_terms, dim=1)
        return torch.nan_to_num(log_p_theta.clamp(max=50.0), nan=0.0)

    def _log_jacobian_exp(self, mp, mk_src, s_adv_of, eta_pm, phi_pm, b_pm,
                          pt_obs, m_obs):
        """Frozen-score continuous-flow log-Jacobian:
        ``log G' = log(1 + s_adv'(m')) − V·∂²_m log p₀(m')/2``.

        Both ``s_adv'(m')`` and ``∂²_m log p₀(m')`` are taken via autograd at
        the SOURCE mass m' (frozen along the smear trajectory). Always finite
        — no floor — but a different operator approximation than the autograd
        Jacobian of the N-step Euler forward map (see _continuity_logp
        docstring). Live grads to θ are preserved when ``mp.requires_grad``.

        Both factors are guarded with ``nan_to_num``: ``∂²_m log p₀`` exposes
        the flow's tail singularities (where the Jacobian floor in the GF
        monotonic transform produces a step) more directly than the softlog
        path's autograd Gp (which combines everything through the chain rule
        — sharp pieces of the score and its derivative tend to cancel). A
        small minority of events at the m-window boundary can return inf/NaN
        d2, and even one such event poisons the batch mean. We mask those
        contributions to zero rather than dropping the events — the model is
        at the edge of its validity there anyway."""
        # ∂s_adv/∂m' (scale-advection Jacobian contribution). With the BACKWARD
        # scale the operators handle the scale Jacobian explicitly at the
        # observed mass (s_adv_of=None here → smear-only Jacobian).
        if self.scale_enabled and s_adv_of is not None:
            if mp.requires_grad:
                s = s_adv_of(mp)
                s_prime = torch.autograd.grad(
                    s.sum(), mp, create_graph=True, retain_graph=True)[0]
            else:
                with torch.enable_grad():
                    mpj = mp.detach().requires_grad_(True)
                    s_prime = torch.autograd.grad(
                        s_adv_of(mpj).sum(), mpj)[0].detach()
            log_J_scale = torch.nan_to_num(
                torch.log1p(s_prime), nan=0.0, posinf=0.0, neginf=0.0)
        else:
            log_J_scale = mp.new_zeros(mp.shape)
        # −V·∂²_m log p₀(m')/2 (smear contribution, =0 if smear disabled).
        if self.smearing_enabled:
            V = self._smear_mass_var(eta_pm, phi_pm, b_pm, pt_obs, mp, m_obs)
            if mp.requires_grad:
                score = self._flow_score(mp, mk_src)
                d2 = torch.autograd.grad(
                    score.sum(), mp, create_graph=True, retain_graph=True)[0]
            else:
                with torch.enable_grad():
                    mpj = mp.detach().requires_grad_(True)
                    score = self._flow_score(mpj, mk_src)
                    d2 = torch.autograd.grad(score.sum(), mpj)[0].detach()
            log_J_smear = torch.nan_to_num(
                -0.5 * V * d2, nan=0.0, posinf=0.0, neginf=0.0)
        else:
            log_J_smear = mp.new_zeros(mp.shape)
        return log_J_scale + log_J_smear

    def _flow_normal_base_z(self, m_std: torch.Tensor,
                            mk: torch.Tensor) -> torch.Tensor:
        """Base-space latent ``z = f(m_std | c)`` for the normal-base archs
        (CDF F = Φ(z)). ``mixture`` uses its analytic ``transform`` (no autograd,
        compile-traceable); gf/nsf use zuko's monotonic ``dist.transform``."""
        if getattr(self, "flow_is_mixture", False):
            return self.flow.transform(m_std, mk)
        dist = self.flow.flow(mk)
        return dist.transform(m_std.unsqueeze(-1)).squeeze(-1)

    def _flow_log_cdf(self, m: torch.Tensor, mk: torch.Tensor) -> torch.Tensor:
        """``log F_0(m | mk)`` — the FLOW's CDF at observed mass ``m``.

        For the GF flow with standard-normal base, ``F_0(m|c) =
        Φ(T_mono(m_std|c))`` where ``T_mono`` is the conditional Gaussianisation
        flow's monotonic transform (zuko's ``dist.transform``) acting on the
        standardised mass. The standardisation is monotone increasing so the
        CDF transforms trivially: ``F_data(m) = F_data_std(m_std)``. Returns
        ``[B]`` log-CDF values clamped from below for log safety."""
        m_std = self._standardise_mll(m).clamp(
            -MLL_STD_FLOW_CLAMP, MLL_STD_FLOW_CLAMP)
        if getattr(self, "flow_is_compact", False):
            # Compact flow: the cumulative is built in (and includes the matched
            # tails), so the window-integral differences F(unkick_hi)-F(unkick_lo)
            # are exact. No Φ — the transform already IS the (renormalised) CDF.
            return self.flow.log_cdf(m_std, mk)
        if getattr(self, "flow_is_nce", False):
            # NCE density: no analytic CDF — per-event Gauss-Legendre quadrature
            # F₀(x) = 1 + ∫_a^x p̂₀. The +1 offset (positivity for the small
            # below-window probes) makes this valid for CDF DIFFERENCES only,
            # which is how every consumer uses it (window Z, stage-1 norm,
            # display norm).
            return self.flow.log_cdf(m_std, mk)
        if getattr(self, "flow_is_dcb", False) or getattr(self, "flow_is_ege", False):
            # DCB/EGE: fully analytic piecewise CDF (Φ core + power-law /
            # exponential tail antiderivatives), normalised so
            # F₀(b)−F₀(a)=1 exactly.
            return self.flow.log_cdf(m_std, mk)
        z = self._flow_normal_base_z(m_std, mk)
        F = 0.5 * (1.0 + torch.erf(z / math.sqrt(2.0)))
        return F.clamp(min=1e-30).log()

    @staticmethod
    def _log1mexp(d: torch.Tensor) -> torch.Tensor:
        """``log(1 − e^d)`` for d ≤ 0, stable for both d→0⁻ and d→−∞
        (Mächler's two-branch form)."""
        d = d.clamp(max=-1e-12)          # coincident bounds → ~log(1e-12), not −inf
        return torch.where(d > -0.6931471805599453,          # −log 2
                           torch.log(-torch.expm1(d)),
                           torch.log1p(-torch.exp(d)))

    def _log_phi_window(self, z_lo: torch.Tensor, z_hi: torch.Tensor) -> torch.Tensor:
        """``log(Φ(z_hi) − Φ(z_lo))`` (z_hi ≥ z_lo), fp64-stable in BOTH tails
        AND for near-coincident bounds:

        * wide bounds: reflect to the side where Φ is small (log_ndtr is
          accurate there) and subtract in log space — exact down to window
          masses ~e−300 (far-tail windows), no exp-difference cancellation.
        * narrow bounds (z_hi − z_lo < 1e−6, the collapsed-transform regime —
          the log-space subtraction amplifies the log_ndtr ulp by 1/width):
          midpoint form log Z = log φ(z̄) + log(width), exact to O(width²),
          with a finite floor and a NONZERO width gradient — so a training
          run whose transform collapses the window gets a correct restoring
          gradient instead of a frozen clamp."""
        zl = z_lo.double(); zh = z_hi.double()
        width = (zh - zl).clamp_min(0.0)
        # wide branch: log-space subtraction on the accurate tail side
        flip = (zl + zh) > 0             # window in the upper Φ tail → reflect
        a = torch.where(flip, -zh, zl)
        b = torch.where(flip, -zl, zh)
        log_Fb = torch.special.log_ndtr(b)
        log_Fa = torch.special.log_ndtr(a)
        log_wide = log_Fb + self._log1mexp(log_Fa - log_Fb)
        # narrow branch: midpoint × width
        zbar = 0.5 * (zl + zh)
        log_narrow = (-0.5 * zbar * zbar - 0.5 * math.log(2.0 * math.pi)
                      + width.clamp_min(1e-290).log())
        out = torch.where(width < 1e-6, log_narrow, log_wide)
        return out.to(z_lo.dtype)

    def _flow_log_window_Z(self, m_lo_t: torch.Tensor, m_hi_t: torch.Tensor,
                           mk: torch.Tensor) -> torch.Tensor:
        """STABLE per-row ``log ∫_{m_lo_t}^{m_hi_t} p₀`` = log(F₀(hi) − F₀(lo)).

        For the normal-base archs (gf/nsf) the naive ``F.exp() − F.exp()``
        difference collapses once the window mass Z(c) falls below the fp
        cancellation/underflow floor — which the truncated stage-1 training
        ALLOWS, since the out-of-window gauge is free and Z(c) can drift
        arbitrarily small (observed at e−40-scale for forward-η conditioning).
        The collapsed Z then (a) explodes the display normalisation by
        e⁶⁹·Z_true and (b) freezes the training/fit logZ at the clamp with
        ZERO gradient, silently switching those events to unnormalised-density
        maximisation. Here the window mass is computed from the base-space z
        values via fp64 ``log_ndtr`` on the accurate tail side — exact to
        Z ~ e−300, with correct gradients throughout.

        compact/nce/dcb keep their native O(1)-conditioned CDF differences
        (exact-Z construction / quadrature integral / analytic CDF)."""
        if (getattr(self, "flow_is_compact", False)
                or getattr(self, "flow_is_nce", False)
                or getattr(self, "flow_is_dcb", False)
                or getattr(self, "flow_is_ege", False)):
            lo = self._flow_log_cdf(m_lo_t, mk)
            hi = self._flow_log_cdf(m_hi_t, mk)
            return (hi.exp() - lo.exp()).clamp_min(1e-30).log()
        m_std_lo = self._standardise_mll(m_lo_t).clamp(
            -MLL_STD_FLOW_CLAMP, MLL_STD_FLOW_CLAMP)
        m_std_hi = self._standardise_mll(m_hi_t).clamp(
            -MLL_STD_FLOW_CLAMP, MLL_STD_FLOW_CLAMP)
        z_lo = self._flow_normal_base_z(m_std_lo, mk)
        z_hi = self._flow_normal_base_z(m_std_hi, mk)
        return self._log_phi_window(z_lo, z_hi)

    def _flow_log_window_Z_chunked(self, m_lo_flat, m_hi_flat, mk_flat):
        """Chunked ``_flow_log_window_Z`` (same chunking as _flow_eval_chunked,
        incl. the torch.compile bypass — see there for the recompile story)."""
        n = m_lo_flat.shape[0]
        if torch.compiler.is_compiling() or n <= _FLOW_EVAL_CHUNK:
            return self._flow_log_window_Z(m_lo_flat, m_hi_flat, mk_flat)
        return torch.cat([
            self._flow_log_window_Z(m_lo_flat[i:i + _FLOW_EVAL_CHUNK],
                                    m_hi_flat[i:i + _FLOW_EVAL_CHUNK],
                                    mk_flat[i:i + _FLOW_EVAL_CHUNK])
            for i in range(0, n, _FLOW_EVAL_CHUNK)])

    def _fit_cut_m_min(self, m_obs, pt_obs, eta_pm, phi_pm) -> torch.Tensor:
        """Per-event LOWER mass edge ``[B]`` set by the fit-time reco pt cuts.

        At fixed conditioning every pt magnitude scales with a single overall
        factor α (pt ∝ α), so a cut ``C`` on a pt quantity ``q`` is met exactly
        at ``α = C/q_obs``; the binding scale is ``α_min = max_c C_c/q_obs,c``
        and the forbidden region is ``m < m_min = mll(α_min · pt_obs)`` —
        evaluated with the FULL muon-mass kinematics (``_event_mll``), so the
        edge is exact (no massless pt∝m approximation). Floored at the window
        ``m_lo`` and capped just below ``m_hi``. With no fit-time cuts this
        returns the constant ``m_lo`` (exact legacy path).

        Events that pass the cuts have ``m_obs ≥ m_min`` by construction, so the
        per-event window ``[m_min, m_hi]`` is non-empty; the clamp only guards
        numerical edge cases."""
        if not (self._fit_ptll_min or self._fit_pt_lead_min
                or self._fit_pt_both_min):
            return m_obs.new_full(m_obs.shape, self._m_lo_f)
        alpha = m_obs.new_zeros(m_obs.shape)
        if self._fit_ptll_min:
            px = pt_obs[:, 0] * torch.cos(phi_pm[:, 0]) + \
                pt_obs[:, 1] * torch.cos(phi_pm[:, 1])
            py = pt_obs[:, 0] * torch.sin(phi_pm[:, 0]) + \
                pt_obs[:, 1] * torch.sin(phi_pm[:, 1])
            ptll = torch.sqrt((px * px + py * py).clamp_min(1e-12))
            alpha = torch.maximum(alpha, self._fit_ptll_min / ptll)
        if self._fit_pt_lead_min:
            lead = torch.maximum(pt_obs[:, 0], pt_obs[:, 1])
            alpha = torch.maximum(alpha, self._fit_pt_lead_min / lead)
        if self._fit_pt_both_min:
            soft = torch.minimum(pt_obs[:, 0], pt_obs[:, 1])
            alpha = torch.maximum(alpha, self._fit_pt_both_min / soft)
        m_min = _event_mll(alpha.unsqueeze(-1) * pt_obs, eta_pm, phi_pm)
        return m_min.clamp(self._m_lo_f, self._m_hi_f - 1e-6)

    def _norm_correction_log_Z(self, m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm,
                               b_pm, n_iter: int = 2) -> torch.Tensor:
        """Per-event ``log Z(θ;c) = log ∫_{m_lo}^{m_hi} p_θ(x|c) dx`` — the
        normalisation of the transformed-flow density on the m-window. Bare
        ``log p_θ(x)`` is normalised on ``T(window)``, NOT on the window; this
        correction restores the proper truncated-likelihood NLL by subtracting
        ``log Z`` from the per-event signal log-density. Dispatches on
        ``self.norm_correction``:

        * ``"none"``: return 0 (no correction).
        * ``"linear"``: leading-order boundary-leakage estimate
          ``1 − Z ≈ p_0(m_lo|c)·(m_lo − T(m_lo))_+ + p_0(m_hi|c)·(T(m_hi) − m_hi)_+``
          — only the positive parts of the boundary shift contribute (a
          narrowing T pushes mass INTO the window, no leakage). Cheap: 2
          boundary forward-map evals + 2 flow density evals per event.
        * ``"flow_cdf"``: EXACT via the flow's CDF
          ``Z(θ;c) = F_0(T⁻¹(m_hi)|c) − F_0(T⁻¹(m_lo)|c)``. The scale part of
          T⁻¹ is the EXPLICIT backward (data → MC) boundary transport (defining
          direction — no inversion); only the smear is inverted (bisection),
          and not at all when smearing is disabled. ~2× the per-event work of
          the bare density.

        When ``self.smear_operator == "gh_convolution"`` both 'linear' and
        'flow_cdf' route to the GH-specific exact formula
        ``Z = Σ_i W_i [F_0(m'_hi(ξ_i)) − F_0(m'_lo(ξ_i))]`` (per-GH-node
        boundary inversion + flow CDF) — the PF-ODE forward map used by the
        modes below is the wrong operator for GH, so the simple boundary
        formulae would give the wrong V-scaling. See
        ``_norm_correction_log_Z_gh``.
        """
        if self.norm_correction == "none":
            return m_obs.new_zeros(m_obs.shape)
        if self.smear_operator == "gh_convolution":
            return self._norm_correction_log_Z_gh(
                m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm)
        if self.smear_operator == "gh_convolution_qop":
            return self._norm_correction_log_Z_gh_qop(
                m_obs, mk, pt_obs, eta_pm, phi_pm, q_pm, b_pm)

        B = m_obs.shape[0]
        theta_scale_pm = (self._scale_per_event(eta_pm, phi_pm, b_pm)
                          if self.scale_enabled
                          else m_obs.new_zeros((B, N_THETA_SCALE_PM)))
        n_step = max(1, int(self.smear_flow_steps))

        mk_src = mk
        if self.scale_enabled:
            if self.cond_basis == "event_level":
                pt1 = self._scale_unapply_pt(pt_obs, eta_pm, phi_pm, q_pm, b_pm)
                mk_src = self._cond_from_muons(pt1, eta_pm, phi_pm, q_pm)
            else:
                mk_src = mk.clone()
                mk_src[..., N_MUON_KIN - 1] = self._scale_source_rho_std(
                    pt_obs, eta_pm, phi_pm, q_pm, b_pm)

        def forward_sm(mp):
            # smear-only forward map (the scale is the explicit backward
            # boundary transport below).
            y = mp
            if self.smearing_enabled:
                V = self._smear_mass_var(eta_pm, phi_pm, b_pm, pt_obs, mp, m_obs)
                dt = V / (2.0 * n_step)
                for _ in range(n_step):
                    y = y - dt * self._flow_score(y, mk_src)
            return y

        m_lo = self._fit_cut_m_min(m_obs, pt_obs, eta_pm, phi_pm)
        m_hi = m_obs.new_full(m_obs.shape, self._m_hi_f)
        # BACKWARD (data → MC) scale images of the window boundaries along the
        # event's observed ray — explicit and exactly θ-differentiable (the
        # defining direction): no boundary fixed point / bisection for the
        # scale, one of the backward-direction simplifications.
        m_s_lo, _ = self._scale_backward_mass_linear(
            m_lo, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)
        m_s_hi, _ = self._scale_backward_mass_linear(
            m_hi, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)

        if self.norm_correction == "linear":
            # Boundary source preimages u = SmearInv(m_s(boundary)) by a short
            # fixed point (precision is secondary for the linear estimate);
            # 1 − Z ≈ p₀(m_lo)·(u_lo − m_lo)₊ + p₀(m_hi)·(m_hi − u_hi)₊ — only
            # sources pushed OUT of the window leak; narrowing contributes 0.
            u_lo, u_hi = m_s_lo, m_s_hi
            if self.smearing_enabled:
                for _ in range(n_iter):
                    u_lo = m_s_lo - (forward_sm(u_lo) - u_lo)
                    u_hi = m_s_hi - (forward_sm(u_hi) - u_hi)
            left_leak = (u_lo - m_lo).clamp(min=0.0)
            right_leak = (m_hi - u_hi).clamp(min=0.0)
            log_p_lo = self.log_p_nominal(m_lo, mk_src)
            log_p_hi = self.log_p_nominal(m_hi, mk_src)
            leakage = (log_p_lo.exp() * left_leak
                       + log_p_hi.exp() * right_leak)
            # clamp at a sub-1 ceiling so log(1−leakage) stays finite even if
            # the linear estimate overshoots in pathological cases
            return torch.log1p(-leakage.clamp(0.0, 1.0 - 1e-6))

        # "flow_cdf": Z = F₀(SmearInv(m_s(m_hi))) − F₀(SmearInv(m_s(m_lo))).
        # Scale-only: the preimage IS the explicit backward image — done. With
        # smearing: invert the SMEAR-ONLY map toward the backward-scaled
        # targets by BISECTION on a wider-than-window bracket (monotonic —
        # G' = 1 + V/(2σ²) > 0 at the peak; the per-event fixed point becomes
        # expansive at the boundaries when V/(2σ²) is large). 24 iterations
        # give ~22-bit precision — far below grid noise. The bisection runs
        # detached; the exact θ-gradient of the SCALE boundary transport is
        # re-attached through m_s (the smear preimage displacement stays
        # detached, as before).
        if not self.smearing_enabled:
            return self._flow_log_window_Z(m_s_lo, m_s_hi, mk_src)
        n_bisect = 26
        # Bracket WIDER than the window: a smear preimage of a (transported)
        # boundary can lie OUTSIDE [m_lo, m_hi]; a window-tight bracket (the
        # old behaviour) silently clamps it to the edge — a first-order logZ
        # error whenever the transport pushes the preimage out.
        marg = 0.25 * (self._m_hi_f - self._m_lo_f)
        @torch.no_grad()
        def bisect(target):
            lo = m_obs.new_full(m_obs.shape, self._m_lo_f - marg)
            hi = m_obs.new_full(m_obs.shape, self._m_hi_f + marg)
            for _ in range(n_bisect):
                mid = 0.5 * (lo + hi)
                t_mid = forward_sm(mid)
                go_right = t_mid < target
                lo = torch.where(go_right, mid, lo)
                hi = torch.where(go_right, hi, mid)
            return 0.5 * (lo + hi)
        mp_lo = m_s_lo - (m_s_lo.detach() - bisect(m_s_lo.detach()))
        mp_hi = m_s_hi - (m_s_hi.detach() - bisect(m_s_hi.detach()))
        # Stable log window mass (no exp-difference cancellation; coincident
        # preimages — T's image missing the window — give a large-negative
        # finite logZ via the _log1mexp guard rather than a hard clamp).
        return self._flow_log_window_Z(mp_lo, mp_hi, mk_src)

    def _norm_correction_log_Z_gh(self, m_obs, mk, pt_obs, eta_pm, phi_pm,
                                  q_pm, b_pm) -> torch.Tensor:
        """``log Z(θ;c)`` for the ``smear_operator='gh_convolution'`` path.

        The scale is the explicit BACKWARD (data → MC) boundary transport
        ``m_s(boundary) = boundary − s_adv(boundary)`` along the event's
        observed ray (defining direction — closed form, exactly
        θ-differentiable, no inversion). The smear keeps its forward map
        ``m_s = m' + √V(m')·ε``, ε ~ N(0,1), inverted toward the transported
        boundaries. Z over the m-window factors per GH node:

            ``Z = Σ_i W_i · [F_0(m'_hi(ξ_i)|c_src,i) − F_0(m'_lo(ξ_i)|c_src,i)]``

        where each m'_{lo,hi}(ξ_i) is the source mass at GH node ξ_i whose
        smear image is the scale-transported boundary, found by BISECTION on
        a wider-than-window bracket (a boundary fixed point oscillates at
        large V — same
        reason as in `_norm_correction_log_Z`'s flow_cdf branch; the scale's
        exact θ-gradient is re-attached through m_s, the smear displacement
        stays detached). With smearing disabled the preimage IS the explicit
        backward image — no bisection at all. Per-GH-node CDF evaluation uses
        the source ρ from `_source_rho_std_gh` for consistency with
        `_continuity_logp_gh`. Cost: 2·n_gh boundary bisections + 2·n_gh
        flow CDF evals per event — roughly 2× the GH density itself.
        """
        B = m_obs.shape[0]
        theta_scale_pm = (self._scale_per_event(eta_pm, phi_pm, b_pm)
                          if self.scale_enabled
                          else m_obs.new_zeros((B, N_THETA_SCALE_PM)))
        xi, logW = _gh_nodes(self.n_gh_nodes, m_obs.device, m_obs.dtype)
        G = xi.shape[0]

        mo = m_obs.unsqueeze(1)
        pto = pt_obs.unsqueeze(1)
        etao = eta_pm.unsqueeze(1)
        phio = phi_pm.unsqueeze(1)
        qo = q_pm.unsqueeze(1)
        bpo = b_pm.unsqueeze(1).expand(B, G, b_pm.shape[-1])
        xig = xi.view(1, G)

        # BACKWARD scale images of the boundaries (per event, node-independent).
        m_lo_t = self._fit_cut_m_min(m_obs, pt_obs, eta_pm, phi_pm)
        m_hi_t = m_obs.new_full(m_obs.shape, self._m_hi_f)
        m_s_lo, _ = self._scale_backward_mass_linear(
            m_lo_t, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)
        m_s_hi, _ = self._scale_backward_mass_linear(
            m_hi_t, m_obs, pt_obs, eta_pm, q_pm, theta_scale_pm)

        def forward_at_eps_sm(me):
            """Smear-only T(me; ξ) per (event, GH node). me [B, G] → [B, G].
            The +EPS inside the sqrt matches the density's `_smear_disp` so the
            boundary preimages are consistent with `_continuity_logp_gh`."""
            V = self._smear_mass_var(
                etao, phio, bpo, pto, me, mo).clamp_min(0.0)
            return me + (V + 1e-12).sqrt() * xig

        if not self.smearing_enabled:
            mp_lo = m_s_lo.unsqueeze(1).expand(B, G)       # [B, G]
            mp_hi = m_s_hi.unsqueeze(1).expand(B, G)
        else:
            n_bisect = 26
            # Wider-than-window bracket: smear preimages of the transported
            # boundaries can lie outside [m_lo, m_hi] (window-tight brackets
            # clamp them to the edge — first-order logZ error; see the
            # _norm_correction_log_Z flow_cdf branch).
            marg = 0.25 * (self._m_hi_f - self._m_lo_f)

            @torch.no_grad()
            def bisect(target):                            # target [B, 1]
                lo = m_obs.new_full((B, G), self._m_lo_f - marg)
                hi = m_obs.new_full((B, G), self._m_hi_f + marg)
                for _ in range(n_bisect):
                    mid = 0.5 * (lo + hi)
                    t_mid = forward_at_eps_sm(mid)
                    go_right = t_mid < target
                    lo = torch.where(go_right, mid, lo)
                    hi = torch.where(go_right, hi, mid)
                return 0.5 * (lo + hi)

            # Re-attach the scale's exact θ-gradient through m_s; the smear
            # preimage displacement (bisected) stays detached, as before.
            mp_lo = (m_s_lo.unsqueeze(1)
                     - (m_s_lo.detach().unsqueeze(1)
                        - bisect(m_s_lo.detach().unsqueeze(1))))   # [B, G]
            mp_hi = (m_s_hi.unsqueeze(1)
                     - (m_s_hi.detach().unsqueeze(1)
                        - bisect(m_s_hi.detach().unsqueeze(1))))   # [B, G]

        # Per-GH-node source ρ (same construction as _continuity_logp_gh).
        mk_g = mk.unsqueeze(1).expand(B, G, mk.shape[-1]).clone()
        if self.scale_enabled or self.smearing_enabled:
            if self.cond_basis == "event_level":
                # smear is guarded to the qop operator, so this is scale-only:
                # the full event-level conditioning from the un-scaled pt.
                pt1 = self._scale_unapply_pt(pt_obs, eta_pm, phi_pm, q_pm, b_pm)
                cond1 = self._cond_from_muons(pt1, eta_pm, phi_pm, q_pm)
                mk_g = cond1.unsqueeze(1).expand(B, G, cond1.shape[-1]).contiguous()
            else:
                mk_g[..., N_MUON_KIN - 1] = self._source_rho_std_gh(
                    m_obs, pt_obs, eta_pm, phi_pm, q_pm, b_pm, xi)
        # Per-node STABLE log window mass, then the GH-weighted sum in log
        # space (logsumexp) — same maths as Σ W·(F_hi − F_lo) but immune to
        # the exp-difference cancellation when the window mass is tiny.
        log_nw = self._flow_log_window_Z_chunked(
            mp_lo.reshape(-1), mp_hi.reshape(-1),
            mk_g.reshape(B * G, -1)).reshape(B, G)                      # [B, G]
        return torch.logsumexp(logW.view(1, G) + log_nw, dim=1)         # [B]

    def _pt_lambda_to_mass(self, pto, etao, phio, m_obs, m_target,
                           n_iter: int = 3):
        """Per-event pt scale ``λ`` such that ``_event_mll(pt·λ) = m_target`` at
        fixed (η, φ) — the muon-mass-EXACT inverse of pt∝m.

        Newton from the massless guess ``λ₀ = m_target/m_obs``; the closed-form
        ``∂m_ll/∂λ = (∂m/∂pt · pt)/λ`` uses ``_dm_dpt_analytic``. Shapes: pto/
        etao/phio ``[...,2]``; m_obs / m_target broadcastable to the leading dims.
        3 iters → median ~3e-5 GeV (the muon-mass correction is ~1e-3, Newton is
        quadratic). Replaces the naive massless ``m_target/m_obs`` scaling whose
        boundary config lands at an observed mass off by ~(m_μ/m)², shrinking the
        norm window ~0.035% → signal renorm ~0.5% high → bkg-fraction bias."""
        lam0 = m_target / m_obs
        # the true correction is < 0.5%; bound the iterate to ±5% of the massless
        # guess so rare near-collinear events (∂m/∂pt → 0, Newton overshoots)
        # stay bounded instead of diverging — a no-op for the bulk.
        lo, hi = lam0 * 0.95, lam0 * 1.05
        lam = lam0
        for _ in range(n_iter):
            ptl = pto * lam.unsqueeze(-1)
            mll = _event_mll(ptl, etao, phio)
            dmdl = ((_dm_dpt_analytic(ptl, etao, phio) * ptl).sum(-1)
                    / lam.clamp_min(1e-12))
            lam = torch.maximum(torch.minimum(
                lam - (mll - m_target) / dmdl.clamp_min(1e-12), hi), lo)
        # guarantee a STRICT improvement over the massless guess for every event
        # (the rare overshooters fall back to ≈massless, never worse):
        err_n = (_event_mll(pto * lam.unsqueeze(-1), etao, phio) - m_target).abs()
        err_0 = (_event_mll(pto * lam0.unsqueeze(-1), etao, phio) - m_target).abs()
        return torch.where(err_n <= err_0, lam, lam0)

    def _norm_correction_log_Z_gh_qop(self, m_obs, mk, pt_obs, eta_pm, phi_pm,
                                       q_pm, b_pm) -> torch.Tensor:
        """``log Z(θ;c)`` for ``smear_operator='gh_convolution_qop'``, in the
        DIRECT un-kick formulation (matching ``_continuity_logp_gh_qop``).

        Swapping the window integral and the GH expectation and substituting
        ``u = m_t(x,ξ)`` (the un-kick map):

            ``Z = Σ_{i,j} W_i W_j [F_0(m_t,hi(ξ) | ρ_t) − F_0(m_t,lo(ξ) | ρ_t)]``

        where ``m_t,lo/hi(ξ)`` are the nominal masses obtained by un-kicking the
        observed config SCALED to the window boundaries m_lo / m_hi (the event's
        pt scaled by m_lo/m_obs and m_hi/m_obs), and ``ρ_t`` is the per-node
        nominal ρ at the event mass — the same un-kick used by the density, so
        no bisection or fixed point is needed. Cost: 3·n_gh² un-kicks (event +
        two boundaries) + 2·n_gh² flow-CDF evals per event."""
        B = m_obs.shape[0]
        # --disable-smearing: single GH node (see _continuity_logp_gh_qop) so the
        # norm correction does pure scale and shapes stay [B, G²=1]-consistent.
        ng = self.n_gh_nodes if self.smearing_enabled else 1
        xi, logW = _gh_nodes(ng, m_obs.device, m_obs.dtype)
        G = xi.shape[0]
        G2 = G * G
        xi_p = xi.view(G, 1).expand(G, G).reshape(G2)
        xi_m = xi.view(1, G).expand(G, G).reshape(G2)
        logW2 = (logW.view(G, 1) + logW.view(1, G)).reshape(G2)
        eps = torch.stack([xi_p, xi_m], dim=-1).view(1, G2, 2)
        mo = m_obs.unsqueeze(1)                                         # [B,1]
        pto = pt_obs.unsqueeze(1)                                       # [B,1,2]
        etao = eta_pm.unsqueeze(1)
        phio = phi_pm.unsqueeze(1)
        qo = q_pm.unsqueeze(1)
        bpo = b_pm.unsqueeze(1)
        # The observed config is self-consistent (mll = _event_mll(pt_obs)), so
        # _event_mll(pt_obs) = m_obs; scale pt_obs to the window boundaries
        # directly (pt∝m at fixed observed ρ). No rescaling to m_obs needed.

        # Per-node nominal conditioning at the event mass (for the CDF). The
        # window integral is over the mass at FIXED event conditioning, so the
        # event-mass un-kicked config gives the right ρ / event-level vector.
        _, pt_truth_evt = self._gh_qop_unsmear(pto, etao, phio, qo, bpo, eps)
        # Nominal masses at the two window boundaries (observed config scaled to
        # the per-event lower edge m_min(c) / m_hi along pt∝m, then un-kicked).
        m_lo_pe = self._fit_cut_m_min(
            m_obs, pt_obs, eta_pm, phi_pm).unsqueeze(1)  # [B,1]
        # muon-mass-EXACT pt scale to the window boundaries (observed mass at the
        # boundary config = m_edge exactly), not the massless pt∝m approximation.
        lam_lo = self._pt_lambda_to_mass(pto, etao, phio, mo, m_lo_pe)
        lam_hi = self._pt_lambda_to_mass(pto, etao, phio, mo, self._m_hi_f)
        m_t_lo, _ = self._gh_qop_unsmear(
            pto * lam_lo.unsqueeze(-1), etao, phio, qo, bpo, eps)
        m_t_hi, _ = self._gh_qop_unsmear(
            pto * lam_hi.unsqueeze(-1), etao, phio, qo, bpo, eps)
        mk_g = mk.unsqueeze(1).expand(B, G2, mk.shape[-1]).clone()
        if self.scale_enabled or self.smearing_enabled:
            mk_g = self._node_cond(mk_g, pt_truth_evt, etao, phio, qo)
        mkf = mk_g.reshape(B * G2, -1)
        # Per-node STABLE log window mass + GH-weighted logsumexp (see
        # _norm_correction_log_Z_gh) — replaces the cancellation-prone
        # Σ W·(F_hi.exp() − F_lo.exp()).
        log_nw = self._flow_log_window_Z_chunked(
            m_t_lo.reshape(-1), m_t_hi.reshape(-1), mkf).reshape(B, G2)
        return torch.logsumexp(logW2.view(1, G2) + log_nw, dim=1)       # [B]

    def data_nll_continuity(
        self,
        mll: torch.Tensor,
        pt_pm: torch.Tensor,
        eta_pm: torch.Tensor,
        phi_pm: torch.Tensor,
        q_pm: torch.Tensor,
        b_pm: torch.Tensor,
        muon_kin_std: torch.Tensor,
        is_data_mask: torch.Tensor,
        eps: float = 1e-30,
        n_iter: int = 2,
    ) -> torch.Tensor:
        """Per-event data NLL (``[B]``, unweighted) for stage 2.

        Signal density from the #2 direct evaluation (``_continuity_logp``):
        the frozen nominal flow forward-folded (advection + smear) by evaluating
        it at the source pre-images with the change-of-variables Jacobian — no
        flow derivatives, normalised by construction. Mixed with a degree-1
        Bernstein background via the MLP ``f(c)``. MC rows are ignored.
        """
        B = mll.shape[0]
        per = torch.zeros(B, dtype=mll.dtype, device=mll.device)
        data_idx = is_data_mask.nonzero(as_tuple=True)[0]
        if data_idx.numel() == 0:
            return per
        m = mll[data_idx]
        mk = muon_kin_std[data_idx]
        pt, eta, phi, q, b = (pt_pm[data_idx], eta_pm[data_idx], phi_pm[data_idx],
                              q_pm[data_idx], b_pm[data_idx])
        log_ps = self._continuity_logp(m, mk, pt, eta, phi, q, b, n_iter=n_iter)
        if self.norm_correction != "none":
            # Truncated-likelihood correction: the transformed-flow density is
            # naturally normalised on T(window), not on the window itself, so
            # log p_θ(x) carries a -log Z(θ;c) bias. Subtract per-event log Z
            # to restore the correct normalisation over the observation window.
            log_Z = self._norm_correction_log_Z(
                m, mk, pt, eta, phi, q, b, n_iter=n_iter)
            log_ps = log_ps - log_Z
        if self.background_enabled:
            # Normalise the background over the SAME per-event window as the
            # signal: [max(m_lo, m_min(c)), m_hi]. With no fit-time cuts this is
            # the scalar m_lo (exact legacy path); else a per-event [B] edge.
            bkg_m_lo = (self._fit_cut_m_min(m, pt, eta, phi)
                        if (self._fit_ptll_min or self._fit_pt_lead_min
                            or self._fit_pt_both_min) else self._m_lo_f)
            if self.bkg_model == "exp":
                # fractions + conditioning-dependent slope in one MLP pass
                f, s_bkg = self.mlp.forward_with_slope(mk)
                p_bkg = f[:, 0] * exp_bkg_density(
                    m, bkg_m_lo, self._m_hi_f, s_bkg)
            else:
                f = self.f_data(mk)                   # [n_d, n_bkg+1], signal LAST
                basis = bernstein_basis_n(
                    m, bkg_m_lo, self._m_hi_f, self.bkg_degree)
                p_bkg = (f[:, :-1] * basis).sum(-1)
            p_mix = p_bkg + f[:, -1] * log_ps.exp()
            per_data = -torch.log(p_mix.clamp_min(eps))
        else:
            # Background disabled (validation closure mode): the data branch
            # is pure signal — NLL = -log p_signal directly, no MLP / no
            # Bernstein. This removes the f_bkg ↔ smear degeneracy where the
            # MLP would otherwise absorb forward-|η| tails the signal can't
            # broaden into.
            per_data = -log_ps
        per = per.index_put((data_idx,), per_data, accumulate=False)
        return per
