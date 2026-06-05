"""Double-sided Crystal Ball conditional density (parametric stage-1 template).

The classic J/ψ line shape — Gaussian core with power-law tails on both sides —
with all six parameters (μ, σ, α_L, n_L, α_R, n_R) CONDITIONAL on the per-event
kinematics c through a single MLP. A fully ANALYTIC alternative to the flow
templates:

  * WINDOW-NORMALISED BY CONSTRUCTION: the density is divided by the analytic
    window integral I_W(c), so Z ≡ 1 per event — no out-of-window gauge
    freedom (the gf failure mode), nothing for --flow-gauge-penalty to do.
  * Analytic CDF (Φ for the core, closed-form power-law antiderivatives for
    the tails) → the window-Z machinery (un-kicked boundaries, decompose,
    display normalisation) is exact and perfectly conditioned.
  * The unnormalised shape extends smoothly to all of ℝ (power-law tails), so
    the operator's just-outside-window probes (un-kick, smear, boundary
    preimages) evaluate the REAL parametric tail — no extrapolation guesswork,
    no matched-tail machinery.
  * Pure closed-form elementwise ops → torch.compile(fullgraph)-friendly,
    vmap-friendly (Fisher), and cross-device deterministic at the ulp level.

Honest limitation: the DSCB is C¹ at the two core/tail junctions t = −α_L and
t = +α_R — the density and SCORE (∂_m log p) are continuous, but the curvature
∂²_m log p jumps there (cf. the C∞ flows). The stage-2 smear d² term sees two
isolated jump points; the scale-only fit is unaffected.

Parameter maps (raw MLP outputs → physical, all smooth, floored):
    μ   = raw                       (standardised-mass units; init 0 = peak)
    σ   = 0.05 + softplus(raw)      (core width, std units)
    α_X = 0.10 + softplus(raw)      (junction positions, in σ units)
    n_X = 1.01 + softplus(raw)      (tail powers; >1 keeps the tail
                                     antiderivatives finite/regular)
The final MLP layer is zero-initialised with biases at (μ, σ, α, n) =
(0, 0.9, 1.4, 3.0): the standardised mass is ~unit-variance with a narrower
core, so the init is already a credible J/ψ shape.

Exposes the interface the mass-fit model needs:
    forward(x_std, c) -> log p₀(x_std | c)   (window-normalised log density,
                                              evaluable on all of ℝ)
    log_cdf(x_std, c) -> log F₀(x_std | c)   (analytic cumulative with
                                              F₀(b) − F₀(a) = 1 exactly)
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

_SQRT2 = math.sqrt(2.0)
_SQRT2PI = math.sqrt(2.0 * math.pi)


def _phi_cdf(t: torch.Tensor) -> torch.Tensor:
    """Standard-normal CDF."""
    return 0.5 * (1.0 + torch.erf(t / _SQRT2))


def _softplus_inv(y: float) -> float:
    return math.log(math.expm1(y))


class DCBDensity(nn.Module):
    """Conditional double-sided Crystal Ball, window-normalised on [a, b]."""

    def __init__(
        self,
        n_cond: int,
        a: float,
        b: float,
        hidden_features: int = 128,
        n_layers: int = 3,
        activation: type[nn.Module] = nn.GELU,
        sigma_floor: float = 0.05,
        alpha_floor: float = 0.10,
        n_floor: float = 1.01,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        self._af = float(a)
        self._bf = float(b)
        self.sigma_floor = float(sigma_floor)
        self.alpha_floor = float(alpha_floor)
        self.n_floor = float(n_floor)

        seq: list[nn.Module] = []
        d = int(n_cond)
        for _ in range(int(n_layers)):
            seq += [nn.Linear(d, hidden_features), activation()]
            d = hidden_features
        final = nn.Linear(d, 6)
        with torch.no_grad():
            final.weight.zero_()
            # (μ, σ, α_L, n_L, α_R, n_R) init = (0, 0.9, 1.4, 3.0, 1.4, 3.0)
            final.bias[0] = 0.0
            final.bias[1] = _softplus_inv(0.9 - self.sigma_floor)
            final.bias[2] = _softplus_inv(1.4 - self.alpha_floor)
            final.bias[3] = _softplus_inv(3.0 - self.n_floor)
            final.bias[4] = _softplus_inv(1.4 - self.alpha_floor)
            final.bias[5] = _softplus_inv(3.0 - self.n_floor)
        seq.append(final)
        self.net = nn.Sequential(*seq)

    # ---- conditional parameters ---------------------------------------------
    def _params(self, c: torch.Tensor):
        """Raw MLP outputs → (μ, σ, α_L, n_L, α_R, n_R), each [B]."""
        h = self.net(c)
        mu = h[:, 0]
        sig = self.sigma_floor + F.softplus(h[:, 1])
        aL = self.alpha_floor + F.softplus(h[:, 2])
        nL = self.n_floor + F.softplus(h[:, 3])
        aR = self.alpha_floor + F.softplus(h[:, 4])
        nR = self.n_floor + F.softplus(h[:, 5])
        return mu, sig, aL, nL, aR, nR

    # ---- unnormalised shape in t = (x − μ)/σ units (peak value 1) -----------
    @staticmethod
    def _log_f_t(t, aL, nL, aR, nR):
        """log f(t): Gaussian core on [−α_L, α_R], CB power-law tails outside.
        The tail arguments are clamped before log so the UNSELECTED branch of
        torch.where never produces NaN (where's backward then stays clean)."""
        BL = nL / aL - aL
        BR = nR / aR - aR
        log_core = -0.5 * t * t
        log_L = (nL * torch.log(nL / aL) - 0.5 * aL * aL
                 - nL * torch.log((BL - t).clamp_min(1e-12)))
        log_R = (nR * torch.log(nR / aR) - 0.5 * aR * aR
                 - nR * torch.log((BR + t).clamp_min(1e-12)))
        return torch.where(t < -aL, log_L, torch.where(t > aR, log_R, log_core))

    @staticmethod
    def _F_t(t, aL, nL, aR, nR):
        """Unnormalised CDF ∫_{−∞}^t f(t')dt' in t units. Piecewise analytic:
        power-law antiderivatives for the tails (n > 1 keeps them finite) and
        √(2π)·Φ for the core. Continuous by construction.

        The tail pieces A·(B∓t)^{1−n}/(n−1) are evaluated in LOG space —
        A = (n/α)^n e^{−α²/2} overflows plain arithmetic already at n ~ 100 —
        with the exponent clamped at +30: the SELECTED region's exponent is
        bounded by ~log(n/α) − α²/2 − log(n−1) ≲ 15, so the clamp only tames
        the UNSELECTED ``where`` branch (whose clamped log argument would
        otherwise drive exp → inf and poison the backward with 0·inf)."""
        BL = nL / aL - aL
        BR = nR / aR - aR
        # total tail masses: e^{−α²/2}·n/(α(n−1)) — always moderate
        IL_tot = torch.exp(-0.5 * aL * aL) * nL / (aL * (nL - 1.0))
        IR_tot = torch.exp(-0.5 * aR * aR) * nR / (aR * (nR - 1.0))
        # left piece: A_L (B_L − t)^{1−n_L}/(n_L − 1), valid t ≤ −α_L
        expo_L = (nL * torch.log(nL / aL) - 0.5 * aL * aL
                  + (1.0 - nL) * torch.log((BL - t).clamp_min(1e-12))
                  - torch.log(nL - 1.0)).clamp(max=30.0)
        F_L = expo_L.exp()
        # core piece: I_L + √(2π)(Φ(t) − Φ(−α_L)), valid −α_L ≤ t ≤ α_R
        F_C = IL_tot + _SQRT2PI * (_phi_cdf(t) - _phi_cdf(-aL))
        core_tot = _SQRT2PI * (_phi_cdf(aR) - _phi_cdf(-aL))
        # right piece: everything below α_R + the partial right-tail integral
        expo_R = (nR * torch.log(nR / aR) - 0.5 * aR * aR
                  + (1.0 - nR) * torch.log((BR + t).clamp_min(1e-12))
                  - torch.log(nR - 1.0)).clamp(max=30.0)
        F_R = IL_tot + core_tot + (IR_tot - expo_R.exp())
        return torch.where(t < -aL, F_L, torch.where(t > aR, F_R, F_C))

    def _t_and_params(self, x, c):
        mu, sig, aL, nL, aR, nR = self._params(c)
        t = (x.reshape(-1) - mu) / sig
        return t, mu, sig, aL, nL, aR, nR

    def _log_window(self, mu, sig, aL, nL, aR, nR):
        """log I_W: the window integral of f in t units (σ-free; the σ factor
        cancels between density and CDF normalisation)."""
        t_a = (self._af - mu) / sig
        t_b = (self._bf - mu) / sig
        IW = (self._F_t(t_b, aL, nL, aR, nR)
              - self._F_t(t_a, aL, nL, aR, nR)).clamp_min(1e-30)
        return IW.log()

    # ---- public API (flow drop-in) ------------------------------------------
    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B] — WINDOW-NORMALISED (∫_a^b p dx ≡ 1):
        p = f(t)/(σ · I_W). Evaluable on all of ℝ (the operator's just-outside
        probes get the real parametric tail). ``x`` may be [B] or [B,1]."""
        t, mu, sig, aL, nL, aR, nR = self._t_and_params(x, c)
        return (self._log_f_t(t, aL, nL, aR, nR) - sig.log()
                - self._log_window(mu, sig, aL, nL, aR, nR))

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B], with F₀ = (cumulative of f)/I_W — so
        F₀(b) − F₀(a) = 1 EXACTLY and the window-Z differences at un-kicked
        boundaries are exact (the σ factors cancel in the ratio)."""
        t, mu, sig, aL, nL, aR, nR = self._t_and_params(x, c)
        logF = self._F_t(t, aL, nL, aR, nR).clamp_min(1e-300).log()
        return logF - self._log_window(mu, sig, aL, nL, aR, nR)


class EGEDensity(nn.Module):
    """Gaussian core + C¹-MATCHED exponential tails (ExpGaussExp), conditional
    on c through an MLP and window-normalised on [a, b].

    Per tail the exponential has two dof (amplitude, slope) and BOTH are spent
    on the matching: C⁰ fixes the amplitude and C¹ forces the slope λ = α (the
    junction position in σ units) — "matched derivatives as much as possible"
    (C² is impossible: a log-linear tail cannot match the core's log-curvature
    −1). So the shape has only FOUR conditional parameters (μ, σ, α_L, α_R):

        f(t) = exp(α_L²/2 + α_L t)   t < −α_L      (log-linear, slope +α_L)
             = exp(−t²/2)            −α_L ≤ t ≤ α_R
             = exp(α_R²/2 − α_R t)   t > α_R

    vs the DCB this trades the power-law tails for exponentials: no n
    parameters, no n>1 integrability constraint, an even simpler analytic CDF
    (pure exp + Φ), and the slope matching is exact by construction. Same
    window-normalised design (Z ≡ 1, no gauge), real tails on all of ℝ for
    the operator probes, fullgraph-compileable, and the same honest C¹
    limitation (score continuous EXACTLY; curvature jumps at two points)."""

    def __init__(
        self,
        n_cond: int,
        a: float,
        b: float,
        hidden_features: int = 128,
        n_layers: int = 3,
        activation: type[nn.Module] = nn.GELU,
        sigma_floor: float = 0.05,
        alpha_floor: float = 0.10,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        self._af = float(a)
        self._bf = float(b)
        self.sigma_floor = float(sigma_floor)
        self.alpha_floor = float(alpha_floor)

        seq: list[nn.Module] = []
        d = int(n_cond)
        for _ in range(int(n_layers)):
            seq += [nn.Linear(d, hidden_features), activation()]
            d = hidden_features
        final = nn.Linear(d, 4)
        with torch.no_grad():
            final.weight.zero_()
            # (μ, σ, α_L, α_R) init = (0, 0.9, 1.4, 1.4)
            final.bias[0] = 0.0
            final.bias[1] = _softplus_inv(0.9 - self.sigma_floor)
            final.bias[2] = _softplus_inv(1.4 - self.alpha_floor)
            final.bias[3] = _softplus_inv(1.4 - self.alpha_floor)
        seq.append(final)
        self.net = nn.Sequential(*seq)

    def _params(self, c: torch.Tensor):
        h = self.net(c)
        mu = h[:, 0]
        sig = self.sigma_floor + F.softplus(h[:, 1])
        aL = self.alpha_floor + F.softplus(h[:, 2])
        aR = self.alpha_floor + F.softplus(h[:, 3])
        return mu, sig, aL, aR

    @staticmethod
    def _log_f_t(t, aL, aR):
        """log f(t): Gaussian core, C¹-matched exponential tails. No clamps
        needed — the tail logs are linear in t (finite everywhere), and the
        unselected ``where`` branches never pass through an exp."""
        log_core = -0.5 * t * t
        log_L = 0.5 * aL * aL + aL * t
        log_R = 0.5 * aR * aR - aR * t
        return torch.where(t < -aL, log_L, torch.where(t > aR, log_R, log_core))

    @staticmethod
    def _F_t(t, aL, aR):
        """Unnormalised CDF in t units: exp antiderivatives for the tails, Φ
        for the core. The tail exponents are clamped at 0 so the UNSELECTED
        ``where`` branch's exp stays finite (the selected region's exponent is
        ≤ −α²/2 < 0, so the clamp never bites there)."""
        IL_tot = torch.exp(-0.5 * aL * aL) / aL
        IR_tot = torch.exp(-0.5 * aR * aR) / aR
        F_L = (0.5 * aL * aL + aL * t).clamp(max=0.0).exp() / aL
        F_C = IL_tot + _SQRT2PI * (_phi_cdf(t) - _phi_cdf(-aL))
        core_tot = _SQRT2PI * (_phi_cdf(aR) - _phi_cdf(-aL))
        F_R = (IL_tot + core_tot
               + (IR_tot - (0.5 * aR * aR - aR * t).clamp(max=0.0).exp() / aR))
        return torch.where(t < -aL, F_L, torch.where(t > aR, F_R, F_C))

    def _t_and_params(self, x, c):
        mu, sig, aL, aR = self._params(c)
        t = (x.reshape(-1) - mu) / sig
        return t, mu, sig, aL, aR

    def _log_window(self, mu, sig, aL, aR):
        t_a = (self._af - mu) / sig
        t_b = (self._bf - mu) / sig
        IW = (self._F_t(t_b, aL, aR) - self._F_t(t_a, aL, aR)).clamp_min(1e-30)
        return IW.log()

    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B] — window-normalised; evaluable on all of ℝ."""
        t, mu, sig, aL, aR = self._t_and_params(x, c)
        return (self._log_f_t(t, aL, aR) - sig.log()
                - self._log_window(mu, sig, aL, aR))

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B], with F₀(b) − F₀(a) = 1 exactly."""
        t, mu, sig, aL, aR = self._t_and_params(x, c)
        logF = self._F_t(t, aL, aR).clamp_min(1e-300).log()
        return logF - self._log_window(mu, sig, aL, aR)


# ---------------------------------------------------------------------------
# Self-test: normalisation, CDF↔density consistency, C¹ junctions, tails,
# gradient flow, extreme-parameter stability, quick MLE recovery.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    a, b = -3.7, 3.9
    net = DCBDensity(7, a, b, hidden_features=64, n_layers=2).double()

    B = 4
    c = torch.randn(B, 7, dtype=torch.float64) * 0.5
    # 1) zero-init (biases only): proper window normalisation + CDF exactness
    xg = torch.linspace(a, b, 40001, dtype=torch.float64)
    for e in range(B):
        ce = c[e:e + 1].expand(xg.shape[0], -1)
        I = torch.trapezoid(net(xg, ce).exp(), xg).item()
        Fa = net.log_cdf(torch.tensor([a], dtype=torch.float64), c[e:e + 1]).exp().item()
        Fb = net.log_cdf(torch.tensor([b], dtype=torch.float64), c[e:e + 1]).exp().item()
        assert abs(I - 1.0) < 1e-6 and abs(Fb - Fa - 1.0) < 1e-12, (I, Fb - Fa)
    print("zero-init: ∫_window p = 1 (trapz 1e-6); F(b)−F(a) = 1 exactly")

    # 2) perturbed MLP: dF/dx ≡ p including OUTSIDE the window (real tails)
    with torch.no_grad():
        for p_ in net.parameters():
            p_.add_(0.3 * torch.randn_like(p_))
    xt = torch.tensor([a - 0.5, a - 0.01, a + 0.01, -1.0, 0.0, 1.0,
                       b - 0.01, b + 0.01, b + 0.5], dtype=torch.float64)
    ce = c[:1].expand(xt.shape[0], -1)
    h = 1e-6
    dF = (net.log_cdf(xt + h, ce).exp() - net.log_cdf(xt - h, ce).exp()) / (2 * h)
    p = net(xt, ce).exp()
    rel = float(((dF - p).abs() / p.clamp_min(1e-12)).max())
    assert rel < 1e-6, rel
    print(f"perturbed: max rel diff dF/dx vs p = {rel:.1e} (incl. out-of-window)")
    # window normalisation survives the perturbation
    for e in range(B):
        ce_ = c[e:e + 1].expand(xg.shape[0], -1)
        I = torch.trapezoid(net(xg, ce_).exp(), xg).item()
        assert abs(I - 1.0) < 1e-6, I
    print("perturbed: ∫_window p = 1 for all conditionings")

    # 3) C¹ at the junctions: density and score continuous (curvature jumps)
    mu, sig, aL, nL, aR, nR = net._params(c[:1])
    for tj in (-(aL[0]), aR[0]):
        xj = (mu[0] + sig[0] * tj).item()
        eps = 1e-7
        lp = net(torch.tensor([xj - eps, xj + eps], dtype=torch.float64),
                 c[:1].expand(2, -1))
        score = (net(torch.tensor([xj - 2 * eps, xj - eps], dtype=torch.float64),
                     c[:1].expand(2, -1)).diff() / eps,
                 net(torch.tensor([xj + eps, xj + 2 * eps], dtype=torch.float64),
                     c[:1].expand(2, -1)).diff() / eps)
        assert abs(float(lp.diff())) < 1e-5
        assert abs(float(score[0] - score[1])) < 1e-3
    print("C¹ at both core/tail junctions (density + score continuous)")

    # 4) gradients reach all parameters through forward and log_cdf
    for fn in (lambda: net(xt, ce).sum(), lambda: net.log_cdf(xt, ce).sum()):
        net.zero_grad()
        fn().backward()
        assert all(p_.grad is not None and torch.isfinite(p_.grad).all()
                   for p_ in net.parameters())
    print("gradients finite through forward and log_cdf")

    # 5) extreme conditioners: finite everywhere on a wide grid
    f2 = DCBDensity(7, a, b, hidden_features=64, n_layers=2).double()
    with torch.no_grad():
        for p_ in f2.parameters():
            p_.normal_(0, 3.0)
    cc = torch.randn(64, 7, dtype=torch.float64)
    xs = torch.linspace(-8, 8, 200, dtype=torch.float64)
    ok = all(torch.isfinite(f2(xs, cc[e:e + 1].expand(200, -1))).all()
             and torch.isfinite(f2.log_cdf(xs, cc[e:e + 1].expand(200, -1))).all()
             for e in range(64))
    assert ok
    print("overflow stress (64 extreme conditioners): all finite")

    # 6) quick MLE recovery of a known shape (truncated DSCB samples via
    #    rejection from the init density)
    torch.manual_seed(1)
    target = DCBDensity(1, a, b, hidden_features=8, n_layers=1).double()
    with torch.no_grad():                      # target: fixed, c-independent
        for p_ in target.parameters():
            p_.zero_()
        target.net[-1].bias.copy_(torch.tensor(
            [0.3, _softplus_inv(0.55), _softplus_inv(1.2 - 0.10),
             _softplus_inv(2.5 - 1.01), _softplus_inv(1.9 - 0.10),
             _softplus_inv(4.0 - 1.01)], dtype=torch.float64))
    xs = torch.linspace(a, b, 200001, dtype=torch.float64)
    w = target(xs, torch.zeros(xs.shape[0], 1, dtype=torch.float64)).exp()
    idx = torch.multinomial(w, 200000, replacement=True)
    sample = xs[idx] + (xs[1] - xs[0]) * (torch.rand(200000, dtype=torch.float64) - 0.5)
    fit = DCBDensity(1, a, b, hidden_features=8, n_layers=1).double()
    opt = torch.optim.Adam(fit.parameters(), lr=2e-2)
    cz = torch.zeros(200000, 1, dtype=torch.float64)
    for it in range(400):
        opt.zero_grad()
        (-fit(sample, cz).mean()).backward()
        opt.step()
    mu_t, sig_t, *_ = target._params(torch.zeros(1, 1, dtype=torch.float64))
    mu_f, sig_f, *_ = fit._params(torch.zeros(1, 1, dtype=torch.float64))
    print(f"MLE recovery: μ {float(mu_f):+.4f} (truth {float(mu_t):+.4f})  "
          f"σ {float(sig_f):.4f} (truth {float(sig_t):.4f})")
    assert abs(float(mu_f - mu_t)) < 0.01 and abs(float(sig_f - sig_t)) < 0.02
    print("ALL DCB SELF-TESTS PASSED")

    # ---- EGE (Gaussian core + C¹-matched exponential tails) ----------------
    torch.manual_seed(2)
    ege = EGEDensity(7, a, b, hidden_features=64, n_layers=2).double()
    with torch.no_grad():
        for p_ in ege.parameters():
            p_.add_(0.3 * torch.randn_like(p_))
    c = torch.randn(4, 7, dtype=torch.float64) * 0.5
    for e in range(4):
        ce_ = c[e:e + 1].expand(xg.shape[0], -1)
        I = torch.trapezoid(ege(xg, ce_).exp(), xg).item()
        Fa = ege.log_cdf(torch.tensor([a], dtype=torch.float64), c[e:e + 1]).exp().item()
        Fb = ege.log_cdf(torch.tensor([b], dtype=torch.float64), c[e:e + 1]).exp().item()
        assert abs(I - 1.0) < 1e-6 and abs(Fb - Fa - 1.0) < 1e-12, (I, Fb - Fa)
    print("EGE: ∫_window p = 1; F(b)−F(a) = 1 exactly (perturbed MLP)")
    ce = c[:1].expand(xt.shape[0], -1)
    dF = (ege.log_cdf(xt + h, ce).exp() - ege.log_cdf(xt - h, ce).exp()) / (2 * h)
    p = ege(xt, ce).exp()
    rel = float(((dF - p).abs() / p.clamp_min(1e-12)).max())
    assert rel < 1e-6, rel
    print(f"EGE: max rel diff dF/dx vs p = {rel:.1e} (incl. out-of-window)")
    # C¹ junctions: the slope matching is ALGEBRAIC for EGE — test tightly
    mu, sig, aL, aR = ege._params(c[:1])
    for tj, slope in ((-aL[0], aL[0]), (aR[0], -aR[0])):
        xj = (mu[0] + sig[0] * tj).item()
        eps = 1e-6
        sc_in = float(ege(torch.tensor([xj - 2 * eps, xj - eps], dtype=torch.float64),
                          c[:1].expand(2, -1)).diff()) / eps
        sc_out = float(ege(torch.tensor([xj + eps, xj + 2 * eps], dtype=torch.float64),
                           c[:1].expand(2, -1)).diff()) / eps
        assert abs(sc_in - sc_out) < 1e-4 * max(abs(sc_in), 1.0), (sc_in, sc_out)
    print("EGE: score continuous at both junctions (C¹ algebraic)")
    for fn in (lambda: ege(xt, ce).sum(), lambda: ege.log_cdf(xt, ce).sum()):
        ege.zero_grad(); fn().backward()
        assert all(p_.grad is not None and torch.isfinite(p_.grad).all()
                   for p_ in ege.parameters())
    e2 = EGEDensity(7, a, b, hidden_features=64, n_layers=2).double()
    with torch.no_grad():
        for p_ in e2.parameters():
            p_.normal_(0, 3.0)
    xs_st = torch.linspace(-8, 8, 200, dtype=torch.float64)
    ok = all(torch.isfinite(e2(xs_st, cc[e:e + 1].expand(200, -1))).all()
             and torch.isfinite(e2.log_cdf(xs_st, cc[e:e + 1].expand(200, -1))).all()
             for e in range(64))
    assert ok
    print("EGE: gradients finite; overflow stress (64 extreme conditioners) OK")
    # quick MLE recovery
    torch.manual_seed(3)
    tgt = EGEDensity(1, a, b, hidden_features=8, n_layers=1).double()
    with torch.no_grad():
        for p_ in tgt.parameters():
            p_.zero_()
        tgt.net[-1].bias.copy_(torch.tensor(
            [0.25, _softplus_inv(0.6 - 0.05), _softplus_inv(1.1 - 0.10),
             _softplus_inv(2.0 - 0.10)], dtype=torch.float64))
    w = tgt(xs2 := torch.linspace(a, b, 200001, dtype=torch.float64),
            torch.zeros(200001, 1, dtype=torch.float64)).exp()
    idx = torch.multinomial(w, 200000, replacement=True)
    sample = xs2[idx] + (xs2[1] - xs2[0]) * (torch.rand(200000, dtype=torch.float64) - 0.5)
    fit2 = EGEDensity(1, a, b, hidden_features=8, n_layers=1).double()
    opt = torch.optim.Adam(fit2.parameters(), lr=2e-2)
    cz = torch.zeros(200000, 1, dtype=torch.float64)
    for it in range(400):
        opt.zero_grad(); (-fit2(sample, cz).mean()).backward(); opt.step()
    mu_t, sig_t, *_ = tgt._params(torch.zeros(1, 1, dtype=torch.float64))
    mu_f, sig_f, *_ = fit2._params(torch.zeros(1, 1, dtype=torch.float64))
    print(f"EGE MLE recovery: μ {float(mu_f):+.4f} (truth {float(mu_t):+.4f})  "
          f"σ {float(sig_f):.4f} (truth {float(sig_t):.4f})")
    assert abs(float(mu_f - mu_t)) < 0.01 and abs(float(sig_f - sig_t)) < 0.02
    print("ALL EGE SELF-TESTS PASSED")


if __name__ == "__main__":
    _selftest()
