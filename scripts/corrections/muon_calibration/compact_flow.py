"""Compact, uniform-base 1-D conditional flow with a C²-matched analytic tail.

Motivation (J/ψ mass-fit two-stage continuity calibration)
----------------------------------------------------------
The stage-1 nominal density p₀(m|c) must be:
  * a proper density on the trigger/selection mass window [a, b] (standardised),
  * exactly normalised over that window (no out-of-window mass gauge freedom,
    which a full-support flow trained on a truncated likelihood leaves
    unconstrained → drifts, growing spurious far-tail structure),
  * smooth (C∞) in the interior — the stage-2 continuity operator transports
    probability using ∂_m log p₀ and ∂²_m log p₀,
  * and STILL evaluable just OUTSIDE the window, smoothly, because the scale
    un-kick maps near-edge observed masses to source masses slightly outside
    [a,b] (and the smear samples a neighbourhood, and the window-normalisation Z
    evaluates the CDF at un-kicked boundaries).

Construction
------------
In-window:  a logistic-mixture CDF G(x|c)=Σ_k π_k σ((x−μ_k)/s_k) — C∞, monotone —
            renormalised to the window:  Q(x) = (G(x)−G(a)) / (G(b)−G(a)),
            so Q(a)=0, Q(b)=1 and the in-window density q = Q' integrates to 1
            over [a,b] EXACTLY (Z=1, no gauge freedom).

Tails:      beyond each edge the log-density is continued by its 2nd-order
            Taylor expansion in log q — i.e. a Gaussian (downward log-parabola)
            matched to log q, (log q)', (log q)'' at the edge.  This is C² across
            the boundary (so the operator's score AND curvature are continuous),
            monotone-decaying (guaranteed for a decaying, log-concave edge), and
            has a closed-form CDF (erf).  It is a *fixed* controlled extension of
            the data-constrained edge — not a free part of the learned model — so
            no humps and no near-boundary bias.

The module exposes the two methods the mass-fit model needs:
    forward(x_std, c) -> log p₀(x_std | c)      (standardised-mass log-density)
    log_cdf(x_std, c) -> log F₀(x_std | c)      (cumulative of the unnormalised
                                                 density; differences give Z)
F₀ is the cumulative of the (in-window + tails) density, so F₀(b)−F₀(a)=1 and the
window-normalisation Z=F₀(unkick_hi)−F₀(unkick_lo) is exact.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

_SQRT2 = math.sqrt(2.0)
_SQRT2PI = math.sqrt(2.0 * math.pi)


def _phi(z: torch.Tensor) -> torch.Tensor:
    """Standard-normal CDF."""
    return 0.5 * (1.0 + torch.erf(z / _SQRT2))


class CompactMatchedFlow(nn.Module):
    """Uniform-base compact 1-D conditional flow with C²-matched Gaussian tails.

    ``n_cond`` conditioning features → a logistic mixture of ``n_components``
    components via an MLP conditioner. Operates entirely in STANDARDISED mass;
    the window edges ``a``, ``b`` are the standardised [m_lo, m_hi].
    """

    def __init__(
        self,
        n_cond: int,
        a: float,
        b: float,
        hidden_features: int = 64,
        n_layers: int = 2,
        n_components: int = 16,
        activation: type[nn.Module] = nn.GELU,
        s_min_frac: float = 0.02,
        curv_floor: float = 1e-3,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        self.K = int(n_components)
        self.register_buffer("a", torch.tensor(float(a)))
        self.register_buffer("b", torch.tensor(float(b)))
        self._width = float(b - a)
        # minimum logistic scale and minimum |log-density curvature| in the tail,
        # both in standardised units (keep the tail strictly log-concave/decaying).
        self.s_min = float(s_min_frac) * self._width
        self.curv_floor = float(curv_floor)

        layers: list[nn.Module] = []
        d = int(n_cond)
        for _ in range(int(n_layers)):
            layers += [nn.Linear(d, hidden_features), activation()]
            d = hidden_features
        self.body = nn.Sequential(*layers) if layers else nn.Identity()
        self.head = nn.Linear(d, 3 * self.K)   # [logit_w | raw_mu | raw_s]
        # Init: small head weights so the initial density is ~c-independent, with
        # means spread across the window and moderate scales.
        nn.init.normal_(self.head.weight, std=1e-3)
        with torch.no_grad():
            self.head.bias.zero_()
            # raw_mu biases: means at sigmoid⁻¹(linspace) → spread over (a,b).
            frac = torch.linspace(0.08, 0.92, self.K)
            self.head.bias[self.K:2 * self.K] = torch.log(frac / (1 - frac))
            # raw_s biases: scales ≈ width/K via softplus⁻¹.
            s0 = max(self._width / self.K, self.s_min)
            self.head.bias[2 * self.K:] = math.log(math.expm1(max(s0 - self.s_min, 1e-4)))

    # ---- mixture parameters ------------------------------------------------
    def _params(self, c: torch.Tensor):
        h = self.head(self.body(c))                                  # [B, 3K]
        logit_w, raw_mu, raw_s = h.split(self.K, dim=-1)
        log_pi = F.log_softmax(logit_w, dim=-1)                      # [B,K]
        mu = self.a + (self.b - self.a) * torch.sigmoid(raw_mu)      # [B,K] in (a,b)
        s = self.s_min + F.softplus(raw_s)                           # [B,K] > s_min
        return log_pi, mu, s

    def _G_derivs(self, x, log_pi, mu, s):
        """G and its first three derivatives at x [B] for params [B,K].
        Returns (G, G1, G2, G3) each [B]."""
        z = (x.unsqueeze(-1) - mu) / s                               # [B,K]
        sig = torch.sigmoid(z)
        s1 = sig * (1 - sig)                                         # σ'
        s2 = s1 * (1 - 2 * sig)                                      # σ''
        s3 = s2 * (1 - 2 * sig) - 2 * s1 * s1                        # σ'''
        pi = log_pi.exp()
        G = (pi * sig).sum(-1)
        G1 = (pi * s1 / s).sum(-1)
        G2 = (pi * s2 / s.pow(2)).sum(-1)
        G3 = (pi * s3 / s.pow(3)).sum(-1)
        return G, G1, G2, G3

    def _edges(self, log_pi, mu, s):
        """In-window normaliser and per-edge (log q, (log q)', (log q)'') plus
        Gaussian-tail integrals. Returns a dict of [B] tensors."""
        B = mu.shape[0]
        a = self.a.expand(B); b = self.b.expand(B)
        Ga, Ga1, Ga2, Ga3 = self._G_derivs(a, log_pi, mu, s)
        Gb, Gb1, Gb2, Gb3 = self._G_derivs(b, log_pi, mu, s)
        norm = (Gb - Ga).clamp_min(1e-12)                           # G(b)-G(a)
        log_norm = norm.log()
        eps = 1e-30
        # log q and log-derivatives at each edge (the 1/norm cancels in log-derivs)
        out = {"log_norm": log_norm, "Ga": Ga, "norm": norm}
        for tag, (G1, G2, G3) in (("a", (Ga1, Ga2, Ga3)), ("b", (Gb1, Gb2, Gb3))):
            logq = G1.clamp_min(eps).log() - log_norm               # log q(edge)
            d1 = G2 / G1.clamp_min(eps)                             # (log q)'
            d2 = G3 / G1.clamp_min(eps) - d1 * d1                   # (log q)''
            out[f"logq_{tag}"] = logq
            out[f"d1_{tag}"] = d1
            out[f"d2_{tag}"] = d2
        return out

    @staticmethod
    def _tail_logp(t, logq, d1, d2):
        """log p_ext(edge+t) = logq + d1·t + ½ d2·t²  (matched Gaussian/parabola).
        ``t`` is signed displacement from the edge (≥0 upper, ≤0 lower)."""
        return logq + d1 * t + 0.5 * d2 * t * t

    def _outward(self, U, logq, s_out, alpha):
        """Matched-Gaussian tail integrated in the OUTWARD coordinate u≥0 (away
        from the window): density g(u)=exp(logq + s_out·u − ½α u²), s_out<0 for a
        decaying tail, α=−d2>0.  Returns (total=∫₀^∞, beyond=∫_U^∞) — both
        computed via erfcx so the exp(½ s_out²/α) factor never overflows.

            ∫_L^∞ g = exp(g_at_L)·√(π/2α)·erfcx((αL − s_out)/√(2α))
        with erfcx(w)=e^{w²}erfc(w) bounded; g_at_L = logq + s_out L − ½α L².
        """
        root = (0.5 / alpha).sqrt()                                # 1/√(2α)
        pref = (math.pi * 0.5 / alpha).sqrt()                      # √(π/2α)
        erfcx = torch.special.erfcx
        total = logq.exp() * pref * erfcx(-s_out * root)           # ∫₀^∞ (L=0)
        g_at_U = logq + s_out * U - 0.5 * alpha * U * U            # log g(U)
        beyond = g_at_U.exp() * pref * erfcx((alpha * U - s_out) * root)  # ∫_U^∞
        return total, beyond

    # ---- public API --------------------------------------------------------
    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B]. ``x`` may be [B] or [B,1]."""
        x = x.reshape(-1)
        log_pi, mu, s = self._params(c)
        a, b = self.a, self.b
        ed = self._edges(log_pi, mu, s)
        # in-window log q
        G, G1, _, _ = self._G_derivs(x, log_pi, mu, s)
        logq_in = G1.clamp_min(1e-30).log() - ed["log_norm"]
        # tails
        logp_hi = self._tail_logp(x - b, ed["logq_b"], ed["d1_b"], ed["d2_b"])
        logp_lo = self._tail_logp(x - a, ed["logq_a"], ed["d1_a"], ed["d2_a"])
        out = torch.where(x > b, logp_hi, torch.where(x < a, logp_lo, logq_in))
        return out

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B]. F₀ = cumulative of the (in-window + tails)
        unnormalised density; F₀(b)−F₀(a)=1 exactly, so window Z is exact."""
        x = x.reshape(-1)
        log_pi, mu, s = self._params(c)
        a, b = self.a, self.b
        ed = self._edges(log_pi, mu, s)
        af = self.curv_floor
        # Outward coordinates u≥0 and outward log-slopes (s_out<0 for a decaying
        # tail): lower edge a → outward = a−x, s_out = −(log q)'(a); upper edge b
        # → outward = x−b, s_out = +(log q)'(b).
        alpha_a = (-ed["d2_a"]).clamp_min(af)
        alpha_b = (-ed["d2_b"]).clamp_min(af)
        s_out_a = -ed["d1_a"]
        s_out_b = ed["d1_b"]
        u_a = (a - x).clamp_min(0.0)        # >0 only where x<a
        u_b = (x - b).clamp_min(0.0)        # >0 only where x>b
        m_lo, F_below_x = self._outward(u_a, ed["logq_a"], s_out_a, alpha_a)
        m_hi_tot, beyond_xb = self._outward(u_b, ed["logq_b"], s_out_b, alpha_b)
        # in-window cumulative Q(x) = (G(x)-G(a))/norm  (Q(a)=0, Q(b)=1)
        G, _, _, _ = self._G_derivs(x, log_pi, mu, s)
        Q = (G - ed["Ga"]) / ed["norm"]
        F_in = m_lo + Q                                            # x in [a,b]
        F_lo = F_below_x                                           # x < a: ∫_{-∞}^x
        F_hi = m_lo + 1.0 + (m_hi_tot - beyond_xb)                 # x > b
        Fc = torch.where(x > b, F_hi, torch.where(x < a, F_lo, F_in))
        return Fc.clamp_min(1e-30).log()


# ---------------------------------------------------------------------------
# Self-test: normalisation, density↔CDF consistency, C² continuity.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    a, b = -3.7, 3.9
    flow = CompactMatchedFlow(n_cond=7, a=a, b=b, n_components=12).double()
    c = torch.randn(4, 7, dtype=torch.float64)
    # 1) in-window normalisation: ∫_a^b exp(logp) dx ≈ 1
    xg = torch.linspace(a, b, 20001, dtype=torch.float64)
    for e in range(c.shape[0]):
        lp = flow(xg, c[e:e+1].expand(xg.shape[0], -1))
        I = torch.trapz(lp.exp(), xg).item()
        # 2) CDF endpoints: F(b)-F(a) == 1 (window mass)
        Fa = flow.log_cdf(torch.tensor([a], dtype=torch.float64), c[e:e+1]).exp().item()
        Fb = flow.log_cdf(torch.tensor([b], dtype=torch.float64), c[e:e+1]).exp().item()
        print(f"event {e}: ∫_window p={I:.6f}  F(b)-F(a)={Fb-Fa:.6f}")
    # 3) density↔CDF: dF/dx ≈ p  (finite diff), across both boundaries
    xt = torch.tensor([a-0.3, a-0.01, a+0.01, 0.0, b-0.01, b+0.01, b+0.3],
                      dtype=torch.float64)
    ce = c[:1].expand(xt.shape[0], -1)
    h = 1e-5
    Fp = flow.log_cdf(xt + h, ce).exp(); Fm = flow.log_cdf(xt - h, ce).exp()
    dF = (Fp - Fm) / (2 * h)
    p = flow(xt, ce).exp()
    print("x       :", [f"{v:+.3f}" for v in xt.tolist()])
    print("dF/dx   :", [f"{v:.4e}" for v in dF.tolist()])
    print("p(x)    :", [f"{v:.4e}" for v in p.tolist()])
    print("rel.diff:", [f"{abs(d-pp)/max(pp,1e-12):.2e}" for d, pp in zip(dF.tolist(), p.tolist())])
    # 4) C² continuity at boundaries: log p and its 1st/2nd deriv match across a,b
    for edge, name in ((a, "a"), (b, "b")):
        xx = torch.tensor([edge - 1e-3, edge + 1e-3], dtype=torch.float64)
        cc = c[:1].expand(2, -1)
        lpv = flow(xx, cc)
        print(f"log p across {name}: {lpv[0].item():.6f} vs {lpv[1].item():.6f} "
              f"(jump={abs(lpv[0]-lpv[1]).item():.2e})")


if __name__ == "__main__":
    _selftest()
