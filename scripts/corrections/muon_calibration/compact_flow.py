"""Compact, COMPOSED uniform-base 1-D conditional flow with C²-matched tails.

Motivation (J/ψ mass-fit two-stage continuity calibration)
----------------------------------------------------------
The stage-1 nominal density p₀(m|c) must be:
  * a proper density on the trigger/selection mass window [a, b] (standardised),
  * EXACTLY normalised over that window (no out-of-window mass gauge freedom →
    no spurious far-tail structure, no wasted capacity outside the window),
  * smooth (C∞) in the interior — the stage-2 continuity operator transports
    probability using ∂_m log p₀ and ∂²_m log p₀,
  * and still evaluable just OUTSIDE the window, smoothly, because the scale
    un-kick maps near-edge masses to source masses slightly outside [a,b] (and
    the smear samples a neighbourhood, and the window-Z evaluates the CDF at
    un-kicked boundaries).

Construction
------------
In-window:  a COMPOSITION of ``n_transforms`` monotone bijections [0,1]→[0,1]
            (same depth idea as the Gaussianization flow, but on a BOUNDED domain
            with a uniform base, so it is intrinsically compact). Each layer is a
            logistic-mixture CDF renormalised to [0,1]:
                S_l(u) = (G_l(u) − G_l(0)) / (G_l(1) − G_l(0)),
                G_l(u) = Σ_k π_k σ((u−μ_k)/s_k),  μ_k∈(0,1), s_k>0.
            With x↦u₀=(x−a)/(b−a) and the uniform base, the composition
            F = S_L∘…∘S₁ : [0,1]→[0,1] IS the in-window CDF (F(a)=0, F(b)=1) and
            the density q = F' integrates to 1 over [a,b] EXACTLY. C∞ interior,
            expressive via depth (few components per layer, like gf).

Tails:      beyond each edge the log-density is continued by its 2nd-order Taylor
            in log q — a Gaussian (downward log-parabola) matched to log q,
            (log q)', (log q)'' at the edge (computed by AUTOGRAD through the
            composition). C² across the boundary, monotone-decaying, closed-form
            CDF (erf, via erfcx for overflow safety). A FIXED controlled extension
            of the data-constrained edge — no humps, no near-boundary bias.

Exposes the interface the mass-fit model needs:
    forward(x_std, c) -> log p₀(x_std | c)      (standardised-mass log-density)
    log_cdf(x_std, c) -> log F₀(x_std | c)      (cumulative of the unnormalised
                                                 in-window+tails density; F₀(b)−
                                                 F₀(a)=1, so window Z is exact)
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

_SQRT2 = math.sqrt(2.0)


def _phi(z: torch.Tensor) -> torch.Tensor:
    """Standard-normal CDF."""
    return 0.5 * (1.0 + torch.erf(z / _SQRT2))


class CompactMatchedFlow(nn.Module):
    """Composed uniform-base compact 1-D conditional flow with matched tails."""

    def __init__(
        self,
        n_cond: int,
        a: float,
        b: float,
        hidden_features: int = 64,
        n_layers: int = 2,
        n_components: int = 8,
        n_transforms: int = 5,
        activation: type[nn.Module] = nn.GELU,
        s_min_frac: float = 0.01,
        curv_floor: float = 1e-3,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        self.K = int(n_components)
        self.L = int(n_transforms)
        self.register_buffer("a", torch.tensor(float(a)))
        self.register_buffer("b", torch.tensor(float(b)))
        self._width = float(b - a)
        self.s_min = float(s_min_frac)          # min logistic scale in [0,1] units
        self.curv_floor = float(curv_floor)     # min |(log p)''| in std-mass units

        # One INDEPENDENT conditioner MLP per layer (c → [logit_w|raw_mu|raw_s]),
        # mirroring gf's per-transform hyper-network (so the conditioning capacity
        # and parameter count match gf), rather than a shared body + linear head.
        frac = torch.linspace(0.08, 0.92, self.K)
        mu_b = torch.log(frac / (1 - frac))                          # sigmoid⁻¹
        s0 = max(1.0 / self.K, self.s_min)
        s_b = math.log(math.expm1(max(s0 - self.s_min, 1e-4)))       # softplus⁻¹
        self.conditioners = nn.ModuleList()
        for _ in range(self.L):
            seq: list[nn.Module] = []
            d = int(n_cond)
            for _ in range(int(n_layers)):
                seq += [nn.Linear(d, hidden_features), activation()]
                d = hidden_features
            final = nn.Linear(d, 3 * self.K)
            nn.init.normal_(final.weight, std=1e-3)                  # ~c-independent init
            with torch.no_grad():
                final.bias.zero_()
                b = final.bias.view(3, self.K)                       # [logit_w|raw_mu|raw_s]
                b[1, :] = mu_b                                       # spread means in (0,1)
                b[2, :] = s_b                                        # moderate scales
            seq.append(final)
            self.conditioners.append(nn.Sequential(*seq))

    # ---- per-layer mixture parameters --------------------------------------
    def _layer_params(self, c: torch.Tensor):
        h = torch.stack([cond(c) for cond in self.conditioners], dim=1)  # [B,L,3K]
        h = h.view(-1, self.L, 3, self.K)                            # [B,L,3,K]
        log_pi = F.log_softmax(h[:, :, 0, :], dim=-1)                # [B,L,K]
        mu = torch.sigmoid(h[:, :, 1, :])                            # [B,L,K] in (0,1)
        s = self.s_min + F.softplus(h[:, :, 2, :])                   # [B,L,K] > s_min
        return log_pi, mu, s

    @staticmethod
    def _mixture(u, lp_l, mu_l, s_l):
        """G(u) and G'(u) for a [0,1] logistic mixture. u [B], params [B,K]."""
        z = (u.unsqueeze(-1) - mu_l) / s_l                           # [B,K]
        sig = torch.sigmoid(z)
        pi = lp_l.exp()
        G = (pi * sig).sum(-1)
        Gp = (pi * sig * (1 - sig) / s_l).sum(-1)
        return G, Gp

    def _layer_S_logSp(self, u, lp_l, mu_l, s_l):
        """One renormalised-mixture-CDF layer [0,1]→[0,1]: returns S(u), log S'(u)."""
        G_u, Gp_u = self._mixture(u, lp_l, mu_l, s_l)
        z0 = torch.zeros_like(u); z1 = torch.ones_like(u)
        G0, _ = self._mixture(z0, lp_l, mu_l, s_l)
        G1, _ = self._mixture(z1, lp_l, mu_l, s_l)
        denom = (G1 - G0).clamp_min(1e-12)
        S = ((G_u - G0) / denom).clamp(0.0, 1.0)
        logSp = Gp_u.clamp_min(1e-30).log() - denom.log()
        return S, logSp

    def _compose(self, x, log_pi, mu, s):
        """Map x∈[a,b]→u₀∈[0,1], compose the L layers. Returns (u_L, log p_std),
        where log p_std = Σ_l log S_l' − log(b−a) is the standardised-mass log
        density and u_L is the in-window CDF (∈[0,1])."""
        u = (x - self.a) / (self.b - self.a)
        logp = -math.log(self._width) + torch.zeros_like(x)
        for l in range(self.L):
            S, logSp = self._layer_S_logSp(u, log_pi[:, l], mu[:, l], s[:, l])
            logp = logp + logSp
            u = S
        return u, logp

    def _edges(self, log_pi, mu, s):
        """(log q, (log q)', (log q)'') at the window edges via autograd through
        the composition. Differentiable w.r.t. the conditioning when grad is on."""
        B = mu.shape[0]
        train_grad = torch.is_grad_enabled()
        out = {}
        with torch.enable_grad():
            for tag, edge in (("a", self.a), ("b", self.b)):
                xe = edge.expand(B).clone().requires_grad_(True)
                _, lp = self._compose(xe, log_pi, mu, s)
                d1 = torch.autograd.grad(lp.sum(), xe, create_graph=True,
                                         retain_graph=True)[0]
                d2 = torch.autograd.grad(d1.sum(), xe, create_graph=train_grad,
                                         retain_graph=True)[0]
                # Force a log-concave (monotone-decaying) tail: d2 ≤ −curv_floor.
                # Used by BOTH forward (_tail_logp) and log_cdf (α=−d2), so the
                # displayed density and its CDF stay exactly consistent. Exact C²
                # match when the true edge is already log-concave (the normal
                # case); a tiny relaxation only for a flat/convex edge.
                d2 = d2.clamp(max=-self.curv_floor)
                out[f"logq_{tag}"] = lp if train_grad else lp.detach()
                out[f"d1_{tag}"] = d1 if train_grad else d1.detach()
                out[f"d2_{tag}"] = d2 if train_grad else d2.detach()
        return out

    # ---- matched-Gaussian tail (erfcx-stable) ------------------------------
    @staticmethod
    def _tail_logp(t, logq, d1, d2):
        """log p_ext(edge+t) = logq + d1·t + ½ d2·t² (2nd-order Taylor of log q)."""
        return logq + d1 * t + 0.5 * d2 * t * t

    def _outward(self, U, logq, s_out, alpha):
        """Matched tail in the OUTWARD coordinate u≥0: g(u)=exp(logq+s_out·u−½α u²),
        s_out<0 (decaying). Returns (total=∫₀^∞, beyond=∫_U^∞), erfcx-stable."""
        root = (0.5 / alpha).sqrt()
        pref = (math.pi * 0.5 / alpha).sqrt()
        erfcx = torch.special.erfcx
        total = logq.exp() * pref * erfcx(-s_out * root)
        g_at_U = logq + s_out * U - 0.5 * alpha * U * U
        beyond = g_at_U.exp() * pref * erfcx((alpha * U - s_out) * root)
        return total, beyond

    # ---- public API --------------------------------------------------------
    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B]. ``x`` may be [B] or [B,1]."""
        x = x.reshape(-1)
        log_pi, mu, s = self._layer_params(c)
        a, b = self.a, self.b
        _, logq_in = self._compose(x, log_pi, mu, s)
        if not (bool((x < a).any()) or bool((x > b).any())):
            return logq_in                                  # all in-window (training)
        ed = self._edges(log_pi, mu, s)
        logp_hi = self._tail_logp(x - b, ed["logq_b"], ed["d1_b"], ed["d2_b"])
        logp_lo = self._tail_logp(x - a, ed["logq_a"], ed["d1_a"], ed["d2_a"])
        return torch.where(x > b, logp_hi, torch.where(x < a, logp_lo, logq_in))

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B]. F₀ = cumulative of the (in-window + tails)
        unnormalised density; F₀(b)−F₀(a)=1 exactly → window Z is exact."""
        x = x.reshape(-1)
        log_pi, mu, s = self._layer_params(c)
        a, b = self.a, self.b
        ed = self._edges(log_pi, mu, s)
        # d2 is already floored ≤ −curv_floor in _edges, so α=−d2 ≥ curv_floor > 0.
        alpha_a = -ed["d2_a"]
        alpha_b = -ed["d2_b"]
        s_out_a = -ed["d1_a"]                       # outward (decreasing x) slope
        s_out_b = ed["d1_b"]                        # outward (increasing x) slope
        u_a = (a - x).clamp_min(0.0)
        u_b = (x - b).clamp_min(0.0)
        m_lo, F_below_x = self._outward(u_a, ed["logq_a"], s_out_a, alpha_a)
        m_hi_tot, beyond_xb = self._outward(u_b, ed["logq_b"], s_out_b, alpha_b)
        Q, _ = self._compose(x, log_pi, mu, s)      # in-window CDF u_L ∈[0,1]
        Q = Q.clamp(0.0, 1.0)
        F_in = m_lo + Q
        F_hi = m_lo + 1.0 + (m_hi_tot - beyond_xb)
        Fc = torch.where(x > b, F_hi, torch.where(x < a, F_below_x, F_in))
        return Fc.clamp_min(1e-30).log()


# ---------------------------------------------------------------------------
# Self-test: normalisation, density↔CDF consistency, C² continuity, no-overflow.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    a, b = -3.7, 3.9
    flow = CompactMatchedFlow(n_cond=7, a=a, b=b, n_components=8, n_transforms=5).double()
    c = torch.randn(4, 7, dtype=torch.float64)
    xg = torch.linspace(a, b, 20001, dtype=torch.float64)
    for e in range(c.shape[0]):
        lp = flow(xg, c[e:e + 1].expand(xg.shape[0], -1))
        I = torch.trapz(lp.exp(), xg).item()
        Fa = flow.log_cdf(torch.tensor([a], dtype=torch.float64), c[e:e + 1]).exp().item()
        Fb = flow.log_cdf(torch.tensor([b], dtype=torch.float64), c[e:e + 1]).exp().item()
        print(f"event {e}: ∫_window p={I:.6f}  F(b)-F(a)={Fb-Fa:.6f}")
    xt = torch.tensor([a - 0.3, a - 0.01, a + 0.01, 0.0, b - 0.01, b + 0.01, b + 0.3],
                      dtype=torch.float64)
    ce = c[:1].expand(xt.shape[0], -1)
    h = 1e-5
    dF = (flow.log_cdf(xt + h, ce).exp() - flow.log_cdf(xt - h, ce).exp()) / (2 * h)
    p = flow(xt, ce).exp()
    print("rel.diff dF/dx vs p:",
          [f"{abs(d - pp) / max(pp, 1e-12):.1e}" for d, pp in zip(dF.tolist(), p.tolist())])
    # overflow stress: extreme conditioners
    f2 = CompactMatchedFlow(7, a, b, n_components=8, n_transforms=5).double()
    with torch.no_grad():
        for p in f2.parameters():
            p.normal_(0, 3.0)         # extreme conditioners (stress the tails)
    cc = torch.randn(64, 7, dtype=torch.float64)
    xs = torch.linspace(-6, 6, 200, dtype=torch.float64)
    ok = all(torch.isfinite(f2(xs, cc[e:e + 1].expand(200, -1))).all()
             and torch.isfinite(f2.log_cdf(xs, cc[e:e + 1].expand(200, -1))).all()
             for e in range(64))
    print(f"overflow stress (64 extreme conditioners): all finite = {ok}")
    # gradient through θ-like conditioning works (operator differentiability)
    cg = torch.randn(8, 7, dtype=torch.float64, requires_grad=True)
    out = flow(torch.full((8,), b + 0.05, dtype=torch.float64), cg).sum()
    g = torch.autograd.grad(out, cg)[0]
    print(f"tail grad wrt conditioning finite = {bool(torch.isfinite(g).all())}")


if __name__ == "__main__":
    _selftest()
