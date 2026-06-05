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
            with a uniform base, so it is intrinsically compact). Two layer types:

            ``logistic`` (default): a logistic-mixture CDF renormalised to [0,1]:
                S_l(u) = (G_l(u) − G_l(0)) / (G_l(1) − G_l(0)),
                G_l(u) = Σ_k π_k σ((u−μ_k)/s_k),  μ_k∈(0,1), s_k>0,
            with x↦u₀=(x−a)/(b−a) (uniform input warp).

            ``bernstein``: a monotone Bernstein polynomial with PINNED endpoints,
                S_l(u) = Σ_{k=1}^{M} θ_k b_{k,M}(u),  θ_k = Σ_{j≤k} δ_j,
                δ = (1−ε)·softmax(conditioner) + ε/M  (so S'_l ≥ ε exactly),
            i.e. θ_0=0, θ_M=1 ALGEBRAICALLY — no renormalisation division, no
            scale floors, positivity by construction. The input warp is the
            FIXED truncated-N(0,1) CDF u₀ = (Φ(x)−Φ(a))/(Φ(b)−Φ(a)) — exactly
            equivalent to a frozen extra layer, giving a PEAKED initial density
            (the standardised mass has μ≈0, σ≈1 by construction, so the base is
            parameter-free) so the fixed-position polynomial basis only has to
            fit the smooth residual warp, not localise the J/ψ peak from a flat
            start. Edge derivatives of log q are CLOSED FORM (Bernstein
            endpoint derivatives + chain rule; the prewarp contributes exactly
            (log w')' = −x, (log w')'' = −1), so the matched tails need NO
            autograd — cheaper, and vmap/compile-friendly.

            With the uniform base the composition F = S_L∘…∘S₁∘u₀ : [a,b]→[0,1]
            IS the in-window CDF (F(a)=0, F(b)=1) and the density q = F'
            integrates to 1 over [a,b] EXACTLY. C∞ interior, expressive via
            depth (few components / moderate degree per layer, like gf).

Tails:      beyond each edge the log-density is continued by its 2nd-order Taylor
            in log q — a Gaussian (downward log-parabola) matched to log q,
            (log q)', (log q)'' at the edge (AUTOGRAD through the composition for
            logistic; ANALYTIC for bernstein). C² across the boundary, monotone-
            decaying, closed-form CDF (erf, via erfcx for overflow safety). A
            FIXED controlled extension of the data-constrained edge — no humps,
            no near-boundary bias.

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


def infer_learn_weights(state_dict, n_components: int,
                        prefix: str = "flow.conditioners.0.",
                        default: bool = False) -> bool:
    """Infer whether a compact-flow checkpoint used LEARNABLE mixture weights
    (3K-output conditioner heads) or fixed equal weights (2K), from the final
    conditioner layer's output dimension. For checkpoints predating the
    ``learn_weights`` option (which were always 3K)."""
    cand = [(k, v) for k, v in state_dict.items()
            if k.startswith(prefix) and k.endswith(".weight")]
    if not cand:
        return default

    def _idx(k):
        try:
            return int(k[len(prefix):].split(".")[0])
        except (ValueError, IndexError):
            return -1

    _, v = max(cand, key=lambda kv: _idx(kv[0]))
    out = int(v.shape[0])
    if out == 3 * int(n_components):
        return True
    if out == 2 * int(n_components):
        return False
    return default


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
        learn_weights: bool = False,
        layer_type: str = "logistic",
        bernstein_degree: int = 16,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        if layer_type not in ("logistic", "bernstein"):
            raise ValueError(f"layer_type must be 'logistic' or 'bernstein'; "
                             f"got {layer_type!r}")
        self.layer_type = str(layer_type)
        self.K = int(n_components)
        self.M = int(bernstein_degree)
        self.L = int(n_transforms)
        # learn_weights=False (default) fixes π_k = 1/K — the per-layer mixture
        # has only means+scales (2K conditioner outputs), exactly matching gf's
        # equal-weight Gaussianization layers. True adds learnable per-component
        # weights (3K outputs, the original behaviour). Logistic-only.
        self.learn_weights = bool(learn_weights)
        if self.layer_type == "bernstein":
            if self.learn_weights:
                raise ValueError("learn_weights applies to logistic layers only")
            if self.M < 3:
                raise ValueError(f"bernstein_degree must be ≥ 3 (need S''' at "
                                 f"the edges for the matched tails); got {self.M}")
        # Buffers kept ONLY for checkpoint compatibility (they are persistent
        # and present in existing state_dicts). All maths/comparisons use the
        # EXACT python-float scalars below: the buffers were created at the
        # default dtype (fp32), so under model.double() they upcast the
        # fp32-ROUNDED edges — an exact-double edge evaluation (decompose at
        # A=0, fit norm at θ→0) then lands 1e-7 INSIDE the flow's internal
        # window, taking the in-window branch whose x-gradient the bernstein
        # pow-safety clamp zeroes. Scalars are exact doubles (and round
        # consistently with the tensor dtype in fp32 ops).
        self.register_buffer("a", torch.tensor(float(a)))
        self.register_buffer("b", torch.tensor(float(b)))
        self._af = float(a)
        self._bf = float(b)
        self._width = float(b - a)
        # logistic: min component scale in [0,1] units; bernstein: EXACT lower
        # bound on every layer slope S' ≥ s_min (via the δ floor) — both keep
        # log S' bounded so the density never collapses to −inf in-window.
        self.s_min = float(s_min_frac)
        self.curv_floor = float(curv_floor)     # min |(log p)''| in std-mass units

        if self.layer_type == "bernstein":
            # Bernstein basis binomials (exact ints, float64 so an fp64 model
            # keeps them exact; non-persistent — rebuilt at __init__ and self-
            # healed in _cond_params after any model-wide fp32 cast).
            self.register_buffer("_binM", torch.tensor(
                [math.comb(self.M, k) for k in range(self.M + 1)],
                dtype=torch.float64), persistent=False)
            self.register_buffer("_binM1", torch.tensor(
                [math.comb(self.M - 1, k) for k in range(self.M)],
                dtype=torch.float64), persistent=False)
            # Fixed truncated-N(0,1) prewarp constants (python floats — dtype-
            # independent). The truncation range is the window itself (forced:
            # the prewarp must be a [a,b]→[0,1] bijection with pinned endpoints
            # or Z≡1 breaks); μ=0, σ=1 because the standardisation stats come
            # from the same windowed sample, so the standardised mass is
            # moment-matched to N(0,1) for free.
            phi_a = 0.5 * (1.0 + math.erf(float(a) / _SQRT2))
            phi_b = 0.5 * (1.0 + math.erf(float(b) / _SQRT2))
            self._phi_a = float(phi_a)
            self._D = float(phi_b - phi_a)
            self._log_norm = 0.5 * math.log(2.0 * math.pi) + math.log(self._D)

        # One INDEPENDENT conditioner MLP per layer, mirroring gf's per-transform
        # hyper-network (so conditioning capacity and parameter count match gf),
        # rather than a shared body + linear head. Per-layer outputs:
        # logistic: [raw_mu|raw_s] (2K, equal weights — default) or
        # [logit_w|raw_mu|raw_s] (3K, learnable weights);
        # bernstein: M raw increment logits (δ = softmax → θ = cumsum; zero bias
        # → δ uniform → S = identity → init density = the truncated N(0,1)).
        if self.layer_type == "bernstein":
            out_dim = self.M
        else:
            nblk = 3 if self.learn_weights else 2
            out_dim = nblk * self.K
            frac = torch.linspace(0.08, 0.92, self.K)
            mu_b = torch.log(frac / (1 - frac))                      # sigmoid⁻¹
            s0 = max(1.0 / self.K, self.s_min)
            s_b = math.log(math.expm1(max(s0 - self.s_min, 1e-4)))   # softplus⁻¹
        self.conditioners = nn.ModuleList()
        for _ in range(self.L):
            seq: list[nn.Module] = []
            d = int(n_cond)
            for _ in range(int(n_layers)):
                seq += [nn.Linear(d, hidden_features), activation()]
                d = hidden_features
            final = nn.Linear(d, out_dim)
            nn.init.normal_(final.weight, std=1e-3)                  # ~c-independent init
            with torch.no_grad():
                final.bias.zero_()
                if self.layer_type == "logistic":
                    b = final.bias.view(nblk, self.K)
                    b[-2, :] = mu_b                                  # spread means in (0,1)
                    b[-1, :] = s_b                                   # moderate scales
            seq.append(final)
            self.conditioners.append(nn.Sequential(*seq))

    # ---- per-layer parameters (arch-specific tuple) -------------------------
    def _cond_params(self, c: torch.Tensor):
        """Conditioner outputs → the per-layer transform parameters, as a tuple
        consumed by _compose/_edges. logistic: (log_pi, mu, s) each [B,L,K];
        bernstein: (delta,) [B,L,M] — positive increments with Σ_k δ_k = 1 and
        δ_k ≥ s_min/M, so θ = cumsum(δ) is increasing with θ_M = 1 and every
        layer slope S' = M·Σ δ_{k+1} b_{k,M−1} ≥ s_min."""
        if self.layer_type == "bernstein":
            if self._binM.dtype != torch.float64:
                # Self-heal after a model-wide fp32 cast (cf. nce GL buffers):
                # the binomials are exact integers; keep them float64.
                self._binM = torch.tensor(
                    [math.comb(self.M, k) for k in range(self.M + 1)],
                    dtype=torch.float64, device=self._binM.device)
                self._binM1 = torch.tensor(
                    [math.comb(self.M - 1, k) for k in range(self.M)],
                    dtype=torch.float64, device=self._binM1.device)
            h = torch.stack([cond(c) for cond in self.conditioners], dim=1)
            delta = ((1.0 - self.s_min) * torch.softmax(h, dim=-1)
                     + self.s_min / self.M)                          # [B,L,M]
            delta = delta / delta.sum(-1, keepdim=True)              # θ_M = 1 (ulp)
            return (delta,)
        nblk = 3 if self.learn_weights else 2
        h = torch.stack([cond(c) for cond in self.conditioners], dim=1)  # [B,L,nblk·K]
        h = h.view(-1, self.L, nblk, self.K)                         # [B,L,nblk,K]
        if self.learn_weights:
            log_pi = F.log_softmax(h[:, :, 0, :], dim=-1)            # [B,L,K]
        else:
            log_pi = h.new_full(h.shape[:2] + (self.K,), -math.log(self.K))
        mu = torch.sigmoid(h[:, :, -2, :])                           # [B,L,K] in (0,1)
        s = self.s_min + F.softplus(h[:, :, -1, :])                  # [B,L,K] > s_min
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

    @staticmethod
    def _bern_basis(u, M, binom):
        """Bernstein basis b_{k,M}(u) = C(M,k) u^k (1−u)^{M−k}, [B] → [B, M+1].
        ``u`` must be strictly inside (0,1): pow's backward for the k=0 / k=M
        terms is 0·u^{−1}, which is 0 for u>0 but NaN at exactly 0/1 — callers
        clamp by an ulp."""
        k = torch.arange(M + 1, device=u.device, dtype=u.dtype)
        uu = u.unsqueeze(-1)
        return binom.to(u.dtype) * uu.pow(k) * (1.0 - uu).pow(M - k)

    def _bern_S_logSp(self, u, delta_l):
        """One pinned monotone Bernstein layer [0,1]→[0,1]: S = Σ_k θ_k b_{k,M},
        θ = cumsum(δ) (θ_0=0, θ_M=1 algebraically); S' = M Σ_k δ_{k+1} b_{k,M−1}
        ≥ s_min by the δ floor. Returns S(u), log S'(u)."""
        eps = torch.finfo(u.dtype).eps
        uc = u.clamp(eps, 1.0 - eps)                 # pow-backward safety (b_0/b_M)
        theta = torch.cumsum(delta_l, dim=-1)        # [B,M] = θ_1..θ_M (θ_0=0)
        bM = self._bern_basis(uc, self.M, self._binM)        # [B,M+1]
        bM1 = self._bern_basis(uc, self.M - 1, self._binM1)  # [B,M]
        S = (theta * bM[:, 1:]).sum(-1).clamp(0.0, 1.0)
        Sp = float(self.M) * (delta_l * bM1).sum(-1)
        return S, Sp.clamp_min(1e-30).log()

    def _compose(self, x, params):
        """Map x∈[a,b]→u₀∈[0,1] (affine for logistic; the fixed truncated-N(0,1)
        CDF prewarp for bernstein), compose the L layers. Returns (u_L, log
        p_std), where log p_std = log u₀'(x) + Σ_l log S_l' is the standardised-
        mass log density and u_L is the in-window CDF (∈[0,1])."""
        if self.layer_type == "bernstein":
            (delta,) = params
            u = ((_phi(x) - self._phi_a) / self._D).clamp(0.0, 1.0)
            # log w'(x) = log φ(x) − log D  (φ = standard-normal pdf)
            logp = -0.5 * x * x - self._log_norm
            for l in range(self.L):
                S, logSp = self._bern_S_logSp(u, delta[:, l])
                logp = logp + logSp
                u = S
            return u, logp
        log_pi, mu, s = params
        u = (x - self._af) / (self._bf - self._af)
        logp = -math.log(self._width) + torch.zeros_like(x)
        for l in range(self.L):
            S, logSp = self._layer_S_logSp(u, log_pi[:, l], mu[:, l], s[:, l])
            logp = logp + logSp
            u = S
        return u, logp

    def _edges(self, params):
        """(log q, (log q)', (log q)'') at the window edges — analytic for
        bernstein, autograd through the composition for logistic. The d2 ≤
        −curv_floor clamp is shared (see _edges_autograd).

        ALWAYS computed in fp64: when the trained edge density is tiny
        (log q ~ −100s, routine for sharp peaks far from a window edge), the
        per-layer S' factors underflow to 0 in fp32 and the autograd second
        derivative hits 0/0 → NaN d2 (observed on a real fp32 checkpoint for
        ~half the conditioning points → non-finite stage-2 loss). The upcast
        is differentiable, the cost is per-event (not per-node/grid), and the
        results are cast back to the input dtype (d1/d2 range-clamped first so
        the downcast itself cannot create ±inf → 0·inf=NaN downstream)."""
        dt = params[0].dtype
        if dt != torch.float64:
            params = tuple(p.double() for p in params)
        out = (self._edges_bernstein(params) if self.layer_type == "bernstein"
               else self._edges_autograd(params))
        if dt != torch.float64:
            fmax = float(torch.finfo(dt).max) / 16.0
            out = {k: (v if k.startswith("logq") else v.clamp(-fmax, fmax)).to(dt)
                   for k, v in out.items()}
        return out

    def _edges_autograd(self, params):
        """(log q, (log q)', (log q)'') at the window edges via autograd through
        the composition. Differentiable w.r.t. the conditioning when grad is on."""
        B = params[0].shape[0]
        train_grad = torch.is_grad_enabled()
        out = {}
        with torch.enable_grad():
            for tag, edge in (("a", self._af), ("b", self._bf)):
                xe = params[0].new_full((B,), edge).requires_grad_(True)
                _, lp = self._compose(xe, params)
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

    def _edges_bernstein(self, params):
        """ANALYTIC (log q, (log q)', (log q)'') at the window edges for the
        bernstein layers — no autograd (cheaper; vmap/compile-friendly; still
        differentiable w.r.t. the conditioning, since everything is a smooth
        closed-form function of the δ).

        At the edges every layer input is exactly 0 (x=a) or 1 (x=b), where the
        Bernstein endpoint derivatives are closed form:
            S'(0)   = M δ_1                S'(1)   = M δ_M
            S''(0)  = M(M−1)(δ_2 − δ_1)    S''(1)  = M(M−1)(δ_M − δ_{M−1})
            S'''(0) = M(M−1)(M−2)(δ_3 − 2δ_2 + δ_1)
            S'''(1) = M(M−1)(M−2)(δ_M − 2δ_{M−1} + δ_{M−2})
        and the truncated-normal prewarp contributes exactly
            (log w')'(x) = −x,  (log w')''(x) = −1,
            w'(x) = φ(x)/D,  w''(x) = −x·w'(x).
        Chain rule through t_l = S_l(t_{l−1}) with log q = log w' + Σ log S_l':
            (log q)'  += (S''/S')(t)·t',
            (log q)'' += (S'''/S' − (S''/S')²)(t)·t'² + (S''/S')(t)·t'',
            t''_new = S''·t'² + S'·t'',  t'_new = S'·t'."""
        (delta,) = params
        B = delta.shape[0]
        Mf = float(self.M)
        c2 = Mf * (Mf - 1.0)
        c3 = Mf * (Mf - 1.0) * (Mf - 2.0)
        out = {}
        for tag, xe in (("a", self._af), ("b", self._bf)):
            at_b = (tag == "b")
            lw = -0.5 * xe * xe - self._log_norm     # log w'(xe)
            lq = delta.new_full((B,), lw)
            d1 = delta.new_full((B,), -xe)           # (log w')' at the edge
            d2 = delta.new_full((B,), -1.0)          # (log w')'' at the edge
            tp = delta.new_full((B,), math.exp(lw))  # t' = w'(xe)
            tpp = tp * (-xe)                         # t'' = w''(xe) = −xe·w'(xe)
            for l in range(self.L):
                dl = delta[:, l]
                if at_b:
                    s1 = Mf * dl[:, -1]
                    s2 = c2 * (dl[:, -1] - dl[:, -2])
                    s3 = c3 * (dl[:, -1] - 2.0 * dl[:, -2] + dl[:, -3])
                else:
                    s1 = Mf * dl[:, 0]
                    s2 = c2 * (dl[:, 1] - dl[:, 0])
                    s3 = c3 * (dl[:, 2] - 2.0 * dl[:, 1] + dl[:, 0])
                r1 = s2 / s1                          # s1 ≥ s_min > 0 by the floor
                r2 = s3 / s1
                lq = lq + s1.log()
                d2 = d2 + (r2 - r1 * r1) * tp * tp + r1 * tpp
                d1 = d1 + r1 * tp
                tpp = s2 * tp * tp + s1 * tpp         # before tp is overwritten
                tp = s1 * tp
            # Shared log-concavity clamp (see _edges_autograd).
            d2 = d2.clamp(max=-self.curv_floor)
            out[f"logq_{tag}"] = lq
            out[f"d1_{tag}"] = d1
            out[f"d2_{tag}"] = d2
        return out

    # ---- matched-Gaussian tail (erfcx-stable) ------------------------------
    @staticmethod
    def _tail_logp(t, logq, d1, d2):
        """log p_ext(edge+t) = logq + d1·t + ½ d2·t² (2nd-order Taylor of log q)."""
        return logq + d1 * t + 0.5 * d2 * t * t

    @staticmethod
    def _log_erfcx(z):
        """log erfcx(z) for any real z. erfcx(z) = 2e^{z²} − erfcx(−z), so for
        z ≤ −6 it is 2e^{z²} to relative error < e^{−36} — but erfcx itself
        overflows there (fp32 already at z ≈ −9.4). Branch in log space; the
        clamp keeps the unselected branch NaN-free in forward AND backward."""
        safe = torch.special.erfcx(z.clamp_min(-6.0)).log()
        return torch.where(z < -6.0, z * z + math.log(2.0), safe)

    def _outward(self, U, logq, s_out, alpha):
        """Matched tail in the OUTWARD coordinate u≥0: g(u)=exp(logq+s_out·u−½α u²).
        Returns (total=∫₀^∞, beyond=∫_U^∞). erfcx is only directly evaluable for
        DECAYING outward slopes (s_out<0); a trained edge can come out RISING
        (s_out>0), where the closed form ∝ e^{s²/2α} overflows — assemble both
        factors in log space and cap the exponent: an e^{80} tail mass is equally
        pathological either way, but the NLL stays finite (and huge) so the
        optimizer steers away instead of crashing on NaN."""
        root = (0.5 / alpha).sqrt()
        log_pref = 0.5 * (math.log(math.pi * 0.5) - alpha.log())
        total = (logq + log_pref + self._log_erfcx(-s_out * root)
                 ).clamp(max=80.0).exp()
        g_at_U = logq + s_out * U - 0.5 * alpha * U * U
        beyond = (g_at_U + log_pref + self._log_erfcx((alpha * U - s_out) * root)
                  ).clamp(max=80.0).exp()
        return total, beyond

    # ---- public API --------------------------------------------------------
    def forward_inwindow(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c) for x GUARANTEED ∈ [a, b] — the stage-1 training
        case (the masses are window-selected, so standardised values lie in
        [a, b] by construction and ``forward``'s tail branch is dead code).
        Numerically identical to ``forward`` there, but skips the
        data-dependent short-circuit (``bool((x<a).any())`` graph-breaks
        dynamo) and the matched tails (logistic ``_edges``' inner
        ``autograd.grad`` cannot be traced), so it is
        torch.compile(fullgraph=True)-friendly. Do NOT use where x can stray
        outside the window (operator un-kick, CDF boundaries) — the tail
        extension would be silently wrong."""
        x = x.reshape(-1)
        params = self._cond_params(c)
        _, logp = self._compose(x, params)
        return logp

    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B]. ``x`` may be [B] or [B,1]."""
        x = x.reshape(-1)
        params = self._cond_params(c)
        a, b = self._af, self._bf
        _, logq_in = self._compose(x, params)
        # STRICT interior fast path / INCLUSIVE tail branches: at exactly x=a/b
        # the tail value equals the in-window value (C² match), but the tail's
        # x-derivative is clamp-free — the bernstein layers' pow-safety clamp
        # sits exactly at its boundary for u=0/1 and silently ZEROES the
        # autograd derivative of the in-window branch there (the decompose
        # g_norm bug: dZ/ds at the un-shifted window edges came out 0).
        if not (bool((x <= a).any()) or bool((x >= b).any())):
            return logq_in                                  # all in-window (training)
        ed = self._edges(params)
        logp_hi = self._tail_logp(x - b, ed["logq_b"], ed["d1_b"], ed["d2_b"])
        logp_lo = self._tail_logp(x - a, ed["logq_a"], ed["d1_a"], ed["d2_a"])
        return torch.where(x >= b, logp_hi, torch.where(x <= a, logp_lo, logq_in))

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B]. F₀ = cumulative of the (in-window + tails)
        unnormalised density; F₀(b)−F₀(a)=1 exactly → window Z is exact."""
        x = x.reshape(-1)
        params = self._cond_params(c)
        a, b = self._af, self._bf
        ed = self._edges(params)
        # d2 is already floored ≤ −curv_floor in _edges, so α=−d2 ≥ curv_floor > 0.
        alpha_a = -ed["d2_a"]
        alpha_b = -ed["d2_b"]
        s_out_a = -ed["d1_a"]                       # outward (decreasing x) slope
        s_out_b = ed["d1_b"]                        # outward (increasing x) slope
        u_a = (a - x).clamp_min(0.0)
        u_b = (x - b).clamp_min(0.0)
        m_lo, F_below_x = self._outward(u_a, ed["logq_a"], s_out_a, alpha_a)
        m_hi_tot, beyond_xb = self._outward(u_b, ed["logq_b"], s_out_b, alpha_b)
        Q, _ = self._compose(x, params)             # in-window CDF u_L ∈[0,1]
        Q = Q.clamp(0.0, 1.0)
        F_in = m_lo + Q
        F_hi = m_lo + 1.0 + (m_hi_tot - beyond_xb)
        # INCLUSIVE tail branches (see forward): at exactly x=a/b the tail
        # F is both exact in value (F_hi(b) = m_lo + 1 algebraically, vs an
        # O(M·ε) clamp residue in Q) and has the correct clamp-free derivative
        # p(edge) — the in-window branch's derivative is zeroed there by the
        # bernstein pow-safety clamp (decompose-g_norm / fit-norm-gradient at
        # θ≈0 would silently lose the p(b)−p(a) term).
        Fc = torch.where(x >= b, F_hi, torch.where(x <= a, F_below_x, F_in))
        return Fc.clamp_min(1e-30).log()


# ---------------------------------------------------------------------------
# Self-test: normalisation, density↔CDF consistency, C² continuity, no-overflow.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    a, b = -3.7, 3.9
    for lw in (False, True):
        flow = CompactMatchedFlow(n_cond=7, a=a, b=b, n_components=8,
                                  n_transforms=5, learn_weights=lw).double()
        c = torch.randn(4, 7, dtype=torch.float64)
        xg = torch.linspace(a, b, 20001, dtype=torch.float64)
        for e in range(c.shape[0]):
            lp = flow(xg, c[e:e + 1].expand(xg.shape[0], -1))
            I = torch.trapz(lp.exp(), xg).item()
            Fa = flow.log_cdf(torch.tensor([a], dtype=torch.float64), c[e:e + 1]).exp().item()
            Fb = flow.log_cdf(torch.tensor([b], dtype=torch.float64), c[e:e + 1]).exp().item()
            print(f"learn_weights={lw} event {e}: ∫_window p={I:.6f}  F(b)-F(a)={Fb-Fa:.6f}")
        xt = torch.tensor([a - 0.3, a - 0.01, a + 0.01, 0.0, b - 0.01, b + 0.01, b + 0.3],
                          dtype=torch.float64)
        ce = c[:1].expand(xt.shape[0], -1)
        h = 1e-5
        dF = (flow.log_cdf(xt + h, ce).exp() - flow.log_cdf(xt - h, ce).exp()) / (2 * h)
        p = flow(xt, ce).exp()
        print(f"learn_weights={lw} rel.diff dF/dx vs p:",
              [f"{abs(d - pp) / max(pp, 1e-12):.1e}" for d, pp in zip(dF.tolist(), p.tolist())])
        # forward_inwindow ≡ forward for in-window x (the compile fast path)
        xin = a + (b - a) * torch.rand(257, dtype=torch.float64)
        cin = c[:1].expand(257, -1)
        assert torch.equal(flow.forward_inwindow(xin, cin), flow(xin, cin)), \
            "forward_inwindow disagrees with forward in-window"
        # checkpoint-shape inference round-trip
        sd = {f"flow.{k}": v for k, v in flow.state_dict().items()}
        assert infer_learn_weights(sd, 8) == lw, "infer_learn_weights round-trip failed"
    print("forward_inwindow ≡ forward (in-window); infer_learn_weights "
          "round-trip OK for both modes")

    # ---- bernstein layers -------------------------------------------------
    fb = CompactMatchedFlow(n_cond=7, a=a, b=b, n_transforms=5,
                            layer_type="bernstein", bernstein_degree=16).double()
    c = torch.randn(4, 7, dtype=torch.float64)
    # 1) zero-init head → identity layers → density == truncated N(0,1)
    with torch.no_grad():
        for cond in fb.conditioners:
            cond[-1].weight.zero_(); cond[-1].bias.zero_()
    xg = torch.linspace(a + 1e-9, b - 1e-9, 20001, dtype=torch.float64)
    ce = c[:1].expand(xg.shape[0], -1)
    lp = fb(xg, ce)
    D = 0.5 * (math.erf(b / _SQRT2) - math.erf(a / _SQRT2))
    lp_tn = -0.5 * xg**2 - 0.5 * math.log(2 * math.pi) - math.log(D)
    err_tn = float((lp - lp_tn).abs().max())
    print(f"bernstein zero-init: max |logp − logTN(0,1)| = {err_tn:.2e}")
    assert err_tn < 1e-9
    # 2) perturbed: normalisation, CDF exactness, dF/dx ≡ p incl. tails
    with torch.no_grad():
        for p_ in fb.parameters():
            p_.add_(0.5 * torch.randn_like(p_))
    for e in range(c.shape[0]):
        ci = c[e:e + 1]
        lp = fb(xg, ci.expand(xg.shape[0], -1))
        I = torch.trapz(lp.exp(), xg).item()
        Fa = fb.log_cdf(torch.tensor([a], dtype=torch.float64), ci).exp().item()
        Fb = fb.log_cdf(torch.tensor([b], dtype=torch.float64), ci).exp().item()
        print(f"bernstein event {e}: ∫_window p={I:.6f}  F(b)-F(a)={Fb-Fa:.6f}")
        assert abs(I - 1.0) < 5e-4 and abs(Fb - Fa - 1.0) < 1e-12
    xt = torch.tensor([a - 0.3, a - 0.01, a + 0.01, 0.0, b - 0.01, b + 0.01, b + 0.3],
                      dtype=torch.float64)
    ce = c[:1].expand(xt.shape[0], -1)
    h = 1e-5
    dF = (fb.log_cdf(xt + h, ce).exp() - fb.log_cdf(xt - h, ce).exp()) / (2 * h)
    p = fb(xt, ce).exp()
    rel = [abs(d - pp) / max(pp, 1e-12) for d, pp in zip(dF.tolist(), p.tolist())]
    print("bernstein rel.diff dF/dx vs p:", [f"{r:.1e}" for r in rel])
    assert max(rel) < 1e-6
    # 3) ANALYTIC edges == autograd reference (the headline closed-form check).
    #    NB _edges_autograd at exactly x=a/b is NOT a valid reference here: the
    #    layer input is exactly 0/1 there, where _bern_S_logSp's pow-safety
    #    clamp zeroes the gradient flow. Instead autograd the LAYER CHAIN
    #    λ(u₀) = Σ_l log S_l' at u₀ = ε (inside the clamp; error O(ε)) and add
    #    the prewarp chain rule analytically:
    #      d1 = −x + λ'(0)·w'(x),  d2 = −1 + λ''(0)·w'(x)² + λ'(0)·w''(x).
    prm = fb._cond_params(c)
    (delta,) = prm
    ana = fb._edges_bernstein(prm)
    # u₀ offset: the reference is evaluated at ε inside the edge, so it carries
    # O(ε·next-derivative) error — the chain amplifies derivatives by ~(M·δ̄)^L,
    # so keep ε tiny and compare RELATIVELY.
    for tag, xe, u0val in (("a", a, 1e-12), ("b", b, 1.0 - 1e-12)):
        u0 = torch.full((c.shape[0],), u0val, dtype=torch.float64,
                        requires_grad=True)
        u, lam = u0, torch.zeros_like(u0)
        for l in range(fb.L):
            S, logSp = fb._bern_S_logSp(u, delta[:, l])
            lam = lam + logSp
            u = S
        l1 = torch.autograd.grad(lam.sum(), u0, create_graph=True)[0]
        l2 = torch.autograd.grad(l1.sum(), u0)[0]
        wp = math.exp(-0.5 * xe * xe - fb._log_norm)
        wpp = -xe * wp
        lq_ref = (-0.5 * xe * xe - fb._log_norm) + lam.detach()
        d1_ref = -xe + l1.detach() * wp
        d2_ref = (-1.0 + l2.detach() * wp * wp
                  + l1.detach() * wpp).clamp(max=-fb.curv_floor)
        for k, ref in (("logq", lq_ref), ("d1", d1_ref), ("d2", d2_ref)):
            d = float(((ana[f"{k}_{tag}"] - ref).abs()
                       / ref.abs().clamp_min(1.0)).max())
            # The REFERENCE carries O(u₀·λ'') truncation (verified: it
            # converges to the analytic values as u₀ → 0 until fp64 round-off
            # floors it) — 1e-4 still rules out any sign/order/coefficient
            # error in the closed forms, which would be O(1).
            assert d < 1e-4, (k, tag, d)
    print("bernstein analytic edges == autograd-chain reference "
          "(logq, d1, d2 at both edges; rel <1e-4, reference-limited)")
    # 3b) EXACT-EDGE shift gradient: dZ/ds at s=0 with Z = F(b+s) − F(a+s)
    #     must equal p(b) − p(a) (the decompose g_norm path; the bernstein
    #     pow-safety clamp used to zero it through the in-window branch).
    for lt in ("bernstein", "logistic"):
        fz = (fb if lt == "bernstein" else
              CompactMatchedFlow(n_cond=7, a=a, b=b, n_transforms=5).double())
        cz = c[:3]
        s0 = torch.zeros((), dtype=torch.float64, requires_grad=True)
        Z = (fz.log_cdf(torch.full((3,), b, dtype=torch.float64) + s0, cz).exp()
             - fz.log_cdf(torch.full((3,), a, dtype=torch.float64) + s0, cz).exp())
        gz = torch.autograd.grad(Z.sum(), s0)[0]
        pe = (fz(torch.full((3,), b, dtype=torch.float64), cz).exp()
              - fz(torch.full((3,), a, dtype=torch.float64), cz).exp()).sum()
        rel = float((gz - pe).abs() / pe.abs().clamp_min(1e-12))
        assert rel < 1e-9, (lt, float(gz), float(pe))
        print(f"{lt} exact-edge dZ/ds == p(b)-p(a) (rel {rel:.1e})")
    # 4) forward_inwindow ≡ forward in-window; conditioning grads through the
    #    ANALYTIC tails (the operator path)
    xin = a + (b - a) * torch.rand(257, dtype=torch.float64)
    cin = c[:1].expand(257, -1)
    assert torch.equal(fb.forward_inwindow(xin, cin), fb(xin, cin))
    cg = torch.randn(8, 7, dtype=torch.float64, requires_grad=True)
    out = fb(torch.full((8,), b + 0.05, dtype=torch.float64), cg).sum()
    out = out + fb.log_cdf(torch.full((8,), a - 0.05, dtype=torch.float64), cg).sum()
    g = torch.autograd.grad(out, cg)[0]
    assert torch.isfinite(g).all()
    print("bernstein forward_inwindow ≡ forward; analytic-tail grads finite")
    # 5) overflow stress (extreme conditioners) + fp32 self-heal of binomials
    f3 = CompactMatchedFlow(7, a, b, n_transforms=5, layer_type="bernstein",
                            bernstein_degree=16).double()
    with torch.no_grad():
        for p_ in f3.parameters():
            p_.normal_(0, 3.0)
    cc = torch.randn(64, 7, dtype=torch.float64)
    xs = torch.linspace(-6, 6, 200, dtype=torch.float64)
    ok = all(torch.isfinite(f3(xs, cc[e:e + 1].expand(200, -1))).all()
             and torch.isfinite(f3.log_cdf(xs, cc[e:e + 1].expand(200, -1))).all()
             for e in range(64))
    print(f"bernstein overflow stress (64 extreme conditioners): all finite = {ok}")
    assert ok
    f3.float(); f3.double()                          # binomials downcast → heal
    _ = f3(xs, cc[:1].expand(200, -1))
    assert f3._binM.dtype == torch.float64
    print("bernstein binomial buffers self-heal after fp32 round-trip")
    flow = CompactMatchedFlow(n_cond=7, a=a, b=b, n_components=8, n_transforms=5).double()
    c = torch.randn(4, 7, dtype=torch.float64)
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

    # 6) fp32 stress — the REAL fit-stage precision. Two failure modes guarded:
    #    (a) NaN d2 from _edges_autograd: fp32 underflow of the per-layer S'
    #        factors when the edge density is tiny (a trained checkpoint hit
    #        this for ~44% of conditioning points → non-finite stage-2 loss);
    #        _edges now computes in fp64.
    #    (b) erfcx overflow in _outward for outward-RISING trained edge slopes
    #        (erfcx(−y) ≈ 2e^{y²} overflows fp32 at y ≳ 9.4); now log-space.
    for lt in ("logistic", "bernstein"):
        f32 = CompactMatchedFlow(7, a, b, n_components=8, n_transforms=4,
                                 learn_weights=(lt == "logistic"),
                                 layer_type=lt, bernstein_degree=16)
        with torch.no_grad():
            for p_ in f32.parameters():
                p_.normal_(0, 3.0)
        cc32 = torch.randn(512, 7)
        xs32 = torch.linspace(a - 1.0, b + 1.0, 101)
        nbad = 0
        for e in range(512):
            ce = cc32[e:e + 1].expand(101, -1)
            nbad += int((~torch.isfinite(f32(xs32, ce))).sum())
            nbad += int((~torch.isfinite(f32.log_cdf(xs32, ce))).sum())
        ed32 = f32._edges(f32._cond_params(cc32))
        edfin = all(bool(torch.isfinite(v).all()) for v in ed32.values())
        assert all(v.dtype == torch.float32 for v in ed32.values())
        print(f"fp32 stress [{lt}]: non-finite forward/log_cdf values = {nbad}, "
              f"edges all finite = {edfin}")
        assert nbad == 0 and edfin
    # _outward: finite for a rising edge slope with curvature at the floor
    for dt in (torch.float32, torch.float64):
        tot, bey = flow._outward(torch.tensor([0.0, 0.5], dtype=dt),
                                 torch.tensor([-2.0, -2.0], dtype=dt),
                                 torch.tensor([5.0, 5.0], dtype=dt),
                                 torch.tensor([1e-3, 1e-3], dtype=dt))
        assert bool(torch.isfinite(tot).all() and torch.isfinite(bey).all()), dt
    # _log_erfcx: matches log∘erfcx where evaluable, continuous at the branch
    z = torch.linspace(-5.9, 6.0, 1001, dtype=torch.float64)
    r = float((flow._log_erfcx(z) - torch.special.erfcx(z).log()).abs().max())
    zb = torch.tensor([-6.0 - 1e-9, -6.0 + 1e-9], dtype=torch.float64)
    jump = float((flow._log_erfcx(zb)[1] - flow._log_erfcx(zb)[0]).abs())
    print(f"_outward rising-slope finite (fp32+fp64); _log_erfcx max err {r:.1e}, "
          f"branch jump {jump:.1e}")
    assert r < 1e-12 and jump < 1e-7


if __name__ == "__main__":
    _selftest()
