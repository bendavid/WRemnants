"""Full-support, normal-base 1-D conditional flow with mixture-CDF transforms.

Motivation (J/ψ mass-fit two-stage continuity calibration)
----------------------------------------------------------
A normal-base alternative to the compact flow: support over ALL of ℝ (a standard
Gaussian base, no fixed matched tails), built from a stack of conditional
Gaussianization layers whose per-layer transform is a Gaussian- OR logistic-
mixture CDF. It keeps the analytic benefits that make the compact flow the
production default — closed-form density, closed-form CDF (so the window-Z is
exact and cheap), and an ANALYTIC Jacobian (no inner ``autograd.grad`` →
torch.compile-traceable) — but trades compact support for ℝ support with LEARNED
tails (Gaussian-mixture → Gaussian tails; logistic-mixture → heavier exponential
tails), instead of the compact flow's bolted-on C²-matched Gaussian tail.

Construction
------------
Base z ~ N(0,1). The inverse map f : x → z is a composition of ``n_transforms``
conditional Gaussianization layers:

    u_0 = x_std
    u_l = Φ⁻¹( C_l(u_{l-1} | c) ),   C_l(u|c) = Σ_k π_k(c)·K_cdf((u−μ_k(c))/s_k(c))

where ``K_cdf`` is the standard normal CDF (``layer_type='gaussian'``) or the
logistic CDF σ (``layer_type='logistic'``), and (π, μ, s) come from a per-layer
conditioner MLP of c only (a CONDITIONAL, not autoregressive, flow). Each layer
maps ℝ→ℝ (the mixture CDF maps ℝ→(0,1), then Φ⁻¹ maps back to ℝ), so the support
is all of ℝ.

Density (analytic, no autograd):
    log p(x) = log φ(z) + Σ_l log T_l'(u_{l-1}),
    T_l(u) = Φ⁻¹(C_l(u)),  log T_l'(u) = log C_l'(u) − log φ(T_l(u)),
where C_l' is the mixture pdf. (For L=1 the φ terms cancel and p = the mixture
pdf — i.e. a plain conditional mixture density; depth adds tail/multi-scale
flexibility, see the design notes.)

CDF / window-Z (exact, cheap): for ANY normal-base flow the data CDF is
    F(x|c) = Φ( f(x|c) ),
so ``transform`` returns the base-space z and the model computes the window mass
log(Φ(z_hi) − Φ(z_lo)) with its fp64-stable ``_log_phi_window`` (the same path as
the gf/nsf normal-base archs — NOT the compact ``exp−exp`` path, since the window
mass of a full-support flow can legitimately be < 1 and drift small).

Kernel choice / a note on the logistic case
-------------------------------------------
``gaussian`` (default) Gaussianization layers surject ℝ→ℝ (Φ saturates fast:
Φ(±8)≈1−1e−15), so the global density integrates to 1 exactly. ``logistic``
kernels have heavier (exponential) tails but saturate slowly (σ(±8)≈0.9997); once
``Φ⁻¹`` caps the intermediate latents at ≈±8, a STACK of logistic layers no longer
spans the full (0,1) per layer, so the composed map's range shrinks and the GLOBAL
integral is a c-dependent factor < 1. This is HARMLESS here: the model only ever
uses the WINDOW-renormalised density ``p / [Φ(z_hi) − Φ(z_lo)]``, which divides out
that factor exactly (both stage-1 training, via the flow-window-Z correction, and
the stage-2 fit window). Prefer ``gaussian`` unless you specifically want the
heavier logistic tails; a single logistic layer (L=1) also surjects exactly.

Exposes the interface the mass-fit model needs:
    forward(x_std, c)   -> log p₀(x_std | c)      (standardised-mass log-density)
    transform(x_std, c) -> z = f(x_std | c)       (base-space latent, for window-Z)
    log_cdf(x_std, c)   -> log F₀(x_std | c) = log Φ(z)
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

_SQRT2 = math.sqrt(2.0)
_LOG_SQRT_2PI = 0.5 * math.log(2.0 * math.pi)
# erfinv argument clamp (in fp64): keeps Φ⁻¹ finite while still reaching
# |z|≈8 (Φ(±8)≈6e-16) so the window mass can be resolved down to ~1e-15 —
# ample for any sane flow⊋fit window (Z is O(0.1–1) there).
_C_EPS = 1e-15


def _normal_logpdf(z: torch.Tensor) -> torch.Tensor:
    return -0.5 * z * z - _LOG_SQRT_2PI


def _normal_icdf(p: torch.Tensor) -> torch.Tensor:
    """Φ⁻¹(p) = √2·erfinv(2p−1), evaluated in fp64 with a clamp so erfinv stays
    finite (its argument never reaches ±1), then cast back to the input dtype."""
    dt = p.dtype
    pp = p.double().clamp(_C_EPS, 1.0 - _C_EPS)
    z = _SQRT2 * torch.erfinv(2.0 * pp - 1.0)
    return z.to(dt)


class MixtureFlow(nn.Module):
    """Normal-base, full-support 1-D conditional flow: a stack of Gaussian- or
    logistic-mixture-CDF Gaussianization layers with an analytic Jacobian/CDF."""

    def __init__(
        self,
        n_cond: int,
        hidden_features: int = 64,
        n_layers: int = 2,
        n_components: int = 8,
        n_transforms: int = 4,
        activation: type[nn.Module] = nn.GELU,
        s_min: float = 1e-3,
        layer_type: str = "gaussian",
    ):
        super().__init__()
        if layer_type not in ("gaussian", "logistic"):
            raise ValueError(
                f"layer_type must be 'gaussian' or 'logistic'; got {layer_type!r}")
        if int(n_components) < 1:
            raise ValueError(f"n_components must be ≥ 1; got {n_components}")
        self.layer_type = str(layer_type)
        self.K = int(n_components)
        self.L = int(n_transforms)
        self.s_min = float(s_min)
        # One INDEPENDENT conditioner MLP per layer (c → 3K = [logit_π | μ | raw_s]),
        # mirroring the compact flow / gf per-transform hyper-network. Init: tiny
        # weights (≈ c-independent), biases giving μ spread over a small range,
        # s≈1, equal weights → each layer ≈ a mild Gaussianization of the
        # standardised (μ≈0, σ≈1) mass, so the initial density ≈ N(0,1) and the
        # composition is a stable near-identity (no erfinv blow-up at init).
        out_dim = 3 * self.K
        mu_b = (torch.linspace(-1.0, 1.0, self.K) if self.K > 1
                else torch.zeros(1))
        s_b = math.log(math.expm1(max(1.0 - self.s_min, 1e-4)))   # softplus⁻¹(≈1)
        self.conditioners = nn.ModuleList()
        for _ in range(self.L):
            seq: list[nn.Module] = []
            d = int(n_cond)
            for _ in range(int(n_layers)):
                seq += [nn.Linear(d, hidden_features), activation()]
                d = hidden_features
            final = nn.Linear(d, out_dim)
            nn.init.normal_(final.weight, std=1e-3)
            with torch.no_grad():
                b = final.bias.view(3, self.K)
                b.zero_()
                b[0, :] = 0.0          # logit_π → uniform
                b[1, :] = mu_b         # μ spread
                b[2, :] = s_b          # s ≈ 1
            seq.append(final)
            self.conditioners.append(nn.Sequential(*seq))

    # ---- per-layer mixture parameters --------------------------------------
    def _cond_params(self, c: torch.Tensor):
        """Conditioner outputs → (log_π, μ, s), each [B, L, K]. π via softmax,
        s = s_min + softplus(raw) > 0 (strictly positive component scales)."""
        h = torch.stack([cond(c) for cond in self.conditioners], dim=1)  # [B,L,3K]
        h = h.view(h.shape[0], self.L, 3, self.K)
        log_pi = F.log_softmax(h[:, :, 0, :], dim=-1)
        mu = h[:, :, 1, :]
        s = self.s_min + F.softplus(h[:, :, 2, :])
        return log_pi, mu, s

    # ---- one mixture layer: C(u), log C'(u) --------------------------------
    def _layer_C_logCp(self, u, lp_l, mu_l, s_l):
        """Mixture CDF C(u) ∈ (0,1) and log of its derivative (the mixture pdf),
        for u [B] and per-layer params [B,K]. Gaussian or logistic kernels."""
        z = (u.unsqueeze(-1) - mu_l) / s_l                       # [B,K]
        log_s = s_l.log()
        if self.layer_type == "gaussian":
            cdf_k = 0.5 * (1.0 + torch.erf(z / _SQRT2))          # Φ(z)
            log_pdf_k = -0.5 * z * z - _LOG_SQRT_2PI             # log φ(z)
        else:                                                    # logistic
            cdf_k = torch.sigmoid(z)                             # σ(z)
            # log[σ(z)(1−σ(z))] = −softplus(z) − softplus(−z)
            log_pdf_k = -F.softplus(z) - F.softplus(-z)
        C = (lp_l.exp() * cdf_k).sum(-1)                         # Σ π_k CDF_k
        # log C'(u) = logsumexp_k( log π_k − log s_k + log pdf_k )  (stable)
        log_cp = torch.logsumexp(lp_l - log_s + log_pdf_k, dim=-1)
        return C, log_cp

    def _forward_zlogp(self, x: torch.Tensor, c: torch.Tensor):
        """Run the stack: returns (z, log p) where z = f(x|c) (base latent) and
        log p = log φ(z) + Σ_l log T_l'. Fully analytic (no autograd)."""
        log_pi, mu, s = self._cond_params(c)
        u = x
        logdet = torch.zeros_like(x)
        for l in range(self.L):
            C, log_cp = self._layer_C_logCp(u, log_pi[:, l], mu[:, l], s[:, l])
            u_next = _normal_icdf(C)
            # log T_l'(u) = log C'(u) − log φ(T_l(u))
            logdet = logdet + log_cp - _normal_logpdf(u_next)
            u = u_next
        z = u
        logp = _normal_logpdf(z) + logdet
        return z, logp

    # ---- public interface (mirrors compact_flow) ---------------------------
    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B]. ``x`` may be [B] or [B,1]."""
        x = x.reshape(-1)
        _, logp = self._forward_zlogp(x, c)
        return logp

    def transform(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """z = f(x_std | c), [B] — the base-space (N(0,1)) latent. The model uses
        this for the fp64-stable window mass log(Φ(z_hi) − Φ(z_lo))."""
        x = x.reshape(-1)
        z, _ = self._forward_zlogp(x, c)
        return z

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c) = log Φ(f(x_std | c)), [B]. Exact (normal-base
        CDF), via fp64 log_ndtr for tail accuracy."""
        z = self.transform(x, c)
        return torch.special.log_ndtr(z.double()).to(z.dtype)


# ---------------------------------------------------------------------------
# Self-test: normalisation, density↔CDF consistency, full-ℝ support, compile.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    xs = torch.linspace(-60, 60, 120001, dtype=torch.float64)
    for layer_type in ("gaussian", "logistic"):
        for L in (1, 3):
            m = MixtureFlow(n_cond=4, hidden_features=16, n_layers=2,
                            n_components=6, n_transforms=L, layer_type=layer_type)
            m.double()
            B = 5
            c = torch.randn(B, 4, dtype=torch.float64)
            for i in range(B):
                ci_g = c[i:i + 1].expand(xs.shape[0], -1)
                with torch.no_grad():
                    p = m.forward(xs, ci_g).exp()
                # (1) THE model-relevant property: the WINDOW-renormalised density
                #     ∫_w p / (F(w_hi) − F(w_lo)) = 1 EXACTLY for both kernels
                #     (this is what the model always uses; cancels any global
                #     sub-normalisation). Window [-1, 1.5].
                wlo, whi = -1.0, 1.5
                msk = (xs >= wlo) & (xs <= whi)
                ci1 = c[i:i + 1]
                Z = float(m.log_cdf(torch.tensor([whi], dtype=torch.float64), ci1).exp()
                          - m.log_cdf(torch.tensor([wlo], dtype=torch.float64), ci1).exp())
                win = float(torch.trapz(p[msk], xs[msk])) / Z
                assert abs(win - 1.0) < 2e-3, \
                    f"{layer_type} L={L}: window-renorm ∫={win}"
                # (2) density ↔ CDF consistency: dF/dx ≈ p
                x0 = torch.tensor([0.3], dtype=torch.float64); h = 1e-4
                dF = float((m.log_cdf(x0 + h, ci1).exp()
                            - m.log_cdf(x0 - h, ci1).exp()) / (2 * h))
                p0 = float(m.forward(x0, ci1).exp())
                assert abs(dF - p0) < 1e-3 * max(p0, 1e-6), \
                    f"{layer_type} L={L}: dF/dx={dF} vs p={p0}"
            # (1b) GAUSSIAN kernel additionally surjects ℝ→ℝ → ∫_ℝ p = 1 exactly
            #      (logistic stacks saturate → harmless global factor < 1, see
            #      module docstring; the model window-renormalises either way).
            if layer_type == "gaussian":
                tot = float(torch.trapz(m.forward(xs, c[0:1].expand(xs.shape[0], -1))
                                        .exp(), xs))
                assert abs(tot - 1.0) < 2e-3, f"gaussian L={L}: ∫_ℝ p={tot}"
            # (3) backward works (analytic graph; no inner autograd)
            mf = MixtureFlow(n_cond=4, n_components=6, n_transforms=L,
                             layer_type=layer_type)
            cf = torch.randn(8, 4, requires_grad=True)
            xf = torch.randn(8)
            mf.forward(xf, cf).sum().backward()
            assert cf.grad is not None and torch.isfinite(cf.grad).all()
            print(f"  OK {layer_type} L={L}: window-renorm ∫=1, dF/dx↔p, "
                  f"backward finite")
    # 4) torch.compile traces the analytic forward (no graph-breaking autograd)
    mc = MixtureFlow(n_cond=4, n_components=6, n_transforms=3,
                     layer_type="gaussian")
    fc = torch.compile(mc.forward, fullgraph=True)
    out = fc(torch.randn(8), torch.randn(8, 4))
    assert out.shape == (8,) and torch.isfinite(out).all()
    print("  OK torch.compile(fullgraph=True) traced the forward")
    print("mixture_flow self-test PASSED")


if __name__ == "__main__":
    _selftest()
