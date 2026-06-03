"""NCE (classifier-vs-uniform) conditional density on the mass window.

Motivation (J/ψ mass-fit two-stage continuity calibration)
----------------------------------------------------------
Noise-contrastive estimation with UNIFORM noise on the (standardised) trigger
mass window [a, b]: a plain MLP classifier is trained with binary cross-entropy
to separate simulation events (m, c) from uniform-mass twins (u, c) sharing the
SAME conditioning c. With per-event class totals balanced (prior odds 1), the
optimal logit is

    logit*(m, c) = log p₀(m | c) − log u(m),     u = 1/(b − a) on [a, b],

so ``log p₀ = logit − log(b−a)`` is an ABSOLUTE, statistically window-normalised
conditional density — a drop-in replacement for the stage-1 flow with:

  * truncation exact by construction — the classifier is trained on the
    windowed population against a window-uniform reference, so the learned
    object IS the truncated density (no full-support-vs-truncated mismatch,
    no out-of-window Z gauge drift);
  * the conditional (not joint) ratio guaranteed by PAIRING — every uniform
    twin carries its event's exact c, so the two classes have identical
    c-marginals event-by-event and the classifier cannot learn p(c);
  * an unconstrained C∞ architecture (no invertibility/monotonicity gymnastics)
    whose logit extrapolates smoothly just outside the window (the scale
    un-kick / smear / window-Z all probe ~O(MeV) beyond the edges).

What the flow had and the NCE density does not: an analytic CDF. The window
normalisation Z = F₀(hi) − F₀(lo) is supplied instead by per-event Gauss-
Legendre quadrature of the density (``log_cdf``). NOTE: ``log_cdf`` carries an
arbitrary additive offset in F (a constant +1) chosen so F stays positive for
the small out-of-window probes — it is valid ONLY for CDF *differences*
F(x₂) − F(x₁), which is the only way the mass-fit model consumes it.

Calibration: a properly trained classifier has Z₀(c) = ∫_window p̂₀ ≈ 1, but
finite capacity/statistics leave Z₀ = 1 + δ(c). δ is a per-event,
m-independent, θ-independent factor: it cancels EXACTLY from the window-
normalised stage-2 density (log p − log Z computed from the same p̂₀ — the same
gauge invariance that made the flow's Z drift harmless) and is pure gauge in
the pure-signal NLL; only the mixture fractions see it at first order. The
trainer reports the post-training Z₀ distribution as a quality flag.

Exposes the interface the mass-fit model needs:
    forward(x_std, c) -> log p₀(x_std | c)   (standardised-mass log-density)
    log_cdf(x_std, c) -> log F₀(x_std | c)   (quadrature cumulative, offset by
                                              +1; DIFFERENCES only — window Z
                                              etc.)
    nce_loss(x_std, c, n_noise) -> [B]       (per-event paired BCE, stage 1)
"""
from __future__ import annotations

import math

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class NCEDensity(nn.Module):
    """Classifier-vs-uniform (NCE) conditional density on [a, b]."""

    def __init__(
        self,
        n_cond: int,
        a: float,
        b: float,
        hidden_features: int = 128,
        n_layers: int = 3,
        activation: type[nn.Module] = nn.GELU,
        quad_nodes: int = 64,
    ):
        super().__init__()
        if not (b > a):
            raise ValueError(f"need b>a; got a={a}, b={b}")
        self.register_buffer("a", torch.tensor(float(a)))
        self.register_buffer("b", torch.tensor(float(b)))
        self._log_width = math.log(float(b) - float(a))

        # Plain MLP logit net: [x_std | c] → 1. GELU keeps it C∞ — the stage-2
        # continuity operator consumes ∂_m log p₀ and ∂²_m log p₀. Final layer
        # zero-init → logit ≡ 0 → p̂₀ = uniform on the window at init (the NCE
        # reference itself), so training starts from a calibrated classifier.
        seq: list[nn.Module] = []
        d = int(n_cond) + 1
        for _ in range(int(n_layers)):
            seq += [nn.Linear(d, hidden_features), activation()]
            d = hidden_features
        final = nn.Linear(d, 1)
        with torch.no_grad():
            final.weight.zero_()
            final.bias.zero_()
        seq.append(final)
        self.net = nn.Sequential(*seq)

        # Gauss-Legendre rule for the quadrature CDF. Non-persistent (rebuilt
        # at __init__ from quad_nodes) so checkpoints don't pin the node count.
        # Stored float64 (downcast on use): the rule must stay exact when the
        # model is evaluated in fp64 (diagnostics/scan).
        t, w = np.polynomial.legendre.leggauss(int(quad_nodes))
        self.register_buffer("_gl_t", torch.from_numpy(t), persistent=False)
        self.register_buffer("_gl_w", torch.from_numpy(w), persistent=False)

    # ---- core ---------------------------------------------------------------
    def logit(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x|c) − log u(x): the classifier logit. x [B], c [B, n_cond]."""
        return self.net(torch.cat([x.unsqueeze(-1), c], dim=-1)).squeeze(-1)

    # ---- public API (flow drop-in) ------------------------------------------
    def forward(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log p₀(x_std | c), [B]. ``x`` may be [B] or [B,1]. Outside [a, b]
        this is the smooth MLP extrapolation of the logit (the un-kick/smear/Z
        probes reach a few MeV beyond the edges)."""
        return self.logit(x.reshape(-1), c) - self._log_width

    def log_cdf(self, x: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        """log F₀(x_std | c), [B], with F₀(x) = 1 + ∫_a^x p̂₀(u|c) du by
        per-event Gauss-Legendre quadrature on [a, x] (signed for x < a; the
        nodes/weights depend on x, so autograd gradients w.r.t. x — the
        θ-dependent un-kicked boundaries — are those of the quadrature itself,
        consistent and accurate).

        The +1 offset keeps F positive for evaluation points slightly BELOW
        ``a`` (a small negative signed integral): this log-CDF is meaningful
        ONLY in differences F(x₂)−F(x₁), which is how all model consumers use
        it (window Z, stage-1 window norm, display normalisation)."""
        x = x.reshape(-1)
        B = x.shape[0]
        half = 0.5 * (x - self.a)              # [B] signed half-length of [a, x]
        mid = 0.5 * (x + self.a)
        t = self._gl_t.to(x.dtype)
        w = self._gl_w.to(x.dtype)
        u = mid.unsqueeze(-1) + half.unsqueeze(-1) * t            # [B, N]
        N = u.shape[1]
        cu = c.unsqueeze(1).expand(B, N, c.shape[-1]).reshape(B * N, -1)
        lp = self.forward(u.reshape(-1), cu).reshape(B, N)
        integral = (lp.exp() * w).sum(-1) * half                   # signed
        return (1.0 + integral).clamp_min(1e-30).log()

    # ---- stage-1 training objective ------------------------------------------
    def nce_loss(self, x: torch.Tensor, c: torch.Tensor,
                 n_noise: int = 8) -> torch.Tensor:
        """Per-event paired NCE binary cross-entropy, [B].

        Each event (x, c) is classified against ``n_noise`` uniform-mass twins
        (u_k, c) sharing its exact conditioning (identical class c-marginals
        event-by-event → the learned ratio is the CONDITIONAL p(m|c)/u, with no
        joint-p(c) leakage). Twin losses are averaged (weight 1/n_noise each)
        so the per-event class totals are balanced → prior odds 1 → the optimal
        logit is log[p₀(x|c)·(b−a)] pointwise. At zero-init (logit ≡ 0) the
        loss is exactly 2·log 2."""
        x = x.reshape(-1)
        B = x.shape[0]
        k = int(n_noise)
        u = self.a + (self.b - self.a) * torch.rand(
            B, k, device=x.device, dtype=x.dtype)
        xs = torch.cat([x.unsqueeze(-1), u], dim=1)                # [B, 1+k]
        cs = c.unsqueeze(1).expand(B, 1 + k, c.shape[-1]).reshape(-1, c.shape[-1])
        lg = self.logit(xs.reshape(-1), cs).reshape(B, 1 + k)
        loss_real = -F.logsigmoid(lg[:, 0])
        loss_noise = -F.logsigmoid(-lg[:, 1:]).mean(dim=1)
        return loss_real + loss_noise


# ---------------------------------------------------------------------------
# Self-test: init calibration, quadrature CDF↔density consistency, signed
# below-window behaviour, loss value/gradients, and a small functional
# recovery test on a known conditional density.
# ---------------------------------------------------------------------------
def _selftest():
    torch.manual_seed(0)
    a, b = -2.5, 3.0
    n_cond = 7
    net = NCEDensity(n_cond, a, b, hidden_features=64, n_layers=2,
                     quad_nodes=64).double()

    B = 5
    c = torch.randn(B, n_cond, dtype=torch.float64)
    x = torch.linspace(a + 0.1, b - 0.1, B, dtype=torch.float64)

    # 1) zero-init: p̂ = uniform exactly; F(b)−F(a) = 1 exactly (GL is exact
    #    for constants); nce_loss = 2 log 2 exactly.
    lp = net(x, c)
    assert torch.allclose(lp, torch.full_like(lp, -math.log(b - a)))
    Fb = net.log_cdf(torch.full((B,), b, dtype=torch.float64), c).exp()
    Fa = net.log_cdf(torch.full((B,), a, dtype=torch.float64), c).exp()
    print(f"zero-init: F(b)-F(a) = {(Fb - Fa).max().item():.12f} (expect 1)")
    assert torch.allclose(Fb - Fa, torch.ones_like(Fb))
    loss = net.nce_loss(x, c, n_noise=4)
    assert torch.allclose(loss, torch.full_like(loss, 2 * math.log(2.0)))
    print(f"zero-init: nce_loss = {loss[0].item():.6f} (expect {2*math.log(2):.6f})")

    # 2) perturb the net → non-trivial density; check the quadrature CDF is
    #    consistent with the density: dF/dx ≈ p (central FD) and F(b)−F(a)
    #    matches a dense trapezoid integral.
    with torch.no_grad():
        for p in net.parameters():
            p.add_(0.3 * torch.randn_like(p))
    eps = 1e-5
    Fp = net.log_cdf(x + eps, c).exp()
    Fm = net.log_cdf(x - eps, c).exp()
    dFdx = (Fp - Fm) / (2 * eps)
    p_x = net(x, c).exp()
    rel = ((dFdx - p_x).abs() / p_x).max().item()
    print(f"perturbed: max rel diff dF/dx vs p = {rel:.2e}")
    assert rel < 1e-6
    grid = torch.linspace(a, b, 20001, dtype=torch.float64)
    for i in range(B):
        ci = c[i].expand(grid.shape[0], -1)
        pg = net(grid, ci).exp()
        z_trap = torch.trapezoid(pg, grid).item()
        z_quad = (net.log_cdf(torch.tensor([b], dtype=torch.float64), c[i:i+1]).exp()
                  - net.log_cdf(torch.tensor([a], dtype=torch.float64), c[i:i+1]).exp()).item()
        assert abs(z_trap - z_quad) / max(abs(z_trap), 1.0) < 1e-6, (z_trap, z_quad)
    print(f"perturbed: quadrature Z matches trapezoid to <1e-6 ({B} events)")

    # 3) signed behaviour just below/above the window (the un-kick probes).
    xb = torch.tensor([a - 0.05, a, b, b + 0.05], dtype=torch.float64)
    cb = c[:1].expand(4, -1)
    Fv = net.log_cdf(xb, cb).exp()
    assert Fv[0] < Fv[1] < Fv[2] < Fv[3]
    print(f"out-of-window probes monotone: {[f'{v:.6f}' for v in Fv.tolist()]}")

    # 4) gradients reach every parameter through both the loss and log_cdf.
    loss = net.nce_loss(x, c, n_noise=4).sum()
    loss.backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all()
               for p in net.parameters())
    net.zero_grad()
    net.log_cdf(x, c).sum().backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all()
               for p in net.parameters())
    net.zero_grad()
    print("gradients finite through nce_loss and log_cdf")

    # 5) functional recovery: 1-cond truncated Gaussian, mean depends on c.
    #    Train briefly; check the learned density tracks the truth and Z₀≈1.
    torch.manual_seed(1)
    net2 = NCEDensity(1, a, b, hidden_features=64, n_layers=2, quad_nodes=64)
    opt = torch.optim.Adam(net2.parameters(), lr=3e-3)
    sigma = 0.5
    for step in range(600):
        cc = torch.rand(4096, 1) * 2 - 1                  # c ∈ [−1, 1]
        mu = 0.8 * cc[:, 0]
        xx = (mu + sigma * torch.randn(4096)).clamp(a + 1e-3, b - 1e-3)
        opt.zero_grad()
        net2.nce_loss(xx, cc, n_noise=8).mean().backward()
        opt.step()
    with torch.no_grad():
        cc = torch.tensor([[-0.5], [0.0], [0.5]])
        gx = torch.linspace(a, b, 401)
        for i in range(3):
            mu = 0.8 * cc[i, 0]
            zt = 0.5 * (math.erf((b - mu) / (sigma * math.sqrt(2)))
                        - math.erf((a - mu) / (sigma * math.sqrt(2))))
            pt = (torch.exp(-0.5 * ((gx - mu) / sigma) ** 2)
                  / (sigma * math.sqrt(2 * math.pi)) / zt)
            ph = net2(gx, cc[i].expand(gx.shape[0], -1)).exp()
            core = pt > 0.05 * pt.max()
            relmax = ((ph - pt).abs() / pt)[core].max().item()
            z0 = (net2.log_cdf(torch.tensor([float(b)]), cc[i:i+1]).exp()
                  - net2.log_cdf(torch.tensor([float(a)]), cc[i:i+1]).exp()).item()
            print(f"recovery c={cc[i,0]:+.1f}: max rel err (core) = {relmax:.3f}, "
                  f"Z0 = {z0:.4f}")
            assert relmax < 0.10, relmax
            assert abs(z0 - 1.0) < 0.05, z0
    print("functional recovery OK (≤10% core density error, |Z0−1| < 5%)")
    print("ALL NCE SELF-TESTS PASSED")


if __name__ == "__main__":
    _selftest()
