"""Diagnostic plots for a trained muon-response normalizing flow.

Loads a ``checkpoint.pt`` (or a fully-trained ``flow.pt``) written by
``train_muon_response_flow.py``, evaluates the flow on a held-out MC
sample (typically the same snapshot used for training), and produces
a small set of diagnostic plots:

  marginals.png             MC vs flow 1D histograms for r_kappa,
                            dlambda, dphi (unconditioned), with ratio.
  slices_log_pt.png         Same three marginals, but stratified into
                            quartile bins of log(pt_gen).
  slices_lambda.png         Same, stratified into quartile bins of
                            lambda_gen.
  slices_charge.png         Same, stratified by muon charge (q = +1 vs -1).
  log_p_per_event.png       Distribution of log p(x|c) across events
                            (sanity: finite, unimodal, reasonable scale).

All plots are produced for a single checkpoint and a single validation
sample. Use multiple invocations (different --checkpoint) to compare
snapshots across epochs.
"""

import argparse
import os
import sys
from dataclasses import asdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

# Import helpers from the sibling training script.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from train_muon_response_flow import (  # noqa: E402
    ACTIVATIONS,
    PreprocStats,
    _sigma_pack_outer,
    apply_preproc,
    build_flow,
    compute_targets_and_conditioning,
    load_ntuples,
)


TARGET_NAMES = ["r_kappa", "dlambda", "dphi"]


def parse_args():
    # Combined formatter: keep the raw module-docstring description
    # (preserves newlines / paragraphs) AND append ``(default: ...)``
    # to every option's help text.
    class _Fmt(
        argparse.RawDescriptionHelpFormatter,
        argparse.ArgumentDefaultsHelpFormatter,
    ):
        pass

    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_Fmt,
    )
    p.add_argument(
        "--checkpoint",
        required=True, default=argparse.SUPPRESS,
        help="(required) Path to checkpoint.pt (or flow.pt) from the "
        "training run.",
    )
    p.add_argument(
        "--input-files",
        nargs="+",
        required=True, default=argparse.SUPPRESS,
        help="(required) Snapshot ROOT file(s) (from "
        "flow_training_snapshot.py).",
    )
    p.add_argument(
        "--tree", default="tree",
        help="ROOT TTree name in --input-files.",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="Directory for plots. None = directory containing "
        "--checkpoint.",
    )
    p.add_argument(
        "--n-events",
        type=int,
        default=500_000,
        help="Number of *muon rows* to use for diagnostics. Applied "
        "post-filter: subsamples the rows that survive the quality "
        "cuts. Independent of --max-events.",
    )
    p.add_argument(
        "--max-events",
        type=int,
        default=-1,
        help="Cap on the number of raw J/psi events kept via an "
        "RDataFrame ``Filter(\"rdfentry_ < N\")`` — keeps the first "
        "N rows of the underlying TTree. ImplicitMT-compatible. "
        "Reduces post-filter pipeline cost (quality cut, "
        "materialization, downstream ops) but not raw disk I/O. "
        "Default -1 keeps all events.",
    )
    p.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device for flow evaluation + sampling.",
    )
    p.add_argument(
        "--batch-size", type=int, default=65536,
        help="Batch size for flow forward / log-prob evaluation.",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="RDataFrame ImplicitMT threads for loading. 0 = ROOT auto.",
    )
    p.add_argument(
        "--pt-min", type=float, default=2.0,
        help="Min gen pt (GeV) per muon; matches snapshot script.",
    )
    p.add_argument(
        "--pt-max", type=float, default=200.0,
        help="Max gen pt (GeV) per muon.",
    )
    p.add_argument(
        "--eta-max", type=float, default=2.5,
        help="Max |gen eta| per muon.",
    )
    p.add_argument(
        "--n-bins",
        type=int,
        default=80,
        help="Histogram bin count for marginal plots.",
    )
    p.add_argument(
        "--range-percentile",
        type=float,
        default=0.05,
        help="Percentile (both tails) used to set histogram x-range on each target.",
    )
    p.add_argument(
        "--shift-factors",
        nargs="+",
        type=float,
        default=[0.05, 0.1, 0.2],
        help="Per-component shift magnitudes (in units of target_std) "
        "for the shift-reweight closure plot.",
    )
    p.add_argument(
        "--smear-factors",
        nargs="+",
        type=float,
        default=[0.1, 0.25, 0.5],
        help="Per-component smear magnitudes (in units of target_std) "
        "for the smear-reweight closure plot.",
    )
    p.add_argument(
        "--reweight-chunk",
        type=int,
        default=4096,
        help="Chunk size for batched score + diag-Hessian evaluation "
        "during reweight closure tests. Lower if memory-constrained.",
    )
    p.add_argument(
        "--skip-reweight",
        action="store_true",
        help="Skip the shift/smear reweight closure plots.",
    )
    p.add_argument(
        "--skip-polyhead-plots",
        action="store_true",
        help="Skip the polyhead validation plots even if a polyhead "
        "is present in the checkpoint.",
    )
    p.add_argument(
        "--polyhead-validate-n",
        type=int,
        default=10000,
        help="Number of MC events used for polyhead validation "
        "(scatter + log-error histogram). Smaller = faster.",
    )
    p.add_argument(
        "--polyhead-shift-deltas",
        nargs="+",
        type=float,
        default=[0.1, 0.3, 0.5, 1.0],
        help="Shift magnitudes |δ| (in standardized cond units) for "
        "the polyhead reweight closure. Each is applied along each "
        "of the n_cond unit-basis directions in turn. Matches the "
        "direct shift-reweight diagnostic's --shift-factors default "
        "for side-by-side comparability.",
    )
    p.add_argument(
        "--polyhead-smear-sigmas",
        nargs="+",
        type=float,
        default=[0.1, 0.3, 0.5, 1.0],
        help="Smear scales σ (in standardized cond units) for the "
        "polyhead smear-reweight closure. Matches the direct "
        "shift-reweight diagnostic's --smear-factors default.",
    )
    p.add_argument(
        "--polyhead-closure-percentile",
        type=float,
        default=0.5,
        help="Percentile (each tail) for setting the histogram x-range "
        "in the polyhead shift / smear closure plots. Defaults to 0.5 "
        "to match the direct shift-reweight diagnostic's "
        "--range-percentile default (so the two scripts produce "
        "side-by-side-comparable closure plots). The tighter "
        "0.05-percentile default of --range-percentile is kept for "
        "the marginal / slice plots which need wider tails.",
    )
    return p.parse_args()


# -----------------------------------------------------------------------------
# Checkpoint loading
# -----------------------------------------------------------------------------

def load_polyhead_from_checkpoint(ckpt, n_features, n_cond, device):
    """Build the reweight head from ``ckpt`` and load its weights;
    return ``None`` if no head is present.

    Both ``--head-arch polyhead`` and ``--head-arch mlp`` are
    supported — the saved ``head_config["head_arch"]`` selects which
    nn.Module to construct. The diagnostic plot suite consumes the
    head's per-event predicted W via :func:`_polyhead_pred_W_at`,
    which dispatches internally on isinstance, so the same plots
    work for both arches.

    For MLP heads, a few attributes the plot code inspects via
    ``getattr(head, ...)`` are attached after construction:
    ``positivity``, ``basis_scale_u``/``basis_scale_sigma`` (used by
    the plots only as the perturbation-magnitude bound for the
    diagnostic-time u/σ sampling — same convention as the polyhead's
    ``oversample · delta_max`` / ``oversample · sigma_max``), and the
    ``sigma_pack_iu/ju`` buffers needed to build σ_pack on the fly
    inside the dual-forward predictor.
    """
    if "polyhead_state_dict" not in ckpt:
        return None
    # ``head_config`` is the new name; older checkpoints used
    # ``polyhead_config``. Try both.
    cfg = ckpt.get("head_config") or ckpt.get("polyhead_config", {})
    head_arch = str(cfg.get("head_arch", "polyhead")).lower()
    activation_cls = ACTIVATIONS.get(
        str(cfg.get("activation", "gelu")).lower(),
        torch.nn.GELU,
    )

    if head_arch == "polyhead":
        from flow_polyhead import PolyHead
        # ``include_smear`` was introduced with the shift-only-by-
        # default polyhead. Older checkpoints pre-date it and were
        # always full (smear-included), so default to True for
        # backwards compat — fresh checkpoints always carry the
        # explicit value.
        include_smear = bool(cfg.get("include_smear", True))
        head = PolyHead(
            n_features=n_features,
            n_cond=n_cond,
            hidden_features=int(
                cfg.get("trunk_hidden", cfg.get("hidden_features", 64)),
            ),
            n_hidden_layers=int(
                cfg.get("trunk_layers", cfg.get("n_hidden_layers", 2)),
            ),
            max_deg_u=int(cfg.get("max_deg_u", cfg.get("max_deg", 3))),
            max_deg_sigma=int(cfg.get("max_deg_sigma", 4)),
            max_cross_deg=int(cfg.get("max_cross_deg", 3)),
            activation=activation_cls,
            smear_K=int(cfg.get("smear_K", 1)),
            smear_residual=bool(cfg.get("smear_residual", False)),
            include_smear=include_smear,
            positivity=str(cfg.get("positivity", "softplus")),
            # ``basis`` was introduced with the Chebyshev option;
            # older checkpoints pre-date it and used monomial.
            basis=str(cfg.get("basis", "monomial")),
            basis_scale_u=float(cfg.get("basis_scale_u", 1.0)),
            basis_scale_sigma=float(cfg.get("basis_scale_sigma", 1.0)),
        )
    elif head_arch == "mlp":
        # Lazy import to avoid the module-level cycle between
        # train_shift_smear_reweight ↔ train_muon_response_flow.
        from train_shift_smear_reweight import (
            ReweightMLP_B, _sigma_pack_indices,
        )
        include_smear = bool(cfg.get("include_smear", True))
        head = ReweightMLP_B(
            n_features=n_features,
            n_cond=n_cond,
            d_emb=int(cfg.get("d_emb", 32)),
            trunk_hidden=int(
                cfg.get("trunk_hidden", cfg.get("hidden_features", 64)),
            ),
            trunk_layers=int(
                cfg.get("trunk_layers", cfg.get("n_hidden_layers", 2)),
            ),
            head_hidden=int(cfg.get("head_hidden", 32)),
            head_layers=int(cfg.get("head_layers", 2)),
            activation=activation_cls,
            shift_only=not include_smear,
        )
        # Attach the metadata the polyhead plot suite reads via
        # ``getattr``: positivity wrap + perturbation magnitude
        # bounds (oversample·delta_max / oversample·sigma_max from
        # the head's training).
        head.positivity = str(cfg.get("positivity", "softplus"))
        head.basis_scale_u = float(cfg.get("basis_scale_u", 1.0))
        head.basis_scale_sigma = float(cfg.get("basis_scale_sigma", 1.0))
        head.include_smear = include_smear
        # σ_pack indices for the dual-forward predictor.
        iu, ju = _sigma_pack_indices(int(n_features))
        head.register_buffer(
            "sigma_pack_iu", iu.to(device), persistent=False,
        )
        head.register_buffer(
            "sigma_pack_ju", ju.to(device), persistent=False,
        )
    else:
        raise SystemExit(
            f"unknown head_arch {head_arch!r} in checkpoint config"
        )

    head.load_state_dict(ckpt["polyhead_state_dict"])
    return head.to(device).eval()


def load_flow_from_checkpoint(path, device):
    """Load a zuko flow + PreprocStats from a checkpoint produced by
    train_muon_response_flow.py. Handles both the per-epoch
    ``checkpoint.pt`` layout (``state_dict`` keyed on the inner flow)
    and the final ``flow.pt`` layout (``wrapper_state_dict`` keyed on
    the FlowWrapper). Architecture is inferred from ``flow_config``."""
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    flow_config = ckpt["flow_config"]
    stats = PreprocStats(
        **(ckpt["stats"] if "stats" in ckpt else ckpt["preproc"])
    )

    flow_type = flow_config.get("flow_type", "RealNVP")
    architecture = flow_config.get(
        "architecture",
        {
            "RealNVP": "realnvp", "Glow": "glow", "MAF": "maf",
            "GF": "gf", "SOSPF": "sospf",
        }.get(flow_type, flow_type.lower()),
    )
    _supported = ("realnvp", "glow", "maf", "gf", "sospf")
    if architecture not in _supported:
        raise ValueError(
            f"unsupported architecture '{architecture}' (flow_type "
            f"'{flow_type}') — supported: {_supported}."
        )

    flow = build_flow(
        n_features=flow_config["n_features"],
        n_cond=flow_config["n_cond"],
        n_transforms=flow_config["n_transforms"],
        hidden_features=flow_config["hidden_features"],
        n_hidden_layers=flow_config["n_hidden_layers"],
        activation=flow_config.get("activation", "gelu"),
        architecture=architecture,
        randmask=flow_config.get("randmask", False),
        gf_components=flow_config.get("gf_components", 8),
        sospf_degree=flow_config.get("sospf_degree", 4),
        sospf_polynomials=flow_config.get("sospf_polynomials", 3),
        sospf_quad_n=flow_config.get("sospf_quad_n", -1),
    )

    if "state_dict" in ckpt:
        raw = ckpt["state_dict"]
    else:
        raw = ckpt["wrapper_state_dict"]

    # Strip any "flow." prefix introduced by the FlowWrapper.
    if any(k.startswith("flow.") for k in raw.keys()):
        raw = {
            k[len("flow."):]: v
            for k, v in raw.items()
            if k.startswith("flow.")
        }
    # Compat: zuko has used both ``base.loc/base.scale`` (older) and
    # ``base._0/base._1`` (newer). A checkpoint produced by one
    # version may not load directly into a flow built with the other,
    # so remap bidirectionally based on what the local zuko expects.
    from train_muon_response_flow import _zuko_base_remap
    raw = _zuko_base_remap(raw, flow)
    flow.load_state_dict(raw)
    flow = flow.to(device).eval()

    return flow, stats, ckpt


# -----------------------------------------------------------------------------
# Evaluation: batched log_prob + sample
# -----------------------------------------------------------------------------

@torch.no_grad()
def evaluate_flow(flow, target_std_t, cond_t, stats, device, batch_size):
    """Return (log_p [N], samples_raw [N, 3]).

    log_p is the flow log-density evaluated at each event's
    standardized target. samples_raw is one flow sample per event
    (conditioned on the event's c), converted back to raw coordinates.
    """
    n = target_std_t.shape[0]
    target_mean_t = torch.tensor(
        stats.target_mean, dtype=torch.float32, device=device
    )
    target_std_vec_t = torch.tensor(
        stats.target_std, dtype=torch.float32, device=device
    )

    log_ps = np.empty(n, dtype=np.float32)
    samples_raw = np.empty((n, target_std_t.shape[1]), dtype=np.float32)

    for s in range(0, n, batch_size):
        e = min(s + batch_size, n)
        xs = target_std_t[s:e].to(device, non_blocking=True)
        cs = cond_t[s:e].to(device, non_blocking=True)
        dist = flow(cs)
        log_p = dist.log_prob(xs)
        z = dist.sample(torch.Size([]))
        x_raw = z * target_std_vec_t + target_mean_t
        log_ps[s:e] = log_p.cpu().numpy()
        samples_raw[s:e] = x_raw.cpu().numpy()

    # Jacobian correction so log p is in raw-space: +/- log(std) sum.
    log_ps = log_ps - float(np.log(np.asarray(stats.target_std)).sum())
    return log_ps, samples_raw


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def _step(ax, edges, y, **kwargs):
    y_plot = np.r_[y, y[-1]]
    ax.step(edges, y_plot, where="post", **kwargs)


def _hist_range(mc, pct):
    lo = np.percentile(mc, pct)
    hi = np.percentile(mc, 100.0 - pct)
    return lo, hi


def _safe_ratio(numer, denom, rel_floor=1e-3):
    """Element-wise ``numer / denom`` with small-denom bins masked to NaN."""
    numer = np.asarray(numer, dtype=np.float64)
    denom = np.asarray(denom, dtype=np.float64)
    thr = rel_floor * np.nanmax(denom) if denom.size else 0.0
    safe = np.where(denom > thr, denom, 1.0)
    return np.where(denom > thr, numer / safe, np.nan)


def _safe_weights(w, factor, label=None):
    """Return ``w * factor`` with non-finite entries replaced by ``w``.

    Extreme-tail events in the diagnostic can produce ``inf`` /
    ``NaN`` reweight factors (autograd-noisy score/Hessian hitting
    pathological values in a resummed exp formula, or strict Taylor
    going to ±inf on overflow). Leaving them in the weight array
    poisons every bin in ``np.histogram(..., density=True)``, so
    the dashed line can disappear entirely. Fall back to the
    nominal per-event weight (no reweight) for those events.
    ``label`` prints a one-line diagnostic noting what fraction of
    events were substituted.
    """
    wf = w * factor
    bad = ~np.isfinite(wf)
    if bad.any():
        frac = float(bad.mean())
        if label is not None:
            print(
                f"  {label}: {frac:.2%} of events had non-finite "
                f"reweight factor → substituted nominal weight"
            )
        wf = np.where(bad, w, wf)
    return wf


def _save_png_pdf(fig, png_path):
    """Save ``fig`` as both PNG (inspection) and PDF (publication),
    then close it. ``png_path`` should end in ``.png``; the PDF is
    written alongside with the same basename.

    Closing here avoids the "20 open figures" warning when the
    diagnostic produces many plots in one run — pyplot retains
    figures by default until explicitly closed.
    """
    base, ext = os.path.splitext(png_path)
    pdf_path = base + ".pdf"
    fig.savefig(png_path, dpi=110, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def _weighted_hist_err(values, bins, weights):
    """Weighted histogram with per-bin sqrt(sum w^2) statistical
    uncertainty. Returns ``(h, sigma)``, both shape ``(n_bins,)``.
    Mirrors the helper in shift_smear_reweight_diagnostics.py for
    cross-diagnostic consistency.
    """
    w = np.asarray(weights, dtype=np.float64)
    h, _ = np.histogram(values, bins=bins, weights=w)
    h2, _ = np.histogram(values, bins=bins, weights=w * w)
    sigma = np.sqrt(np.maximum(h2, 0.0))
    return h, sigma


def _ratio_with_err(h_num, sigma_num, h_den, sigma_den):
    """Bin-wise ratio with propagated uncertainty assuming uncorrelated
    bins. Holds for the closure ratios here (h_pred / h_ref): ``h_pred``
    is filled at ``y_orig`` with reweight, ``h_ref`` at ``y_orig + dy``
    so per-bin counts are statistically independent (modulo bin-edge
    effects, neglected). Returns ``(ratio, sigma_ratio)`` with NaN
    where ``h_den ≤ 0``.
    """
    h_num = np.asarray(h_num, dtype=np.float64)
    h_den = np.asarray(h_den, dtype=np.float64)
    den_safe = np.where(h_den > 0, h_den, np.nan)
    num_safe = np.where(h_num > 0, h_num, np.nan)
    ratio = h_num / den_safe
    rel = np.sqrt(
        (sigma_num / num_safe) ** 2 + (sigma_den / den_safe) ** 2
    )
    return ratio, ratio * rel


def _stepped_errorbar(ax, centers, h, sigma, color, **kwargs):
    """Vertical sqrt(Σw²) error bars at each bin center, matched to
    a step histogram. ``kwargs`` are forwarded to errorbar."""
    ax.errorbar(
        centers, h, yerr=sigma, fmt="none", ecolor=color,
        elinewidth=0.7, capsize=0, alpha=0.6, **kwargs,
    )


def _flow_log_prob_batched(
    flow, target_std_t, cond_std_t, device, chunk,
):
    """Return ``flow(c).log_prob(x)`` in standardized space, batched.

    The Jacobian correction from standardized→raw cancels in lp
    differences, so we keep values in standardized space.
    """
    N = target_std_t.shape[0]
    out = torch.empty(N, dtype=torch.float32)
    with torch.no_grad():
        for s in range(0, N, chunk):
            e = min(s + chunk, N)
            x = target_std_t[s:e].to(device, non_blocking=True)
            c = cond_std_t[s:e].to(device, non_blocking=True)
            out[s:e] = flow(c).log_prob(x).detach().cpu()
    return out


def _flow_score_and_hessian_diag(
    flow, target_std_t, cond_std_t, device, chunk,
):
    """Return (score_std [N, d], H_diag_std [N, d]) in standardized space.

    ``score_std[:, i] = ∂ log p_std / ∂ x_std_i``; ``H_diag_std[:, i]
    = ∂² log p_std / ∂ x_std_i²``. Computed per-chunk via autograd:
    one first backward gives the full score vector; d second-backward
    calls give the diagonal of the Hessian.
    """
    N, d = target_std_t.shape
    score_out = torch.empty(N, d, dtype=torch.float32)
    hess_out = torch.empty(N, d, dtype=torch.float32)
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        x = (
            target_std_t[s:e]
            .to(device, non_blocking=True)
            .detach()
            .requires_grad_(True)
        )
        c = cond_std_t[s:e].to(device, non_blocking=True)
        log_p = flow(c).log_prob(x)
        (grad_x,) = torch.autograd.grad(
            log_p.sum(), x, create_graph=True,
        )
        score_out[s:e] = grad_x.detach().cpu()
        for i in range(d):
            (row_i,) = torch.autograd.grad(
                grad_x[..., i].sum(), x,
                create_graph=False,
                retain_graph=(i < d - 1),
            )
            hess_out[s:e, i] = row_i[:, i].detach().cpu()
    return score_out.numpy(), hess_out.numpy()


def _flow_base_quantities(
    flow, target_std_t, cond_std_t, device, chunk,
):
    """Return ``(z, J, G, K)`` in standardized space, per event.

    ``z[N, d]``: base-space latent ``z = transform(y_std; c_std)``,
    with the flow's ``DiagNormal(0, 1)`` base.
    ``J[N, d, d]``: inverse Jacobian ``J[n, j, i] = ∂z_j / ∂y_std_i``.
    ``G[N, d]``: ``G[n, i] = ∂L / ∂y_std_i`` with
    ``L = log|det ∂z/∂y_std|``.
    ``K[N, d]``: ``K[n, i] = ∂² L / ∂y_std_i²``.

    Building blocks of the base-space factorized reweight: since
    ``log p = log φ(z) + L``, shift/smear ratios factor as
    ``exp(Δ log φ)`` — closed form in ``z`` under a linearized
    ``Δz = J · δ`` — times ``exp(ΔL)`` expanded as a polynomial in
    ``δ`` via ``(G, K)``. Computed per-chunk via autograd on
    ``flow(c).transform.call_and_ladj(y)``: ``d`` reverse-mode
    backwards give the full Jacobian; one first-backward and ``d``
    second-backwards give ``(G, K)`` in the same style as
    ``_flow_score_and_hessian_diag``.
    """
    N, d = target_std_t.shape
    z_out = torch.empty(N, d, dtype=torch.float32)
    J_out = torch.empty(N, d, d, dtype=torch.float32)
    G_out = torch.empty(N, d, dtype=torch.float32)
    K_out = torch.empty(N, d, dtype=torch.float32)
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        x = (
            target_std_t[s:e]
            .to(device, non_blocking=True)
            .detach()
            .requires_grad_(True)
        )
        c = cond_std_t[s:e].to(device, non_blocking=True)
        z, ladj = flow(c).transform.call_and_ladj(x)
        # Full inverse Jacobian: d reverse-mode passes on z[:, j].
        for j in range(d):
            (row_j,) = torch.autograd.grad(
                z[:, j].sum(), x,
                create_graph=False,
                retain_graph=True,
            )
            J_out[s:e, j, :] = row_j.detach().cpu()
        # Gradient + diagonal Hessian of ladj wrt y.
        (grad_ladj,) = torch.autograd.grad(
            ladj.sum(), x,
            create_graph=True,
            retain_graph=True,
        )
        G_out[s:e] = grad_ladj.detach().cpu()
        for i in range(d):
            (row_i,) = torch.autograd.grad(
                grad_ladj[..., i].sum(), x,
                create_graph=False,
                retain_graph=(i < d - 1),
            )
            K_out[s:e, i] = row_i[:, i].detach().cpu()
        z_out[s:e] = z.detach().cpu()
    return (
        z_out.numpy(), J_out.numpy(),
        G_out.numpy(), K_out.numpy(),
    )


def logp_derivatives_on_grid(flow, stats, cond_ref_raw, grid_ranges,
                             n_grid, device, n_sigma_pack=0):
    """Evaluate ``log p(x)``, ``∂log p/∂x``, ``∂²log p/∂x²`` along
    1D slices of each target axis.

    Same slicing convention as ``pdf_derivatives_on_grid``: for each
    target dimension we vary the corresponding component across its
    ``grid_ranges`` window (raw units) while the other components
    sit at their raw mean and the conditioning vector is fixed at
    ``cond_ref_raw``. Derivatives come from autograd on the flow's
    log-density.

    These are exactly the ``(s, H_ii)`` quantities that enter the
    reweight formulas (``exp(−δ·s + ½δ²·H_ii)`` for shifts and the
    Gaussian-MGF analog for smears). Plotting them across each axis
    gives a direct read on how well-behaved the flow's first and
    second log-derivatives are — spikes or sign changes in ``H`` in
    the tails are the same pathologies that drive ``_safe_weights``
    substitutions in the shift/smear closures.

    Returns ``{target_name: (x_grid_raw, log_p, s_i, H_ii)}`` where
    all quantities are in raw coordinates.
    """
    cond_mean = np.asarray(stats.cond_mean, dtype=np.float32)
    cond_std = np.asarray(stats.cond_std, dtype=np.float32)
    cond_std_vec_np = (
        (cond_ref_raw.astype(np.float32) - cond_mean) / cond_std
    )
    cond_std_t = torch.from_numpy(cond_std_vec_np).to(device)
    if n_sigma_pack > 0:
        # Smear-conditioned flow expects (c, σ_pack); evaluate at the
        # σ=0 (unsmeared baseline) slice for diagnostics.
        cond_std_t = torch.cat(
            [
                cond_std_t,
                torch.zeros(
                    n_sigma_pack, device=device, dtype=cond_std_t.dtype,
                ),
            ],
            dim=0,
        )

    target_mean = np.asarray(stats.target_mean, dtype=np.float32)
    target_std = np.asarray(stats.target_std, dtype=np.float32)
    log_jac_raw = float(np.log(target_std).sum())

    results = {}
    for dim in range(len(TARGET_NAMES)):
        lo, hi = grid_ranges[dim]
        x_grid_raw = np.linspace(lo, hi, n_grid, dtype=np.float32)
        x_std_np = np.zeros((n_grid, len(TARGET_NAMES)), dtype=np.float32)
        x_std_np[:, dim] = (x_grid_raw - target_mean[dim]) / target_std[dim]

        x_t = torch.from_numpy(x_std_np).to(device).requires_grad_(True)
        c_t = cond_std_t.unsqueeze(0).expand(n_grid, -1)

        log_p_std = flow(c_t).log_prob(x_t)                       # (n_grid,)
        grad = torch.autograd.grad(
            log_p_std.sum(), x_t, create_graph=True,
        )[0]                                                      # (n_grid, 3)
        score_i = grad[:, dim]
        d2 = torch.autograd.grad(score_i.sum(), x_t)[0]           # (n_grid, 3)
        hess_i = d2[:, dim]

        # Raw-space conversions. log p_raw differs from log p_std by
        # the (constant in x) Jacobian term; derivatives of log p
        # rescale by 1/target_std (first) and 1/target_std² (second).
        log_p_raw = log_p_std.detach().cpu().numpy() - log_jac_raw
        s_raw = score_i.detach().cpu().numpy() / target_std[dim]
        H_raw = hess_i.detach().cpu().numpy() / (target_std[dim] ** 2)

        results[TARGET_NAMES[dim]] = (x_grid_raw, log_p_raw, s_raw, H_raw)

    return results


def plot_logp_derivatives(results, outpath, cond_ref_raw):
    """Three-row panel plot of ``log p(x | c_ref)``, ``s_i = ∂log p/∂x``,
    and ``H_ii = ∂²log p/∂x²`` along each target axis.

    All three rows use linear y-scale (``log p`` is already a log,
    ``s`` and ``H`` are signed so log-scale is meaningless).
    """
    n_dim = len(TARGET_NAMES)
    fig, axes = plt.subplots(3, n_dim, figsize=(4.2 * n_dim, 9.5),
                             sharex="col")
    for j, name in enumerate(TARGET_NAMES):
        x_grid, log_p, s, H = results[name]

        ax = axes[0, j]
        ax.plot(x_grid, log_p, color="C0")
        ax.set_title(name)
        ax.set_ylabel("log p(x | c_ref)")
        ax.grid(alpha=0.3)

        ax = axes[1, j]
        ax.plot(x_grid, s, color="C2")
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylabel("s = ∂log p/∂x")
        ax.grid(alpha=0.3)

        ax = axes[2, j]
        ax.plot(x_grid, H, color="C3")
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylabel("H = ∂²log p/∂x²")
        ax.set_xlabel(name)
        ax.grid(alpha=0.3)

    cond_str = " ".join(f"{n}={v:+.3g}" for n, v in
                        zip(["logpt", "q", "lam", "sin_phi", "cos_phi"],
                            cond_ref_raw))
    fig.suptitle(f"log p + derivatives at reference c_ref: {cond_str}",
                 y=1.01, fontsize=10)
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


def pdf_derivatives_on_grid(flow, stats, cond_ref_raw, grid_ranges,
                            n_grid, device, n_sigma_pack=0):
    """Evaluate p(x), dp/dx, d²p/dx² along 1D slices of the target.

    For each target dimension ``d``, we vary ``x_d`` over a uniform
    grid inside ``grid_ranges[d]`` (raw units) while holding the other
    target components fixed at their raw mean. The conditioning
    vector is fixed at ``cond_ref_raw`` (raw-space 5-vector) for the
    entire slice. Derivatives are taken along the varied dimension
    via autograd on the log-density, then converted from standardized
    to raw space using the stored preprocessing.

    Returns ``{target_name: (x_grid_raw, p, dp_dx, d2p_dx2)}``.
    """
    cond_mean = np.asarray(stats.cond_mean, dtype=np.float32)
    cond_std = np.asarray(stats.cond_std, dtype=np.float32)
    cond_std_vec_np = (
        (cond_ref_raw.astype(np.float32) - cond_mean) / cond_std
    )
    cond_std_t = torch.from_numpy(cond_std_vec_np).to(device)
    if n_sigma_pack > 0:
        cond_std_t = torch.cat(
            [
                cond_std_t,
                torch.zeros(
                    n_sigma_pack, device=device, dtype=cond_std_t.dtype,
                ),
            ],
            dim=0,
        )

    target_mean = np.asarray(stats.target_mean, dtype=np.float32)
    target_std = np.asarray(stats.target_std, dtype=np.float32)
    log_jac_raw = float(np.log(target_std).sum())

    results = {}
    for dim in range(len(TARGET_NAMES)):
        lo, hi = grid_ranges[dim]
        x_grid_raw = np.linspace(lo, hi, n_grid, dtype=np.float32)
        x_std_np = np.zeros((n_grid, len(TARGET_NAMES)), dtype=np.float32)
        x_std_np[:, dim] = (x_grid_raw - target_mean[dim]) / target_std[dim]

        x_t = torch.from_numpy(x_std_np).to(device).requires_grad_(True)
        c_t = cond_std_t.unsqueeze(0).expand(n_grid, -1)

        log_p_std = flow(c_t).log_prob(x_t)                       # (n_grid,)
        grad = torch.autograd.grad(
            log_p_std.sum(), x_t, create_graph=True
        )[0]                                                      # (n_grid, 3)
        score_i = grad[:, dim]
        d2 = torch.autograd.grad(score_i.sum(), x_t)[0]           # (n_grid, 3)
        hess_i = d2[:, dim]

        # Raw-space conversions.
        log_p_raw = (
            log_p_std.detach().cpu().numpy() - log_jac_raw
        )
        p_raw = np.exp(log_p_raw)
        score_raw = score_i.detach().cpu().numpy() / target_std[dim]
        hess_raw = hess_i.detach().cpu().numpy() / (target_std[dim] ** 2)

        dp_dx = p_raw * score_raw
        d2p_dx2 = p_raw * (score_raw ** 2 + hess_raw)

        results[TARGET_NAMES[dim]] = (x_grid_raw, p_raw, dp_dx, d2p_dx2)

    return results


def plot_marginals(
    target, samples_raw, outpath, n_bins=80, pct=0.05, yscale="log"
):
    fig, axes = plt.subplots(
        2, 3, figsize=(14, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex="col",
    )
    for i, name in enumerate(TARGET_NAMES):
        mc = target[:, i]
        fl = samples_raw[:, i]
        lo, hi = _hist_range(mc, pct)
        bins = np.linspace(lo, hi, n_bins + 1)
        h_mc, _ = np.histogram(mc, bins=bins, density=True)
        h_fl, _ = np.histogram(fl, bins=bins, density=True)

        ax = axes[0, i]
        _step(ax, bins, h_mc, label="MC", color="C0")
        _step(ax, bins, h_fl, label="flow", color="C1", linestyle="--")
        ax.set_ylabel("pdf")
        ax.set_yscale(yscale)
        ax.set_title(name)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(loc="best")

        ax2 = axes[1, i]
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(h_mc > 0, h_fl / h_mc, np.nan)
        _step(ax2, bins, r, color="C1")
        ax2.axhline(1.0, color="k", linestyle=":", linewidth=0.8)
        ax2.set_ylim(0.7, 1.3)
        ax2.set_ylabel("flow/MC")
        ax2.set_xlabel(name)
        ax2.grid(alpha=0.3)

    fig.suptitle("Marginal distributions: MC vs flow samples", y=1.02)
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


def plot_slices(target, samples_raw, slice_var, edges, slice_label,
                outpath, n_bins=60, pct=0.5, yscale="log"):
    """Stratify events by a continuous ``slice_var`` into the bins
    defined by ``edges`` and plot MC vs flow-sampled marginals for
    each target, in each slice.
    """
    n_slices = len(edges) - 1
    fig, axes = plt.subplots(
        3, n_slices, figsize=(3.3 * n_slices, 8.5), sharex="row"
    )
    if n_slices == 1:
        axes = axes.reshape(3, 1)
    for j in range(n_slices):
        mask = (slice_var >= edges[j]) & (slice_var < edges[j + 1])
        for i, name in enumerate(TARGET_NAMES):
            mc = target[mask, i]
            fl = samples_raw[mask, i]
            if mc.size == 0:
                continue
            lo, hi = _hist_range(mc, pct)
            bins = np.linspace(lo, hi, n_bins + 1)
            h_mc, _ = np.histogram(mc, bins=bins, density=True)
            h_fl, _ = np.histogram(fl, bins=bins, density=True)
            ax = axes[i, j]
            _step(ax, bins, h_mc, label="MC", color="C0")
            _step(ax, bins, h_fl, label="flow", color="C1", linestyle="--")
            ax.set_yscale(yscale)
            ax.grid(alpha=0.3)
            if j == 0:
                ax.set_ylabel(f"{name}  pdf")
            if i == 0:
                ax.set_title(
                    f"{slice_label}: [{edges[j]:.3g}, {edges[j+1]:.3g}]  "
                    f"n={mask.sum()}"
                )
            if i == len(TARGET_NAMES) - 1:
                ax.set_xlabel(name)
            if i == 0 and j == 0:
                ax.legend(loc="best")
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


def plot_slices_discrete(target, samples_raw, slice_var, values,
                         labels, slice_label, outpath,
                         n_bins=60, pct=0.5, yscale="log"):
    """Like plot_slices but for a discrete conditioning (e.g. charge)."""
    n_slices = len(values)
    fig, axes = plt.subplots(
        3, n_slices, figsize=(3.3 * n_slices, 8.5), sharex="row"
    )
    if n_slices == 1:
        axes = axes.reshape(3, 1)
    for j, (v, lbl) in enumerate(zip(values, labels)):
        mask = slice_var == v
        for i, name in enumerate(TARGET_NAMES):
            mc = target[mask, i]
            fl = samples_raw[mask, i]
            if mc.size == 0:
                continue
            lo, hi = _hist_range(mc, pct)
            bins = np.linspace(lo, hi, n_bins + 1)
            h_mc, _ = np.histogram(mc, bins=bins, density=True)
            h_fl, _ = np.histogram(fl, bins=bins, density=True)
            ax = axes[i, j]
            _step(ax, bins, h_mc, label="MC", color="C0")
            _step(ax, bins, h_fl, label="flow", color="C1", linestyle="--")
            ax.set_yscale(yscale)
            ax.grid(alpha=0.3)
            if j == 0:
                ax.set_ylabel(f"{name}  pdf")
            if i == 0:
                ax.set_title(f"{slice_label}: {lbl}  n={mask.sum()}")
            if i == len(TARGET_NAMES) - 1:
                ax.set_xlabel(name)
            if i == 0 and j == 0:
                ax.legend(loc="best")
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


def plot_pdf_derivatives(results, outpath, cond_ref_raw, yscale="linear"):
    """Grid plot: 3 target dims (cols) × {p, dp/dx, d²p/dx²} (rows).

    ``yscale`` applies to the p(x) row only — derivatives (which can
    be negative) always use linear scale.
    """
    n_dim = len(TARGET_NAMES)
    fig, axes = plt.subplots(3, n_dim, figsize=(4.2 * n_dim, 9.5),
                             sharex="col")
    for j, name in enumerate(TARGET_NAMES):
        x_grid, p, dp, d2p = results[name]

        ax = axes[0, j]
        ax.plot(x_grid, p, color="C0")
        ax.set_title(name)
        ax.set_ylabel("p(x | c_ref)")
        ax.set_yscale(yscale)
        ax.grid(alpha=0.3)

        ax = axes[1, j]
        ax.plot(x_grid, dp, color="C2")
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylabel("dp/dx")
        ax.grid(alpha=0.3)

        ax = axes[2, j]
        ax.plot(x_grid, d2p, color="C3")
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylabel("d²p/dx²")
        ax.set_xlabel(name)
        ax.grid(alpha=0.3)

    # Annotate the reference conditioning.
    cond_str = " ".join(f"{n}={v:+.3g}" for n, v in
                        zip(["logpt", "q", "lam", "sin_phi", "cos_phi"],
                            cond_ref_raw))
    fig.suptitle(f"pdf + derivatives at reference c_ref: {cond_str}",
                 y=1.01, fontsize=10)
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


def plot_mc_shift_reweight(
    flow, target, cond_std_np, stats, w, out_dir, device, chunk,
    shift_factors=(0.05, 0.1, 0.2),
    n_bins=80, n_std=5.0,
):
    """Validate the flow's log-density via MC shift reweighting.

    For each target axis ``i`` and magnitude ``δ = factor · target_std[i]``
    (raw units) produce histograms that all estimate the same target
    density — ``p`` translated by ``δ·e_i``:

      * **unshifted MC** (gray dotted): ``{y_m}`` reference baseline.
      * **shifted data** (red shaded): ``{y_m + δ·e_i}`` — explicit shift.
      * **1st-order Taylor** (dashed): weights
        ``w_m · (1 − δ·s_i)``. Strict linear Taylor of the ratio
        ``p(y − δ·e_i)/p(y)`` in δ. Uses only the score; can go
        negative at large δ but cannot diverge.
      * **2nd-order Taylor** (dot-dashed): weights
        ``w_m · [1 − δ·s_i + ½·δ²·(s_i² + H_ii)]``. Strict
        polynomial Taylor of the ratio in δ. Uses ``(s, H)``; can
        go negative for pathological ``(s, H)`` but cannot
        diverge numerically.
      * **base-space** (densely dash-dotted): weights
        ``w_m · exp(δ·zᵀv_i − ½δ²·‖v_i‖²)
              · [1 − δ·G_i + ½δ²·(G_i² + K − zᵀw)]``,
        with ``zᵀw ≡ K − H − ‖v‖²`` derived from existing
        quantities (no extra autograd). ``log p`` is split as
        ``log φ(z) + L`` with ``L = log|det ∂z/∂y|``. The
        Gaussian piece uses the linearized ``Δz = −δ·v_i`` —
        kept inside the exp so the exponent remains a
        concave-down quadratic (`−‖v‖² ≤ 0`) and bounded above
        per event. The `O(δ²)` contribution ``−½δ²·zᵀw`` from
        the second-order ``Δz`` correction is added as a
        polynomial modification ``K → K − zᵀw`` in the Jac
        Taylor — reproduces strict 2nd-order Taylor at ``O(δ²)``
        without the exp-blowup pathology of embedding it in the
        exponent. Negative resulting weights are clipped to
        zero; non-finite events fall back to ``_safe_weights``.
      * **exact** (solid): weights
        ``w_m · exp(log p(y_m − δ·e_i) − log p(y_m))``.

    Each panel has a ratio sub-panel plotting reweight / shifted-data
    with a ``y = 1`` reference. Because the flow gives exact log-p,
    the "exact" curve is ground truth up to MC noise — the 1st/2nd-
    order curves converging to it is a validation of the learned
    score / diagonal Hessian.

    ``cond_std_np`` is the pre-standardized conditioning (shape [N,
    n_cond]), as returned by ``apply_preproc``.
    """
    std = np.asarray(stats.target_std, dtype=np.float64)
    mean = np.asarray(stats.target_mean, dtype=np.float64)
    n_mag = len(shift_factors)
    N, d = target.shape

    target_std_np = ((target - mean) / std).astype(np.float32)
    target_std_t = torch.from_numpy(target_std_np)
    cond_std_t = torch.from_numpy(cond_std_np.astype(np.float32))

    print(
        f"  shift reweight: score + diag-Hessian at {N} events"
    )
    score_std, hess_diag_std = _flow_score_and_hessian_diag(
        flow, target_std_t, cond_std_t, device, chunk,
    )
    # Chain rule: raw-space derivatives of log p_raw (which differ
    # from log p_std by a c-dependent but x-independent constant).
    #   ∂ log p_raw / ∂ x_raw_i   = (1/std_i) ∂ log p_std / ∂ x_std_i
    #   ∂² log p_raw / ∂ x_raw_i² = (1/std_i²) ∂² log p_std / ∂ x_std_i²
    score_raw = score_std / std.astype(np.float32)
    hess_raw = hess_diag_std / (std.astype(np.float32) ** 2)

    print(
        f"  shift reweight: base-space (z, J, ∇L, diag ∇²L) at "
        f"{N} events"
    )
    z_std, J_std, G_std, K_std = _flow_base_quantities(
        flow, target_std_t, cond_std_t, device, chunk,
    )

    lp0 = _flow_log_prob_batched(
        flow, target_std_t, cond_std_t, device, chunk,
    ).numpy()

    fig, axes = plt.subplots(
        2 * n_mag, d, figsize=(5 * d, 4.0 * n_mag),
        gridspec_kw={"height_ratios": [3, 1] * n_mag},
        sharex="col",
    )
    if n_mag == 1 and axes.ndim == 1:
        axes = axes[None, :]
    for col in range(d):
        lo, hi = mean[col] - n_std * std[col], mean[col] + n_std * std[col]
        for row, factor in enumerate(shift_factors):
            ax = axes[2 * row, col]
            ax_r = axes[2 * row + 1, col]
            delta = factor * std[col]

            shifted = target[:, col] + delta
            hist_d, edges = np.histogram(
                shifted, bins=n_bins, range=(lo, hi),
                weights=w, density=True,
            )
            hist_orig, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w, density=True,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])

            ax.plot(
                centers, hist_orig, color="gray", lw=0.8,
                linestyle=":", label="unshifted MC",
            )
            ax.stairs(
                hist_d, edges, color="C3", fill=True, alpha=0.35,
                label="shifted data",
            )

            s_col = score_raw[:, col]
            H_col = hess_raw[:, col]
            # 1st-order Taylor of the ratio p(y−δ)/p(y) in δ:
            #   1 − δ·s.
            # Strict linear truncation; can go negative at large δ.
            w1 = _safe_weights(
                w, 1.0 - delta * s_col,
                label=f"shift 1st-order Taylor col={col} factor={factor}",
            )
            # 2nd-order Taylor of the ratio p(y−δ)/p(y) in δ:
            #   1 − δ·s + ½δ²·(s² + H).
            # Polynomial truncation using (s, H). Can go negative
            # for events with pathological (s, H) but cannot
            # diverge numerically.
            w2 = _safe_weights(
                w,
                1.0 - delta * s_col
                + 0.5 * delta * delta * (s_col ** 2 + H_col),
                label=f"shift 2nd-order Taylor col={col} factor={factor}",
            )

            # Base-space decomposition. Work in standardized
            # coordinates (where the flow's base is N(0, I)), so
            # the standardized shift magnitude is simply ``factor``.
            #   log p(y) = log φ(z(y)) + L(y)
            # with L = log|det ∂z/∂y|. For y → y − factor·e_col,
            # the Gaussian log-ratio under the linearized
            # Δz = −factor · v (v = ∂z/∂y_col) is
            #   log φ(z + Δz) − log φ(z)
            #     = factor · zᵀv − ½·factor²·‖v‖².
            # The log-det log-ratio is 2nd-order Taylored:
            #   L(y − factor·e_col) − L(y)
            #     ≈ −factor·G + ½factor²·K_ii
            # ⇒  exp(·) ≈ 1 − factor·G + ½factor²·(G² + K_ii).
            # The exp(Gaussian) factor is bounded above per
            # event (concave-down quadratic in factor with
            # leading coefficient −‖v‖² ≤ 0). Linearising Δz
            # drops the ``−½δ²·zᵀw`` piece of the full Gaussian
            # log-ratio (``w = ∂²z/∂y_col²``); including it in
            # the exponent via ``K − H`` reintroduces the
            # weight-spike pathology. Instead we fold it into
            # the jac polynomial via the identity
            #   zᵀw ≡ K − H − ‖v‖²,
            # i.e. replace ``K`` with ``K − zᵀw = H + ‖v‖²`` in
            # the polynomial coefficient. At O(δ²) this
            # reproduces strict 2nd-order Taylor exactly
            # (``(zv − G)² − ‖v‖² + K − zw = (zv − G)² + H``),
            # while keeping the exponent's concavity. The
            # polynomial can become negative for extreme events
            # — clipped to zero below so weights stay ≥ 0.
            v_b = J_std[:, :, col]                 # [N, d] = ∂z/∂y_col
            zv_b = np.einsum("nj,nj->n", z_std, v_b)
            vv_b = np.einsum("nj,nj->n", v_b, v_b)
            G_b = G_std[:, col]
            K_b = K_std[:, col]
            zw_b = K_b - hess_diag_std[:, col] - vv_b
            gauss_ratio_b = np.exp(
                factor * zv_b - 0.5 * factor ** 2 * vv_b
            )
            # "Linearised" jac polynomial (no zᵀw correction) —
            # retained as an intermediate validation curve so the
            # contribution of the zᵀw term is visible separately
            # from the linearisation-of-Δz approximation.
            jac_taylor_lin = (
                1.0
                - factor * G_b
                + 0.5 * factor ** 2 * (G_b ** 2 + K_b)
            )
            jac_taylor_full = (
                1.0
                - factor * G_b
                + 0.5 * factor ** 2 * (G_b ** 2 + K_b - zw_b)
            )
            # Clip negative weight factors to zero (NaN/inf pass
            # through and are handled by ``_safe_weights``).
            raw_lin = gauss_ratio_b * jac_taylor_lin
            raw_lin = np.where(raw_lin < 0, 0.0, raw_lin)
            raw_full = gauss_ratio_b * jac_taylor_full
            raw_full = np.where(raw_full < 0, 0.0, raw_full)
            w_base_lin = _safe_weights(
                w, raw_lin,
                label=f"shift base-space lin col={col} factor={factor}",
            )
            w_base = _safe_weights(
                w, raw_full,
                label=f"shift base-space col={col} factor={factor}",
            )
            hist_r1, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w1, density=True,
            )
            hist_r2, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w2, density=True,
            )
            hist_rblin, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_base_lin, density=True,
            )
            hist_rb, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_base, density=True,
            )

            # Exact: evaluate flow at y − δ·e_col.
            y_shift_raw = target.copy()
            y_shift_raw[:, col] -= delta
            y_shift_std_np = ((y_shift_raw - mean) / std).astype(np.float32)
            y_shift_t = torch.from_numpy(y_shift_std_np)
            lp_shift = _flow_log_prob_batched(
                flow, y_shift_t, cond_std_t, device, chunk,
            ).numpy()
            w_exact = _safe_weights(
                w, np.exp(lp_shift - lp0),
                label=f"shift exact col={col} factor={factor}",
            )
            hist_ex, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_exact, density=True,
            )

            show_labels = (row == 0 and col == 0)
            # zorder: exact solid (C2) keeps default (=2); overlay
            # Taylor curves use zorder=3 so they stay visible on
            # top of the exact line when very close.
            ax.plot(
                centers, hist_r1, color="C0", lw=1.0, linestyle="--",
                zorder=3,
                label=(
                    "1st-order Taylor (1 − δ·s)"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_r2, color="C1", lw=1.0, linestyle="-.",
                zorder=3,
                label=(
                    "2nd-order Taylor (1 − δs + ½δ²·(s² + H))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_rblin, color="C5", lw=1.0,
                linestyle=(0, (1, 1)),
                zorder=3,
                label=(
                    "base-space lin (exp(δzv − ½δ²‖v‖²)·(1 − δG + ½δ²(G²+K)))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_rb, color="C4", lw=1.0,
                linestyle=(0, (3, 1, 1, 1, 1, 1)),
                zorder=3,
                label=(
                    "base-space full (·(1 − δG + ½δ²(G²+K−zw)))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_ex, color="C2", lw=1.5, linestyle="-",
                label="exact (lp-diff)" if show_labels else None,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_r1, hist_d),
                color="C0", lw=1.0, linestyle="--", zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_r2, hist_d),
                color="C1", lw=1.0, linestyle="-.", zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_rblin, hist_d),
                color="C5", lw=1.0,
                linestyle=(0, (1, 1)), zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_rb, hist_d),
                color="C4", lw=1.0,
                linestyle=(0, (3, 1, 1, 1, 1, 1)), zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_ex, hist_d),
                color="C2", lw=1.5, linestyle="-",
            )
            ax.grid(alpha=0.3)
            ax_r.axhline(1.0, color="C3", lw=0.8, alpha=0.6)
            ax_r.grid(alpha=0.3)
            ax_r.set_ylim(0.7, 1.3)
            if col == 0:
                ax.set_ylabel(f"density (δ={factor:g}·σ)")
                ax_r.set_ylabel("rw / data")
            if row == 0:
                ax.set_title(TARGET_NAMES[col])
            if row == n_mag - 1:
                ax_r.set_xlabel(TARGET_NAMES[col])
            if show_labels:
                ax.legend(fontsize=6, loc="best")
    fig.suptitle(
        "MC shift vs Taylor reweighting (flow)  "
        "(δ in units of target_std; Taylor of p(y−δ)/p(y), "
        "base-space split, exact lp-diff)"
    )
    fig.tight_layout()
    _save_png_pdf(fig, os.path.join(out_dir, "mc_shift_reweight.png"))


def plot_mc_smear_reweight(
    flow, target, cond_std_np, stats, w, out_dir, device, chunk,
    smear_factors=(0.1, 0.25, 0.5),
    n_bins=80, n_std=5.0, seed=0,
    n_sigma_pack=0,
):
    """Validate the flow via rank-1 Gaussian smearing.

    For each axis ``i`` and magnitude ``σ_smear = factor · target_std[i]``:

      * **unsmeared MC** (gray dotted): ``{y_m}``.
      * **smeared data** (red shaded): ``{y_m + σ_smear · ε_m · e_i}``,
        ``ε_m ~ N(0, 1)``.
      * **2nd-order Taylor** (dot-dashed): weights
        ``w_m · [1 + ½·σ²·(s_i² + H_ii)]`` — the Gaussian
        expectation over ε of the 2nd-order Taylor of the ratio
        ``p(y − σε·e_i)/p(y)`` in σε (linear-in-σ² truncation).
        Uses ``(s, H)``; can go negative for events with very
        negative ``s² + H``. Cannot diverge numerically.
      * **base-space** (densely dash-dotted): weights
        ``w_m · MGF_base · J_exp``. ``log p`` is split as
        ``log φ(z) + L`` with ``L = log|det ∂z/∂y|`` and the
        flow's ``N(0, I)`` base. Under the linearized
        ``Δz = −σε · v`` (``v = ∂z/∂y_i``), the Gaussian
        log-ratio is exactly quadratic in ε, so its expectation
        over ε ~ N(0, 1) is the closed-form MGF
        ``MGF_base = exp(σ²(zᵀv)² / (2D)) / √D``,
        ``D = 1 + σ²‖v‖²``. Because ``D ≥ 1`` always,
        ``MGF_base`` is bounded — no pathological tail where the
        resummed data-space formula blows up. The log-det piece
        is 2nd-order Taylored with the ``zᵀw`` curvature
        correction ``K → K − zᵀw`` (identity
        ``zᵀw ≡ K − H − ‖v‖²``, no extra autograd), and
        expectation-matched against the tilted Gaussian in ε
        (mean μ = σ(zᵀv)/D, variance τ² = 1/D):
        ``J_exp = 1 − σ·G_i·μ
                  + ½σ²·(G_i² + K − zᵀw)·(μ² + τ²)``.
        At O(σ²) this reproduces strict 2nd-order Taylor
        exactly. Negative resulting weights are clipped to zero.
      * **exact** (solid): ``w_m · exp(log p(y_m − σ_smear·ε_m·e_i) −
        log p(y_m))``, using the *same* ``ε_m`` as the explicit smear
        (one-sample MC identity for the convolution).

    Any event producing a non-finite reweight factor (``inf`` /
    ``NaN``) in either formula falls back to its nominal per-
    event weight. A one-line diagnostic print reports the
    substituted fraction.
    """
    std = np.asarray(stats.target_std, dtype=np.float64)
    mean = np.asarray(stats.target_mean, dtype=np.float64)
    n_mag = len(smear_factors)
    N, d = target.shape
    rng = np.random.default_rng(seed)
    eps = rng.standard_normal(N).astype(np.float32)

    target_std_np = ((target - mean) / std).astype(np.float32)
    target_std_t = torch.from_numpy(target_std_np)
    cond_std_t = torch.from_numpy(cond_std_np.astype(np.float32))

    print(
        f"  smear reweight: score + diag-Hessian at {N} events"
    )
    score_std, hess_diag_std = _flow_score_and_hessian_diag(
        flow, target_std_t, cond_std_t, device, chunk,
    )
    score_raw = score_std / std.astype(np.float32)
    hess_raw = hess_diag_std / (std.astype(np.float32) ** 2)

    print(
        f"  smear reweight: base-space (z, J, ∇L, diag ∇²L) at "
        f"{N} events"
    )
    z_std, J_std, G_std, K_std = _flow_base_quantities(
        flow, target_std_t, cond_std_t, device, chunk,
    )

    lp0 = _flow_log_prob_batched(
        flow, target_std_t, cond_std_t, device, chunk,
    ).numpy()

    fig, axes = plt.subplots(
        2 * n_mag, d, figsize=(5 * d, 4.0 * n_mag),
        gridspec_kw={"height_ratios": [3, 1] * n_mag},
        sharex="col",
    )
    if n_mag == 1 and axes.ndim == 1:
        axes = axes[None, :]
    for col in range(d):
        lo, hi = mean[col] - n_std * std[col], mean[col] + n_std * std[col]
        for row, factor in enumerate(smear_factors):
            ax = axes[2 * row, col]
            ax_r = axes[2 * row + 1, col]
            sigma_smear = factor * std[col]

            smeared = target[:, col] + sigma_smear * eps
            hist_d, edges = np.histogram(
                smeared, bins=n_bins, range=(lo, hi),
                weights=w, density=True,
            )
            hist_orig, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w, density=True,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])

            s_col = score_raw[:, col]
            H_col = hess_raw[:, col]

            # 2nd-order Taylor weight for rank-1 Gaussian smearing:
            # Taylor-expand p(y − σε·e_i)/p(y) in σε to 2nd order
            # and take E_ε with ε ~ N(0, 1) (using E[ε] = 0,
            # E[ε²] = 1):
            #   E_ε[1 − σε·s + ½σ²ε²·(s² + H)]
            #       = 1 + ½σ²·(s² + H).
            # Uses (s, H). Can go negative for events with very
            # negative ``s² + H``; cannot diverge numerically.
            taylor_factor = 1.0 + 0.5 * sigma_smear ** 2 * (
                s_col ** 2 + H_col
            )
            w_taylor = _safe_weights(
                w, taylor_factor,
                label=f"smear 2nd-order Taylor col={col} factor={factor}",
            )
            hist_r, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_taylor, density=True,
            )

            # Base-space E_ε of the reweight ratio. Standardized
            # coordinates (base = N(0, I)), so the standardized
            # smear is just ``factor``.
            #   log p(y) = log φ(z) + L(y),  L = log|det ∂z/∂y|.
            # Under linearized Δz = −σε·v (v = ∂z/∂y_col):
            #   Gauss log-ratio(ε) = σε·zᵀv
            #                       − ½·σ²ε²·‖v‖²
            #                      ≡ a·ε + b·ε²,
            # with a = σ·zᵀv, b = −½·σ²·‖v‖². Because b ≤ 0,
            # ``exp(Gauss)`` is an integrable Gaussian in ε; the
            # ratio MGF is
            #   M = (1 − 2b)^{−1/2} · exp(a²/(2(1−2b)))
            # with D ≡ 1 − 2b = 1 + σ²·‖v‖² ≥ 1. The log-det
            # log-ratio is Taylored to O(σ²):
            #   L(y − σε·e_col) − L(y)
            #       ≈ −σε·G_col + ½·σ²ε²·K_col,
            # so exp(·) ≈ 1 − σε·G + ½·σ²ε²·(G² + K). Under the
            # tilted measure ``exp(Gauss)/M · φ(ε)`` (Gaussian
            # with mean μ = a/D, variance τ² = 1/D), moments of
            # ε give
            #   E_ε[exp(Gauss)·1]  = M,
            #   E_ε[exp(Gauss)·ε]  = M · μ,
            #   E_ε[exp(Gauss)·ε²] = M · (μ² + τ²).
            # ⇒ E_ε[exp(Gauss)·Jac-Taylor] =
            #     M · [1 − σ·G·μ
            #          + ½·σ²·(G² + K)·(μ² + τ²)].
            # Using the always-positive ``‖v‖²`` keeps D ≥ 1
            # and M bounded. The full 2nd-order Gaussian
            # log-ratio picks up an extra ``−½σ²ε²·zᵀw``
            # (``w = ∂²z/∂y_col²``) that linearisation drops;
            # putting it in the exp (``‖v‖² → K − H``) can flip
            # sign and re-introduces the weight-spike pathology,
            # so we instead fold it into the Jac polynomial via
            # ``K → K − zᵀw`` (identity: ``zᵀw ≡ K − H − ‖v‖²``,
            # so ``K − zᵀw = H + ‖v‖²``). Adding ``E_ε`` over the
            # extra ``−½σ²ε²·zᵀw·exp(Gauss)`` yields
            # ``−½σ²·zᵀw·(μ² + τ²)``, absorbed into the existing
            # ``(μ² + τ²)`` coefficient. At O(σ²) this reproduces
            # strict 2nd-order Taylor exactly; negative resulting
            # weights are clipped to zero.
            v_b = J_std[:, :, col]                 # [N, d] = ∂z/∂y_col
            zv_b = np.einsum("nj,nj->n", z_std, v_b)
            vv_b = np.einsum("nj,nj->n", v_b, v_b)
            G_b = G_std[:, col]
            K_b = K_std[:, col]
            zw_b = K_b - hess_diag_std[:, col] - vv_b
            D_b = 1.0 + factor ** 2 * vv_b
            mgf_base = (
                np.exp(factor ** 2 * zv_b ** 2 / (2.0 * D_b))
                / np.sqrt(D_b)
            )
            mu_b = factor * zv_b / D_b
            tau2_b = 1.0 / D_b
            # "Linearised" jac polynomial (no zᵀw correction) —
            # retained as an intermediate validation curve so the
            # contribution of the zᵀw term is visible separately
            # from the linearisation-of-Δz approximation.
            jac_expect_lin = (
                1.0
                - factor * G_b * mu_b
                + 0.5 * factor ** 2 * (G_b ** 2 + K_b)
                * (mu_b ** 2 + tau2_b)
            )
            jac_expect_full = (
                1.0
                - factor * G_b * mu_b
                + 0.5 * factor ** 2 * (G_b ** 2 + K_b - zw_b)
                * (mu_b ** 2 + tau2_b)
            )
            # Clip negative weight factors to zero (NaN/inf pass
            # through and are handled by ``_safe_weights``).
            raw_lin = mgf_base * jac_expect_lin
            raw_lin = np.where(raw_lin < 0, 0.0, raw_lin)
            raw_full = mgf_base * jac_expect_full
            raw_full = np.where(raw_full < 0, 0.0, raw_full)
            w_base_lin = _safe_weights(
                w, raw_lin,
                label=f"smear base-space lin col={col} factor={factor}",
            )
            w_base = _safe_weights(
                w, raw_full,
                label=f"smear base-space col={col} factor={factor}",
            )
            hist_rblin, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_base_lin, density=True,
            )
            hist_rb, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_base, density=True,
            )

            # Exact, using the same ε as the explicit smear.
            y_pert_raw = target.copy()
            y_pert_raw[:, col] -= sigma_smear * eps
            y_pert_std_np = ((y_pert_raw - mean) / std).astype(np.float32)
            y_pert_t = torch.from_numpy(y_pert_std_np)
            lp_pert = _flow_log_prob_batched(
                flow, y_pert_t, cond_std_t, device, chunk,
            ).numpy()
            exact_factor = np.exp(lp_pert - lp0)
            w_exact = _safe_weights(
                w, exact_factor,
                label=f"smear exact col={col} factor={factor}",
            )
            hist_ex, _ = np.histogram(
                target[:, col], bins=n_bins, range=(lo, hi),
                weights=w_exact, density=True,
            )

            # σ-conditioned-flow direct: the flow trained with
            # ``--cond-on-smear`` *is* the smear-integrated density,
            # so the rank-1 smear weight at this column reduces to
            # ``p_flow(y | c, σσᵀ_pack) / p_flow(y | c, σ=0)`` —
            # one perturbed flow forward, no per-event ε draw, no
            # quadrature truncation. The σσᵀ pack is one-hot at
            # diagonal index ``col`` with value ``factor²`` (since
            # ``factor`` is already in standardized target units).
            # Skipped when the flow wasn't trained σ-conditioned
            # (``n_sigma_pack=0``).
            hist_cond = None
            if n_sigma_pack > 0:
                n_cond_total = cond_std_np.shape[1]
                n_cond_base = n_cond_total - n_sigma_pack
                # Build σ_pack from a one-hot σ_std at column ``col``
                # with magnitude ``factor`` (standardized target
                # units).
                sigma_one = torch.zeros(d, dtype=torch.float32)
                sigma_one[col] = float(factor)
                sigma_pack_one = _sigma_pack_outer(
                    sigma_one.unsqueeze(0)
                ).squeeze(0)                         # [n_sigma_pack]
                cond_aug_np = cond_std_np.astype(np.float32).copy()
                cond_aug_np[:, n_cond_base:] = sigma_pack_one.numpy()
                cond_aug_t = torch.from_numpy(cond_aug_np)
                lp_cond = _flow_log_prob_batched(
                    flow, target_std_t, cond_aug_t, device, chunk,
                ).numpy()
                cond_factor = np.exp(lp_cond - lp0)
                w_cond = _safe_weights(
                    w, cond_factor,
                    label=(
                        f"smear σ-cond direct col={col} "
                        f"factor={factor}"
                    ),
                )
                hist_cond, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w_cond, density=True,
                )

            show_labels = (row == 0 and col == 0)
            ax.plot(
                centers, hist_orig, color="gray", lw=0.8,
                linestyle=":", label="unsmeared MC" if show_labels else None,
            )
            ax.stairs(
                hist_d, edges, color="C3", fill=True, alpha=0.35,
                label="smeared data" if show_labels else None,
            )
            # zorder: exact solid (C2) keeps default (=2); the
            # Taylor overlay uses zorder=3 so it stays visible
            # on top of the exact line when close.
            ax.plot(
                centers, hist_r, color="C1", lw=1.0, linestyle="-.",
                zorder=3,
                label=(
                    "2nd-order Taylor (1 + ½σ²·(s² + H))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_rblin, color="C5", lw=1.0,
                linestyle=(0, (1, 1)),
                zorder=3,
                label=(
                    "base-space lin (·(G²+K)(μ²+τ²))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_rb, color="C4", lw=1.0,
                linestyle=(0, (3, 1, 1, 1, 1, 1)),
                zorder=3,
                label=(
                    "base-space full (·(G²+K−zw)(μ²+τ²))"
                    if show_labels else None
                ),
            )
            ax.plot(
                centers, hist_ex, color="C2", lw=1.5, linestyle="-",
                label="exact (lp-diff)" if show_labels else None,
            )
            if hist_cond is not None:
                ax.plot(
                    centers, hist_cond, color="C0", lw=1.3,
                    linestyle=(0, (5, 2)), zorder=4,
                    label=(
                        "σ-cond flow direct (1 forward)"
                        if show_labels else None
                    ),
                )
            ax_r.plot(
                centers, _safe_ratio(hist_r, hist_d),
                color="C1", lw=1.0, linestyle="-.", zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_rblin, hist_d),
                color="C5", lw=1.0,
                linestyle=(0, (1, 1)), zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_rb, hist_d),
                color="C4", lw=1.0,
                linestyle=(0, (3, 1, 1, 1, 1, 1)), zorder=3,
            )
            ax_r.plot(
                centers, _safe_ratio(hist_ex, hist_d),
                color="C2", lw=1.5, linestyle="-",
            )
            if hist_cond is not None:
                ax_r.plot(
                    centers, _safe_ratio(hist_cond, hist_d),
                    color="C0", lw=1.3,
                    linestyle=(0, (5, 2)), zorder=4,
                )
            ax.grid(alpha=0.3)
            ax_r.axhline(1.0, color="C3", lw=0.8, alpha=0.6)
            ax_r.grid(alpha=0.3)
            ax_r.set_ylim(0.7, 1.3)
            if col == 0:
                ax.set_ylabel(f"density (σ_smear={factor:g}·σ)")
                ax_r.set_ylabel("rw / data")
            if row == 0:
                ax.set_title(TARGET_NAMES[col])
            if row == n_mag - 1:
                ax_r.set_xlabel(TARGET_NAMES[col])
            if show_labels:
                ax.legend(fontsize=6, loc="best")
    fig.suptitle(
        "MC smear vs Taylor reweighting (flow)  "
        "(σ_smear in units of target_std; "
        "E_ε of Taylor, base-space split, exact lp-diff)"
    )
    fig.tight_layout()
    _save_png_pdf(fig, os.path.join(out_dir, "mc_smear_reweight.png"))


# ---------------------------------------------------------------------------
# Polyhead validation
# ---------------------------------------------------------------------------

def _flow_log_w_at(flow, y_std_t, c_std_t, u_std_t, batch_size, device):
    """Evaluate true log w = log p(y - u | c) − log p(y | c) for each
    event (y-space shift weight), batched, no_grad.

    All three input arrays are shape ``[N, *]`` (same N) and live in
    standardized space. The perturbation ``u`` is in target space
    (R^{n_features}), matching the polyhead's input convention.
    """
    N = y_std_t.shape[0]
    out = np.empty(N, dtype=np.float32)
    flow.eval()
    with torch.no_grad():
        for s in range(0, N, batch_size):
            e = min(s + batch_size, N)
            y = y_std_t[s:e].to(device, non_blocking=True)
            c = c_std_t[s:e].to(device, non_blocking=True)
            u = u_std_t[s:e].to(device, non_blocking=True)
            z, ladj = flow(c).transform.call_and_ladj(y)
            z_p, ladj_p = flow(c).transform.call_and_ladj(y - u)
            log_w = -0.5 * (
                (z_p * z_p).sum(dim=-1) - (z * z).sum(dim=-1)
            ) + (ladj_p - ladj)
            out[s:e] = log_w.detach().cpu().numpy().astype(np.float32)
    return out


def _polyhead_pred_W_at(
    head, y_std_t, c_std_t, u_shift_t, sigma_vec_t,
    batch_size, device,
):
    """Arch-agnostic per-event predicted W = positivity(d), batched,
    no_grad. Dispatches on ``isinstance(head, PolyHead)``:

    * polyhead → :func:`predicted_W` on the head's basis-coef output
      (``predicted_W(coefs, u, σ_vec, joint_indices, positivity,
      basis, scale_u, scale_sigma)``).
    * MLP → dual-forward construction
      ``d = head_forward(e, u, σ_pack) − head_forward(e, 0, 0)``,
      followed by ``positivity_W(d)``. ``e = trunk_forward(y, c)`` is
      shared between the two head calls so the trunk runs once per
      batch.

    Returns CPU numpy array of shape [N].
    """
    from train_muon_response_flow import (
        PolyHead, predicted_W, _apply_positivity_W,
    )
    N = y_std_t.shape[0]
    out = np.empty(N, dtype=np.float32)
    head.eval()
    positivity = getattr(head, "positivity", "softplus")
    is_poly = isinstance(head, PolyHead)
    if is_poly:
        basis = getattr(head, "basis", "monomial")
        scale_u = float(getattr(head, "basis_scale_u", 1.0))
        scale_sigma = float(getattr(head, "basis_scale_sigma", 1.0))

    with torch.no_grad():
        for s in range(0, N, batch_size):
            e = min(s + batch_size, N)
            y = y_std_t[s:e].to(device, non_blocking=True)
            c = c_std_t[s:e].to(device, non_blocking=True)
            u = u_shift_t[s:e].to(device, non_blocking=True)
            sv = sigma_vec_t[s:e].to(device, non_blocking=True)
            if is_poly:
                joint = head(y, c)
                W = predicted_W(
                    joint, u, sv, head.joint_indices,
                    positivity=positivity,
                    basis=basis,
                    scale_u=scale_u,
                    scale_sigma=scale_sigma,
                )
            else:
                # MLP arch dual-forward.
                emb = head.trunk_forward(y, c)
                if head.head_layer1_sigma is not None:
                    sigma_pack = (
                        sv[..., head.sigma_pack_iu]
                        * sv[..., head.sigma_pack_ju]
                    )
                else:
                    sigma_pack = torch.zeros(
                        y.shape[0], 0,
                        device=y.device, dtype=y.dtype,
                    )
                f_pert = head.head_forward(emb, u, sigma_pack)
                f_zero = head.head_forward(
                    emb,
                    torch.zeros_like(u),
                    torch.zeros_like(sigma_pack),
                )
                d = f_pert - f_zero
                W = _apply_positivity_W(d, positivity)
            out[s:e] = W.detach().cpu().numpy().astype(np.float32)
    return out


def plot_polyhead_pred_vs_true(
    flow, polyhead, target_std, cond, w_event, stats, args, out_dir,
):
    """Per-mode scatter / 2D-density of polyhead-predicted W vs the
    flow's true W. One row of three panels: SHIFT, SMEAR, JOINT.
    A perfect polyhead would have all points on the diagonal.
    """
    n_features = target_std.shape[1]
    device = args.device
    n_use = min(args.polyhead_validate_n, target_std.shape[0])
    rng = np.random.default_rng(0)
    idx = rng.choice(target_std.shape[0], size=n_use, replace=False)
    y_std_t = torch.from_numpy(target_std[idx])
    c_std_t = torch.from_numpy(cond[idx])
    w_sub = w_event[idx].astype(np.float32)

    # Sample perturbations in standardized target / y-space, matching
    # the training-time shift range. ``polyhead.basis_scale_u`` equals
    # ``oversample · delta_max`` from training, so it's the right
    # magnitude bound for shifts; same for sigma.
    g = torch.Generator().manual_seed(0)
    half = float(getattr(polyhead, "basis_scale_u", 1.3))
    half_sig = float(getattr(polyhead, "basis_scale_sigma", 1.3))
    delta_shift = (torch.rand(n_use, generator=g) * 2.0 - 1.0) * half
    v_shift = torch.randn(n_use, n_features, generator=g)
    v_shift = v_shift / v_shift.norm(dim=-1, keepdim=True).clamp_min(1e-30)
    sigma_smear = torch.rand(n_use, generator=g) * half_sig
    v_smear = torch.randn(n_use, n_features, generator=g)
    v_smear = v_smear / v_smear.norm(dim=-1, keepdim=True).clamp_min(1e-30)
    delta_smear = sigma_smear * torch.randn(n_use, generator=g)
    u_shift_full = delta_shift.unsqueeze(-1) * v_shift
    sigma_vec_full = sigma_smear.unsqueeze(-1) * v_smear

    titles = ("SHIFT (σ=0)", "SMEAR (u=0)", "JOINT (both)")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), sharey=False)
    summary_lines = []
    for col, title in enumerate(titles):
        if col == 0:           # SHIFT: u=u_shift_full, σ=0, δ_smear=0
            u_eval = u_shift_full.clone()
            u_shift = u_shift_full.clone()
            sigma_vec = torch.zeros_like(sigma_vec_full)
        elif col == 1:         # SMEAR: u=0, σ=σ_vec_full, δ_smear∼N
            u_eval = delta_smear.unsqueeze(-1) * v_smear
            u_shift = torch.zeros_like(u_shift_full)
            sigma_vec = sigma_vec_full.clone()
        else:                  # JOINT
            u_eval = u_shift_full + delta_smear.unsqueeze(-1) * v_smear
            u_shift = u_shift_full.clone()
            sigma_vec = sigma_vec_full.clone()

        # True log_w from flow forwards at the y-perturbed point.
        true_lw = _flow_log_w_at(
            flow, y_std_t, c_std_t, u_eval,
            args.batch_size, device,
        )
        # Polyhead prediction.
        pred_W = _polyhead_pred_W_at(
            polyhead, y_std_t, c_std_t, u_shift, sigma_vec,
            args.batch_size, device,
        )
        true_W = np.exp(true_lw)
        # Trim to a sensible range for the scatter so a few extreme
        # tail events don't dominate the axes; report what was clipped.
        finite = np.isfinite(pred_W) & np.isfinite(true_W)
        n_inf = int((~finite).sum())
        pred_W = pred_W[finite]; true_W = true_W[finite]
        w_finite = w_sub[finite]
        lo = float(np.quantile(np.minimum(pred_W, true_W), 0.005))
        hi = float(np.quantile(np.maximum(pred_W, true_W), 0.995))
        lo = max(lo, 1e-6); hi = max(hi, lo * 10)
        # Weighted RMS log-error for the panel summary.
        log_err = np.log(np.clip(pred_W, 1e-30, None)) - np.log(
            np.clip(true_W, 1e-30, None))
        rms = float(np.sqrt(np.average(log_err ** 2, weights=w_finite)))
        bias = float(np.average(log_err, weights=w_finite))
        ax = axes[col]
        ax.hist2d(
            np.log10(np.clip(true_W, 1e-30, None)),
            np.log10(np.clip(pred_W, 1e-30, None)),
            bins=80,
            range=[
                [np.log10(lo), np.log10(hi)],
                [np.log10(lo), np.log10(hi)],
            ],
            cmin=1, cmap="viridis",
            weights=w_finite,
        )
        ax.plot(
            [np.log10(lo), np.log10(hi)],
            [np.log10(lo), np.log10(hi)],
            "r-", lw=1.0, alpha=0.7, label="diagonal",
        )
        ax.set_xlabel("log10  true W")
        ax.set_ylabel("log10  pred W")
        ax.set_title(
            f"{title}\nlog-W rms={rms:.3f}  bias={bias:+.3f}"
            + (f"  (skip {n_inf} non-finite)" if n_inf else "")
        )
        ax.legend(loc="upper left", fontsize=8)
        summary_lines.append(
            f"  {title:24s}: rms(log W err)={rms:.4f}  "
            f"bias={bias:+.4f}  (n={len(pred_W)})"
        )
    fig.suptitle("polyhead pred W  vs  true W (flow forward)")
    fig.tight_layout()
    _save_png_pdf(fig, os.path.join(out_dir, "polyhead_pred_vs_true.png"))
    print("polyhead pred-vs-true:")
    for ln in summary_lines:
        print(ln)


def plot_polyhead_logw_error(
    flow, polyhead, target_std, cond, w_event, stats, args, out_dir,
):
    """Per-event polyhead-prediction error histograms, weighted by
    ``w_event``, stratified by mode. Saves four files:

      * ``polyhead_logw_error.{png,pdf}`` (linear y)
      * ``polyhead_logw_error_log.{png,pdf}`` (log y)
        — histograms of ``log(pred_W) − log W_true``.
      * ``polyhead_relW_error.{png,pdf}`` (linear y)
      * ``polyhead_relW_error_log.{png,pdf}`` (log y)
        — histograms of ``(pred_W − W_true) / W_true``.

    Both metrics measure the same per-event error but penalize
    extreme over- vs underestimates differently. ``log(pred/true)``
    is symmetric and finite for all positive W; ``(pred−true)/true``
    is the natural "relative bias" for downstream binned-fit use, but
    can blow up when ``true ≈ 0``. The log-y versions expose tails.
    """
    n_features = target_std.shape[1]
    device = args.device
    n_use = min(args.polyhead_validate_n, target_std.shape[0])
    rng = np.random.default_rng(1)
    idx = rng.choice(target_std.shape[0], size=n_use, replace=False)
    y_std_t = torch.from_numpy(target_std[idx])
    c_std_t = torch.from_numpy(cond[idx])
    w_sub = w_event[idx].astype(np.float32)

    # Sample perturbations in standardized target / y-space, matching
    # the training-time shift range via ``polyhead.basis_scale_u`` (=
    # ``oversample · delta_max`` at training).
    g = torch.Generator().manual_seed(1)
    half = float(getattr(polyhead, "basis_scale_u", 1.3))
    half_sig = float(getattr(polyhead, "basis_scale_sigma", 1.3))
    delta_shift = (torch.rand(n_use, generator=g) * 2.0 - 1.0) * half
    v_shift = torch.randn(n_use, n_features, generator=g)
    v_shift = v_shift / v_shift.norm(dim=-1, keepdim=True).clamp_min(1e-30)
    sigma_smear = torch.rand(n_use, generator=g) * half_sig
    v_smear = torch.randn(n_use, n_features, generator=g)
    v_smear = v_smear / v_smear.norm(dim=-1, keepdim=True).clamp_min(1e-30)
    delta_smear = sigma_smear * torch.randn(n_use, generator=g)
    u_shift_full = delta_shift.unsqueeze(-1) * v_shift
    sigma_vec_full = sigma_smear.unsqueeze(-1) * v_smear

    titles = ("SHIFT", "SMEAR", "JOINT")
    # Compute pred_W and true_lw once per mode; both metrics derive
    # from these.
    panel_log = []
    panel_rel = []
    for col in range(3):
        if col == 0:
            u_eval = u_shift_full.clone()
            u_shift = u_shift_full.clone()
            sigma_vec = torch.zeros_like(sigma_vec_full)
        elif col == 1:
            u_eval = delta_smear.unsqueeze(-1) * v_smear
            u_shift = torch.zeros_like(u_shift_full)
            sigma_vec = sigma_vec_full.clone()
        else:
            u_eval = u_shift_full + delta_smear.unsqueeze(-1) * v_smear
            u_shift = u_shift_full.clone()
            sigma_vec = sigma_vec_full.clone()
        true_lw = _flow_log_w_at(
            flow, y_std_t, c_std_t, u_eval, args.batch_size, device,
        )
        pred_W = _polyhead_pred_W_at(
            polyhead, y_std_t, c_std_t, u_shift, sigma_vec,
            args.batch_size, device,
        )
        true_W = np.exp(true_lw)

        # log(pred) - log(true).
        log_err = np.log(np.clip(pred_W, 1e-30, None)) - true_lw
        finite = np.isfinite(log_err)
        le = log_err[finite]
        wle = w_sub[finite]
        lo, hi = np.quantile(le, [0.005, 0.995])
        spread = max(abs(lo), abs(hi))
        bins_log = np.linspace(-spread, spread, 81)
        rms_log = float(np.sqrt(np.average(le ** 2, weights=wle)))
        bias_log = float(np.average(le, weights=wle))
        panel_log.append((le, wle, bins_log, rms_log, bias_log))

        # (pred - true) / true. Drop events with |true| below a tiny
        # floor — those would give astronomically large rel errors
        # that aren't representative.
        true_floor = max(1e-12, 1e-6 * float(np.median(true_W)))
        rel_err = (pred_W - true_W) / np.where(
            true_W > true_floor, true_W, np.nan,
        )
        finite_r = np.isfinite(rel_err)
        re = rel_err[finite_r]
        wre = w_sub[finite_r]
        lo, hi = np.quantile(re, [0.005, 0.995])
        spread = max(abs(lo), abs(hi))
        bins_rel = np.linspace(-spread, spread, 81)
        rms_rel = float(np.sqrt(np.average(re ** 2, weights=wre)))
        bias_rel = float(np.average(re, weights=wre))
        panel_rel.append((re, wre, bins_rel, rms_rel, bias_rel))

    def _render(panel, out_basename, xlabel, suptitle):
        for yscale, suffix in (("linear", ""), ("log", "_log")):
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
            for col, title in enumerate(titles):
                e, w_f, bins, rms, bias = panel[col]
                ax = axes[col]
                ax.hist(e, bins=bins, weights=w_f, histtype="step",
                        color="C0", lw=1.5)
                ax.axvline(0.0, color="k", lw=0.8, alpha=0.5)
                ax.set_xlabel(xlabel)
                ax.set_ylabel("weighted events")
                ax.set_title(
                    f"{title}: rms={rms:.4f}  bias={bias:+.4f}"
                )
                ax.set_yscale(yscale)
            fig.suptitle(suptitle)
            fig.tight_layout()
            _save_png_pdf(
                fig,
                os.path.join(out_dir, f"{out_basename}{suffix}.png"),
            )

    _render(
        panel_log, "polyhead_logw_error",
        r"$\log(\hat W) - \log W_{\rm true}$",
        "polyhead reconstruction error per event  (log-W)",
    )
    _render(
        panel_rel, "polyhead_relW_error",
        r"$(\hat W - W_{\rm true}) / W_{\rm true}$",
        "polyhead reconstruction error per event  (rel-W)",
    )


def plot_polyhead_axis_logw_error(
    flow, polyhead, target_std, cond, w_event, stats, args, out_dir,
):
    """Per-axis polyhead reconstruction-error histograms at the
    closure shift magnitudes.

    Lays out a grid of ``log(pred_W) − log(W_true)`` histograms with:
      * rows: shift magnitudes from ``args.polyhead_shift_deltas``
      * cols: target axes (r_kappa, dlambda, dphi, ...)

    Each panel evaluates the polyhead at axis-aligned shifts
    ``u = δ · ê_axis`` for every event and shows the weighted
    distribution of the log-W error. Weighted RMS and bias are
    annotated per panel. Complements :func:`plot_polyhead_logw_error`
    (which uses random-direction shifts averaged over the unit
    sphere) by isolating axis-aligned behavior — useful for
    diagnosing per-axis polynomial-fit quality and distinguishing
    real polyhead error from statistical noise in the closure tails.

    Saves four files:
      * ``polyhead_axis_logw_error.{png,pdf}`` (linear y)
      * ``polyhead_axis_logw_error_log.{png,pdf}`` (log y)
    """
    n_features = target_std.shape[1]
    device = args.device
    n_use = min(args.polyhead_validate_n, target_std.shape[0])
    rng = np.random.default_rng(2)
    idx = rng.choice(target_std.shape[0], size=n_use, replace=False)
    y_std_t = torch.from_numpy(target_std[idx])
    c_std_t = torch.from_numpy(cond[idx])
    w_sub = w_event[idx].astype(np.float32)

    deltas = list(args.polyhead_shift_deltas)
    target_components = list(range(min(n_features, 3)))
    target_names = ["r_kappa", "dlambda", "dphi"][: len(target_components)]

    n_rows = len(deltas)
    n_cols = len(target_components)

    # Compute panel data once (independent of the y-scale).
    panel_data = [[None] * n_cols for _ in range(n_rows)]
    for r, factor in enumerate(deltas):
        for c, axis in enumerate(target_components):
            v = torch.zeros(n_use, n_features, dtype=torch.float32)
            v[:, axis] = 1.0
            u = float(factor) * v
            sigma_vec = torch.zeros_like(v)
            true_lw = _flow_log_w_at(
                flow, y_std_t, c_std_t, u, args.batch_size, device,
            )
            pred_W = _polyhead_pred_W_at(
                polyhead, y_std_t, c_std_t, u, sigma_vec,
                args.batch_size, device,
            )
            log_err = np.log(np.clip(pred_W, 1e-30, None)) - true_lw
            finite = np.isfinite(log_err)
            le = log_err[finite]
            wle = w_sub[finite]
            rms = float(np.sqrt(np.average(le ** 2, weights=wle)))
            bias = float(np.average(le, weights=wle))
            panel_data[r][c] = (le, wle, rms, bias)

    for yscale, suffix in (("linear", ""), ("log", "_log")):
        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(4.0 * n_cols, 3.0 * n_rows),
            sharex="row", sharey="row", squeeze=False,
        )
        for r, factor in enumerate(deltas):
            # Common bin range across columns at this row, so the
            # histograms are visually comparable across axes.
            spreads = []
            for c in range(n_cols):
                le, *_ = panel_data[r][c]
                lo, hi = np.quantile(le, [0.005, 0.995])
                spreads.append(max(abs(lo), abs(hi)))
            spread = max(spreads) if spreads else 1.0
            bins = np.linspace(-spread, spread, 81)
            for c, name in enumerate(target_names):
                ax = axes[r][c]
                le, wle, rms, bias = panel_data[r][c]
                ax.hist(
                    le, bins=bins, weights=wle, histtype="step",
                    color="C0", lw=1.5,
                )
                ax.axvline(0.0, color="k", lw=0.8, alpha=0.5)
                ax.set_title(
                    f"|δ|={factor:g} along {name}: "
                    f"rms={rms:.4f} bias={bias:+.4f}"
                )
                ax.set_yscale(yscale)
                if r == n_rows - 1:
                    ax.set_xlabel(
                        r"$\log(\hat W) - \log W_{\rm true}$"
                    )
                if c == 0:
                    ax.set_ylabel("weighted events")
        fig.suptitle(
            "polyhead reconstruction error, axis-aligned shifts"
        )
        fig.tight_layout()
        _save_png_pdf(
            fig,
            os.path.join(
                out_dir, f"polyhead_axis_logw_error{suffix}.png",
            ),
        )


def plot_polyhead_reweight_closure(
    flow, polyhead, target, target_std, cond, w_event, stats,
    args, out_dir, mode, n_sigma_pack=0,
):
    """Closure plot: histogram each target component four ways and
    compare, with a ratio sub-panel below each main panel. The
    polyhead predicts y-space shift/smear weights, so the
    explicitly y-shifted/smeared raw MC is the literal closure target
    here — `flow-true`, `polyhead`, and `shifted/smeared MC` should
    all overlap when the polyhead is well-trained.

    ``polyhead`` may be ``None`` (no head trained / available). The
    head curves are then omitted but the rest of the panel — raw MC,
    shifted/smeared MC reference, flow-true K=1 reweight, and (in
    smear mode with σ-cond) the σ-cond direct curve — are still
    produced as a flow-only closure check.

      * ``raw MC`` (gray dotted): unweighted MC.
      * ``shifted/smeared MC`` (black solid): raw MC y values with a
        deterministic y-space shift (or per-event Gaussian smear)
        applied to target component ``tcol``. Magnitude ``factor *
        stats.target_std[tcol]``. **Literal closure target.**
      * ``flow-true W`` (red): original MC y reweighted by
        ``w = p(y - u | c) / p(y | c)`` computed from the flow.
        For ``mode='smear'`` this uses K=1 stochastic ε per event,
        so it carries per-event MC noise on top of the flow's
        density-ratio noise.
      * ``σ-cond flow direct`` (green dot-dashed, *smear mode only*,
        only when the flow was trained with ``--cond-on-smear``):
        single perturbed forward of the σ-conditioned flow,
        ``p_flow(y | c, σσᵀ_pack) / p_flow(y | c, σ=0)``. No per-
        event ε draw, no quadrature truncation — the smoothest
        flow-side reference for the smear closure.
      * ``polyhead pred W`` (blue dashed): original MC y reweighted by
        the polyhead's prediction.

    Ratio sub-panel: the flow-true and polyhead-pred reweights
    divided by the shifted/smeared MC (the literal closure target).
    Both should sit on 1 when well-trained. The raw (unshifted/un-
    smeared) MC is shown only in the main panel — its ratio against
    the shifted/smeared MC is just the inverse of the perturbation
    and would dominate the y-range without carrying closure info.

    Iteration is over the first three target components (one panel
    column each: r_kappa, dlambda, dphi), with shift/smear applied
    along that single target axis.

    ``mode`` is one of ``"shift"`` (u_shift = δ·v_y, σ_vec = 0) or
    ``"smear"`` (u_shift = 0, σ_vec = σ·v_y with stochastic K=1
    δ_smear ~ N(0, σ²) per event).
    """
    assert mode in ("shift", "smear")
    n_features = target.shape[1]
    device = args.device
    factors = (
        args.polyhead_shift_deltas if mode == "shift"
        else args.polyhead_smear_sigmas
    )
    # Iterate over the first n target components; each gets a column
    # with the perturbation applied along its standardized basis
    # direction in y-space.
    target_components = list(range(min(n_features, 3)))

    n_rows = len(factors)
    n_cols = len(target_components)
    fig, axes = plt.subplots(
        2 * n_rows, n_cols,
        figsize=(4.0 * n_cols, 3.5 * n_rows),
        sharex="col",
        gridspec_kw={
            "height_ratios": [3, 1] * n_rows,
            "hspace": 0.05,
        },
        squeeze=False,
        layout="constrained",
    )
    target_std_t = torch.from_numpy(target_std)
    cond_t = torch.from_numpy(cond)
    target_std_per_dim = np.asarray(stats.target_std, dtype=np.float64)

    # σ-conditioned-flow direct curve (smear mode only). The σ=0
    # baseline log p is the same for every (factor, tcol), so compute
    # it once. The cond passed in already has the trailing σ_pack
    # columns set to zero (σ=0 baseline) when the flow was trained
    # with ``--cond-on-smear``.
    use_cond_direct = (mode == "smear" and n_sigma_pack > 0)
    if use_cond_direct:
        n_cond_total = cond.shape[1]
        n_cond_base = n_cond_total - n_sigma_pack
        lp0_cond = _flow_log_prob_batched(
            flow, target_std_t, cond_t, device, args.batch_size,
        ).numpy()

    for r, factor in enumerate(factors):
        for cidx, tcol in enumerate(target_components):
            # Unit basis vector in standardized target space along
            # component tcol.
            v = torch.zeros(
                target_std.shape[0], n_features, dtype=torch.float32,
            )
            v[:, tcol] = 1.0
            if mode == "shift":
                u_eval = factor * v
                u_shift = factor * v
                sigma_vec = torch.zeros_like(v)
            else:
                # SMEAR: stochastic K=1 — δ_smear ~ N(0, σ²) per event.
                g = torch.Generator().manual_seed(100 * r + tcol)
                d = factor * torch.randn(
                    target_std.shape[0], generator=g,
                )
                u_eval = d.unsqueeze(-1) * v
                u_shift = torch.zeros_like(v)
                sigma_vec = factor * v

            true_lw = _flow_log_w_at(
                flow, target_std_t, cond_t, u_eval,
                args.batch_size, device,
            )
            if polyhead is not None:
                pred_W = _polyhead_pred_W_at(
                    polyhead, target_std_t, cond_t, u_shift, sigma_vec,
                    args.batch_size, device,
                )
            else:
                pred_W = None
            true_W = np.exp(true_lw)
            # σ-conditioned-flow direct: replace the σ_pack tail of
            # cond with the one-hot σ_pack at column tcol with
            # magnitude ``factor`` (standardized target units). The
            # smear weight is then ``p(y | c, σσᵀ) / p(y | c, 0)`` —
            # one perturbed forward, no per-event ε draw. Skipped
            # outside smear mode and when the flow wasn't σ-cond.
            cond_W = None
            if use_cond_direct:
                sigma_one = torch.zeros(n_features, dtype=torch.float32)
                sigma_one[tcol] = float(factor)
                sigma_pack_one = _sigma_pack_outer(
                    sigma_one.unsqueeze(0)
                ).squeeze(0)
                cond_aug_np = cond.astype(np.float32).copy()
                cond_aug_np[:, n_cond_base:] = sigma_pack_one.numpy()
                cond_aug_t = torch.from_numpy(cond_aug_np)
                lp_cond = _flow_log_prob_batched(
                    flow, target_std_t, cond_aug_t, device,
                    args.batch_size,
                ).numpy()
                cond_W = np.exp(lp_cond - lp0_cond)
            # Use the closure-specific percentile so this plot's x-
            # range and binning match the direct shift-reweight
            # diagnostic's ``shift_closure.png`` style for side-by-
            # side comparability. Falls back to ``range_percentile``
            # if the closure-specific arg isn't set (older configs).
            closure_pct = float(getattr(
                args, "polyhead_closure_percentile",
                args.range_percentile,
            ))
            lo, hi = np.percentile(
                target[:, tcol],
                [closure_pct, 100 - closure_pct],
            )
            bins = np.linspace(lo, hi, args.n_bins + 1)
            centers = 0.5 * (bins[:-1] + bins[1:])

            # Reference: y-space shifted / smeared raw MC. Magnitude
            # ``factor * target_std[tcol]`` matches the polyhead's
            # standardized-space input ``factor * v``, so this is a
            # literal closure target in this design (perturbations
            # are in y-space).
            dy = float(factor) * float(target_std_per_dim[tcol])
            if mode == "shift":
                # Shift y → y + dy; reweight is
                # p(y - dy | c) / p(y | c) so the shifted MC sits
                # at original y values + dy. Sign matches polyhead.
                y_perturbed = target[:, tcol] + dy
                ref_label = f"raw MC + {factor:g}·σ_y"
            else:
                rng_y = np.random.default_rng(1000 * r + tcol + 1)
                y_perturbed = (
                    target[:, tcol] + dy * rng_y.standard_normal(
                        target.shape[0]
                    ).astype(target.dtype)
                )
                ref_label = f"raw MC + N(0, ({factor:g}·σ_y)²)"

            w64 = np.asarray(w_event, dtype=np.float64)
            h_raw, e_raw = _weighted_hist_err(
                target[:, tcol], bins, w64,
            )
            h_ref, e_ref = _weighted_hist_err(
                y_perturbed, bins, w64,
            )
            h_true, e_true = _weighted_hist_err(
                target[:, tcol], bins,
                w64 * np.asarray(true_W, dtype=np.float64),
            )
            if pred_W is not None:
                h_pred, e_pred = _weighted_hist_err(
                    target[:, tcol], bins,
                    w64 * np.asarray(pred_W, dtype=np.float64),
                )
            else:
                h_pred = e_pred = None
            if cond_W is not None:
                h_cond, e_cond = _weighted_hist_err(
                    target[:, tcol], bins,
                    w64 * np.asarray(cond_W, dtype=np.float64),
                )
            else:
                h_cond = e_cond = None

            # Main panel.
            ax_main = axes[2 * r][cidx]
            ax_main.step(centers, h_raw, where="mid",
                         color="0.4", linestyle=":",
                         label="raw MC", lw=1.0)
            ax_main.step(centers, h_ref, where="mid",
                         color="k", label=ref_label, lw=1.0)
            _stepped_errorbar(ax_main, centers, h_ref, e_ref, "k")
            ax_main.step(centers, h_true, where="mid",
                         color="C3", label="flow-true W reweight",
                         lw=1.0)
            _stepped_errorbar(ax_main, centers, h_true, e_true, "C3")
            if h_cond is not None:
                ax_main.step(
                    centers, h_cond, where="mid",
                    color="C2", linestyle="-.",
                    label="σ-cond flow direct (1 forward)",
                    lw=1.0,
                )
                _stepped_errorbar(
                    ax_main, centers, h_cond, e_cond, "C2",
                )
            if h_pred is not None:
                ax_main.step(centers, h_pred, where="mid",
                             color="C0", linestyle="--",
                             label="polyhead pred W reweight", lw=1.0)
                _stepped_errorbar(
                    ax_main, centers, h_pred, e_pred, "C0",
                )
            ax_main.set_ylabel("events")
            ax_main.set_title(
                f"{mode}: |δ|={factor}  along {TARGET_NAMES[tcol] if tcol < len(TARGET_NAMES) else f'target[{tcol}]'}"
            )
            if r == 0 and cidx == 0:
                ax_main.legend(loc="best", fontsize=7)

            # Ratio panel: reweight curves / shifted-or-smeared MC
            # (the literal closure target). The raw (unshifted/un-
            # smeared) MC ratio is excluded — it just plots the
            # inverse of the perturbation and dominates the y-range
            # without carrying closure information.
            ax_ratio = axes[2 * r + 1][cidx]
            ratio_true, eratio_true = _ratio_with_err(
                h_true, e_true, h_ref, e_ref,
            )
            if h_pred is not None:
                ratio_pred, eratio_pred = _ratio_with_err(
                    h_pred, e_pred, h_ref, e_ref,
                )
            else:
                ratio_pred = eratio_pred = None
            if h_cond is not None:
                ratio_cond, eratio_cond = _ratio_with_err(
                    h_cond, e_cond, h_ref, e_ref,
                )
            else:
                ratio_cond = eratio_cond = None
            ax_ratio.step(centers, ratio_true, where="mid",
                          color="C3", lw=1.0)
            _stepped_errorbar(
                ax_ratio, centers, ratio_true, eratio_true, "C3",
            )
            if ratio_cond is not None:
                ax_ratio.step(
                    centers, ratio_cond, where="mid",
                    color="C2", linestyle="-.", lw=1.0,
                )
                _stepped_errorbar(
                    ax_ratio, centers, ratio_cond, eratio_cond, "C2",
                )
            if ratio_pred is not None:
                ax_ratio.step(centers, ratio_pred, where="mid",
                              color="C0", linestyle="--", lw=1.0)
                _stepped_errorbar(
                    ax_ratio, centers, ratio_pred, eratio_pred, "C0",
                )
            ax_ratio.axhline(1.0, color="k", lw=0.5, alpha=0.5)
            ax_ratio.set_ylabel(
                f"/ {'shifted' if mode == 'shift' else 'smeared'} MC",
                fontsize=8,
            )
            ax_ratio.set_xlabel(
                TARGET_NAMES[tcol] if tcol < len(TARGET_NAMES)
                else f"target[{tcol}]"
            )
            ratio_arrays = [np.atleast_1d(ratio_true)]
            if ratio_pred is not None:
                ratio_arrays.append(np.atleast_1d(ratio_pred))
            if ratio_cond is not None:
                ratio_arrays.append(np.atleast_1d(ratio_cond))
            ratios = np.concatenate(ratio_arrays)
            ratios = ratios[np.isfinite(ratios)]
            if ratios.size:
                lo_r, hi_r = np.quantile(ratios, [0.05, 0.95])
                pad = max(0.05, 0.5 * (hi_r - lo_r))
                ax_ratio.set_ylim(
                    max(0.0, min(lo_r - pad, 1.0 - pad)),
                    max(hi_r + pad, 1.0 + pad),
                )

    title = (
        f"polyhead {mode}-reweight closure"
        if polyhead is not None
        else f"flow {mode}-reweight closure (no head)"
    )
    fig.suptitle(title)
    _save_png_pdf(
        fig,
        os.path.join(out_dir, f"polyhead_{mode}_closure.png"),
    )


def plot_polyhead_reweight_error_ratio(
    flow, polyhead, target, target_std, cond, w_event, stats,
    args, out_dir, mode,
):
    """Per-bin √Σw² error comparison: flow-true reweight, polyhead
    reweight, and the variance-optimal constant-per-bin reweight as
    the achievable floor. Parallel to
    :func:`plot_polyhead_reweight_closure` but plotting per-bin
    statistical errors rather than bin contents.

    Per (factor, axis) cell, two stacked panels:
      * Top (log y): four overlaid step lines.
          - σ_ref      = √Σ wᵢ² at perturbed y (literal shifted /
            smeared MC).
          - σ_true     = √Σ (wᵢ · W_true)² at unperturbed y
            (flow-true reweight).
          - σ_pred     = √Σ (wᵢ · W_pred)² at unperturbed y
            (polyhead reweight).
          - σ_optimal  = (h_ref / h_raw) · σ_raw — the variance-
            minimal per-bin error achievable by **any** reweight
            that matches the perturbed-MC bin totals (constant
            per-bin reweight). Lower bound on what either flow or
            polyhead reweight can possibly reach.
      * Bottom (linear y): three ratio lines, all relative to
        ``σ_optimal``.
          - ``σ_true / σ_optimal``  — flow-true reweight inflation
            over the floor.
          - ``σ_pred / σ_optimal``  — polyhead reweight inflation.
          - ``σ_ref  / σ_optimal``  — literal shifted / smeared
            sample compared to the floor.

    Both ``σ_true / σ_optimal`` and ``σ_pred / σ_optimal`` are
    bounded below by 1; values close to 1 mean the corresponding
    reweight is achieving near-uniform per-bin reweight (intra-bin
    reweight values vary little). Values ≫ 1 mean the reweight is
    assigning highly variable weights within bin — visible signal
    of intra-bin reweight roughness. ``σ_pred / σ_true`` (read off
    by eye) tells you how much per-bin variance the polyhead adds
    on top of the flow's already-imperfect reweight.

    ``mode = 'shift'`` uses the deterministic shifted-MC reference;
    ``mode = 'smear'`` uses the K=1 stochastic smeared-MC reference.
    Layout, binning, and perturbation scheme match
    :func:`plot_polyhead_reweight_closure`.
    """
    assert mode in ("shift", "smear")
    n_features = target.shape[1]
    device = args.device
    factors = (
        args.polyhead_shift_deltas if mode == "shift"
        else args.polyhead_smear_sigmas
    )
    target_components = list(range(min(n_features, 3)))

    n_rows = len(factors)
    n_cols = len(target_components)
    fig, axes = plt.subplots(
        2 * n_rows, n_cols,
        figsize=(4.0 * n_cols, 3.5 * n_rows),
        sharex="col",
        gridspec_kw={
            "height_ratios": [3, 1] * n_rows,
            "hspace": 0.05,
        },
        squeeze=False,
        layout="constrained",
    )
    target_std_t = torch.from_numpy(target_std)
    cond_t = torch.from_numpy(cond)
    target_std_per_dim = np.asarray(stats.target_std, dtype=np.float64)

    for r, factor in enumerate(factors):
        for cidx, tcol in enumerate(target_components):
            v = torch.zeros(
                target_std.shape[0], n_features, dtype=torch.float32,
            )
            v[:, tcol] = 1.0
            if mode == "shift":
                u_eval = factor * v
                u_shift = factor * v
                sigma_vec = torch.zeros_like(v)
            else:
                # SMEAR: stochastic K=1 — δ_smear ~ N(0, σ²) per event.
                # Same RNG seed convention as plot_polyhead_reweight_closure
                # so the two plots use the same perturbed sample (the
                # closure and error-ratio panels are interpreted
                # together).
                g = torch.Generator().manual_seed(100 * r + tcol)
                d = factor * torch.randn(
                    target_std.shape[0], generator=g,
                )
                u_eval = d.unsqueeze(-1) * v
                u_shift = torch.zeros_like(v)
                sigma_vec = factor * v

            true_lw = _flow_log_w_at(
                flow, target_std_t, cond_t, u_eval,
                args.batch_size, device,
            )
            pred_W = _polyhead_pred_W_at(
                polyhead, target_std_t, cond_t, u_shift, sigma_vec,
                args.batch_size, device,
            )
            true_W = np.exp(true_lw)
            closure_pct = float(getattr(
                args, "polyhead_closure_percentile",
                args.range_percentile,
            ))
            lo, hi = np.percentile(
                target[:, tcol],
                [closure_pct, 100 - closure_pct],
            )
            bins = np.linspace(lo, hi, args.n_bins + 1)
            centers = 0.5 * (bins[:-1] + bins[1:])

            # Reference: y-space shifted / smeared raw MC. Same
            # construction as plot_polyhead_reweight_closure (and
            # the same RNG seeds for the smear sample).
            dy = float(factor) * float(target_std_per_dim[tcol])
            if mode == "shift":
                y_perturbed = target[:, tcol] + dy
                ref_label = f"shifted MC ({factor:g}·σ_y)"
            else:
                rng_y = np.random.default_rng(1000 * r + tcol + 1)
                y_perturbed = (
                    target[:, tcol] + dy * rng_y.standard_normal(
                        target.shape[0]
                    ).astype(target.dtype)
                )
                ref_label = (
                    f"smeared MC ({factor:g}·σ_y, K=1 stochastic)"
                )

            # Four σ vectors:
            #   raw      — unperturbed MC at original y (bin-fill).
            #   ref      — perturbed MC at shifted/smeared y.
            #   true     — flow-true reweighted MC at original y.
            #   pred     — polyhead reweighted MC at original y.
            w64 = np.asarray(w_event, dtype=np.float64)
            h_raw, e_raw = _weighted_hist_err(
                target[:, tcol], bins, w64,
            )
            h_ref, e_ref = _weighted_hist_err(y_perturbed, bins, w64)
            _, e_true = _weighted_hist_err(
                target[:, tcol], bins,
                w64 * np.asarray(true_W, dtype=np.float64),
            )
            _, e_pred = _weighted_hist_err(
                target[:, tcol], bins,
                w64 * np.asarray(pred_W, dtype=np.float64),
            )
            # Optimal: constant-per-bin reweight = h_ref / h_raw.
            # σ_opt = (h_ref / h_raw) · σ_raw. NaN where h_raw == 0
            # (no events to reweight) so the line drops out cleanly
            # in those bins.
            scale = h_ref / np.where(h_raw > 0, h_raw, np.nan)
            e_optimal = scale * e_raw

            ax_main = axes[2 * r][cidx]
            ax_main.step(
                centers, e_ref, where="mid",
                color="k", lw=1.0, label=ref_label,
            )
            ax_main.step(
                centers, e_true, where="mid",
                color="C3", lw=1.0,
                label=r"flow-true: $\sqrt{\sum (w W_{\rm true})^2}$",
            )
            ax_main.step(
                centers, e_pred, where="mid",
                color="C0", linestyle="--", lw=1.0,
                label=r"polyhead: $\sqrt{\sum (w \hat W)^2}$",
            )
            ax_main.step(
                centers, e_optimal, where="mid",
                color="C2", linestyle=":", lw=1.0,
                label=(
                    r"optimal: $(h_{\rm ref}/h_{\rm raw})"
                    r" \sqrt{\sum w^2}$"
                ),
            )
            ax_main.set_yscale("log")
            ax_main.set_ylabel(r"$\sigma$ per bin")
            tname = (
                TARGET_NAMES[tcol] if tcol < len(TARGET_NAMES)
                else f"target[{tcol}]"
            )
            ax_main.set_title(
                f"{mode}: |{'δ' if mode == 'shift' else 'σ'}|="
                f"{factor:g} along {tname}",
                fontsize=10,
            )
            if r == 0 and cidx == 0:
                ax_main.legend(loc="best", fontsize=7)

            # Ratio panel — three step lines vs the σ_optimal floor.
            # Bins with h_raw == 0 give σ_optimal = NaN and drop all
            # lines cleanly.
            denom = np.where(e_optimal > 0, e_optimal, np.nan)
            ratio_true = e_true / denom
            ratio_pred = e_pred / denom
            ratio_ref = e_ref / denom
            ax_ratio = axes[2 * r + 1][cidx]
            ax_ratio.step(
                centers, ratio_true, where="mid",
                color="C3", lw=1.0,
                label=r"flow-true / optimal",
            )
            ax_ratio.step(
                centers, ratio_pred, where="mid",
                color="C0", linestyle="--", lw=1.0,
                label=r"polyhead / optimal",
            )
            ax_ratio.step(
                centers, ratio_ref, where="mid",
                color="k", lw=1.0,
                label="ref / optimal",
            )
            ax_ratio.axhline(1.0, color="k", lw=0.5, alpha=0.5)
            ax_ratio.set_ylabel(
                r"$\sigma / \sigma_{\rm optimal}$", fontsize=8,
            )
            ax_ratio.set_xlabel(tname)
            if r == 0 and cidx == 0:
                ax_ratio.legend(loc="best", fontsize=7)
            ratios_all = np.concatenate([
                ratio_true[np.isfinite(ratio_true)],
                ratio_pred[np.isfinite(ratio_pred)],
                ratio_ref[np.isfinite(ratio_ref)],
            ])
            if ratios_all.size:
                lo_r, hi_r = np.quantile(ratios_all, [0.05, 0.95])
                pad = max(0.05, 0.5 * (hi_r - lo_r))
                ax_ratio.set_ylim(
                    max(0.0, min(lo_r - pad, 0.9)),
                    max(hi_r + pad, 1.1),
                )

    fig.suptitle(
        f"polyhead {mode}-reweight per-bin error vs {mode} MC: "
        r"ratio relative to the optimal-reweight floor"
    )
    _save_png_pdf(
        fig,
        os.path.join(out_dir, f"polyhead_{mode}_closure_error_ratio.png"),
    )


def plot_logp(log_p, outpath, n_bins=100, yscale="linear"):
    lo, hi = np.percentile(log_p, [0.1, 99.9])
    bins = np.linspace(lo, hi, n_bins + 1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(log_p, bins=bins, histtype="step", color="C0")
    ax.set_yscale(yscale)
    ax.set_xlabel("log p(x|c)  [raw-space]")
    ax.set_ylabel("events")
    ax.set_title(
        f"mean={log_p.mean():+.4f}  std={log_p.std():.4f}  "
        f"min={log_p.min():+.2f}  max={log_p.max():+.2f}"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    _save_png_pdf(fig, outpath)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    output_dir = args.output_dir or os.path.dirname(
        os.path.abspath(args.checkpoint)
    )
    os.makedirs(output_dir, exist_ok=True)

    print(f"loading checkpoint {args.checkpoint}")
    flow, stats, ckpt = load_flow_from_checkpoint(args.checkpoint, args.device)
    if "epoch" in ckpt:
        print(
            f"  epoch {ckpt['epoch']}  "
            f"train_nll {ckpt.get('train_nll', float('nan')):+.4f}  "
            f"val_nll {ckpt.get('val_nll', float('nan')):+.4f}"
        )

    print(f"loading MC events from {len(args.input_files)} file(s)")
    (
        pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q, _w,
    ) = load_ntuples(
        args.input_files,
        args.tree,
        max_muons=args.n_events,
        pt_min=args.pt_min,
        pt_max=args.pt_max,
        eta_max=args.eta_max,
        threads=args.threads,
        max_events=args.max_events,
    )
    target, cond_raw = compute_targets_and_conditioning(
        pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q
    )
    target_std, cond = apply_preproc(target, cond_raw, stats)

    # If the flow was trained with --cond-on-smear, its conditioner
    # expects an augmented input ``(c_orig, σ_pack)`` where ``σ_pack``
    # is the F·(F+1)/2-component upper-triangular pack of ``σσᵀ``.
    # Diagnostics evaluate the flow at the unsmeared baseline (σ=0),
    # so pad ``cond`` with the corresponding zeros once here — every
    # downstream consumer that takes ``cond`` as an array picks up the
    # right shape.
    flow_config = ckpt.get("flow_config", {})
    if flow_config.get("cond_on_smear", False):
        n_sigma_pack = int(flow_config.get("n_sigma_pack", 0))
        if n_sigma_pack > 0:
            cond = np.concatenate(
                [
                    cond,
                    np.zeros(
                        (cond.shape[0], n_sigma_pack),
                        dtype=cond.dtype,
                    ),
                ],
                axis=-1,
            )
            print(
                f"flow trained with --cond-on-smear: padding cond "
                f"with {n_sigma_pack} σ-pack zeros (σ=0 baseline) "
                f"→ n_cond_total={cond.shape[1]}"
            )

    target_std_t = torch.from_numpy(target_std)
    cond_t = torch.from_numpy(cond)

    print(
        f"evaluating flow on {target_std_t.shape[0]} events "
        f"(device {args.device}, batch {args.batch_size})"
    )
    log_p, samples_raw = evaluate_flow(
        flow, target_std_t, cond_t, stats, args.device, args.batch_size
    )

    # Print raw-space NLL statistics
    nll_mean = -float(log_p.mean())
    print(f"raw-space NLL per event: mean {nll_mean:+.4f}  "
          f"std {log_p.std():.4f}")

    # --- Plots ---
    # Each family produces one "_log" (log-y) and one "_linear" file so
    # the core of the distribution and the tails can both be inspected.
    log_pt_gen = cond_raw["log_pt_gen"]
    edges_pt = np.quantile(log_pt_gen, [0.0, 0.25, 0.5, 0.75, 1.0])
    lam = cond_raw["lambda_gen"]
    edges_lam = np.quantile(lam, [0.0, 0.25, 0.5, 0.75, 1.0])
    charge = cond_raw["charge"]

    # Reference conditioning for the 1D pdf-derivative slices:
    # per-component median of the raw conditioning, matching the flow's
    # bulk population.
    cond_ref_raw = np.array(
        [np.median(cond_raw[name]) for name in stats.cond_names],
        dtype=np.float32,
    )
    # Project (sin_phi_gen, cos_phi_gen) back onto the unit circle.
    # The per-component median is NOT on the manifold the flow was
    # trained on — for a nearly-uniform phi distribution each median
    # is ≈ 0 independently, giving a reference at (≈ 0, ≈ 0), far
    # from the unit circle where every real muon lives. Extrapolating
    # the flow to that point yields spurious shifts in the conditional
    # peak. Re-derive phi_ref from the (sin, cos) medians, then set
    # the pair back to (sin(phi_ref), cos(phi_ref)) which is on the
    # circle by construction. If only one of the two names is
    # present, or neither, skip silently.
    try:
        i_sin = stats.cond_names.index("sin_phi_gen")
        i_cos = stats.cond_names.index("cos_phi_gen")
        phi_ref = float(
            np.arctan2(cond_ref_raw[i_sin], cond_ref_raw[i_cos])
        )
        cond_ref_raw[i_sin] = np.float32(np.sin(phi_ref))
        cond_ref_raw[i_cos] = np.float32(np.cos(phi_ref))
    except ValueError:
        pass
    grid_ranges = [
        _hist_range(target[:, i], args.range_percentile)
        for i in range(target.shape[1])
    ]
    print(
        "evaluating p, dp/dx, d^2p/dx^2 on 1D slices "
        f"at c_ref = {dict(zip(stats.cond_names, cond_ref_raw.tolist()))}"
    )
    n_sigma_pack = int(flow_config.get("n_sigma_pack", 0))
    deriv_results = pdf_derivatives_on_grid(
        flow, stats, cond_ref_raw, grid_ranges,
        n_grid=500, device=args.device, n_sigma_pack=n_sigma_pack,
    )
    print(
        "evaluating log p, s = ∂log p/∂x, H = ∂²log p/∂x² on 1D slices"
    )
    logp_deriv_results = logp_derivatives_on_grid(
        flow, stats, cond_ref_raw, grid_ranges,
        n_grid=500, device=args.device, n_sigma_pack=n_sigma_pack,
    )
    plot_logp_derivatives(
        logp_deriv_results,
        os.path.join(output_dir, "logp_derivatives.png"),
        cond_ref_raw=cond_ref_raw,
    )

    for yscale in ("log", "linear"):
        suffix = f"_{yscale}"
        plot_logp(
            log_p,
            os.path.join(output_dir, f"log_p_per_event{suffix}.png"),
            n_bins=args.n_bins,
            yscale=yscale,
        )
        plot_marginals(
            target, samples_raw,
            os.path.join(output_dir, f"marginals{suffix}.png"),
            n_bins=args.n_bins,
            pct=args.range_percentile,
            yscale=yscale,
        )
        plot_slices(
            target, samples_raw, log_pt_gen, edges_pt,
            slice_label="log(pt_gen)",
            outpath=os.path.join(output_dir, f"slices_log_pt{suffix}.png"),
            n_bins=max(40, args.n_bins // 2),
            yscale=yscale,
        )
        plot_slices(
            target, samples_raw, lam, edges_lam,
            slice_label="lambda_gen",
            outpath=os.path.join(output_dir, f"slices_lambda{suffix}.png"),
            n_bins=max(40, args.n_bins // 2),
            yscale=yscale,
        )
        plot_slices_discrete(
            target, samples_raw, charge,
            values=[+1.0, -1.0],
            labels=["mu+ (q=+1)", "mu- (q=-1)"],
            slice_label="charge",
            outpath=os.path.join(output_dir, f"slices_charge{suffix}.png"),
            n_bins=max(40, args.n_bins // 2),
            yscale=yscale,
        )
        plot_pdf_derivatives(
            deriv_results,
            os.path.join(output_dir, f"pdf_derivatives{suffix}.png"),
            cond_ref_raw=cond_ref_raw,
            yscale=yscale,
        )

    if not args.skip_reweight:
        # Unit weights for the reweight closure: per-event MC weights
        # don't change the reweight meaning (we're comparing data-
        # shifted-explicitly vs data-reweighted-by-model).
        w_unit = np.ones(target.shape[0], dtype=np.float64)
        print("shift reweight closure:")
        plot_mc_shift_reweight(
            flow, target, cond, stats, w_unit,
            output_dir, args.device, args.reweight_chunk,
            shift_factors=tuple(args.shift_factors),
            n_bins=args.n_bins,
        )
        print("smear reweight closure:")
        plot_mc_smear_reweight(
            flow, target, cond, stats, w_unit,
            output_dir, args.device, args.reweight_chunk,
            smear_factors=tuple(args.smear_factors),
            n_bins=args.n_bins,
            n_sigma_pack=int(flow_config.get("n_sigma_pack", 0)),
        )

    # Reweight-head validation (works for both polyhead and MLP
    # arches; the per-event predicted W is computed via an arch-
    # agnostic dispatcher in :func:`_polyhead_pred_W_at`).
    #
    # The shift / smear closure plots ALSO run when no head is
    # present, with the polyhead curve omitted — they still show the
    # raw MC, shifted/smeared MC reference, flow-true K=1 reweight,
    # and (for smear with σ-cond) the σ-cond direct curve, which is a
    # useful flow-only sanity check independent of any head.
    if not args.skip_polyhead_plots:
        polyhead = load_polyhead_from_checkpoint(
            ckpt,
            n_features=int(target_std.shape[1]),
            n_cond=int(cond.shape[1]),
            device=args.device,
        )
        w_unit = np.ones(target.shape[0], dtype=np.float64)
        if polyhead is not None:
            n_polyhead = sum(p.numel() for p in polyhead.parameters())
            head_arch = polyhead.__class__.__name__
            extra = (
                f", n_joint_basis={polyhead.n_joint_basis}"
                if hasattr(polyhead, "n_joint_basis") else ""
            )
            print(
                f"head present ({head_arch}, {n_polyhead:,} params"
                f"{extra}); running validation plots:"
            )
            print("  pred-vs-true scatter:")
            plot_polyhead_pred_vs_true(
                flow, polyhead, target_std, cond, w_unit, stats,
                args, output_dir,
            )
            print("  log-W error histogram:")
            plot_polyhead_logw_error(
                flow, polyhead, target_std, cond, w_unit, stats,
                args, output_dir,
            )
            print("  axis-aligned log-W error histogram:")
            plot_polyhead_axis_logw_error(
                flow, polyhead, target_std, cond, w_unit, stats,
                args, output_dir,
            )
        else:
            print(
                "no head in checkpoint — running flow-only closure "
                "plots (head curve omitted from shift/smear closures, "
                "head-only plots skipped)."
            )
        print(
            "  shift-reweight closure"
            f"{' (flow-only)' if polyhead is None else ''}:"
        )
        plot_polyhead_reweight_closure(
            flow, polyhead, target, target_std, cond, w_unit,
            stats, args, output_dir, mode="shift",
        )
        if polyhead is not None:
            print("  shift-reweight per-bin error ratio:")
            plot_polyhead_reweight_error_ratio(
                flow, polyhead, target, target_std, cond, w_unit,
                stats, args, output_dir, mode="shift",
            )
        print(
            "  smear-reweight closure"
            f"{' (flow-only)' if polyhead is None else ''}:"
        )
        plot_polyhead_reweight_closure(
            flow, polyhead, target, target_std, cond, w_unit,
            stats, args, output_dir, mode="smear",
            n_sigma_pack=int(flow_config.get("n_sigma_pack", 0)),
        )
        if polyhead is not None:
            print("  smear-reweight per-bin error ratio:")
            plot_polyhead_reweight_error_ratio(
                flow, polyhead, target, target_std, cond, w_unit,
                stats, args, output_dir, mode="smear",
            )

    print(f"done. plots in {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
