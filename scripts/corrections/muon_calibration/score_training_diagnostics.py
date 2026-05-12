"""Diagnostic plots for a trained energy-based muon-response score model.

Loads a ``score_model.pt`` (or ``checkpoint.pt``) written by
``train_score_model.py`` and produces a fixed set of diagnostics
focused on the score and its derivative — which is what the model
is actually expected to deliver downstream. Everything the energy
head produces beyond an additive constant per ``c`` is a derivative,
so the plots center on derivatives:

  training_curve.png    train/val DSM losses vs epoch (DSM1 always,
                        DSM2 overlaid when present). Shows convergence
                        behavior.
  log_p_profiles.png    Unnormalized log-density ``-E(y, c, σ)`` as a
                        1D profile along each target axis, for several
                        conditioning points. Others held at target
                        mean. Y-axis shifted so each curve's max is 0.
  log_p_2d.png          2D heatmaps of ``-E`` on each target pair for
                        one canonical context. Visualizes the shape
                        of the peak in each target plane.
  score_profiles.png    3×3 grid of 1D score profiles. Row i, col j:
                        ``s_i(y)`` with axis j varying and others at
                        mean. Diagonal panels (i = j) should cross
                        zero cleanly through the mode; off-diagonals
                        show how the score components couple.
  hessian_profiles.png  3×3 grid of diagonal-Hessian components along
                        each axis. Diagonal panels (i = j) should be
                        negative near the mode (log-concavity) and
                        trend smoothly away from zero in the tails.
  sigma_sweep.png       At a fixed ``y = target_mean`` (presumed
                        mode), plot score components, diagonal
                        Hessian components, and ``trace(H)`` versus σ
                        on a log σ axis. Exposes the σ-conditioning
                        behavior end-to-end; useful to spot training
                        pathologies outside the sampled σ range.
  eig_spectrum.png      Hessian eigenvalues at the presumed mode of
                        each context, as a grouped bar chart. All
                        three eigenvalues should be negative for a
                        well-trained model.
  eig_spectrum_log.png  Same, but plotting |λ| on a log y-axis so the
                        (typically multi-decade) curvature range is
                        legible at a glance. Positive eigenvalues —
                        pathological at a mode — are drawn red.
  score_at_data.png     (optional, --snapshot) Distribution of
                        ``‖s(y, c)‖`` and per-component scores at
                        data points. Sanity check that the model is
                        finite and well-calibrated under data.
  mc_marginals.png      (optional, --snapshot) Pooled 1D MC histograms
                        vs a model marginal reconstructed from the
                        precomputed per-event scores:
                        ``d log p_marg / dy_i = ⟨s_i | y_i⟩`` gives
                        the log-density slope from the binned score,
                        which is cumulatively integrated and
                        exponentiated to produce a normalized density
                        at bin centers. Ratio panel under each axis.
                        No per-point forward density evaluations.
  mc_derivatives.png    (optional, --snapshot) 2×3 grid comparing the
                        first and second derivatives of the 1D log
                        marginal along each target axis between data
                        (finite differences on a weighted histogram)
                        and the model's per-event derivative
                        identities aggregated at data
                        (``E[s_i|y_i]`` for the first derivative,
                        ``Var[s_i|y_i] + E[H_ii|y_i]`` for the
                        second). Residual panels show ``model −
                        data``.
  mc_slices_{pt,lambda,charge}.png
                        (optional, --snapshot) Same MC-vs-model
                        comparison, stratified by log(pt_gen) quartile,
                        lambda_gen quartile, or charge. The model
                        curve per stratum is reconstructed from
                        binned scores *within* that stratum, so every
                        curve is matched to its own data histogram by
                        construction.
  mc_shift_reweight.png (optional, --snapshot) 3×3 grid. For each
                        target axis (column) and each shift magnitude
                        δ = factor·target_std[i] (row), compares the
                        histogram of explicitly shifted MC events
                        ``{y_m + δ·e_i}`` with originals reweighted by
                        the Taylor expansion
                        ``exp(-δ·s_i + 0.5·δ²·H_ii)`` (2nd-order) and
                        its 1st-order truncation. Agreement validates
                        the model's score/Hessian outputs for finite
                        perturbations; divergence at large δ shows
                        where the Taylor series breaks down.
  mc_smear_reweight.png (optional, --snapshot) 3×3 grid. Same layout
                        but the perturbation is a Gaussian smear along
                        one axis. Compares the explicitly smeared
                        histogram to originals reweighted by
                        ``exp(0.5·σ_smear²·(s_i² + H_ii))`` — the
                        identity ``∂²p/∂y² / p = s² + H_ii`` applied
                        to the leading-order convolution. Validates
                        the ``s² + H_ii`` combination pointwise.
  mc_hessian_stratified_{pt,lambda,charge}.png
                        (optional, --snapshot) Stratified analogue of
                        ``mc_derivatives.png``'s second-derivative
                        row. Per-stratum data ``FD² log hist`` vs
                        model ``FD² log p_IS_marginal`` along each
                        target axis, with ``model − data`` residual
                        panel.

Plot shape by construction. The model trained under denoising score
matching approximates ``∇log p_σ`` at whatever σ you query with, and
with the default σ-conditioning mode the same network represents a
whole family. The diagnostics query at ``sigma_inference`` (baked
into the wrapper) unless ``--sigma`` is passed to override; the σ
sweep visualizes the full family so you can confirm that your
inference-σ choice sits comfortably inside the trained range.
"""

import argparse
import os
import sys
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

# Import sibling training module to reuse the model classes and the
# preprocessing dataclass. No module-load side effects (unlike
# flow_training_diagnostics which monkey-patches zuko).
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from train_score_model import (  # noqa: E402
    ACTIVATIONS,
    BRANCHES_PER_MUON,
    EnergyMLP,
    PreprocStats,
    ScoreWrapper,
    WEIGHT_BRANCH,
    compute_targets_and_conditioning,
    load_ntuples,
)


TARGET_NAMES = ["r_kappa", "dlambda", "dphi"]


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--checkpoint",
        required=True,
        help="Path to score_model.pt or checkpoint.pt from "
        "train_score_model.py.",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output directory for PNG diagnostics. "
        "Default: same directory as the --checkpoint file.",
    )
    p.add_argument(
        "--training-log",
        default=None,
        help="Path to training.log for the training-curve plot. "
        "Default: training.log next to the checkpoint; skipped if "
        "that file does not exist.",
    )
    p.add_argument(
        "--snapshot",
        nargs="+",
        default=None,
        help="Optional ROOT snapshot file(s) (same format as the ones "
        "used for training). When provided, produces score_at_data, "
        "mc_marginals, and mc_slices_* plots. Multiple paths are read "
        "together via RDataFrame with ImplicitMT — same multithreaded "
        "loader as train_score_model.py.",
    )
    p.add_argument(
        "--tree",
        default="tree",
        help="TTree name inside the snapshot file(s).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="RDataFrame ImplicitMT threads for snapshot loading. "
        "0 = ROOT's default (all cores); 1 disables MT.",
    )
    p.add_argument(
        "--pt-min",
        type=float,
        default=2.0,
        help="Minimum gen pt (GeV) — same default as the training "
        "script; keeps the MC comparison on the same event set.",
    )
    p.add_argument(
        "--pt-max",
        type=float,
        default=200.0,
        help="Maximum gen pt (GeV), matching the training cut.",
    )
    p.add_argument(
        "--eta-max",
        type=float,
        default=2.4,
        help="Gen/reco |eta| cut, matching the training cut.",
    )
    p.add_argument(
        "--max-data-samples",
        type=int,
        default=200000,
        help="Subsample cap applied after RDataFrame loads. Score "
        "evaluation builds a grad graph per sample, so keep this "
        "moderate.",
    )
    p.add_argument(
        "--n-points",
        type=int,
        default=201,
        help="Grid points per 1D profile axis.",
    )
    p.add_argument(
        "--n-grid-2d",
        type=int,
        default=80,
        help="Grid points per axis in the 2D heatmaps. Cost is N²; "
        "80 is visibly smooth without being slow.",
    )
    p.add_argument(
        "--n-sigma",
        type=int,
        default=16,
        help="Number of σ points in the sigma-sweep plot.",
    )
    p.add_argument(
        "--sigma",
        type=float,
        default=None,
        help="σ value for single-σ plots (profiles, 2D, eigen). "
        "Default: the sigma_inference baked into the checkpoint.",
    )
    p.add_argument(
        "--sigma-grid",
        nargs="+",
        type=float,
        default=None,
        help="Explicit σ values for the sweep plot. Default: "
        "--n-sigma points geometrically spaced between sigma_min "
        "and sigma_max from the checkpoint config.",
    )
    p.add_argument(
        "--reweight-sigmas",
        nargs="+",
        type=float,
        default=None,
        help="σ values at which to evaluate the shift/smear reweight "
        "plots. Each σ gets its own color-coded set of Taylor + "
        "exact reweight curves. Default: for the shift plot "
        "[0.01] + shift_factors; for the smear plot "
        "[0.01] + smear_factors. An explicit list here overrides "
        "both plots.",
    )
    p.add_argument(
        "--n-range-std",
        type=float,
        default=5.0,
        help="Half-range of 1D profiles, in units of per-axis "
        "target standard deviation. 5 covers the peak and the tails.",
    )
    p.add_argument(
        "--device",
        default="cpu",
        help="Device for plotting (cpu is fine; most plots are <10k "
        "total forwards). Set to cuda:0 for the data-score plot only "
        "if the snapshot is large.",
    )
    return p.parse_args()


# -----------------------------------------------------------------------------
# Model loading
# -----------------------------------------------------------------------------

def load_wrapper(
    ckpt_path: str, device: str
) -> Tuple[ScoreWrapper, PreprocStats, dict]:
    """Rebuild a ScoreWrapper from a saved checkpoint.

    Handles both the final ``score_model.pt`` (wrapper_state_dict
    under the top-level key) and ``checkpoint.pt`` (bare inner
    joint-model state_dict with stats in ``stats`` key). If the
    checkpoint includes reweight-head weights (post-refactor), the
    head is built and loaded; otherwise ``wrapper.head`` is ``None``
    and the diagnostic plots use only the base energy + score/Hessian.
    """
    from train_score_model import HeadMLP
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    mc = ckpt["model_config"]

    if "preproc" in ckpt:
        preproc = PreprocStats(**ckpt["preproc"])
    elif "stats" in ckpt and ckpt["stats"] is not None:
        preproc = PreprocStats(**ckpt["stats"])
    else:
        raise RuntimeError(
            f"{ckpt_path}: no 'preproc' or 'stats' key; can't rebuild "
            f"preprocessing stats."
        )

    act_cls = ACTIVATIONS[mc.get("activation", "gelu").lower()]
    energy = EnergyMLP(
        n_features=int(mc["n_features"]),
        n_cond=int(mc["n_cond"]),
        hidden_features=int(mc["hidden_features"]),
        n_hidden_layers=int(mc["n_hidden_layers"]),
        activation=act_cls,
    )
    # Build head only if the config advertises one.
    head = None
    if "head_hidden_features" in mc and "head_n_hidden_layers" in mc:
        head = HeadMLP(
            n_features=int(mc["n_features"]),
            n_cond=int(mc["n_cond"]),
            hidden_features=int(mc["head_hidden_features"]),
            n_hidden_layers=int(mc["head_n_hidden_layers"]),
            activation=act_cls,
        )

    sigma_inf = float(
        ckpt.get("sigma_inference", mc.get("sigma_inference", 0.01))
    )
    wrapper = ScoreWrapper(
        energy=energy,
        head=head,
        target_mean=torch.tensor(preproc.target_mean, dtype=torch.float32),
        target_std=torch.tensor(preproc.target_std, dtype=torch.float32),
        cond_mean=torch.tensor(preproc.cond_mean, dtype=torch.float32),
        cond_std=torch.tensor(preproc.cond_std, dtype=torch.float32),
        sigma_inference=sigma_inf,
    )

    if "wrapper_state_dict" in ckpt:
        # score_model.pt final export; may or may not include head.
        wrapper.load_state_dict(ckpt["wrapper_state_dict"], strict=head is not None)
    elif "state_dict" in ckpt:
        # checkpoint.pt: inner joint-model state → extract energy and
        # head submodules separately (keys prefixed with 'energy.'/'head.').
        inner_state = ckpt["state_dict"]
        energy_state = {
            k[len("energy."):]: v
            for k, v in inner_state.items()
            if k.startswith("energy.")
        }
        energy.load_state_dict(energy_state)
        if head is not None:
            head_state = {
                k[len("head."):]: v
                for k, v in inner_state.items()
                if k.startswith("head.")
            }
            if head_state:
                head.load_state_dict(head_state, strict=False)
    else:
        raise RuntimeError(
            f"{ckpt_path}: no recognized state_dict key "
            f"('wrapper_state_dict' or 'state_dict')."
        )

    wrapper.eval().to(device)
    return wrapper, preproc, mc


# -----------------------------------------------------------------------------
# Contexts
# -----------------------------------------------------------------------------

def build_contexts(preproc: PreprocStats) -> List[Tuple[str, torch.Tensor]]:
    """Return a compact grid of representative conditioning points.

    Covers low/mid/high pt, both charges, and central/barrel/endcap η.
    Each entry is ``(label, c_raw [1, n_cond])`` in the *raw* (pre-
    standardization) ordering matching ``preproc.cond_names``.
    """
    cond_names = list(preproc.cond_names)
    expected = [
        "log_pt_gen", "charge", "lambda_gen",
        "sin_phi_gen", "cos_phi_gen",
    ]
    if cond_names != expected:
        raise RuntimeError(
            f"unexpected cond_names {cond_names}; expected {expected}"
        )

    specs = [
        ("pt=5, q=+, η=0",   5.0,  +1.0, 0.0,  0.0),
        ("pt=20, q=+, η=0",  20.0, +1.0, 0.0,  0.0),
        ("pt=80, q=+, η=0",  80.0, +1.0, 0.0,  0.0),
        ("pt=20, q=-, η=0",  20.0, -1.0, 0.0,  0.0),
        ("pt=20, q=+, η=1.5", 20.0, +1.0, float(np.arctan(np.sinh(1.5))), 0.0),
        ("pt=20, q=+, η=2.2", 20.0, +1.0, float(np.arctan(np.sinh(2.2))), 0.0),
    ]
    ctxs: List[Tuple[str, torch.Tensor]] = []
    for label, pt, q, lam, phi in specs:
        c = np.array(
            [np.log(pt), q, lam, np.sin(phi), np.cos(phi)],
            dtype=np.float32,
        )
        ctxs.append((label, torch.from_numpy(c)[None]))
    return ctxs


def target_grid_along_axis(
    axis_idx: int, preproc: PreprocStats, n_points: int, n_std: float
) -> torch.Tensor:
    """Return a ``[n_points, d]`` grid with axis_idx sweeping ±n_std
    standard deviations around the mean, and other axes fixed at the
    per-axis target mean."""
    mean = np.asarray(preproc.target_mean, dtype=np.float32)
    std = np.asarray(preproc.target_std, dtype=np.float32)
    lo = mean[axis_idx] - n_std * std[axis_idx]
    hi = mean[axis_idx] + n_std * std[axis_idx]
    ys = np.tile(mean[None, :], (n_points, 1))
    ys[:, axis_idx] = np.linspace(lo, hi, n_points)
    return torch.from_numpy(ys)


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

def _save(fig, path: str):
    """Save ``fig`` as both PNG (quick inspection) and PDF
    (publication-quality). ``path`` should end in ``.png``; the PDF
    is written alongside with the same basename."""
    base, ext = os.path.splitext(path)
    png_path = path if ext.lower() == ".png" else base + ".png"
    pdf_path = base + ".pdf"
    fig.savefig(png_path, dpi=110, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png_path} and {pdf_path}")


def _chunked_score(wrapper, y, c, sigma, chunk=4096):
    """Evaluate scores in CPU-friendly chunks, returning a detached
    ``[N, d]`` tensor. Each chunk rebuilds the autograd graph afresh
    (score() uses a leaf tensor internally)."""
    out = []
    for s in range(0, y.shape[0], chunk):
        e = min(s + chunk, y.shape[0])
        out.append(
            wrapper.score(y[s:e], c[s:e], sigma=sigma).detach()
        )
    return torch.cat(out, dim=0)


def _chunked_unnorm_log_density(wrapper, y, c, sigma, chunk=16384):
    """Pure-forward ``wrapper.unnormalized_log_density`` in chunks.

    No autograd — a batched inference loop. Used by the shift / smear
    closure tests to evaluate ``log p(y' | c)`` at perturbed positions
    so the *exact* model reweight ``exp(lp(y') − lp(y))`` can be
    compared against the linear Taylor truncations.
    """
    out = []
    with torch.no_grad():
        for s in range(0, y.shape[0], chunk):
            e = min(s + chunk, y.shape[0])
            out.append(
                wrapper.unnormalized_log_density(
                    y[s:e], c[s:e], sigma=sigma,
                ).detach()
            )
    return torch.cat(out, dim=0).numpy()


def _chunked_hessian(wrapper, y, c, sigma, chunk=1024):
    out = []
    for s in range(0, y.shape[0], chunk):
        e = min(s + chunk, y.shape[0])
        out.append(
            wrapper.hessian(y[s:e], c[s:e], sigma=sigma).detach()
        )
    return torch.cat(out, dim=0)


# -----------------------------------------------------------------------------
# Training-curve plot
# -----------------------------------------------------------------------------

def parse_training_log(path: str):
    """Extract per-epoch metrics from the log file written by train().

    Recognizes all of the historical schema variants:
      * DSM1-only:  ``train_dsm ... val_dsm ...``
      * DSM1+DSM2:  ``train_dsm1 ... train_dsm2 ... val_dsm1 ... val_dsm2``
      * Joint head: ``train_dsm ... train_shift ... train_smear ...
                     val_dsm ... val_shift ... val_smear ...``
    """
    recs = []
    known_keys = {
        "train_dsm", "train_dsm1", "train_dsm2",
        "val_dsm", "val_dsm1", "val_dsm2",
        "train_shift", "val_shift",
        "train_smear", "val_smear",
    }
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln.startswith("epoch "):
                continue
            toks = ln.split()
            try:
                ep = int(toks[1])
            except (ValueError, IndexError):
                continue
            d = {"epoch": ep}
            for i, t in enumerate(toks):
                if t in known_keys and i + 1 < len(toks):
                    try:
                        d[t] = float(toks[i + 1])
                    except ValueError:
                        pass
            if len(d) > 1:
                recs.append(d)
    return recs


def plot_training_curve(log_path: str, out_dir: str):
    recs = parse_training_log(log_path)
    if not recs:
        print(f"[skip] training log at {log_path} has no parseable lines")
        return
    epochs = [r["epoch"] for r in recs]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    style = {
        "train_dsm":   dict(color="C0", ls="-",  marker="o", ms=3, label="train DSM"),
        "val_dsm":     dict(color="C0", ls="--", marker="o", ms=3, label="val DSM"),
        "train_dsm1":  dict(color="C0", ls="-",  marker="o", ms=3, label="train DSM1"),
        "val_dsm1":    dict(color="C0", ls="--", marker="o", ms=3, label="val DSM1"),
        "train_dsm2":  dict(color="C3", ls="-",  marker="s", ms=3, label="train DSM2"),
        "val_dsm2":    dict(color="C3", ls="--", marker="s", ms=3, label="val DSM2"),
        "train_shift": dict(color="C2", ls="-",  marker="^", ms=3, label="train shift"),
        "val_shift":   dict(color="C2", ls="--", marker="^", ms=3, label="val shift"),
        "train_smear": dict(color="C4", ls="-",  marker="v", ms=3, label="train smear"),
        "val_smear":   dict(color="C4", ls="--", marker="v", ms=3, label="val smear"),
    }
    secondary_keys = {
        "train_dsm2", "val_dsm2",
        "train_shift", "val_shift",
        "train_smear", "val_smear",
    }
    has_secondary = any(
        any(k in r for r in recs) for k in secondary_keys
    )
    if has_secondary:
        ax2 = ax.twinx()
    else:
        ax2 = None
    for key, st in style.items():
        if not any(key in r for r in recs):
            continue
        vals = [r.get(key, np.nan) for r in recs]
        target_ax = ax2 if (key in secondary_keys and ax2 is not None) else ax
        target_ax.plot(epochs, vals, **st)
    ax.set_xlabel("epoch")
    ax.set_ylabel("DSM loss")
    ax.grid(alpha=0.3)
    if ax2 is not None:
        ax2.set_ylabel("head / DSM2 loss")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)
    else:
        ax.legend(loc="best")
    ax.set_title("training curve")
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "training_curve.png"))


# -----------------------------------------------------------------------------
# Profile plots
# -----------------------------------------------------------------------------

def plot_log_p_profiles(
    wrapper, contexts, preproc, out_dir, sigma, n_points, n_std
):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax_idx, ax in enumerate(axes):
        y_grid = target_grid_along_axis(
            ax_idx, preproc, n_points, n_std
        )
        for label, c_raw in contexts:
            c = c_raw.expand(n_points, -1)
            with torch.no_grad():
                lp = wrapper.unnormalized_log_density(
                    y_grid, c, sigma=sigma
                ).detach()
            lp = lp - lp.max()
            ax.plot(
                y_grid[:, ax_idx].numpy(), lp.numpy(),
                label=label, alpha=0.85,
            )
        ax.set_xlabel(TARGET_NAMES[ax_idx])
        ax.set_ylabel("log p − max (others at mean)")
        ax.grid(alpha=0.3)
        ax.axvline(preproc.target_mean[ax_idx], color="k",
                   ls=":", lw=0.6)
    axes[0].legend(fontsize=7, loc="lower center")
    fig.suptitle(f"unnormalized log-density 1D profiles, σ = {sigma:.4g}")
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "log_p_profiles.png"))


def plot_log_p_2d(
    wrapper, contexts, preproc, out_dir, sigma, n_grid
):
    # Pick the canonical central context (index 1: pt=20, q=+, η=0).
    label, c_raw = contexts[1] if len(contexts) > 1 else contexts[0]
    mean = np.asarray(preproc.target_mean, dtype=np.float32)
    std = np.asarray(preproc.target_std, dtype=np.float32)
    pairs = [(0, 1), (0, 2), (1, 2)]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (i, j) in zip(axes, pairs):
        lo_i, hi_i = mean[i] - 4 * std[i], mean[i] + 4 * std[i]
        lo_j, hi_j = mean[j] - 4 * std[j], mean[j] + 4 * std[j]
        yi = np.linspace(lo_i, hi_i, n_grid)
        yj = np.linspace(lo_j, hi_j, n_grid)
        YI, YJ = np.meshgrid(yi, yj, indexing="xy")
        Y = np.tile(mean[None, None, :], (n_grid, n_grid, 1))
        Y[..., i] = YI
        Y[..., j] = YJ
        Y_flat = torch.from_numpy(Y.reshape(-1, 3).astype(np.float32))
        c_flat = c_raw.expand(Y_flat.shape[0], -1)
        with torch.no_grad():
            lp = wrapper.unnormalized_log_density(
                Y_flat, c_flat, sigma=sigma
            ).detach()
        lp = lp.numpy().reshape(n_grid, n_grid) - lp.max().item()
        im = ax.imshow(
            lp, origin="lower",
            extent=[lo_i, hi_i, lo_j, hi_j],
            aspect="auto", cmap="viridis",
            vmin=-15, vmax=0,
        )
        ax.set_xlabel(TARGET_NAMES[i])
        ax.set_ylabel(TARGET_NAMES[j])
        fig.colorbar(im, ax=ax, label="log p − max")
    fig.suptitle(
        f"log-density 2D slices at [{label}], σ = {sigma:.4g} "
        f"(remaining axis at its mean)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "log_p_2d.png"))


def plot_score_profiles(
    wrapper, contexts, preproc, out_dir, sigma, n_points, n_std
):
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex="col")
    for col, ax_idx in enumerate(range(3)):
        y_grid = target_grid_along_axis(
            ax_idx, preproc, n_points, n_std
        )
        for label, c_raw in contexts:
            c = c_raw.expand(n_points, -1)
            s = _chunked_score(wrapper, y_grid, c, sigma)
            for row in range(3):
                axes[row, col].plot(
                    y_grid[:, ax_idx].numpy(),
                    s[:, row].numpy(),
                    label=label, alpha=0.75,
                )
        axes[2, col].set_xlabel(TARGET_NAMES[ax_idx])
        for row in range(3):
            ax = axes[row, col]
            ax.axhline(0, color="k", lw=0.5)
            ax.axvline(preproc.target_mean[ax_idx], color="k",
                       ls=":", lw=0.5)
            ax.grid(alpha=0.3)
            if col == 0:
                ax.set_ylabel(f"s_{TARGET_NAMES[row]}")
    axes[0, 0].legend(fontsize=6, loc="best")
    fig.suptitle(
        f"score profiles at σ = {sigma:.4g}  "
        f"(rows: score component; columns: target axis varying)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "score_profiles.png"))


def plot_hessian_profiles(
    wrapper, contexts, preproc, out_dir, sigma, n_points, n_std
):
    """Plot diagonal Hessian components along each target axis.

    Only diagonals: off-diagonals are much smaller and hard to read
    in a shared panel. Each panel shows H_{ii} as axis col varies
    while other axes sit at the target mean.
    """
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex="col")
    for col, ax_idx in enumerate(range(3)):
        y_grid = target_grid_along_axis(
            ax_idx, preproc, n_points, n_std
        )
        for label, c_raw in contexts:
            c = c_raw.expand(n_points, -1)
            H = _chunked_hessian(wrapper, y_grid, c, sigma)
            for row in range(3):
                axes[row, col].plot(
                    y_grid[:, ax_idx].numpy(),
                    H[:, row, row].numpy(),
                    label=label, alpha=0.75,
                )
        axes[2, col].set_xlabel(TARGET_NAMES[ax_idx])
        for row in range(3):
            ax = axes[row, col]
            ax.axhline(0, color="k", lw=0.5)
            ax.axvline(preproc.target_mean[ax_idx], color="k",
                       ls=":", lw=0.5)
            ax.grid(alpha=0.3)
            if col == 0:
                ax.set_ylabel(f"H_{TARGET_NAMES[row]}{TARGET_NAMES[row]}")
    axes[0, 0].legend(fontsize=6, loc="best")
    fig.suptitle(
        f"Hessian diagonal profiles at σ = {sigma:.4g}  "
        f"(rows: diag element; columns: target axis varying)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "hessian_profiles.png"))


# -----------------------------------------------------------------------------
# σ sweep
# -----------------------------------------------------------------------------

def plot_sigma_sweep(
    wrapper, contexts, preproc, out_dir, sigma_grid
):
    """At y = target_mean (a good proxy for the mode for most c),
    trace the score, the Hessian diagonal, and ``trace(H)`` across
    the σ grid. Linear σ-axis on a log scale because σ is sampled
    log-uniform during training."""
    target_mean = torch.tensor(preproc.target_mean, dtype=torch.float32)
    y = target_mean[None, :]
    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    for (label, c_raw) in contexts:
        scores = np.empty((len(sigma_grid), 3))
        diag_H = np.empty((len(sigma_grid), 3))
        trace_H = np.empty(len(sigma_grid))
        for k, sig in enumerate(sigma_grid):
            s = wrapper.score(y, c_raw, sigma=float(sig)).detach().squeeze(0)
            H = wrapper.hessian(y, c_raw, sigma=float(sig)).detach().squeeze(0)
            scores[k] = s.numpy()
            diag_H[k] = torch.diag(H).numpy()
            trace_H[k] = torch.trace(H).item()
        for comp in range(3):
            axes[0].plot(
                sigma_grid, scores[:, comp],
                label=f"{label} / s_{TARGET_NAMES[comp]}", alpha=0.6,
            )
            axes[1].plot(
                sigma_grid, diag_H[:, comp],
                label=f"{label} / H_{TARGET_NAMES[comp]}{TARGET_NAMES[comp]}",
                alpha=0.6,
            )
        axes[2].plot(sigma_grid, trace_H, label=label,
                     marker="o", ms=3, alpha=0.8)
    axes[0].set_ylabel("score components at y=mean")
    axes[1].set_ylabel("H_ii diagonal components")
    axes[2].set_ylabel("trace(H)")
    axes[2].set_xlabel("σ")
    for ax in axes:
        ax.set_xscale("log")
        ax.axhline(0, color="k", lw=0.5)
        ax.grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=5, ncol=3, loc="upper right")
    axes[2].legend(fontsize=7, loc="best")
    fig.suptitle("σ sweep at y = target_mean")
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "sigma_sweep.png"))


# -----------------------------------------------------------------------------
# Eigenvalue spectrum
# -----------------------------------------------------------------------------

def plot_eigenvalues(
    wrapper, contexts, preproc, out_dir, sigma
):
    target_mean = torch.tensor(preproc.target_mean, dtype=torch.float32)
    y = target_mean[None, :]
    eigs_rows = []
    labels = []
    for label, c_raw in contexts:
        H = wrapper.hessian(y, c_raw, sigma=sigma).detach().squeeze(0).numpy()
        H = 0.5 * (H + H.T)  # symmetrize away the ~1e-9 asymmetry
        eigs = np.linalg.eigvalsh(H)
        eigs_rows.append(eigs)
        labels.append(label)
    eigs = np.asarray(eigs_rows)  # [N_ctx, 3]
    x = np.arange(len(labels))
    width = 0.26

    # Linear y-axis: signed eigenvalues, ref line at zero.
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        ax.bar(
            x + (i - 1) * width, eigs[:, i], width,
            label=f"λ_{i+1} (ordered)",
        )
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(f"Hessian eigenvalues at y=mean, σ={sigma:.4g}")
    ax.set_title(
        "Hessian eigenvalues (all should be < 0 at a mode)"
    )
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "eig_spectrum.png"))

    # Log y-axis companion: plots |λ| so the full dynamic range
    # (often several decades between r_kappa and dlambda/dphi
    # curvatures) is legible. Color by sign — red bars flag any
    # spurious positive eigenvalues, which should not be present at
    # a true mode.
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [
        ("C0", "C1", "C2"),  # λ_1, λ_2, λ_3 when negative (OK)
    ][0]
    for i in range(3):
        abs_vals = np.abs(eigs[:, i])
        bar_colors = [
            colors[i] if e < 0 else "C3" for e in eigs[:, i]
        ]
        ax.bar(
            x + (i - 1) * width, abs_vals, width,
            color=bar_colors, label=f"λ_{i+1} (ordered)",
        )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(f"|Hessian eigenvalues| at y=mean, σ={sigma:.4g}")
    ax.set_title(
        "Hessian eigenvalue magnitudes (red = positive = pathological)"
    )
    ax.legend()
    ax.grid(alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "eig_spectrum_log.png"))


# -----------------------------------------------------------------------------
# MC data loader + shared helpers
# -----------------------------------------------------------------------------

def load_snapshot_data(
    snapshot_paths, preproc, max_samples,
    tree="tree", threads=0,
    pt_min=2.0, pt_max=200.0, eta_max=2.4,
):
    """Load snapshot file(s) via the training script's ``load_ntuples``
    (RDataFrame + ImplicitMT + C++ quality filter), pool mu+/mu-, and
    build the arrays needed by every MC-comparison plot.

    Returns
    -------
    target : np.ndarray [N, 3]
        Raw-coordinate targets (r_kappa, dlambda, dphi).
    c_raw : np.ndarray [N, n_cond]
        Raw conditioning in the order of ``preproc.cond_names``.
    w : np.ndarray [N]
        Per-muon event weights.

    Returns ``None`` on ROOT import failure.
    """
    try:
        import ROOT  # noqa: F401
    except ImportError:
        print("[skip] ROOT not available; no MC comparisons possible")
        return None

    paths = (
        [snapshot_paths]
        if isinstance(snapshot_paths, str)
        else list(snapshot_paths)
    )

    pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q, w = load_ntuples(
        paths=paths,
        tree_name=tree,
        max_muons=max_samples,
        pt_min=pt_min,
        pt_max=pt_max,
        eta_max=eta_max,
        threads=threads,
    )

    target, cond_raw = compute_targets_and_conditioning(
        pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q
    )
    c_cols = [cond_raw[name] for name in preproc.cond_names]
    c_raw = np.stack(c_cols, axis=1).astype(np.float32)
    return target.astype(np.float32), c_raw, w.astype(np.float64)


def _data_1d_histogram(target, weights, axis_idx, preproc, n_bins, n_std):
    """Weighted, range-normalized 1D histogram of ``target[:, axis_idx]``."""
    mean = preproc.target_mean[axis_idx]
    std = preproc.target_std[axis_idx]
    lo, hi = mean - n_std * std, mean + n_std * std
    hist, edges = np.histogram(
        target[:, axis_idx], bins=n_bins, range=(lo, hi),
        weights=weights, density=True,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    return hist, edges, centers


def _reconstruct_marginal_from_score(centers, mean_s):
    """Reconstruct the normalized 1D marginal ``p_marg(y_i)`` by
    integrating the binned marginal score ``E[s_i | y_i] ≈ mean_s``.

    Uses the identity ``d log p_marg / dy_i = E[s_i | y_i]`` — the
    cumulative trapezoidal integral of mean_s is ``log p_marg`` up to
    an additive constant. We exponentiate, renormalize to unit area,
    and return the density sampled at ``centers``.

    NaN bins (empty strata) in ``mean_s`` are treated as zero-slope
    to avoid poisoning the cumulative sum; the surrounding non-NaN
    bins still carry the correct shape. The computation is a couple
    of O(n_bins) numpy ops — no per-point network evaluations.
    """
    mean_s = np.asarray(mean_s, dtype=np.float64)
    s = np.where(np.isnan(mean_s), 0.0, mean_s)
    dx = np.diff(centers)
    mid_s = 0.5 * (s[:-1] + s[1:])
    log_p = np.concatenate([[0.0], np.cumsum(mid_s * dx)])
    log_p -= log_p.max()  # numerical stability before exp
    p = np.exp(log_p)
    trap = getattr(np, "trapezoid", np.trapz)
    norm = trap(p, centers)
    return p / max(norm, 1e-30)


# -----------------------------------------------------------------------------
# MC comparison plots (optional)
# -----------------------------------------------------------------------------

def plot_score_at_data(
    wrapper, target, c_raw, w, preproc, out_dir, sigma,
):
    """Distribution of ‖s(y, c)‖ and per-component scores at data
    points. E[s] ≈ 0 under p (Stein's identity), so the per-component
    histograms should be tight around zero."""
    y_t = torch.from_numpy(target).contiguous()
    c_t = torch.from_numpy(c_raw).contiguous()
    s = _chunked_score(wrapper, y_t, c_t, sigma, chunk=8192).numpy()
    norms = np.linalg.norm(s, axis=1)

    fig, axes = plt.subplots(1, 4, figsize=(17, 4))
    axes[0].hist(norms, bins=120, weights=w)
    axes[0].set_xlabel("‖s(y, c)‖")
    axes[0].set_ylabel("weighted count")
    axes[0].set_yscale("log")
    axes[0].set_title("score norm at data")
    axes[0].grid(alpha=0.3)
    for i, name in enumerate(TARGET_NAMES):
        w_mean = float(np.average(s[:, i], weights=w))
        axes[i + 1].hist(
            s[:, i], bins=120, weights=w, density=True, alpha=0.7
        )
        axes[i + 1].axvline(0, color="k", lw=0.5)
        axes[i + 1].axvline(
            w_mean, color="C3", ls="--", lw=0.8,
            label=f"⟨s⟩ = {w_mean:+.2e}",
        )
        axes[i + 1].set_xlabel(f"s_{name} at data")
        axes[i + 1].grid(alpha=0.3)
        axes[i + 1].legend(fontsize=8)
    fig.suptitle(
        f"score distribution at data (σ={sigma:.4g}); "
        f"Stein's identity says E[s] ≈ 0 under p."
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "score_at_data.png"))


def plot_mc_marginals(
    target, c_raw, w, s_data, preproc, out_dir, sigma,
    n_bins=80, n_std=5.0,
):
    """Pooled data histograms vs model 1D marginals reconstructed
    from the model's per-event scores.

    For each target axis i,
      - data: histogram of y_i across all events;
      - model: ``p_marg(y_i)`` obtained by binning ``s_model_i(y_m,
        c_m)`` by y_i, averaging, and integrating
        (``d log p_marg / dy_i = E[s_i | y_i]``).
    No per-point forward density evaluations — the only network work
    was the one-pass ``s_data`` computation already in hand.
    """
    fig, axes = plt.subplots(
        2, 3, figsize=(15, 6),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex="col",
    )
    for i in range(3):
        ax, ax_r = axes[0, i], axes[1, i]
        hist, edges, centers = _data_1d_histogram(
            target, w, i, preproc, n_bins, n_std
        )
        lo, hi = edges[0], edges[-1]
        in_range = (target[:, i] >= lo) & (target[:, i] < hi)
        bin_idx = np.clip(
            np.digitize(target[in_range, i], edges) - 1, 0, n_bins - 1,
        )
        mean_s = _weighted_bin_mean(
            s_data[in_range, i], bin_idx, n_bins, w[in_range],
        )
        p_model = _reconstruct_marginal_from_score(centers, mean_s)

        ax.stairs(
            hist, edges, label="data", color="C3",
            fill=True, alpha=0.35,
        )
        ax.plot(
            centers, p_model,
            label=f"model (score-integrated, N={target.shape[0]})",
            color="C0", lw=1.5,
        )
        ax.set_ylabel("density")
        ax.grid(alpha=0.3)
        ax.set_title(TARGET_NAMES[i])
        if i == 0:
            ax.legend(fontsize=8, loc="best")

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(p_model > 0, hist / p_model, np.nan)
        ax_r.step(edges[:-1], ratio, where="post", color="C3")
        ax_r.axhline(1.0, color="k", lw=0.5)
        ax_r.set_ylabel("data / model")
        ax_r.set_xlabel(TARGET_NAMES[i])
        ax_r.set_ylim(0, 2)
        ax_r.grid(alpha=0.3)
    fig.suptitle(
        f"MC vs model 1D marginals (score-integrated), σ = {sigma:.4g}"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "mc_marginals.png"))


def _stratify_groups(c_raw, w, preproc, stratify_by):
    """Return ``list[(label, mask)]`` for the requested stratification.

    ``stratify_by`` is one of 'pt', 'lambda', 'charge'. For pt and
    lambda, produces weighted quartiles. For charge, produces two ±
    groups. Masks are boolean arrays over the data rows.
    """
    cond_map = {
        "pt": ("log_pt_gen", "log(pt_gen)"),
        "lambda": ("lambda_gen", "lambda_gen"),
        "charge": ("charge", "charge"),
    }
    col_name, human = cond_map[stratify_by]
    col_idx = preproc.cond_names.index(col_name)
    var = c_raw[:, col_idx]

    if stratify_by == "charge":
        return [("q = +", var > 0), ("q = −", var < 0)]

    # Weighted quartile cut points.
    order = np.argsort(var)
    cumw = np.cumsum(w[order])
    total_w = cumw[-1]
    cuts = [0.25, 0.5, 0.75]
    edges = [
        var[order[np.searchsorted(cumw, q * total_w)]] for q in cuts
    ]
    if stratify_by == "pt":
        disp = [f"pt<{np.exp(edges[0]):.1f}",
                f"{np.exp(edges[0]):.1f}<pt<{np.exp(edges[1]):.1f}",
                f"{np.exp(edges[1]):.1f}<pt<{np.exp(edges[2]):.1f}",
                f"pt>{np.exp(edges[2]):.1f}"]
    else:
        disp = [f"{human}<{edges[0]:+.2f}",
                f"{edges[0]:+.2f}<{human}<{edges[1]:+.2f}",
                f"{edges[1]:+.2f}<{human}<{edges[2]:+.2f}",
                f"{human}>{edges[2]:+.2f}"]
    return [
        (disp[0], var < edges[0]),
        (disp[1], (var >= edges[0]) & (var < edges[1])),
        (disp[2], (var >= edges[1]) & (var < edges[2])),
        (disp[3], var >= edges[2]),
    ]


def _fd_log(y, vals, floor_factor=1e3):
    """Return first and second derivatives of ``log(vals)`` at ``y`` via
    central differences (``np.gradient``). Points whose ``vals`` is
    below ``max(vals)/floor_factor`` are NaN-masked so derivatives in
    empty-tail regions don't blow up the y-axis."""
    eps = np.finfo(np.float64).tiny
    vals = np.asarray(vals, dtype=np.float64)
    log_v = np.log(np.maximum(vals, eps))
    d1 = np.gradient(log_v, y)
    d2 = np.gradient(d1, y)
    mask = vals < (vals.max() / floor_factor)
    d1[mask] = np.nan
    d2[mask] = np.nan
    return d1, d2


def _weighted_bin_mean(values, bin_idx, n_bins, weights):
    """``values`` mean in each of ``n_bins`` bins, weighted by ``weights``.

    Uses ``np.add.at`` for the weighted accumulation (handles repeated
    bin indices correctly). NaN for empty bins.
    """
    sums = np.zeros(n_bins, dtype=np.float64)
    wsums = np.zeros(n_bins, dtype=np.float64)
    np.add.at(sums, bin_idx, values.astype(np.float64) * weights)
    np.add.at(wsums, bin_idx, weights)
    return np.where(wsums > 0, sums / np.maximum(wsums, 1e-30), np.nan)


def plot_mc_derivatives(
    target, c_raw, w, s_data, H_diag, preproc, out_dir, sigma,
    n_bins=50, n_std=5.0,
):
    """Data-vs-model comparison of first and second derivatives of
    the 1D log marginal along each target axis.

    Two curves per panel (no more grid-based IS marginal FD):

      * ``data FD``: finite difference of ``log(hist(y_i))`` on a
        weighted histogram.
      * ``model binned at data``:
          - first derivative → ``⟨s_i(y_m, c_m) | y_i⟩``, which
            equals ``d log p_marg / dy_i`` in expectation.
          - second derivative → ``⟨H_ii + s_i²⟩ − ⟨s_i⟩² | y_i``, the
            conditional-expectation identity for
            ``d²log p_marg / dy_i²``.

    Both model curves are aggregated from the precomputed ``s_data``
    and ``H_diag`` arrays — no forward density evaluations here.
    """
    n_events = target.shape[0]

    fig, axes = plt.subplots(
        4, 3, figsize=(16, 12), sharex="col",
        gridspec_kw={"height_ratios": [3, 1, 3, 1], "hspace": 0.05},
    )
    for i in range(3):
        mean_i = preproc.target_mean[i]
        std_i = preproc.target_std[i]
        lo, hi = mean_i - n_std * std_i, mean_i + n_std * std_i

        # --- data histogram + FD on log hist -----------------------
        hist, edges = np.histogram(
            target[:, i], bins=n_bins, range=(lo, hi),
            weights=w, density=True,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])
        d1_data, d2_data = _fd_log(centers, hist)

        # --- binned model derivatives at data events --------------
        in_range = (target[:, i] >= lo) & (target[:, i] < hi)
        bin_idx = np.clip(
            np.digitize(target[in_range, i], edges) - 1,
            0, n_bins - 1,
        )
        s_i = s_data[in_range, i]
        H_ii = H_diag[in_range, i]
        w_in = w[in_range]
        mean_s = _weighted_bin_mean(s_i, bin_idx, n_bins, w_in)
        mean_s_sq = _weighted_bin_mean(s_i * s_i, bin_idx, n_bins, w_in)
        mean_H = _weighted_bin_mean(H_ii, bin_idx, n_bins, w_in)
        # Identity: d²log p_marg/dy_i² = Var[s_i|y_i] + E[H_ii|y_i]
        #                              = E[s_i²|y_i] - E[s_i|y_i]² + E[H_ii|y_i]
        d2_binned = mean_s_sq - mean_s ** 2 + mean_H

        r1_binned = mean_s - d1_data
        r2_binned = d2_binned - d2_data

        ax1, axr1 = axes[0, i], axes[1, i]
        ax2, axr2 = axes[2, i], axes[3, i]

        # --- 1st-derivative panel ----------------------------------
        ax1.plot(
            centers, d1_data, "o", color="C3", ms=3,
            label="data (FD)",
        )
        ax1.plot(
            centers, mean_s, "s-", color="C2", ms=3, alpha=0.85,
            label="model ⟨s_i|y_i⟩ at data",
        )
        ax1.axhline(0, color="k", lw=0.5)
        ax1.set_ylabel(f"∂/∂{TARGET_NAMES[i]} log p_marg")
        ax1.grid(alpha=0.3)
        if i == 0:
            ax1.legend(fontsize=7, loc="best")
        ax1.set_title(TARGET_NAMES[i])

        axr1.axhline(0.0, color="k", lw=0.5)
        axr1.plot(
            centers, r1_binned, "s-", color="C2", ms=3, alpha=0.85,
        )
        axr1.set_ylabel("model − data")
        axr1.grid(alpha=0.3)

        # --- 2nd-derivative panel ----------------------------------
        ax2.plot(
            centers, d2_data, "o", color="C3", ms=3,
            label="data (FD²)",
        )
        ax2.plot(
            centers, d2_binned, "s-", color="C2", ms=3, alpha=0.85,
            label="model ⟨Var[s]+⟨H⟩|y⟩ at data",
        )
        ax2.axhline(0, color="k", lw=0.5)
        ax2.set_ylabel(f"∂²/∂{TARGET_NAMES[i]}² log p_marg")
        ax2.grid(alpha=0.3)
        if i == 0:
            ax2.legend(fontsize=7, loc="best")

        axr2.axhline(0.0, color="k", lw=0.5)
        axr2.plot(
            centers, r2_binned, "s-", color="C2", ms=3, alpha=0.85,
        )
        axr2.set_ylabel("model − data")
        axr2.set_xlabel(TARGET_NAMES[i])
        axr2.grid(alpha=0.3)

    fig.suptitle(
        f"1D marginal log-density derivatives: data vs model "
        f"(σ = {sigma:.4g}, {n_bins} FD bins, {n_events} events); "
        f"residuals below each panel"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "mc_derivatives.png"))


def _safe_ratio(numer, denom, rel_floor=1e-3):
    """Element-wise ``numer / denom`` with small-denom bins masked.

    Bins where ``denom`` is below ``rel_floor · max(denom)`` are set
    to NaN so they don't blow up the ratio plot.
    """
    numer = np.asarray(numer, dtype=np.float64)
    denom = np.asarray(denom, dtype=np.float64)
    thr = rel_floor * np.nanmax(denom) if denom.size else 0.0
    safe = np.where(denom > thr, denom, 1.0)
    return np.where(denom > thr, numer / safe, np.nan)


def plot_mc_shift_reweight(
    wrapper, target, c_raw, w, preproc, out_dir, sigmas=None,
    shift_factors=(0.05, 0.1, 0.2),
    n_bins=80, n_std=5.0,
):
    """Validate the model's score and Hessian by shifting MC events.

    For each target axis ``i`` and each magnitude ``δ = factor ·
    target_std[i]`` we produce histograms that all estimate the same
    target density — ``p`` translated by ``δ·e_i``:

      * **data** (red shaded): events ``{y_m + δ·e_i}`` — explicit shift.
      * **unshifted MC** (gray dotted): ``{y_m}`` — reference baseline.
      * For each σ in ``sigmas`` (color-coded):
          - **1st-order reweight** (dashed): weights
            ``w_m · (1 − δ·s_i(σ))`` — strict linear Taylor of
            ``p_σ(y − δ·e_i)/p_σ(y)``.
          - **2nd-order reweight** (dot-dashed): weights
            ``w_m · [1 − δ·s_i(σ) + 0.5·δ²·(s_i(σ)² + H_ii(σ))]``.
          - **exact** (solid): weights
            ``w_m · exp(lp_σ(y_m − δ·e_i) − lp_σ(y_m))``.

    Each panel has a ratio sub-panel below it plotting the reweight
    curves divided by the explicit-shift data histogram; ratios of 1
    indicate a matching reweight.

    If ``sigmas`` is None, defaults to ``[0.01] + list(shift_factors)``
    — the minimal model σ plus values matching the δ scales tested.
    """
    if sigmas is None:
        sigmas = [0.01] + list(shift_factors)
    std = np.asarray(preproc.target_std, dtype=np.float64)
    mean = np.asarray(preproc.target_mean, dtype=np.float64)
    n_mag = len(shift_factors)
    sigmas = list(sigmas)

    target_t = torch.from_numpy(target.astype(np.float32))
    c_t = torch.from_numpy(c_raw.astype(np.float32))

    # Per-σ: score, diag-Hessian, lp at originals. Each takes one
    # chunked forward pass over all events; lp at the shifted points
    # is done inside the panel loop (it depends on δ as well).
    print(
        f"  shift reweight: precomputing s/H/lp0 at "
        f"{len(sigmas)} σ values × {target.shape[0]} events"
    )
    per_sigma = []
    for sig in sigmas:
        s_sig = _chunked_score(
            wrapper, target_t, c_t, sig, chunk=8192,
        ).numpy()
        H_full = _chunked_hessian(
            wrapper, target_t, c_t, sig, chunk=1024,
        ).numpy()
        H_sig = np.stack([H_full[:, k, k] for k in range(3)], axis=1)
        del H_full
        lp0_sig = _chunked_unnorm_log_density(
            wrapper, target_t, c_t, sig,
        )
        per_sigma.append((sig, s_sig, H_sig, lp0_sig))

    cmap = plt.get_cmap("viridis")
    if len(sigmas) == 1:
        sig_colors = [cmap(0.5)]
    else:
        sig_colors = [
            cmap(i / (len(sigmas) - 1)) for i in range(len(sigmas))
        ]

    # 2 rows per magnitude: main panel + ratio panel.
    fig, axes = plt.subplots(
        2 * n_mag, 3, figsize=(15, 4.0 * n_mag),
        gridspec_kw={"height_ratios": [3, 1] * n_mag},
        sharex="col",
    )
    if n_mag == 1 and axes.ndim == 1:
        axes = axes[None, :]
    for col in range(3):
        lo, hi = mean[col] - n_std * std[col], mean[col] + n_std * std[col]
        for row, factor in enumerate(shift_factors):
            ax = axes[2 * row, col]
            ax_r = axes[2 * row + 1, col]
            delta = factor * std[col]

            # Explicit shift + unshifted reference (σ-independent).
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

            # Shifted positions shared across σ for the exact reweight.
            y_shift_np = target.copy()
            y_shift_np[:, col] -= delta
            y_shift_t = torch.from_numpy(y_shift_np.astype(np.float32))

            show_labels = (row == 0 and col == 0)
            for (sig, s_sig, H_sig, lp0_sig), sc in zip(
                per_sigma, sig_colors,
            ):
                s_col = s_sig[:, col]
                H_col = H_sig[:, col]

                w1_factor = 1.0 - delta * s_col
                w2_factor = (
                    w1_factor + 0.5 * delta * delta * (s_col ** 2 + H_col)
                )
                hist_r1, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w * w1_factor, density=True,
                )
                hist_r2, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w * w2_factor, density=True,
                )

                lp_shift = _chunked_unnorm_log_density(
                    wrapper, y_shift_t, c_t, sig,
                )
                w_exact = w * np.exp(lp_shift - lp0_sig)
                hist_ex, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w_exact, density=True,
                )

                # zorder: exact solid keeps default (=2); Taylor
                # 1st/2nd dashed / dash-dot overlays use zorder=3
                # so they remain visible above the exact curve.
                ax.plot(
                    centers, hist_r1, color=sc, lw=1.0, linestyle="--",
                    zorder=3,
                    label=(
                        f"1st (σ={sig:.3g})" if show_labels else None
                    ),
                )
                ax.plot(
                    centers, hist_r2, color=sc, lw=1.0, linestyle="-.",
                    zorder=3,
                    label=(
                        f"2nd (σ={sig:.3g})" if show_labels else None
                    ),
                )
                ax.plot(
                    centers, hist_ex, color=sc, lw=1.5, linestyle="-",
                    label=(
                        f"exact (σ={sig:.3g})" if show_labels else None
                    ),
                )

                # Ratio panel: reweight / explicit-shift data.
                ax_r.plot(
                    centers, _safe_ratio(hist_r1, hist_d),
                    color=sc, lw=1.0, linestyle="--", zorder=3,
                )
                ax_r.plot(
                    centers, _safe_ratio(hist_r2, hist_d),
                    color=sc, lw=1.0, linestyle="-.", zorder=3,
                )
                ax_r.plot(
                    centers, _safe_ratio(hist_ex, hist_d),
                    color=sc, lw=1.5, linestyle="-",
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
                ax.legend(fontsize=6, loc="best", ncol=1)
    fig.suptitle(
        "MC shift vs score/Hessian reweighting  "
        f"(σ ∈ {[f'{s:.3g}' for s in sigmas]}; "
        "δ in units of target_std; strict Taylor + exact lp-diff)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "mc_shift_reweight.png"))


def plot_mc_smear_reweight(
    wrapper, target, c_raw, w, preproc, out_dir, sigmas=None,
    smear_factors=(0.1, 0.25, 0.5),
    n_bins=80, n_std=5.0, seed=0,
):
    """Validate the model's ``s² + H_ii`` combination by adding
    Gaussian noise to MC events along each target axis.

    For each axis ``i`` and magnitude
    ``σ_smear = factor · target_std[i]`` we compare:

      * **data** (red shaded): ``{y_m + σ_smear·ε_m·e_i}``,
        ``ε_m ~ N(0,1)``.
      * **unsmeared MC** (gray dotted): ``{y_m}`` reference.
      * For each σ in ``sigmas`` (color-coded):
          - **Taylor reweight** (dashed): weights
            ``w_m · [1 + 0.5·σ_smear²·(s_i(σ)² + H_ii(σ))]``, strict
            linear truncation in σ_smear² using
            ``∂²p_σ/∂y² / p_σ = s² + H_ii``.
          - **exact** (solid): weights
            ``w_m · exp(lp_σ(y_m − σ_smear·ε_m·e_i) − lp_σ(y_m))``
            using the *same* ε as the explicit smear (one-sample MC
            identity for the convolution).

    Each panel has a ratio sub-panel below it plotting the reweight
    curves divided by the explicit-smear data histogram; ratios of 1
    indicate a matching reweight.

    If ``sigmas`` is None, defaults to ``[0.01] + list(smear_factors)``
    — the minimal model σ plus values matching the σ_smear scales
    tested.
    """
    if sigmas is None:
        sigmas = [0.01] + list(smear_factors)
    std = np.asarray(preproc.target_std, dtype=np.float64)
    mean = np.asarray(preproc.target_mean, dtype=np.float64)
    n_mag = len(smear_factors)
    sigmas = list(sigmas)
    rng = np.random.default_rng(seed)
    # One ε draw per event, shared across factors and σ values and
    # reused for explicit smear + exact reweight (same-ε identity).
    eps = rng.standard_normal(target.shape[0])

    target_t = torch.from_numpy(target.astype(np.float32))
    c_t = torch.from_numpy(c_raw.astype(np.float32))

    print(
        f"  smear reweight: precomputing s/H/lp0 at "
        f"{len(sigmas)} σ values × {target.shape[0]} events"
    )
    per_sigma = []
    for sig in sigmas:
        s_sig = _chunked_score(
            wrapper, target_t, c_t, sig, chunk=8192,
        ).numpy()
        H_full = _chunked_hessian(
            wrapper, target_t, c_t, sig, chunk=1024,
        ).numpy()
        H_sig = np.stack([H_full[:, k, k] for k in range(3)], axis=1)
        del H_full
        lp0_sig = _chunked_unnorm_log_density(
            wrapper, target_t, c_t, sig,
        )
        per_sigma.append((sig, s_sig, H_sig, lp0_sig))

    cmap = plt.get_cmap("viridis")
    if len(sigmas) == 1:
        sig_colors = [cmap(0.5)]
    else:
        sig_colors = [
            cmap(i / (len(sigmas) - 1)) for i in range(len(sigmas))
        ]

    fig, axes = plt.subplots(
        2 * n_mag, 3, figsize=(15, 4.0 * n_mag),
        gridspec_kw={"height_ratios": [3, 1] * n_mag},
        sharex="col",
    )
    if n_mag == 1 and axes.ndim == 1:
        axes = axes[None, :]
    for col in range(3):
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

            ax.plot(
                centers, hist_orig, color="gray", lw=0.8,
                linestyle=":", label="unsmeared MC",
            )
            ax.stairs(
                hist_d, edges, color="C3", fill=True, alpha=0.35,
                label="smeared data",
            )

            # Perturbed positions shared across σ for the exact reweight.
            y_perturbed = target.copy()
            y_perturbed[:, col] -= sigma_smear * eps
            y_perturbed_t = torch.from_numpy(
                y_perturbed.astype(np.float32),
            )

            show_labels = (row == 0 and col == 0)
            for (sig, s_sig, H_sig, lp0_sig), sc in zip(
                per_sigma, sig_colors,
            ):
                s_col = s_sig[:, col]
                H_col = H_sig[:, col]

                w_rw_factor = 1.0 + 0.5 * sigma_smear * sigma_smear * (
                    s_col ** 2 + H_col
                )
                hist_r, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w * w_rw_factor, density=True,
                )

                lp_pert = _chunked_unnorm_log_density(
                    wrapper, y_perturbed_t, c_t, sig,
                )
                w_exact = w * np.exp(lp_pert - lp0_sig)
                hist_ex, _ = np.histogram(
                    target[:, col], bins=n_bins, range=(lo, hi),
                    weights=w_exact, density=True,
                )

                # zorder: exact solid keeps default (=2); Taylor
                # dashed overlay uses zorder=3 so it remains
                # visible above the exact curve when close.
                ax.plot(
                    centers, hist_r, color=sc, lw=1.0, linestyle="--",
                    zorder=3,
                    label=(
                        f"Taylor (σ={sig:.3g})" if show_labels else None
                    ),
                )
                ax.plot(
                    centers, hist_ex, color=sc, lw=1.5, linestyle="-",
                    label=(
                        f"exact (σ={sig:.3g})" if show_labels else None
                    ),
                )

                # Ratio panel: reweight / explicit-smear data.
                ax_r.plot(
                    centers, _safe_ratio(hist_r, hist_d),
                    color=sc, lw=1.0, linestyle="--", zorder=3,
                )
                ax_r.plot(
                    centers, _safe_ratio(hist_ex, hist_d),
                    color=sc, lw=1.5, linestyle="-",
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
                ax.legend(fontsize=6, loc="best", ncol=1)
    fig.suptitle(
        "MC smear vs (s² + H_ii) reweighting  "
        f"(σ ∈ {[f'{s:.3g}' for s in sigmas]}; "
        "σ_smear in units of target_std; strict Taylor + exact lp-diff)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, "mc_smear_reweight.png"))


def plot_mc_hessian_stratified(
    target, c_raw, w, s_data, H_diag, preproc, out_dir, sigma, stratify_by,
    n_bins=50, n_std=5.0,
):
    """Stratified data-vs-model comparison of the 1D log-marginal
    second derivative ``d²log p_marg / dy_i²``.

    Per-stratum curves:
      - data: second finite difference of ``log(hist(y_i))`` on a
        weighted histogram of stratum events.
      - model: the per-bin identity
        ``Var[s_i | y_i] + E[H_ii | y_i]``,
        evaluated on the stratum's events. Both quantities are
        aggregated from the precomputed ``s_data`` and ``H_diag``
        arrays — no per-point forward density evaluations.

    Residual panel: ``model − data`` absolute difference (ratios
    blow up at zero crossings of the second derivative).
    """
    groups = _stratify_groups(c_raw, w, preproc, stratify_by)
    fig, axes = plt.subplots(
        2, 3, figsize=(15, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex="col",
    )
    for i in range(3):
        ax, ax_r = axes[0, i], axes[1, i]
        mean_i = preproc.target_mean[i]
        std_i = preproc.target_std[i]
        lo, hi = mean_i - n_std * std_i, mean_i + n_std * std_i

        for g, (label, mask) in enumerate(groups):
            if mask.sum() == 0:
                continue
            color = f"C{g}"
            hist, edges = np.histogram(
                target[mask, i], bins=n_bins, range=(lo, hi),
                weights=w[mask], density=True,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])
            _, d2_data = _fd_log(centers, hist)

            tgt_m = target[mask]
            s_m = s_data[mask]
            H_m = H_diag[mask]
            w_m = w[mask]
            in_range = (tgt_m[:, i] >= lo) & (tgt_m[:, i] < hi)
            bin_idx = np.clip(
                np.digitize(tgt_m[in_range, i], edges) - 1,
                0, n_bins - 1,
            )
            s_i = s_m[in_range, i]
            H_ii = H_m[in_range, i]
            w_in = w_m[in_range]
            mean_s = _weighted_bin_mean(s_i, bin_idx, n_bins, w_in)
            mean_s_sq = _weighted_bin_mean(
                s_i * s_i, bin_idx, n_bins, w_in
            )
            mean_H = _weighted_bin_mean(H_ii, bin_idx, n_bins, w_in)
            d2_model = mean_s_sq - mean_s ** 2 + mean_H

            ax.plot(
                centers, d2_data, "o", color=color, ms=3, alpha=0.7,
                label=f"{label}  data",
            )
            ax.plot(
                centers, d2_model, color=color, lw=1.5,
                label=f"{label}  model",
            )
            ax_r.plot(
                centers, d2_model - d2_data,
                color=color, lw=1.2, alpha=0.85,
            )

        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(f"∂²/∂{TARGET_NAMES[i]}² log p_marg")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=6, loc="best", ncol=2)
        ax_r.axhline(0.0, color="k", lw=0.5)
        ax_r.set_ylabel("model − data")
        ax_r.set_xlabel(TARGET_NAMES[i])
        ax_r.grid(alpha=0.3)

    fig.suptitle(
        f"2nd-derivative (diagonal Hessian) comparison, stratified by "
        f"{stratify_by}, σ = {sigma:.4g}  "
        f"(data = FD² of log hist; model = ⟨Var[s|y] + ⟨H|y⟩⟩)"
    )
    fig.tight_layout()
    _save(
        fig,
        os.path.join(out_dir, f"mc_hessian_stratified_{stratify_by}.png"),
    )


def plot_mc_slices_stratified(
    target, c_raw, w, s_data, preproc, out_dir, sigma, stratify_by,
    n_bins=80, n_std=5.0,
):
    """Stratified MC comparison. Per-stratum data histogram and model
    marginal reconstructed from the stratum's binned scores
    (``d log p_marg / dy_i = E[s_i | y_i]``). Same one-pass
    ``s_data`` array used by ``plot_mc_marginals`` — no extra forward
    passes for any stratum.
    """
    groups = _stratify_groups(c_raw, w, preproc, stratify_by)
    fig, axes = plt.subplots(
        2, 3, figsize=(15, 6.5),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex="col",
    )
    for i in range(3):
        ax, ax_r = axes[0, i], axes[1, i]
        for g, (label, mask) in enumerate(groups):
            if mask.sum() == 0:
                continue
            color = f"C{g}"
            hist, edges, centers = _data_1d_histogram(
                target[mask], w[mask], i, preproc, n_bins, n_std,
            )
            lo, hi = edges[0], edges[-1]
            tgt_m = target[mask]
            s_m = s_data[mask]
            w_m = w[mask]
            in_range = (tgt_m[:, i] >= lo) & (tgt_m[:, i] < hi)
            bin_idx = np.clip(
                np.digitize(tgt_m[in_range, i], edges) - 1,
                0, n_bins - 1,
            )
            mean_s = _weighted_bin_mean(
                s_m[in_range, i], bin_idx, n_bins, w_m[in_range],
            )
            p_model = _reconstruct_marginal_from_score(centers, mean_s)

            ax.stairs(
                hist, edges, color=color, alpha=0.32, fill=True,
            )
            ax.plot(
                centers, p_model, color=color, lw=1.5,
                label=f"{label}  (N={int(mask.sum())})",
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(p_model > 0, hist / p_model, np.nan)
            ax_r.step(
                edges[:-1], ratio, where="post",
                color=color, alpha=0.85,
            )
        ax.set_ylabel("density")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=7, loc="best")
        ax_r.axhline(1.0, color="k", lw=0.5)
        ax_r.set_ylabel("data / model")
        ax_r.set_xlabel(TARGET_NAMES[i])
        ax_r.set_ylim(0, 2)
        ax_r.grid(alpha=0.3)
    fig.suptitle(
        f"MC vs model 1D marginals, stratified by "
        f"{stratify_by}, σ = {sigma:.4g}  "
        f"(shaded = data; line = model score-integrated over stratum)"
    )
    fig.tight_layout()
    _save(fig, os.path.join(out_dir, f"mc_slices_{stratify_by}.png"))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    if args.output is None:
        args.output = os.path.dirname(os.path.abspath(args.checkpoint))
    os.makedirs(args.output, exist_ok=True)

    print(f"loading wrapper from {args.checkpoint}")
    wrapper, preproc, model_config = load_wrapper(
        args.checkpoint, args.device
    )

    sigma = (
        float(args.sigma)
        if args.sigma is not None
        else float(wrapper.sigma_inference.item())
    )
    print(f"single-σ plots at σ = {sigma:.4g}")

    smin = float(model_config.get("sigma_min", 0.01))
    smax = float(model_config.get("sigma_max", 0.3))
    if args.sigma_grid:
        sigma_grid = [float(s) for s in args.sigma_grid]
    else:
        sigma_grid = np.geomspace(smin, smax, args.n_sigma).tolist()
    print(
        f"sigma sweep: {len(sigma_grid)} points in "
        f"[{sigma_grid[0]:.4g}, {sigma_grid[-1]:.4g}]"
    )

    if args.reweight_sigmas:
        reweight_sigmas = [float(s) for s in args.reweight_sigmas]
        print(
            f"reweight σ values (override): "
            f"{[f'{s:.4g}' for s in reweight_sigmas]}"
        )
    else:
        # Each plot picks its own default: [0.01] + its factor list.
        reweight_sigmas = None

    # Training curve (cheap, do first so it's available even if the
    # score plots error).
    log_path = args.training_log or os.path.join(
        os.path.dirname(os.path.abspath(args.checkpoint)),
        "training.log",
    )
    if os.path.exists(log_path):
        plot_training_curve(log_path, args.output)
    else:
        print(f"[skip] training log not found at {log_path}")

    contexts = build_contexts(preproc)
    print(f"evaluating at {len(contexts)} conditioning points")

    plot_log_p_profiles(
        wrapper, contexts, preproc, args.output,
        sigma, args.n_points, args.n_range_std,
    )
    plot_log_p_2d(
        wrapper, contexts, preproc, args.output,
        sigma, args.n_grid_2d,
    )
    plot_score_profiles(
        wrapper, contexts, preproc, args.output,
        sigma, args.n_points, args.n_range_std,
    )
    plot_hessian_profiles(
        wrapper, contexts, preproc, args.output,
        sigma, args.n_points, args.n_range_std,
    )
    plot_sigma_sweep(
        wrapper, contexts, preproc, args.output, sigma_grid,
    )
    plot_eigenvalues(
        wrapper, contexts, preproc, args.output, sigma,
    )

    if args.snapshot:
        data = load_snapshot_data(
            args.snapshot, preproc, args.max_data_samples,
            tree=args.tree, threads=args.threads,
            pt_min=args.pt_min, pt_max=args.pt_max, eta_max=args.eta_max,
        )
        if data is not None:
            target, c_raw_np, w_np = data
            # Precompute the model score and diagonal Hessian at every
            # data event once. All MC plots below read these arrays —
            # no plot does any further forward-density evaluation,
            # which is what used to dominate the cost.
            print(
                f"  evaluating model score + diag-Hessian "
                f"at {target.shape[0]} data events"
            )
            target_t = torch.from_numpy(target.astype(np.float32))
            c_t = torch.from_numpy(c_raw_np.astype(np.float32))
            s_data = _chunked_score(
                wrapper, target_t, c_t, sigma, chunk=8192,
            ).numpy()
            H_full = _chunked_hessian(
                wrapper, target_t, c_t, sigma, chunk=1024,
            ).numpy()
            H_diag = np.stack(
                [H_full[:, k, k] for k in range(3)], axis=1,
            )
            del H_full

            plot_score_at_data(
                wrapper, target, c_raw_np, w_np, preproc,
                args.output, sigma,
            )
            plot_mc_marginals(
                target, c_raw_np, w_np, s_data, preproc,
                args.output, sigma,
            )
            plot_mc_derivatives(
                target, c_raw_np, w_np, s_data, H_diag, preproc,
                args.output, sigma,
            )
            plot_mc_shift_reweight(
                wrapper, target, c_raw_np, w_np,
                preproc, args.output, reweight_sigmas,
            )
            plot_mc_smear_reweight(
                wrapper, target, c_raw_np, w_np,
                preproc, args.output, reweight_sigmas,
            )
            for by in ("pt", "lambda", "charge"):
                plot_mc_slices_stratified(
                    target, c_raw_np, w_np, s_data, preproc,
                    args.output, sigma, by,
                )
                plot_mc_hessian_stratified(
                    target, c_raw_np, w_np, s_data, H_diag, preproc,
                    args.output, sigma, by,
                )
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
