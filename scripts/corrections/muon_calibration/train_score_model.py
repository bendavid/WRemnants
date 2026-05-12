"""Train an energy base + structured reweight head for the muon response.

Two sub-networks, trained jointly in a single optimizer step:

  1. **Energy base** ``E_θ(y, c, σ)``: scalar MLP. Score via autograd,
     ``s_θ = -∇_y E_θ``. Curl-free by construction. Trained with 1st-
     order Denoising Score Matching (DSM, Vincent 2011):
         L_DSM = E[ ‖σ·s_θ(y+σε, c, σ) + ε‖² ]
     with optional σ²-weighting for the Song et al. 2021 maximum-
     likelihood variational bound on ``log p_0`` (``--mle-weight``,
     default off). Second-order DSM (Hessian matching) has been
     removed — the Hessian is no longer a training product.

  2. **Reweight head** ``H_φ(y, c, σ, s, δ, v) → (R_shift, S_smear)``:
     structured form with context-only quadratic tensors plus free-
     form residuals:
         R_shift = ½·δᵀ·M_φ(context)·δ  +  R_res_φ(context, δ, v)
         S_smear = ½·(sᵀv)²  +  ½·vᵀ·N_φ(context)·v
                   + S_res_φ(context, δ, v)
     where ``context = (y, c, log σ, s)`` with ``s`` the base's
     detached score. ``M`` and ``N`` are 3×3 symmetric tensors (6
     coefficients each) predicted from context alone; they recover
     the Hessian-level Taylor behavior analytically:
         R_shift ≈ ½·δᵀH·δ     (leading O(|δ|²))
         S_smear ≈ ½·vᵀ(ssᵀ+H)·v    (leading O(|v|²),
                                     deterministic smear)
     The residuals ``R_res, S_res`` absorb all higher-order terms.
     The full shift reweight in lp-space is ``Δlp_shift = -δ·s +
     R_shift`` — linear-in-δ dominant (from the analytic score
     prefix), quadratic from the structured M term, cubic+ from
     the residual.

  **Smear target — deterministic.** ``S_smear`` is trained against
  ``logsumexp_k[lp(y-ε_k·v) - lp(y)] - log K`` with K MC samples of
  ``ε ~ N(0,1)``; the trained output is the log of the expectation-
  over-ε weight, i.e. deterministic in (y, c, σ, v).

  **Head inputs.** The head sees both δ and v simultaneously but the
  supervision is decoupled (shift target doesn't use v, smear
  target doesn't use δ), so at convergence R_shift and S_smear are
  almost entirely functions of only their relevant argument. For a
  combined shift+smear reweight, the cheap default is the naive
  additive sum; the exact-factorization path evaluates S_smear at
  the shifted position ``y - δ``.

  **Joint training with clean separation.** Base-model outputs used
  as teacher targets for the head (``lp``, ``s``) are detached, so
  gradients from head losses do not flow back into the base. The
  base is trained only by DSM; the head is trained only by its
  distillation losses. Both in a single backward pass per batch.

Noise-level conditioning. The base takes σ as an extra input:
``E_θ(y, c, σ)``, trained on a log-uniform schedule
``σ ~ LogUniform(σ_min, σ_max)`` (defaults 0.01..0.3 in
standardized-target units). Queryable at any σ at inference.

Inference (ScoreWrapper methods):
    unnormalized_log_density(y_raw, c_raw[, σ])
        → -E_θ(y, c, σ) up to a ``c, σ``-constant.
    score(y_raw, c_raw[, σ])
        → ∇_y log p_σ. One forward + one backward.
    apply_shift(y_raw, c_raw, δ_raw[, σ])
        → Δlp for shift δ. One base forward + score + one head call.
    apply_smear(y_raw, c_raw, v_raw[, σ])
        → Δlp for rank-1 smearing v. One base forward + score + one
        head call. Deterministic.
    apply_combined(y_raw, c_raw, δ_raw, v_raw[, σ], exact=False)
        → Δlp for combined shift + smear. Default naive additive;
        ``exact=True`` uses the factorization
        Δlp_combined = Δlp_shift(y, δ) + Δlp_smear(y-δ, v).

Input: the snapshot produced by ``flow_training_snapshot.py``.

Outputs (to ``--output``):
    score_model.pt    ScoreWrapper state_dict (includes base + head).
    score_scripted.pt TorchScript-traced base energy (head not traced).
    preproc.json      Preprocessing stats + model config.
    training.log      Per-epoch loss components (DSM, shift, smear).
    checkpoint.pt     Per-epoch snapshot (overwritten) for resume.

DDP: ``--num-gpus N`` spawns N worker processes via
torch.multiprocessing.spawn. Shared CPU-tensor preprocessing, per-rank
index sharding, and metric all-reduce.

Requires: torch, numpy, and ROOT for snapshot loading.
"""

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
from tqdm import tqdm


ACTIVATIONS = {
    "gelu": nn.GELU,
    "silu": nn.SiLU,
    "mish": nn.Mish,
    "softplus": nn.Softplus,
    "elu": nn.ELU,
    "tanh": nn.Tanh,
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="Explicit list of ROOT snapshot files (from "
        "flow_training_snapshot.py). Paths are passed directly to "
        "RDataFrame — no recursive directory search.",
    )
    p.add_argument(
        "--tree", default="tree", help="TTree name inside the input files"
    )
    p.add_argument(
        "--output",
        default="./muon_response_score/",
        help="Output directory",
    )
    p.add_argument(
        "--max-muons",
        type=int,
        default=-1,
        help="Cap on number of muon rows loaded (mu+ and mu- each count "
        "as one row). -1 loads all.",
    )
    p.add_argument(
        "--val-fraction",
        type=float,
        default=0.1,
        help="Fraction of muons held out for validation.",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Training batch size. Defaults: 32768 on CUDA, 16384 on CPU. "
        "DSM requires a double-backward per step, so reduce if OOM.",
    )
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight-decay", type=float, default=0.0)
    p.add_argument(
        "--patience",
        type=int,
        default=5,
        help="Early-stopping patience on validation DSM loss.",
    )
    p.add_argument(
        "--hidden-features",
        type=int,
        default=256,
        help="Hidden-layer width of the energy MLP.",
    )
    p.add_argument(
        "--n-hidden-layers",
        type=int,
        default=4,
        help="Number of hidden layers in the energy MLP.",
    )
    p.add_argument(
        "--activation",
        default="gelu",
        choices=sorted(ACTIVATIONS),
        help="Energy-MLP activation. Default 'gelu' is C^∞, so the "
        "score field (first derivative of E) is C^∞ and the score "
        "derivatives (second derivative of E) are at least C^0. "
        "Avoid 'elu' here — its second derivative has a jump at 0, "
        "which shows up as a discontinuity in ∂s/∂y. 'softplus' "
        "and 'silu'/'mish' are also smooth choices.",
    )
    p.add_argument(
        "--sigma-schedule",
        default="log-uniform",
        choices=["log-uniform", "fixed"],
        help="'log-uniform' (default) trains with σ-conditioning: "
        "the network takes log(σ) as an input and σ is sampled per-"
        "sample from LogUniform(--sigma-min, --sigma-max). 'fixed' "
        "trains DSM at --sigma-dsm for every sample.",
    )
    p.add_argument(
        "--sigma-min",
        type=float,
        default=0.01,
        help="Lower bound of the log-uniform σ schedule, in "
        "standardized-target units. Also the default --sigma-"
        "inference.",
    )
    p.add_argument(
        "--sigma-max",
        type=float,
        default=0.3,
        help="Upper bound of the log-uniform σ schedule.",
    )
    p.add_argument(
        "--sigma-inference",
        type=float,
        default=None,
        help="Default σ baked into ScoreWrapper for inference calls "
        "without an explicit σ. Unset → uses --sigma-min.",
    )
    p.add_argument(
        "--sigma-dsm",
        type=float,
        default=0.02,
        help="(--sigma-schedule fixed only) Noise level used for "
        "DSM at every sample.",
    )
    p.add_argument(
        "--mle-weight",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use σ²-weighted DSM (Song et al. 2021 MLE weighting), "
        "which yields a variational lower bound on log p_0. Default "
        "off: per-sample DSM loss has uniform σ weighting, giving "
        "balanced per-σ score accuracy across the schedule.",
    )
    p.add_argument(
        "--head-hidden-features",
        type=int,
        default=128,
        help="Hidden-layer width of the reweight head MLP. Head is "
        "typically much smaller than the base energy MLP since it "
        "only has to fit the Taylor residual.",
    )
    p.add_argument(
        "--head-n-hidden-layers",
        type=int,
        default=3,
        help="Number of hidden layers in the reweight head MLP.",
    )
    p.add_argument(
        "--shift-scale",
        type=float,
        default=0.3,
        help="Half-width of the uniform distribution used to sample "
        "training shift vectors δ per-component (in standardized-"
        "target units). Values outside the range ±scale are "
        "extrapolation at inference.",
    )
    p.add_argument(
        "--smear-scale",
        type=float,
        default=0.3,
        help="Half-width of the uniform distribution used to sample "
        "training smearing vectors v per-component (in standardized-"
        "target units).",
    )
    p.add_argument(
        "--shift-loss-weight",
        type=float,
        default=1.0,
        help="Relative weight on the shift-head loss in the joint "
        "objective. 0 disables head-shift training entirely.",
    )
    p.add_argument(
        "--smear-loss-weight",
        type=float,
        default=1.0,
        help="Relative weight on the smear-head loss in the joint "
        "objective. 0 disables head-smear training entirely.",
    )
    p.add_argument(
        "--smear-mc-samples",
        type=int,
        default=4,
        help="Number of Gaussian ε samples used to estimate the "
        "deterministic smear target ``log E_ε[exp(Δlp_ε)]`` per "
        "training event. Higher K reduces Jensen bias but costs "
        "K base forward passes per step.",
    )
    p.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="cuda or cpu. Multi-GPU requires cuda (see --num-gpus).",
    )
    p.add_argument(
        "--amp",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable bfloat16 autocast for the training/validation "
        "forward passes. Default: off (FP32). bf16 is risky for DSM "
        "because the create_graph=True double backward amplifies "
        "rounding error; turn it on only if FP32 training is the "
        "bottleneck and accept a higher noise floor on the score.",
    )
    p.add_argument(
        "--num-gpus",
        type=int,
        default=-1,
        help="Number of GPUs to use when --device=cuda. -1 "
        "(default) auto-detects via torch.cuda.device_count().",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="RDataFrame ImplicitMT threads during the snapshot load. "
        "0 = ROOT's default; 1 disables MT.",
    )
    p.add_argument(
        "--pt-min",
        type=float,
        default=2.0,
        help="Minimum gen pt (GeV) for each muon.",
    )
    p.add_argument(
        "--pt-max",
        type=float,
        default=200.0,
        help="Maximum gen pt (GeV).",
    )
    p.add_argument(
        "--eta-max",
        type=float,
        default=2.4,
        help="Absolute value of gen/reco eta cut.",
    )
    return p.parse_args()


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

BRANCHES_PER_MUON = {
    "plus": {
        "pt_reco": "Mupluscor_pt",
        "eta_reco": "Mupluscor_eta",
        "phi_reco": "Mupluscor_phi",
        "pt_gen": "Muplusgen_pt",
        "eta_gen": "Muplusgen_eta",
        "phi_gen": "Muplusgen_phi",
    },
    "minus": {
        "pt_reco": "Muminuscor_pt",
        "eta_reco": "Muminuscor_eta",
        "phi_reco": "Muminuscor_phi",
        "pt_gen": "Muminusgen_pt",
        "eta_gen": "Muminusgen_eta",
        "phi_gen": "Muminusgen_phi",
    },
}

WEIGHT_BRANCH = "nominal_weight"


def load_ntuples(
    paths: List[str],
    tree_name: str,
    max_muons: int,
    pt_min: float,
    pt_max: float,
    eta_max: float,
    threads: int = 0,
):
    """Return per-muon arrays. Same schema as train_muon_response_flow.py."""
    import ROOT

    if not ROOT.ROOT.IsImplicitMTEnabled():
        if threads == 0:
            ROOT.ROOT.EnableImplicitMT()
        elif threads > 1:
            ROOT.ROOT.EnableImplicitMT(threads)

    files_vec = ROOT.std.vector("string")()
    for p in paths:
        files_vec.push_back(p)
    df = ROOT.ROOT.RDataFrame(tree_name, files_vec)
    ROOT.ROOT.RDF.Experimental.AddProgressBar(df)

    filt = (
        f"Muplusgen_pt > {pt_min} && Muminusgen_pt > {pt_min} && "
        f"Muplusgen_pt < {pt_max} && Muminusgen_pt < {pt_max} && "
        f"std::fabs(Muplusgen_eta) < {eta_max} && "
        f"std::fabs(Muminusgen_eta) < {eta_max} && "
        f"Mupluscor_pt > 0. && Muminuscor_pt > 0. && "
        f"{WEIGHT_BRANCH} > 0."
    )
    df = df.Filter(filt, "quality")

    # If ``max_muons`` is set well below the dataset size, drop most
    # rows inside the RDF graph via an ``rdfentry_ % stride == 0``
    # filter so the subsequent AsNumpy pass only materializes ~the
    # rows we will keep. This cuts the event-loop cost, the post-
    # event-loop numpy copies (astype + concat of mu+/mu- into float64),
    # and the per-event trig in compute_targets_and_conditioning all
    # by the same factor. A short Count() first tells us the stride;
    # the subsequent AsNumpy reads are usually page-cache-warm.
    if max_muons > 0:
        print("counting filtered events to size the stride subsample...")
        n_pass = int(df.Count().GetValue())
        # Each event yields two muons (mu+ and mu-) after pooling;
        # pick enough events to reach max_muons with a small buffer.
        n_needed = max(1, max_muons // 2)
        if n_pass > 2 * n_needed:
            stride = max(2, n_pass // n_needed)
            df = df.Filter(
                f"rdfentry_ % {stride} == 0", "stride_subsample",
            )
            print(
                f"  quality-passing events: {n_pass}; "
                f"stride-subsampling 1/{stride} → "
                f"~{n_pass // stride} events"
            )
        else:
            print(
                f"  quality-passing events: {n_pass}; no stride "
                f"subsample needed (max_muons={max_muons})"
            )

    all_branches = {WEIGHT_BRANCH}
    for m in BRANCHES_PER_MUON.values():
        all_branches.update(m.values())

    arrs = df.AsNumpy(columns=list(all_branches))
    w_event = arrs[WEIGHT_BRANCH].astype(np.float64)

    per_muon_rows = []
    for sign, charge in (("plus", +1.0), ("minus", -1.0)):
        b = BRANCHES_PER_MUON[sign]
        n = w_event.shape[0]
        charge_arr = np.full(n, charge, dtype=np.float64)
        per_muon_rows.append((
            arrs[b["pt_reco"]].astype(np.float64),
            arrs[b["eta_reco"]].astype(np.float64),
            arrs[b["phi_reco"]].astype(np.float64),
            arrs[b["pt_gen"]].astype(np.float64),
            arrs[b["eta_gen"]].astype(np.float64),
            arrs[b["phi_gen"]].astype(np.float64),
            charge_arr,
            w_event,
        ))

    pt_r = np.concatenate([p[0] for p in per_muon_rows])
    eta_r = np.concatenate([p[1] for p in per_muon_rows])
    phi_r = np.concatenate([p[2] for p in per_muon_rows])
    pt_g = np.concatenate([p[3] for p in per_muon_rows])
    eta_g = np.concatenate([p[4] for p in per_muon_rows])
    phi_g = np.concatenate([p[5] for p in per_muon_rows])
    q = np.concatenate([p[6] for p in per_muon_rows])
    w = np.concatenate([p[7] for p in per_muon_rows])

    arrs = (pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q, w)
    n = arrs[0].shape[0]
    print(
        f"loaded {n} muons after filters "
        f"({pt_min} < pt_gen < {pt_max}, |eta| < {eta_max}, w > 0)"
    )
    print(
        f"  weight: mean {arrs[7].mean():.4f}  std {arrs[7].std():.4f}  "
        f"min {arrs[7].min():.4f}  max {arrs[7].max():.4f}"
    )

    if max_muons > 0 and n > max_muons:
        rng = np.random.default_rng(0)
        idx = rng.choice(n, size=max_muons, replace=False)
        arrs = tuple(a[idx] for a in arrs)
        print(f"subsampled to {max_muons} muons for training")

    return arrs


def compute_targets_and_conditioning(
    pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q
):
    """Return (target [N,3], cond_raw dict). Identical to flow script
    so the two approaches use the same pre-network representation."""
    lam_r = np.arctan(np.sinh(eta_r))
    lam_g = np.arctan(np.sinh(eta_g))

    kappa_r = q * np.cos(lam_r) / pt_r
    kappa_g = q * np.cos(lam_g) / pt_g

    r_kappa = kappa_r / kappa_g - 1.0

    dphi = np.arctan2(
        np.sin(phi_r - phi_g), np.cos(phi_r - phi_g)
    )
    dlambda = lam_r - lam_g

    target = np.stack([r_kappa, dlambda, dphi], axis=1).astype(np.float32)

    cond_raw = {
        "log_pt_gen": np.log(pt_g).astype(np.float32),
        "charge": q.astype(np.float32),
        "lambda_gen": lam_g.astype(np.float32),
        "sin_phi_gen": np.sin(phi_g).astype(np.float32),
        "cos_phi_gen": np.cos(phi_g).astype(np.float32),
    }

    return target, cond_raw


# -----------------------------------------------------------------------------
# Preprocessing
# -----------------------------------------------------------------------------

@dataclass
class PreprocStats:
    target_names: List[str]
    target_mean: List[float]
    target_std: List[float]

    cond_names: List[str]
    cond_mean: List[float]
    cond_std: List[float]


def build_preproc(target: np.ndarray, cond_raw: dict) -> PreprocStats:
    target_names = ["r_kappa", "dlambda", "dphi"]
    target_mean = target.mean(axis=0).tolist()
    target_std = target.std(axis=0).tolist()

    cond_names = [
        "log_pt_gen",
        "charge",
        "lambda_gen",
        "sin_phi_gen",
        "cos_phi_gen",
    ]
    cond_mean, cond_std = [], []
    for name in cond_names:
        arr = cond_raw[name]
        cond_mean.append(float(arr.mean()))
        cond_std.append(float(arr.std()) if arr.std() > 1e-6 else 1.0)

    return PreprocStats(
        target_names=target_names,
        target_mean=target_mean,
        target_std=target_std,
        cond_names=cond_names,
        cond_mean=cond_mean,
        cond_std=cond_std,
    )


def apply_preproc(
    target: np.ndarray, cond_raw: dict, stats: PreprocStats
) -> Tuple[np.ndarray, np.ndarray]:
    tmean = np.asarray(stats.target_mean, dtype=np.float32)
    tstd = np.asarray(stats.target_std, dtype=np.float32)
    target_std = (target - tmean) / tstd

    cond_cols = []
    for name, mean, std in zip(stats.cond_names, stats.cond_mean, stats.cond_std):
        cond_cols.append((cond_raw[name] - mean) / std)
    cond = np.stack(cond_cols, axis=1).astype(np.float32)

    return target_std, cond


# -----------------------------------------------------------------------------
# Energy-based model
# -----------------------------------------------------------------------------

class EnergyMLP(nn.Module):
    """Scalar energy ``E_θ(y, c, σ)``; score is ``-∇_y E``.

    Plain concat-input MLP: input ``[y, c, log σ]`` passes through
    ``n_hidden_layers`` Linear+activation pairs ending in a
    Linear-to-1 readout. σ enters as its log so the geometric spread
    of a log-uniform schedule maps to a linear span. The extra input
    channel costs <1% runtime (only the first linear layer grows by
    one column) in exchange for learning a whole family of σ-smoothed
    densities instead of one — see module docstring.

    Autograd derivatives are taken w.r.t. ``y`` only; ``c`` and
    ``log σ`` are treated as constants at each forward, which is
    exactly what we want: the Hessian we want is ``∂²log p(y|c, σ) /
    ∂y²``, not ``∂²/∂σ²`` or mixed partials.
    """

    def __init__(
        self,
        n_features: int,
        n_cond: int,
        hidden_features: int,
        n_hidden_layers: int,
        activation: type = nn.GELU,
    ):
        super().__init__()
        self.n_features = n_features
        self.n_cond = n_cond
        layers: List[nn.Module] = []
        # +1 input column for log σ (noise-level conditioning).
        in_dim = n_features + n_cond + 1
        for _ in range(n_hidden_layers):
            layers.append(nn.Linear(in_dim, hidden_features))
            layers.append(activation())
            in_dim = hidden_features
        layers.append(nn.Linear(in_dim, 1))
        self.net = nn.Sequential(*layers)

    @staticmethod
    def _broadcast_log_sigma(
        log_sigma: torch.Tensor, y: torch.Tensor
    ) -> torch.Tensor:
        """Turn a per-sample or scalar log-σ into a trailing-dim
        tensor compatible with concatenation against ``y`` on the
        last axis."""
        if log_sigma.dim() == 0:
            log_sigma = log_sigma.expand(y.shape[:-1])
        if log_sigma.dim() == y.dim() - 1:
            log_sigma = log_sigma.unsqueeze(-1)
        return log_sigma

    def forward(
        self,
        y: torch.Tensor,
        c: torch.Tensor,
        log_sigma: torch.Tensor,
    ) -> torch.Tensor:
        log_sigma = self._broadcast_log_sigma(log_sigma, y)
        h = torch.cat([y, c, log_sigma], dim=-1)
        return self.net(h).squeeze(-1)

    def score(
        self,
        y: torch.Tensor,
        c: torch.Tensor,
        log_sigma: torch.Tensor,
        create_graph: bool = False,
    ) -> torch.Tensor:
        """Return the conservative score field ``s = -∇_y E``.

        ``create_graph=True`` during training so the loss backward
        can propagate through the gradient computation to the
        network parameters (classic DSM double-backward).
        """
        y = y.requires_grad_(True)
        E = self.forward(y, c, log_sigma)
        (grad_y,) = torch.autograd.grad(
            E.sum(), y, create_graph=create_graph
        )
        return -grad_y

class HeadMLP(nn.Module):
    """Structured reweight head producing ``R_shift`` and ``S_smear``.

    Output structure:

      ``R_shift = ½·δᵀ·M(context)·δ + R_res(context, δ, v)``
      ``S_smear = ½·(sᵀv)² + ½·vᵀ·N(context)·v + S_res(context, δ, v)``

    where ``context = (y, c, log σ, s)`` (s detached), and ``M`` and
    ``N`` are symmetric ``d × d`` tensors predicted from context
    alone — so the quadratic forms are manifestly quadratic in δ
    and v respectively. The ``(sᵀv)²`` term in S_smear uses the
    input score directly (no learning required). ``R_res`` and
    ``S_res`` are free-form scalar residuals absorbing higher-order
    corrections.

    Two MLP branches share no weights (cheap; tiny vs the base):

      *context branch*: input ``[y, c, log σ, s]`` → 2·d(d+1)/2
        output coefficients for the two symmetric tensors.
      *residual branch*: input ``[y, c, log σ, s, δ, v]`` → 2 scalars.

    All training supervision uses ``s`` and teacher lp values
    detached from the base, so head gradients do not perturb the
    base training.
    """

    def __init__(
        self,
        n_features: int,
        n_cond: int,
        hidden_features: int,
        n_hidden_layers: int,
        activation: type = nn.GELU,
    ):
        super().__init__()
        d = int(n_features)
        self.n_features = d
        self.n_cond = int(n_cond)
        self.n_sym = d * (d + 1) // 2

        # +1 for log σ, +d for score.
        ctx_dim = d + n_cond + 1 + d
        res_dim = ctx_dim + d + d  # + δ + v

        def mlp(in_dim: int, out_dim: int) -> nn.Sequential:
            layers: List[nn.Module] = []
            cur = in_dim
            for _ in range(n_hidden_layers):
                layers.append(nn.Linear(cur, hidden_features))
                layers.append(activation())
                cur = hidden_features
            layers.append(nn.Linear(cur, out_dim))
            return nn.Sequential(*layers)

        # Two symmetric tensors: M (shift Hessian) and N (smear
        # anti-score quadratic). Each has d(d+1)/2 free params.
        self.ctx_net = mlp(ctx_dim, 2 * self.n_sym)
        # Two scalar residuals (R_res, S_res).
        self.res_net = mlp(res_dim, 2)

        # Pre-computed (row, col) indices for the upper triangle of
        # a symmetric d × d matrix, used to scatter the predicted
        # coefficients into a full tensor.
        idx_i, idx_j = [], []
        for i in range(d):
            for j in range(i, d):
                idx_i.append(i)
                idx_j.append(j)
        self.register_buffer(
            "_tri_i", torch.tensor(idx_i, dtype=torch.long), persistent=False
        )
        self.register_buffer(
            "_tri_j", torch.tensor(idx_j, dtype=torch.long), persistent=False
        )

    @staticmethod
    def _broadcast_log_sigma(
        log_sigma: torch.Tensor, y: torch.Tensor
    ) -> torch.Tensor:
        if log_sigma.dim() == 0:
            log_sigma = log_sigma.expand(y.shape[:-1])
        if log_sigma.dim() == y.dim() - 1:
            log_sigma = log_sigma.unsqueeze(-1)
        return log_sigma

    def _pack_symmetric(self, coeffs: torch.Tensor) -> torch.Tensor:
        """``[..., d*(d+1)/2]`` → symmetric ``[..., d, d]`` matrix."""
        d = self.n_features
        batch_shape = coeffs.shape[:-1]
        M = coeffs.new_zeros(*batch_shape, d, d)
        M[..., self._tri_i, self._tri_j] = coeffs
        M[..., self._tri_j, self._tri_i] = coeffs  # symmetrize
        return M

    def forward(
        self,
        y: torch.Tensor,
        c: torch.Tensor,
        log_sigma: torch.Tensor,
        score: torch.Tensor,
        delta: torch.Tensor,
        v: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Return ``(R_shift, S_smear)`` scalars per sample.

        All tensor inputs must broadcast on the same leading batch
        shape; ``y, score, δ, v`` carry the trailing ``d`` axis,
        ``c`` the trailing ``n_cond`` axis, and ``log_sigma`` is
        either scalar, per-sample, or per-sample-with-trailing-1.
        """
        log_sigma = self._broadcast_log_sigma(log_sigma, y)
        ctx = torch.cat([y, c, log_sigma, score], dim=-1)

        tensor_coeffs = self.ctx_net(ctx)
        n_sym = self.n_sym
        M = self._pack_symmetric(tensor_coeffs[..., :n_sym])
        N = self._pack_symmetric(tensor_coeffs[..., n_sym:])

        res_in = torch.cat([ctx, delta, v], dim=-1)
        res_out = self.res_net(res_in)
        R_res = res_out[..., 0]
        S_res = res_out[..., 1]

        # ½ δᵀ M δ, handled as a sum of outer products to keep the
        # graph simple and TorchScript-friendly.
        quad_shift = 0.5 * (delta.unsqueeze(-2) @ M @ delta.unsqueeze(-1))
        quad_shift = quad_shift.squeeze(-1).squeeze(-1)

        # ½ (sᵀv)² + ½ vᵀ N v.
        sv = (score * v).sum(dim=-1)
        quad_smear_n = 0.5 * (v.unsqueeze(-2) @ N @ v.unsqueeze(-1))
        quad_smear_n = quad_smear_n.squeeze(-1).squeeze(-1)
        quad_smear = 0.5 * sv.pow(2) + quad_smear_n

        R_shift = quad_shift + R_res
        S_smear = quad_smear + S_res
        return R_shift, S_smear


class JointScoreModel(nn.Module):
    """DDP-friendly wrapper holding the base energy + reweight head.

    forward() runs one training step and returns a scalar combined
    loss (with gradient) plus detached per-component sums for logging.
    The base is trained only by DSM; the head is trained only by its
    distillation losses (teacher lp / score detached from the base).

    Loss components:
      - ``L_DSM = ‖σ·s_θ(y+σε) + ε‖²`` per sample, optionally
        multiplied by ``σ²`` for MLE-DSM weighting.
      - ``L_shift = (R_shift − target_shift)²`` with target
        ``lp(y−δ) − lp(y) + δ·s`` (the 1st-order-Taylor residual).
      - ``L_smear = (S_smear − target_smear)²`` with target
        ``logsumexp_k[lp(y−ε_k·v) − lp(y)] − log K`` (deterministic
        rank-1 smear, K MC samples per event).
    """

    def __init__(
        self,
        energy: EnergyMLP,
        head: HeadMLP,
        sigma_schedule: str,
        sigma_min: float,
        sigma_max: float,
        sigma_dsm: float,
        mle_weight: bool = False,
        shift_scale: float = 0.3,
        smear_scale: float = 0.3,
        shift_loss_weight: float = 1.0,
        smear_loss_weight: float = 1.0,
        smear_mc_samples: int = 4,
    ):
        super().__init__()
        if sigma_schedule not in ("log-uniform", "fixed"):
            raise ValueError(
                f"unknown sigma_schedule '{sigma_schedule}'"
            )
        self.energy = energy
        self.head = head
        self.sigma_schedule = sigma_schedule
        self.mle_weight = bool(mle_weight)
        self.shift_scale = float(shift_scale)
        self.smear_scale = float(smear_scale)
        self.shift_loss_weight = float(shift_loss_weight)
        self.smear_loss_weight = float(smear_loss_weight)
        self.smear_mc_samples = int(smear_mc_samples)
        self.register_buffer(
            "sigma_min", torch.tensor(float(sigma_min), dtype=torch.float32)
        )
        self.register_buffer(
            "sigma_max", torch.tensor(float(sigma_max), dtype=torch.float32)
        )
        self.register_buffer(
            "sigma_dsm", torch.tensor(float(sigma_dsm), dtype=torch.float32)
        )

    def _sample_sigma(
        self, batch_shape: Tuple[int, ...]
    ) -> torch.Tensor:
        """Per-sample σ from the configured schedule."""
        device = self.sigma_min.device
        dtype = self.sigma_min.dtype
        if self.sigma_schedule == "log-uniform":
            log_smin = torch.log(self.sigma_min)
            log_smax = torch.log(self.sigma_max)
            u = torch.rand(batch_shape, device=device, dtype=dtype)
            return torch.exp(log_smin + u * (log_smax - log_smin))
        return self.sigma_dsm.expand(batch_shape)

    def dsm_loss(
        self, y: torch.Tensor, c: torch.Tensor, create_graph: bool,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """1st-order DSM residual ``‖σ·s_θ(y+σε) + ε‖²`` per sample.

        Returns ``(per_sample_loss, sigma)``. The caller applies the
        ``σ²`` MLE weighting if requested — this keeps the raw
        unweighted residual available for logging.
        """
        B = y.shape[:-1]
        sigma = self._sample_sigma(B).to(y.dtype)
        log_sigma = torch.log(sigma)
        eps = torch.randn_like(y)
        y_tilde = y + sigma.unsqueeze(-1) * eps
        score = self.energy.score(
            y_tilde, c, log_sigma, create_graph=create_graph
        )
        r = sigma.unsqueeze(-1) * score + eps
        return r.pow(2).sum(dim=-1), sigma

    def head_losses(
        self, y: torch.Tensor, c: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Shift and smear losses for the head, one sample of (δ, v, {ε}) each.

        All base-model evaluations used as teacher targets are done
        under ``torch.no_grad()`` so gradients from the head loss do
        not reach the base. The head's score input comes from a
        fresh no-grad score evaluation — detached.
        """
        B = y.shape[:-1]
        d = y.shape[-1]
        device = y.device
        dtype = y.dtype

        sigma = self._sample_sigma(B).to(dtype)
        log_sigma = torch.log(sigma)

        # Sample δ, v per-sample uniformly in a per-component box.
        delta = (2.0 * torch.rand(*B, d, device=device, dtype=dtype) - 1.0) \
            * self.shift_scale
        v_vec = (2.0 * torch.rand(*B, d, device=device, dtype=dtype) - 1.0) \
            * self.smear_scale

        # Score at y, for the shift target AND the head's score
        # input. Score requires autograd internally, but we detach
        # the result so head gradients don't propagate back into
        # the base through this channel.
        score_input = self.energy.score(
            y, c, log_sigma, create_graph=False,
        ).detach()

        # Teacher lp values (simple forwards — no autograd needed).
        with torch.no_grad():
            E0 = self.energy(y, c, log_sigma)
            E_shift = self.energy(y - delta, c, log_sigma)
            lp0 = -E0
            lp_shift = -E_shift

            # Shift target: lp(y−δ) − lp(y) + δ·s — residual after
            # subtracting the 1st-order Taylor prefix.
            target_shift = (lp_shift - lp0) + (delta * score_input).sum(dim=-1)

            # Smear target: deterministic log-E-exp with K MC samples.
            K = max(1, self.smear_mc_samples)
            eps_k = torch.randn(K, *B, device=device, dtype=dtype)
            v_exp = v_vec.unsqueeze(0).expand(K, *([-1] * (v_vec.dim())))
            y_pert = y.unsqueeze(0) - eps_k.unsqueeze(-1) * v_exp
            c_exp = c.unsqueeze(0).expand(K, *([-1] * (c.dim())))
            ls_exp = log_sigma.unsqueeze(0).expand(
                K, *([-1] * (log_sigma.dim()))
            )
            flat_y = y_pert.reshape(-1, d)
            flat_c = c_exp.reshape(-1, c.shape[-1])
            flat_ls = ls_exp.reshape(-1)
            flat_E = self.energy(flat_y, flat_c, flat_ls)
            E_pert = flat_E.reshape(K, *B)
            lp_pert = -E_pert
            lp0_bcast = lp0.unsqueeze(0).expand_as(lp_pert)
            dlp_eps = lp_pert - lp0_bcast  # [K, *B]
            target_smear = torch.logsumexp(dlp_eps, dim=0) \
                - torch.log(torch.tensor(float(K), device=device, dtype=dtype))

        # Head forward — this is where gradient flows.
        R_shift_pred, S_smear_pred = self.head(
            y, c, log_sigma, score_input, delta, v_vec,
        )
        shift_loss_per = (R_shift_pred - target_shift).pow(2)
        smear_loss_per = (S_smear_pred - target_smear).pow(2)
        return shift_loss_per, smear_loss_per

    def forward(
        self, y: torch.Tensor, c: torch.Tensor, w: torch.Tensor,
    ) -> Tuple[
        torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor
    ]:
        """Return ``(combined_wsum, dsm_wsum, shift_wsum, smear_wsum, wsum)``.

        Only ``combined_wsum`` carries a gradient; the per-component
        sums are detached for logging. Caller divides by ``wsum`` at
        epoch end to get per-sample averages.
        """
        dsm_per, sigma = self.dsm_loss(
            y, c, create_graph=self.training,
        )
        dsm_weighted = dsm_per
        if self.mle_weight:
            dsm_weighted = sigma.pow(2) * dsm_weighted

        do_shift = self.shift_loss_weight > 0.0
        do_smear = self.smear_loss_weight > 0.0
        if do_shift or do_smear:
            shift_per, smear_per = self.head_losses(y, c)
        else:
            shift_per = torch.zeros_like(dsm_per)
            smear_per = torch.zeros_like(dsm_per)

        combined_per = (
            dsm_weighted
            + self.shift_loss_weight * shift_per
            + self.smear_loss_weight * smear_per
        )
        wsum = w.sum()
        combined_wsum = (w * combined_per).sum()
        dsm_wsum = (w * dsm_per).sum().detach()
        shift_wsum = (w * shift_per).sum().detach()
        smear_wsum = (w * smear_per).sum().detach()
        return combined_wsum, dsm_wsum, shift_wsum, smear_wsum, wsum


def build_joint_model(
    n_features: int,
    n_cond: int,
    hidden_features: int,
    n_hidden_layers: int,
    head_hidden_features: int,
    head_n_hidden_layers: int,
    activation: str,
    sigma_schedule: str,
    sigma_min: float,
    sigma_max: float,
    sigma_dsm: float,
    mle_weight: bool = False,
    shift_scale: float = 0.3,
    smear_scale: float = 0.3,
    shift_loss_weight: float = 1.0,
    smear_loss_weight: float = 1.0,
    smear_mc_samples: int = 4,
) -> JointScoreModel:
    act_cls = ACTIVATIONS.get(activation.lower())
    if act_cls is None:
        raise ValueError(
            f"unknown activation '{activation}'; "
            f"available: {sorted(ACTIVATIONS)}"
        )
    energy = EnergyMLP(
        n_features=n_features,
        n_cond=n_cond,
        hidden_features=hidden_features,
        n_hidden_layers=n_hidden_layers,
        activation=act_cls,
    )
    head = HeadMLP(
        n_features=n_features,
        n_cond=n_cond,
        hidden_features=head_hidden_features,
        n_hidden_layers=head_n_hidden_layers,
        activation=act_cls,
    )
    return JointScoreModel(
        energy=energy,
        head=head,
        sigma_schedule=sigma_schedule,
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        sigma_dsm=sigma_dsm,
        mle_weight=mle_weight,
        shift_scale=shift_scale,
        smear_scale=smear_scale,
        shift_loss_weight=shift_loss_weight,
        smear_loss_weight=smear_loss_weight,
        smear_mc_samples=smear_mc_samples,
    )


def _all_reduce_sum_(t: torch.Tensor, is_dist: bool) -> torch.Tensor:
    if is_dist:
        import torch.distributed as dist
        dist.all_reduce(t, op=dist.ReduceOp.SUM)
    return t


# -----------------------------------------------------------------------------
# Training loop
# -----------------------------------------------------------------------------

def train(
    model: nn.Module,
    inner_model: "JointScoreModel",
    train_loader,
    val_loader,
    epochs: int,
    lr: float,
    weight_decay: float,
    patience: int,
    device: str,
    log_lines: List[str],
    is_dist: bool = False,
    is_rank0: bool = True,
    checkpoint_path=None,
    stats=None,
    model_config=None,
    use_amp: bool = False,
):
    """Train the energy-based score model with DSM.

    Mirrors ``train_muon_response_flow.train`` except that the per-batch
    forward returns ``(weighted_residual_sum, weight_sum)`` directly
    (no NLL conversion needed) and the inner forward already includes
    the double-backward via ``score(..., create_graph=True)``.
    """
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=lr, weight_decay=weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=max(1, patience // 2)
    )

    best_val = float("inf")
    best_state = None
    no_improve = 0

    postfix_every = 50

    amp_device_type = "cuda" if str(device).startswith("cuda") else "cpu"
    amp_ctx = lambda: torch.amp.autocast(
        device_type=amp_device_type,
        dtype=torch.bfloat16,
        enabled=use_amp,
    )

    shift_w = inner_model.shift_loss_weight
    smear_w = inner_model.smear_loss_weight
    mle_weight = inner_model.mle_weight

    def run_val(epoch_idx: int):
        # Evaluate all three loss components for validation. No
        # parameter gradients needed — scores via create_graph=False.
        model.eval()
        total_dsm = torch.zeros((), device=device)
        total_shift = torch.zeros((), device=device)
        total_smear = torch.zeros((), device=device)
        wsum = torch.zeros((), device=device)
        val_bar = tqdm(
            val_loader,
            desc=f"epoch {epoch_idx:03d} val  ",
            leave=False,
            dynamic_ncols=True,
            disable=not is_rank0,
        )
        for x, c, w in val_bar:
            x = x.to(device, non_blocking=True)
            c = c.to(device, non_blocking=True)
            w = w.to(device, non_blocking=True)
            with amp_ctx():
                dsm_per, _ = inner_model.dsm_loss(
                    x, c, create_graph=False,
                )
                if shift_w > 0.0 or smear_w > 0.0:
                    shift_per, smear_per = inner_model.head_losses(x, c)
                else:
                    shift_per = torch.zeros_like(dsm_per)
                    smear_per = torch.zeros_like(dsm_per)
            total_dsm = total_dsm + (w * dsm_per.float()).sum().detach()
            total_shift = total_shift + (w * shift_per.float()).sum().detach()
            total_smear = total_smear + (w * smear_per.float()).sum().detach()
            wsum = wsum + w.sum()
        _all_reduce_sum_(total_dsm, is_dist)
        _all_reduce_sum_(total_shift, is_dist)
        _all_reduce_sum_(total_smear, is_dist)
        _all_reduce_sum_(wsum, is_dist)
        wsum_py = max(wsum.item(), 1e-30)
        return (
            total_dsm.item() / wsum_py,
            total_shift.item() / wsum_py,
            total_smear.item() / wsum_py,
        )

    for epoch in range(1, epochs + 1):
        if is_rank0:
            print(f"epoch {epoch:03d} starting", flush=True)
        model.train()
        t0 = time.time()
        total_combined = torch.zeros((), device=device)
        total_dsm = torch.zeros((), device=device)
        total_shift = torch.zeros((), device=device)
        total_smear = torch.zeros((), device=device)
        wsum = torch.zeros((), device=device)
        bar = tqdm(
            train_loader,
            desc=f"epoch {epoch:03d} train",
            leave=False,
            dynamic_ncols=True,
            disable=not is_rank0,
        )
        for i, (x, c, w) in enumerate(bar):
            x = x.to(device, non_blocking=True)
            c = c.to(device, non_blocking=True)
            w = w.to(device, non_blocking=True)

            optimizer.zero_grad(set_to_none=True)
            with amp_ctx():
                (
                    combined_wsum,
                    dsm_wsum,
                    shift_wsum,
                    smear_wsum,
                    weight_sum,
                ) = model(x, c, w)
            # Keep the loss in FP32 for the backward.
            combined_wsum = combined_wsum.float()
            weight_sum = weight_sum.float().clamp_min(1e-30)
            loss = combined_wsum / weight_sum
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            total_combined = total_combined + combined_wsum.detach()
            total_dsm = total_dsm + dsm_wsum
            total_shift = total_shift + shift_wsum
            total_smear = total_smear + smear_wsum
            wsum = wsum + weight_sum.detach()
            if is_rank0 and (i + 1) % postfix_every == 0:
                denom = max(wsum.item(), 1e-30)
                bar.set_postfix(
                    dsm=f"{total_dsm.item() / denom:.4f}",
                    shft=f"{total_shift.item() / denom:.4f}",
                    smr=f"{total_smear.item() / denom:.4f}",
                )

        _all_reduce_sum_(total_combined, is_dist)
        _all_reduce_sum_(total_dsm, is_dist)
        _all_reduce_sum_(total_shift, is_dist)
        _all_reduce_sum_(total_smear, is_dist)
        _all_reduce_sum_(wsum, is_dist)
        wsum_py = max(wsum.item(), 1e-30)
        train_combined = total_combined.item() / wsum_py
        train_dsm = total_dsm.item() / wsum_py
        train_shift = total_shift.item() / wsum_py
        train_smear = total_smear.item() / wsum_py
        val_dsm, val_shift, val_smear = run_val(epoch)
        val_combined = (
            val_dsm + shift_w * val_shift + smear_w * val_smear
        )
        scheduler.step(val_combined)
        lr_now = optimizer.param_groups[0]["lr"]
        line = (
            f"epoch {epoch:03d}  "
            f"train_dsm {train_dsm:.4f} "
            f"train_shift {train_shift:.4f} "
            f"train_smear {train_smear:.4f}  "
            f"val_dsm {val_dsm:.4f} "
            f"val_shift {val_shift:.4f} "
            f"val_smear {val_smear:.4f}  "
            f"lr {lr_now:.2e}  dt {time.time()-t0:.1f}s"
        )
        if is_rank0:
            print(line, flush=True)
        log_lines.append(line)

        if val_combined < best_val - 1e-4:
            best_val = val_combined
            if is_rank0:
                best_state = {
                    k: v.detach().cpu().clone()
                    for k, v in inner_model.state_dict().items()
                }
            no_improve = 0
        else:
            no_improve += 1

        if is_rank0 and checkpoint_path is not None:
            current_state = {
                k: v.detach().cpu().clone()
                for k, v in inner_model.state_dict().items()
            }
            ckpt = {
                "epoch": epoch,
                "train_dsm": train_dsm,
                "train_shift": train_shift,
                "train_smear": train_smear,
                "val_dsm": val_dsm,
                "val_shift": val_shift,
                "val_smear": val_smear,
                "best_val": best_val,
                "no_improve": no_improve,
                "state_dict": current_state,
                "stats": asdict(stats) if stats is not None else None,
                "model_config": model_config,
            }
            torch.save(ckpt, checkpoint_path)

        if no_improve >= patience:
            line = (
                f"early stopping at epoch {epoch} "
                f"(best val {best_val:.4f})"
            )
            if is_rank0:
                print(line)
            log_lines.append(line)
            break

    if is_rank0 and best_state is not None:
        inner_model.load_state_dict(best_state)
    return inner_model, best_val


# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

class ScoreWrapper(nn.Module):
    """Raw-coordinate inference interface for base + reweight head.

    All methods take an optional ``sigma`` argument — either a scalar
    tensor or a per-sample tensor. When omitted, the wrapper uses its
    baked-in ``sigma_inference`` buffer.

    Methods:
      ``energy_raw(y_raw, c_raw[, sigma])``
          → ``E_θ(y, c, σ)`` in standardized space.
      ``unnormalized_log_density(y_raw, c_raw[, sigma])``
          → ``-E_θ - Σ log target_std``. Up to a ``(c, σ)`` constant.
      ``score(y_raw, c_raw[, sigma])``
          → ``∇_{y_raw} log p_σ(y_raw | c_raw)``. One forward + one
            backward through the base.
      ``hessian(y_raw, c_raw[, sigma])``
          → ``∂² log p_σ / ∂y_raw²``, shape ``[..., d, d]``. Kept for
            diagnostics; not used by ``apply_shift`` / ``apply_smear``.
      ``apply_shift(y_raw, c_raw, delta_raw[, sigma])``
          → ``Δlp_shift = -δ·s + R_shift``. Expects δ in raw y-units.
      ``apply_smear(y_raw, c_raw, v_raw[, sigma])``
          → ``Δlp_smear = S_smear``. Deterministic rank-1 smear along
            ``v``. Expects v in raw y-units.
      ``apply_combined(y_raw, c_raw, delta_raw, v_raw[, sigma], exact=False)``
          → Combined shift+smear Δlp. ``exact=False`` (default) does
            the naive additive composition in one head call. ``exact
            =True`` evaluates the smear at the shifted position
            ``y − δ`` (two head calls per systematic) to respect the
            full factorization
            ``Δlp = Δlp_shift(y, δ) + Δlp_smear(y − δ, v)``.
    """

    def __init__(
        self,
        energy: EnergyMLP,
        head,
        target_mean: torch.Tensor,
        target_std: torch.Tensor,
        cond_mean: torch.Tensor,
        cond_std: torch.Tensor,
        sigma_inference: float,
    ):
        """``head`` may be ``None`` to load an old checkpoint that
        predates the reweight head. In that case the ``apply_*``
        methods raise; ``score`` / ``hessian`` / ``unnormalized_log_
        density`` still work."""
        super().__init__()
        self.energy = energy
        # Keep head accessible as self.head (None if absent). Assigning
        # nn.Module to an attribute registers it; None is stored as a
        # plain Python attribute.
        self.head = head
        self.register_buffer("target_mean", target_mean)
        self.register_buffer("target_std", target_std)
        self.register_buffer("cond_mean", cond_mean)
        self.register_buffer("cond_std", cond_std)
        self.register_buffer(
            "sigma_inference",
            torch.tensor(float(sigma_inference), dtype=torch.float32),
        )

    def _standardize_target(self, y_raw: torch.Tensor) -> torch.Tensor:
        return (y_raw - self.target_mean) / self.target_std

    def _standardize_cond(self, c_raw: torch.Tensor) -> torch.Tensor:
        return (c_raw - self.cond_mean) / self.cond_std

    def _resolve_log_sigma(
        self, sigma, y_raw: torch.Tensor
    ) -> torch.Tensor:
        if sigma is None:
            sigma_t = self.sigma_inference
        elif isinstance(sigma, torch.Tensor):
            sigma_t = sigma.to(y_raw.device, dtype=y_raw.dtype)
        else:
            sigma_t = torch.tensor(
                float(sigma), device=y_raw.device, dtype=y_raw.dtype
            )
        return torch.log(sigma_t)

    def energy_raw(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        y = self._standardize_target(y_raw)
        c = self._standardize_cond(c_raw)
        log_sigma = self._resolve_log_sigma(sigma, y_raw)
        return self.energy(y, c, log_sigma)

    def unnormalized_log_density(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        return -self.energy_raw(y_raw, c_raw, sigma=sigma) - torch.log(
            self.target_std
        ).sum()

    def score(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        """Return ``∇_{y_raw} log p_σ(y_raw | c_raw)``."""
        y_raw = y_raw.detach().requires_grad_(True)
        log_p = self.unnormalized_log_density(y_raw, c_raw, sigma=sigma)
        (grad,) = torch.autograd.grad(
            log_p.sum(), y_raw, create_graph=False
        )
        return grad

    def hessian(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        """Return ``∂²log p_σ(y_raw | c_raw) / ∂y_raw_i ∂y_raw_j``
        with shape ``[..., d, d]``.

        Differentiates the unnormalized log-density twice wrt
        ``y_raw`` directly — the standardization chain rule and the
        additive constant from the Jacobian correction drop out of
        the second derivative. The result is symmetric by construction
        (second derivatives of a scalar).

        Obtained by autograd from the learned energy. Kept available
        for diagnostics; reweight inference uses the head instead.
        """
        y_raw = y_raw.detach().requires_grad_(True)
        log_p = self.unnormalized_log_density(y_raw, c_raw, sigma=sigma)
        (grad,) = torch.autograd.grad(
            log_p.sum(), y_raw, create_graph=True
        )
        d = y_raw.shape[-1]
        H_rows: List[torch.Tensor] = []
        for i in range(d):
            (row,) = torch.autograd.grad(
                grad[..., i].sum(),
                y_raw,
                retain_graph=(i < d - 1),
                create_graph=False,
            )
            H_rows.append(row)
        return torch.stack(H_rows, dim=-2)

    def _head_outputs(
        self,
        y: torch.Tensor,
        c: torch.Tensor,
        log_sigma: torch.Tensor,
        score_std: torch.Tensor,
        delta_std: torch.Tensor,
        v_std: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Head forward in standardized space (no-grad-friendly)."""
        if self.head is None:
            raise RuntimeError(
                "ScoreWrapper was built without a head; "
                "apply_shift/apply_smear/apply_combined are unavailable."
            )
        return self.head(y, c, log_sigma, score_std, delta_std, v_std)

    def _score_standardized(
        self, y_raw: torch.Tensor, c_raw: torch.Tensor, sigma,
    ) -> Tuple[
        torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor
    ]:
        """Compute the standardized-space score at ``y_raw`` plus the
        standardized ``(y, c, log σ)`` needed by the head.

        The head was trained with standardized inputs, so we convert
        here: ``y_std = (y_raw − mean)/std`` and the score in that
        space is ``s_std = std ⊙ s_raw``.
        """
        y_std0 = self._standardize_target(y_raw)
        c_std = self._standardize_cond(c_raw)
        log_sigma = self._resolve_log_sigma(sigma, y_raw)
        # Autograd score in raw space, then rescale to std space.
        y_inp = y_raw.detach().requires_grad_(True)
        E = self.energy_raw(y_inp, c_raw, sigma=sigma)
        (grad_raw,) = torch.autograd.grad(
            E.sum(), y_inp, create_graph=False
        )
        score_raw = -grad_raw
        score_std = score_raw * self.target_std
        return y_std0.detach(), c_std.detach(), log_sigma, score_std.detach()

    def apply_shift(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        delta_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        """``Δlp = log p_σ(y − δ | c) − log p_σ(y | c)`` via the head.

        Expects ``delta_raw`` in *raw* y-units — it is standardized
        internally. Returns a scalar per sample (no trailing dim).
        """
        y_std, c_std, log_sigma, score_std = self._score_standardized(
            y_raw, c_raw, sigma,
        )
        delta_std = delta_raw / self.target_std
        zero_v = torch.zeros_like(delta_std)
        R_shift, _ = self._head_outputs(
            y_std, c_std, log_sigma, score_std, delta_std, zero_v,
        )
        taylor = -(delta_std * score_std).sum(dim=-1)
        return taylor + R_shift

    def apply_smear(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        v_raw: torch.Tensor,
        sigma=None,
    ) -> torch.Tensor:
        """Deterministic rank-1 smear ``Δlp = log E_ε[p(y − εv)/p(y)]``.

        ``v_raw`` is in raw y-units; standardized internally.
        """
        y_std, c_std, log_sigma, score_std = self._score_standardized(
            y_raw, c_raw, sigma,
        )
        v_std = v_raw / self.target_std
        zero_d = torch.zeros_like(v_std)
        _, S_smear = self._head_outputs(
            y_std, c_std, log_sigma, score_std, zero_d, v_std,
        )
        return S_smear

    def apply_combined(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        delta_raw: torch.Tensor,
        v_raw: torch.Tensor,
        sigma=None,
        exact: bool = False,
    ) -> torch.Tensor:
        """Combined shift + rank-1 smear reweight in lp space.

        ``exact=False`` (default): naive additive composition. One
        head call with both δ and v nonzero; cheapest. Composition
        error is ``O(|δ|·|v|²)``, negligible for small systematics.

        ``exact=True``: two head calls per systematic. Evaluates
        ``Δlp_smear`` at the shifted position ``y − δ`` so the full
        factorization ``Δlp = Δlp_shift(y, δ) + Δlp_smear(y − δ, v)``
        is respected. Use for large-δ systematics.
        """
        if exact:
            shift_part = self.apply_shift(y_raw, c_raw, delta_raw, sigma)
            smear_part = self.apply_smear(
                y_raw - delta_raw, c_raw, v_raw, sigma,
            )
            return shift_part + smear_part
        y_std, c_std, log_sigma, score_std = self._score_standardized(
            y_raw, c_raw, sigma,
        )
        delta_std = delta_raw / self.target_std
        v_std = v_raw / self.target_std
        R_shift, S_smear = self._head_outputs(
            y_std, c_std, log_sigma, score_std, delta_std, v_std,
        )
        taylor = -(delta_std * score_std).sum(dim=-1)
        return taylor + R_shift + S_smear


def export_model(
    energy: EnergyMLP,
    head: HeadMLP,
    stats: PreprocStats,
    outpath_pt: str,
    outpath_ts: str,
    model_config: dict,
    sigma_inference: float,
):
    """Save the trained wrapper (base + head) two ways.

    1. ``<outpath_pt>``: ``torch.save`` of ScoreWrapper state_dict +
       config. Loading requires Python + torch.
    2. ``<outpath_ts>``: ``torch.jit.trace`` of the energy forward
       only. The head is saved in state_dict but not traced.
    """
    energy = energy.cpu().eval()
    head = head.cpu().eval()
    wrapper = ScoreWrapper(
        energy=energy,
        head=head,
        target_mean=torch.tensor(stats.target_mean, dtype=torch.float32),
        target_std=torch.tensor(stats.target_std, dtype=torch.float32),
        cond_mean=torch.tensor(stats.cond_mean, dtype=torch.float32),
        cond_std=torch.tensor(stats.cond_std, dtype=torch.float32),
        sigma_inference=float(sigma_inference),
    )

    torch.save(
        {
            "wrapper_state_dict": wrapper.state_dict(),
            "model_config": model_config,
            "preproc": asdict(stats),
            "sigma_inference": float(sigma_inference),
        },
        outpath_pt,
    )
    print(
        f"saved wrapper state_dict + config to {outpath_pt} "
        f"(sigma_inference={sigma_inference:.4g})"
    )

    n_features = len(stats.target_mean)
    n_cond = len(stats.cond_mean)
    x_example = torch.randn(1, n_features, dtype=torch.float32)
    c_example = torch.randn(1, n_cond, dtype=torch.float32)
    sigma_example = torch.tensor([sigma_inference], dtype=torch.float32)

    class _TracedEnergy(nn.Module):
        def __init__(self, inner: ScoreWrapper):
            super().__init__()
            self.inner = inner

        def forward(self, y_raw, c_raw, sigma):
            return self.inner.energy_raw(y_raw, c_raw, sigma=sigma)

    try:
        with torch.no_grad():
            traced = torch.jit.trace(
                _TracedEnergy(wrapper),
                (x_example, c_example, sigma_example),
                check_trace=False,
                strict=False,
            )
        traced.save(outpath_ts)
        print(
            f"saved TorchScript (traced, energy only) to {outpath_ts}"
        )
    except Exception as e:
        print(
            f"[note] torch.jit.trace export skipped "
            f"({type(e).__name__}: {e}). "
            f"The .pt wrapper file is the portable alternative."
        )


# -----------------------------------------------------------------------------
# In-memory loader
# -----------------------------------------------------------------------------

class InMemoryLoader:
    """Bulk-indexed per-batch iterator over in-RAM CPU tensors."""

    def __init__(self, x, c, w, batch_size, shuffle, drop_last):
        self.x, self.c, self.w = x, c, w
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.drop_last = drop_last
        self.n = x.shape[0]

    def __len__(self):
        if self.drop_last:
            return self.n // self.batch_size
        return (self.n + self.batch_size - 1) // self.batch_size

    def __iter__(self):
        bs = self.batch_size
        if self.shuffle:
            perm = torch.randperm(self.n)
            n_batches = self.n // bs if self.drop_last else (
                (self.n + bs - 1) // bs
            )
            for i in range(n_batches):
                idx = perm[i * bs:min((i + 1) * bs, self.n)]
                yield (
                    self.x.index_select(0, idx),
                    self.c.index_select(0, idx),
                    self.w.index_select(0, idx),
                )
        else:
            n_batches = self.n // bs if self.drop_last else (
                (self.n + bs - 1) // bs
            )
            for i in range(n_batches):
                s = i * bs
                e = min(s + bs, self.n)
                yield self.x[s:e], self.c[s:e], self.w[s:e]


# -----------------------------------------------------------------------------
# Worker
# -----------------------------------------------------------------------------

def main_worker(
    rank,
    args,
    world_size,
    master_port,
    stats,
    model_config,
    target_std_t,
    cond_t,
    w_t,
    train_sel,
    val_sel,
):
    is_dist = world_size > 1
    is_rank0 = rank == 0

    if is_dist:
        import torch.distributed as dist
        os.environ["MASTER_ADDR"] = "localhost"
        os.environ["MASTER_PORT"] = str(master_port)
        if args.device.startswith("cuda"):
            torch.cuda.set_device(rank)
            device = f"cuda:{rank}"
            backend = "nccl"
        else:
            device = "cpu"
            backend = "gloo"
        dist.init_process_group(
            backend=backend, rank=rank, world_size=world_size
        )
    else:
        if args.device.startswith("cuda"):
            torch.cuda.set_device(0)
            device = "cuda:0"
        else:
            device = args.device

    if is_dist:
        train_sel_rank = torch.chunk(
            train_sel, world_size
        )[rank].contiguous()
        val_sel_rank = torch.chunk(
            val_sel, world_size
        )[rank].contiguous()
    else:
        train_sel_rank = train_sel
        val_sel_rank = val_sel

    train_x_cpu = target_std_t.index_select(0, train_sel_rank).contiguous()
    train_c_cpu = cond_t.index_select(0, train_sel_rank).contiguous()
    train_w_cpu = w_t.index_select(0, train_sel_rank).contiguous()
    val_x_cpu = target_std_t.index_select(0, val_sel_rank).contiguous()
    val_c_cpu = cond_t.index_select(0, val_sel_rank).contiguous()
    val_w_cpu = w_t.index_select(0, val_sel_rank).contiguous()
    n_train_rank = train_x_cpu.shape[0]
    n_val_rank = val_x_cpu.shape[0]

    train_loader = InMemoryLoader(
        train_x_cpu, train_c_cpu, train_w_cpu,
        batch_size=args.batch_size, shuffle=True, drop_last=True,
    )
    val_loader = InMemoryLoader(
        val_x_cpu, val_c_cpu, val_w_cpu,
        batch_size=args.batch_size, shuffle=False, drop_last=False,
    )
    if is_rank0:
        print(
            f"rank {rank}/{world_size}  "
            f"train {n_train_rank}  val {n_val_rank}  "
            f"batch {args.batch_size}  device {device}"
        )

    inner_model = build_joint_model(
        n_features=model_config["n_features"],
        n_cond=model_config["n_cond"],
        hidden_features=model_config["hidden_features"],
        n_hidden_layers=model_config["n_hidden_layers"],
        head_hidden_features=model_config["head_hidden_features"],
        head_n_hidden_layers=model_config["head_n_hidden_layers"],
        activation=model_config.get("activation", "gelu"),
        sigma_schedule=model_config.get("sigma_schedule", "log-uniform"),
        sigma_min=model_config.get("sigma_min", 0.01),
        sigma_max=model_config.get("sigma_max", 0.3),
        sigma_dsm=model_config["sigma_dsm"],
        mle_weight=bool(model_config.get("mle_weight", False)),
        shift_scale=float(model_config.get("shift_scale", 0.3)),
        smear_scale=float(model_config.get("smear_scale", 0.3)),
        shift_loss_weight=float(model_config.get("shift_loss_weight", 1.0)),
        smear_loss_weight=float(model_config.get("smear_loss_weight", 1.0)),
        smear_mc_samples=int(model_config.get("smear_mc_samples", 4)),
    ).to(device)

    if is_dist:
        if args.device.startswith("cuda"):
            model = nn.parallel.DistributedDataParallel(
                inner_model, device_ids=[rank]
            )
        else:
            model = nn.parallel.DistributedDataParallel(inner_model)
    else:
        model = inner_model

    n_energy_params = sum(p.numel() for p in inner_model.energy.parameters())
    n_head_params = sum(p.numel() for p in inner_model.head.parameters())
    n_params = n_energy_params + n_head_params
    if is_rank0:
        print(
            f"parameters: energy {n_energy_params:,}  "
            f"head {n_head_params:,}  total {n_params:,}"
        )

    log_lines: List[str] = []
    if is_rank0:
        log_lines.append(
            f"parameters: energy {n_energy_params} head {n_head_params} "
            f"total {n_params}"
        )
        log_lines.append(
            f"train {n_train_rank * world_size} "
            f"val {n_val_rank * world_size} "
            f"batch {args.batch_size}  world_size {world_size}"
        )

    checkpoint_path = (
        os.path.join(args.output, "checkpoint.pt") if is_rank0 else None
    )
    use_amp = bool(args.amp)
    if is_rank0:
        print(f"bf16 autocast: {'on' if use_amp else 'off'}")
        schedule = model_config.get("sigma_schedule", "log-uniform")
        print(f"sigma schedule: {schedule}")
        if schedule == "log-uniform":
            print(
                f"  sigma range: [{model_config['sigma_min']:.4g}, "
                f"{model_config['sigma_max']:.4g}]"
            )
            print(
                f"  sigma_inference: {model_config['sigma_inference']:.4g}"
            )
        else:
            print(f"  sigma_dsm:  {model_config['sigma_dsm']:.4g}")
        print(
            f"mle weight: "
            f"{'on' if model_config.get('mle_weight', False) else 'off'}"
        )
        print(
            f"head: {model_config['head_hidden_features']} × "
            f"{model_config['head_n_hidden_layers']}  "
            f"shift_w={model_config['shift_loss_weight']:.2g}  "
            f"smear_w={model_config['smear_loss_weight']:.2g}  "
            f"K={model_config['smear_mc_samples']}"
        )
        print(
            f"training δ ~ U(±{model_config['shift_scale']:.2g}), "
            f"v ~ U(±{model_config['smear_scale']:.2g})"
        )

    trained_model, best_val = train(
        model,
        inner_model,
        train_loader,
        val_loader,
        args.epochs,
        args.lr,
        args.weight_decay,
        args.patience,
        device,
        log_lines,
        is_dist=is_dist,
        is_rank0=is_rank0,
        checkpoint_path=checkpoint_path,
        stats=stats,
        model_config=model_config,
        use_amp=use_amp,
    )

    if is_rank0:
        log_lines.append(f"best val_combined: {best_val:.4f}")
        with open(os.path.join(args.output, "training.log"), "w") as f:
            f.write("\n".join(log_lines) + "\n")
        export_model(
            trained_model.energy,
            trained_model.head,
            stats,
            outpath_pt=os.path.join(args.output, "score_model.pt"),
            outpath_ts=os.path.join(args.output, "score_scripted.pt"),
            model_config=model_config,
            sigma_inference=model_config["sigma_inference"],
        )
        print(f"best val_combined: {best_val:.4f}")
        print("done")

    if is_dist:
        import torch.distributed as dist
        dist.destroy_process_group()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    print(f"loading ntuples from {len(args.input_files)} file(s)")
    (
        pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q, w,
    ) = load_ntuples(
        args.input_files,
        args.tree,
        args.max_muons,
        args.pt_min,
        args.pt_max,
        args.eta_max,
        threads=args.threads,
    )

    target, cond_raw = compute_targets_and_conditioning(
        pt_r, eta_r, phi_r, pt_g, eta_g, phi_g, q
    )

    w = (w / w.mean()).astype(np.float32)
    print(
        f"target  r_kappa mean={target[:,0].mean():+.3e} "
        f"std={target[:,0].std():.3e}"
    )
    print(
        f"         dlambda mean={target[:,1].mean():+.3e} "
        f"std={target[:,1].std():.3e}"
    )
    print(
        f"         dphi    mean={target[:,2].mean():+.3e} "
        f"std={target[:,2].std():.3e}"
    )

    stats = build_preproc(target, cond_raw)
    target_std, cond = apply_preproc(target, cond_raw, stats)

    with open(os.path.join(args.output, "preproc.json"), "w") as f:
        json.dump(asdict(stats), f, indent=2)

    n_total = target_std.shape[0]
    n_features = target_std.shape[1]
    n_cond = cond.shape[1]
    n_val = int(n_total * args.val_fraction)
    n_train = n_total - n_val
    gen = torch.Generator().manual_seed(args.seed)
    perm_all = torch.randperm(n_total, generator=gen)
    val_sel = perm_all[:n_val].contiguous()
    train_sel = perm_all[n_val:].contiguous()

    target_std_t = torch.from_numpy(target_std).contiguous()
    cond_t = torch.from_numpy(cond).contiguous()
    w_t = torch.from_numpy(w).contiguous()
    del target_std, cond, w

    if args.device.startswith("cuda"):
        if args.num_gpus == -1:
            world_size = max(1, torch.cuda.device_count())
        else:
            world_size = max(1, args.num_gpus)
    else:
        world_size = 1

    if args.batch_size is None:
        args.batch_size = 32768 if args.device != "cpu" else 16384

    sigma_inference = (
        float(args.sigma_inference)
        if args.sigma_inference is not None
        else float(args.sigma_min)
    )
    model_config = {
        "model_type": "EnergyMLP+HeadMLP",
        "n_features": int(n_features),
        "n_cond": int(n_cond),
        "hidden_features": args.hidden_features,
        "n_hidden_layers": args.n_hidden_layers,
        "head_hidden_features": args.head_hidden_features,
        "head_n_hidden_layers": args.head_n_hidden_layers,
        "activation": args.activation,
        "sigma_schedule": args.sigma_schedule,
        "sigma_min": float(args.sigma_min),
        "sigma_max": float(args.sigma_max),
        "sigma_inference": sigma_inference,
        "sigma_dsm": float(args.sigma_dsm),
        "mle_weight": bool(args.mle_weight),
        "shift_scale": float(args.shift_scale),
        "smear_scale": float(args.smear_scale),
        "shift_loss_weight": float(args.shift_loss_weight),
        "smear_loss_weight": float(args.smear_loss_weight),
        "smear_mc_samples": int(args.smear_mc_samples),
        "loss": "dsm+shift+smear",
    }

    print(
        f"world_size {world_size}  batch_size {args.batch_size}  "
        f"device {args.device}"
    )
    print(f"train {n_train}  val {n_val}  total {n_total}")

    if world_size == 1:
        main_worker(
            0, args, 1, 0, stats, model_config,
            target_std_t, cond_t, w_t, train_sel, val_sel,
        )
    else:
        target_std_t.share_memory_()
        cond_t.share_memory_()
        w_t.share_memory_()
        train_sel.share_memory_()
        val_sel.share_memory_()
        import socket
        sock = socket.socket()
        sock.bind(("", 0))
        master_port = sock.getsockname()[1]
        sock.close()
        import torch.multiprocessing as mp
        mp.spawn(
            main_worker,
            args=(
                args, world_size, master_port, stats, model_config,
                target_std_t, cond_t, w_t, train_sel, val_sel,
            ),
            nprocs=world_size,
            join=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
