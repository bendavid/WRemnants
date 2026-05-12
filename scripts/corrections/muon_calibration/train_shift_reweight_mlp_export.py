"""Export a trained shift-reweight MLP for C++ inference.

Loads a checkpoint produced by ``train_shift_reweight_mlp.py`` and
packages the inference forward in two formats:

  * **ONNX** (single dynamic-batch file) — runtime: onnxruntime.
  * **AOTInductor** ``.pt2`` (one per requested static batch size) —
    runtime: libtorch's ``AOTIModelPackageLoader``.

Inference signature (both formats):

    forward(y_raw[B, n_features], c_raw[B, n_cond], u_raw[B, n_features])
        -> log_r[B]

where ``y_raw`` and ``c_raw`` are in physical units (the wrapper
standardizes internally using the saved PreprocStats), and ``u_raw``
is a y-space shift in the same physical units as ``y_raw`` (also
standardized internally as ``u_std = u_raw / target_std``).

The wrapper:
  1. Standardizes (y, c) into (y_std, c_std) using baked buffers.
  2. Standardizes u into u_std (no mean subtraction — it's a delta).
  3. Runs the vector-output MLP once: g(y_std, c_std, u_std) ∈
     ℝ^{n_features}.
  4. Projects: ``d = u_std · g``.
  5. Applies the positivity wrap (exp clamp or softplus/log 2)
     consistent with how the model was trained, and returns ``log r``.

Single MLP forward per query — the ``uᵀ·g`` parameterization enforces
``log r = 0`` at ``u = 0`` structurally, with no anchor-subtraction
overhead.

Static-batch AOTI (the per-event narf use case): mkldnn weight-prepack
fusion + freezing + max_autotune produce the lowest-latency
single-event path. Compile multiple .pt2 files via
``--aoti-batches 1 8 32 128 512 2048`` and benchmark each with the
companion C++ bench. ONNX gets a single dynamic-batch file usable at
any runtime size — typically faster than AOTI at large batch but
slower at B=1.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import List

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from train_muon_response_flow import PreprocStats  # noqa: E402
from train_shift_reweight_mlp import (  # noqa: E402
    _ACTIVATIONS,
    _LOG_LOG2,
    _LOG_W_CLAMP,
    ShiftReweightMLP,
)


# -----------------------------------------------------------------------------
# Inference wrapper — bakes preprocessing + positivity, returns log r
# -----------------------------------------------------------------------------

class ShiftReweightInference(nn.Module):
    """``forward(y_raw, c_raw, u_raw) → log_r``.

    Single MLP forward per query. The shift weight is computed as

        log r = positivity( uᵀ_std · g(y_std, c_std, u_std) )

    where ``g`` is the vector-output MLP (n_features outputs). The
    ``u·g`` projection makes ``log r = 0`` exact at ``u = 0`` while
    leaving the score ``∂_u log r |_{u=0} = g(y, c, 0)`` unconstrained.

    Pure tensor ops; ``torch.export.export`` and ``torch.onnx.export``
    both trace this directly.

    Convention:
      * ``y_raw``, ``c_raw`` are in the physical units of the input
        ntuples (the same units the trainer's ``compute_targets_and_
        conditioning`` produces). The wrapper standardizes via baked
        buffers `target_mean / target_std / cond_mean / cond_std`.
      * ``u_raw`` is a y-space shift in the same units as ``y_raw``.
        Standardized as ``u_std = u_raw / target_std`` (no mean
        subtraction — it's a delta).
    """

    def __init__(
        self,
        mlp: ShiftReweightMLP,
        target_mean: torch.Tensor,
        target_std: torch.Tensor,
        cond_mean: torch.Tensor,
        cond_std: torch.Tensor,
        positivity: str = "exp",
        clamp: float = _LOG_W_CLAMP,
    ):
        super().__init__()
        if positivity not in ("exp", "softplus"):
            raise ValueError(f"unknown positivity {positivity!r}")
        self.mlp = mlp
        self.positivity = str(positivity)
        self.clamp = float(clamp)
        self.register_buffer("target_mean", target_mean)
        self.register_buffer("target_std", target_std)
        self.register_buffer("cond_mean", cond_mean)
        self.register_buffer("cond_std", cond_std)

    def forward(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
        u_raw: torch.Tensor,
    ) -> torch.Tensor:
        y_std = (y_raw - self.target_mean) / self.target_std
        c_std = (c_raw - self.cond_mean) / self.cond_std
        u_std = u_raw / self.target_std
        g = self.mlp(y_std, c_std, u_std)        # [B, n_features]
        d = (u_std * g).sum(dim=-1)              # [B]
        if self.positivity == "exp":
            return d.clamp(min=-self.clamp, max=self.clamp)
        d_clamp = d.clamp(min=-self.clamp, max=self.clamp)
        return torch.log(F.softplus(d_clamp).clamp_min(1e-38)) - _LOG_LOG2


# -----------------------------------------------------------------------------
# Checkpoint -> wrapper
# -----------------------------------------------------------------------------

def _load_stats_from_ckpt_or_dir(ckpt: dict, ckpt_path: str) -> PreprocStats:
    """Same logic as the diagnostic script: prefer ``ckpt['stats']``,
    fall back to ``preproc.json`` next to the checkpoint."""
    s = ckpt.get("stats")
    if s is not None:
        return PreprocStats(**s)
    pj = os.path.join(
        os.path.dirname(os.path.abspath(ckpt_path)), "preproc.json",
    )
    if not os.path.exists(pj):
        raise SystemExit(
            f"Checkpoint has no PreprocStats and no preproc.json at "
            f"{pj}. Cannot apply training-time preprocessing."
        )
    with open(pj) as f:
        return PreprocStats(**json.load(f))


def build_wrapper(checkpoint_path: str) -> tuple:
    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    cfg = ckpt["model_config"]
    train_cfg = ckpt.get("train_config", {})
    activation_name = str(cfg.get("activation", "gelu")).lower()
    mlp = ShiftReweightMLP(
        n_features=int(cfg["n_features"]),
        n_cond=int(cfg["n_cond"]),
        hidden_features=int(cfg["hidden_features"]),
        n_hidden_layers=int(cfg["n_hidden_layers"]),
        activation=_ACTIVATIONS.get(activation_name, nn.GELU),
    )
    mlp.load_state_dict(ckpt["state_dict"])
    mlp.eval()

    stats = _load_stats_from_ckpt_or_dir(ckpt, checkpoint_path)
    positivity = str(train_cfg.get("positivity", "exp"))

    wrapper = ShiftReweightInference(
        mlp=mlp,
        target_mean=torch.tensor(stats.target_mean, dtype=torch.float32),
        target_std=torch.tensor(stats.target_std, dtype=torch.float32),
        cond_mean=torch.tensor(stats.cond_mean, dtype=torch.float32),
        cond_std=torch.tensor(stats.cond_std, dtype=torch.float32),
        positivity=positivity,
    ).eval()
    return wrapper, stats, train_cfg


# -----------------------------------------------------------------------------
# ONNX export
# -----------------------------------------------------------------------------

def export_onnx(
    wrapper: ShiftReweightInference,
    output_path: str,
    n_features: int,
    n_cond: int,
    batch_example: int = 1,
    dynamic_batch: bool = True,
    opset: int = 17,
    validate: bool = True,
    validate_n: int = 256,
):
    print(f"\n=== ONNX export -> {output_path} ===")
    y = torch.randn(batch_example, n_features)
    c = torch.randn(batch_example, n_cond)
    u = torch.randn(batch_example, n_features) * 0.1
    with torch.no_grad():
        ref = wrapper(y, c, u)
    print(f"  eager log_r shape: {tuple(ref.shape)}")

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    dynamic_axes = None
    if dynamic_batch:
        dynamic_axes = {
            "y_raw": {0: "batch"},
            "c_raw": {0: "batch"},
            "u_raw": {0: "batch"},
            "log_r": {0: "batch"},
        }
    torch.onnx.export(
        wrapper,
        (y, c, u),
        output_path,
        input_names=["y_raw", "c_raw", "u_raw"],
        output_names=["log_r"],
        dynamic_axes=dynamic_axes,
        opset_version=opset,
    )
    size_kb = os.path.getsize(output_path) / 1024.0
    print(f"  wrote {output_path} ({size_kb:.1f} KB)")

    if not validate:
        return
    try:
        import onnxruntime as ort
    except ImportError:
        print("  [validate] onnxruntime not importable; skipping")
        return
    sess = ort.InferenceSession(
        output_path, providers=["CPUExecutionProvider"],
    )
    rng = np.random.default_rng(0)
    n = max(batch_example, validate_n)
    yv = rng.standard_normal((n, n_features)).astype(np.float32)
    cv = rng.standard_normal((n, n_cond)).astype(np.float32)
    uv = (rng.standard_normal((n, n_features)) * 0.1).astype(np.float32)
    with torch.no_grad():
        ref_v = wrapper(
            torch.from_numpy(yv),
            torch.from_numpy(cv),
            torch.from_numpy(uv),
        ).cpu().numpy()
    out = sess.run(
        ["log_r"], {"y_raw": yv, "c_raw": cv, "u_raw": uv},
    )[0]
    max_abs = float(np.abs(ref_v - out).max())
    denom = max(float(np.abs(ref_v).max()), 1e-30)
    rel = max_abs / denom
    status = "OK" if rel < 1e-4 else "MISMATCH"
    print(
        f"  [validate] max |eager-onnx|={max_abs:.3e} rel={rel:.2e} "
        f"{status}"
    )


# -----------------------------------------------------------------------------
# AOTI export (one .pt2 per static batch size)
# -----------------------------------------------------------------------------

def _to_tuple(out):
    if isinstance(out, (tuple, list)):
        return tuple(out)
    return (out,)


def export_aoti(
    wrapper: ShiftReweightInference,
    output_path: str,
    n_features: int,
    n_cond: int,
    B: int,
    validate: bool = True,
    validate_n: int = 64,
    bench: bool = False,
    bench_duration: float = 2.0,
):
    print(f"\n=== AOTI export -> {output_path} (B={B}) ===")
    y = torch.randn(B, n_features)
    c = torch.randn(B, n_cond)
    u = torch.randn(B, n_features) * 0.1

    with torch.no_grad():
        ref = wrapper(y, c, u)
    print(f"  eager log_r shape: {tuple(ref.shape)}")

    print(f"  torch.export.export (static B={B}) ...")
    ep = torch.export.export(wrapper, (y, c, u))
    print("    OK")

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"  aoti_compile_and_package -> {output_path}")
    with torch.no_grad():
        package_path = torch._inductor.aoti_compile_and_package(
            ep, package_path=output_path,
        )
    size_mb = os.path.getsize(package_path) / 1e6
    print(f"    OK ({size_mb:.2f} MB)")

    runner = torch._inductor.aoti_load_package(package_path)
    out = _to_tuple(runner(y, c, u))[0]
    diff = (out - ref).abs().max().item()
    print(f"  loaded; max |aoti − eager| log_r = {diff:.2e}")

    if validate and validate_n > 0:
        rng = np.random.default_rng(1)
        n_events = max(B, ((validate_n + B - 1) // B) * B)
        yy = torch.from_numpy(
            rng.standard_normal((n_events, n_features)).astype(np.float32))
        cc = torch.from_numpy(
            rng.standard_normal((n_events, n_cond)).astype(np.float32))
        uu = torch.from_numpy(
            (rng.standard_normal((n_events, n_features)) * 0.1)
            .astype(np.float32))
        eager_chunks, aoti_chunks = [], []
        with torch.no_grad():
            for s in range(0, n_events, B):
                yi, ci, ui = yy[s:s + B], cc[s:s + B], uu[s:s + B]
                eager_chunks.append(_to_tuple(wrapper(yi, ci, ui))[0])
                aoti_chunks.append(_to_tuple(runner(yi, ci, ui))[0])
        et = torch.cat(eager_chunks)
        at = torch.cat(aoti_chunks)
        max_abs = (et - at).abs().max().item()
        rel = max_abs / max(et.abs().max().item(), 1e-30)
        status = "OK" if rel < 1e-4 else "MISMATCH"
        print(
            f"  [validate] log_r: max |eager-aoti|={max_abs:.3e} "
            f"rel={rel:.2e} {status}"
        )

    if bench:
        try:
            os.sched_setaffinity(0, {0})
        except Exception:
            pass
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        for _ in range(20):
            runner(y, c, u)
        n_iters = 0
        t0 = time.perf_counter()
        deadline = t0 + bench_duration
        while True:
            runner(y, c, u)
            n_iters += 1
            if n_iters >= 5 and time.perf_counter() >= deadline:
                break
        elapsed = time.perf_counter() - t0
        ms_b = 1000.0 * elapsed / n_iters
        us_e = 1e6 * elapsed / n_iters / B
        ev_s = n_iters * B / elapsed
        print(
            f"  bench (1 core): {ms_b:.3f} ms/call  {us_e:.3f} us/event  "
            f"{ev_s:.1f} ev/s  ({n_iters} iters)"
        )


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--checkpoint", required=True)
    p.add_argument(
        "--output-dir", default=None,
        help="Output directory for ONNX + AOTI files. Default: "
        "directory containing --checkpoint.",
    )
    p.add_argument(
        "--export-onnx", action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to write a dynamic-batch ONNX file. Default: yes.",
    )
    p.add_argument(
        "--onnx-name", default="shift_reweight_mlp.onnx",
        help="ONNX filename (within --output-dir).",
    )
    p.add_argument("--onnx-opset", type=int, default=17)
    p.add_argument("--onnx-batch-example", type=int, default=1)
    p.add_argument(
        "--dynamic-batch", action=argparse.BooleanOptionalAction,
        default=True,
        help="Mark batch dim as dynamic in ONNX. Default: yes.",
    )
    p.add_argument(
        "--export-aoti", action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to write AOTI .pt2 files at the requested static "
        "batch sizes. Default: yes.",
    )
    p.add_argument(
        "--aoti-batches", nargs="+", type=int,
        default=[1, 8, 32, 128, 512, 2048],
        help="Static batch sizes to compile AOTI .pt2 packages for.",
    )
    p.add_argument(
        "--aoti-name-template", default="shift_reweight_mlp_b{batch}.pt2",
        help="AOTI filename template (within --output-dir). "
        "{batch} is substituted with the static batch size.",
    )
    p.add_argument(
        "--max-autotune", action=argparse.BooleanOptionalAction,
        default=True,
        help="Inductor max-autotune for GEMMs (slower compile, "
        "fastest runtime).",
    )
    p.add_argument(
        "--freeze", action=argparse.BooleanOptionalAction, default=True,
        help="Constant-fold weights (required for oneDNN paths).",
    )
    p.add_argument(
        "--prepack", action=argparse.BooleanOptionalAction, default=True,
        help="Prepack Linear weights into oneDNN's preferred layout. "
        "Only fires for static shapes.",
    )
    p.add_argument(
        "--mkldnn-fp32-bf16", action="store_true",
        help="Route fp32 matmul through mkldnn's BF16 path (set_"
        "float32_matmul_precision='medium').",
    )
    p.add_argument(
        "--validate", action=argparse.BooleanOptionalAction, default=True,
        help="Compare ONNX/AOTI outputs against eager on random "
        "inputs after export.",
    )
    p.add_argument("--validate-n", type=int, default=256)
    p.add_argument(
        "--bench", action="store_true",
        help="Microbenchmark each AOTI .pt2 from Python (single core).",
    )
    p.add_argument("--bench-duration", type=float, default=2.0)
    return p.parse_args()


def main():
    args = parse_args()
    if args.output_dir is None:
        args.output_dir = os.path.dirname(os.path.abspath(args.checkpoint))
    os.makedirs(args.output_dir, exist_ok=True)

    if args.mkldnn_fp32_bf16:
        torch.set_float32_matmul_precision("medium")

    from torch._inductor import config as ind_cfg
    if args.freeze:
        ind_cfg.freezing = True
    if args.prepack:
        ind_cfg.cpp.weight_prepack = True
    if args.max_autotune:
        ind_cfg.max_autotune = True
        ind_cfg.max_autotune_gemm = True
    print(
        f"inductor: freezing={ind_cfg.freezing} "
        f"weight_prepack={ind_cfg.cpp.weight_prepack} "
        f"max_autotune={getattr(ind_cfg, 'max_autotune', False)}"
    )

    print(f"loading checkpoint {args.checkpoint}")
    wrapper, stats, train_cfg = build_wrapper(args.checkpoint)
    n_features = len(stats.target_mean)
    n_cond = len(stats.cond_mean)
    print(
        f"  n_features={n_features} n_cond={n_cond} "
        f"positivity={wrapper.positivity}"
    )

    if args.export_onnx:
        export_onnx(
            wrapper,
            output_path=os.path.join(args.output_dir, args.onnx_name),
            n_features=n_features,
            n_cond=n_cond,
            batch_example=args.onnx_batch_example,
            dynamic_batch=args.dynamic_batch,
            opset=args.onnx_opset,
            validate=args.validate,
            validate_n=args.validate_n,
        )

    if args.export_aoti:
        for B in args.aoti_batches:
            out = os.path.join(
                args.output_dir,
                args.aoti_name_template.format(batch=B),
            )
            export_aoti(
                wrapper,
                output_path=out,
                n_features=n_features,
                n_cond=n_cond,
                B=B,
                validate=args.validate,
                validate_n=max(B, args.validate_n),
                bench=args.bench,
                bench_duration=args.bench_duration,
            )

    print("\ndone.")


if __name__ == "__main__":
    main()
