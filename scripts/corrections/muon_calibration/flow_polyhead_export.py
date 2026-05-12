"""AOTI export of a (flow + polyhead) checkpoint.

Loads a checkpoint produced by ``train_muon_response_flow.py
--polyhead`` and packages the inference forward as a self-contained
``.pt2`` (AOT-Inductor) file. The packaged forward signature is

    forward(y_raw, c_raw)
        -> joint_coefs

i.e., per-event polyhead joint-polynomial coefficients only — the
flow is **not** part of the inference graph (the polyhead's
prediction has no flow-state dependence, so a flow forward at
inference would be wasted compute). The downstream forms y-space
perturbations in standardized target units —
``u_shift_std = u_shift_raw / target_std`` and
``σ_vec_std = σ_vec_raw / target_std`` (either / both can be zero) —
and evaluates

    W_pred(u_shift, σ_vec) = softplus( joint(u_shift, σ_vec) ) / log 2

via :func:`flow_polyhead.predicted_W` (or its hand-coded C++
equivalent). One model call per event covers pure shift, pure
smear, and joint shift+smear evaluations for any number of
``(u_shift, σ_vec)`` pairs.

Perturbations are in target / y-space (the natural space for muon
calibration: a momentum-scale shift on reco pt directly affects
r_kappa, and a detector resolution smear is just per-event Gaussian
noise on the y-targets).

If you need ``z``, ``ladj``, or ``log p`` from the flow, export
the flow separately via ``flow_export_onnx.py``. This script
produces a polyhead-only inference package on purpose to minimize
per-event latency for the reweight use case.

Why batch=1 + AOTI by default: muon-calibration reweighting in narf
runs event-by-event in RDataFrame's per-row C++ helper. Static batch
plus ``cpp.weight_prepack=True`` + ``freezing=True`` + ``max_autotune``
gives the lowest-latency single-event path. Dynamic batch is also
supported (``--dynamic-batch``); the trade is that mkldnn weight-
prepack fusion does not fire, so per-event latency is ~2x worse.

The wrapper's forward is pure tensor ops (no functorch, no
``autograd.grad``), so ``torch.export.export`` traces it directly —
no ``make_fx`` fallback needed.
"""
import argparse
import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn
from torch._inductor import config as ind_cfg

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from flow_export_onnx import bake_masked_linears  # noqa: E402
from flow_polyhead import (  # noqa: E402
    FlowPolyheadInference,
    predicted_W,
)
from flow_training_diagnostics import (  # noqa: E402
    load_flow_from_checkpoint,
    load_polyhead_from_checkpoint,
)
from bake_couplings import bake_coupling_indices  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--checkpoint", required=True,
                   help="Path to checkpoint.pt or flow.pt with a "
                   "polyhead_state_dict entry.")
    p.add_argument("--output", required=True,
                   help="Path to .pt2 AOTI package output.")
    p.add_argument(
        "--batch", type=int, default=1,
        help="Static batch size baked into the package when "
        "--no-dynamic-batch is in effect. Default 1 — narf calls the "
        "model once per event, and static batch=1 is what makes "
        "Inductor's mkldnn weight-prepack fusion fire.",
    )
    p.add_argument(
        "--dynamic-batch", action=argparse.BooleanOptionalAction,
        default=False,
        help="Compile a single dynamic-batch package usable at any "
        "call size. Off by default (static batch wins on per-event "
        "latency at B=1; mkldnn weight-prepack does not fire on "
        "dynamic shapes).",
    )
    p.add_argument(
        "--max-autotune", action=argparse.BooleanOptionalAction,
        default=True,
        help="Inductor max-autotune for GEMMs. Slower compile, "
        "fastest runtime. On by default.",
    )
    p.add_argument(
        "--freeze", action=argparse.BooleanOptionalAction, default=True,
        help="Constant-fold weights so they can be prepacked / fused. "
        "Required for the oneDNN paths.",
    )
    p.add_argument(
        "--prepack", action=argparse.BooleanOptionalAction, default=True,
        help="Prepack Linear weights into oneDNN's preferred layout "
        "(cpp.weight_prepack=True). Only fires for static shapes.",
    )
    p.add_argument(
        "--mkldnn-fp32-bf16", action="store_true",
        help="Route fp32 matmul through mkldnn's BF16-accelerated "
        "path (set_float32_matmul_precision='medium'). On hardware "
        "without BF16 acceleration this is ignored or slower.",
    )
    p.add_argument(
        "--validate-n", type=int, default=64,
        help="Number of random events to validate eager vs AOTI on.",
    )
    p.add_argument(
        "--validate-tol", type=float, default=1e-4,
        help="Relative tolerance for validation.",
    )
    p.add_argument(
        "--bench", action="store_true",
        help="After export, microbenchmark the .pt2 forward at the "
        "static batch (pinned to one core).",
    )
    p.add_argument(
        "--bench-duration", type=float, default=2.0,
        help="Seconds to bench for.",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Build wrapper from checkpoint
# ---------------------------------------------------------------------------

def build_wrapper(checkpoint_path: str):
    """Load checkpoint, build flow + polyhead, return an
    :class:`FlowPolyheadInference` ready to export. Bakes masked
    linears and coupling indices so the FX graph has no
    ``aten.nonzero`` and no bool-mask-prim-device cast (those are
    what blocks Inductor freezing).
    """
    flow, stats, ckpt = load_flow_from_checkpoint(checkpoint_path, "cpu")
    n_features = int(len(stats.target_mean))
    n_cond = int(len(stats.cond_mean))
    polyhead = load_polyhead_from_checkpoint(
        ckpt, n_features=n_features, n_cond=n_cond, device="cpu",
    )
    if polyhead is None:
        raise RuntimeError(
            f"checkpoint {checkpoint_path!r} has no polyhead_state_dict; "
            "rerun training with --polyhead"
        )

    wrapper = FlowPolyheadInference(
        flow=flow,
        polyhead=polyhead,
        target_mean=torch.tensor(stats.target_mean, dtype=torch.float32),
        target_std=torch.tensor(stats.target_std, dtype=torch.float32),
        cond_mean=torch.tensor(stats.cond_mean, dtype=torch.float32),
        cond_std=torch.tensor(stats.cond_std, dtype=torch.float32),
    ).eval()

    n_baked_lin = bake_masked_linears(wrapper)
    n_baked_idx = bake_coupling_indices(wrapper)
    print(
        f"baked {n_baked_lin} MaskedLinear -> nn.Linear, "
        f"{n_baked_idx} coupling-mask index pair(s)"
    )
    return wrapper, stats


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def _decompose_linear(ep):
    """Replace ``aten.linear`` with explicit ``aten.mm`` / ``aten.addmm``
    so the mkldnn fusion can match — same trick as :file:`aoti_flow.py`.
    """
    aten = torch.ops.aten

    def _linear_decomp(input, weight, bias=None):
        if bias is None:
            return aten.mm.default(input, weight.permute(1, 0))
        return aten.addmm.default(bias, input, weight.permute(1, 0))

    decomp_table = {aten.linear.default: _linear_decomp}
    return ep.run_decompositions(decomp_table)


def _print_op_counts(ep, label):
    counts = {}
    for n in ep.graph_module.graph.nodes:
        if n.op == "call_function":
            t = str(n.target)
            counts[t] = counts.get(t, 0) + 1
    interesting = (
        "aten.linear.default", "aten.addmm.default", "aten.mm.default",
        "aten.nonzero.default",
    )
    print(f"  {label}:")
    for k in interesting:
        print(f"    {k}: {counts.get(k, 0)}")


def _to_tuple(out):
    """Normalize wrapper / runner output to a tuple — AOTI sometimes
    returns a single tensor, sometimes a 1-tuple.
    """
    if isinstance(out, (tuple, list)):
        return tuple(out)
    return (out,)


def _validate(wrapper, runner, n_features, n_cond, B, n_events, tol):
    """Compare eager vs AOTI on random inputs at the package's batch
    shape. ``n_events`` is rounded up to a multiple of ``B``.
    """
    rng = np.random.default_rng(0)
    n_events = max(B, ((n_events + B - 1) // B) * B)
    y = torch.from_numpy(
        rng.standard_normal((n_events, n_features)).astype(np.float32))
    c = torch.from_numpy(
        rng.standard_normal((n_events, n_cond)).astype(np.float32))

    eager_chunks = []
    aoti_chunks = []
    with torch.no_grad():
        for s in range(0, n_events, B):
            yi, ci = y[s:s+B], c[s:s+B]
            e = _to_tuple(wrapper(yi, ci))[0]
            a = _to_tuple(runner(yi, ci))[0]
            eager_chunks.append(e)
            aoti_chunks.append(a)
    et = torch.cat(eager_chunks)
    at = torch.cat(aoti_chunks)
    max_abs = (et - at).abs().max().item()
    denom = max(et.abs().max().item(), 1e-30)
    rel = max_abs / denom
    status = "OK" if rel < tol else "MISMATCH"
    print(
        f"[validate] joint_coefs: max |eager-aoti|={max_abs:.3e} "
        f"rel={rel:.2e}  {status}"
    )
    return rel < tol


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    if args.mkldnn_fp32_bf16:
        torch.set_float32_matmul_precision("medium")
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

    wrapper, stats = build_wrapper(args.checkpoint)
    n_features = int(len(stats.target_mean))
    n_cond = int(len(stats.cond_mean))

    B = args.batch
    # torch.export.export specializes a size-1 example dim as a
    # constant even when marked dynamic, so dynamic-batch tracing
    # needs an example with B >= 2.
    B_trace = max(B, 2) if args.dynamic_batch else B
    y = torch.randn(B_trace, n_features)
    c = torch.randn(B_trace, n_cond)

    # Eager smoke check on the (post-baking) wrapper.
    with torch.no_grad():
        ref = wrapper(y, c)
    print(f"eager: joint_coefs={tuple(ref.shape)}")

    if args.dynamic_batch:
        batch_dim = torch.export.Dim("batch", min=1, max=2**20)
        dynamic_shapes = {
            "y_raw": {0: batch_dim},
            "c_raw": {0: batch_dim},
        }
        print(
            f"\ntorch.export.export (dynamic batch, "
            f"B_example={B_trace}) ..."
        )
    else:
        dynamic_shapes = None
        print(f"\ntorch.export.export (static B={B}) ...")
    ep = torch.export.export(
        wrapper, (y, c), dynamic_shapes=dynamic_shapes,
    )
    print(" OK")

    ep = _decompose_linear(ep)
    _print_op_counts(ep, "post-decomp op counts")

    print(f"\naoti_compile_and_package -> {args.output}")
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with torch.no_grad():
        package_path = torch._inductor.aoti_compile_and_package(
            ep, package_path=args.output,
        )
    size_mb = os.path.getsize(package_path) / 1e6
    print(f" OK ({size_mb:.2f} MB)")

    runner = torch._inductor.aoti_load_package(package_path)
    out = _to_tuple(runner(y, c))[0]
    diff = (out - ref).abs().max().item()
    print(f"loaded; max |aoti - eager| joint_coefs = {diff:.2e}")

    if args.validate_n > 0:
        print()
        all_ok = _validate(
            wrapper, runner, n_features, n_cond, B,
            args.validate_n, args.validate_tol,
        )
        if not all_ok:
            print("[validate] WARNING: tolerance exceeded.")

    if args.bench:
        try:
            os.sched_setaffinity(0, {0})
        except Exception:
            pass
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        for _ in range(20):
            runner(y, c)
        n = 0
        t0 = time.perf_counter()
        deadline = t0 + args.bench_duration
        while True:
            runner(y, c)
            n += 1
            if n >= 5 and time.perf_counter() >= deadline:
                break
        elapsed = time.perf_counter() - t0
        ms_b = 1000.0 * elapsed / n
        us_e = 1e6 * elapsed / n / B
        ev_s = n * B / elapsed
        print(
            f"\nbench (B={B}, 1 core): "
            f"{ms_b:.3f} ms/call  {us_e:.3f} us/event  "
            f"{ev_s:.1f} ev/s  ({n} iters)"
        )


if __name__ == "__main__":
    main()
