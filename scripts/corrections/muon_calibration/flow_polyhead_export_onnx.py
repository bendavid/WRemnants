"""ONNX export of a (flow + polyhead) checkpoint — polyhead only.

Loads a checkpoint produced by ``train_muon_response_flow.py
--polyhead`` and emits an ONNX model with signature

    forward(y_raw, c_raw) -> joint_coefs

i.e., the same wrapper :class:`FlowPolyheadInference` used by the AOTI
exporter (:file:`flow_polyhead_export.py`), but serialized to ONNX
instead of an AOT-Inductor ``.pt2`` package.

For the C++ ONNX-Runtime benchmark and for any downstream consumer
that prefers ORT over libtorch+AOTI. The flow itself is *not* part of
the exported graph — the polyhead's prediction has no flow-state
dependence (pure function of ``(y_std, c_std)``).

The caller forms y-space perturbations in standardized target units
outside the model and evaluates ``W_pred`` via
:func:`flow_polyhead.predicted_W` on the emitted ``joint_coefs``,
exactly as in the AOTI path.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import torch

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from flow_polyhead_export import build_wrapper  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--checkpoint", required=True,
        help="[required] Path to checkpoint.pt or flow.pt with a "
        "polyhead_state_dict entry.",
    )
    p.add_argument(
        "--output", required=True,
        help="[required] Output ONNX file path.",
    )
    p.add_argument(
        "--batch-example", type=int, default=1,
        help="[default: %(default)s] Batch size used for the tracing "
        "example input. With --dynamic-batch the runtime batch size "
        "is unconstrained.",
    )
    p.add_argument(
        "--dynamic-batch",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="[default: on] Mark the batch dimension as dynamic so "
        "the exported model accepts variable batch sizes at runtime.",
    )
    p.add_argument(
        "--opset", type=int, default=17,
        help="[default: %(default)s] ONNX opset version.",
    )
    p.add_argument(
        "--validate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="[default: on] Run the exported ONNX through onnxruntime "
        "and compare to PyTorch on random inputs. Skipped if "
        "onnxruntime isn't importable.",
    )
    p.add_argument(
        "--validate-n", type=int, default=256,
        help="[default: %(default)s] Number of random events for "
        "validation.",
    )
    p.add_argument(
        "--validate-tol", type=float, default=1e-4,
        help="[default: %(default)s] Relative tolerance for "
        "validation.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    wrapper, stats = build_wrapper(args.checkpoint)
    n_features = int(len(stats.target_mean))
    n_cond = int(len(stats.cond_mean))

    B = args.batch_example
    y = torch.randn(B, n_features)
    c = torch.randn(B, n_cond)

    with torch.no_grad():
        ref = wrapper(y, c)
    print(f"eager: joint_coefs={tuple(ref.shape)}")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    dynamic_axes = None
    if args.dynamic_batch:
        dynamic_axes = {
            "y_raw": {0: "batch"},
            "c_raw": {0: "batch"},
            "joint_coefs": {0: "batch"},
        }

    print(
        f"\ntorch.onnx.export ("
        f"{'dynamic' if args.dynamic_batch else 'static'} batch, "
        f"B_example={B}, opset={args.opset}) -> {args.output}"
    )
    torch.onnx.export(
        wrapper,
        (y, c),
        args.output,
        input_names=["y_raw", "c_raw"],
        output_names=["joint_coefs"],
        dynamic_axes=dynamic_axes,
        opset_version=args.opset,
        do_constant_folding=True,
        export_params=True,
    )
    size_mb = os.path.getsize(args.output) / 1e6
    print(f" OK ({size_mb:.2f} MB)")

    if args.validate:
        try:
            import onnxruntime as ort
        except ImportError:
            print("[validate] onnxruntime not installed — skipping.")
            return

        sess = ort.InferenceSession(
            args.output, providers=["CPUExecutionProvider"],
        )
        rng = np.random.default_rng(0)
        N = max(args.validate_n, B)
        y_v = rng.standard_normal((N, n_features)).astype(np.float32)
        c_v = rng.standard_normal((N, n_cond)).astype(np.float32)
        with torch.no_grad():
            ref = wrapper(
                torch.from_numpy(y_v), torch.from_numpy(c_v),
            ).numpy()
        ort_out = sess.run(
            ["joint_coefs"],
            {"y_raw": y_v, "c_raw": c_v},
        )[0]
        max_abs = float(np.abs(ref - ort_out).max())
        denom = max(float(np.abs(ref).max()), 1e-30)
        rel = max_abs / denom
        status = "OK" if rel < args.validate_tol else "MISMATCH"
        print(
            f"[validate] joint_coefs over N={N}: "
            f"max |eager-ort|={max_abs:.3e} rel={rel:.2e}  {status}"
        )


if __name__ == "__main__":
    main()
