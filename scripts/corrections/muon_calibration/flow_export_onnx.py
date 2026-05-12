"""Export a trained muon-response flow checkpoint to ONNX.

Loads ``checkpoint.pt`` (or ``flow.pt``) produced by
``train_muon_response_flow.py`` and exports the flow as an ONNX model
suitable for C++/Python inference via ``onnxruntime``.

The exported model accepts raw-space inputs
``(x_raw, c_raw)`` — i.e. the same per-event targets
``(r_kappa, dlambda, dphi)`` and raw conditioning ``(log_pt_gen,
charge, lambda_gen, sin_phi_gen, cos_phi_gen)`` you would use in the
diagnostic script — and returns raw-space log-density ``log_p`` and
optionally the raw-space score ``∇_x log p``. Preprocessing
(standardization) and the Jacobian correction are baked into the
exported graph, so the caller doesn't need the ``preproc.json``.

Precision. The default export is FP32, which runs anywhere.
``--dtype bfloat16`` casts weights and exposes bf16 I/O. On CPUs with
hardware BF16 acceleration (Zen 4+ EPYC, Sapphire Rapids+ Intel, ARM
NEON-BF16) this gives ~1.3-1.8x inference speedup. On older CPUs the
bf16 path is typically the same or slightly slower than FP32, so
only turn it on when the deployment target is known to be BF16-
accelerated.

Score export. ``--include-score`` adds a second output that is the
gradient of ``log_p`` wrt ``x_raw``, computed analytically via
``torch.autograd.grad`` at export time. Requires the dynamo ONNX
backend (``--dynamo``).

Base-space derivative export. ``--include-base={lin,full,both}`` adds
further outputs that give the caller everything needed to evaluate the
base-space reweight formulas (see ``plot_mc_{shift,smear}_reweight`` in
``flow_training_diagnostics.py``), working in standardized coordinates
where the flow's base is ``N(0, I)``:

  * ``z`` (``[B, d]``): latent, ``z = transform(x_std; c_std)``.
  * ``J`` (``[B, d, d]``): inverse Jacobian ``J[n, j, i] = ∂z_j/∂x_std_i``.
  * ``s_std`` (``[B, d]``, full/both): standardized score
    ``∂ log p_std / ∂ x_std``.
  * ``H_std`` (``[B, d]``, full/both): diag of standardized Hessian
    ``∂² log p_std / ∂ x_std_i²``.
  * ``G_std`` (``[B, d]``, lin/both): ``∂ L / ∂ x_std`` with
    ``L = log|det ∂z/∂x_std|``.
  * ``K_std`` (``[B, d]``, lin/both): ``∂² L / ∂ x_std_i²``.

With ``full`` the downstream computes ``G_i = s_std_i + zᵀv_i``
(``v_i = J[:, :, i]``) and assembles the FULL weight as
``exp(δ·zᵀv − ½δ²‖v‖²) · (1 − δG + ½δ²·(G² + H + ‖v‖²))``. With
``lin`` the downstream uses ``(G_std, K_std)`` directly in the LIN form
``(1 − δG + ½δ²·(G² + K))``. ``both`` emits the union so either form
can be chosen at runtime. Requires ``--dynamo`` (higher-order
``torch.autograd.grad`` inside the traced graph).

Validation. ``--validate`` (default on) runs the exported ONNX model
through ``onnxruntime`` on random inputs and compares to the
PyTorch forward. Skipped silently if ``onnxruntime`` isn't importable.
"""

import argparse
import math
import os
import sys
from typing import Optional

import numpy as np
import torch
import torch.nn as nn

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Importing flow_training_diagnostics triggers the zuko patches that
# make gauss_legendre twice-differentiable AND polynomial basis safe at
# x=0. We only need first-derivative autograd for --include-score (the
# score), but the patches are harmless for export.
import flow_training_diagnostics  # noqa: F401
from flow_training_diagnostics import (  # noqa: E402
    load_flow_from_checkpoint,
)


def bake_masked_linears(module: nn.Module) -> int:
    """Replace each ``zuko.nn.MaskedLinear`` with a plain ``nn.Linear``
    whose weight has been pre-multiplied by the mask.

    zuko's ``MaskedLinear.forward`` does ``F.linear(x, self.mask *
    self.weight, self.bias)``. The BoolTensor-buffer-times-float-
    parameter implicit cast produces a ``prim.device`` op that
    torch.onnx.export (both the dynamo-based and legacy paths) fails
    to decompose. Pre-multiplying the weights offline turns this into
    an ordinary ``nn.Linear`` call that exports cleanly.

    In-place replacement, returns the number of layers replaced.
    Safe only for a frozen/eval model: gradient structure w.r.t. the
    original unmasked weight is lost.
    """
    from zuko.nn import MaskedLinear

    replaced = 0
    for name, submodule in list(module.named_children()):
        if isinstance(submodule, MaskedLinear):
            out_f, in_f = submodule.weight.shape
            use_bias = submodule.bias is not None
            linear = nn.Linear(in_f, out_f, bias=use_bias)
            with torch.no_grad():
                masked_w = (
                    submodule.weight
                    * submodule.mask.to(submodule.weight.dtype)
                )
                linear.weight.copy_(masked_w)
                if use_bias:
                    linear.bias.copy_(submodule.bias)
            # Match dtype of the source.
            linear = linear.to(submodule.weight.dtype)
            setattr(module, name, linear)
            replaced += 1
        else:
            replaced += bake_masked_linears(submodule)
    return replaced


# -----------------------------------------------------------------------------
# ONNX-friendly wrapper
# -----------------------------------------------------------------------------

_BASE_MODES = ("none", "lin", "full", "both")


def _output_names_for(include_score: bool, include_base: str) -> list:
    """Canonical order of output names for (include_score, include_base)."""
    names = ["log_p"]
    if include_score:
        names.append("score")
    if include_base != "none":
        names.extend(["z", "J"])
    if include_base in ("full", "both"):
        names.extend(["s_std", "H_std"])
    if include_base in ("lin", "both"):
        names.extend(["G_std", "K_std"])
    return names


class FlowONNXWrapper(nn.Module):
    """Wrap a zuko flow so its ``forward`` returns plain tensors.

    Minimal mode: ``forward(x_raw, c_raw) -> log_p`` (raw-space).
    With ``include_score`` adds the raw-space score
    ``∂ log_p / ∂ x_raw``. With ``include_base ∈ {lin, full, both}``
    adds the standardized-space derivatives needed for base-space
    reweighting — see module docstring for the per-mode output list.

    All standardization constants are registered as buffers so they
    travel with the exported graph; the caller works purely in raw
    units for ``(x, c)`` and ``log_p``.
    """

    def __init__(
        self, flow, stats,
        include_score: bool = False,
        include_base: str = "none",
    ):
        super().__init__()
        if include_base not in _BASE_MODES:
            raise ValueError(
                f"include_base must be one of {_BASE_MODES}, got "
                f"{include_base!r}"
            )
        self.flow = flow
        self.include_score = include_score
        self.include_base = include_base

        target_mean = torch.tensor(stats.target_mean, dtype=torch.float32)
        target_std = torch.tensor(stats.target_std, dtype=torch.float32)
        cond_mean = torch.tensor(stats.cond_mean, dtype=torch.float32)
        cond_std = torch.tensor(stats.cond_std, dtype=torch.float32)

        self.register_buffer("target_mean", target_mean)
        self.register_buffer("target_std", target_std)
        self.register_buffer("cond_mean", cond_mean)
        self.register_buffer("cond_std", cond_std)
        # Log-Jacobian of the target standardization (additive constant).
        self.register_buffer(
            "log_target_std_sum",
            torch.log(target_std).sum(),
        )

    def _standardize(self, x_raw, c_raw):
        x = (x_raw - self.target_mean) / self.target_std
        c = (c_raw - self.cond_mean) / self.cond_std
        return x, c

    def log_density(self, x_raw, c_raw):
        x, c = self._standardize(x_raw, c_raw)
        log_p_std = self.flow(c).log_prob(x)
        return log_p_std - self.log_target_std_sum

    def output_names(self) -> list:
        return _output_names_for(self.include_score, self.include_base)

    def _derivatives_forward(self, x_raw, c_raw):
        """Forward path that computes log_p plus every requested derivative.

        Uses ``flow.transform.call_and_ladj`` so ``z`` and the
        log-det-Jacobian ``L = log|det ∂z/∂x_std|`` are available
        separately, then autograd on ``x_std`` gives J, s_std, H_std,
        G_std, K_std. Derivatives are in standardized coordinates (the
        flow's natural space with base ``N(0, I)``); ``log_p`` and
        ``score`` are kept in raw space to match the existing contract.
        """
        # Detach and standardize; make x_std the leaf that autograd
        # differentiates against.
        x_std = ((x_raw.detach() - self.target_mean) / self.target_std)
        x_std = x_std.requires_grad_(True)
        c_std = (c_raw.detach() - self.cond_mean) / self.cond_std

        # Flow forward: z and L = log|det ∂z/∂x_std|.
        z, ladj = self.flow(c_std).transform.call_and_ladj(x_std)
        d = z.shape[-1]

        # log p in standardized space = log φ(z) + L, with φ = N(0, I).
        log_phi = -0.5 * (z * z).sum(dim=-1) - 0.5 * float(d) * math.log(
            2.0 * math.pi
        )
        log_p_std = log_phi + ladj
        log_p_raw = log_p_std - self.log_target_std_sum

        outputs = {"log_p": log_p_raw}

        # Need Jacobian J and z whenever include_base is on.
        need_base = self.include_base != "none"
        need_full = self.include_base in ("full", "both")
        need_lin = self.include_base in ("lin", "both")
        need_score = self.include_score

        # ---- J = ∂z/∂x_std via d vector-jacobian products (rows). ----
        if need_base:
            j_rows = []
            for j in range(d):
                (row_j,) = torch.autograd.grad(
                    z[..., j].sum(), x_std,
                    create_graph=False,
                    retain_graph=True,
                )
                j_rows.append(row_j)
            # Stack as J[n, j, i] = ∂z_j/∂x_std_i.
            outputs["z"] = z
            outputs["J"] = torch.stack(j_rows, dim=-2)

        # ---- Score and diag Hessian of log_p in std-space. ----
        if need_score or need_full:
            (grad_lp,) = torch.autograd.grad(
                log_p_std.sum(), x_std,
                create_graph=need_full,
                retain_graph=True,
            )
            if need_score:
                # Chain rule: ∂log_p_raw/∂x_raw = grad_lp / target_std.
                outputs["score"] = grad_lp.detach() / self.target_std
            if need_full:
                outputs["s_std"] = grad_lp.detach()
                h_cols = []
                for i in range(d):
                    retain = (i < d - 1) or need_lin
                    (row_i,) = torch.autograd.grad(
                        grad_lp[..., i].sum(), x_std,
                        create_graph=False,
                        retain_graph=retain,
                    )
                    h_cols.append(row_i[..., i])
                outputs["H_std"] = torch.stack(h_cols, dim=-1)

        # ---- Grad and diag Hessian of L = log|det ∂z/∂x_std|. ----
        if need_lin:
            (grad_ladj,) = torch.autograd.grad(
                ladj.sum(), x_std,
                create_graph=True,
                retain_graph=True,
            )
            outputs["G_std"] = grad_ladj.detach()
            k_cols = []
            for i in range(d):
                retain = i < d - 1
                (row_i,) = torch.autograd.grad(
                    grad_ladj[..., i].sum(), x_std,
                    create_graph=False,
                    retain_graph=retain,
                )
                k_cols.append(row_i[..., i])
            outputs["K_std"] = torch.stack(k_cols, dim=-1)

        return outputs

    def forward(self, x_raw, c_raw):
        if not self.include_score and self.include_base == "none":
            return self.log_density(x_raw, c_raw)

        outputs = self._derivatives_forward(x_raw, c_raw)
        names = self.output_names()
        return tuple(outputs[n] for n in names)


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
        help="Path to checkpoint.pt / flow.pt from train_muon_response_flow.py.",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output ONNX file path.",
    )
    p.add_argument(
        "--dtypes",
        nargs="+",
        choices=["float32", "bfloat16"],
        default=["float32", "bfloat16"],
        help="One or both of float32/bfloat16. Default: both. "
        "Writes <output> for float32 and <output_stem>_bf16.onnx "
        "for bfloat16. The bfloat16 file needs ORT CPU BF16 kernels "
        "(Zen 4+/Sapphire Rapids+, ORT 1.17+ with oneDNN, etc.) — "
        "the FP32 file is the portable fallback. The suggested "
        "caller pattern: try BF16 first, catch NotImplemented, fall "
        "back to FP32.",
    )
    p.add_argument(
        "--include-score",
        action="store_true",
        help="Also export the score (grad of log_p wrt x_raw). "
        "Requires --dynamo.",
    )
    p.add_argument(
        "--include-base",
        choices=_BASE_MODES,
        default="none",
        help="Also export standardized-space derivatives for base-"
        "space reweight. 'lin' → adds (z, J, G_std, K_std); 'full' "
        "→ adds (z, J, s_std, H_std); 'both' → union. Requires "
        "--dynamo. Default: none.",
    )
    p.add_argument(
        "--dynamo",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Use torch.onnx.export(dynamo=True) — the TorchDynamo-based "
        "exporter. Default: on when --include-score or --include-base "
        "is set (needed for the autograd.grad calls), off otherwise.",
    )
    p.add_argument(
        "--opset",
        type=int,
        default=17,
        help="ONNX opset version.",
    )
    p.add_argument(
        "--batch-example",
        type=int,
        default=1,
        help="Batch size used for the tracing example input. With "
        "dynamic axes enabled, runtime batch size is unconstrained.",
    )
    p.add_argument(
        "--dynamic-batch",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to mark the batch dimension as dynamic so the "
        "exported model accepts variable batch sizes at runtime.",
    )
    p.add_argument(
        "--validate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run the exported ONNX through onnxruntime on random "
        "inputs and compare to PyTorch. Skipped if onnxruntime isn't "
        "importable.",
    )
    p.add_argument(
        "--validate-n",
        type=int,
        default=256,
        help="Number of random events used for validation.",
    )
    p.add_argument(
        "--validate-tol",
        type=float,
        default=1e-4,
        help="Relative tolerance for validation (per-event max abs "
        "diff divided by PyTorch result's abs max). Warns if exceeded.",
    )
    return p.parse_args()


# -----------------------------------------------------------------------------
# Validation: run ORT vs PyTorch and report discrepancies
# -----------------------------------------------------------------------------

def _to_np(tensor: torch.Tensor) -> np.ndarray:
    """Convert a tensor to numpy, downgrading any tensor subclass."""
    t = tensor.detach().cpu()
    if type(t) is not torch.Tensor:
        print(f"[_to_np] downgrading {type(t).__name__}")
        t = t.as_subclass(torch.Tensor)
    if t.dtype == torch.bfloat16:
        t = t.float()
    try:
        return t.numpy()
    except RuntimeError:
        # Last resort: copy through Python list.
        return np.asarray(t.tolist())


def _compute_torch_reference(
    wrapper: nn.Module,
    n_features: int,
    n_cond: int,
    n_events: int,
    needs_grad: bool,
    dtype: torch.dtype,
):
    """Precompute PyTorch reference outputs + inputs BEFORE export.

    torch.onnx.export / torch.export leave FakeTensor tracing state
    behind that contaminates later tensor->numpy conversions; do the
    reference run up front and stash the plain numpy arrays.
    """
    rng = np.random.default_rng(0)
    x_np = rng.standard_normal((n_events, n_features)).astype(np.float32)
    c_np = rng.standard_normal((n_events, n_cond)).astype(np.float32)

    x_t = torch.from_numpy(x_np).to(dtype)
    c_t = torch.from_numpy(c_np).to(dtype)

    if needs_grad:
        out_torch = wrapper(x_t, c_t)
        if not isinstance(out_torch, (tuple, list)):
            out_torch = (out_torch,)
        outputs_np = [_to_np(o) for o in out_torch]
    else:
        with torch.no_grad():
            out_torch = wrapper.log_density(x_t, c_t)
        outputs_np = [_to_np(out_torch)]

    return x_np, c_np, outputs_np


def _validate_with_onnxruntime(
    onnx_path: str,
    torch_inputs: tuple,
    torch_outputs: list,
    tol: float,
) -> bool:
    """Compare pre-computed PyTorch outputs against ORT on the same
    inputs. No torch forward pass here — avoids post-export subclass
    contamination.
    """
    try:
        import onnxruntime as ort
    except ImportError:
        print("[validate] onnxruntime not available; skipping.")
        return True

    x_np, c_np = torch_inputs
    sess = ort.InferenceSession(
        onnx_path, providers=["CPUExecutionProvider"],
    )
    input_names = [i.name for i in sess.get_inputs()]
    output_names = [o.name for o in sess.get_outputs()]
    print(f"[validate] ORT inputs={input_names} outputs={output_names}")

    feed = dict(zip(input_names, [x_np, c_np]))
    outputs_ort = sess.run(None, feed)

    ok = True
    for name, ort_val, torch_val in zip(
        output_names, outputs_ort, torch_outputs
    ):
        abs_max_torch = np.abs(torch_val).max() + 1e-30
        max_abs_diff = np.abs(ort_val - torch_val).max()
        rel = max_abs_diff / abs_max_torch
        status = "OK" if rel < tol else "MISMATCH"
        print(
            f"[validate] {name}: max |ort - pytorch| = {max_abs_diff:.3e} "
            f"(rel {rel:.2e})  {status}"
        )
        if rel >= tol:
            ok = False
    return ok


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def _output_path_for_dtype(base_path: str, dtype_name: str) -> str:
    """FP32 uses ``base_path`` verbatim; BF16 inserts ``_bf16`` before
    the extension so both files can coexist in one directory.
    """
    if dtype_name == "float32":
        return base_path
    stem, ext = os.path.splitext(base_path)
    if not ext:
        ext = ".onnx"
    return f"{stem}_bf16{ext}"


def _build_wrapper(
    checkpoint_path: str,
    include_score: bool,
    include_base: str,
    dtype: torch.dtype,
):
    """Load a fresh wrapper at the requested dtype. Fresh per-dtype so
    we don't share state (bake_masked_linears + dtype cast) across
    iterations."""
    flow, stats, _ = load_flow_from_checkpoint(checkpoint_path, "cpu")
    wrapper = FlowONNXWrapper(
        flow=flow, stats=stats,
        include_score=include_score, include_base=include_base,
    ).eval()
    n_baked = bake_masked_linears(wrapper)
    if dtype != torch.float32:
        wrapper = wrapper.to(dtype)
    return wrapper, stats, n_baked


def _export_one(
    *,
    checkpoint_path: str,
    dtype_name: str,
    output_path: str,
    include_score: bool,
    include_base: str,
    use_dynamo: bool,
    opset: int,
    dynamic_batch: bool,
    batch_example: int,
    validate: bool,
    validate_n: int,
    validate_tol: float,
) -> bool:
    """Export one dtype to one ONNX file, with its own fresh wrapper.
    Returns True on success (validation passing or skipped).
    """
    dtype = torch.bfloat16 if dtype_name == "bfloat16" else torch.float32
    print(
        f"\n=== exporting {dtype_name} → {output_path} "
        f"(opset {opset}, dynamo {use_dynamo}, include_base={include_base}) ==="
    )

    wrapper, stats, n_baked = _build_wrapper(
        checkpoint_path, include_score, include_base, dtype,
    )
    if n_baked:
        print(f"baked {n_baked} MaskedLinear layer(s) into nn.Linear")

    n_features = int(len(stats.target_mean))
    n_cond = int(len(stats.cond_mean))

    x_example = torch.randn(batch_example, n_features, dtype=dtype)
    c_example = torch.randn(batch_example, n_cond, dtype=dtype)

    needs_grad = include_score or include_base != "none"
    ref_inputs = ref_outputs = None
    if validate:
        ref_x, ref_c, ref_outputs = _compute_torch_reference(
            wrapper,
            n_features=n_features,
            n_cond=n_cond,
            n_events=validate_n,
            needs_grad=needs_grad,
            dtype=dtype,
        )
        ref_inputs = (ref_x, ref_c)

    output_names = _output_names_for(include_score, include_base)

    dynamic_axes = None
    if dynamic_batch and not use_dynamo:
        dynamic_axes = {
            "x_raw": {0: "batch"},
            "c_raw": {0: "batch"},
        }
        for name in output_names:
            dynamic_axes[name] = {0: "batch"}

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if use_dynamo:
        torch.onnx.export(
            wrapper,
            (x_example, c_example),
            output_path,
            input_names=["x_raw", "c_raw"],
            output_names=output_names,
            opset_version=opset,
            dynamo=True,
        )
    else:
        torch.onnx.export(
            wrapper,
            (x_example, c_example),
            output_path,
            input_names=["x_raw", "c_raw"],
            output_names=output_names,
            opset_version=opset,
            dynamic_axes=dynamic_axes,
        )
    print(
        f"saved {output_path}  "
        f"({os.path.getsize(output_path) / 1e6:.1f} MB)"
    )

    # Reset dynamo state left behind by the export pass.
    try:
        import torch._dynamo as _tdyn
        _tdyn.reset()
    except Exception:
        pass

    if validate:
        try:
            ok = _validate_with_onnxruntime(
                output_path,
                torch_inputs=ref_inputs,
                torch_outputs=ref_outputs,
                tol=validate_tol,
            )
        except Exception as e:
            print(
                f"[validate] skipped ({type(e).__name__}: {e}). "
                f"This is expected for BF16 on CPUs without hardware "
                f"BF16 kernels (Pre-Zen4, pre-Sapphire-Rapids)."
            )
            ok = True
        if not ok:
            print("[validate] WARNING: max relative error exceeds tolerance.")
            return False
    return True


def main():
    args = parse_args()

    needs_autograd = args.include_score or args.include_base != "none"
    if needs_autograd and args.dynamo is False:
        raise SystemExit(
            "--include-score / --include-base require the dynamo "
            "exporter (pass --dynamo or drop --no-dynamo)."
        )
    use_dynamo = args.dynamo if args.dynamo is not None else needs_autograd

    device = "cpu"
    print(f"loading checkpoint {args.checkpoint}")
    _, _, ckpt = load_flow_from_checkpoint(args.checkpoint, device)
    if "epoch" in ckpt:
        print(
            f"  epoch {ckpt['epoch']}  "
            f"val_nll {ckpt.get('val_nll', float('nan')):+.4f}"
        )

    # Deduplicate while preserving order; write FP32 first if requested
    # so the base --output path gets the portable model.
    dtype_order = sorted(
        set(args.dtypes),
        key=lambda d: 0 if d == "float32" else 1,
    )

    all_ok = True
    for dtype_name in dtype_order:
        out_path = _output_path_for_dtype(args.output, dtype_name)
        ok = _export_one(
            checkpoint_path=args.checkpoint,
            dtype_name=dtype_name,
            output_path=out_path,
            include_score=args.include_score,
            include_base=args.include_base,
            use_dynamo=use_dynamo,
            opset=args.opset,
            dynamic_batch=args.dynamic_batch,
            batch_example=args.batch_example,
            validate=args.validate,
            validate_n=args.validate_n,
            validate_tol=(
                args.validate_tol if dtype_name == "float32"
                else max(args.validate_tol, 1e-2)
            ),
        )
        all_ok = all_ok and ok

    print("\ndone.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
