"""Inference-side wrapper + re-exports for the polynomial-correction head.

The polyhead module class (:class:`PolyHead`), the joint training step
(:func:`joint_loss_step`), the polynomial-evaluation helpers
(:func:`joint_monomial_basis`, :func:`evaluate_joint`,
:func:`predicted_W`, :func:`predicted_log_W`), and the perturbation
sampler all live in :mod:`train_muon_response_flow` so the training
script is fully self-contained when copied around without sibling
files.

This module re-exports them under their canonical name and adds
:class:`FlowPolyheadInference`, the AOTI-export-time wrapper that
bakes standardization into buffers and exposes a pure-tensor-op
``forward(y_raw, c_raw) → (z, joint_coefs)``. The downstream weight
evaluation uses :func:`predicted_W` (re-exported here) on the
package's outputs:

    W_pred(u_shift, σ_vec) = softplus( joint(u_shift, σ_vec) ) / log 2

— a pure function of the head's joint-polynomial coefs. The flow's
``z`` is exposed alongside for downstream uses (e.g., ``log p`` via
``log_phi(z) + ladj`` if needed) but is not consumed by ``predicted_W``
itself.
"""
from __future__ import annotations

import os
import sys

import torch
import torch.nn as nn

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from train_muon_response_flow import (  # noqa: E402
    PolyHead,
    evaluate_joint,
    joint_loss_step,
    joint_monomial_basis,
    predicted_W,
    predicted_log_W,
    sample_perturbations,
)


__all__ = [
    "PolyHead",
    "FlowPolyheadInference",
    "evaluate_joint",
    "joint_loss_step",
    "joint_monomial_basis",
    "predicted_W",
    "predicted_log_W",
    "sample_perturbations",
]


class FlowPolyheadInference(nn.Module):
    """``forward(y_raw, c_raw) → joint_coefs``.

    Inference-time wrapper that runs **only the polyhead**, with
    standardization baked into buffers. The flow is not part of the
    inference graph — the polyhead's prediction has no flow-state
    dependence (it's a pure function of ``(y_std, c_std)``), so the
    flow forward would be wasted compute.

    Pure tensor ops, ``torch.export.export`` traces directly,
    dynamic-batch supported. The caller forms y-space perturbations
    in standardized target units::

        u_shift_std = u_shift_raw / target_std
        σ_vec_std   = σ_vec_raw / target_std

    outside the package and calls :func:`predicted_W` on
    ``joint_coefs`` for each ``(u_shift, σ_vec)`` pair — pure shift,
    pure smear, or joint shift+smear all use the same formula. The
    polynomial inputs live in target / y-space (R^{n_features}), so
    a "1% pt-scale shift" maps directly to a deterministic y-space
    shift of the corresponding ``r_kappa`` component.

    If you also need the flow's ``z`` or ``log p`` (for diagnostics
    or a separate density-evaluation use case), export the flow
    separately via ``flow_export_onnx.py``.

    The ``flow`` argument is accepted (and stored) for API
    compatibility but is not used in the forward — keeping the
    constructor signature stable means downstream tooling that
    builds the wrapper from a checkpoint doesn't have to branch.
    """

    def __init__(
        self,
        flow: nn.Module,
        polyhead: PolyHead,
        target_mean: torch.Tensor,
        target_std: torch.Tensor,
        cond_mean: torch.Tensor,
        cond_std: torch.Tensor,
    ):
        super().__init__()
        # Flow is stored for API compat but not used in the inference
        # forward — this keeps the AOTI graph polyhead-only.
        self.flow = flow
        self.polyhead = polyhead
        self.register_buffer("target_mean", target_mean)
        self.register_buffer("target_std", target_std)
        self.register_buffer("cond_mean", cond_mean)
        self.register_buffer("cond_std", cond_std)

    def forward(
        self,
        y_raw: torch.Tensor,
        c_raw: torch.Tensor,
    ) -> torch.Tensor:
        y_std = (y_raw - self.target_mean) / self.target_std
        c_std = (c_raw - self.cond_mean) / self.cond_std
        return self.polyhead(y_std, c_std)
