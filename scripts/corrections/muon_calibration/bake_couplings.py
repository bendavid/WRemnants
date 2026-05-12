"""Bake zuko's masked-coupling-layer indices as constants.

zuko's ``GeneralCouplingTransform.forward(c)`` constructs a fresh
``CouplingTransform`` on every call, whose ``__init__`` does

    self.idx_a = mask.nonzero().squeeze(-1)
    self.idx_b = (~mask).nonzero().squeeze(-1)

so each forward pass through the flow runs ``aten.nonzero()`` per
coupling layer. That is what blocks Inductor freezing (and ExecuTorch's
view-spec validator) — ``nonzero()`` returns a data-dependent shape,
which the constant-folder cannot track through.

This module:

1. Defines :class:`FastCouplingTransform`, a subclass that takes
   pre-computed ``idx_a`` / ``idx_b`` directly (no nonzero call).
2. Provides :func:`bake_coupling_indices` which walks a flow, finds
   every module with a ``mask`` buffer (the GeneralCouplingTransform-
   like modules), pre-computes ``idx_a``/``idx_b`` once and stores them
   as buffers, then patches the module's ``forward`` to instantiate
   :class:`FastCouplingTransform` with those buffers.

After baking, the exported FX graph contains only the gather/scatter
indexing using constant index tensors — no ``aten.nonzero`` and no
data-dependent shapes. Same numerics as the original.
"""
from __future__ import annotations

from functools import partial

import torch
from zuko.transforms import CouplingTransform


class FastCouplingTransform(CouplingTransform):
    """``CouplingTransform`` with pre-computed split indices.

    Unlike the parent, this does not call ``mask.nonzero()`` in
    ``__init__``; the indices are passed in directly. Everything else
    (split/merge/_call/log_abs_det_jacobian/call_and_ladj) is inherited
    unchanged.
    """

    def __init__(self, meta, idx_a, idx_b, **kwargs):
        # Skip CouplingTransform.__init__ entirely; jump to the
        # ancestor Transform's __init__ via super().__init__ chain.
        # CouplingTransform.__mro__ is
        # (CouplingTransform, Transform, object).
        from torch.distributions.transforms import Transform
        Transform.__init__(self, **kwargs)
        self.meta = meta
        self.idx_a = idx_a
        self.idx_b = idx_b


def _make_baked_forward(module, idx_a, idx_b):
    """Build a replacement ``forward(c=None)`` that returns a
    :class:`FastCouplingTransform` with the pre-computed indices.

    Mirrors ``GeneralCouplingTransform.forward`` in
    ``zuko/flows/coupling.py``::

        return CouplingTransform(partial(self.meta, c), self.mask)
    """
    def baked_forward(c=None):
        meta_fn = partial(module.meta, c)
        return FastCouplingTransform(meta_fn, idx_a, idx_b)
    return baked_forward


def bake_coupling_indices(model: torch.nn.Module) -> int:
    """Walk ``model``, find masked-coupling modules, bake their indices.

    A module is treated as a coupling layer if it has all of:
      * a ``mask`` buffer of dtype bool/int and ``ndim == 1``
      * a ``meta(c, x)`` method (the conditioner-builder)
      * a ``forward(c=None)`` returning a Transform

    Returns the number of modules patched.
    """
    n = 0
    for module in model.modules():
        mask = getattr(module, "mask", None)
        if not isinstance(mask, torch.Tensor):
            continue
        if mask.ndim != 1:
            continue
        if not hasattr(module, "meta") or not callable(module.meta):
            continue
        # Skip if already baked.
        if getattr(module, "_coupling_indices_baked", False):
            continue

        idx_a = mask.nonzero().squeeze(-1).contiguous()
        idx_b = (~mask).nonzero().squeeze(-1).contiguous()

        module.register_buffer("_baked_idx_a", idx_a, persistent=False)
        module.register_buffer("_baked_idx_b", idx_b, persistent=False)
        module._coupling_indices_baked = True

        module.forward = _make_baked_forward(
            module, module._baked_idx_a, module._baked_idx_b,
        )
        n += 1
    return n
