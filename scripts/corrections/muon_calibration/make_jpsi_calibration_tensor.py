"""Build a rabbit tensor for the J/psi muon-calibration fit.

Inputs (from /home/b/bendavid/cvhdev4/work by default):

  jpsi_module_corrections_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr.hdf5
      Output of scripts/histmakers/jpsi_module_corrections.py. Contains:
        - jpsi_ideal/output/hmuplus, hmuminus     (4D dense templates)
        - jpsi_ideal/output/<uuid> (x2)            (5D sparse templates with
                                                    extra "corparms" axis,
                                                    delta wrt nominal)
        - jpsi_nom/output/hmuplus, hmuminus       (treated as observed data)

  combinedgrads.hdf5
      Output of scripts/corrections/muon_calibration/aggregategrads.py.
      Contains an external gradient (1D) and hessian (2D scipy CSR) for the
      same set of correction parameters.

The tensor encodes a fit with two channels (positive and negative muon
kinematics), one signal process (jpsi), one uniform-distribution background
process (~10% of jpsi) with a single fully-correlated unconstrained
normalization uncertainty, one unconstrained nuisance per correction
parameter built from the per-corparm SparseHist deltas (interpreted as
differences relative to the nominal), and an external additive likelihood
term (g^T x + 0.5 x^T H x) on the corparms.
"""

import argparse
import os
import time

import h5py
import hist
import numpy as np

import wums.ioutils
from rabbit import tensorwriter
from wums import logging as wums_logging
from wums.sparse_hist import SparseHist


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input-jpsi",
        default="/home/b/bendavid/cvhdev4/work/jpsi_module_corrections_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr.hdf5",
        help="jpsi_module_corrections.py output",
    )
    p.add_argument(
        "--input-grads",
        default="/home/b/bendavid/cvhdev4/work/combinedgrads.hdf5",
        help="aggregategrads.py output",
    )
    p.add_argument(
        "--bkg-fraction",
        type=float,
        default=0.01,
        help="background total normalization as fraction of jpsi",
    )
    p.add_argument(
        "--bkg-norm-uncertainty",
        type=float,
        default=1.5,
        help="prefit relative uncertainty for the unconstrained bkg normalization (size only sets the prefit step)",
    )
    p.add_argument(
        "--systematic-type",
        default="log_normal",
        help="systematic type for the tensor writer (e.g. normal, log_normal)",
    )
    p.add_argument(
        "--reg-sigma",
        type=float,
        default=None,
        help="if set, add an L2 regularization term on the corparms with "
        "equivalent 1-sigma uncertainty of this value; the diagonal "
        "Hessian contribution is 1/sigma^2. Disabled by default.",
    )
    p.add_argument(
        "--reg-zero-diag-only",
        action="store_true",
        help="only regularize corparm indices whose diagonal entry in the "
        "external Hessian is exactly zero (or absent from the sparse "
        "Hessian). Useful for closing the null subspace of a PSD "
        "external Hessian without perturbing constrained directions. "
        "Pass --reg-diag-threshold to broaden to near-null pivots.",
    )
    p.add_argument(
        "--reg-diag-threshold",
        type=float,
        default=None,
        help="threshold on the external-Hessian diagonal used by "
        "--reg-zero-diag-only: indices with diag < threshold are "
        "regularized. Default is strictly zero (diag == 0).",
    )
    p.add_argument(
        "--parameterized",
        action="store_true",
        help="Use parameterized A, e, M scale shift variations instead of "
        "per-module corparm variations.  Skips external likelihood terms "
        "and uses dense tensor mode.",
    )
    p.add_argument(
        "-o", "--output", default="./", help="output directory"
    )
    p.add_argument(
        "--outname", default="jpsi_calibration_tensor", help="output filename (without extension)"
    )
    p.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity (-v for INFO, -vv for DEBUG, "
        "-vvv for NOTSET). Controls rabbit/wums loggers, including the "
        "detailed preconditioner factorization output.",
    )
    return p.parse_args()


def load_jpsi_inputs(filename, parameterized=False):
    """Return (data_plus, data_minus, ideal_plus, ideal_minus, syst_plus, syst_minus, volume_hist).

    When ``parameterized`` is False (default), the two SparseHist syst
    histograms are identified by being the entries in the jpsi_ideal output
    dict that have an extra "corparms" axis.

    When ``parameterized`` is True, the syst histograms are the A/e/M scale
    shift variations (hmuplus_scale, hmuminus_scale) with extra "scale_eta"
    and "param" axes.

    ``volume_hist`` is the quantile-transform bin-volume histogram written by
    jpsi_module_corrections.py under the ``quantile_info`` key. It is
    ``None`` when the histmaker was run with ``--rawAxes`` (no quantile
    transform), in which case ``make_volume_bkg`` falls back to computing
    volumes from the template histogram's axis widths directly.
    """
    with h5py.File(filename, "r") as f:
        jideal = wums.ioutils.pickle_load_h5py(f["jpsi_ideal"])
        jnom = wums.ioutils.pickle_load_h5py(f["jpsi_nom"])

        data_plus = jnom["output"]["hmuplus"].get()
        data_minus = jnom["output"]["hmuminus"].get()
        ideal_plus = jideal["output"]["hmuplus"].get()
        ideal_minus = jideal["output"]["hmuminus"].get()

        if parameterized:
            syst_plus = jideal["output"]["hmuplus_scale"].get()
            syst_minus = jideal["output"]["hmuminus_scale"].get()
        else:
            syst_plus = None
            syst_minus = None
            for key, val in jideal["output"].items():
                obj = val.get()
                ax_names = [a.name for a in obj.axes] if hasattr(obj, "axes") else []
                if "corparms" in ax_names:
                    if syst_plus is None:
                        syst_plus = obj
                    else:
                        syst_minus = obj

        if "quantile_info" in f:
            quantile_info = wums.ioutils.pickle_load_h5py(f["quantile_info"])
            volume_hist = quantile_info["volume"].get()
        else:
            # Raw-axes mode (no quantile transform): no volume histogram
            # was written. make_volume_bkg will fall back to the template
            # hist's own axis widths.
            volume_hist = None

    if syst_plus is None or syst_minus is None:
        kind = "scale shift" if parameterized else "per-corparm SparseHist"
        raise RuntimeError(
            f"Did not find both {kind} objects in {filename}"
        )
    return data_plus, data_minus, ideal_plus, ideal_minus, syst_plus, syst_minus, volume_hist


def make_volume_bkg(template_hist, total, volume_hist):
    """Build a background hist proportional to the bin volume.

    The ``volume_hist`` axes are a subset of ``template_hist`` axes (e.g.
    eta, pt_quant, mass_quant but not phi). The volume is broadcast across
    any axes present in the template but missing from the volume histogram,
    then normalized so the total yield equals ``total``. Per-bin variance is
    set to zero (smooth, statistically perfect prediction).

    When ``volume_hist`` is ``None`` (raw-axes mode), the per-bin volume is
    built as the outer product of each template axis's bin widths (with
    categorical axes contributing unit weight).
    """
    h = hist.Hist(*template_hist.axes, storage=hist.storage.Weight())

    if volume_hist is None:
        # Raw-axes fallback: outer product of axis widths.
        shape = h.values().shape
        vol_broadcast = np.ones(shape, dtype=np.float64)
        for i, ax in enumerate(template_hist.axes):
            widths = getattr(ax, "widths", None)
            if widths is None:
                continue
            reshape = [1] * len(shape)
            reshape[i] = len(widths)
            vol_broadcast = vol_broadcast * np.asarray(widths).reshape(reshape)
    else:
        # Build a mapping from volume_hist axis names to template axis indices
        # so the volume can be broadcast into the full shape.
        vol = volume_hist.values()
        template_names = [a.name for a in template_hist.axes]
        vol_names = [a.name for a in volume_hist.axes]

        # Reshape vol so that matched axes align and unmatched axes get size-1
        # (broadcast) dimensions.
        shape = []
        vol_axis = 0
        for tname in template_names:
            if vol_axis < len(vol_names) and tname == vol_names[vol_axis]:
                shape.append(vol.shape[vol_axis])
                vol_axis += 1
            else:
                shape.append(1)
        vol_broadcast = vol.reshape(shape)

    vol_sum = vol_broadcast.sum()
    if vol_sum > 0:
        h.values()[...] = total * vol_broadcast / vol_sum
    else:
        n = int(np.prod(h.values().shape))
        h.values()[...] = total / n
    h.variances()[...] = 0.0
    return h


def add_external_grad_hess(writer, grad_np, hess_csr, param_names):
    """Wrap the external gradient and hessian and book the term."""
    # Sanity checks
    n = len(param_names)
    if grad_np.shape != (n,):
        raise RuntimeError(
            f"grad shape {grad_np.shape} does not match nparms {n}"
        )
    if hess_csr.shape != (n, n):
        raise RuntimeError(
            f"hess shape {hess_csr.shape} does not match nparms {n}"
        )

    # 1D StrCategory axis (overflow=False since the labels are exhaustive)
    ax_grad = hist.axis.StrCategory(
        list(param_names), name="params", overflow=False
    )
    grad_hist = hist.Hist(ax_grad, storage=hist.storage.Double())
    grad_hist.values()[...] = grad_np.astype(np.float64)

    # 2D StrCategory x StrCategory axis for the hessian, wrapped in SparseHist
    ax_h0 = hist.axis.StrCategory(
        list(param_names), name="params0", overflow=False
    )
    ax_h1 = hist.axis.StrCategory(
        list(param_names), name="params1", overflow=False
    )
    # SparseHist requires shape == product of extents; with overflow=False the
    # extent equals the size, so the original (n, n) CSR fits exactly.
    hess_sh = SparseHist(hess_csr, [ax_h0, ax_h1])

    writer.add_external_likelihood_term(grad=grad_hist, hess=hess_sh, name="corparm_grad_hess")


def add_regularization_term(writer, param_names, sigma=1.0, mask=None):
    """Add an L2 regularization term on the corparms as an extra external
    likelihood contribution: zero gradient and a diagonal Hessian with
    ``1/sigma^2`` on the selected diagonal entries.

    If ``mask`` is provided, it must be a boolean array of length
    ``len(param_names)``; only indices where ``mask[i]`` is True receive
    the ``1/sigma^2`` diagonal contribution (the rest are zero).
    """
    import scipy.sparse

    n = len(param_names)

    # Zero gradient
    ax_grad = hist.axis.StrCategory(
        list(param_names), name="params", overflow=False
    )
    grad_hist = hist.Hist(ax_grad, storage=hist.storage.Double())

    # Diagonal hessian = (1/sigma^2) on selected indices (all indices by default)
    strength = 1.0 / (float(sigma) ** 2)
    diag_values = np.zeros(n, dtype=np.float64)
    if mask is None:
        diag_values[...] = strength
    else:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (n,):
            raise ValueError(
                f"mask shape {mask.shape} does not match nparms {n}"
            )
        diag_values[mask] = strength
    diag_csr = scipy.sparse.diags(diag_values, format="csr")

    ax_h0 = hist.axis.StrCategory(
        list(param_names), name="params0", overflow=False
    )
    ax_h1 = hist.axis.StrCategory(
        list(param_names), name="params1", overflow=False
    )
    hess_sh = SparseHist(diag_csr, [ax_h0, ax_h1])

    writer.add_external_likelihood_term(
        grad=grad_hist, hess=hess_sh, name="corparm_regularization"
    )


def main():
    args = parse_args()

    # Configure logging. Default verbosity=3 (INFO). -v bumps to DEBUG, etc.
    # Using setup_logger with initName="wums" covers both wums.* and rabbit.*
    # child loggers, including rabbit.tensorwriter (preconditioner debug).
    wums_logging.setup_logger(
        __file__, verbosity=3 + args.verbose
    )

    print(f"Loading jpsi inputs from {args.input_jpsi}")
    (
        data_plus,
        data_minus,
        ideal_plus,
        ideal_minus,
        syst_plus,
        syst_minus,
        volume_hist,
    ) = load_jpsi_inputs(args.input_jpsi, parameterized=args.parameterized)

    if args.parameterized:
        # Dense mode, no external likelihood terms.
        print("  parameterized mode (A, e, M scale shift variations)")
        syst_ax_names = [
            a.name for a in syst_plus.axes
            if a.name not in [a.name for a in ideal_plus.axes]
        ]
        nsyst = int(np.prod([syst_plus.axes[n].size for n in syst_ax_names]))
        print(f"  syst axes: {syst_ax_names}, total variations: {nsyst}")
    else:
        nparms = len(syst_plus.axes[-1])
        print(f"  nparms = {nparms}")
        print(f"  syst_plus  nnz = {syst_plus.nnz}")
        print(f"  syst_minus nnz = {syst_minus.nnz}")

        print(f"Loading external gradient/hessian from {args.input_grads}")
        with h5py.File(args.input_grads, "r") as f:
            grad_np = wums.ioutils.pickle_load_h5py(f["grad"])
            hess_csr = wums.ioutils.pickle_load_h5py(f["hess"])
        print(f"  grad shape = {grad_np.shape}")
        print(f"  hess shape = {hess_csr.shape}, nnz = {hess_csr.nnz}")

        if grad_np.shape[0] != nparms or hess_csr.shape != (nparms, nparms):
            raise RuntimeError(
                f"corparm count mismatch: jpsi has {nparms}, grads has "
                f"grad={grad_np.shape}, hess={hess_csr.shape}"
            )

        # The sparse corparm histograms store the variations evaluated at a step
        # size of ``dparm`` in parameter space (set by CVHCorrectorUncertainty).
        # With kfactor=reg_sigma on the systematic, one unit of the fit
        # parameter corresponds to reg_sigma * dparm in physical space.
        # The external gradient and hessian (defined in the unscaled parameter
        # space) must be scaled to match: grad by (reg_sigma * dparm), hess
        # by (reg_sigma * dparm)^2.
        meta_plus = getattr(syst_plus, "metadata", None) or {}
        meta_minus = getattr(syst_minus, "metadata", None) or {}
        dparm_plus = meta_plus.get("dparm")
        dparm_minus = meta_minus.get("dparm")
        if dparm_plus is None or dparm_minus is None:
            raise RuntimeError(
                "dparm metadata missing on syst histograms; re-run "
                "jpsi_module_corrections.py after the metadata support was added"
            )
        if dparm_plus != dparm_minus:
            raise RuntimeError(
                f"dparm mismatch between plus ({dparm_plus}) and minus "
                f"({dparm_minus}) syst histograms"
            )
        dparm = float(dparm_plus)
        # Capture the raw (unscaled) external-Hessian diagonal for the
        # reg-zero-diag-only mask; the threshold is applied in raw
        # external-parameter units.
        ext_diag_raw = np.asarray(hess_csr.diagonal())
        print(f"  scaling external grad by dparm={dparm} and hess by dparm^2")
        grad_np = grad_np * dparm
        hess_csr = hess_csr * (dparm * dparm)

        # Parameter names that label both the per-corparm systematics and the
        # rows/columns of the external gradient/hessian.
        param_names = [f"corparm_{i}" for i in range(nparms)]

    # Build the writer: dense for parameterized, sparse for corparms.
    writer = tensorwriter.TensorWriter(
        sparse=not args.parameterized,
        systematic_type=args.systematic_type,
    )

    # Channels: positive and negative muon kinematics
    writer.add_channel(ideal_plus.axes, "ch_plus")
    writer.add_channel(ideal_minus.axes, "ch_minus")

    # Observed data: jpsi_nom
    writer.add_data(data_plus, "ch_plus")
    writer.add_data(data_minus, "ch_minus")

    # Signal process: jpsi (from jpsi_ideal nominal)
    writer.add_process(ideal_plus, "jpsi", "ch_plus", signal=False)
    writer.add_process(ideal_minus, "jpsi", "ch_minus", signal=False)

    # Background: distributed proportionally to quantile bin volume, ~10% of jpsi total
    bkg_total_plus = args.bkg_fraction * float(ideal_plus.values().sum())
    bkg_total_minus = args.bkg_fraction * float(ideal_minus.values().sum())
    bkg_plus = make_volume_bkg(ideal_plus, bkg_total_plus, volume_hist)
    bkg_minus = make_volume_bkg(ideal_minus, bkg_total_minus, volume_hist)
    # writer.add_process(bkg_plus, "bkg", "ch_plus")
    # writer.add_process(bkg_minus, "bkg", "ch_minus")

    # # Single fully-correlated, unconstrained background normalization syst
    # writer.add_norm_systematic(
    #     "bkg_norm",
    #     "bkg",
    #     "ch_plus",
    #     args.bkg_norm_uncertainty,
    #     constrained=False,
    # )
    # writer.add_norm_systematic(
    #     "bkg_norm",
    #     "bkg",
    #     "ch_minus",
    #     args.bkg_norm_uncertainty,
    #     constrained=False,
    # )

    for chan in ["ch_plus", "ch_minus"]:
        writer.add_norm_systematic(
            "jpsi_norm",
            "jpsi",
            chan,
            1. + 2e-3,
            constrained=False,
            noi = False,
        )

    if args.parameterized:
        # A, e, M scale shift shape systematics.  The multi-systematic
        # dispatch in add_systematic auto-detects the extra tensor axes
        # (scale_eta, param) and books one sub-systematic per combination.
        print("Booking parameterized scale shift systematics for ch_plus")
        t0 = time.time()
        writer.add_systematic(
            syst_plus,
            "scale",
            "jpsi",
            "ch_plus",
            mirror=True,
            as_difference=False,
            constrained=False,
        )
        print(f"  done in {time.time() - t0:.1f}s")

        print("Booking parameterized scale shift systematics for ch_minus")
        t0 = time.time()
        writer.add_systematic(
            syst_minus,
            "scale",
            "jpsi",
            "ch_minus",
            mirror=True,
            as_difference=False,
            constrained=False,
        )
        print(f"  done in {time.time() - t0:.1f}s")

    else:
        # Per-corparm shape systematics built from the SparseHist deltas. Each
        # SparseHist has an extra "corparms" axis that the writer's multi-systematic
        # dispatch detects and iterates: it produces one sub-systematic per
        # corparm (named "corparm_<i>" using the bin label) for *every* bin on the
        # corparms axis, including those with no nonzero variation, so that all
        # nparms nuisances appear in the fit parameter list and can be constrained
        # by the external term. The same name is used in both channels which makes
        # the resulting nuisances fully correlated.
        print("Booking corparm systematics for ch_plus")
        t0 = time.time()
        writer.add_systematic(
            syst_plus,
            "corparm",
            "jpsi",
            "ch_plus",
            mirror=True,
            as_difference=True,
            constrained=False,
            groups=["corparms"],
        )
        print(f"  done in {time.time() - t0:.1f}s")

        print("Booking corparm systematics for ch_minus")
        t0 = time.time()
        writer.add_systematic(
            syst_minus,
            "corparm",
            "jpsi",
            "ch_minus",
            mirror=True,
            as_difference=True,
            constrained=False,
            groups=["corparms"],
        )
        print(f"  done in {time.time() - t0:.1f}s")

        # External additive likelihood term on the corparms
        print("Adding external gradient + hessian likelihood term")
        add_external_grad_hess(writer, grad_np, hess_csr, param_names)

        # L2 regularization on the corparms (diagonal hessian = 1/sigma^2 * I,
        # zero gradient). Disabled unless --reg-sigma is provided. Can be
        # restricted to indices whose raw external-Hessian diagonal is
        # below --reg-diag-threshold, so well-constrained directions
        # aren't perturbed.
        if args.reg_sigma is not None:
            reg_mask = None
            if args.reg_zero_diag_only:
                if args.reg_diag_threshold is None:
                    # Default: regularize only diagonals that are exactly
                    # zero (or absent from the sparse Hessian). Indices
                    # missing from the COO representation have zero diag
                    # by construction, so ``== 0`` catches both.
                    reg_mask = ext_diag_raw == 0.0
                    desc = "raw diag == 0"
                else:
                    reg_mask = ext_diag_raw < args.reg_diag_threshold
                    desc = f"raw diag < {args.reg_diag_threshold}"
                n_reg = int(reg_mask.sum())
                print(
                    f"Adding corparm regularization term on "
                    f"{n_reg}/{len(param_names)} indices with {desc} "
                    f"(sigma={args.reg_sigma})"
                )
            else:
                print(
                    f"Adding corparm regularization term (sigma={args.reg_sigma})"
                )
            add_regularization_term(
                writer, param_names, sigma=args.reg_sigma, mask=reg_mask
            )
        else:
            print("Corparm regularization disabled (no --reg-sigma)")

    # Collect metadata to propagate through the fit to the results.
    if args.parameterized:
        # Read dA, de, dM from the param axis metadata on the syst histogram.
        param_ax = [a for a in syst_plus.axes if a.name == "param"][0]
        scale_eta_ax = [a for a in syst_plus.axes if a.name == "scale_eta"][0]
        calibration_meta = param_ax.metadata.copy()
        calibration_meta["parameterized"] = True
        calibration_meta["scale_eta_axis"] = scale_eta_ax
        calibration_meta["param_axis"] = param_ax
    else:
        calibration_meta = {"dparm": dparm, "parameterized": False}

    # Write the tensor
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    print(f"Writing tensor to {args.output}/{args.outname}.hdf5")
    writer.write(
        outfolder=args.output,
        outfilename=args.outname,
        meta_data_dict={"calibration": calibration_meta},
    )
    print("Done")


if __name__ == "__main__":
    main()
