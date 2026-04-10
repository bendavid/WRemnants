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
        default=0.1,
        help="background total normalization as fraction of jpsi",
    )
    p.add_argument(
        "--bkg-norm-uncertainty",
        type=float,
        default=1.5,
        help="prefit relative uncertainty for the unconstrained bkg normalization (size only sets the prefit step)",
    )
    p.add_argument(
        "-o", "--output", default="./", help="output directory"
    )
    p.add_argument(
        "--outname", default="jpsi_calibration_tensor", help="output filename (without extension)"
    )
    return p.parse_args()


def load_jpsi_inputs(filename):
    """Return (data_plus, data_minus, ideal_plus, ideal_minus, syst_plus, syst_minus).

    The two SparseHist syst histograms are identified by being the entries in
    the jpsi_ideal output dict that have an extra "corparms" axis. The order
    matches the order in which jpsi_module_corrections.py appends them
    (plus then minus).
    """
    with h5py.File(filename, "r") as f:
        jideal = wums.ioutils.pickle_load_h5py(f["jpsi_ideal"])
        jnom = wums.ioutils.pickle_load_h5py(f["jpsi_nom"])

        data_plus = jnom["output"]["hmuplus"].get()
        data_minus = jnom["output"]["hmuminus"].get()
        ideal_plus = jideal["output"]["hmuplus"].get()
        ideal_minus = jideal["output"]["hmuminus"].get()

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

    if syst_plus is None or syst_minus is None:
        raise RuntimeError(
            f"Did not find both per-corparm SparseHist objects in {filename}"
        )
    return data_plus, data_minus, ideal_plus, ideal_minus, syst_plus, syst_minus


def make_uniform_bkg(template_hist, total):
    """Build a flat-distribution background hist matching the template's axes.

    Total normalization is set to ``total``; per-bin variance is set to zero
    (treating the background as a smooth, statistically perfect prediction).
    """
    h = hist.Hist(*template_hist.axes, storage=hist.storage.Weight())
    n = int(np.prod(h.values().shape))
    per_bin = total / n
    h.values()[...] = per_bin
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


def main():
    args = parse_args()

    print(f"Loading jpsi inputs from {args.input_jpsi}")
    (
        data_plus,
        data_minus,
        ideal_plus,
        ideal_minus,
        syst_plus,
        syst_minus,
    ) = load_jpsi_inputs(args.input_jpsi)

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

    # Parameter names that label both the per-corparm systematics and the
    # rows/columns of the external gradient/hessian.
    param_names = [f"corparm_{i}" for i in range(nparms)]

    # Build the writer in sparse mode with normal systematic type so the
    # corparm deltas can be applied with as_difference=True.
    writer = tensorwriter.TensorWriter(
        sparse=True,
        systematic_type="log_normal",
    )

    # Channels: positive and negative muon kinematics
    writer.add_channel(ideal_plus.axes, "ch_plus")
    writer.add_channel(ideal_minus.axes, "ch_minus")

    # Observed data: jpsi_nom
    writer.add_data(data_plus, "ch_plus")
    writer.add_data(data_minus, "ch_minus")

    # Signal process: jpsi (from jpsi_ideal nominal)
    writer.add_process(ideal_plus, "jpsi", "ch_plus", signal=True)
    writer.add_process(ideal_minus, "jpsi", "ch_minus", signal=True)

    # Background: uniform distribution, ~10% of jpsi total
    bkg_total_plus = args.bkg_fraction * float(ideal_plus.values().sum())
    bkg_total_minus = args.bkg_fraction * float(ideal_minus.values().sum())
    bkg_plus = make_uniform_bkg(ideal_plus, bkg_total_plus)
    bkg_minus = make_uniform_bkg(ideal_minus, bkg_total_minus)
    writer.add_process(bkg_plus, "bkg", "ch_plus")
    writer.add_process(bkg_minus, "bkg", "ch_minus")

    # Single fully-correlated, unconstrained background normalization syst
    writer.add_norm_systematic(
        "bkg_norm",
        "bkg",
        "ch_plus",
        args.bkg_norm_uncertainty,
        constrained=False,
    )
    writer.add_norm_systematic(
        "bkg_norm",
        "bkg",
        "ch_minus",
        args.bkg_norm_uncertainty,
        constrained=False,
    )

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
        constrained=True,
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
        constrained=True,
        groups=["corparms"],
    )
    print(f"  done in {time.time() - t0:.1f}s")

    # External additive likelihood term on the corparms
    print("Adding external gradient + hessian likelihood term")
    add_external_grad_hess(writer, grad_np, hess_csr, param_names)

    # Write the tensor
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    print(f"Writing tensor to {args.output}/{args.outname}.hdf5")
    writer.write(outfolder=args.output, outfilename=args.outname)
    print("Done")


if __name__ == "__main__":
    main()
