"""Convert rabbit fitresults to correction output files.

Supports two modes, detected automatically from the calibration metadata
embedded in the fitresults:

**Corparm mode** (``parameterized=False``):
    Reads corparm fit values, scales by ``dparm``, and writes a
    correctionResults ROOT file with ``idxmaptree`` and ``parmtree``
    TTrees (the standard module-level format).

**Parameterized mode** (``parameterized=True``):
    Reads the fitted A, e, M scale parameters, scales each by its
    corresponding step size (dA, de, dM), and writes the results as
    an HDF5 file containing a ``hist.Hist`` with ``scale_eta`` and
    ``param`` axes (values and variances).

Usage:
    python fitresults_to_correctionResults.py \\
        --fitresults fitresults.hdf5 \\
        -o correctionResults.root   # or .hdf5 for parameterized
"""

import argparse
import array
import os

import h5py
import hist
import numpy as np
import scipy.linalg
import scipy.stats

import wums.ioutils


def _compat_chi2(values, cov, label):
    """Print a chi2 compatibility test of ``values`` with zero.

    ``chi2 = values^T cov^{-1} values`` under the null hypothesis
    ``values = 0`` is distributed as ``chi2_ndf`` with ``ndf = len(values)``.
    The chi2 is invariant under diagonal rescaling of ``values`` (since
    ``cov`` is rescaled by the same factor on both sides), so it is the
    same whether computed in raw or dparm-scaled units.

    Falls back to a diagonal-only chi2 (ignoring correlations) when
    ``cov`` is ``None`` or the full-matrix solve fails.
    """
    values = np.asarray(values, dtype=np.float64).ravel()
    ndf = int(values.size)
    if ndf == 0:
        print(f"  [{label}] no parameters — skipped")
        return

    diag_var = None
    if cov is not None:
        diag_var = np.diag(cov)

    chi2 = None
    mode = None
    if cov is not None:
        try:
            # Symmetrize numerical asymmetries and solve. cho_solve
            # requires PD; fall back to lstsq via pinvh on failure.
            C = 0.5 * (cov + cov.T)
            cf = scipy.linalg.cho_factor(C, lower=True)
            y = scipy.linalg.cho_solve(cf, values)
            chi2 = float(values @ y)
            mode = "full-cov"
        except Exception:
            try:
                C = 0.5 * (cov + cov.T)
                C_inv = scipy.linalg.pinvh(C)
                chi2 = float(values @ (C_inv @ values))
                mode = "full-cov (pinv)"
            except Exception:
                chi2 = None

    if chi2 is None:
        # Diagonal-only fallback.
        if diag_var is None:
            print(
                f"  [{label}] no covariance available — compatibility "
                f"test skipped"
            )
            return
        var_safe = np.where(np.isfinite(diag_var) & (diag_var > 0),
                            diag_var, np.inf)
        chi2 = float(np.sum(values * values / var_safe))
        mode = "diagonal-only"

    pval = float(scipy.stats.chi2.sf(chi2, ndf))
    # Equivalent z-score for a two-sided normal test.
    if pval > 0.0:
        z = float(scipy.stats.norm.isf(0.5 * pval))
    else:
        z = float("inf")
    print(
        f"  [{label}] chi2/ndf = {chi2:.2f} / {ndf} = {chi2/ndf:.3f}, "
        f"p = {pval:.3e}, z = {z:.2f}  ({mode})"
    )


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--fitresults",
        required=True,
        help="rabbit fitresults HDF5 file",
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="output file path (default: correctionResults.root for "
        "corparm mode, correctionResults.hdf5 for parameterized mode)",
    )
    return p.parse_args()


def write_corparm_root(output, corparm_values, corparm_errors):
    """Write corparm results as a correctionResults ROOT file."""
    import ROOT

    nparms = len(corparm_values)

    outdir = os.path.dirname(output)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)

    fout = ROOT.TFile.Open(output, "RECREATE")

    # idxmaptree: identity mapping (idx[i] = i).
    idxmaptree = ROOT.TTree("idxmaptree", "")
    idx = array.array("I", [0])  # UInt_t
    idxmaptree.Branch("idx", idx, "idx/i")
    for i in range(nparms):
        idx[0] = i
        idxmaptree.Fill()
    idxmaptree.Write()

    # parmtree: correction values and uncertainties.
    parmtree = ROOT.TTree("parmtree", "")
    x = array.array("f", [0.0])  # Float_t
    err = array.array("f", [0.0])  # Float_t
    parmtree.Branch("x", x, "x/F")
    parmtree.Branch("err", err, "err/F")
    for i in range(nparms):
        x[0] = float(corparm_values[i])
        err[0] = float(corparm_errors[i])
        parmtree.Fill()
    parmtree.Write()

    fout.Close()
    print(f"Done: {nparms} parameters written to {output}")


def write_parameterized_hdf5(output, scale_hist):
    """Write parameterized A, e, M results as an HDF5 file."""
    outdir = os.path.dirname(output)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)

    with h5py.File(output, "w") as f:
        wums.ioutils.pickle_dump_h5py("scale_corrections", scale_hist, f)

    print(f"Done: scale corrections written to {output}")


def main():
    args = parse_args()

    print(f"Reading fitresults from {args.fitresults}")
    with h5py.File(args.fitresults, "r") as f:
        results = wums.ioutils.pickle_load_h5py(f["results"])
        meta = wums.ioutils.pickle_load_h5py(f["meta"])
        parms_proxy = results["parms"]
        parms_hist = parms_proxy.get() if hasattr(parms_proxy, "get") else parms_proxy
        cov_hist = None
        if "cov" in results:
            cov_proxy = results["cov"]
            cov_hist = cov_proxy.get() if hasattr(cov_proxy, "get") else cov_proxy

    # Read calibration metadata propagated from the input tensor.
    cal_meta = meta.get("meta_info_input", {}).get("calibration", {})
    if not cal_meta:
        raise RuntimeError(
            "No calibration metadata found in fitresults.  Re-run "
            "make_jpsi_calibration_tensor.py to embed it."
        )
    print(f"  calibration metadata: {cal_meta}")

    parameterized = cal_meta.get("parameterized", False)

    # Extract parameter names, values, and variances.
    parm_names = list(parms_hist.axes[0])
    parm_values = parms_hist.values()
    parm_variances = parms_hist.variances()

    if parameterized:
        # --- Parameterized A, e, M mode ---
        dA = float(cal_meta["dA"])
        de = float(cal_meta["de"])
        dM = float(cal_meta["dM"])
        scale_eta_axis = cal_meta["scale_eta_axis"]
        param_axis = cal_meta["param_axis"]
        neta = scale_eta_axis.size
        param_labels = list(param_axis)
        dparms = {"A": dA, "e": de, "M": dM}
        print(f"  parameterized mode: dA={dA}, de={de}, dM={dM}, neta={neta}")

        # Build output histogram with Weight storage (value + variance).
        out_hist = hist.Hist(
            scale_eta_axis, param_axis, storage=hist.storage.Weight()
        )

        # Match fit parameters to (eta_bin, param) bins.
        # Parameter names are "scale_{eta_idx}_{A|e|M}".
        prefix = "scale_"
        # Collect fit-parameter indices per scale (eta_idx, param_label)
        # pair for the subsequent compatibility test.
        scale_fit_idx_by_label = {lbl: [] for lbl in dparms}
        scale_fit_idx_all = []
        scale_fit_labels_all = []
        for iparm, name in enumerate(parm_names):
            if not name.startswith(prefix):
                continue
            # Parse "scale_{eta_idx}_{param}" from the name.
            rest = name[len(prefix):]
            parts = rest.rsplit("_", 1)
            if len(parts) != 2:
                continue
            eta_idx_str, param_label = parts
            if param_label not in dparms:
                continue
            try:
                eta_idx = int(eta_idx_str)
            except ValueError:
                continue

            dp = dparms[param_label]
            val = parm_values[iparm] * dp
            var = parm_variances[iparm] if parm_variances is not None else 0.0
            if not (np.isfinite(var) and var > 0):
                var = 0.0
            var = var * dp * dp

            iparam = param_labels.index(param_label)
            out_hist.values()[eta_idx, iparam] = val
            out_hist.variances()[eta_idx, iparam] = var

            scale_fit_idx_by_label[param_label].append(iparm)
            scale_fit_idx_all.append(iparm)
            scale_fit_labels_all.append(param_label)

        print(
            f"  A values: min={out_hist.values()[:, 0].min():.6e}, "
            f"max={out_hist.values()[:, 0].max():.6e}"
        )
        print(
            f"  e values: min={out_hist.values()[:, 1].min():.6e}, "
            f"max={out_hist.values()[:, 1].max():.6e}"
        )
        print(
            f"  M values: min={out_hist.values()[:, 2].min():.6e}, "
            f"max={out_hist.values()[:, 2].max():.6e}"
        )

        # ------------------------------------------------------------------
        # Statistical compatibility of the scale parameters with zero.
        # Uses the postfit parameter covariance (stored as a hist with
        # parms_x/parms_y string axes). chi2 = p^T C^{-1} p under the
        # null p=0 follows chi2_ndf, where ndf = number of params in
        # the block. dparm rescaling cancels between p and cov so the
        # test is identical whether computed in raw or scaled units —
        # we use the raw (fitted) values.
        # ------------------------------------------------------------------
        print("Compatibility of scale parameters with zero:")
        scale_fit_idx_all = np.asarray(scale_fit_idx_all, dtype=int)
        if cov_hist is not None:
            cov_values_full = cov_hist.values()
            cov_all = cov_values_full[
                np.ix_(scale_fit_idx_all, scale_fit_idx_all)
            ]
        else:
            cov_all = None
            print(
                "  WARNING: no 'cov' hist in fitresults — did the fit run "
                "with --noHessian?  Falling back to diagonal-only chi2 "
                "(ignores correlations)."
            )
        _compat_chi2(
            parm_values[scale_fit_idx_all], cov_all,
            f"all scale (A+e+M, {len(scale_fit_idx_all)} params)",
        )
        for lbl in ("A", "e", "M"):
            idx = np.asarray(scale_fit_idx_by_label[lbl], dtype=int)
            if idx.size == 0:
                continue
            if cov_hist is not None:
                cov_sub = cov_hist.values()[np.ix_(idx, idx)]
            else:
                cov_sub = None
            _compat_chi2(
                parm_values[idx], cov_sub,
                f"{lbl} only ({idx.size} params)",
            )

        output = args.output or "correctionResults.hdf5"
        print(f"Writing {output}")
        write_parameterized_hdf5(output, out_hist)

    else:
        # --- Corparm mode ---
        dparm = float(cal_meta["dparm"])

        # Select only the corparm parameters.
        corparm_mask = np.array(
            [name.startswith("corparm_") for name in parm_names]
        )
        corparm_indices = np.where(corparm_mask)[0]
        nparms = len(corparm_indices)
        print(f"  found {nparms} corparm parameters")

        if nparms == 0:
            raise RuntimeError(
                f"No corparm parameters found in fitresults. "
                f"Available: {parm_names[:10]}..."
            )

        # Extract corparm values and uncertainties, scaled by dparm.
        corparm_values = parm_values[corparm_indices] * dparm
        if parm_variances is not None:
            var = parm_variances[corparm_indices]
            corparm_errors = np.sqrt(np.where(
                np.isfinite(var) & (var > 0), var, 0.0
            )) * dparm
        else:
            corparm_errors = np.zeros(nparms)

        print(f"  dparm = {dparm}, scaling fitted values by dparm")
        print(
            f"  corparm values: min={corparm_values.min():.6e}, "
            f"max={corparm_values.max():.6e}, "
            f"mean={corparm_values.mean():.6e}"
        )
        print(
            f"  corparm errors: min={corparm_errors.min():.6e}, "
            f"max={corparm_errors.max():.6e}"
        )

        output = args.output or "correctionResults.root"
        print(f"Writing {output}")
        write_corparm_root(output, corparm_values, corparm_errors)


if __name__ == "__main__":
    main()
