"""Plot fitted A, e, M parameterized scale corrections as a function of eta.

Usage:
    python scripts/plotting/muon_calibration/plot_parameterized_corrections.py \
        --input /home/b/bendavid/cvhdev4/work/correctionResults_param.hdf5 \
        -o plots/
"""

import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np

import wums.ioutils


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input",
        default="/home/b/bendavid/cvhdev4/work/correctionResults_param.hdf5",
        help="correctionResults HDF5 from fitresults_to_correctionResults.py",
    )
    p.add_argument(
        "-o",
        "--output",
        default="./",
        help="output directory for plots",
    )
    return p.parse_args()


def main():
    args = parse_args()

    with h5py.File(args.input, "r") as f:
        h = wums.ioutils.pickle_load_h5py(f["scale_corrections"])

    eta_axis = h.axes["scale_eta"]
    param_axis = h.axes["param"]
    eta_centers = eta_axis.centers
    eta_edges = eta_axis.edges
    param_labels = list(param_axis)

    values = h.values()
    variances = h.variances()
    errors = np.sqrt(np.where(np.isfinite(variances) & (variances > 0), variances, 0.0))

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    # Individual parameter plots.
    for iparam, label in enumerate(param_labels):
        fig, ax = plt.subplots(figsize=(8, 5))
        vals = values[:, iparam]
        errs = errors[:, iparam]

        ax.errorbar(
            eta_centers, vals, yerr=errs,
            fmt="o", markersize=4, capsize=2, label=label,
        )
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_xlabel(r"$\eta$")
        ax.set_ylabel(f"Fitted {label}")
        ax.set_xlim(eta_edges[0], eta_edges[-1])
        ax.legend()
        ax.set_title(f"Parameterized scale correction: {label}")
        fig.tight_layout()

        outpath = os.path.join(args.output, f"scale_correction_{label}.pdf")
        fig.savefig(outpath)
        plt.close(fig)
        print(f"Saved {outpath}")

    # Combined plot with all three parameters on separate panels.
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    for iparam, (label, ax) in enumerate(zip(param_labels, axes)):
        vals = values[:, iparam]
        errs = errors[:, iparam]

        ax.errorbar(
            eta_centers, vals, yerr=errs,
            fmt="o", markersize=4, capsize=2, color=f"C{iparam}",
        )
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_ylabel(f"{label}")
        ax.set_xlim(eta_edges[0], eta_edges[-1])

    axes[-1].set_xlabel(r"$\eta$")
    axes[0].set_title("Parameterized scale corrections")
    fig.tight_layout()

    outpath = os.path.join(args.output, "scale_corrections_combined.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"Saved {outpath}")


if __name__ == "__main__":
    main()
