"""Visualize per-module correction parameters across the CMS tracker.

Produces r-phi, r-z, and 3D scatter views for each parameter type, with
module positions colored by the fitted correction value.

Parameter types:
  0: dx  (local x translation)     5: rz  (local z rotation)
  1: dy  (local y translation)     6: dBz (B-field z offset)
  2: dz  (local z translation)     7: dxi (material scaling)
  3: rx  (local x rotation)
  4: ry  (local y rotation)

Usage:
    python scripts/plotting/muon_calibration/plot_module_corrections.py \\
        --corrections /home/b/bendavid/cvhdev4/work/correctionResults_module.root \\
        --runtree /home/b/bendavid/cvhdev4/work/globalcor_0_1.root \\
        -o plots/
"""

import argparse
import os

import numpy as np
import ROOT


PARM_LABELS = {
    0: r"$\delta x$ (local x transl.) [cm]",
    1: r"$\delta y$ (local y transl.) [cm]",
    2: r"$\delta z$ (local z transl.) [cm]",
    3: r"$\delta\theta_x$ (local x rot.) [rad]",
    4: r"$\delta\theta_y$ (local y rot.) [rad]",
    5: r"$\delta\theta_z$ (local z rot.) [rad]",
    6: r"$\delta B_z$ (B-field offset) [T]",
    7: r"$\delta\xi$ (material scaling)",
}

PARM_SHORT = {
    0: "dx", 1: "dy", 2: "dz",
    3: "rx", 4: "ry", 5: "rz",
    6: "dBz", 7: "dxi",
}

SUBDET_NAMES = {0: "PXB-L1", 1: "PXB", 2: "PXF", 3: "TIB", 4: "TID", 5: "TOB"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corrections",
                   default="/home/b/bendavid/cvhdev4/work/correctionResults_module.root")
    p.add_argument("--runtree",
                   default="/home/b/bendavid/cvhdev4/work/globalcor_0_1.root")
    p.add_argument("-o", "--output", default="./")
    p.add_argument("--parmtypes", type=int, nargs="*", default=None,
                   help="parameter types to plot (default: all)")
    return p.parse_args()


def load_data(corrections_file, runtree_file):
    """Load correction values and module geometry into numpy arrays."""
    f1 = ROOT.TFile.Open(corrections_file)
    f2 = ROOT.TFile.Open(runtree_file)

    parmtree = f1.Get("parmtree")
    runtree = f2.Get("runtree")
    n = runtree.GetEntries()

    corr = np.zeros(n, dtype=np.float64)
    for i, entry in enumerate(parmtree):
        corr[i] = entry.x

    fields = {}
    for name in ["parmtype", "subdet", "layer", "x", "y", "z",
                  "rho", "eta", "phi", "bz", "b0"]:
        fields[name] = np.zeros(n, dtype=np.float64)

    for i, entry in enumerate(runtree):
        for name in fields:
            fields[name][i] = getattr(entry, name)

    fields["parmtype"] = fields["parmtype"].astype(int)
    fields["subdet"] = fields["subdet"].astype(int)
    fields["layer"] = fields["layer"].astype(int)

    f1.Close()
    f2.Close()
    return corr, fields


def plot_parmtype(corr, fields, pt, output_dir):
    """Generate r-phi, r-z, and 3D views for one parameter type."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    mask = fields["parmtype"] == pt
    vals = corr[mask]
    xx = fields["x"][mask]
    yy = fields["y"][mask]
    zz = fields["z"][mask]
    rr = fields["rho"][mask]
    pp = fields["phi"][mask]

    label = PARM_LABELS[pt]
    short = PARM_SHORT[pt]

    # Symmetric color range centered on zero
    vmax = np.percentile(np.abs(vals), 99)
    if vmax == 0:
        vmax = 1e-10
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Sort by absolute value so extreme points are drawn on top
    order = np.argsort(np.abs(vals))

    # --- r-phi view ---
    fig, ax = plt.subplots(figsize=(8, 8))
    sc = ax.scatter(xx[order], yy[order], c=vals[order], s=0.3,
                    cmap="RdBu_r", norm=norm, rasterized=True)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label(label)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_title(f"{short}: r-$\\phi$ view")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{short}_rphi.pdf"))
    plt.close(fig)

    # --- r-z view ---
    fig, ax = plt.subplots(figsize=(12, 5))
    sc = ax.scatter(zz[order], rr[order], c=vals[order], s=0.3,
                    cmap="RdBu_r", norm=norm, rasterized=True)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label(label)
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("r [cm]")
    ax.set_title(f"{short}: r-z view")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{short}_rz.pdf"))
    plt.close(fig)

    # --- 3D view ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(xx[order], yy[order], zz[order], c=vals[order], s=0.3,
                    cmap="RdBu_r", norm=norm, rasterized=True)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
    cb.set_label(label)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_zlabel("z [cm]")
    ax.set_title(f"{short}: 3D view")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{short}_3d.pdf"))
    plt.close(fig)

    print(f"  {short}: {mask.sum()} modules, "
          f"range [{vals.min():.3e}, {vals.max():.3e}], "
          f"color range [-{vmax:.3e}, {vmax:.3e}]")


def main():
    args = parse_args()

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    print("Loading data...")
    corr, fields = load_data(args.corrections, args.runtree)
    print(f"  {len(corr)} total parameters")

    parmtypes = args.parmtypes if args.parmtypes is not None else list(range(8))

    print("Generating plots...")
    for pt in parmtypes:
        if pt not in PARM_LABELS:
            print(f"  skipping unknown parmtype {pt}")
            continue
        plot_parmtype(corr, fields, pt, args.output)

    print("Done")


if __name__ == "__main__":
    main()
