import numpy as np
import lz4.frame
import pickle
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools
import pathlib
import hist
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minnlo_file", type=str, default="w_z_gen_dists.pkl.lz4", help="MiNNLO gen file, denominator in ratio") 
parser.add_argument("-c", "--corr_files", type=str, nargs='+', required=True, help="Reference files for the corrections (both W+ and W- for the W)") 
parser.add_argument("-g", "--generator", type=str, choices=["dyturbo", "scetlib", "scetlib_dyturbo", "matrix_radish"], 
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")
parser.add_argument("--axes", nargs="*", type=str, default=None, help="Use only specified axes in hist")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--no-color-logger", action="store_false", dest="color_logger", default=False, 
                    help="Do not use logging with colors")

args = parser.parse_args()

logger = common.setup_logger("make_theory_corr", 4 if args.debug else 3, args.color_logger)

def read_corr(procName, generator, corr_file):
    charge = 0 if procName[0] == "Z" else (1 if "Wplus" in procName else -1)
    if "scetlib" in generator:
        if "dyturbo" in generator == "scetlib_dyturbo":
            scetlib_files = [x for x in corr_files if pathlib.Path(x).suffix == ".pkl"]
            if len(scetlib_files) != 2:
                raise ValueError(f"scetlib_dyturbo correction requires two SCETlib files (resummed and FO nonsingular). Found {len(scetlib_files)}")
            nlo_nons_idx = 0 if "nlo" in scetlib_files[0] else 1
            resumf = scetlib_files[~nlo_nons_idx]
            nlo_nonsf = scetlib_files[nlo_nons_idx]

            dyturbo_files = [x for x in corr_files if pathlib.Path(x).suffix == ".txt"]
            if len(dyturbo_files) != 1:
                raise ValueError("scetlib_dyturbo correction requires one DYTurbo file (fixed order contribution)")

            numh = input_tools.read_matched_scetlib_dyturbo_hist(resumf, nlo_nonsf, dyturbo_files[0], args.axes, charge=charge)
        else:
            nons = "auto"
            if not os.path.isfile(corr_file.replace(".", "_nons.")):
                nons = ""
                logger.warning("No nonsingular file found for SCETlib. Will skip it.")
            numh = input_tools.read_scetlib_hist(corr_file, charge=charge, nonsing=nons)
        if "Y" in numh.axes.name:
            numh = hh.makeAbsHist(numh, "Y")
        return numh 
    else:
        if args.generator == "matrix_radish":
            h = input_tools.read_matrixRadish_hist(corr_file, "ptVgen")
        else:
            axnames = args.axes
            if not axnames:
                axnames = ("y", "pt") if "2d" in corr_file else ("pt")
            h = input_tools.read_dyturbo_hist(corr_file.split(":"), axes=axnames, charge=charge)
            if "y" in h.axes.name:
                h = hh.makeAbsHist(h, "y")

        vars_ax = h.axes["vars"] if "vars" in h.axes.name else hist.axis.Integer(0, 1, flow=False, name="vars") 
        axes.append(vars_ax)

        h5D = hist.Hist(*axes)
        # Leave off the overflow, we won't use it anyway
        h5D[...] = np.reshape(h.values(), h5D.shape)
        numh = h5D
    return numh

if args.proc == "z":
    filesByProc = { "ZmumuPostVFP" : args.corr_files }
elif args.proc == "w":
    wpfiles = filter(lambda x: "wp" in x.lower(), args.corr_files)
    wmfiles = filter(lambda x: "wm" in x.lower(), args.corr_files)
    if len(wpfiles) != len(wmfiles):
        raise ValueError(f"Expected equal number of files for W+ and W-, found {len(wpfiles)} (Wp) and {len(wmfiles)} (Wm)")
    filesByProc = { "WplusmunuPostVFP" : wpfiles,
        "WminusmunuPostVFP" : wmfiles}

minnloh = input_tools.read_all_and_scale(args.minnlo_file, list(filesByProc.keys()), ["nominal_gen"])[0]

if "y" in minnloh.axes.name:
    minnloh = hh.makeAbsHist(minnloh, "y")

# Get more stats in the correction
add_taus = True 
if add_taus:
    logger.info("Combining muon and tau decay samples for increased stats")
    from wremnants.datasets.datasetDict_v9 import Z_TAU_TO_LEP_RATIO,BR_TAUToMU
    taus = ["WplustaunuPostVFP", "WminustaunuPostVFP"] if args.proc == "w" else ["ZtautauPostVFP"]
    taush = input_tools.read_all_and_scale(args.minnlo_file, taus, ["nominal_gen"])[0]
    if "y" in taush.axes.name:
        taush = hh.makeAbsHist(taush, "y")
    # Rescale taus to mu to effectively have more stats
    minnloh = 0.5*(minnloh + taush/(Z_TAU_TO_LEP_RATIO if args.proc == "z" else BR_TAUToMU))

numh = hh.sumHists([read_corr(procName, args.generator, corr_file) for procName, corr_file in filesByProc.items()])

if numh.ndim-1 < minnloh.ndim:
    ax_map = {
        "ptVgen" : "qT",
        "absYVgen" : "absY",
        "massVgen" : "Q",
        "chargeVgen" : "charge",
    }
    axes = []
    for ax in minnloh.axes.name:
        axname = ax_map[ax]
        if axname in numh.axes.name:
            axes.append(numh.axes[axname])
        else:
            minnlo_ax = minnloh.axes[ax]
            # TODO: Should be a little careful because this won't include overflow, as long as the
            # axis range is large enough, it shouldn't matter much
            axes.append(hist.axis.Regular(1, minnlo_ax.edges[0], minnlo_ax.edges[-1], 
                underflow=minnlo_ax.traits.underflow, overflow=minnlo_ax.traits.overflow, name=ax))
    # NOTE: This leaves out the underflow, but there can't be any underflow anyway
    numh = hist.Hist(*axes, storage=numh._storage_type(), data=np.reshape(numh.view(), [ax.size for ax in axes]))


corrh_unc  = theory_corrections.make_corr_from_ratio(minnloh, numh)
corrh = hist.Hist(*corrh_unc.axes, name=corrh_unc.name, storage=hist.storage.Double(), data=corrh_unc.values(flow=True))

if args.postfix:
    args.generator += args.postfix
outfile = f"{args.outpath}/{args.generator}Corr{args.proc.upper()}.pkl.lz4"

with lz4.frame.open(outfile, "wb") as f:
    pickle.dump({
            args.proc.upper() : {
                f"{args.generator}_minnlo_ratio" : corrh,
                f"{args.generator}_hist" : numh,
                "minnlo_ref_hist" : minnloh,
            },
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

logger.info("Correction binning is")
for ax in corrh.axes:
    logger.info(f"Axis {ax.name}: {ax.edges}")
logger.info(f"Average correction is {corrh.sum()/np.ones_like(corrh).sum()}")
logger.info(f"Wrote file {outfile}")
