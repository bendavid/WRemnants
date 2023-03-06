import numpy as np
import lz4.frame
import pickle
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools, logging
import hist
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minnlo_file", type=str, default="w_z_gen_dists.pkl.lz4", help="MiNNLO gen file, denominator in ratio") 
parser.add_argument("-ul", "--corr_ul", type=str, nargs='+', required=True, help="Reference file of the corrected sigma_UL") 
parser.add_argument("-a4", "--corr_a4", type=str, nargs='+', required=True, help="Reference file of the corrected sigma_4") 
parser.add_argument("-g", "--generator", type=str, choices=["scetlib", ], default="scetlib",
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")

args = parser.parse_args()

logger = logging.setup_base_logger("makeDataMCStackPlot", args.debug)

def read_corr(procName, generator, corr, isA4=False):
    charge = 0 if procName[0] == "Z" else (1 if "Wplus" in procName else -1)
    if args.generator == "scetlib":
        h = input_tools.read_scetlib_hist(corr, charge=charge, flip_y_sign=isA4)
        h = hh.makeAbsHist(h, "y")
    else:
        raise NotImplementedError("Only SCETlib supported so far")

    return h

if len(args.corr_ul) != len(args.corr_a4):
    raise ValueError("Must pass the same number of UL and sigma4 files")
if args.proc == "z":
    if len(args.corr_ul) != 1:
        raise ValueError("Only one file expected for Z")
    filesByProc = { "ZmumuPostVFP" : (args.corr_ul[0], args.corr_a4[0]) }
elif args.proc == "w":
    if len(args.corr_ul) != 2:
        raise ValueError("Requires two files for W (W+ and W-)")
    plus_idx = 0 if "Wp" in args.corr_ul[0] else 1
    plus_idx4 = 0 if "Wp" in args.corr_a4[0] else 1
    filesByProc = { "WplusmunuPostVFP" : (args.corr_ul[plus_idx], args.corr_a4[plus_idx4]),
        "WminusmunuPostVFP" : (args.corr_ul[not plus_idx], args.corr_a4[not plus_idx4])}

# TODO: Read taus here already
minnloh = input_tools.read_all_and_scale(args.minnlo_file, list(filesByProc.keys()), ["helicity_moments_scale"])[0]
# # TODO: Really need to reuse code here and for the inclusive corr
# Get more stats in the correction
add_taus = True 
if add_taus:
    logger.info("Combine muon and tau decay samples for increased stats")
    from wremnants.datasets.datasetDict_v9 import Z_TAU_TO_LEP_RATIO,BR_TAUToMU
    taus = ["WplustaunuPostVFP", "WminustaunuPostVFP"] if args.proc == "w" else ["ZtautauPostVFP"]
    taush = input_tools.read_all_and_scale(args.minnlo_file, taus, ["helicity_moments_scale"])[0]
    # Rescale taus to mu to effectively have more stats
    minnloh = 0.5*(minnloh + taush/(Z_TAU_TO_LEP_RATIO if args.proc == "z" else BR_TAUToMU))
minnloh = minnloh[{"muRfact" : 1.j, "muFfact" : 1.j}]

sigma_ulh = hh.sumHists([read_corr(procName, args.generator, corr_files[0]) for procName, corr_files in filesByProc.items()])
sigma4h = hh.sumHists([read_corr(procName, args.generator, corr_files[1], isA4=True) for procName, corr_files in filesByProc.items()])

corrh = theory_corrections.make_corr_by_helicity(minnloh, sigma_ulh, sigma4h)

outName = "Z" if args.proc == "z" else "W"
outfile = f"{args.outpath}/{args.generator}HelicityCorr{outName}.pkl.lz4"
if args.postfix:
    outfile = outfile.replace(".pkl.lz4", f"_{args.postfix}.pkl.lz4")

out_dict = {
    f"{args.generator}_minnlo_coeffs" : corrh,
    f"{args.generator}_sigma4_hist" : sigma4h,
    f"{args.generator}_sigmaUL_hist" : sigma_ulh,
    "minnlo_helicity_hist" : minnloh,
}
with lz4.frame.open(outfile, "wb") as f:
    pickle.dump({ outName : out_dict,
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

logger.info("Correction binning is")
for ax in corrh.axes:
    logger.info(f"Axis {ax.name}: {ax.edges}")
helcorr = np.average(corrh[{'helicity' : -1.j, 'corr' : True}])
logger.info(f"Average correction for sigmaUL is {helcorr}")
logger.info(f"Wrote file {outfile}")

