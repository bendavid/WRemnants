import numpy as np
import lz4.frame
import pickle
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools, logging
import hist
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minnlo_file", type=str, required=True, help="MiNNLO gen file, denominator in ratio") 
parser.add_argument("-ul", "--corr-ul", type=str, nargs='+', required=False, help="Reference file of the corrected sigma_UL") 
parser.add_argument("-a4", "--corr-a4", type=str, nargs='+', required=False, help="Reference file of the corrected sigma_4") 
parser.add_argument("-f", "--other-coeff-files", type=str, nargs='+', required=False, help="Files with other Ai coefficients") 
parser.add_argument("-c", "--other-coeffs", type=int, nargs='+', default=[0,1,2,3], help="Files with other Ai coefficients") 
parser.add_argument("-g", "--generator", type=str, choices=["scetlib", ], default="scetlib",
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--verbose", choices=[0,1,2,3,4], default=3, help="Output verbosity level")
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")

args = parser.parse_args()

logger = logging.setup_base_logger("makeDataMCStackPlot", args.verbose)

if all(getattr(args, x) is None for x in ["corr_ul", "corr_a4", "other_coeffs"]):
    raise ValueError("Must specify at least one correction file")

def read_corr(procName, generator, corr, isA4=False):
    charge = 0 if procName[0] == "Z" else (1 if "Wplus" in procName else -1)
    if args.generator == "scetlib":
        h = input_tools.read_scetlib_hist(corr, charge=charge, flip_y_sign=isA4)
        h = hh.makeAbsHist(h, "y")
    else:
        raise NotImplementedError("Only SCETlib supported so far")

    return h

filesByProc = {"ZmumuPostVFP" : []} if args.proc == "z" else {"WplusmunuPostVFP" : [], "WminusmunuPostVFP" : []}

if args.corr_ul:
    if args.corr_a4 and len(args.corr_ul) != len(args.corr_a4):
        raise ValueError("Must pass the same number of UL and sigma4 files")

    def readCorrs(args, idx):
        return [getattr(args, corr)[idx] for corr in ["corr_ul", "corr_a4"] if hasattr(args, corr)]
    
    if args.proc == "z":
        if len(args.corr_ul) != 1:
            raise ValueError("Only one file expected for Z")
        filesByProc["ZmumuPostVFP"] = readCorrs(args, 0)
    elif args.proc == "w":
        if len(args.corr_ul) != 2:
            raise ValueError("Requires two files for W (W+ and W-)")
        plus_idx = 0 if args.corr_ul and "Wp" in args.corr_ul[0] else 1
        plus_idx4 = 0 if args.corr_a4 and "Wp" in args.corr_a4[0] else 1
        if plus_idx != plus_idx4:
            raise ValueError("Expected a common order between UL and sigma4 files!")
        filesByProc["WplusmunuPostVFP"] = readCorrs(args, plusidx)
        filesBuProc["WminusmunuPostVFP"] = readCorrs(args, not plusidx)

minnloh = hh.sumHists([input_tools.read_mu_hist_combine_tau(args.minnlo_file, proc, "nominal_gen_helicity_moments_scale") for proc in filesByProc.keys()])
minnloh = minnloh[{"muRfact" : 1.j, "muFfact" : 1.j}]

sigma_ulh = hh.sumHists([read_corr(procName, args.generator, corr_files[0]) for procName, corr_files in filesByProc.items()]) if args.corr_ul else None
sigma4h = hh.sumHists([read_corr(procName, args.generator, corr_files[1], isA4=True) for procName, corr_files in filesByProc.items()]) if args.corr_a4 else None

dyturbo_coeffs = None
ptV_binning = common.ptV_binning[:-4]+list(range(30,110,10))
binning = {"qT" : ptV_binning, "y" : [-4.5+.5*i for i in range(19)]}
if args.other_coeff_files and args.other_coeffs:
    dyturbo_coeffs = hh.sumHists([input_tools.read_dyturbo_angular_coeffs(f, rebin=binning, 
                                        absy=True, add_axes=[minnloh.axes["massVgen"]]) for f in args.other_coeff_files])

# Workaround for the different axis names
binning["ptVgen"] = ptV_binning
corrh = theory_corrections.make_corr_by_helicity(minnloh, sigma_ulh, sigma4h, coeff_hist=dyturbo_coeffs, coeffs_from_hist=args.other_coeffs, binning=binning)

outName = "Z" if args.proc == "z" else "W"
outfile = f"{args.outpath}/{args.generator}HelicityCorr{outName}.pkl.lz4"
if args.postfix:
    outfile = outfile.replace(".pkl.lz4", f"_{args.postfix}.pkl.lz4")

out_dict = {
    f"{args.generator}_minnlo_coeffs" : corrh,
    f"{args.generator}_sigma4_hist" : sigma4h,
    f"{args.generator}_sigmaUL_hist" : sigma_ulh,
    "minnlo_helicity_hist" : minnloh,
    f"other_coeffs" : dyturbo_coeffs,
}
with lz4.frame.open(outfile, "wb") as f:
    pickle.dump({ outName : out_dict,
            "meta_data" : output_tools.metaInfoDict(args=args),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

logger.info("Correction binning is")
for ax in corrh.axes:
    logger.info(f"Axis {ax.name}: {ax.edges}")
for i in range(-1, 8):
    h = corrh[{'helicity' : complex(0, i)}]
    helcorr = np.average(hh.divideHists(h[{'corr' : True}], h[{'corr' : False}]) )
    logger.info(f"Average correction for sigma_{i} is {helcorr}")
logger.info(f"Wrote file {outfile}")

