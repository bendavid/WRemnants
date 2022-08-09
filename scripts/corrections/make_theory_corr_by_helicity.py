import numpy as np
import lz4.frame
import pickle
import logging
from wremnants import plot_tools, theory_corrections, theory_tools
from wremnants import boostHistHelpers as hh
from wremnants import common, input_tools, output_tools
import hist
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minnlo_file", type=str, default="w_z_gen_dists.pkl.lz4", help="MiNNLO gen file, denominator in ratio") 
parser.add_argument("-ul", "--corr_ul", type=str, required=True, help="Reference file of the corrections") 
parser.add_argument("-a4", "--corr_a4", type=str, required=True, help="Reference file of the corrections") 
parser.add_argument("-g", "--generator", type=str, choices=["scetlib", ], default="scetlib",
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--proc", type=str, required=True, choices=["z", "wm", "wp"], help="Process")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

if args.proc == "z":
    procName = "ZmumuPostVFP" 
    charge = 0
elif args.proc == "wm":
    procName = "WminusmunuPostVFP" 
    charge = -1
elif args.proc == "wp":
    procName = "WplusmunuPostVFP" 
    charge = 1

minnloh_all = input_tools.read_and_scale(args.minnlo_file, procName, "helicity_moments_scale")
minnloh = minnloh_all[{"muRfact" : 1.j, "muFfact" : 1.j}]

if args.generator == "scetlib":
    # TODO: Fix this, excluding nonsingular is just for now
    sigma_ulh = theory_corrections.read_scetlib_hist(args.corr_ul, charge=charge, nonsing="" if args.proc=="z" else "auto")
    sigma_ulh = hh.makeAbsHist(sigma_ulh, "y")
    sigma4h = theory_corrections.read_scetlib_hist(args.corr_a4, charge=charge, nonsing="" if args.proc=="z" else "auto")
    sigma4h = hh.makeAbsHist(sigma4h, "y")
    a4h = theory_corrections.make_a4_coeff(sigma4h, sigma_ulh)
else:
    raise NotImplementedError("Only SCETlib supported so far")

corrh  = theory_corrections.make_corr_by_helicity(minnloh, sigma_ulh, a4h)

outfile = f"{args.outpath}/{args.generator}HelicityCorr{procName}.pkl.lz4"
if args.postfix:
    outfile = outfile.replace(".pkl.lz4", f"_{args.postfix}.pkl.lz4")

with lz4.frame.open(outfile, "wb") as f:
    pickle.dump({
            procName : {
                f"{args.generator}_minnlo_coeffs" : corrh,
                f"{args.generator}_sigma4_hist" : sigma4h,
                f"{args.generator}_sigmaUL_hist" : sigma_ulh,
                "minnlo_helicity_hist" : minnloh,
            },
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

logging.info("Correction binning is")
for ax in corrh.axes:
    logging.info(f"Axis {ax.name}: {ax.edges}")
logging.info(f"Average correction is {corrh.sum()/np.ones_like(corrh).sum()}")
logging.info(f"Wrote file {outfile}")

