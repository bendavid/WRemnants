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
parser.add_argument("-c", "--corr_file", type=str, required=True, help="Reference file of the corrections") 
parser.add_argument("-g", "--generator", type=str, choices=["dyturbo", "scetlib", "matrix_radish"], 
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

if args.proc == "z":
    procNames = ["ZmumuPostVFP"]
    charges = [0]
elif args.proc == "w":
    procName = ["WminusmunuPostVFP", "WplusmunuPostVFP"]
    charges = [-1, 1]
elif args.proc == "wp":
    procName = "WplusmunuPostVFP" 
    charge = 1

minnloh = input_tools.read_and_scale(args.minnlo_file, procName, "nominal_gen")

if args.generator == "scetlib":
    numh = theory_corrections.read_scetlib_hist(args.corr_file, charge=charge)
    numh = hh.makeAbsHist(numh, "y")
else:
    if args.generator == "matrix_radish":
        h = theory_corrections.read_matrixRadish_hist(args.corr_file, "ptVgen")
    else:
        h = theory_corrections.read_dyturbo_hist(args.corr_file.split(":"), axis="ptVgen")

    charge_ax = minnloh.axes["chargeVgen"]
    absy_ax = hist.axis.Regular(1, 0, 10, flow=True, name="absYVgen")
    mass_ax = hist.axis.Regular(1, 60, 120, flow=True, name="massVgen")
    vars_ax = hist.axis.Integer(0, 1, flow=False, name="vars")
    h5D = hist.Hist(mass_ax, absy_ax, *h.axes, charge_ax, vars_ax)
    h5D[...] = h.values(flow=True)[np.newaxis, np.newaxis,:,np.newaxis,np.newaxis]
    numh = h5D

corrh_unc  = theory_corrections.make_corr_from_ratio(minnloh, numh)
corrh = hist.Hist(*corrh_unc.axes, name=corrh_unc.name, storage=hist.storage.Double(), data=corrh_unc.values(flow=True))

procName = "Z" if args.proc == "z" else "W"
outfile = f"{args.outpath}/{args.generator}Corr{procName}.pkl.lz4"
if args.postfix:
    outfile = outfile.replace(".pkl.lz4", f"_{args.postfix}.pkl.lz4")

with lz4.frame.open(outfile, "wb") as f:
    pickle.dump({
            procName : {
                f"{args.generator}_minnlo_ratio" : corrh,
                f"{args.generator}_hist" : numh,
                "minnlo_ref_hist" : minnloh,
            },
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

logging.info("Correction binning is")
for ax in corrh.axes:
    logging.info(f"Axis {ax.name}: {ax.edges}")
logging.info(f"Wrote file {outfile}")
