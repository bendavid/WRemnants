import numpy as np
import lz4.frame
import pickle
import logging
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools
import hist
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--minnlo_file", type=str, default="w_z_gen_dists.pkl.lz4", help="MiNNLO gen file, denominator in ratio") 
parser.add_argument("-c", "--corr_files", type=str, nargs='+', required=True, help="Reference files for the corrections (both W+ and W- for the W)") 
parser.add_argument("-g", "--generator", type=str, choices=["dyturbo", "scetlib", "matrix_radish"], 
    required=True, help="Generator used to produce correction hist")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

def read_corr(procName, generator, corr_file):
    charge = 0 if procName[0] == "Z" else (1 if "Wplus" in procName else -1)
    if generator == "scetlib":
        nons = "auto"
        if not os.path.isfile(corr_file.replace(".", "_nons.")):
            nons = ""
            logging.warning("Non nonsingular file found for SCETlib. Will skip it. This is an approximation!")
        numh = input_tools.read_scetlib_hist(corr_file, charge=charge, nonsing=nons)
        numh = hh.makeAbsHist(numh, "y")
    else:
        if args.generator == "matrix_radish":
            h = input_tools.read_matrixRadish_hist(corr_file, "ptVgen")
        else:
            h = input_tools.read_dyturbo_hist(corr_file.split(":"), axis="ptVgen")

        charge_ax = minnloh.axes["chargeVgen"]
        #absy_ax = hist.axis.Regular(1, 0, 10, flow=True, name="absYVgen")
        # TODO: This is needed to match axes currently, but it's strictly correct
        absy_ax = hist.axis.Regular(1, 0, 5, underflow=False, overflow=True, name="absYVgen")
        mass_ax = hist.axis.Regular(1, 60, 120, flow=True, name="massVgen")
        vars_ax = hist.axis.Integer(0, 1, flow=False, name="vars")
        h5D = hist.Hist(mass_ax, absy_ax, *h.axes, charge_ax, vars_ax)
        h5D[...] = h.values(flow=True)[np.newaxis, np.newaxis,:,np.newaxis,np.newaxis]
        numh = h5D
    return numh

if args.proc == "z":
    if len(args.corr_files) != 1:
        raise ValueError("Only one file expected for Z")
    filesByProc = { "ZmumuPostVFP" : args.corr_files[0] }
elif args.proc == "w":
    if len(args.corr_files) != 2:
        raise ValueError("Requires two files for W (W+ and W-)")
    plus_idx = 0 if "Wp" in args.corr_files[0] else 1
    filesByProc = { "WplusmunuPostVFP" : args.corr_files[plus_idx],
        "WminusmunuPostVFP" : args.corr_files[0 if plus_idx else 1]}

minnloh = input_tools.read_all_and_scale(args.minnlo_file, list(filesByProc.keys()), ["nominal_gen"])[0]
# Get more stats in the correction
add_taus = True 
if add_taus:
    logging.info("Combining muon and tau decay samples for increased stats")
    from wremnants.datasets.datasetDict_v9 import Z_TAU_TO_LEP_RATIO,BR_TAUToMU
    taus = ["WplustaunuPostVFP", "WminustaunuPostVFP"] if args.proc == "w" else ["ZtautauPostVFP"]
    taush = input_tools.read_all_and_scale(args.minnlo_file, taus, ["nominal_gen"])[0]
    # Rescale taus to mu to effectively have more stats
    minnloh = 0.5*(minnloh + taush/(Z_TAU_TO_LEP_RATIO if args.proc == "z" else BR_TAUToMU))

numh = hh.sumHists([read_corr(procName, args.generator, corr_file) for procName, corr_file in filesByProc.items()])

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

logging.info("Correction binning is")
for ax in corrh.axes:
    logging.info(f"Axis {ax.name}: {ax.edges}")
logging.info(f"Average correction is {corrh.sum()/np.ones_like(corrh).sum()}")
logging.info(f"Wrote file {outfile}")
