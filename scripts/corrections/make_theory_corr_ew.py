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
parser.add_argument("-i", "--input", type=str, default="w_z_gen_dists_RawPFMET.pkl.lz4", help="File containing EW hists")
parser.add_argument("--num", type=str, default="horace-nlo", help="Numerator")
parser.add_argument("--den", type=str, default="horace-lo-photos", help="Denominator")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/EWCorrections", help="Output path")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']

# file created with `python WRemnants/scripts/histmakers/w_z_gen_dists.py --skipAngularCoeffs --filter horace --ewHists`
with lz4.frame.open(args.input) as f:
    res = pickle.load(f)

for proc in procs:
    hnum = hh.normalize(res[f'{proc}_{args.num}']['output']['nominal_ew'])
    hden = hh.normalize(res[f'{proc}_{args.den}']['output']['nominal_ew'])
    corrh = hh.divideHists(hnum, hden)
    corrh = hh.smoothenTowardsOne(corrh)

    outfile = f"{args.outpath}/{proc}_{args.num}_over_{args.den}.pkl.lz4"

    with lz4.frame.open(outfile, "wb") as f:
        pickle.dump({
                proc : {
                    "ewratio" : corrh,
                },
                "meta_data" : output_tools.metaInfoDict(),
            }, f, protocol = pickle.HIGHEST_PROTOCOL)

    logging.info("Correction binning is")
    for ax in corrh.axes:
        logging.info(f"Axis {ax.name}: {ax.edges}")
    logging.info(f"Wrote file {outfile}")
