import argparse
import h5py

from utilities import logging
from utilities.io_tools import input_tools
import narf

import pdb

parser = argparse.ArgumentParser()
parser.add_argument("infiles", type=str, nargs="+", help="Input hdf5 files")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("-o", "--outfolder", type=str, default="", help="Output folder")
args = parser.parse_args()
logger = logging.setup_logger(__file__)

results = {}
for i, infile in enumerate(args.infiles):
    h5file = h5py.File(infile, "r")
    result = narf.ioutils.pickle_load_h5py(h5file["results"])

    for key, value in result.items():
        if key in results.keys():
            if key == "meta_info":
                results[f"meta_info_{i}"] = value
            else:
                raise NotImplementedError(f"The object with key {key} from file {infile} was already present in a previous file; No implementation to solve this conflict")
        else:
            results[key] = value

outfile = args.infiles[0]
if args.postfix:
    outfile = outfile.replace(".hdf5", f"_{args.postfix}.hdf5")

with h5py.File(f"{args.outfolder}/{outfile}", 'w') as f:
    narf.ioutils.pickle_dump_h5py("results", results, f)