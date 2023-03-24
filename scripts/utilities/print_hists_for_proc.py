import argparse
from utilities import input_tools

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="Input hdf5 file")
parser.add_argument("-s", "--sample", type=str, required=True, help="Sample name (e.g., ZmumuPostVFP)")

args = parser.parse_args()

names = input_tools.read_hist_names(args.infile, args.sample)
print("Valid names for process {args.sample} are:")
print(names)
