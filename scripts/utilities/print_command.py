import lz4.frame
import pickle
import argparse
import os
import pathlib
from utilities import input_tools

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help=".pkl.lz4 from with meta_info")

args = parser.parse_args()

if not os.path.isfile(args.infile):
    raise ValueError(f"{args.infile} is not a valid file!")

exts = pathlib.Path(args.infile).suffixes

def print_command_from_root(rtfile_name):
    import ROOT
    rtfile = ROOT.TFile.Open(rtfile_name)
    command = rtfile.Get("meta_info/command")
    print(command.GetTitle())

def print_command_from_dict(infile):
    meta_data = input_tools.get_metadata(infile)
    print(meta_data["command"])

if args.infile.endswith(".root"):
    print_command_from_root(args.infile)
else:
    print_command_from_dict(args.infile)
