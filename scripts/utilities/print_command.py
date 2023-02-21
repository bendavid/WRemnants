import lz4.frame
import pickle
import argparse
import os
import pathlib
from wremnants.datasets.datagroups import datagroups2016

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

def print_command_from_dict(res):
    meta_data = res["meta_data"] if "meta_data" in res else res["meta_info"]
    print(meta_data["command"])


datagroups = datagroups2016(args.infile)
if datagroups.rtfile:
    print_command_from_root(args.infile)
else:
    print_command_from_dict(datagroups.results)
