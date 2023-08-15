import mplhep as hep
import matplotlib.pyplot as plt

import uproot
import argparse
import os
import numpy as np
import pandas as pd
import hist
import pdb

from utilities import logging, output_tools
from wremnants import plot_tools

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult root file")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["covariance"], choices=["correlation", "covariance"], help="Define which plots to make")
parser.add_argument("--normalize", action='store_true', help="Plot normalized distributions")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("-c", "--channels", type=str, nargs="+", choices=["plus", "minus", "all"], default=["plus", "minus", "all"], help="Select channel to plot")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("plotFitresult", 4 if args.debug else 3, False)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

rfile = uproot.open(args.infile)
if args.asimov:
    asimov = uproot.open(args.asimov)
else:
    asimov=None

input_subdir = args.infile.split("/")[-2]

cms_decor = "Preliminary"


def get_bin(name, var):
    name_split = name.split(var)
    if len(name_split) == 1:
        return -1
    else:
        return int(name_split[-1].split("_")[0])

def plot_matrix_poi(matrix="covariance_matrix_channelmu", base_process=None, axes=None, keys=None):
    hist_cov = matrix_poi(matrix, base_process=base_process, axes=axes, keys=keys)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot() 

    hep.hist2dplot(hist_cov)#, labels=(xlabels,ylabels))

    # calculate condition number
    cond = np.linalg.cond(hist_cov.values())
    logger.info(f"Condition number: {cond}")
    plt.text(0.2, 0.9, round(cond), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

    outfile = "covariance" if "covariance" in matrix else "correlation"

    if base_process is not None:
        outfile += f"_{base_process}"

    if keys is not None:
        outfile += "_" + "_".join(keys)

    if axes is not None:
        outfile += "_" + "_".join(axes)

    outfile += (f"_{args.postfix}_" if args.postfix else "_") + matrix.split("_")[-1].replace("channel","")

    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Values" : hist_cov.values()}, nround=2 if "correlation" in matrix else 10,
        analysis_meta_info=None,
        args=args,
    )

if "correlation" in args.plots:
    plot_matrix_poi("correlation_matrix_channelmu")
    plot_matrix_poi("correlation_matrix_channelpmaskedexp")
    plot_matrix_poi("correlation_matrix_channelpmaskedexpnorm")

if "covariance" in args.plots:
    plot_matrix_poi("covariance_matrix_channelmu")
    plot_matrix_poi("covariance_matrix_channelpmaskedexp")
    plot_matrix_poi("covariance_matrix_channelpmaskedexpnorm")

if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)