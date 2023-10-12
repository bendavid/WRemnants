import mplhep as hep
import matplotlib.pyplot as plt
from matplotlib import cm

import uproot
import argparse
import os
import numpy as np
import pandas as pd
import hist
import pdb
from utilities import logging, input_tools, output_tools
from wremnants import plot_tools
from utilities.input_tools_combinetf import get_fitresult, load_covariance_pois, select_covariance_poi

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult .root or .hdf5 file")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./test", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--cmsDecor", default="Preliminary", type=str, choices=[None,"Preliminary", "Work in progress", "Internal"], help="Name to append to file name")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--flow", action='store_true', help="Show overflow/underflow pois")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["covariance"], choices=["correlation", "covariance"], help="Define which plots to make")
parser.add_argument("-c", "--channels", type=str, nargs="+", choices=["plus", "minus", "all"], default=["plus", "minus", "all"], help="Select channel to plot")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("unfolding_covariance", 4 if args.debug else 3, False)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

fitresult = get_fitresult(args.infile)
meta_info = input_tools.get_metadata(args.infile)

label_dict = {
    "qGen": "charge",
    "ptGen": r"$p_\mathrm{T}$",
    "absEtaGen": r"$|\eta|$",
    "ptVGen": r"$p_\mathrm{T}^{V}$",
    "absYVGen": r"$|Y^{V}|$",
}

selection_dict = {
    "Z": r"$\mathrm{Z}\rightarrow\mu\mu$",
    "W": r"$\mathrm{W}^{\pm}\rightarrow\mu\nu$",
    "W_qGen0": r"$\mathrm{W}^-\rightarrow\mu\nu$",
    "W_qGen1": r"$\mathrm{W}^+\rightarrow\mu\nu$",
}

def plot_matrix_poi(h2d, names, poi_type="mu", gen_axes=None, base_processes=[], selections=None, covariance=False, 
    xlabel=None, ylabel=None, cms_decor="Preliminary", flow=False
    ):

    hist_cov = select_covariance_pois(h2d, names, gen_axes=gen_axes, base_processes=base_processes, selections=selections, flow=flow)
    if not covariance:
        std_devs = np.sqrt(np.diag(hist_cov.values()))

        # Calculate the correlation matrix
        hist_cov[...] = hist_cov.values() / np.outer(std_devs, std_devs)

    if xlabel is None:
        if gen_axes is not None:
            xlabel = "gen "+ "-".join([label_dict.get(a, a) for a in gen_axes if a not in selections.keys()])+" bin"
        else:
            xlabel = "POI bin"
    if ylabel is None:
        ylabel = xlabel

    fig, ax = plot_tools.figure(hist_cov, xlabel=xlabel, ylabel=ylabel, cms_label=cms_decor, automatic_scale=False, width_scale=1.2)

    if covariance:
        cmin, cmax = None, None
    else:
        cmin, cmax = -1.0, 1.0

    # Create a mask to None out the upper triangle
    values = hist_cov.values().copy()
    mask = np.tril(np.ones_like(values, dtype=bool))
    values[~mask] = None

    hep.hist2dplot(values, ax=ax)#, cmap=cm.RdBu)#, clim=(cmin, cmax))#, labels=(xlabels,ylabels))

    # calculate condition number
    cond = np.linalg.cond(hist_cov.values())
    logger.info(f"Condition number: {cond}")
    plt.text(0.2, 0.9, round(cond), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

    key = base_process
    for k,v in selections.items():
        key += f"_{k}{v}"

    ax.text(1.0, 1.003, selection_dict.get(key, key), transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")

    outfile = "covariance" if covariance else "correlation"
    if base_process is not None:
        outfile += f"_{key}"
    if gen_axes is not None:
        outfile += "_" + "_".join([a for a in gen_axes if a not in selections.keys()])
    outfile += (f"_{args.postfix}_" if args.postfix else "_") + poi_type.split("_")[-1].replace("channel","")

    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Values" : hist_cov.values()}, nround=10 if covariance else 2,
        analysis_meta_info={args.infile : datagroups.getMetaInfo()},
        args=args,
    )

for plot_type in args.plots:
    
    covariance = False if plot_type == "correlation" else True

    for poi_type in (
        # "mu", 
        # "pmaskedexp", 
        # "pmaskedexpnorm", 
        "sumpois", 
        "sumpoisnorm",
        ):
        hist_cov, names = load_covariance_pois(fitresult, poi_type)

        for i, base_process in enumerate(map(lambda x: "W" if x.split("/")[-1].startswith("mw") else "Z", meta_info["args"]["inputFile"])):
            logger.debug(f"Now at {i}, base process {base_process}")

            # TODO, get axis information from meta data
            if base_process == "W":
                selections=({"qGen":0}, {"qGen":1})

                if poi_type.startswith("sum"):
                    axes_combinations = (("qGen","absEtaGen"), ("qGen","ptGen"))
                else:
                    axes_combinations = (("qGen","ptGen","absEtaGen"),)

            elif base_process == "Z":
                selections = ({}, )
                if poi_type.startswith("sum"):
                    axes_combinations = ("ptVGen", "absYVGen")
                else:
                    axes_combinations = (("ptVGen", "absYVGen"),)

            for axes in axes_combinations:
                if isinstance(axes, str):
                    axes = [axes]
                for selection in selections:
                    plot_matrix_poi(hist_cov, names, poi_type, cms_decor=args.cmsDecor, 
                        gen_axes=axes, selections=selection, base_processes=[base_process,], covariance=covariance, flow=args.flow)


if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)