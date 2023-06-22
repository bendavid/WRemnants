import argparse
import os
import mplhep as hep
import matplotlib.pyplot as plt

from wremnants import logging
from wremnants import plot_tools
from wremnants.datasets.datagroups2016 import make_datagroups_2016
from utilities import boostHistHelpers as hh, output_tools

import pdb

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histogrdams")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="test", help="Subfolder for output")
parser.add_argument("--procFilters", type=str, nargs="*", default="Zmumu", help="Filter to plot (default no filter, only specify if you want a subset")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--axes", type=str, nargs="+", default=["pt-ptGen","abs(eta)-absEtaGen"], help="Define for which axes the response matrix to be plotted")
parser.add_argument("-c", "--channels", type=str, nargs="+", choices=["plus", "minus", "all"], default=["all"], help="Select channel to plot")

args = parser.parse_args()

logger = logging.setup_logger("makeDataMCStackPlot", 4 if args.debug else 3)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

groups = make_datagroups_2016(args.infile, filterGroups=args.procFilters)
datasets = groups.getNames()
logger.info(f"Will plot datasets {datasets}")

groups.setGenAxes([]) # set gen axes empty to not integrate over when loading

groups.loadHistsForDatagroups("nominal", syst="nominal", procsToRead=datasets)
datagroups = groups.getDatagroups()

hist = datagroups["Zmumu"].hists["nominal"]

for channel in args.channels:
    select = {} if channel == "all" else {"charge" : -1.j if channel == "minus" else 1.j}

    for axes_string in args.axes:
        axes = axes_string.split("-")

        if axes[0].startswith("abs("):
            # mirror axis at half
            hist2d = hh.projectNoFlow(hist[select], [axes[0][4:-1], *axes[1:]])
            nbins = len(hist2d.axes.edges[0])-1
            values = hist2d.values()[:int(nbins/2)][::-1] + hist2d.values()[int(nbins/2):]
            xbins = hist2d.axes[0].edges[int(nbins/2):]
        else:
            hist2d = hh.projectNoFlow(hist[select], axes)
            values = hist2d.values()
            xbins = hist2d.axes[0].edges

        ybins = hist2d.axes[1].edges

        fig = plt.figure()#figsize=(8*width,8))
        ax = fig.add_subplot() 

        hep.hist2dplot(values, xbins=xbins, ybins=ybins)#, labels=(xlabels,ylabels))

        # ax.set_xticks(np.arange(len(xlabels))+0.5)
        # ax.set_yticks(np.arange(len(xlabels))+0.5)
        # ax.set_xticklabels(xlabels, rotation = 90)
        # ax.set_yticklabels(xlabels)

        outfile = "responce_matrix_"+"_".join([a.replace("(","").replace(")","") for a in axes])

        outfile += (f"_{channel}" if channel != "all" else "") 

        plot_tools.save_pdf_and_png(outdir, outfile)

        # plot_tools.write_index_and_log(outdir, outfile, 
        #     yield_tables={"Values" : cov_mat}, nround=2 if "correlation" in matrix else 10,
        #     analysis_meta_info=None,
        #     args=args,
        # )