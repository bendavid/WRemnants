import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import hist
import numpy as np
import argparse

from utilities import logging, boostHistHelpers as hh
from utilities.styles import styles
from wremnants import plot_tools
from utilities.io_tools import output_tools, combinetf2_input

import pdb

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="hdf5 file from combinetf2")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./test", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--cmsDecor", default="Preliminary", type=str, choices=[None,"Preliminary", "Work in progress", "Internal"], help="CMS label")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9,1.1], help="y range for ratio plot")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--logy", action='store_true', help="Make the yscale logarithmic")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--noRatio", action='store_true', help="Don't make the ratio in the plot")
parser.add_argument("--noData", action='store_true', help="Don't plot the data")
parser.add_argument("--prefit", action='store_true', help="Make prefit plot, else postfit")
parser.add_argument("--selectionAxes", type=str, default=["charge"], help="List of axes where for each bin a seperate plot is created")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger(__file__, 4 if args.debug else 3)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)

fittype = "prefit" if args.prefit else "postfit"
ratio = not args.noRatio
data = not args.noData
fitresult = combinetf2_input.get_fitresult(args.infile)

meta = fitresult["meta"].get()
axes = meta["axes"]
procs = meta["procs"].astype(str)

labels = [styles.process_labels.get(p, p) for p in procs]
colors = [styles.process_colors.get(p, "grey") for p in procs]

hists_data_obs = combinetf2_input.load_histograms(fitresult, "hist_data_obs", axes)
hists_stack = combinetf2_input.load_histograms(fitresult, f"hist_{fittype}", axes)
hists_inclusive = combinetf2_input.load_histograms(fitresult, f"hist_{fittype}_inclusive", axes)

translate_selection = {
    "charge": {
        0 : "minus",
        1 : "plus"
    }
}

def make_plot(h_data, h_inclusive, h_stack, axis_name, suffix):

    if ratio:
        fig, ax1, ax2 = plot_tools.figureWithRatio(h_data, styles.xlabels.get(axis_name, "Bin number"), "Entries/bin", args.ylim, "Data/Pred.", args.rrange)
    else:
        fig, ax1 = plot_tools.figure(h_data, styles.xlabels.get(axis_name, "Bin number"), "Entries/bin", args.ylim)

    hep.histplot(
        h_stack,
        xerr=False,
        yerr=False,
        histtype="fill",
        color=colors,
        label=labels,
        stack=True,
        density=False,
        binwnorm=1.0,
        ax=ax1,
        zorder=1,
    )

    if data:
        hep.histplot(
            h_data,
            yerr=True,
            histtype="errorbar",
            color="black",
            label="Data",
            binwnorm=1.0,
            ax=ax1,
            alpha=1.,
            zorder=2,
        )    

    if ratio:

        hep.histplot(
            hh.divideHists(h_inclusive, h_inclusive, cutoff=1e-8, rel_unc=True, flow=False, by_ax_name=False),
            histtype="step",
            color="grey",
            alpha=0.5,
            yerr=True,
            ax=ax2,
            linewidth=2,
            flow='none',
        )

        if data:
            hep.histplot(
                hh.divideHists(h_data, h_inclusive, cutoff=0.01, rel_unc=False),
                histtype="errorbar",
                color="black",
                label="Data",
                yerr=False,
                linewidth=2,
                ax=ax2
            )

            # for uncertaity bands
            edges = h_inclusive.axes[0].edges

            # need to divide by bin width
            binwidth = edges[1:]-edges[:-1]
            if h_inclusive.storage_type != hist.storage.Weight:
                raise ValueError(f"Did not find uncertainties in {fittype} hist. Make sure you run combinetf with --computeHistErrors!")
            nom = h_inclusive.values() / binwidth
            std = np.sqrt(h_inclusive.variances()) / binwidth

            hatchstyle = '///'
            ax1.fill_between(edges, 
                    np.append(nom+std, (nom+std)[-1]), 
                    np.append(nom-std, (nom-std)[-1]),
                step='post',facecolor="none", zorder=2, hatch=hatchstyle, edgecolor="k", linewidth=0.0, label="Uncertainty")

            ax2.fill_between(edges, 
                    np.append((nom+std)/nom, ((nom+std)/nom)[-1]), 
                    np.append((nom-std)/nom, ((nom-std)/nom)[-1]),
                step='post',facecolor="none", zorder=2, hatch=hatchstyle, edgecolor="k", linewidth=0.0)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data=data)

    plot_tools.addLegend(ax1, ncols=2, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    to_join = [fittype, args.postfix, axis_name, suffix]
    outfile = "_".join(filter(lambda x: x, to_join))

    plot_tools.save_pdf_and_png(outdir, outfile)

    stack_yields = None
    unstacked_yields = None
    plot_tools.write_index_and_log(outdir, outfile, 
        # yield_tables={"Stacked processes" : stack_yields, "Unstacked processes" : unstacked_yields},
        analysis_meta_info={"AnalysisOutput" : meta["meta_info"]},
        args=args,
    )


# loop over channels (channels here correspond to the different histograms e.g. wmass & dilepton combined fit)
nchannels=len(hists_data_obs)
for i, (hist_data, hist_inclusive, hist_stack) in enumerate(zip(hists_data_obs, hists_inclusive, hists_stack)):
    logger.info(f"Now at channel {i}")
    suffix = f"channel{i}"  if nchannels>1 else ""

    axes_names = [a.name for a in axes[i]]

    # loop over all axes, make a plot for each axis
    for axis in [n for n in axes_names if n not in args.selectionAxes]:

        # make inclusive plot
        h_data = hist_data.project(axis)
        h_inclusive = hist_inclusive.project(axis)
        h_stack = [h.project(axis) for h in hist_stack]

        make_plot(h_data, h_inclusive, h_stack, axis, suffix=f"{suffix}")

        # loop over selections, make a plot for each bin
        for sa in [a for a in axes[i] if a.name in args.selectionAxes]:
            for idx in range(len(sa)):
                sel = translate_selection[sa.name][idx]
                logger.info(f"Make plot for axis {axis}, selection {sel}")

                h_data = hist_data[{sa.name: idx}].project(axis)
                h_inclusive = hist_inclusive[{sa.name: idx}].project(axis)
                h_stack = [h[{sa.name: idx}].project(axis) for h in hist_stack]

                make_plot(h_data, h_inclusive, h_stack, axis, suffix=f"{suffix}{sel}")
        
