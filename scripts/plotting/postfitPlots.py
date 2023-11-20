import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import hist
import numpy as np
import argparse
import pandas as pd

from utilities import logging, boostHistHelpers as hh
from utilities.styles import styles
from wremnants import plot_tools, histselections as sel
from utilities.io_tools import output_tools, combinetf_input, combinetf2_input

import pdb

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="hdf5 file from combinetf2 or root file from combinetf1")
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

if args.infile.endswith(".hdf5"):
    fitresult = combinetf2_input.get_fitresult(args.infile)
    combinetf2 = True
elif args.infile.endswith(".root"):
    fitresult = combinetf_input.get_fitresult(args.infile)
    combinetf2 = False
else:
    raise IOError("Unknown format of input file")

translate_selection = {
    "charge": {
        0 : "minus",
        1 : "plus"
    }
}

def make_plot(h_data, h_inclusive, h_stack, axis_name, colors=None, labels=None, suffix="", chi2=None, meta=None):

    if ratio:
        fig, ax1, ax2 = plot_tools.figureWithRatio(h_data, styles.xlabels.get(axis_name.replace("_","-"), "Bin number"), "Entries/bin", args.ylim, "Data/Pred.", args.rrange)
    else:
        fig, ax1 = plot_tools.figure(h_data, styles.xlabels.get(axis_name.replace("_","-"), "Bin number"), "Entries/bin", args.ylim)

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
        flow='none',
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
            flow='none',
        )    

    if ratio:

        hep.histplot(
            hh.divideHists(h_inclusive, h_inclusive, cutoff=1e-8, rel_unc=True, flow=False, by_ax_name=False),
            histtype="step",
            color="grey",
            alpha=0.5,
            yerr=False,
            ax=ax2,
            linewidth=2,
            flow='none',
        )

        if data:
            hep.histplot(
                hh.divideHists(h_data, h_inclusive, cutoff=0.01, rel_unc=True),
                histtype="errorbar",
                color="black",
                label="Data",
                yerr=True,
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

    if chi2 is not None:
        plt.text(0.05, 0.94, f"$\chi^2/ndf = {round(chi2[0],1)}/{chi2[1]}$", horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes,
            fontsize=20*args.scaleleg*scale)

    plot_tools.redo_axis_ticks(ax1, "x")
    plot_tools.redo_axis_ticks(ax2, "x")

    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data=data)

    plot_tools.addLegend(ax1, ncols=2, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    to_join = [fittype, args.postfix, axis_name, suffix]
    outfile = "_".join(filter(lambda x: x, to_join))

    plot_tools.save_pdf_and_png(outdir, outfile)

    stack_yields = None
    unstacked_yields = None
    if meta is not None:
        kwargs=dict(analysis_meta_info={"AnalysisOutput" : meta["meta_info"]})
    else:
        kwargs=dict()

    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={
            "Stacked processes" : pd.DataFrame([(k, sum(h.values()), sum(h.variances())**0.5) for k,h in zip(labels, h_stack)], columns=["Process", "Yield", "Uncertainty"]), 
            "Unstacked processes" : pd.DataFrame([(k, sum(h.values()), sum(h.variances())**0.5) for k,h in zip(["Data", "Inclusive"], [h_data, h_inclusive])], columns=["Process", "Yield", "Uncertainty"])},
        args=args, **kwargs
    )

if combinetf2:
    meta = fitresult["meta"].get()
    procs = meta["procs"].astype(str)

    labels = [styles.process_labels.get(p, p) for p in procs]
    colors = [styles.process_colors.get(p, "grey") for p in procs]

    chi2=None
    if f"chi2_{fittype}" in fitresult:
        chi2 = fitresult[f"chi2_{fittype}"], fitresult[f"ndf"]

    for channel, axes in meta["channel_axes"].items():

        hist_data = fitresult["hist_data_obs"][channel].get()
        hist_inclusive = fitresult[f"hist_{fittype}_inclusive"][channel].get()
        hist_stack = fitresult[f"hist_{fittype}"][channel].get()
        hist_stack = [hist_stack[{"processes" : p}] for p in procs]

        # make unrolled 1D histograms
        if len(axes) > 1:
            hist_data = sel.unrolledHist(hist_data, binwnorm=1, obs=None)
            hist_inclusive = sel.unrolledHist(hist_inclusive, binwnorm=1, obs=None)
            hist_stack = [sel.unrolledHist(h, binwnorm=1, obs=None) for h in hist_stack]

            axis_name = "_".join([a.name for a in axes])
        else:
            axis_name = axes[0].name

        make_plot(hist_data, hist_inclusive, hist_stack, axis_name, colors=colors, labels=labels, suffix=f"{channel}", chi2=chi2, meta=meta)

        # loop over selections, make a plot for each bin
        for sa in [a for a in axes if a.name in args.selectionAxes]:
            for idx in range(len(sa)):
                sel = translate_selection[sa.name][idx]
                logger.info(f"Make plot for axes {axis_name}, selection {sel}")

                h_data = hist_data[{sa.name: idx}]
                h_inclusive = hist_inclusive[{sa.name: idx}]
                h_stack = [h[{sa.name: idx}] for h in hist_stack]

                make_plot(h_data, h_inclusive, h_stack, axis_name, colors=colors, labels=labels, suffix=f"{channel}_{sel}", chi2=chi2, meta=meta)
else:
    # combinetf1    
    procs = [k.replace("expproc_","").replace(f"_{fittype};1", "") for k in fitresult.keys() if fittype in k and k.startswith("expproc_") and "hybrid" not in k]

    labels = [styles.process_labels.get(p, p) for p in procs]
    colors = [styles.process_colors.get(p, "grey") for p in procs]

    logger.info(f"Found processes {procs} in fitresult")

    axes_names = [part for part in args.infile.split("/")[-2].split("_") if part in ["pt", "eta", "charge", "ptll", "yll", "mll"]]
    axis_name = "_".join(axes_names)

    # get axes from the directory name
    analysis = args.infile.split("/")[-2].split("_")[0]
    if analysis=="ZMassDilepton":
        all_axes = {
            "mll": hist.axis.Regular(60, 60., 120., name = "mll", overflow=False, underflow=False),
            "yll": hist.axis.Regular(20, -2.5, 2.5, name = "yll", overflow=False, underflow=False),
            "ptll": hist.axis.Variable([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 23, 27, 32, 40, 54, 100], name = "ptll", underflow=False, overflow=False),
        }
    elif analysis=="ZMassWLike":
        all_axes = {
            "pt": hist.axis.Regular(30, 26, 56, name = "pt", overflow=False, underflow=False),
            "eta": hist.axis.Regular(48, -2.4, 2.4, name = "eta", overflow=False, underflow=False),
            "charge": common.axis_charge
        }
    elif analysis=="WMass":
        all_axes = {
            "pt": hist.axis.Regular(34, 26, 60, name = "pt", overflow=False, underflow=False),
            "eta": hist.axis.Regular(48, -2.4, 2.4, name = "eta", overflow=False, underflow=False),
            "charge": common.axis_charge
        }
    axes = [all_axes[a] for a in axes_names]
    shape = [len(a) for a in axes]

    hist_data = fitresult["obs;1"].to_hist()
    values = np.reshape(hist_data.values(), shape)
    hist_data = hist.Hist(*axes, storage=hist.storage.Weight(), data=np.stack((values, values), axis=-1))  

    hist_inclusive = fitresult[f"expfull_{fittype};1"].to_hist()
    hist_inclusive = hist.Hist(*axes, storage=hist.storage.Weight(), 
        data=np.stack((np.reshape(hist_inclusive.values(), shape), np.reshape(hist_inclusive.variances(), shape)), axis=-1))  
    hist_stack = [fitresult[f"expproc_{p}_{fittype};1"].to_hist() for p in procs]
    hist_stack = [hist.Hist(*axes, storage=hist.storage.Weight(), 
        data=np.stack((np.reshape(h.values(), shape), np.reshape(h.variances(), shape)), axis=-1)) for h in hist_stack]

    if len(axes) > 1:
        # make unrolled 1D histograms
        hist_data = sel.unrolledHist(hist_data, binwnorm=1, obs=None)
        hist_inclusive = sel.unrolledHist(hist_inclusive, binwnorm=1, obs=None)
        hist_stack = [sel.unrolledHist(h, binwnorm=1, obs=None) for h in hist_stack]

        axis_name = "_".join([a.name for a in axes])
    else:
        axis_name = axes[0].name

    make_plot(hist_data, hist_inclusive, hist_stack, axis_name, labels=labels, colors=colors)
