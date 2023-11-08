import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt

import itertools
import argparse
import os
import numpy as np
import pandas as pd
import hist
import json

from utilities import boostHistHelpers as hh, logging
from utilities.styles.styles import nuisance_groupings as groupings
from wremnants import plot_tools, histselections as sel
from wremnants.datasets.datagroups import Datagroups
from utilities.io_tools import input_tools, output_tools
from utilities.io_tools.combinetf_input import get_fitresult, read_impacts_pois, select_pois

import pdb

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="Output file of the analysis stage, containing ND boost histogrdams")
parser.add_argument("--fitresult",  type=str, help="Combine fitresult root file")
parser.add_argument("--reference",  type=str, default=None, help="Optional combine fitresult root file from an reference fit for comparison")
parser.add_argument("--refName",  type=str, default="Reference model", help="Name for reference source")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./test", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--cmsDecor", default="Preliminary", type=str, choices=[None,"Preliminary", "Work in progress", "Internal"], help="CMS label")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9,1.1], help="y range for ratio plot")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--logy", action='store_true', help="Make the yscale logarithmic")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--noData", action='store_true', help="Don't plot data")
parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["xsec", "uncertainties"], choices=["xsec", "uncertainties", "ratio"], help="Define which plots to make")
parser.add_argument("--genFlow", action='store_true', help="Show overflow/underflow pois")
parser.add_argument("--poi", action='store_true', help="Plot signal strength parameters (mu)")
parser.add_argument("--noi", action='store_true', help="Plot nuisance of interest parameters (noi)")
parser.add_argument("--poiRef", action='store_true', help="Plot for reference input signal strength parameters (mu)")
parser.add_argument("--noiRef", action='store_true', help="Plot for reference input nuisance of interest parameters (noi)")
parser.add_argument("--normalize", action='store_true', help="Plot normalized distributions")
parser.add_argument("--absolute", action='store_true', help="Plot absolute uncertainties, otherwise relative")
parser.add_argument("--plotSumPOIs", action='store_true', help="Plot xsecs from sum POI groups")
parser.add_argument("--scaleXsec", type=float, default=1.0, help="Scale xsec predictions with this number")
parser.add_argument("--grouping", type=str, default=None, help="Select nuisances by a predefined grouping", choices=groupings.keys())
parser.add_argument("--genAxes", type=str, nargs="+", default=None, help="Gen axes used in unfolding")
parser.add_argument("-t","--translate", type=str, default=None, help="Specify .json file to translate labels")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("unfolding_xsec", 4 if args.debug else 3)

grouping = groupings[args.grouping] if args.grouping else None

translate_label = {}
if args.translate:
    with open(args.translate) as f:
        translate_label = json.load(f)    

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

meta_info = input_tools.get_metadata(args.fitresult)

groups = Datagroups(args.infile)

if groups.mode in ["wmass", "lowpu_w"]:
    process = "Wenu" if groups.flavor == "e" else "Wmunu"
else:
    process = "Zee" if groups.flavor == "ee" else "Zmumu"

base_process = process[0]

gen_axes = groups.gen_axes if args.genAxes is None else args.genAxes

groups.setNominalName(args.baseName)
groups.loadHistsForDatagroups(args.baseName, syst="", procsToRead=[process])

input_subdir = args.fitresult.split("/")[-2]

xlabels = {
    "ptGen" : r"$p_{T}^{\ell}$ [GeV]",
    "absEtaGen" : r"$|\eta^{\ell}|$",
    "ptVGen" : r"$p_{T}^{PROCESS}$ [GeV]",
    "absYVGen" : r"$|Y^{PROCESS}|$",
    "qGen" : r"$q^{PROCESS}$",
}

def get_xlabel(axis, process_label):
    p = process_label.replace("$","")
    l = xlabels.get(axis, axis)
    return l.replace("PROCESS", p)

def make_yields_df(hists, procs, signal=None, per_bin=False, yield_only=False, percentage=True):
    logger.debug(f"Make yield df for {procs}")

    if per_bin:
        def sum_and_unc(h,scale=100 if percentage else 1):
            return (h.values()*scale, np.sqrt(h.variances())*scale)   
    else:
        def sum_and_unc(h,scale=100 if percentage else 1):
            return (sum(h.values())*scale, np.sqrt(sum(h.variances())*scale))

    if per_bin:
        entries = [(i, v[0], v[1]) for i,v in enumerate(zip(*sum_and_unc(hists[0])))]
        index = "Bin"
    else:
        index = "Process"
        if signal is not None:
            entries = [(signal, sum([ sum(v.values()) for k,v in zip(procs, hists) if signal in k]), np.sqrt(sum([ sum(v.variances()) for k,v in zip(procs, hists) if signal in k])))]
        else:
            entries = [(k, *sum_and_unc(v)) for k,v in zip(procs, hists)]

    if yield_only:
        entries = [(e[0], e[1]) for e in entries]
        columns = [index, *procs]
    else:
        columns = [index, "Yield", "Uncertainty"]


    return pd.DataFrame(entries, columns=columns)

def plot_xsec_unfolded(df, edges, poi_type, df_reference=None, bin_widths=None, channel=None, scale=1., normalize=False, process_label="V", axes=None,
    hist_others=[], label_others=[], color_others=[]
):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    if poi_type == "mu":
        yLabel="$\mu\ ("+process_label+")$"
    elif poi_type == "nois":
        yLabel="$\mathrm{NOI} + 1\ ("+process_label+")$"
    elif normalize:
        yLabel="1/$\sigma$ d$\sigma("+process_label+")$"
    else:
        yLabel="d$\sigma ("+process_label+")$ [pb]"

    hist_xsec = hist.Hist(
        hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())

    if bin_widths is None:
        bin_widths = edges[1:] - edges[:-1]

    hist_xsec.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_total"].values/bin_widths)**2], axis=-1)

    unc_ratio = np.sqrt(hist_xsec.variances()) /hist_xsec.values() 

    if "err_stat" in df.keys():
        hist_xsec_stat = hist.Hist(
            hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())
        hist_xsec_stat.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_stat"].values/bin_widths)**2], axis=-1)
        unc_ratio_stat = np.sqrt(hist_xsec_stat.variances()) /hist_xsec.values() 

    if df_reference is not None:
        ha_xsec = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))
        ha_xsec.view(flow=False)[...] = df_reference["value"].values/bin_widths

    # make plots
    if args.ylim is None:
        ylim = (0, 1.1 * max(hist_xsec.values() + np.sqrt(hist_xsec.variances())))
    else:
        ylim = args.ylim

    rrange = args.rrange

    xlabel = "-".join([get_xlabel(a, process_label) for a in axes])
    if len(axes) >= 2:
        xlabel = xlabel.replace("[GeV]","")
        xlabel += " Bin"

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_xsec, xlabel, yLabel, ylim, "Pred./Data", rrange, width_scale=2)

    hep.histplot(
        hist_xsec,
        yerr=np.sqrt(hist_xsec.variances()),
        histtype="errorbar",
        color="black",
        label="Unfolded data",
        ax=ax1,
        alpha=1.,
        zorder=2,
    )    

    if args.genFlow:
        ax1.fill([0,18.5, 18.5, 0,0], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)
        ax1.fill([len(edges)-17.5, len(edges)+0.5, len(edges)+0.5, len(edges)-17.5, len(edges)-17.5], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)

        ax2.fill([0,18.5, 18.5, 0,0], [rrange[0],*rrange,rrange[1],rrange[0]], color="grey", alpha=0.3)
        ax2.fill([len(edges)-17.5, len(edges)+0.5, len(edges)+0.5, len(edges)-17.5, len(edges)-17.5], [rrange[0],*rrange,rrange[1],rrange[0]], color="grey", alpha=0.3)

    centers = hist_xsec.axes.centers[0]

    ax2.bar(centers, height=2*unc_ratio, bottom=1-unc_ratio, width=edges[1:] - edges[:-1], color="silver", label="Total")
    if "err_stat" in df.keys():
        ax2.bar(centers, height=2*unc_ratio_stat, bottom=1-unc_ratio_stat, width=edges[1:] - edges[:-1], color="gold", label="Stat")

    ax2.plot([min(edges), max(edges)], [1,1], color="black", linestyle="-")

    for h, l, c in zip(hist_others, label_others, color_others):
        h_flat = hist.Hist(
            hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())
        h_flat.view(flow=False)[...] = np.stack([h.values(flow=args.genFlow).flatten()/bin_widths, h.variances(flow=args.genFlow).flatten()/bin_widths**2], axis=-1)

        hep.histplot(
            h_flat,
            yerr=False,
            histtype="step",
            color=c,
            label=l,
            ax=ax1,
            alpha=1.,
            zorder=2,
        )            

        hep.histplot(
            hh.divideHists(h_flat, hist_xsec, cutoff=0, rel_unc=True),
            yerr=False,
            histtype="step",
            color=c,
            ax=ax2,
            zorder=2,
        )            

    if df_reference is not None:
        hep.histplot(
            ha_xsec,
            yerr=False,
            histtype="step",
            color="blue",
            label=args.refName,
            ax=ax1,
            alpha=1.,
            zorder=2,
        ) 

        hep.histplot(
            hh.divideHists(ha_xsec, hist_xsec, cutoff=0, rel_unc=True),
            yerr=False,
            histtype="step",
            color="blue",
            # label="Model",
            ax=ax2
        )

    plot_tools.addLegend(ax1, ncols=2, text_size=15*args.scaleleg)
    plot_tools.addLegend(ax2, ncols=2, text_size=15*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_xsec"

    if normalize:
        outfile += "_normalized"

    outfile += f"_{base_process}"
    outfile += "_"+"_".join(axes)
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    if df_reference is not None:
        reference_yields = make_yields_df([ha_xsec], ["Model"], per_bin=True)
        reference_yields["Uncertainty"] *= 0 # artificially set uncertainty on model hard coded to 0
    data_yields = make_yields_df([hist_xsec], ["Data"], per_bin=True)
    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Data" : data_yields, "Model": reference_yields} if df_reference is not None else {"Data" : data_yields},
        analysis_meta_info={args.infile : groups.getMetaInfo(), args.fitresult: meta_info},
        args=args,
    )
    plt.close()

def plot_uncertainties_unfolded(df, poi_type, channel=None, edges=None, scale=1., normalize=False, 
    logy=False, process_label="", axes=None, relative_uncertainty=False, percentage=True,
    error_threshold=0.001,   # only uncertainties are shown with a max error larger than this threshold
):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    # read nominal values and uncertainties from fit result and fill histograms
    logger.debug(f"Produce histograms")

    if poi_type == "mu":
        yLabel="$\mu\ ("+process_label+")$"
    elif poi_type == "nois":
        yLabel="$\mathrm{NOI}\ ("+process_label+")$"
    elif normalize:
        yLabel="1/$\sigma$ d$\sigma("+process_label+")$"
    else:
        yLabel="d$\sigma ("+process_label+")$ [pb]"

    if relative_uncertainty:
        yLabel = "$\delta$ "+ yLabel
        yLabel = yLabel.replace(" [pb]","")
        if percentage:
            yLabel += " [%]"
    else:
        yLabel = "$\Delta$ "+ yLabel
    
    if poi_type in ["mu", "nois"]:
        bin_widths = 1.0
    else:
        # central values
        bin_widths = edges[1:] - edges[:-1]

    errors = df["err_total"].values/bin_widths
    if relative_uncertainty:
        values = df["value"].values/bin_widths
        errors /= values
        if percentage:
            errors *= 100

    hist_xsec = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))
    hist_xsec.view(flow=False)[...] = errors

    # make plots
    if args.ylim is None:
        if logy:
            ylim = (max(errors)/10000., 1000 * max(errors))
        else:
            ylim = (0, 2 * max(errors))
    else:
        ylim = args.ylim

    xlabel = "-".join([get_xlabel(a, process_label) for a in axes])
    if len(axes) >= 2:
        xlabel = xlabel.replace("[GeV]","")
        xlabel += "Bin"

    fig, ax1 = plot_tools.figure(hist_xsec, xlabel, yLabel, ylim, logy=logy, width_scale=2)

    hep.histplot(
        hist_xsec,
        yerr=False,
        histtype="step",
        color="black",
        label="Total",
        ax=ax1,
        alpha=1.,
        zorder=2,
    )
    uncertainties = make_yields_df([hist_xsec], ["Total"], per_bin=True, yield_only=True, percentage=percentage)
    
    if args.genFlow:
        ax1.fill([0,18.5, 18.5, 0,0], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)
        ax1.fill([len(errors)-17.5, len(errors)+0.5, len(errors)+0.5, len(errors)-17.5, len(errors)-17.5], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)

    if grouping:
        sources = ["err_"+g for g in grouping]
    else:
        sources =["err_stat"]
        sources += list(sorted([s for s in filter(lambda x: x.startswith("err"), df.keys()) 
            if s.replace("err_","") not in ["stat", "total"] ])) # total and stat are added first

    NUM_COLORS = len(sources)-1
    cm = mpl.colormaps["gist_rainbow"]
    # add each source of uncertainty
    i=0
    for source in sources:

        name = source.replace("err_","")

        name = translate_label.get(name,name)

        if source =="err_stat":
            color = "grey"
        else:
            color = cm(1.*i/NUM_COLORS)
            i += 1

        if i%3 == 0:
            linestype = "-" 
        elif i%3 == 1:
            linestype = "--" 
        else:
            linestype = ":" 

        hist_unc = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))

        if source not in df:
            logger.warning(f"Source {source} not found in dataframe")
            continue

        errors = df[source].values/bin_widths
        if relative_uncertainty:
            errors /= values

        if max(errors) < error_threshold:
            logger.debug(f"Source {source} smaller than threshold of {error_threshold}")
            continue

        logger.debug(f"Plot source {source}")

        if percentage:
            errors *= 100
        
        hist_unc.view(flow=False)[...] = errors

        hep.histplot(
            hist_unc,
            yerr=False,
            histtype="step",
            color=color,
            linestyle=linestype,
            label=name,
            ax=ax1,
            alpha=1.,
            zorder=2,
        )

        unc_df = make_yields_df([hist_unc], [name], per_bin=True, yield_only=True, percentage=True)
        uncertainties[name] = unc_df[name]

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)

    plot_tools.addLegend(ax1, ncols=4, text_size=15*args.scaleleg*scale)

    if args.yscale:
        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim(ymin, ymax*args.yscale)

    if not logy:
        plot_tools.redo_axis_ticks(ax1, "y")
    plot_tools.redo_axis_ticks(ax1, "x", no_labels=len(axes) >= 2)

    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_uncertainties"

    if poi_type in ["mu","nois"]:
        outfile += f"_{poi_type}"
    if relative_uncertainty:
        outfile += "_relative"   
    if normalize:
        outfile += "_normalized"
    if logy:
        outfile += "_log"
    outfile += f"_{base_process}"
    outfile += "_"+"_".join(axes)
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Unfolded data uncertainty [%]": uncertainties},
        analysis_meta_info={"AnalysisOutput" : groups.getMetaInfo(), 'CardmakerOutput': meta_info},
        args=args,
    )

    plt.close()


def plot_uncertainties_ratio(df, df_ref, poi_type, poi_type_ref, channel=None, edges=None, scale=1., normalize=False, 
    logy=False, process_label="", axes=None, relative_uncertainty=False, 
):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    # read nominal values and uncertainties from fit result and fill histograms
    logger.debug(f"Produce histograms")

    yLabel = f"ref.({poi_type_ref}) / nom.({poi_type}) -1"

    if poi_type in ["mu", "nois"]:
        bin_widths = 1.0
    else:
        # central values
        bin_widths = edges[1:] - edges[:-1]

    errors = df["err_total"].values/bin_widths
    errors_ref = df_ref["err_total"].values/bin_widths
    if relative_uncertainty:
        values = df["value"].values/bin_widths
        values_ref = df_ref["value"].values/bin_widths
        errors /= values
        errors_ref /= values_ref

    xx = errors_ref/errors -1

    hist_ratio = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))
    hist_ratio.view(flow=False)[...] = xx

    # make plots
    if args.ylim is None:
        if logy:
            ylim = (max(xx)/10000., 1000 * max(xx))
        else:
            miny = min(xx)
            maxy = max(xx)
            rangey = maxy-miny
            ylim = (miny-rangey*0.1, maxy+rangey*0.7)
    else:
        ylim = args.ylim

    xlabel = "-".join([get_xlabel(a, process_label) for a in axes])
    if len(axes) >= 2:
        xlabel = xlabel.replace("[GeV]","")
        xlabel += "Bin"

    fig, ax1 = plot_tools.figure(hist_ratio, xlabel, yLabel, ylim, logy=logy, width_scale=2)

    ax1.plot([min(edges), max(edges)], [0.,0.], color="black", linestyle="-")

    hep.histplot(
        hist_ratio,
        yerr=False,
        histtype="step",
        color="black",
        label="Total",
        ax=ax1,
        alpha=1.,
        zorder=2,
    )
    uncertainties = make_yields_df([hist_ratio], ["Total"], per_bin=True, yield_only=True)


    if args.genFlow:
        ax1.fill([0,18.5, 18.5, 0,0], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)
        ax1.fill([len(errors)-17.5, len(errors)+0.5, len(errors)+0.5, len(errors)-17.5, len(errors)-17.5], [ylim[0],*ylim,ylim[1],ylim[0]], color="grey", alpha=0.3)

    if grouping:
        sources = ["err_"+g for g in grouping]
    else:
        sources =["err_stat"]
        sources += list(sorted([s for s in filter(lambda x: x.startswith("err"), df.keys()) 
            if s.replace("err_","") not in ["stat", "total"] ])) # total and stat are added first

    NUM_COLORS = len(sources)-1
    cm = mpl.colormaps["gist_rainbow"]
    # add each source of uncertainty
    i=0
    for source in sources:

        name = source.replace("err_","")

        name = translate_label.get(name,name)

        if source =="err_stat":
            color = "grey"
        else:
            color = cm(1.*i/NUM_COLORS)
            i += 1

        if i%3 == 0:
            linestype = "-" 
        elif i%3 == 1:
            linestype = "--" 
        else:
            linestype = ":" 

        hist_unc = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))

        if source not in df:
            logger.warning(f"Source {source} not found in dataframe")
            continue

        errors = df[source].values/bin_widths
        errors_ref = df_ref[source].values/bin_widths
        if relative_uncertainty:
            errors /= values
            errors_ref /= values_ref

        logger.debug(f"Plot source {source}")

        
        hist_unc.view(flow=False)[...] = errors_ref/errors-1

        hep.histplot(
            hist_unc,
            yerr=False,
            histtype="step",
            color=color,
            linestyle=linestype,
            label=name,
            ax=ax1,
            alpha=1.,
            zorder=2,
        )

        unc_df = make_yields_df([hist_unc], [name], per_bin=True, yield_only=True)
        uncertainties[name] = unc_df[name]

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)

    plot_tools.addLegend(ax1, ncols=4, text_size=15*args.scaleleg*scale)

    if args.yscale:
        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim(ymin, ymax*args.yscale)

    if not logy:
        plot_tools.redo_axis_ticks(ax1, "y")
    plot_tools.redo_axis_ticks(ax1, "x", no_labels=len(axes) >= 2)

    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_uncertainties_ratio"

    if poi_type in ["mu","nois"]:
        outfile += f"_{poi_type}"
    if relative_uncertainty:
        outfile += "_relative"   
    if normalize:
        outfile += "_normalized"
    if logy:
        outfile += "_log"
    outfile += f"_{base_process}"
    outfile += "_"+"_".join(axes)
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Unfolded data uncertainty [%]": uncertainties},
        analysis_meta_info={"AnalysisOutput" : groups.getMetaInfo(), 'CardmakerOutput': meta_info},
        args=args,
    )

    plt.close()

fitresult = get_fitresult(args.fitresult)
meta_info_ref=None
if args.reference:
    meta_info_ref = input_tools.get_metadata(args.reference)
    fitresult_reference = get_fitresult(args.reference)

def get_poi_types(poi, noi, meta):
    if poi or noi:
        poi_types = ["mu",] if poi else ["nois",]
        scale = 1./(meta["args"]["scaleNormXsecHistYields"]*meta["args"]["priorNormXsec"]) if noi else 1
    else:
        poi_types = ["pmaskedexpnorm",] if args.normalize else ["pmaskedexp",]
        if args.plotSumPOIs:
            poi_types += ["sumpoisnorm",] if args.normalize else ["sumpois",]
        scale = 1 if args.normalize else args.lumi * 1000
    return poi_types, scale

poi_types, scale = get_poi_types(args.poi, args.noi, meta_info)
poi_types_ref, scale_ref = get_poi_types(args.poiRef, args.noiRef, meta_info_ref)


    
for poi_type, poi_type_ref in zip(poi_types, poi_types_ref):
    data = read_impacts_pois(fitresult, poi_type=poi_type, scale=scale)
    data_ref = read_impacts_pois(fitresult_reference, poi_type=poi_type_ref, scale=scale_ref) if args.reference else None        

    if poi_type == "nois":
        data["value"] += 1.0
    if poi_type_ref == "nois":
        data_ref["value"] += 1.0

    if poi_type.startswith("sum"):
        # make all possible lower dimensional gen axes combinations; wmass only combinations including qGen
        gen_axes_permutations = [list(k) for n in range(1, len(gen_axes)) for k in itertools.combinations(gen_axes, n)]
    else:
        # gen_axes = ["ptGen", ]
        gen_axes_permutations = [gen_axes[:],]

    for axes in gen_axes_permutations:
        logger.info(f"Make plots for process {base_process} and gen axes {axes}")

        if groups.mode in ["wmass", "lowpu_w"]:
            axes.append("qGen")

        channels = ["plus", "minus"] if "qGen" in axes else ["all"]

        for channel in channels:
            logger.info(f"Now at channel {channel}")
            if channel == "minus":
                channel_sel = {"qGen":0}
                channel_axes = [a for a in axes if a != "qGen"]
                process_label = r"\mathrm{W}^{-}" if base_process == "W" else r"\mathrm{Z}"
            elif channel == "plus":
                channel_sel = {"qGen":1}
                channel_axes = [a for a in axes if a != "qGen"]
                process_label = r"\mathrm{W}^{+}" if base_process == "W" else r"\mathrm{Z}"
            else:
                channel_sel = {}
                process_label = r"\mathrm{W}" if base_process == "W" else r"\mathrm{Z}"
                channel_axes = axes[:]

            data_c = select_pois(data, axes, selections=channel_sel, base_processes=base_process, flow=args.genFlow)
            data_c_ref = select_pois(data_ref, axes, selections=channel_sel, base_processes=base_process, flow=args.genFlow) if args.reference else None

            if len(data_c) == 0:
                logger.info(f"No entries found for axes {axes} in channel {channel}, skip!")
                continue
            
            # find bin widths
            def get_histo(name):
                h = sum([groups.results[m.name]["output"][name].get() for m in groups.groups[process].members 
                    if not m.name.startswith("Bkg") and (base_process=="Z" or channel=="all" or channel in m.name)])
                h = h.project(*channel_axes)
                # for wlike the sample is randomly split in two based on reco charge
                this_scale = 2*scale if groups.mode == "wlike" else scale
                if "xnorm" in name:
                    this_scale /= args.scaleXsec
                h = hh.scaleHist(h, 1./this_scale)
                return h

            histo = get_histo(args.baseName)
            hxnorm = get_histo("xnorm")
            hMiNNLO = get_histo("xnorm_uncorr")

            if args.genFlow:
                edges = plot_tools.extendEdgesByFlow(histo)
            else:
                edges = histo.axes.edges

            binwidths = None
            if len(channel_axes) == 1:
                if channel_axes[0] == "qGen":
                    edges = np.array([-2,0,2])
                else:
                    edges = np.array(edges[0])
            elif len(channel_axes) == 2:
                xbins, ybins = edges
                xbinwidths = np.diff(xbins.flatten())
                ybinwidths = np.diff(ybins.flatten())
                binwidths = np.outer(xbinwidths, ybinwidths).flatten()
                edges = np.arange(0.5, len(binwidths)+1.5, 1.0)
            else:
                bins = np.product([len(e) for e in edges])
                edges = np.arange(0.5, bins+1.5, 1.0)

            if poi_type in ["mu", "nois"]:
                binwidths = None

            if "xsec" in args.plots:
                plot_xsec_unfolded(data_c, edges, poi_type, data_c_ref, bin_widths=binwidths, channel=channel, scale=scale, normalize=args.normalize, axes=channel_axes, 
                    process_label=process_label, 
                    #hist_others=[hxnorm, hMiNNLO], label_others=[r"MiNNLO $\times$ SCETlib+DYTurbo", "MiNNLO"], color_others=["blue", "red"]
                )

            if "uncertainties" in args.plots:
                plot_uncertainties_unfolded(data_c, poi_type, edges=edges, channel=channel, scale=scale, 
                    normalize=args.normalize, relative_uncertainty=not args.absolute, logy=args.logy, process_label = process_label, axes=channel_axes)
                if args.reference:
                    plot_uncertainties_unfolded(data_c_ref, poi_type_ref, edges=edges, channel=channel, scale=scale, 
                        normalize=args.normalize, relative_uncertainty=not args.absolute, logy=args.logy, process_label = process_label, axes=channel_axes)
            
            if "ratio" in args.plots:
                plot_uncertainties_ratio(data_c, data_c_ref, poi_type, poi_type_ref, edges=edges, channel=channel, scale=scale, 
                    normalize=args.normalize, relative_uncertainty=not args.absolute, logy=args.logy, process_label = process_label, axes=channel_axes)


if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)
