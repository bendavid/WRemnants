import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import itertools
import uproot
import argparse
import os
import numpy as np
import matplotlib as mpl
import pandas as pd
import hist
import h5py
from wremnants import histselections as sel

from utilities import boostHistHelpers as hh, logging, input_tools, common, differential, output_tools
from wremnants import plot_tools
from wremnants.datasets.datagroups import Datagroups
from wremnants.unfolding_tools import get_bin, getProcessBins, get_results, matrix_poi

import pdb

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="Output file of the analysis stage, containing ND boost histogrdams")
parser.add_argument("--fitresult",  type=str, help="Combine fitresult root file")
parser.add_argument("--asimov",  type=str, default=None, help="Optional combine fitresult root file from an asimov fit for comparison")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=None, help="y range for ratio plot")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--noData", action='store_true', help="Don't plot data")
parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["xsec", "uncertainties"], choices=["xsec", "uncertainties"], help="Define which plots to make")
parser.add_argument("--normalize", action='store_true', help="Plot normalized distributions")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("--plotSumPOIs", action='store_true', help="Plot xsecs from sum POI groups")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("unfolding_xsec", 4 if args.debug else 3)

cms_decor = "Preliminary"

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

groups = Datagroups(args.infile)

if groups.wmass:
    process = "Wenu" if groups.flavor == "e" else "Wmunu"
else:
    process = "Zee" if groups.flavor == "ee" else "Zmumu"

base_process = process[0]

gen_axes = groups.gen_axes
if base_process == "W":
    gen_axes.append("qGen")

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

def plot_xsec_unfolded(df, df_asimov=None, hreco=None, hxnorm=None, edges=None, bin_widths=None, channel=None, scale=1., normalize=False, process_label="V", axes=None):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    if normalize:
        yLabel="1/$\sigma$ d$\sigma("+process_label+")$"
    else:
        yLabel="d$\sigma ("+process_label+")$ [pb]"

    hist_xsec = hist.Hist(
        hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())

    if bin_widths is None:
        bin_widths = edges[1:] - edges[:-1]

    hist_xsec.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_total"].values/bin_widths)**2], axis=-1)

    if hreco:
        hreco_flat = hist.Hist(
            hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())
        hreco_flat.view(flow=False)[...] = np.stack([hreco.values().flatten()/bin_widths, hreco.variances().flatten()/bin_widths**2], axis=-1)

    if hxnorm:
        hxnorm_flat = hist.Hist(
            hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())
        hxnorm_flat.view(flow=False)[...] = np.stack([hxnorm.values().flatten()/bin_widths, hxnorm.variances().flatten()/bin_widths**2], axis=-1)

    unc_ratio = np.sqrt(hist_xsec.variances()) /hist_xsec.values() 

    if "err_stat" in df.keys():
        hist_xsec_stat = hist.Hist(
            hist.axis.Variable(edges, underflow=False, overflow=False), storage=hist.storage.Weight())
        hist_xsec_stat.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_stat"].values/bin_widths)**2], axis=-1)
        unc_ratio_stat = np.sqrt(hist_xsec_stat.variances()) /hist_xsec.values() 

    if data_asimov is not None:
        ha_xsec = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))
        ha_xsec.view(flow=False)[...] = df_asimov["value"].values/bin_widths

    # make plots
    if args.ylim is None:
        ylim = (0, 1.1 * max(hist_xsec.values() + np.sqrt(hist_xsec.variances())))
    else:
        ylim = args.ylim

    if args.rrange is None:
        rrange = [0.9,1.1]
    else:
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

    if hreco_flat:
        hep.histplot(
            hreco_flat,
            yerr=False,
            histtype="step",
            color="grey",
            label="Reconstructed data",
            ax=ax1,
            alpha=1.,
            zorder=2,
        )            
    if hxnorm_flat:
        hep.histplot(
            hxnorm_flat,
            yerr=False,
            histtype="step",
            color="blue",
            label="Theory model",
            ax=ax1,
            alpha=1.,
            zorder=2,
        )        

    centers = hist_xsec.axes.centers[0]

    ax2.bar(centers, height=2*unc_ratio, bottom=1-unc_ratio, width=edges[1:] - edges[:-1], color="silver", label="Total")
    if "err_stat" in df.keys():
        ax2.bar(centers, height=2*unc_ratio_stat, bottom=1-unc_ratio_stat, width=edges[1:] - edges[:-1], color="gold", label="Stat")

    ax2.plot([min(edges), max(edges)], [1,1], color="black", linestyle="-")

    # hep.histplot(
    #     hh.divideHists(hist_xsec_stat, hist_xsec, cutoff=0, rel_unc=True),
    #     histtype="errorbar",
    #     color="black",
    #     label="Data",
    #     yerr=True,
    #     linewidth=2,
    #     ax=ax2, capsize=2, elinewidth=0, markersize=0
    # )

    # hep.histplot(
    #     hh.divideHists(hist_xsec, hist_xsec, cutoff=0, rel_unc=True),
    #     histtype="errorbar",
    #     color="black",
    #     label="Data",
    #     yerr=True,
    #     linewidth=2,
    #     ax=ax2
    # )

    if data_asimov is not None:
        hep.histplot(
            ha_xsec,
            yerr=False,
            histtype="step",
            color="blue",
            label="Prefit model",
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

    if hreco_flat:
        hep.histplot(
            hh.divideHists(hreco_flat, hist_xsec, cutoff=0, rel_unc=True),
            yerr=False,
            histtype="step",
            color="grey",
            ax=ax2,
        )            
    if hxnorm_flat:
        hep.histplot(
            hh.divideHists(hxnorm_flat, hist_xsec, cutoff=0, rel_unc=True),
            yerr=False,
            histtype="step",
            color="blue",
            ax=ax2,
        )     

    plot_tools.addLegend(ax1, ncols=2, text_size=15*args.scaleleg)
    plot_tools.addLegend(ax2, ncols=2, text_size=15*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_xsec"

    if normalize:
        outfile += "_normalized"

    outfile += f"_{base_process}"
    outfile += "_"+"_".join(axes)
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    if data_asimov is not None:
        asimov_yields = make_yields_df([ha_xsec], ["Model"], per_bin=True)
        asimov_yields["Uncertainty"] *= 0 # artificially set uncertainty on model hard coded to 0
    data_yields = make_yields_df([hist_xsec], ["Data"], per_bin=True)
    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Data" : data_yields, "Model": asimov_yields} if df_asimov else {"Data" : data_yields},
        analysis_meta_info=None,
        args=args,
    )
    plt.close()

translate = {
    "QCDscalePtChargeMiNNLO": "QCDscale",
    "QCDscaleZPtChargeMiNNLO": "QCDscaleZ",
    "QCDscaleWPtChargeMiNNLO": "QCDscaleW",
    "QCDscaleZPtHelicityMiNNLO": "QCDscaleZ",
    "QCDscaleWPtHelicityMiNNLO": "QCDscaleW",
    "QCDscaleZPtChargeHelicityMiNNLO": "QCDscaleZ",
    "QCDscaleWPtChargeHelicityMiNNLO": "QCDscaleW",
    "resumNonpert": "resumNP",
    "resumTransition": "resumT",
    "binByBinStat": "BBlite",
    "CMS_recoil": "recoil",
    "CMS_background": "Bkg.",
    "NonpromptHighMT": "NPhighMT",
    "NonpromptLowMT": "NPlowMT",
    "massShiftZ": "massZ",
    "massShiftW": "massW",
}

def plot_uncertainties_unfolded(df, channel=None, edges=None, scale=1., normalize=False, relative_uncertainty=False, logy=False, process_label="", axes=None):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    # read nominal values and uncertainties from fit result and fill histograms
    logger.debug(f"Produce histograms")

    if normalize:
        yLabel="1/$\sigma$ d$\sigma("+process_label+")$"
    else:
        yLabel="d$\sigma ("+process_label+")$ [pb]"

    if relative_uncertainty:
        yLabel = "$\delta$ "+ yLabel
    else:
        yLabel = "$\Delta$ "+ yLabel
    
    #central values
    bin_widths = edges[1:] - edges[:-1]

    values = df["value"].values/bin_widths

    hist_xsec = hist.Hist(hist.axis.Variable(edges, underflow=False, overflow=False))

    errors = df["err_total"].values/bin_widths
    if relative_uncertainty:
        errors /= values
    
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
    uncertainties = make_yields_df([hist_xsec], ["Total"], per_bin=True, yield_only=True, percentage=True)

    sources =["err_stat"]
    sources += ["err_NonpromptHighMT", "err_NonpromptLowMT", "err_fakerate"]
    sources += list(sorted([s for s in filter(lambda x: x.startswith("err"), df.keys()) 
        if s not in ["err_stat", "err_total", "err_NonpromptHighMT", "err_NonpromptLowMT", "err_fakerate"] 
            and "eff_stat_" not in s and "eff_syst_" not in s]))    # only take eff grouped stat and syst

    NUM_COLORS = len(sources)-1
    cm = mpl.colormaps["gist_rainbow"]
    # add each source of uncertainty
    i=0
    for source in sources:

        name = source.replace("err_","").replace("muon_","")

        name = translate.get(name,name)

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
    plot_tools.redo_axis_ticks(ax1, "x", True)

    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_uncertainties"

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
        analysis_meta_info=None,
        args=args,
    )

    plt.close()

# store unfolded data
# outfile = h5py.File(f'{outdir}/fitresult.hdf5', 'w')

scale = 1 if args.normalize else args.lumi * 1000

poi_type = ["pmaskedexpnorm",] if args.normalize else ["pmaskedexp",]
if args.plotSumPOIs:
    poi_type += ["sumpoisnorm",] if args.normalize else ["sumpois",]

data = get_results(args.fitresult, poi_type=poi_type, scale=scale)
data_asimov = get_results(args.asimov, poi_type=poi_type, scale=scale)

# make all possible gen axes 
gen_axes_permutations = [gen_axes]
if args.plotSumPOIs:
    # Include all combinations of axes
    for n in range(1, len(gen_axes)):
        gen_axes_permutations += [k for k in itertools.combinations(gen_axes, n)]

all_axes = ["qGen", "ptGen",  "absEtaGen",  "ptVGen",  "absYVGen"]

for axes in gen_axes_permutations:
    data_group = data.loc[(sum([data[ax] != -1 for ax in axes]) == len(axes)) 
        & (np.product([data[ax] == -1 for ax in data.keys() if ax in all_axes and ax not in axes], axis=0))
        & (data["Name"].apply(lambda x: x[0] == base_process))]

    if len(data_group) == 0:
        logger.debug(f"No entries found with gen axes {axes}, next one!")
        continue

    if data_asimov:
        data_group_asimov = data_asimov.loc[(sum([data_asimov[ax] != -1 for ax in axes]) == len(axes)) 
            & (np.product([data_asimov[ax] == -1 for ax in data_asimov.keys() if ax in all_axes and ax not in axes], axis=0))
            & (data["Name"].apply(lambda x: x[0] == base_process))]

    if len(set(data_group["qGen"].values)) == 1 or (len(axes) == 1 and axes[0] == "qGen"):
        channels = ["all"]
    else:
        channels = ["plus", "minus"]

    logger.info(f"Make plots for process {base_process} and gen axes {axes}")

    for channel in channels:
        logger.info(f"Now at channel {channel}")

        if channel == "minus":
            data_channel = data_group.loc[data_group["qGen"]==0]
            channel_keys = {"qGen":-1}
            channel_axes = [a for a in axes if a != "qGen"]
            data_channel_asimov = data_group_asimov.loc[data_group_asimov["qGen"]==0] if data_asimov else None
            process_label = r"\mathrm{W}^{-}" if base_process == "W" else r"\mathrm{Z}"
        elif channel == "plus":
            data_channel = data_group.loc[data_group["qGen"]==1]
            channel_keys = {"qGen":1}
            channel_axes = [a for a in axes if a != "qGen"]
            data_channel_asimov = data_group_asimov.loc[data_group_asimov["qGen"]==1] if data_asimov else None
            process_label = r"\mathrm{W}^{+}" if base_process == "W" else r"\mathrm{Z}"
        else:
            process_label = r"\mathrm{W}" if base_process == "W" else r"\mathrm{Z}"
            channel_keys = None
            channel_axes = [*axes]
            data_channel = data_group
            data_channel_asimov = data_group_asimov if data_asimov else None

        if len(data_channel) == 0:
            logger.info(f"No entries found for channel {channel}, skip!")
            continue
        
        # find bin widths
        histo = sum([groups.results[m.name]["output"][args.baseName].get() for m in groups.groups[process].members 
            if channel=="all" and not m.name.startswith("Bkg") or channel in m.name])
        
        hproj = histo.project(*channel_axes)
        hproj = hh.scaleHist(hproj, 1./scale)

        hxnorm = sum([groups.results[m.name]["output"]["xnorm"].get() for m in groups.groups[process].members 
            if channel=="all" and not m.name.startswith("Bkg") or channel in m.name])
        hxproj = hxnorm.project(*channel_axes)
        hxproj = hh.scaleHist(hxproj, 1./scale)

        bins = np.product(hproj.axes.size)

        binwidths = None
        if len(channel_axes) == 1:
            if channel_axes[0] == "qGen":
                edges = np.array([-2,0,2])
            else:
                edges = np.array(histo.axes[channel_axes[0]].edges)
        elif len(channel_axes) == 2:
            xbins, ybins = [a.edges for a in histo.axes if a.name in channel_axes]
            xbinwidths = xbins[1:]-xbins[:-1]
            ybinwidths = ybins[1:]-ybins[:-1]
            binwidths = np.outer(xbinwidths, ybinwidths).flatten()
            edges = np.arange(0.5, len(binwidths)+1.5, 1.0)
        else:
            bins = np.product(hproj.axes.size)
            edges = np.arange(0.5, bins+1.5, 1.0)

        # sort values
        data_c = data_channel.sort_values(by=channel_axes)
        data_c_asimov = data_channel_asimov.sort_values(by=channel_axes) if data_asimov else None

        if "xsec" in args.plots:
            plot_xsec_unfolded(data_c, data_c_asimov, hproj, hxproj, edges=edges, bin_widths=binwidths, channel=channel, scale=scale, normalize=args.normalize, process_label = process_label, axes=channel_axes)

        if "uncertainties" in args.plots:
            # absolute uncertainty
            # plot_uncertainties_unfolded(data_c, channel=channel, scale=scale, normalize=args.normalize, process_label = process_label)
            # plot_uncertainties_unfolded(data_c, channel=channel, scale=scale, normalize=args.normalize, logy=True, process_label = process_label)            
            
            # relative uncertainty
            # plot_uncertainties_unfolded(data_c, channel=channel, scale=scale, normalize=args.normalize, relative_uncertainty=True, process_label = process_label)
            plot_uncertainties_unfolded(data_c, edges=edges, channel=channel, scale=scale, normalize=args.normalize, relative_uncertainty=True, logy=True, process_label = process_label, axes=channel_axes)


        # if not args.normalize:
        #     if set(gen_axes) != set(axes):
        #         continue

        #     # write out 1D distributions
        #     logger.info(f"Save measured differential cross secction distribution")
        #     hist_xsec = hist.Hist(
        #         hist.axis.Regular(bins=len(data_c), start=0.5, stop=len(data_c)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
        #     hist_xsec.view(flow=False)[...] = np.stack([data_c["value"].values, (data_c["err_total"].values)**2], axis=-1)

        #     # write out covariance as 2D hist
        #     hist_cov = matrix_poi(rfile, poi_type[0], base_process=base_process, axes=channel_axes, keys=channel_keys)
            
        #     histname = f"data_{base_process}"
        #     covname = f"covariance_matrix_{base_process}"
        #     if channel != "all":
        #         histname += f"_{channel}_"
        #         covname += f"_{channel}_"

        #     histname += "_".join(axes)
        #     covname += "_".join(axes)

        #     outfile.create_dataset(histname, data=hist_xsec)
        #     outfile.create_dataset(covname, data=hist_cov)

        # if "correlation" in args.plots:
        #     plot_matrix_poi(f"correlation_matrix_channel{poi_type}", base_process=base_process, axes=channel_axes, keys=channel_keys)
        # if "covariance" in args.plots:
        #     plot_matrix_poi(f"covariance_matrix_channel{poi_type}", base_process=base_process, axes=channel_axes, keys=channel_keys)

# outfile.close()

if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)