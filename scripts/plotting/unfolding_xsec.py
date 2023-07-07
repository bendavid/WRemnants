import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import uproot
import argparse
import os
import numpy as np
import matplotlib as mpl
import pandas as pd
import hist
import pdb

from utilities import boostHistHelpers as hh, logging, input_tools, common, common_lowPU, differential, output_tools
from wremnants import plot_tools
from wremnants.unfolding_tools import get_bin, getProcessBins

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult root file")
parser.add_argument("--asimov",  type=str, default=None, help="Optional combine fitresult root file from an asimov fit for comparison")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=None, help="y range for ratio plot")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--noData", action='store_true', help="Don't plot data")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["xsec", "uncertainties"], 
    choices=["xsec", "uncertainties", "correlation", "covariance"], help="Define which plots to make")
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

binwnorm = 1.0

# labels
labelmap = {
    "Zmumu" : r"Z$\to\mu\mu$",
    "Zmumu_bkg" : r"Z$\to\mu\mu$ (bkg)",
    "Ztautau" : r"Z$\to\tau\tau$",
    "Wmunu_bkg" : r"W$^{\pm}\to\mu\nu$ (bkg)",
    "Wmunu" : r"W$^{\pm}\to\mu\nu$ ",
    "Wtaunu" : r"W$^{\pm}\to\tau\nu$",
} 

def get_label(name, gen_axes=["qGen", "ptGen", "absEtaGen", "ptVGen", "absYVGen"]):
    logger.debug(f"Get label for {name}")
    if name in labelmap.keys():
        return labelmap[name]

    if dilepton:
        res = getProcessBins(name, axes=gen_axes)
        idx = res[gen_axes]

        label = r"Z"
        label += f"({idx})"
        return label

    else:
        res = getProcessBins(name)
        eta = res["absEtaGen"]
        pt = res["ptGen"]
        charge = res["qGen"]    

        if name.startswith("Wmunu"):
            label = r"W$^{+}$" if charge else r"W$^{-}$"
            label += f"({eta};{pt})"
            return label

        if name.startswith("Zmumu"):
            label = r"Z$^{+}$" if charge else r"Z$^{-}$"
            label += f"({eta};{pt})"
            return label

        logger.warning(f"No label found for {name}")
        return name

def get_bin_widths(df, dilepton=False, gen_axes=None):
    nPtVBins = len(set(df["ptVGen"].values))
    nQBins = len(set(df["qGen"].values))

    if dilepton:
        if "ptVGen" in gen_axes and "absYVGen" in gen_axes:
            bins_y = np.array(common.ptV_binning)
            bin_widths_y = bins_y[1:] - bins_y[:-1]
            bin_widths_x = np.ones(int(len(df)/len(bin_widths_y)))
            bin_widths = np.tensordot(bin_widths_y, bin_widths_x, axes=0).flatten()
        elif "ptVGen" in gen_axes:
            if nPtVBins == len(common.ptV_binning)-1:
                bins = np.array(common.ptV_binning)
            elif nPtVBins == len(common_lowPU.ptV_binning)-1:
                bins = np.array(common_lowPU.ptV_binning)
            bin_widths = bins[1:] - bins[:-1]
        elif "absYVGen" in gen_axes:
            bin_widths = np.ones(int(len(df)))
    else:
        nEtaBins = len(set(df["absEtaGen"].values))
        if nEtaBins == 18:
            bins_y = np.array(differential.eta_binning)
        else:
            bins_y = np.linspace(0, 2.4, nEtaBins+1)
        bin_widths_y = bins_y[1:] - bins_y[:-1]

        nPtBins = len(set(df["ptGen"].values))

        base_process = df["Name"].apply(lambda x: x.split("_")[0]).values[0]
        if base_process == "Zmumu":
            bin_width_x = 34/nPtBins
        else:
            bin_width_x = 30/nPtBins

        bin_widths_x = bin_width_x * np.ones(int(len(df)/len(bin_widths_y)))

        bin_widths = np.tensordot(bin_widths_x,bin_widths_y, axes=0).flatten()

    return np.tile(bin_widths, nQBins) 

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

def get_results(rtfile, poi_type, scale=1.0, group=True, uncertainties=None):
    results = []

    # pois = input_tools.getPOInames(rfile, poi_type=poi_type)

    fitresult = rtfile["fitresults"]

    histname = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"

    impacts = rtfile[histname].to_hist()

    # process names
    names = [k for k in impacts.axes[0]]

    logger.debug(f"Load ucertainties")
    # pick uncertainties
    if uncertainties is None:
        uncertainties = {f"err_{k}": impacts.values()[:,i] for i, k in enumerate(impacts.axes[1])}
    else:
        uncertainties = {f"err_{k}": impacts.values()[:,i] for i, k in enumerate(impacts.axes[1]) if k in uncertainties}

    # measured central value
    centrals = [fitresult[n].array()[0] for n in names]

    # total uncertainties
    totals = [fitresult[n+"_err"].array()[0] for n in names]

    df = pd.DataFrame({"Name":names, "value":centrals, "err_total":totals, **uncertainties})

    if scale != 1:
        df["value"] /= scale
        df["err_total"] /= scale
        for u in uncertainties.keys():
            df[u] /= scale

    # try to decode the name string into bin number
    for axis in ["qGen", "ptGen", "absEtaGen", "ptVGen", "absYVGen"]:
        df[axis] = df["Name"].apply(lambda x, a=axis: get_bin(x, a))

    df = df.sort_values(["qGen", "ptGen", "absEtaGen", "ptVGen", "absYVGen"], ignore_index=True)
    return df

def plot_xsec_unfolded(df, data_asimov=None, channel=None, poi_type="mu", scale=1., normalize=False,process_label="", dilepton=False):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    # read nominal values and uncertainties from fit result and fill histograms
    logger.debug(f"Read fitresults")
    logger.debug(f"Produce histograms")

    if data_asimov is not None:
        logger.debug(f"Read fitresults for asimov")
        if channel == "minus":
            df_asimov = data_asimov.loc[data_asimov["qGen"]==0]
        elif channel == "plus":
            df_asimov = data_asimov.loc[data_asimov["qGen"]==1]
        else:
            df_asimov = data_asimov

    if normalize:
        yLabel="1/$\sigma$ d$\sigma$("+process_label+")"
    else:
        yLabel="d$\sigma$("+process_label+") [pb]"

    bin_widths = get_bin_widths(df, dilepton, gen_axes)

    hist_xsec = hist.Hist(
        hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
    hist_xsec_stat = hist.Hist(
        hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())

    hist_xsec.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_total"].values/bin_widths)**2], axis=-1)
    hist_xsec_stat.view(flow=False)[...] = np.stack([df["value"].values/bin_widths, (df["err_stat"].values/bin_widths)**2], axis=-1)

    if data_asimov is not None:
        ha_xsec = hist.Hist(
            hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False))

        ha_xsec.view(flow=False)[...] = df_asimov["value"].values/bin_widths

    # make plots
    if args.ylim is None:
        ylim = (0, 1.1 * max((hist_xsec.values() + np.sqrt(hist_xsec.variances()))))
    else:
        ylim = args.ylim

    if args.rrange is None:
        rrange = [0.9,1.1]
    else:
        rrange = args.rrange

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_xsec, "Bin number", yLabel, ylim, "Pred./Data", rrange, width_scale=2)

    hep.histplot(
        hist_xsec,
        yerr=True,
        histtype="errorbar",
        color="black",
        label="Data",
        ax=ax1,
        alpha=1.,
        binwnorm=binwnorm,
        zorder=2,
    )    


    unc_ratio = np.sqrt(hist_xsec.variances()) /hist_xsec.values() 
    unc_ratio_stat = np.sqrt(hist_xsec_stat.variances()) /hist_xsec.values() 

    centers = hist_xsec.axes.centers[0]

    ax2.bar(centers, height=2*unc_ratio, bottom=1-unc_ratio, width=1, color="silver")
    ax2.bar(centers, height=2*unc_ratio_stat, bottom=1-unc_ratio_stat, width=1, color="gold")

    ax2.plot([0, len(df)+1], [1,1], color="black", linestyle="-")

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
            label="Model",
            ax=ax1,
            alpha=1.,
            binwnorm=binwnorm,
            zorder=2,
        ) 

        hep.histplot(
            hh.divideHists(ha_xsec, hist_xsec, cutoff=0, rel_unc=True),
            histtype="step",
            color="blue",
            label="Model",
            yerr=False,
            ax=ax2
        )

    plot_tools.addLegend(ax1, ncols=4, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{input_subdir}_unfolded_xsec"
    if poi_type=="mu":
        outfile += "_mu"
    elif normalize:
        outfile += "_normalized"

    outfile += f"_{base_process}"
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    if data_asimov is not None:
        asimov_yields = make_yields_df([ha_xsec], ["Model"], per_bin=True)
        asimov_yields["Uncertainty"] *= 0 # artificially set uncertainty on model hard coded to 0
    data_yields = make_yields_df([hist_xsec], ["Data"], per_bin=True)
    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Data" : data_yields, "Model": asimov_yields} if asimov else {"Data" : data_yields},
        analysis_meta_info=None,
        args=args,
    )
    plt.close()

translate = {
    "QCDscalePtChargeMiNNLO": "QCD scale",
    "QCDscaleZPtChargeMiNNLO": "QCD scale (Z)",
    "QCDscaleWPtChargeMiNNLO": "QCD scale (W)",
    "QCDscaleZPtHelicityMiNNLO": "QCD scale (Z)",
    "QCDscaleWPtHelicityMiNNLO": "QCD scale (W)",
    "resumNonpert": "resum. NP",
    "resumTransition": "resum. T",
    "binByBinStat": "BB lite",
    "CMS_recoil": "recoil",
}

def plot_uncertainties_unfolded(df, channel=None, poi_type="mu", scale=1., normalize=False, relative_uncertainty=False, logy=False, process_label="", dilepton=False):
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    # read nominal values and uncertainties from fit result and fill histograms
    logger.debug(f"Read fitresults")
    logger.debug(f"Produce histograms")

    if normalize:
        yLabel="1/$\sigma$ d$\sigma$("+process_label+")"
    else:
        yLabel="d$\sigma$("+process_label+") [pb]"

    if relative_uncertainty:
        yLabel = "$\delta$ "+ yLabel
    else:
        yLabel = "$\Delta$ "+ yLabel
    
    bin_widths = get_bin_widths(df, dilepton, gen_axes)

    #central values
    values = df["value"].values/bin_widths

    hist_xsec = hist.Hist(hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False))

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

    fig, ax1 = plot_tools.figure(hist_xsec, "Bin number", yLabel, ylim, logy=logy, width_scale=2)

    hep.histplot(
        hist_xsec,
        yerr=False,
        histtype="step",
        color="black",
        label="Total",
        ax=ax1,
        alpha=1.,
        binwnorm=binwnorm,
        zorder=2,
    )
    uncertainties = make_yields_df([hist_xsec], ["Total"], per_bin=True, yield_only=True, percentage=True)

    sources =["err_stat"]
    sources += [s for s in filter(lambda x: x.startswith("err"), df.keys()) 
        if s not in ["err_stat", "err_total"] 
            and "eff_stat_" not in s and "eff_syst_" not in s]    # only take eff grouped stat and syst

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

        hist_unc = hist.Hist(hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False))

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
            binwnorm=binwnorm,
            zorder=2,
        )

        unc_df = make_yields_df([hist_unc], [name], per_bin=True, yield_only=True, percentage=True)
        uncertainties[name] = unc_df[name]

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)

    plot_tools.addLegend(ax1, ncols=4, text_size=18*args.scaleleg*scale)

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

    if poi_type=="mu":
        outfile += "_mu"
    elif normalize:
        outfile += "_normalized"

    if logy:
        outfile += "_log"

    outfile += f"_{base_process}"
    outfile += (f"_{channel}" if channel else "")

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Unfolded data uncertainty [%]": uncertainties},
        analysis_meta_info=None,
        args=args,
    )

    plt.close()

# # store unfolded data
# outfile = uproot.recreate(f"{outdir}/unfolded_data.root")

poi_types = ["pmaskedexp", "pmaskedexpnorm"]

if args.normalize:
    poi_type = "pmaskedexpnorm"
else:
    poi_type = "pmaskedexp"


if any([key in args.plots for key in ["xsec", "uncertainties"]]):

    scale = 1 if args.normalize else args.lumi * 1000

    data = get_results(rfile, poi_type=poi_type, scale=scale)

    if asimov:
        data_asimov = get_results(asimov, poi_type=poi_type, scale=scale)
    else:
        data_asimov = None

    for dilepton, gen_axes in (
        (True, ["ptVGen"]),
        (True, ["ptVGen", "absYVGen"]),
        (False, ["qGen", "ptGen", "absEtaGen"]),
    ):

        data_group = data.loc[sum([data[ax] != -1 for ax in gen_axes]) == len(gen_axes)]

        if len(data_group) == 0:
            logger.debug(f"No entries found with gen axes {gen_axes}, next one!")
            continue

        if asimov:
            data_group_asimov = data_asimov.loc[sum([data_asimov[ax] != -1 for ax in gen_axes]) == len(gen_axes)]

        if len(set(data_group["qGen"].values)) == 1:
            channels = ["all"]
        else:
            channels = args.channels
        
        base_process = [x for x in set(data_group["Name"].apply(lambda x: x.split("_")[0]))]
        if len(base_process) != 1:
            logger.warning(f"Mutliple base processes are found: {base_process}, take first one")
        else:
            base_process = base_process[0]

        logger.info(f"Make plots for process {base_process} and gen axes {gen_axes}")

        for channel in channels:
            logger.info(f"Now at channel {channel}")

            if channel == "minus":
                data_channel = data_group.loc[data_group["qGen"]==0]
                channel_keys = ["qGen0"]
                channel_axes = [a for a in gen_axes if a != "qGen"]
                data_channel_asimov = data_group_asimov.loc[data_group_asimov["qGen"]==0] if asimov else None
                process_label = r"W$^{-}\to\mu\nu$" if base_process == "Wmunu" else r"Z$\to\mu^{-}$"
            elif channel == "plus":
                data_channel = data_group.loc[data_group["qGen"]==1]
                channel_keys = ["qGen1"]
                channel_axes = [a for a in gen_axes if a != "qGen"]
                data_channel_asimov = data_group_asimov.loc[data_group_asimov["qGen"]==1] if asimov else None
                process_label = r"W$^{+}\to\mu\nu$" if base_process == "Wmunu" else r"Z$\to\mu^{+}$"
            else:
                process_label = r"W$\to\mu\nu$" if base_process == "Wmunu" else r"Z$\to\mu\mu$"
                channel_keys = None
                channel_axes = gen_axes
                data_channel = data_group
                data_channel_asimov = data_group_asimov if asimov else None

            if len(data_channel) == 0:
                logger.info(f"No entries found for channel {channel}, skip!")
                continue
            
            # sort values
            data_c = data_channel.sort_values(by=channel_axes)
            data_c_asimov = data_channel_asimov.sort_values(by=channel_axes) if asimov else None

            # if not args.normalize:
            #     # write out 1D distributions
            #     logger.info(f"Save measured differential cross secction distribution")
            #     hist_xsec = hist.Hist(
            #         hist.axis.Regular(bins=len(data_c), start=0.5, stop=len(data_c)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
            #     hist_xsec.view(flow=False)[...] = np.stack([data_c["value"].values, (data_c["err_total"].values)**2], axis=-1)

            #     histname = f"data_{base_process}"
            #     covname = f"covariance_matrix_{base_process}"
            #     if channel != "all":
            #         histname += f"_{channel}"
            #         covname += f"_{channel}"

            #     outfile[histname] = hist_xsec
                
            #     # write out covariance as 2D hist
            #     hist_cov = matrix_poi(f"covariance_matrix_channel{poi_type}", base_process=base_process, axes=channel_axes, keys=channel_keys)
                
            #     outfile[covname] = hist_cov

            if "xsec" in args.plots:
                plot_xsec_unfolded(data_c, data_c_asimov, channel=channel, poi_type=poi_type, scale=scale, normalize=args.normalize, process_label = process_label, dilepton=dilepton)

            # if "correlation" in args.plots:
            #     plot_matrix_poi(f"correlation_matrix_channel{poi_type}", base_process=base_process, axes=channel_axes, keys=channel_keys)
            # if "covariance" in args.plots:
            #     plot_matrix_poi(f"covariance_matrix_channel{poi_type}", base_process=base_process, axes=channel_axes, keys=channel_keys)

            if "uncertainties" in args.plots:
                # absolute uncertainty
                # plot_uncertainties_unfolded(data_c, channel=channel, poi_type=poi_type, scale=scale, normalize=args.normalize, process_label = process_label, dilepton=dilepton)
                # plot_uncertainties_unfolded(data_c, channel=channel, poi_type=poi_type, scale=scale, normalize=args.normalize, logy=True, process_label = process_label, dilepton=dilepton)            
                
                # relative uncertainty
                # plot_uncertainties_unfolded(data_c, channel=channel, poi_type=poi_type, scale=scale, normalize=args.normalize, relative_uncertainty=True, process_label = process_label, dilepton=dilepton)
                plot_uncertainties_unfolded(data_c, channel=channel, poi_type=poi_type, scale=scale, normalize=args.normalize, relative_uncertainty=True, logy=True, process_label = process_label, dilepton=dilepton)

# if "correlation" in args.plots:
#     plot_matrix_poi("correlation_matrix_channelmu")
#     plot_matrix_poi("correlation_matrix_channelpmaskedexp")
#     plot_matrix_poi("correlation_matrix_channelpmaskedexpnorm")

# if "covariance" in args.plots:
#     plot_matrix_poi("covariance_matrix_channelmu")
#     plot_matrix_poi("covariance_matrix_channelpmaskedexp")
#     plot_matrix_poi("covariance_matrix_channelpmaskedexpnorm")

if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)