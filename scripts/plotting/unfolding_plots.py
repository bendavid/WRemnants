import mplhep as hep
import matplotlib.pyplot as plt
import uproot
import argparse
import os
import numpy as np
import matplotlib as mpl
import pandas as pd
import boost_histogram as bh
import hist
import pdb

from utilities import boostHistHelpers as hh, logging, input_tools, common, output_tools
from wremnants import plot_tools

hep.style.use(hep.style.ROOT)


parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult root file")
# parser.add_argument("--ratioToData", action='store_true', help="Use data as denominator in ratio")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=None, help="y range for ratio plot")
# parser.add_argument("--rebin", type=int, default=1, help="Rebin (for now must be an int)")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
# parser.add_argument("--xlim", type=float, nargs=2, help="min and max for x axis")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--noData", action='store_true', help="Don't plot data")
# parser.add_argument("--noFill", action='store_true', help="Don't fill stack")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--plots", type=str, nargs="+", default=["postfit"], choices=["prefit", "postfit"], help="Define which plots to make")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")
parser.add_argument("-c", "--channels", type=str, nargs="+", choices=["plus", "minus", "all"], default=["plus", "minus"], help="Select channel to plot")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("plotFitresult", 4 if args.debug else 3)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)

rfile = uproot.open(args.infile)

# reco bin settings
nbins_reco_charge = 2
nbins_reco_eta = 48

input_subdir = args.infile.split("/")[-2]

# gen bins
nbins_charge = 2
nbins_pt = 3
nbins_eta = 2

if input_subdir.startswith("WMass"):
    mode="wmass"
    base_process = "Wmunu"
    nbins_reco_pt = 29
    cm = mpl.colormaps["autumn"]
elif input_subdir.startswith("ZMassWLike"):
    mode="wlike"
    base_process = "Zmumu"
    nbins_reco_pt = 34
    cm = mpl.colormaps["winter"]
elif input_subdir.startswith("ZMassDilepton"):
    mode="dilepton"
    base_process = "Zmumu"
    nbins_reco_pt = 20
    cm = mpl.colormaps["winter"]

nbins_reco = nbins_reco_charge * nbins_reco_pt * nbins_reco_eta

cms_decor = "Preliminary" if not args.noData else "Simulation Preliminary"

binwnorm = 1.0

def get_bin(name, var):
    name_split = name.split(var)
    if len(name_split) == 1:
        return 0
    else:
        return int(name_split[-1].split("_")[0])

def getProcessBins(name, axes = ["qGen", "ptGen", "absEtaGen"]):
    res = {
        x: get_bin(name, x) if get_bin(name, x) else 0 for x in axes
    }

    return res

# labels
labelmap = {
    "Zmumu" : r"Z$\to\mu\mu$",
    "Zmumu_bkg" : r"Z$\to\mu\mu$ (bkg)",
    "Ztautau" : r"Z$\to\tau\tau$",
    "Wmunu_bkg" : r"W$^{\pm}\to\mu\nu$ (bkg)",
    "Wmunu" : r"W$^{\pm}\to\mu\nu$ ",
    "Wtaunu" : r"W$^{\pm}\to\tau\nu$",
    # "Top" : "green",
    # "Diboson" : "pink",
    # "QCD" : "grey",
    # "Fake" : "grey",
} 

def get_label(name):
    logger.debug(f"Get label for {name}")
    if name in labelmap.keys():
        return labelmap[name]

    if mode == "dilepton":
        res = getProcessBins(name, axes=["ptVGen"])
        idx = res["ptVGen"]
        label = r"Z$\to\mu\mu$"
        label += f"({idx})"
        return label
    else:
        res = getProcessBins(name)
        eta = res["absEtaGen"]
        pt = res["ptGen"]
        charge = res["qGen"]    

        if name.startswith("Wmunu"):
            label = r"W$^{+}\to\mu\nu$" if charge else r"W$^{-}\to\mu\nu$"
            label += f"({eta};{pt})"
            return label

        if name.startswith("Zmumu"):
            label = r"Z$\to\mu^{+}$" if charge else r"Z$\to\mu^{-}$"
            label += f"({eta};{pt})"
            return label

        logger.warning(f"No label found for {name}")
        return name

# colors
colormap= {
    "Zmumu" : "lightblue",
    "Zmumu_bkg" : "lightblue",
    "Ztautau" : "darkblue",
    "Wmunu_bkg" : "darkred",
    "Wmunu" : "darkred",
    "Wtaunu" : "orange",
    "Top" : "green",
    "Diboson" : "pink",
    "QCD" : "grey",
    "Fake" : "grey",
    "Other" : "grey",
}

def get_bin_number(pt, eta, charge):
    return eta + nbins_eta*pt + nbins_eta*nbins_pt*charge

def get_color(name):
    logger.debug(f"Get color for {name}")
    if name in colormap.keys():
        return colormap[name]
    if name.startswith(base_process):

        if mode == "dilepton":
            res = getProcessBins(name, axes=["ptVGen"])
            idx = res["ptVGen"]
            icol = (1+idx) / (nbins_reco_pt)

        else:
            res = getProcessBins(name)
            eta = res["absEtaGen"]
            pt = res["ptGen"]
            charge = res["qGen"]    

            icol = (1+eta + nbins_eta*pt + nbins_eta*nbins_pt*charge) / (nbins_eta*nbins_pt*nbins_charge)

        return cm(icol) 

    else:
        logger.warning(f"No color found for {name}")
        return "white"

def make_yields_df(hists, procs, signal=None, per_bin=False):
    logger.debug(f"Make yield df for {procs}")

    if per_bin:
        def sum_and_unc(h):
            variances = h.variances() if h.variances() is not None else h.values()
            return (h.values(), np.sqrt(variances))   
    else:
        def sum_and_unc(h):
            variances = h.variances() if h.variances() is not None else h.values()
            return (sum(h.values()), np.sqrt(sum(variances)))

    if per_bin:
        entries = [(i, v[0], v[1]) for i,v in enumerate(zip(*sum_and_unc(hists[0])))]
    else:
        if signal is not None:
            entries = [(signal, sum([ sum(v.values()) for k,v in zip(procs, hists) if signal in k]), np.sqrt(sum([ sum(v.variances()) for k,v in zip(procs, hists) if signal in k])))]
        else:
            entries = [(k, *sum_and_unc(v)) for k,v in zip(procs, hists)]

    return pd.DataFrame(entries, columns=["Process", "Yield", "Uncertainty"])



def plot(fittype, bins=(None, None), channel=None):
    logger.info(f"Make {fittype} plot"+(f" in channel {channel}" if channel else ""))

    procs = [p for p in filter(lambda x: x.startswith("expproc_") and x.endswith(f"_{fittype};1"), rfile.keys())]

    proc_sig = filter(lambda x: f"_{base_process}_qGen" in x, procs)
    proc_bkg = filter(lambda x: f"_{base_process}_qGen" not in x, procs)

    proc_bkg = [s for s in sorted(proc_bkg, key=lambda x: sum(rfile[x].to_hist().values()))]
    proc_sig = [s for s in sorted(proc_sig, reverse=True)]

    bin_lo = bins[0]
    bin_hi = bins[1]

    # pre computation
    if "obs" not in [c.replace(";1","") for c in  rfile.keys()]:
        logger.error(f"Shapes not found in fitresult file, run combine with --saveHists --computeHistErrors to get shapes.")
        return
    hist_data = rfile["obs"].to_hist()[bin_lo:bin_hi]
    nbins_data = len(*hist_data.axes.edges) - 1

    # last bin of exp histograms contain xnorm counts, remove for plotting by taking range from obs histogram
    if bin_lo is None:
        bin_lo = 0
    if bin_hi is None:
        bin_hi = nbins_data
    hist_pred = rfile[f"expfull_{fittype}"].to_hist()[bin_lo:bin_hi]

    #divide histogram by bin widths
    if mode == "dilepton":
        bins = np.array(common.ptV_binning)
        bin_widths = bins[1:] - bins[:-1]
    else: 
        bin_widths = np.ones(nbins_data)    

    hist_data /= bin_widths
    hist_pred /= bin_widths
    
    if args.ylim is None:
        ylim = (0, 1.1 * max(max(hist_data.values()), max(hist_pred.values())))
    else:
        ylim = args.ylim

    if args.rrange is None:
        rrange = [0.90,1.1] if fit_type=="prefit" else [0.95, 1.05]
    else:
        rrange = args.rrange

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_data, "Bin number", "Events/bin", ylim, "Data/Pred.", rrange)

    stack = [rfile[p].to_hist()[bin_lo:bin_hi]/bin_widths for p in procs]
    names = [k.replace("expproc_","").replace(f"_{fittype};1","") for k in procs]
    colors = [get_color(l) for l in names]
    labels = [get_label(l) for l in names]

    for i, p in enumerate(procs):
        logger.debug(f"Process: {p} | Label: {labels[i]} | Color {colors[i]} | Yield: {sum(stack[i].values())}")

    hep.histplot(
        stack,
        histtype="fill",
        color=colors,
        label=labels,
        stack=True,
        ax=ax1,
        binwnorm=binwnorm,
        zorder=1,
    )

    if not args.noData:
        hep.histplot(
            hist_data,
            yerr=True,
            histtype="errorbar",
            color="black",
            label="Data",
            ax=ax1,
            alpha=1.,
            binwnorm=binwnorm,
            zorder=2,
        )    

        hep.histplot(
            hh.divideHists(hist_data, hist_pred, cutoff=0.01, rel_unc=False),
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=False,
            linewidth=2,
            ax=ax2
        )

    # uncertainties
    # need to divide by bin width
    axis = hist_pred.axes[0].edges
    #binwidth =  axis[1:] - axis[:-1]
    nom = hist_pred.values() #/ binwidth
    std = np.sqrt(hist_pred.variances()) #/ binwidth

    hatchstyle = '///'
    ax1.fill_between(axis, 
            np.append(nom+std, (nom+std)[-1]), 
            np.append(nom-std, (nom-std)[-1]),
        step='post',facecolor="none", zorder=2, hatch=hatchstyle, edgecolor="k", linewidth=0.0, label="Uncertainty")

    ax2.fill_between(axis, 
            np.append((nom+std)/nom, ((nom+std)/nom)[-1]), 
            np.append((nom-std)/nom, ((nom-std)/nom)[-1]),
        step='post',facecolor="none", zorder=2, hatch=hatchstyle, edgecolor="k", linewidth=0.0)


    # plot_tools.addLegend(ax1, ncols=4, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{fittype}"
    outfile += "_"+input_subdir
    outfile += (f"_{args.postfix}" if args.postfix else "") 
    outfile += (f"_{channel}" if channel else "")

    plot_tools.save_pdf_and_png(outdir, outfile)

    # make yield tables
    stack = [s*bin_widths for s in stack]
    stack_yields = make_yields_df(stack, names)

    # summed up
    summed_yields = make_yields_df(stack, names, signal=base_process)
    if not args.noData:
        summed_yields.append(make_yields_df([hist_data*bin_widths], ["Data"]))

    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Stacked processes" : stack_yields, "Summed processes": summed_yields},#, "Unstacked processes" : unstacked_yields},
        analysis_meta_info=None,
        args=args,
    )

for fit_type in args.plots:

    if mode == "dilepton":
        channels = ["all"]
    else:
        channels = args.channels

    if "all" in channels:
        plot(fit_type)
    if "minus" in channels:
        plot("prefit", bins=(None,int(nbins_reco/2)), channel="minus")
    if "plus" in channels:
        plot("prefit", bins=(int(nbins_reco/2), nbins_reco), channel="plus")

if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)