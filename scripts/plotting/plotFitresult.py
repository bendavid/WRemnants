import mplhep as hep
import uproot
import argparse
import os
import numpy as np
import matplotlib as mpl
import pandas as pd

import pdb

from utilities import boostHistHelpers as hh, logging
from wremnants import plot_tools

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histogrdams")
# parser.add_argument("--ratioToData", action='store_true', help="Use data as denominator in ratio")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="test", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9, 1.1], help="y range for ratio plot")
# parser.add_argument("--rebin", type=int, default=1, help="Rebin (for now must be an int)")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
# parser.add_argument("--xlim", type=float, nargs=2, help="min and max for x axis")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--noData", action='store_true', help="Don't plot data")
# parser.add_argument("--noFill", action='store_true', help="Don't fill stack")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")
parser.add_argument("--fittypes", type=str, nargs="+", default=["prefit", "postfit"], help="Make prefit or postfit plots, or both")

args = parser.parse_args()

logger = logging.setup_logger("plotFitresult", 4 if args.debug else 3, True)

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

combine_result = uproot.open(args.infile)

# settings
nbins_reco_charge = 2
nbins_reco_pt = 34 #29
nbins_reco_eta = 48

nbins_reco = nbins_reco_charge * nbins_reco_pt * nbins_reco_eta

# gen bins
nbins_charge = 1#2
nbins_pt = 2 #29
nbins_eta = 3 #48

cms_decor = "Preliminary" if not args.noData else "Simulation Preliminary"

lumi=16.8

binwnorm = 1.0

def getProcessPtEtaCharge(name):
    if name.startswith("Zmumu"):
        eta, pt = name.split("_")[1:]
        charge = 0
    else:  
        charge, eta, pt = name.split("_")[1:]
        charge = int(charge.replace("qGen",""))

    eta = int(eta.replace("etaGen",""))
    pt = int(pt.replace("ptGen",""))
    return pt, eta, charge

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
    elif name.startswith("Wmunu"):
        pt, eta, charge = getProcessPtEtaCharge(name)

        label = r"W$^{+}\to\mu\nu$" if charge else r"W$^{-}\to\mu\nu$"
        label += f"({eta};{pt})"

        return label
    elif name.startswith("Zmumu"):
        pt, eta, charge = getProcessPtEtaCharge(name)

        label = r"Z$\to\mu\mu$" #if charge else r"W$^{-}\to\mu\nu$"
        label += f"({eta};{pt})"

        return label
    else:
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
}

cm_w = mpl.colormaps["autumn"]
cm_z = mpl.colormaps["winter"]

def get_color(name):
    logger.debug(f"Get color for {name}")
    if name in colormap.keys():
        return colormap[name]
    elif name.startswith("Wmunu"):
        pt, eta, charge = getProcessPtEtaCharge(name)

        icol = (eta + nbins_eta*pt + nbins_eta*nbins_pt*charge) / (nbins_eta*nbins_pt*nbins_charge-1)

        return cm_w(icol)
    elif name.startswith("Zmumu"):
        pt, eta, charge = getProcessPtEtaCharge(name)

        icol = (eta + nbins_eta*pt + nbins_eta*nbins_pt*charge) / (nbins_eta*nbins_pt*nbins_charge-1)

        return cm_z(icol)        

    else:
        logger.warning(f"No color found for {name}")
        return "white"

def make_yields_df(hists, procs, signal=None):
    logger.debug(f"Make yield df for {procs}")

    def sum_and_unc(h):
        return (sum(h.values()), np.sqrt(sum(h.variances())))

    if signal is not None:
        entries = [(signal, sum([ sum(v.values()) for k,v in zip(procs, hists) if signal in k]), np.sqrt(sum([ sum(v.variances()) for k,v in zip(procs, hists) if signal in k])))]
    else:
        entries = [(k, *sum_and_unc(v)) for k,v in zip(procs, hists)]

    return pd.DataFrame(entries, columns=["Process", "Yield", "Uncertainty"])


def plot(fittype, bins=(None, None), channel=None):
    logger.info(f"Make {fittype} plot"+(f" in channel {channel}" if channel else ""))

    bin_lo = bins[0]
    bin_hi = bins[1]


    # pre computation
    hist_data = combine_result["obs"].to_hist()[bin_lo:bin_hi]
    hist_pred = combine_result[f"expfull_{fittype}"].to_hist()[bin_lo:bin_hi]

    if args.ylim is None:
        ylim = (0, 1.1 * max(max(hist_data.values()), max(hist_pred.values())))
    else:
        ylim = args.ylim

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_data, "Bin number", "Events/bin", ylim, "Data/Pred.", args.rrange)

    procs = [p for p in filter(lambda x: x.startswith("expproc_") and x.endswith(f"_{fittype};1"), combine_result.keys())]

    proc_sig = filter(lambda x: "_Wmunu_qGen" in x, procs)
    proc_bkg = filter(lambda x: "_Wmunu_qGen" not in x, procs)

    proc_bkg = [s for s in sorted(proc_bkg, key=lambda x: sum(combine_result[x].to_hist().values()))]
    proc_sig = [s for s in sorted(proc_sig, reverse=True)]

    procs = proc_bkg + proc_sig

    stack = [combine_result[p].to_hist()[bin_lo:bin_hi] for p in procs]
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


    extra_text=None, 
    plot_tools.addLegend(ax1, ncols=4, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{fittype}" + (f"_{args.postfix}" if args.postfix else "") + (f"_{channel}" if channel else "")
    plot_tools.save_pdf_and_png(outdir, outfile)
    stack_yields = make_yields_df(stack, names)
    signal_yields = make_yields_df(stack, names, signal="Wmunu")
    unstacked_yields = make_yields_df([hist_data], ["Data"])
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Stacked processes" : stack_yields, "Signal process": signal_yields, "Unstacked processes" : unstacked_yields},
        analysis_meta_info=None,
        args=args,
    )


for fittype in args.fittypes:
    plot(fittype)
    plot(fittype, bins=(None,int(nbins_reco/2)), channel="minus")
    plot(fittype, bins=(int(nbins_reco/2),None), channel="plus")
