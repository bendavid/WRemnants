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

from utilities import boostHistHelpers as hh, logging, input_tools
from wremnants import plot_tools

hep.style.use(hep.style.ROOT)


parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult root file")
parser.add_argument("--asimov",  type=str, default=None, help="Optional combine fitresult root file from an asimov fit for comparison")
# parser.add_argument("--ratioToData", action='store_true', help="Use data as denominator in ratio")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./", help="Subfolder for output")
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
parser.add_argument("--plots", type=str, nargs="+", default="xsec", choices=["xsec", "prefit", "postfit", "correlation", "covariance"], help="Define which plots to make")
parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity used in the fit, needed to get the absolute cross section")

args = parser.parse_args()

logger = logging.setup_logger("plotFitresult", 4 if args.debug else 3, False)

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

rfile = uproot.open(args.infile)
if args.asimov:
    asimov = uproot.open(args.asimov)
else:
    asimov=None

# reco bin settings
nbins_reco_charge = 2
nbins_reco_eta = 48

if any([x.startswith("WMass") for x in args.infile.split("/")]):
    base_process = "Wmunu"
    nbins_reco_pt = 29
    cm = mpl.colormaps["autumn"]
else:
    base_process = "Zmumu"
    nbins_reco_pt = 34

    cm = mpl.colormaps["winter"]

nbins_reco = nbins_reco_charge * nbins_reco_pt * nbins_reco_eta

# gen bins
nbins_charge = 2#2
nbins_pt = 3 #29
nbins_eta = 2 #48

cms_decor = "Preliminary" if not args.noData else "Simulation Preliminary"

binwnorm = 1.0

def getProcessPtEtaCharge(name):
    charge, eta, pt = name.split("_")[1:4]
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

        label = r"Z$\to\mu^{+}$" if charge else r"Z$\to\mu^{-}$"
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
    "Other" : "grey",
}

def get_bin_number(pt, eta, charge):
    return eta + nbins_eta*pt + nbins_eta*nbins_pt*charge

def get_color(name):
    logger.debug(f"Get color for {name}")
    if name in colormap.keys():
        return colormap[name]
    elif name.startswith(base_process):
        pt, eta, charge = getProcessPtEtaCharge(name)

        icol = (1+eta + nbins_eta*pt + nbins_eta*nbins_pt*charge) / (nbins_eta*nbins_pt*nbins_charge)

        return cm(icol) 

    else:
        logger.warning(f"No color found for {name}")
        return "white"

def make_yields_df(hists, procs, signal=None, per_bin=False):
    logger.debug(f"Make yield df for {procs}")

    if per_bin:
        def sum_and_unc(h):
            return (h.values(), np.sqrt(h.variances()))   
    else:
        def sum_and_unc(h):
            return (sum(h.values()), np.sqrt(sum(h.variances())))

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
    hist_pred = rfile[f"expfull_{fittype}"].to_hist()[bin_lo:bin_hi]

    if args.ylim is None:
        ylim = (0, 1.1 * max(max(hist_data.values()), max(hist_pred.values())))
    else:
        ylim = args.ylim

    if args.rrange is None:
        rrange = [0.90,1.1]
    else:
        rrange = args.rrange

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_data, "Bin number", "Events/bin", ylim, "Data/Pred.", rrange)

    stack = [rfile[p].to_hist()[bin_lo:bin_hi] for p in procs]
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


    plot_tools.addLegend(ax1, ncols=4, text_size=20*args.scaleleg)
    plot_tools.fix_axes(ax1, ax2, yscale=args.yscale)

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=cms_decor, data=not args.noData)

    outfile = f"{fittype}" + (f"_{args.postfix}" if args.postfix else "") + (f"_{channel}" if channel else "")
    plot_tools.save_pdf_and_png(outdir, outfile)
    stack_yields = make_yields_df(stack, names)
    signal_yields = make_yields_df(stack, names, signal=base_process)
    unstacked_yields = make_yields_df([hist_data], ["Data"])
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Stacked processes" : stack_yields, "Signal process": signal_yields, "Unstacked processes" : unstacked_yields},
        analysis_meta_info=None,
        args=args,
    )

def plot_matrix_poi(matrix="covariance_matrix_channelmu"):

    if matrix not in [c.replace(";1","") for c in rfile.keys()]:
        logger.error(f"Histogram {matrix} was not found in the fit results file!")
        return

    hist2d = rfile[matrix].to_hist()

    # select signal parameters
    key = matrix.split("channel")[-1]
    xentries = [(i, hist2d.axes[0][i]) for i in range(len(hist2d.axes[0])) if hist2d.axes[0][i].endswith(key)]
    xentries = sorted(xentries, key=lambda x: x[1], reverse=True)

    # make matrix between POIs only
    cov_mat = np.zeros((len(xentries), len(xentries)))
    for i, ia in xentries:
        for j, ja in xentries:
            cov_mat[i][j] = hist2d[i,j]

    xlabels = [get_label(x[1]) for x in xentries]

    fig = plt.figure()#figsize=(8*width,8))
    ax = fig.add_subplot() 

    hep.hist2dplot(cov_mat)#, labels=(xlabels,ylabels))
    

    ax.set_xticks(np.arange(len(xlabels))+0.5)
    ax.set_yticks(np.arange(len(xlabels))+0.5)
    ax.set_xticklabels(xlabels, rotation = 90)
    ax.set_yticklabels(xlabels)

    outfile = "covariance" if "covariance" in matrix else "correlation"

    outfile += (f"_{args.postfix}_" if args.postfix else "_") + matrix.split("_")[-1].replace("channel","")
        
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Values" : cov_mat}, nround=2 if "correlation" in matrix else 10,
        analysis_meta_info=None,
        args=args,
    )

def get_results(rootfile, pois, scale=1.0):
    results = []
    for poi in pois:
        impacts, labels, norm = input_tools.readImpacts(rootfile, True, add_total=True, POI=poi, normalize=False)
        res = {l: i/scale for i,l in zip(impacts, labels) if l in ["Total", "stat"]}

        res["norm"] = norm/scale

        charge, eta, pt = getProcessPtEtaCharge(poi)
        res.update({"pt":pt, "eta":eta, "charge": charge})

        res["bin"] = get_bin_number(pt, eta, charge)

        results.append(res)

    df = pd.DataFrame(results)
    df = df.sort_values("bin")
    return df

def plot_xsec_unfolded(bins=(None, None), channel=None, poi_type="mu"):
    normalize = poi_type=="pmaskedexpnorm"
    logger.info(f"Make "+("normalized " if normalize else "")+"unfoled xsec plot"+(f" in channel {channel}" if channel else ""))

    if normalize:
        yLabel="1/$\sigma$ d$\sigma$"
        scale = 1
    else:
        yLabel="d$\sigma$ [pb]"
        scale= args.lumi * 1000

    # read nominal values and uncertainties from fit result and fill histograms
    pois = input_tools.getPOInames(rfile, poi_type=poi_type)

    df = get_results(rfile, pois, scale=scale)

    hist_xsec = hist.Hist(
        hist.axis.Regular(bins=len(df), start=0, stop=max(df["bin"].values)+1, underflow=False, overflow=False), storage=hist.storage.Weight())
    hist_xsec_stat = hist.Hist(
        hist.axis.Regular(bins=len(df), start=0, stop=max(df["bin"].values)+1, underflow=False, overflow=False), storage=hist.storage.Weight())

    hist_xsec.view(flow=False)[...] = np.stack([df["norm"].values, (df["Total"].values)**2], axis=-1)
    hist_xsec_stat.view(flow=False)[...] = np.stack([df["norm"].values, (df["stat"].values)**2], axis=-1)

    if asimov:
        df_asimov = get_results(asimov, pois, scale=scale)

        ha_xsec = hist.Hist(
            hist.axis.Regular(bins=len(df), start=0, stop=max(df["bin"].values)+1, underflow=False, overflow=False))

        ha_xsec.view(flow=False)[...] = df_asimov["norm"].values

    # make plots
    if args.ylim is None:
        ylim = (0, 1.1 * max((df["norm"]+df["Total"]).values))
    else:
        ylim = args.ylim

    if args.rrange is None:
        rrange = [0.99,1.01] if normalize else [0.97,1.03]
    else:
        rrange = args.rrange

    fig, ax1, ax2 = plot_tools.figureWithRatio(hist_xsec, "Bin number", yLabel, ylim, "Data/Pred.", rrange)

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

    hep.histplot(
        hh.divideHists(hist_xsec, hist_xsec, cutoff=0.01, rel_unc=True),
        histtype="errorbar",
        color="black",
        label="Data",
        yerr=True,
        linewidth=2,
        ax=ax2
    )

    hep.histplot(
        hh.divideHists(hist_xsec_stat, hist_xsec, cutoff=0.01, rel_unc=True),
        histtype="errorbar",
        color="black",
        label="Data",
        yerr=True,
        linewidth=2,
        ax=ax2, capsize=2, elinewidth=0, markersize=0
    )

    if asimov:
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
            hh.divideHists(ha_xsec, hist_xsec, cutoff=0.01, rel_unc=True),
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

    if poi_type=="mu":
        outfile = "mu"
    elif normalize:
        outfile = "unfolded_xsec_normalized"
    else: 
        outfile = "unfolded_xsec" 

    outfile += (f"_{args.postfix}" if args.postfix else "")
    plot_tools.save_pdf_and_png(outdir, outfile)

    if asimov:
        asimov_yields = make_yields_df([ha_xsec], ["Model"], per_bin=True)
        asimov_yields["Uncertainty"] *= 0 # artificially set uncertainty on model hard coded to 0
    data_yields = make_yields_df([hist_xsec], ["Data"], per_bin=True)
    plot_tools.write_index_and_log(outdir, outfile, nround=4 if normalize else 2,
        yield_tables={"Data" : data_yields, "Model": asimov_yields} if asimov else {"Data" : data_yields},
        analysis_meta_info=None,
        args=args,
    )



if "xsec" in args.plots:
    plot_xsec_unfolded(poi_type="pmaskedexp")
    plot_xsec_unfolded(poi_type="pmaskedexpnorm")

if "prefit" in args.plots:
    plot("prefit", bins=(None, nbins_reco))
    plot("prefit", bins=(None,int(nbins_reco/2)), channel="minus")
    plot("prefit", bins=(int(nbins_reco/2), nbins_reco), channel="plus")

if "postfit" in args.plots:
    plot("postfit", bins=(None, nbins_reco))
    plot("postfit", bins=(None,int(nbins_reco/2)), channel="minus")
    plot("postfit", bins=(int(nbins_reco/2), nbins_reco), channel="plus")

if "correlation" in args.plots:
    plot_matrix_poi("correlation_matrix_channelmu")
    plot_matrix_poi("correlation_matrix_channelpmaskedexp")
    plot_matrix_poi("correlation_matrix_channelpmaskedexpnorm")

if "covariance" in args.plots:
    plot_matrix_poi("covariance_matrix_channelmu")
    plot_matrix_poi("covariance_matrix_channelpmaskedexp")
    plot_matrix_poi("covariance_matrix_channelpmaskedexpnorm")
