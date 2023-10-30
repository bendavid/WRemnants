import wremnants
from wremnants import plot_tools,theory_tools,histselections as sel
from utilities import boostHistHelpers as hh
from utilities.io_tools import input_tools, output_tools
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import boost_histogram as bh
import argparse
import pathlib
import os
import pickle
import lz4.frame
import logging
import shutil

xlabels = {
    "pt" : r"p$_{T}^{\ell}$ (GeV)",
    "eta" : r"$\eta^{\ell}$",
    "unrolled" : r"(p$_{T}^{\ell}$, $\eta^{\ell}$) bin",
    "unrolled_gen" : r"($|\mathrm{y}^{Z}|$,p$_{T}^{Z}$) bin",
    "unrolled_gen_hel" : r"($|\mathrm{y}^{Z}|$,p$_{T}^{Z}$) bin",
    "ptVgen" : r"p$_{T}^{Z}$ (GeV)",
    "absYVgen" : r"$|\mathrm{y}^{Z}|$",
    "ptll" : r"p$_{\mathrm{T}}^{\ell\ell}$ (GeV)",
    "yll" : r"y$^{\ell\ell}$",
    "mll" : r"m$_{\ell\ell}$ (GeV)",
}

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histograms")
parser.add_argument("--pdfs", type=str, nargs='+', help="List of histograms to plot", choices=theory_tools.pdfMap.keys(), required=True)
parser.add_argument("-c", "--channel", type=str, choices=["plus", "minus", "all"], default="all", help="Select channel to plot")
parser.add_argument("-p", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="test", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9, 1.1], help="y range for ratio plot")
parser.add_argument("-d", "--datasets", type=str, nargs="+", help="Dataset to plot", required=True)
parser.add_argument("--obs", type=str, nargs='+', choices=xlabels.keys(), help="Observable to plot", required=True)
parser.add_argument("--together", action='store_true', help="y range for ratio plot")
parser.add_argument("--baseName", type=str, help="Name of nominal hist")
parser.add_argument("--ymax", type=float, help="Max value for y axis (if not specified, range set automatically)")
args = parser.parse_args()

for pdf in args.pdfs:
    if pdf not in theory_tools.pdfMap:
        raise ValueError(f"pdf {pdf} is not a valid hist (not defined in theory_tools.pdfMap)")

print(args.datasets[0][0])
if "W" in args.datasets[0][0]:
    xlabels["ptVgen"] = xlabels["ptVgen"].replace("Z", "W")
    xlabels["absYVgen"] = xlabels["absYVgen"].replace("Z", "W")
    xlabels["unrolled_gen"] = xlabels["unrolled_gen"].replace("Z", "W")
    xlabels["unrolled_gen_hel"] = xlabels["unrolled_gen_hel"].replace("Z", "W")

pdfInfo = theory_tools.pdfMap 
pdfNames = [pdfInfo[pdf]["name"] for pdf in args.pdfs]
histNames = pdfNames if not args.baseName else [f"{args.baseName}_{pdfName}" for pdfName in pdfNames]

pdfHists = input_tools.read_all_and_scale(args.infile, args.datasets, histNames)
axis_label = "pdfVar"

uncType = [pdfInfo[pdf]["combine"] for pdf in args.pdfs]
uncScale = [pdfInfo[pdf]["scale"] if "scale" in pdfInfo[pdf] else 1. for pdf in args.pdfs]
uncHists = [[h[{axis_label : 0}], *theory_tools.hessianPdfUnc(h, axis_label, unc, scale)] for h,unc,scale in zip(pdfHists, uncType, uncScale)]
names = [[pdfName+" $\pm1\sigma$", "", ""] for pdfName in pdfNames]
cmap = cm.get_cmap("tab10")
colors = [[cmap(i)]*3 for i in range(len(args.pdfs))]


if "unrolled_gen_hel" in args.obs:
    moments = input_tools.read_all_and_scale(args.infile, args.datasets, ["helicity_moments_scale"])
    coeffs =  theory_tools.scale_angular_moments(moments[0].project('helicity','absYVgen','ptVgen','muRfact','muFfact'))
    moments_pdf = input_tools.read_all_and_scale(args.infile, args.datasets, [f"helicity_{args.baseName}_{pdfName}" for pdfName in pdfNames])
    coeffs_pdf = []
    for moments in moments_pdf:
        coeffs_pdf.append(theory_tools.scale_angular_moments(moments.project('helicity','absYVgen','ptVgen',axis_label)))
    uncHists = [[h[{axis_label : 0}], *theory_tools.hessianPdfUnc(h, axis_label, unc, scale)] for h,unc,scale in zip(coeffs_pdf, uncType, uncScale)]

    # add alphaS
    alphaNames = []
    axis_label = "alphasVar"
    for ipdf,pdf in enumerate(args.pdfs):
        if pdfInfo[pdf]["alphasRange"] == "001":
            alphaNames.append(f"helicity_{args.baseName}_{pdfNames[ipdf]}alphaS001")
        else:
            alphaNames.append(f"helicity_{args.baseName}_{pdfNames[ipdf]}alphaS002")
        alphaHists = input_tools.read_all_and_scale(args.infile, args.datasets, alphaNames)
        alphaHists_hel = []
        for alphaHist in alphaHists:
            alphaHists_hel.append(theory_tools.scale_angular_moments(alphaHist.project('helicity','absYVgen','ptVgen',axis_label)))
        uncHists[ipdf].extend([alphaHists_hel[ipdf][...,0],alphaHists_hel[ipdf][...,1]])
        names[ipdf].extend([pdfNames[ipdf]+"alpha $\pm1\sigma$",""])
        colors[ipdf].extend([[cmap(i)]*2 for i in range(len(args.pdfs),2*len(args.pdfs))][0])
    
    # add QCD scales
    uncHists.append([coeffs[{"muRfact" : 2.j, "muFfact" : 2.j}],coeffs[{"muRfact" : 0.5j, "muFfact" : 0.5j}],coeffs[{"muRfact" : 2.j, "muFfact" : 1.j}], coeffs[{"muRfact" : 0.5j, "muFfact" : 1.j}],coeffs[{"muRfact" : 1.j, "muFfact" : 2.j}],coeffs[{"muRfact" : 1.j, "muFfact" : 0.5j}]])
    names.append(["QCDscale_muRmuFUp","QCDscale_muRmuFDown","QCDscale_muRUp","QCDscale_muRDown","QCDscale_muFUp","QCDscale_muFDown"])
    colors.append([[cmap(i)]*6 for i in range(1)][0])
    # uncHists.append([coeffs[{"muRfact" : 1.j, "muFfact" : 1.j}]])
    # names.append(["QCDscale_central"])
    # colors.append([[cmap(i)]*1 for i in range(2*len(args.pdfs),2*len(args.pdfs)+1)][0])

# TODO
#if args.together:

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)
plot_names = args.pdfs
plot_names.append("QCD_scales")

for obs in args.obs:
    all_hists = []
    all_colors = []
    all_names = []
    for name,color,labels,hists in zip(plot_names,colors, names, uncHists):
        # This is the reference
        if not "unrolled" in obs:
            action = lambda x: x.project(obs)
            hists1D = [action(x) for x in hists]
        else:
            obs2unroll = ["absYVgen","ptVgen"] if "unrolled_gen" in obs else ["pt","eta"]
            action = sel.unrolledHist
            if not "hel" in obs:
                hists1D = [action(x,obs2unroll,binwnorm=True) for x in hists]
            else:
                hists1D = [action(x[{'helicity': 3.j}],obs2unroll,binwnorm=True) for x in hists]

        all_hists.extend(hists1D)
        all_colors.extend(color)
        all_names.extend(labels)
    
    # Add the nominal for reference
    if len(uncHists) > 1:
        hists1D = [all_hists[0], *all_hists]
        plot_cols = [colors[0][0], *all_colors]
        plot_labels = ["", *all_names]
        print(all_names)
        fig = plot_tools.makePlotWithRatioToRef(hists1D, colors=plot_cols, labels=plot_labels, alpha=0.7,
            rrange=args.rrange, ylabel="$\sigma$/bin", xlabel=xlabels[obs], rlabel=f"x/{args.pdfs[0].upper()}", binwnorm=None, nlegcols=1)
        outfile = f"{name}Hist_{obs}_{args.channel}_sigma3"
        ax1, ax2 = fig.axes
        ax1.fill_between(hists1D[0].axes[0].centers,np.minimum.reduce([h.values() for h in hists1D]),np.maximum.reduce([h.values() for h in hists1D]),color="grey",alpha=0.5, label="theory agnostic variation")
        ax2.fill_between(hists1D[0].axes[0].centers,np.minimum.reduce([h.values() for h in hists1D])/hists1D[0].values(),np.maximum.reduce([h.values() for h in hists1D])/hists1D[0].values(),color="grey",alpha=0.5, label="theory agnostic variation")
        plot_tools.save_pdf_and_png(outdir, outfile)
        plot_tools.write_index_and_log(outdir, outfile)
