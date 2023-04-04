from wremnants import plot_tools,theory_tools,histselections as sel
from utilities import input_tools
import matplotlib.pyplot as plt
from matplotlib import cm
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
    "ptVgen" : r"p$_{T}^{Z}$ (GeV)",
    "absYVgen" : r"$|\mathrm{y}^{Z}|$",
    "ptll" : r"p$_{\mathrm{T}}^{\ell\ell}$ (GeV)",
    "yll" : r"y$^{\ell\ell}$",
    "mll" : r"m$_{\ell\ell}$ (GeV)",
}

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histograms")
parser.add_argument("--pdfs", type=str, nargs='+', help="List of histograms to plot", choices=theory_tools.pdfMapExtended.keys(), required=True)
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
    if pdf not in theory_tools.pdfMapExtended:
        raise ValueError(f"pdf {pdf} is not a valid hist (not defined in theory_tools.pdfMapExtended)")

if "Z" in args.datasets[0][0]:
    xlabels["ptVgen"] = xlabels["ptVgen"].replace("Z", "W")
    xlabels["absYVgen"] = xlabels["absYVgen"].replace("Z", "W")

pdfInfo = theory_tools.pdfMapExtended 
pdfNames = [pdfInfo[pdf]["name"] for pdf in args.pdfs]
histNames = pdfNames if not args.baseName or "nominal" in args.baseName else [f"{args.baseName}_{pdfName}" for pdfName in pdfNames]
pdfHists = input_tools.read_all_and_scale(args.infile, args.datasets, histNames)
axis_label = "tensor_axis_0"

uncType = [pdfInfo[pdf]["combine"] for pdf in args.pdfs]
uncScale = [pdfInfo[pdf]["scale"] if "scale" in pdfInfo[pdf] else 1. for pdf in args.pdfs]
uncHists = [(h[{axis_label : 0}], *theory_tools.hessianPdfUnc(h, axis_label, unc, scale)) for h,unc in zip(pdfHists, uncType, uncScale)]
names = [(pdfName+" $\pm1\sigma$", "", "") for pdfName in pdfNames]
cmap = cm.get_cmap("tab10")
colors = [[cmap(i)]*3 for i in range(len(args.pdfs))]

# TODO
#if args.together:

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

for obs in args.obs:
    for name,color,labels,hists in zip(args.pdfs,colors, names, uncHists):
        # This is the reference
        action = sel.unrolledHist if obs == "unrolled" else lambda x: x.project(obs) 
        hists1D = [action(x) for x in hists]
        plot_cols = color
        plot_labels = labels
        # Add the nominal for reference
        if len(uncHists) > 1:
            hists1D = [uncHists[0], *hists1D]
            plot_cols = [*colors[0], *color]
            plot_labels = [*names[0], *labels]
        fig = plot_tools.makePlotWithRatioToRef(hists1D, colors=plot_cols, labels=plot_labels, alpha=0.7,
                rrange=args.rrange, ylabel="$\sigma$/bin", xlabel=xlabels[obs], rlabel=f"x/{args.pdfs[0].upper()}", binwnorm=1.0, nlegcols=1)
        outfile = f"{name}Hist_{obs}_{args.channel}"
        plot_tools.save_pdf_and_png(outdir, outfile)
        plot_tools.write_index_and_log(outdir, outfile)
