from wremnants.datasets.datagroups import datagroups2016
from wremnants import histselections as sel
from wremnants import plot_tools,theory_tools
from utilities import boostHistHelpers as hh
import matplotlib.pyplot as plt
import argparse
import os
import shutil
import logging
import pathlib
import hist

logging.basicConfig(level=logging.INFO)

xlabels = {
    "pt" : r"p$_{T}^{\ell}$ (GeV)",
    "eta" : r"$\eta^{\ell}$",
    "unrolled" : r"(p$_{T}^{\ell}$, $\eta^{\ell}$) bin",
    "ptll" : r"p$_{\mathrm{T}}^{\ell\ell}$ (GeV)",
    "yll" : r"y$^{\ell\ell}$",
    "mll" : r"m$_{\ell\ell}$ (GeV)",
    "costhetastarll" : r"$\cos{\phi^{\star}_{\ell\ell}}$",
    "phistarll" : r"$\phi^{\star}_{\ell\ell}$",
    "recoil_MET_pt" : r"p$_{\mathrm{T}}^{miss}$ (recoil corr.)",
}

pdfInfo = theory_tools.pdfMapExtended 
pdfNames = [pdfInfo[k]["name"] for k in pdfInfo.keys()]

def pdfUnc(h, pdfName):
    print(pdfName, list(pdfNames).index(pdfName))
    key =  list(pdfInfo.keys())[list(pdfNames).index(pdfName)]
    unc = pdfInfo[key]["combine"]
    scale = pdfInfo[key]["scale"] if "scale" in pdfInfo[key] else 1.
    print(key, pdfName, scale)
    return theory_tools.hessianPdfUnc(h, "tensor_axis_0", unc, scale)

transforms = {}
transforms.update({pdf+"Up" : lambda h,p=pdf: pdfUnc(h, p)[0] for pdf in pdfNames})
transforms.update({pdf+"Down" : lambda h,p=pdf: pdfUnc(h, p)[1] for pdf in pdfNames})
transforms.update({f"dilepton_{pdf}Up" : transforms[pdf+"Up"] for pdf in pdfNames})
transforms.update({f"dilepton_{pdf}Down" : transforms[pdf+"Down"] for pdf in pdfNames})
transforms.update(
    {"massWeightDown" : lambda h: h[{"tensor_axis_0" : 0}],
    "massWeightUp" : lambda h: h[{"tensor_axis_0" : 20}]}
)

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histograms")
parser.add_argument("--wlike", action='store_true', help="Make W like plots")
parser.add_argument("--ratio_to_data", action='store_true', help="Use data as denominator in ratio")
parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
parser.add_argument("--nominalRef", type=str, help="Specify the nominal his if baseName is a variation hist (for plotting alt hists)")
parser.add_argument("--hists", type=str, nargs='+', required=True, choices=xlabels.keys(), help="List of histograms to plot")
parser.add_argument("-c", "--channel", type=str, choices=["plus", "minus", "all"], default="all", help="Select channel to plot")
parser.add_argument("-p", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="test", help="Subfolder for output")
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9, 1.1], help="y range for ratio plot")
parser.add_argument("--rebin", type=int, default=1, help="Rebin (for now must be an int)")
parser.add_argument("--ylim", type=float, nargs=2, help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--xlim", type=float, nargs=2, help="min and max for x axis")
parser.add_argument("-a", "--name_append", type=str, help="Name to append to file name")

subparsers = parser.add_subparsers()
variation = subparsers.add_parser("variation", help="Arguments for adding variation hists")
variation.add_argument("--varName", type=str, nargs='+', required=True, help="Name of variation hist")
variation.add_argument("--varLabel", type=str, nargs='+', required=True, help="Label(s) of variation hist for plotting")
variation.add_argument("--selectAxis", type=str, nargs='+', help="If you need to select a variation axis")
variation.add_argument("--selectEntries", type=int, nargs='+', help="entries to read from the selected axis")
variation.add_argument("--colors", type=str, nargs='+', help="Variation colors")
variation.add_argument("--transform", action='store_true', help="Apply variation-specific transformation")

args = parser.parse_args()

def padArray(ref, matchLength):
    return ref+ref[-1:]*(len(matchLength)-len(ref))

addVariation = hasattr(args, "varName") and args.varName is not None

if addVariation and (args.selectAxis or args.selectEntries):
    if not (args.selectAxis and args.selectEntries):
        raise ValueError("Must --selectAxis and --selectEntires together")
    if len(args.varLabel) != 1 and len(args.varLabel) != len(args.selectEntries):
        raise ValueError("Must specify the same number of args for --selectEntires, and --varLabel")

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

groups = datagroups2016(args.infile, wlike=args.wlike)
if not args.nominalRef:
    nominalName = args.baseName.rsplit("_", 1)[0]
    groups.setNominalName(nominalName)
    groups.loadHistsForDatagroups(args.baseName, syst="")
else:
    nominalName = args.nominalRef
    groups.setNominalName(nominalName)
    groups.loadHistsForDatagroups(nominalName, syst=args.baseName)

exclude = ["Data"]
colors = args.colors if hasattr(args, "colors") and args.colors else ["red", "green", "blue"]
unstack = exclude[:]

if addVariation:
    logging.info(f"Adding variation {args.varName}")
    varLabels = padArray(args.varLabel, args.varName)
    for i, (label,name,color) in enumerate(zip(varLabels, args.varName, colors)):
        name_toload = name.replace("Up", "").replace("Down", "")
        name = name if name != "" else nominalName
        exclude.append(name)
        # A bit hacky but load the hist without Up/Down, which comes from the transform
        groups.addSummedProc(nominalName, name=name_toload, label=label, exclude=exclude,
            relabel=args.baseName, color=color, reload=name_toload != args.baseName)

        varname = name_toload
        if args.selectAxis or (args.transform and name in transforms):
            exclude.append(name_toload)
            if args.transform:
                action = transforms[name]
                varname = name
            else:
                entries = padArray(args.selectAxis, args.varLabel)
                action = lambda x: x[{ax : entries[i]}]
                varname = name+str(entry)

            groups.copyWithAction(action=action, name=varname, refproc=name_toload, 
                refname=args.baseName, label=label, color=colors[i])
        elif (args.transform and name not in transforms):
            logging.warning(f"No known transformation for variation {name}. No transform applied!")

        exclude.append(varname)
        unstack.append(varname)

histInfo = groups.getDatagroups()

prednames = [x for x in reversed(histInfo.keys()) if x not in exclude]
select = {} if args.channel == "all" else {"select" : -1.j if args.channel == "minus" else 1.j}

def collapseSyst(h):
    for ax in ["systIdx", "tensor_axis_0"]:
        if ax in h.axes.name:
            return h[{ax : 0}]
    return h

for h in args.hists:
    action = (lambda x: sel.unrolledHist(collapseSyst(x))) if "unrolled" in h else lambda x: hh.projectNoFlow(collapseSyst(x), h, ["ptll"])
    fig = plot_tools.makeStackPlotWithRatio(histInfo, prednames, histName=args.baseName, ylim=args.ylim, action=action, unstacked=unstack, 
            xlabel=xlabels[h], ylabel="Events/bin", rrange=args.rrange, select=select, binwnorm=1.0,
            ratio_to_data=args.ratio_to_data, rlabel="Pred./Data" if args.ratio_to_data else "Data/Pred.",
            xlim=args.xlim) 
    outfile = f"{h}_{args.baseName}_{args.channel}"+ (f"_{args.name_append}" if args.name_append else "")
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile)
