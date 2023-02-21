from wremnants.datasets.datagroups import datagroups2016
from wremnants import histselections as sel
from wremnants import plot_tools,theory_tools,syst_tools
from utilities import boostHistHelpers as hh,common
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import os
import shutil
import logging
import pathlib
import hist
import re

xlabels = {
    "pt" : r"p$_{T}^{\ell}$ (GeV)",
    "eta" : r"$\eta^{\ell}$",
    "ptll" : r"p$_{\mathrm{T}}^{\ell\ell}$ (GeV)",
    "yll" : r"y$^{\ell\ell}$",
    "mll" : r"m$_{\ell\ell}$ (GeV)",
    "costhetastarll" : r"$\cos{\phi^{\star}_{\ell\ell}}$",
    "phistarll" : r"$\phi^{\star}_{\ell\ell}$",
    "MET_pt" : r"p$_{\mathrm{T}}^{miss}$ (GeV)",
    "mt" : r"m$_{T}^{\ell\nu}$ (GeV)",
    "etaSum":r"$\eta^{\ell(+)} + \eta^{\ell(-)}$",
    "etaDiff":r"$\eta^{\ell(+)} - \eta^{\ell(-)}$",
    # add 2d unrolled plots 
    "pt-eta" : r"(p$_{T}^{\ell}$, $\eta^{\ell}$) bin",
    "ptll-yll":r"p$_{\mathrm{T}}^{\ell\ell}$, y$^{\ell\ell}$ bin",
    "mll-yll":r"m$_{\ell\ell}$, y$^{\ell\ell}$ bin",
    "mll-ptll":r"m$_{\ell\ell}$, p$_{\mathrm{T}}^{\ell\ell}$ bin",
    "mll-etaPlus":r"m$_{\ell\ell}$, $\eta^{\ell(+)}$ bin",
    "mll-etaMinus":r"m$_{\ell\ell}$, $\eta^{\ell(-)}$ bin",
    "etaPlus-etaMinus":r"$\eta^{\ell(+)}$, $\eta^{\ell(-)}$ bin",
    "etaSum-etaDiff":r"$\eta^{\ell(+)} + \eta^{\ell(-)}$, $\eta^{\ell(+)} - \eta^{\ell(-)}$ bin",
    # add 3d unrolled plots 
    "mll-etaPlus-etaMinus":r"m$_{\ell\ell}$, $\eta^{\ell(+)}$, $\eta^{\ell(-)}$ bin",
    "mll-etaSum-etaDiff":r"m$_{\ell\ell}$, $\eta^{\ell(+)} + \eta^{\ell(-)}$, $\eta^{\ell(+)} - \eta^{\ell(-)}$ bin",
}

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Output file of the analysis stage, containing ND boost histogrdams")
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
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("--xlim", type=float, nargs=2, help="min and max for x axis")
parser.add_argument("-a", "--name_append", type=str, help="Name to append to file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--procFilters", type=str, nargs="*", help="Filter to plot (default no filter, only specify if you want a subset")
parser.add_argument("--no_data", action='store_true', help="Don't plot data")
parser.add_argument("--no_fill", action='store_true', help="Don't fill stack")
parser.add_argument("--scaleleg", type=float, default=1.0, help="Scale legend text")

subparsers = parser.add_subparsers(dest="variation")
variation = subparsers.add_parser("variation", help="Arguments for adding variation hists")
variation.add_argument("--varName", type=str, nargs='+', required=True, help="Name of variation hist")
variation.add_argument("--varLabel", type=str, nargs='+', required=True, help="Label(s) of variation hist for plotting")
variation.add_argument("--selectAxis", type=str, nargs='+', help="If you need to select a variation axis")
variation.add_argument("--selectEntries", type=str, nargs='+', help="entries to read from the selected axis")
variation.add_argument("--colors", type=str, nargs='+', help="Variation colors")
variation.add_argument("--double-colors", action='store_true', help="Auto generate colors in pairs (useful for systematics)")
variation.add_argument("--transform", action='store_true', help="Apply variation-specific transformation")
variation.add_argument("--fill_between", action='store_true', help="Fill between uncertainty hists in ratio")

args = parser.parse_args()

logger = common.setup_logger("makeDataMCStackPlot", 4 if args.debug else 3, True)

def padArray(ref, matchLength):
    return ref+ref[-1:]*(len(matchLength)-len(ref))

addVariation = hasattr(args, "varName") and args.varName is not None

entries = []
if addVariation and (args.selectAxis or args.selectEntries):
    if not (args.selectAxis and args.selectEntries):
        raise ValueError("Must --selectAxis and --selectEntries together")
    if len(args.varLabel) != 1 and len(args.varLabel) != len(args.selectEntries):
        raise ValueError("Must specify the same number of args for --selectEntries, and --varLabel"
                         f" found selectEntries={len(args.selectEntries)} and varLabel={len(args.varLabel)}")
    if len(args.varName) < len(args.selectEntries):
        args.varName = padArray(args.varName, args.selectEntries)
    axes = padArray(args.selectAxis, args.varLabel)
    entries = padArray(args.selectEntries, args.varLabel)

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

groups = datagroups2016(args.infile)
datasets = groups.getNames(args.procFilters, exclude=False)
logger.info(f"Will plot datasets {datasets}")

if not args.nominalRef:
    nominalName = args.baseName.rsplit("-", 1)[0]
    groups.setNominalName(nominalName)
    groups.loadHistsForDatagroups(args.baseName, syst="", procsToRead=datasets)
else:
    nominalName = args.nominalRef
    groups.setNominalName(nominalName)
    groups.loadHistsForDatagroups(nominalName, syst=args.baseName, procsToRead=datasets)

exclude = ["Data"] if not args.no_data else []
unstack = exclude[:]

# TODO: In should select the correct hist for the transform, not just the first
transforms = syst_tools.syst_transform_map(nominalName, args.hists[0])

histInfo = groups.getDatagroups()

if addVariation:
    logger.info(f"Adding variation {args.varName}")
    varLabels = padArray(args.varLabel, args.varName)
    # If none matplotlib will pick a random color
    ncols = len(args.varName) if not args.double_colors else int(len(args.varName)/2)
    colors = args.colors if args.colors else [cm.get_cmap("tab10" if ncols < 10 else "tab20")(int(i/2) if args.double_colors else i) for i in range(len(args.varName))]
    for i, (label,name,color) in enumerate(zip(varLabels,args.varName,colors)):
        entry = entries[i] if entries else None
        do_transform = args.transform and entry in transforms
        name = name if name != "" else nominalName
        load_op = {}
        action=None

        if entry and entry.isdigit():
            entry = int(entry)

        if args.selectAxis or do_transform:
            transform_procs = groups.getProcNames(exclude_group=exclude)
            if do_transform:
                action = transforms[entry]["action"]
                if "procs" in transforms[entry]:
                    transform_procs = transforms[entry]["procs"]
                varname = entry
            else:
                ax = axes[i]
                action = lambda x: x[{ax : entry}] if ax in x.axes.name else x
                varname = name+str(entry)
            load_op = {p : action for p in transform_procs}
        else:
            varname = name

        if (args.transform and entry not in transforms):
            logger.warning(f"No known transformation for variation {entry}. No transform applied!")

        reload = name != args.baseName
        # The action map will only work if reloading, otherwise need to apply some transform
        # to the already loaded hist
        if load_op and reload:
            action = None
        groups.addSummedProc(nominalName, relabel=args.baseName, name=name, label=label, exclude=exclude,
            color=color, reload=reload, rename=varname, procsToRead=datasets,
            preOpMap=load_op, action=action)

        exclude.append(varname)
        unstack.append(varname)

groups.sortByYields(args.baseName, nominalName=nominalName)
histInfo = groups.getDatagroups()

logger.info(f"Unstacked processes are {exclude}")
prednames = list(reversed(groups.getNames([d for d in datasets if d not in exclude], exclude=False)))
logger.info(f"Stacked processes are {prednames}")

select = {} if args.channel == "all" else {"charge" : -1.j if args.channel == "minus" else 1.j}

def collapseSyst(h):
    for ax in ["systIdx", "tensor_axis_0", "vars"]:
        if ax in h.axes.name:
            return h[{ax : 0}].copy()
    return h

overflow_ax = ["ptll", "chargeVgen", "massVgen", "ptVgen"]
for h in args.hists:
    if len(h.split("-")) > 1:
        action = lambda x: sel.unrolledHist(collapseSyst(x[select]), obs=h.split("-"))
    else:
        action = lambda x: hh.projectNoFlow(collapseSyst(x[select]), h, overflow_ax)
    fig = plot_tools.makeStackPlotWithRatio(histInfo, prednames, histName=args.baseName, ylim=args.ylim, yscale=args.yscale,
            fill_between=args.fill_between if hasattr(args, "fill_between") else None, action=action, unstacked=unstack, 
            xlabel=xlabels[h], ylabel="Events/bin", rrange=args.rrange, binwnorm=1.0, lumi=groups.lumi,
            ratio_to_data=args.ratio_to_data, rlabel="Pred./Data" if args.ratio_to_data else "Data/Pred.",
            xlim=args.xlim, no_fill=args.no_fill, cms_decor="Preliminary" if not args.no_data else "Simulation Preliminary",
            legtext_size=20*args.scaleleg)
    outfile = f"{h.replace('-','_')}_{args.baseName}_{args.channel}"+ (f"_{args.name_append}" if args.name_append else "")
    plot_tools.save_pdf_and_png(outdir, outfile)
    stack_yields = groups.make_yields_df(args.baseName, prednames, action)
    unstacked_yields = groups.make_yields_df(args.baseName, unstack, action)
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Stacked processes" : stack_yields, "Unstacked processes" : unstacked_yields},
        analysis_meta_info={"AnalysisOutput" : groups.getMetaInfo()},
        args=args,
    )
