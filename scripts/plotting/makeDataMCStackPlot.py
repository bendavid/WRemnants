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
parser.add_argument("--yscale", type=float, help="Scale the upper y axis by this factor (useful when auto scaling cuts off legend)")
parser.add_argument("--xlim", type=float, nargs=2, help="min and max for x axis")
parser.add_argument("-a", "--name_append", type=str, help="Name to append to file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")

subparsers = parser.add_subparsers()
variation = subparsers.add_parser("variation", help="Arguments for adding variation hists")
variation.add_argument("--varName", type=str, nargs='+', required=True, help="Name of variation hist")
variation.add_argument("--varLabel", type=str, nargs='+', required=True, help="Label(s) of variation hist for plotting")
variation.add_argument("--selectAxis", type=str, nargs='+', help="If you need to select a variation axis")
variation.add_argument("--selectEntries", type=int, nargs='+', help="entries to read from the selected axis")
variation.add_argument("--colors", type=str, nargs='+', help="Variation colors")
variation.add_argument("--transform", action='store_true', help="Apply variation-specific transformation")
variation.add_argument("--fill_between", action='store_true', help="Fill between uncertainty hists in ratio")

args = parser.parse_args()

logger = common.setup_base_logger("makeDataMCStackPlot", args.debug)

def padArray(ref, matchLength):
    return ref+ref[-1:]*(len(matchLength)-len(ref))

addVariation = hasattr(args, "varName") and args.varName is not None

if addVariation and (args.selectAxis or args.selectEntries):
    if not (args.selectAxis and args.selectEntries):
        raise ValueError("Must --selectAxis and --selectEntires together")
    if len(args.varLabel) != 1 and len(args.varLabel) != len(args.selectEntries):
        raise ValueError("Must specify the same number of args for --selectEntires, and --varLabel")
    if len(args.varName) == 1 and len(args.selectEntries):
        args.varName = padArray(args.varName, args.selectEntries)
    axes = padArray(args.selectAxis, args.varLabel)
    entries = padArray(args.selectEntries, args.varLabel)

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
unstack = exclude[:]

# TODO: In should select the correct hist for the transform, not just the first
transforms = syst_tools.syst_transform_map(nominalName, args.hists[0])

histInfo = groups.getDatagroups()

if addVariation:
    logger.info(f"Adding variation {args.varName}")
    varLabels = padArray(args.varLabel, args.varName)
    # If none matplotlib will pick a random color
    colors = args.colors if args.colors else [cm.get_cmap("tab10")(i) for i in range(len(args.varName))]
    for i, (label,name,color) in enumerate(zip(varLabels, args.varName, colors)):
        do_transform = args.transform and name in transforms
        name_toload = name if not do_transform else transforms[name]["hist"]
        name = name if name != "" else nominalName
        exclude.append(name)
        load_op = {}
        action=None

        if args.selectAxis or do_transform:
            transform_procs = groups.getProcNames(exclude_group=exclude)
            if args.selectAxis:
                ax = axes[i]
                entry = entries[i]
                action = lambda x: x[{ax : entry}] if ax in x.axes.name else x
                varname = name+str(entry)
            else:
                action = transforms[name]["action"]
                varname = name
                if "procs" in transforms[name]:
                    transform_procs = transforms[name]["procs"]
            load_op = {p : action for p in transform_procs}
        else:
            varname = name_toload

        if (args.transform and name not in transforms):
            logger.warning(f"No known transformation for variation {name}. No transform applied!")

        reload = name_toload != args.baseName
        # The action map will only work if reloading, otherwise need to apply some transform
        # to the already loaded hist
        if load_op and reload:
            action = None
        groups.addSummedProc(args.baseName, name=name_toload, label=label, exclude=exclude,
            color=color, reload=reload, rename=varname, 
            preOpMap=load_op, action=action)

        exclude.append(varname)
        unstack.append(varname)

groups.sortByYields(args.baseName)
histInfo = groups.getDatagroups()

prednames = list(reversed(groups.getNames(exclude)))
select = {} if args.channel == "all" else {"select" : -1.j if args.channel == "minus" else 1.j}

def collapseSyst(h):
    for ax in ["systIdx", "tensor_axis_0"]:
        if ax in h.axes.name:
            return h[{ax : 0}].copy()
    return h

overflow_ax = ["ptll", "chargeVgen", "massVgen", "ptVgen"]
for h in args.hists:
    action = (lambda x: sel.unrolledHist(collapseSyst(x))) if "unrolled" in h else lambda x: hh.projectNoFlow(collapseSyst(x), h, overflow_ax)
    fig = plot_tools.makeStackPlotWithRatio(histInfo, prednames, histName=args.baseName, ylim=args.ylim, yscale=args.yscale,
            fill_between=args.fill_between if hasattr(args, "fill_between") else None, action=action, unstacked=unstack, 
            xlabel=xlabels[h], ylabel="Events/bin", rrange=args.rrange, select=select, binwnorm=1.0,
            ratio_to_data=args.ratio_to_data, rlabel="Pred./Data" if args.ratio_to_data else "Data/Pred.",
            xlim=args.xlim) 
    outfile = f"{h}_{args.baseName}_{args.channel}"+ (f"_{args.name_append}" if args.name_append else "")
    plot_tools.save_pdf_and_png(outdir, outfile)
    stack_yields = groups.make_yields_df(args.baseName, prednames, action)
    unstacked_yields = groups.make_yields_df(args.baseName, unstack, action)
    plot_tools.write_index_and_log(outdir, outfile, 
        yield_tables={"Stacked processes" : stack_yields, "Unstacked processes" : unstacked_yields})
