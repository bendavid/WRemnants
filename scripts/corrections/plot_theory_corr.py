import argparse
import os
import numpy as np
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LogNorm
import hist
from itertools import combinations

from utilities import common, logging, boostHistHelpers as hh
from utilities.io_tools import output_tools
from utilities.styles import styles
from wremnants import theory_corrections, plot_tools

parser = common.plot_parser()
parser.add_argument("--theoryCorr", nargs="*", default=["scetlib_dyturbo", "horacenloew"], #choices=theory_corrections.valid_theory_corrections(),
    help="Apply corrections from indicated generator. First will be nominal correction.")
parser.add_argument("--idxs", nargs="*", default=None, help="Indexes from systematic axis to be used for plotting.")
parser.add_argument("--datasets", nargs="*", default=["ZmumuPostVFP"], 
    help="Apply corrections from indicated generator. First will be nominal correction.")
parser.add_argument("--baseDir", type=str, default=f"{common.data_dir}/TheoryCorrections/", help="Base directory to the theory corrections")
parser.add_argument("--noFlow", action='store_true', help="Do not show underlfow and overflow bins in plots")
parser.add_argument("--showUncertainties", action='store_true', help="Show uncertainty bands")
parser.add_argument("--axes", type=str, nargs="*", default=None, help="Which axes to plot, if not specified plot all axes")
parser.add_argument("--xlim", type=float, nargs=2, default=[None,None], help="Min and max values for x axis (if not specified, range set automatically)")
parser.add_argument("--ylim", type=float, nargs=2, default=[None,None], help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--clim", type=float, nargs=2, default=[None,None], help="Min and max values for color in 2d plot (if not specified, range set automatically)")
parser.add_argument("--plots", type=str, nargs="+", default=["1d", "2d"], choices=["1d", "2d"], help="Define which plots to make")


args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

colors = mpl.colormaps["tab10"]

corr_dict = theory_corrections.load_corr_helpers(args.datasets, args.theoryCorr, make_tensor=False, base_dir=args.baseDir)

def project(flow=False):
    # average over bins
    if flow:
        return lambda h, axes: h.project(*axes)/np.prod([a.extent for a in h.axes if a.name not in axes])
    else:
        return lambda h, axes: hh.projectNoFlow(h, axes)/np.prod([a.size for a in h.axes if a.name not in axes])

def make_plot_2d(h, name, proc, axes, corr=None, plot_error=False, clim=None, flow=True, density=False, log=False):
    logger.info(f"Make 2d plot {name} with axes {axes[0]}, {axes[1]}")
    
    h2d = project(flow)(h, axes)

    xlabel = styles.axis_labels.get(axes[0],axes[0])
    ylabel = styles.axis_labels.get(axes[1],axes[1])

    if flow:
        xedges, yedges = plot_tools.extendEdgesByFlow(h2d)
    else:
        edges = h2d.axes.edges
        xedges = np.reshape(edges[0], len(edges[0]))
        yedges = edges[1][0]

    if density:
        xbinwidths = np.diff(xedges)
        ybinwidths = np.diff(yedges)
        binwidths = np.outer(xbinwidths, ybinwidths) 
        h2d.values(flow=flow)[...] = h2d.values(flow=flow) / binwidths

    if plot_error:
        # plot relative errors instead
        h2d.values(flow=flow)[...] = np.sqrt(hh.relVariance(h2d.values(flow=flow), h2d.variances(flow=flow), fillOnes=True))

    if args.xlim[0] is None:
        xlim = (xedges[0],xedges[-1])
    else:
        xlim = args.xlim

    if args.ylim[0] is None:
        ylim = (yedges[0],yedges[-1])
    else:
        ylim = args.ylim

    fig, ax = plot_tools.figure(h2d, xlabel=xlabel, ylabel=ylabel, cms_label=args.cmsDecor, automatic_scale=False, width_scale=1.2, xlim=xlim, ylim=ylim)

    if clim is None:
        if log:
            cmin = min(h2d.values(flow=flow)[h2d.values(flow=flow)>0]) # smallest value that is not 0
            cmax = h2d.values(flow=flow).max()
        else:
            cmin = max(0.95,h2d.values(flow=flow).min())
            cmax = min(1.05,h2d.values(flow=flow).max())
        # make symmetric range
        crange = max((cmax-1), (1-cmin))
        clim = [max(0.95,1-crange), min(1.05,1+crange)]
    else:
        colormesh = ax.pcolormesh(xedges, yedges, h2d.values(flow=flow).T, norm=LogNorm(vmin=clim[0], vmax=clim[1]), cmap=cm.RdBu)

    if log:
        colormesh = ax.pcolormesh(xedges, yedges, h2d.values(flow=flow).T, cmap=cm.RdBu, norm=LogNorm(vmin=clim[0], vmax=clim[1]))
    else:
        colormesh = ax.pcolormesh(xedges, yedges, h2d.values(flow=flow).T, cmap=cm.RdBu, vmin=clim[0], vmax=clim[1])

    cbar = fig.colorbar(colormesh, ax=ax)

    ax.text(1.0, 1.003, styles.text_dict.get(proc, proc), transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")

    outfile = f"hist2d_{'_'.join(axes)}_{proc}_{name}"
    if corr:
        outfile += f'_{corr.replace("(","").replace(")","")}'
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, args=args)

def make_plot_1d(hists, names, proc, axis, labels=None, corr=None, 
    ratio=True, normalize=False, xmin=None, xmax=None, ymin=None, ymax=None, flow=True, density=False, uncertainty_bands=False
):
    logger.info(f"Make 1D plot for corr {corr} with {len(names)} entries for axis {axis}")

    if not isinstance(hists, list):
        hists = [hists]
        names = [names]

    h1ds = [project(flow)(h, [axis]) if set(hists[0].axes.name) != set([axis]) else h for h in hists]

    if flow:
        xedges = plot_tools.extendEdgesByFlow(h1ds[0])
    else:
        xedges = h1ds[0].axes.edges[0]

    if normalize:
        h1ds = [h/np.sum(h.values(flow=flow)) for h in h1ds]
    if density:
        for i, h1d in enumerate(h1ds):
            binwidths = xedges[1:]-xedges[:-1]
            hh.scaleHist(h1d, 1./binwidths, createNew=False)

    if xmin is None:
        xmin, xmax = (xedges[0],xedges[-1])

    if ymin is None or ymax is None:
        xmap = (xedges[1:] > xmin) & (xedges[:-1] < xmax)
        ymax = ymax if ymax is not None else max([max(h.values(flow=flow)[xmap]) for h in h1ds])
        ymin = ymin if ymin is not None else min([min(h.values(flow=flow)[xmap]) for h in h1ds])
        yrange = ymax - ymin
        ymin = ymin if ymin == 0 else ymin - yrange*0.3
        ymax = ymax + yrange*0.3

    if ratio:
        ylabel = "correction"# "1/{0}".format(names[0].split("_div_")[-1])
    else:
        ylabel = "a.u."

    fig, ax = plot_tools.figure(h1ds[0], xlabel=styles.axis_labels.get(axis, axis), ylabel=ylabel, cms_label=args.cmsDecor, automatic_scale=False, width_scale=1.2,
        ylim=(ymin, ymax), xlim=(xmin, xmax))

    if ratio:
        ax.plot([min(xedges), max(xedges)], [1,1], color="black", linestyle="--", zorder=-1)

    for i, h1d in enumerate(h1ds):
        y = h1d.values(flow=flow)

        ax.stairs(y, xedges, color=colors(i), label=labels[i])

        if uncertainty_bands:
            err = np.sqrt(h1d.variances(flow=flow))
            ax.bar(x=xedges[:-1], height=2*err, bottom=y - err, width=np.diff(xedges), align='edge', linewidth=0, alpha=0.3, color=colors(i), zorder=-1)

    ax.text(1.0, 1.003, styles.text_dict.get(proc, proc), transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")
    plot_tools.addLegend(ax, ncols=1+int(len(names)/7), text_size=12)

    outfile = f"hist_{axis}_{proc}"
    if corr:
        outfile += f"_{corr}"
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, args=args)


for dataset, corr_hists in corr_dict.items():
    base_proc = dataset.replace("PostVFP","")

    # loop over charge, if multiple charges exist make plots for each charge separately
    for charge_idx in [0, 1]:

        all_hists = []
        all_names = []
        all_labels = []
        all_axes = []
        
        for corr, corrh in corr_hists.items():
            has_charge=False
            sel = {} # bin selections
            if charge_idx == 1 and ("charge" not in corrh.axes.name or corrh.axes["charge"].size <= 1):
                continue # only one charge for this correction

            proc = base_proc
            if "charge" in corrh.axes.name:
                sel.update({"charge": charge_idx})
                if corrh.axes["charge"].size > 1:
                    has_charge=True
                    proc = f"{base_proc[0]}{'minus' if charge_idx==0 else 'plus'}{base_proc[1:]}"

            for systAxName in ["systIdx", "tensor_axis_0", "vars"]:
                if systAxName in corrh.axes.name:
                    syst_axis = systAxName
                    break
            else:
                raise RuntimeError(f"Systematics axis not found, available axes are {corrh.axes.name}")

            corr = corr.split("_")[0]
            # retreive 
            if type(corrh.axes[syst_axis]) == hist.axis.StrCategory:
                idxs = [i for i, idx in enumerate(corrh.axes[syst_axis]) if args.idxs is None or str(idx) in args.idxs or str(i) in args.idxs]
                labels = idxs
            else:
                idxs = [i for i in range(corrh.axes[syst_axis].size) if args.idxs is None or str(i) in args.idxs]
                if corr in styles.syst_labels:
                    label = styles.syst_labels[corr]
                else:
                    label = corr
                if type(label) == dict:
                    labels = [label.get(i, f"{label}({i})") for i in idxs] 
                else:
                    labels = [f"{label}({i})" for i in idxs] 
            
            if len(idxs) == 0:
                raise RuntimeError(f"No index found in systematic axis!")
            elif args.idxs is not None and len(idxs) != len(args.idxs):
                logger.warning(f"Some of the indices ({set(args.idxs) - set(idxs)}) have not been found in the systematic axis!")

            # split hists into systematics
            corrh_systs = {idx: corrh[{**sel, syst_axis:idx}] for idx in idxs}

            hists = [h for h in corrh_systs.values()]
            names = [f"{n}" for n in corrh_systs.keys()]
            all_hists+=hists
            all_names+=names
            all_labels+=labels

            axes = [n for n in corrh.axes.name if n not in [syst_axis, "charge"] and (args.axes is None or n in args.axes)]
            all_axes+=axes
            if "2d" in args.plots and len(axes) >= 2:
                for label, (n, h) in zip(labels, corrh_systs.items()):
                    for ax1, ax2 in list(combinations(axes, 2)): 
                        make_plot_2d(h, n, proc, [ax1, ax2], corr=corr, flow=not args.noFlow, clim=args.clim)                

            # if "1d" in args.plots:
            #     for axis in axes:
            #         make_plot_1d(hists, names, proc, axis, labels=labels, flow=not args.noFlow, corr=corr, 
            #             xmin=args.xlim[0], xmax=args.xlim[1], ymin=args.ylim[0], ymax=args.ylim[1], uncertainty_bands=args.showUncertainties)

        if has_charge:
            proc = f"{base_proc[0]}{'minus' if charge_idx==0 else 'plus'}{base_proc[1:]}"

        if "1d" in args.plots:
            for axis in set(all_axes):
                make_plot_1d(all_hists, all_names, proc, axis, labels=all_labels, flow=not args.noFlow, corr="all",
                    xmin=args.xlim[0], xmax=args.xlim[1], ymin=args.ylim[0], ymax=args.ylim[1], uncertainty_bands=args.showUncertainties)

            
