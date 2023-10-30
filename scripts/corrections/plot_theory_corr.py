import argparse
import os
import numpy as np
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LogNorm

from utilities import logging, boostHistHelpers as hh
from utilities.io_tools import output_tools
from utilities.styles import styles
from wremnants import theory_corrections, plot_tools

parser = argparse.ArgumentParser()
parser.add_argument("--theoryCorr", nargs="*", default=["scetlib_dyturbo", "horacenloew"], choices=theory_corrections.valid_theory_corrections(),
    help="Apply corrections from indicated generator. First will be nominal correction.")
parser.add_argument("--datasets", nargs="*", default=["ZmumuPostVFP"], 
    help="Apply corrections from indicated generator. First will be nominal correction.")
parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4],
                    help="Set verbosity level with logging, the larger the more verbose")
parser.add_argument("--noColorLogger", action="store_true", help="Do not use logging with colors")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./test", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--cmsDecor", default="Preliminary", type=str, choices=[None,"Preliminary", "Work in progress", "Internal"], help="CMS label")
parser.add_argument("--flow", action='store_true', help="Show underlfow and overflow bins in plots")
parser.add_argument("--axes", type=str, nargs="*", default=None, help="Which axes to plot, if not specified plot all axes")
parser.add_argument("--xlim", type=float, nargs=2, default=[None,None], help="Min and max values for x axis (if not specified, range set automatically)")
parser.add_argument("--ylim", type=float, nargs=2, default=[None,None], help="Min and max values for y axis (if not specified, range set automatically)")
parser.add_argument("--plots", type=str, nargs="+", default=["1d", "2d"], choices=["1d", "2d"], help="Define which plots to make")


args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

colors = mpl.colormaps["tab10"]

corr_dict = theory_corrections.load_corr_helpers(args.datasets, args.theoryCorr, make_tensor=False)

def make_plot_2d(h, name, proc, axes, corr=None, plot_error=False, cmin=None, cmax=None, flow=True, density=False, log=False):

    logger.info(f"Make plot {name} with axes {h.axes.name}")
    # average over bins

    h2d = h.project(*axes)

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

    xlim = (xedges[0],xedges[-1])
    ylim = (yedges[0],yedges[-1])

    fig, ax = plot_tools.figure(h2d, xlabel=xlabel, ylabel=ylabel, cms_label=args.cmsDecor, automatic_scale=False, width_scale=1.2, xlim=xlim, ylim=ylim)

    if log:
        cmin = min(h2d.values(flow=flow)[h2d.values(flow=flow)>0]) if cmin is None else cmin # smallest value that is not 0
        cmax = h2d.values(flow=flow).max() if cmax is None else cmax
        colormesh = ax.pcolormesh(xedges, yedges, h2d.values(flow=flow).T, norm=LogNorm(vmin=cmin, vmax=cmax), cmap=cm.RdBu)
    else:
        cmin = max(0.95,h2d.values(flow=flow).min()) if cmin is None else cmin
        cmax = min(1.05,h2d.values(flow=flow).max()) if cmax is None else cmax
        crange = max((cmax-1), (1-cmin))
        colormesh = ax.pcolormesh(xedges, yedges, h2d.values(flow=flow).T, cmap=cm.RdBu, vmin=max(0.95,1-crange), vmax=min(1.05,1+crange))

    cbar = fig.colorbar(colormesh, ax=ax)

    ax.text(1.0, 1.003, styles.text_dict[proc], transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")

    outfile = f"hist2d_{'_'.join(axes)}_{proc}_{name}"
    if corr:
        outfile += f"_{corr}"
    # if ratio:
    #     outfile += "_div"
    # if args.noSmoothing and "_div_" in name[0]:
    #     outfile += "_noSmoothing"
    # if args.normalize:
    #     outfile += "_normalize"
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, args=args)

def make_plot_1d(hists, names, proc, axis, corr=None, ratio=False, normalize=False, xmin=None, xmax=None, ymin=None, ymax=None, flow=True, density=False):
    logger.info(f"Make 1D plot for {names} with axis {axis}")

    if not isinstance(hists, list):
        hists = [hists]
        names = [names]

    if flow:
        project = lambda h, axis: h.project(axis)/np.prod([a.extent for a in h.axes if a.name != axis])
    else:
        project = lambda h, axis: hh.projectNoFlow(h, axis)/np.prod([a.size for a in h.axes if a.name != axis])

    h1ds = [project(h, axis) for h in hists]

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
        ylabel = "1/{0}".format(names[0].split("_div_")[-1])
    else:
        ylabel = "a.u."

    fig, ax = plot_tools.figure(h1ds[0], xlabel=styles.axis_labels.get(axis, axis), ylabel=ylabel, cms_label=args.cmsDecor, automatic_scale=False, width_scale=1.2,
        ylim=(ymin, ymax), xlim=(xmin, xmax))

    if ratio:
        ax.plot([min(xedges), max(xedges)], [1,1], color="black", linestyle="--")

    if corr:
        labels = styles.syst_labels[corr]
    else:
        labels = {}

    for i, h1d in enumerate(h1ds):
        y = h1d.values(flow=flow)

        ax.stairs(y, xedges, color=colors(i), label=labels.get(names[i], names[i]))
        # err = np.sqrt(h1d.variances(flow=flow))
        # ax.bar(x=xedges[:-1], height=2*err, bottom=y - err, width=np.diff(xedges), align='edge', linewidth=0, alpha=0.3, color=colors(i), zorder=-1)

    ax.text(1.0, 1.003, styles.text_dict[proc], transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")
    plot_tools.addLegend(ax, ncols=1+int(len(names)/7), text_size=12)

    outfile = f"hist_{axis}_{proc}"
    if corr:
        outfile += f"_{corr}"
    # if ratio:
    #     outfile += "_div"
    # if args.noSmoothing and "_div_" in name[0]:
    #     outfile += "_noSmoothing"
    # if args.normalize:
    #     outfile += "_normalize"
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, args=args)


for dataset, corr_hists in corr_dict.items():
    for corr, corrh in corr_hists.items():
        if "systIdx" in corrh.axes.name:
            syst_axis = "systIdx" 
        elif "vars" in corrh.axes.name:
            syst_axis = "vars" 
        else:
            raise RuntimeError(f"Systematics axis not found, available axes are {corrh.axes.name}")

        # loop over charge, if exist
        if "charge" in corrh.axes.name:
            corrh_charges = [corrh[{"charge":idx}] for idx in range(corrh.axes["charge"].size)]
        else: 
            corrh_charges = [corrh]

        for corrh_charge in corrh_charges:           
            # split hists into systematics
            corrh_systs = {idx: corrh_charge[{syst_axis:idx}] for idx in range(corrh_charge.axes[syst_axis].size)}

            hists = [h for h in corrh_systs.values()]
            names = [f"{n}" for n in corrh_systs.keys()]
            axes = [n for n in corrh_charge.axes.name if n != syst_axis and (args.axes is None or n in args.axes)]
            
            if len(axes) == 2 and "2d" in args.plots:
                for n, h in corrh_systs.items():
                    make_plot_2d(h, n, dataset.replace("PostVFP",""), axes, corr=corr, flow=args.flow)                

            if "1d" not in args.plots:
                continue

            for axis in axes:
                make_plot_1d(hists, names, dataset.replace("PostVFP",""), axis, corr=corr, flow=args.flow,
                    xmin=args.xlim[0], xmax=args.xlim[1], ymin=args.ylim[0], ymax=args.ylim[1])
        
