import pathlib
import mplhep as hep
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import patches
from matplotlib.ticker import StrMethodFormatter # for setting number of decimal places on tick labels
from utilities import boostHistHelpers as hh,common
from wremnants import histselections as sel
import math
import numpy as np
import re
import os
import logging
import shutil
import sys
import datetime

hep.style.use(hep.style.ROOT)

logger = common.child_logger(__name__)

def figureWithRatio(href, xlabel, ylabel, ylim, rlabel, rrange, xlim=None,
    grid_on_main_plot = False, grid_on_ratio_plot = False, plot_title = None, x_ticks_ndp = None,
    bin_density = 300
):
    if not xlim:
        xlim = [href.axes[0].edges[0], href.axes[0].edges[-1]]
    hax = href.axes[0]
    xlim_range = float(xlim[1] - xlim[0])
    original_xrange = float(hax.edges[-1] - hax.edges[0])
    raw_width = (hax.size/float(bin_density)) * (xlim_range / original_xrange)
    width = math.ceil(raw_width)

    fig = plt.figure(figsize=(8*width,8))
    ax1 = fig.add_subplot(4, 1, (1, 3)) 
    ax2 = fig.add_subplot(4, 1, 4) 

    ax2.set_xlabel(xlabel)
    ax1.set_xlabel(" ")
    ax1.set_ylabel(ylabel)
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    if x_ticks_ndp: ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.' + str(x_ticks_ndp) + 'f}'))
    ax2.set_ylabel(rlabel, fontsize=22)
    ax2.set_ylim(rrange)

    if ylim:
        ax1.set_ylim(ylim)
    else:
        ax1.autoscale(axis='y')

    if grid_on_main_plot:  ax1.grid(which = "both")
    if grid_on_ratio_plot: ax2.grid(which = "both")
    if plot_title: ax1.set_title(plot_title)
    return fig,ax1,ax2

def addLegend(ax, ncols=2, extra_text=None, text_size=20):
    has_extra_text = extra_text is not None
    handles, labels = ax.get_legend_handles_labels()
    
    if has_extra_text:
        #handles.append(patches.Patch(color='none', label=extra_text))
        ax.plot([], [], ' ', ' ')

    shape = np.divide(*ax.get_figure().get_size_inches())
    #TODO: The goal is to leave the data in order, but it should be less hacky
    handles[:] = reversed(handles)
    labels[:] = reversed(labels)
    if len(handles) % 2 and ncols == 2:
        handles.insert(math.floor(len(handles)/2), patches.Patch(color='none', label = ' '))
        labels.insert(math.floor(len(labels)/2), ' ')
    #handles= reversed(handles)
    #labels= reversed(labels)
    ax.legend(handles=handles, labels=labels, prop={'size' : text_size*(0.7 if shape == 1 else 1.3)}, ncol=ncols, loc='upper right')

def makeStackPlotWithRatio(
    histInfo, stackedProcs, histName="nominal", unstacked=None, 
    xlabel="", ylabel="Events/bin", rlabel = "Data/Pred.", rrange=[0.9, 1.1], ylim=None, xlim=None, nlegcols=2,
    binwnorm=None, select={},  action = (lambda x: x), extra_text=None, grid = False, plot_title = None, yscale=None,
    fill_between=False, ratio_to_data=False, baseline=True, legtex_size=20, cms_decor="Preliminary", lumi=16.8,
    bin_density=300,
):
    stack = [action(histInfo[k][histName])[select] for k in stackedProcs if histInfo[k][histName]]
    colors = [histInfo[k]["color"] for k in stackedProcs if histInfo[k][histName]]
    labels = [histInfo[k]["label"] for k in stackedProcs if histInfo[k][histName]]
    fig, ax1, ax2 = figureWithRatio(stack[0], xlabel, ylabel, ylim, rlabel, rrange, xlim=xlim, grid_on_ratio_plot = grid, plot_title = plot_title, bin_density = bin_density)

    hep.histplot(
        stack,
        histtype="fill",
        color=colors,
        label=labels,
        stack=True,
        ax=ax1,
        binwnorm=binwnorm,
    )
    
    data_hist = None
    if "Data" in histInfo and ratio_to_data:
        data_hist = action(histInfo["Data"][histName][select])
        hep.histplot(
            hh.divideHists(sum(stack), data_hist, cutoff=0.01),
            histtype="step",
            color=histInfo[stackedProcs[-1]]["color"],
            label=histInfo[stackedProcs[-1]]["label"],
            yerr=False,
            ax=ax2
        )

    if unstacked:
        if type(unstacked) == str: 
            unstacked = unstacked.split(",")
        ratio_ref = data_hist if data_hist else sum(stack) 
        if baseline:
            hep.histplot(
                hh.divideHists(ratio_ref, ratio_ref, cutoff=1e-8),
                histtype="step",
                color="black",
                yerr=False,
                ax=ax2,
                linewidth=2,
            )

        for proc in unstacked:
            logger.debug(f"Plotting proc {proc}")
            unstack = action(histInfo[proc][histName][select])
            hep.histplot(
                unstack,
                yerr=True if proc == "Data" else False,
                histtype="errorbar" if (proc == "Data" or re.search("^pdf.*_sum", proc)) else "step",
                color=histInfo[proc]["color"],
                label=histInfo[proc]["label"],
                ax=ax1,
                alpha=0.7 if not proc == "Data" else 1.,
                binwnorm=binwnorm,
            )
            # TODO: Add option to leave data off ratio, I guess
            #if proc == "Data":
            #    continue
            hep.histplot(
                hh.divideHists(unstack, ratio_ref, cutoff=0.01),
                histtype="errorbar" if proc == "Data" and not data_hist else "step",
                color=histInfo[proc]["color"],
                label=histInfo[proc]["label"],
                yerr=True if (proc == "Data" and not data_hist) else False,
                linewidth=2,
                ax=ax2
            )

        if fill_between:
            fill_procs = [x for x in unstacked if x != "Data"]
            for up,down in zip(fill_procs[::2], fill_procs[1::2]):
                unstack_up = hh.divideHists(action(histInfo[up][histName][select]), ratio_ref, 1e-6)
                unstack_down = hh.divideHists(action(histInfo[down][histName][select]), ratio_ref, 1e-6)
                ax2.fill_between(unstack_up.axes[0].edges, 
                        np.append(unstack_up.values(), unstack_up.values()[-1]), 
                        np.append(unstack_down.values(), unstack_up.values()[-1]),
                    step='post', color=histInfo[up]["color"], alpha=0.5)

    addLegend(ax1, nlegcols, extra_text)
    fix_axes(ax1, ax2, yscale=yscale)

    if cms_decor:
        scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
        hep.cms.label(ax=ax1, lumi=lumi, fontsize=legtex_size*scale, 
            label=cms_decor, data="Data" in histInfo)

    return fig

def makePlotWithRatioToRef(
    hists, labels, colors, xlabel="", ylabel="Events/bin", rlabel="x/nominal",
    rrange=[0.9, 1.1], ylim=None, xlim=None, nlegcols=2, binwnorm=None, alpha=1.,
    baseline=True, data=False, autorrange=None, grid = False,
    yerr=False, legtext_size=20, plot_title=None, x_ticks_ndp = None, bin_density = 300, yscale=None,
):
    # nominal is always at first, data is always at last, if included
    ratio_hists = [hh.divideHists(h, hists[0], cutoff=0.00001) for h in hists[not baseline:]]
    fig, ax1, ax2 = figureWithRatio(hists[0], xlabel, ylabel, ylim, rlabel, rrange, xlim=xlim, grid_on_ratio_plot = grid, plot_title = plot_title, bin_density = bin_density)
    
    hep.histplot(
        hists[:len(hists) - data],
        histtype="step",
        color=colors[:(len(colors)- data)],
        label=labels[:(len(labels)- data)],
        stack=False,
        ax=ax1,
        yerr=yerr,
        binwnorm=binwnorm,
        alpha=alpha,
    )
    
    if len(hists) > 1:
        hep.histplot(
            ratio_hists[:len(ratio_hists) - data],
            histtype="step",
            color=colors[(not baseline):(len(colors)- data)],
            label=labels[(not baseline):(len(labels)- data)],
            yerr=False,
            stack=False,
            ax=ax2,
            alpha=alpha,
        )
    if data:
        hep.histplot(
            hists[-1],
            histtype="errorbar",
            color=colors[-1],
            label=labels[-1],
            stack=False,
            ax=ax1,
            binwnorm=binwnorm,
            alpha=alpha,
        )
        hep.histplot(
            hh.divideHists(data, hists[0], cutoff=1.e-8),
            histtype="errorbar",
            color=colors[-1],
            label=labels[-1],
            xerr=False,
            yerr=False,
            stack=False,
            ax=ax2,
            alpha=alpha,
        )

    addLegend(ax1, nlegcols, legtext_size)
    
    # This seems like a bug, but it's needed
    if not xlim:
        xlim = [hists[0].axes[0].edges[0], hists[0].axes[0].edges[-1]]
    fix_axes(ax1, ax2)
    if x_ticks_ndp: ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.' + str(x_ticks_ndp) + 'f}'))

    return fig

def fix_axes(ax1, ax2, yscale=None):
    #TODO: Would be good to get this working
    #ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    if yscale:
        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim(ymin, ymax*yscale)
    redo_axis_ticks(ax1, "y")
    redo_axis_ticks(ax2, "x")
    redo_axis_ticks(ax1, "x", True)
    ax1.set_xticklabels([])

def redo_axis_ticks(ax, axlabel, no_labels=False):
    autoloc = ticker.AutoLocator()
    # Need this to avoid a warning when you set the axis values manually
    fixedloc = ticker.FixedLocator(autoloc.tick_values(*getattr(ax, f"get_{axlabel}lim")()))
    getattr(ax, f"{axlabel}axis").set_major_locator(fixedloc)
    ticks = getattr(ax, f"get_{axlabel}ticks")()
    labels = [format_axis_num(x) for x in ticks] if not no_labels else []
    getattr(ax, f"set_{axlabel}ticklabels")(labels)

def format_axis_num(val):
    if type(val) == int or val.is_integer():
        # This is kinda dumb and I might change it
        return f"{val:.0f}" if val > 5 else f"{val:0.1f}"
    return f"{x:0.3g}" if val > 10 else f"{val:0.2g}"

def make_plot_dir(outpath, outfolder):
    full_outpath = "/".join([outpath, outfolder])
    if not os.path.isdir(outpath):
        raise IOError(f"The path {outpath} doesn't not exist. You should create it (and possibly link it to your web area)")
        
    if not os.path.isdir(full_outpath):
        logger.info(f"Creating folder {full_outpath}")
        os.makedirs(full_outpath)

    return full_outpath

def save_pdf_and_png(outdir, basename):
    fname = f"{outdir}/{basename}.pdf"
    plt.savefig(fname, bbox_inches='tight')
    plt.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
    logger.info(f"Wrote file(s) {fname}(.png)")

def write_index_and_log(outpath, logname, indexname="index.php", template_dir=f"{pathlib.Path(__file__).parent}/Templates", yield_tables=None):
    if not os.path.isfile(f"{outpath}/{indexname}"):
        shutil.copyfile(f"{template_dir}/{indexname}", f"{outpath}/{indexname}")

    logdir = f"{outpath}/logs"
    if not os.path.isdir(logdir):
        os.mkdir(logdir)

    with open(f"{logdir}/{logname}.log", "w") as logf:
        meta_info = '-'*80 + '\n' + \
            'Script called at %s\n' % datetime.datetime.now() + \
            'The command was: %s\n' % ' '.join(sys.argv) + \
            '-'*80 + '\n'
        logf.write(meta_info)

        if yield_tables:
            for k,v in yield_tables.items():
                logf.write(f"Yield information for {k}\n")
                logf.write("-"*80+"\n")
                logf.write(str(v.round(2))+"\n\n")
