import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib import cm, ticker
import lz4.frame
import pickle
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools, logging
import hist
import argparse
import os
import h5py
import narf
import pdb

generator_choices = ["horace-nlo", "horace-lo-photos", "horace-qed", "horace-lo", "winhac-nlo", "winhac-lo-photos", "winhac-lo", "MiNNLO"]

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default="w_z_gen_dists_scetlib_dyturboCorr_ewinput.hdf5", help="File containing EW hists")
parser.add_argument("--nums", nargs="+", type=str, default=["horace-nlo"], choices=generator_choices, help="Numerators")
parser.add_argument("--den", type=str, default="horace-lo-photos", choices=generator_choices, help="Denominatos")
parser.add_argument("--normalize", action="store_true", default=False, help="Normalize distributions before computing ratio")
parser.add_argument("--noSmoothing", action="store_true", default=False, help="Disable smoothing of corrections")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--plots", nargs="*", type=str, default=["1D","2D"], choices=["1D","2D", "2Derr"], help="What plots to produce")
parser.add_argument("--showFlow", action='store_true', help="Show underlfow and overflow bins in plots")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-o", "--plotdir", type=str, help="Output directory for plots")

args = parser.parse_args()

logger = logging.setup_logger("make_theory_corr_ew", 4 if args.debug else 3)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']
charge_dict = {'ZToMuMu': 0, 'WplusToMuNu': 1, 'WminusToMuNu': 0}

procs_dict = {
    "ZToMuMu": "ZmumuPostVFP",
    "WminusToMuNu": "WminusmunuPostVFP",
    "WplusToMuNu": "WplusmunuPostVFP",
}

text_dict = {
    "ZToMuMu": r"$\mathrm{Z}\rightarrow\mu\mu$",
    "WplusToMuNu": r"$\mathrm{W}^+\rightarrow\mu\nu$",
    "WminusToMuNu": r"$\mathrm{W}^-\rightarrow\mu\nu$"
}

base_name = 'nominal_ew'
project = ["ewMll", "ewLogDeltaM"]

# base_name = 'nominal_ewMlly'
# project = ["ewMlly"]


# file created with `python WRemnants/scripts/histmakers/w_z_gen_dists.py --skipAngularCoeffs --filter horace -p ewinput`
f = h5py.File(args.input, 'r')
res = narf.ioutils.pickle_load_h5py(f["results"])

corrh = {}

try:
    meta = input_tools.get_metadata(args.input)
except ValueError as e:
    logger.warning(f"No meta data found for file {args.input}")
    pass

labels = {
    "ewMll": "$m_{\ell\ell}$", 
    "ewMlly": "$m_{\ell\ell\gamma}$", 
    "ewLogDeltaM": "$\log_{10}(m_{\ell\ell\gamma} - m_{\ell\ell})$"
}

colors = mpl.colormaps["tab10"]

def make_plot_2d(h, name, proc, plot_error=False, cmin=None, cmax=None, flow=True, density=False, log=False):

    nBinsTotal = h.values(flow=flow).size
    nBins = h.project(*project).values(flow=flow).size
    logger.info(f"Make plot {name} with axes {h.axes.name} and {nBins}/{nBinsTotal} bins")
    # average over bins
    h2d = h.project(*project) #*nBins/nBinsTotal

    xlabel = labels.get(project[0],project[0])
    ylabel = labels.get(project[1],project[1])

    edges = h2d.axes.edges
    xbins = np.reshape(edges[0], len(edges[0]))
    ybins = edges[1][0]

    if flow:
        # add extra bin with bin wdith of 1% of total width
        rangex = xbins[-1] - xbins[0]
        xbins = np.insert(xbins, 0, xbins[0] - rangex*0.02)
        xbins = np.append(xbins, xbins[-1] + rangex*0.02)

        rangey = ybins[-1] - ybins[0]
        ybins = np.insert(ybins, 0, ybins[0] - rangey*0.02)
        ybins = np.append(ybins, ybins[-1] + rangey*0.02)

    # z = h2d.values(flow=flow).T

    if density:
        xbinwidths = xbins[1:]-xbins[:-1]
        ybinwidths = ybins[1:]-ybins[:-1]
        binwidths = np.outer(xbinwidths, ybinwidths) 
        h2d.values(flow=flow)[...] = h2d.values(flow=flow) / binwidths

    if plot_error:
        # plot relative errors instead
        h2d.values(flow=flow)[...] = np.sqrt(hh.relVariance(h2d.values(flow=flow), h2d.variances(flow=flow), fillOnes=True))

    fig, ax = plot_tools.figure(h2d, xlabel=xlabel, ylabel=ylabel, cms_label="Preliminary", automatic_scale=False, width_scale=1.2)

    if log:
        cmin = min(h2d.values(flow=flow)[h2d.values(flow=flow)>0]) if cmin is None else cmin # smallest value that is not 0
        cmax = h2d.values(flow=flow).max() if cmax is None else cmax
        colormesh = ax.pcolormesh(xbins, ybins, h2d.values(flow=flow).T, norm=LogNorm(vmin=cmin, vmax=cmax), cmap=cm.RdBu)
    else:
        cmin = 0 if cmin is None else cmin
        cmax = h2d.values(flow=flow).max() if cmax is None else cmax
        crange = max((cmax-1), (cmin+1))
        colormesh = ax.pcolormesh(xbins, ybins, h2d.values(flow=flow).T, vmin=1-crange, vmax=1+crange, cmap=cm.RdBu)

    cbar = fig.colorbar(colormesh, ax=ax)

    ax.text(1.0, 1.003, text_dict[proc], transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")

    output_tools.make_plot_dir(*args.plotdir.rsplit("/", 1))
    plot_name = f"hist2d_{name}_{proc}"
    if args.noSmoothing and "_div_" in name:
        plot_name += "_noSmoothing"
    if args.normalize:
        plot_name += "_normalize"
    if log:
        plot_name += "_log"
    plot_tools.save_pdf_and_png(args.plotdir, plot_name)
    plot_tools.write_index_and_log(args.plotdir, plot_name, args=args, analysis_meta_info=meta)

def make_plot_1d(hists, name, proc, axis, ratio=False, normalize=False, ymin=None, ymax=None, flow=True, density=False):
    logger.info(f"Make 1D plot for {name} with axes {axis}")

    if not isinstance(hists, list):
        hists = [hists]
        name = [name]
    
    h1ds = [h.project(axis) for h in hists]

    if normalize:
        h1ds = [h/np.sum(h.values(flow=True)) for h in h1ds]
    if density:
        for i, h1d in enumerate(h1ds):
            x = h1d.axes.edges[0]
            if flow:
                # add extra bin with bin wdith of 1% of total width
                rangex = x[-1] - x[0]

                x = np.insert(x, 0, x[0] - rangex*0.02)
                x = np.append(x, x[-1] + rangex*0.02)
            binwidths = x[1:]-x[:-1]
            hh.scaleHist(h1d, 1./binwidths, createNew=False)

    ymax = ymax if ymax is not None else max([max(h.values(flow=args.showFlow)) for h in h1ds])
    ymin = ymin if ymin is not None else min([min(h.values(flow=args.showFlow)) for h in h1ds])
    yrange = ymax - ymin
    ymin = ymin if ymin == 0 else ymin - yrange*0.3
    ymax = ymax + yrange*0.3

    if ratio:
        ylabel = "1/{0}".format(name[0].split("_div_")[-1])
    else:
        ylabel = "a.u."

    fig, ax = plot_tools.figure(h1ds[0], xlabel=labels.get(axis,project[0]), ylabel=ylabel, cms_label="Preliminary", automatic_scale=False, width_scale=1.2,
        ylim=(ymin, ymax))

    if ratio:
        x = h1ds[0].axes.edges[0][:-1]
        ax.plot([min(x), max(x)], [1,1], color="black", linestyle="--")

    outname = ""
    for i, h1d in enumerate(h1ds):
        x = h1d.axes.edges[0][:-1]
        y = h1d.values(flow=flow)
        err = np.sqrt(h1d.variances(flow=flow))

        ax.set_xlim((min(x),max(x)))
        ax.set_ylim((ymin,ymax))

        if flow:
            # add extra bin with bin wdith of 1% of total width
            rangex = x[-1] - x[0]

            x = np.insert(x, 0, x[0] - rangex*0.02)
            x = np.append(x, x[-1] + rangex*0.02)

        if args.showFlow:
            ax.set_xlim((min(x),max(x)))

        ax.step(x, y, color=colors(i), label=name[i].split("_div_")[0], where="post")
        ax.fill_between(x, y - err, y + err, alpha=0.3, color=colors(i), step="post")

        # outname += f"_{name[i]}"

    ax.text(1.0, 1.003, text_dict[proc], transform=ax.transAxes, fontsize=30,
            verticalalignment='bottom', horizontalalignment="right")
    plot_tools.addLegend(ax, text_size=16)

    output_tools.make_plot_dir(*args.plotdir.rsplit("/", 1))
    plot_name = f"hist{outname}_{axis}_{proc}"
    if ratio:
        plot_name += "_div"
    if args.noSmoothing and "_div_" in name[0]:
        plot_name += "_noSmoothing"
    if args.normalize:
        plot_name += "_normalize"
    plot_tools.save_pdf_and_png(args.plotdir, plot_name)
    plot_tools.write_index_and_log(args.plotdir, plot_name, args=args, analysis_meta_info=meta)

for proc in procs:

    # Make 2D ratio
    logger.info(f'Make 2D ratio for {proc}')
    def prepare_hist(name):
        if name == 'MiNNLO':
            proc_name = procs_dict[proc]
        else:
            proc_name = f'{proc}_{name}'

        if proc_name not in res:
            return None

        histo = res[proc_name]['output'][base_name].get()
        logger.info(f'Integrals for {name} {np.sum(histo.values(flow=True))}')

        if args.normalize:
            histo = hh.normalize(histo)
            logger.info(f'Integrals for {name} after normalizing {np.sum(histo.values(flow=True))}')
        else:
            histo = hh.scaleHist(histo, 10e6/res[proc_name]['event_count'], createNew=False)
            logger.info(f'Integrals for {name} after scaling {np.sum(histo.values(flow=True))}')

        return histo
    
    nums = []
    hnums = []
    for num in args.nums:
        hnum = prepare_hist(num)
        if hnum is not None:
            nums.append(num)
            hnums.append(hnum)

    hden = prepare_hist(args.den)
    
    hratios = [hh.divideHists(hnum, hden) for hnum in hnums]

    corrhs = {}
    for num, hratio in zip(nums, hratios):
        if not args.noSmoothing:
            hratio = hh.smoothTowardsOne(hratio)
            scale = np.sum(hden.values(flow=True)) / np.sum(hden.values(flow=True)*hratio.values(flow=True))
            logger.info('Adjustment after smoothing {0}'.format(scale))
            hratio = hh.scaleHist(hratio, scale)
            logger.info('Integrals after adjustment {0} {1}'.format(np.sum(hden.values(flow=True)), np.sum(hden.values(flow=True)*hratio.values(flow=True))))

        # Add dummy axis
        axis_dummy = hist.axis.Regular(1, -10., 10., underflow=False, overflow=False, name = "dummy")
        hdummy = hist.Hist(*hratio.axes, axis_dummy, storage=hist.storage.Weight())
        hdummy.view(flow=True)[...,0] = hratio.view(flow=True)
        hratio = hdummy

        # Add charge axis
        if proc[0] == 'W':
            axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
        elif proc[0] == 'Z':
            axis_charge = hist.axis.Regular(1, -1., 1., underflow=False, overflow=False, name = "charge")
        hcharge = hist.Hist(*hratio.axes, axis_charge, storage=hist.storage.Weight())
        hcharge.view(flow=True)[...,charge_dict[proc]] = hratio.view(flow=True)
        hratio = hcharge

        # Add syst axis
        corrh[proc] = hist.Hist(*hratio.axes, hist.axis.Regular(3, 0, 3, underflow=False, overflow=False, name="systIdx"), storage=hist.storage.Weight())
        # Variations: 0=original MiNNLO, 1=Horace NLO, 2=mirrored
        hones = hist.Hist(*hratio.axes, storage=hist.storage.Weight())
        hones.values(flow=True)[...,charge_dict[proc]] = np.ones(hdummy.values(flow=True).shape)
        hmirror = hist.Hist(*hratio.axes, storage=hist.storage.Weight())
        hmirror.values(flow=True)[...] = 2*hratio.values(flow=True) - hones.values(flow=True)
        mirrorscale = np.sum(hden.values(flow=True)) / np.sum(hden.values(flow=True)*hmirror.values(flow=True)[...,0,charge_dict[proc]])
        # logger.info(f'{proc} mirrorscale = {mirrorscale}')
        hmirror = hh.scaleHist(hmirror, mirrorscale)
        corrh[proc].values(flow=True)[...,0] = hones.values(flow=True)
        corrh[proc].values(flow=True)[...,1] = hratio.values(flow=True)
        corrh[proc].values(flow=True)[...,2] = hmirror.values(flow=True)

        corrhs[num] = corrh

    if "2D" in args.plots:
        logger.info("Make 2D control plots")
        make_plot_2d(hden, args.den, proc, density=True, cmin=10e-6, cmax=10e6, log=True)
        if "2Derr" in args.plots:
            make_plot_2d(hden, f"{args.den}_err", proc, plot_error=True, log=True)

        for i, hnum in enumerate(hnums):
            make_plot_2d(hnum, nums[i], proc, density=True, cmin=10e-6, cmax=10e6, log=True)
            if "2Derr" in args.plots:
                make_plot_2d(hnum, f"{nums[i]}_err", proc, plot_error=True, log=True)

        for i, hratio in enumerate(hratios):
            make_plot_2d(hratio, f"{nums[i]}_div_{args.den}", proc, cmin=0, cmax=3)
            if "2Derr" in args.plots:
                make_plot_2d(hratio, f"{nums[i]}_div_{args.den}_err", proc, plot_error=True, log=True)

    if "1D" in args.plots:
        logger.info("Make 1D control plots")

        for ax in project:
            if ax == "ewMll":
                ymin, ymax = 0.8, 1.2
            else:
                ymin, ymax = 0, 4

            hratios1D = [hh.divideHists(hnum.project(ax), hden.project(ax)) for hnum in hnums]

            make_plot_1d(hratios1D, [f"{n}_div_{args.den}" for n in nums], proc, ax, ymin=ymin, ymax=ymax, ratio=True)
            make_plot_1d([*hnums, hden], [*nums, args.den], proc, ax, ymin=0, density=False)

for num, corrh in corrhs.items():
    outname = num.replace('-', '') + 'ew'
    if 'ZToMuMu' in corrh:
        with lz4.frame.open(f"{args.outpath}/{outname}CorrZ.pkl.lz4", "wb") as f:
            pickle.dump({
                    'Z' : {
                        f"{outname}_minnlo_ratio" : corrh['ZToMuMu'],
                    },
                    "meta_data" : output_tools.metaInfoDict(args=args),
                }, f, protocol = pickle.HIGHEST_PROTOCOL)

    if 'WplusToMuNu' in corrh and "WminusToMuNu" in corrh:
        corrh['W'] = corrh['WplusToMuNu']+corrh['WminusToMuNu']
        with lz4.frame.open(f"{args.outpath}/{outname}CorrW.pkl.lz4", "wb") as f:
            pickle.dump({
                    'W' : {
                        f"{outname}_minnlo_ratio" : corrh['W'],
                    },
                    "meta_data" : output_tools.metaInfoDict(args=args),
                }, f, protocol = pickle.HIGHEST_PROTOCOL)

# logger.info(corrh['ZToMuMu'])