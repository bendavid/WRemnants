import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default="w_z_gen_dists_scetlib_dyturboCorr_ewinput.hdf5", help="File containing EW hists")
parser.add_argument("--num", type=str, default="horace-nlo", help="Numerator")
parser.add_argument("--den", type=str, default="horace-lo-photos", help="Denominator")
parser.add_argument("--noSmoothing", action="store_true", default=False, help="Disable smoothing of corrections")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("-o", "--plotdir", type=str, help="Output directory for plots")
parser.add_argument("--logPlots", action='store_true', help="Plot with logarithmic color code")

args = parser.parse_args()

logger = logging.setup_logger("make_theory_corr_ew", 4 if args.debug else 3)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']
charge_dict = {'ZToMuMu': 0, 'WplusToMuNu': 1, 'WminusToMuNu': 0}

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
    "ewLogDeltaM": "$\log_{10}(m_{\ell\ell} - m_{\ell\ell\gamma})$"
}

for proc in procs:
    # Make 2D ratio
    logger.info(f'Make 2D ratio for {proc}')
    hnum = res[f'{proc}_{args.num}']['output']['nominal_ew'].get()
    hden = res[f'{proc}_{args.den}']['output']['nominal_ew'].get()
    logger.info('Integrals {0} {1}'.format(np.sum(hnum.values(flow=True)), np.sum(hden.values(flow=True))))
    hnum = hh.normalize(hnum)
    hden = hh.normalize(hden)
    logger.info('Integrals after normalizing {0} {1}'.format(np.sum(hnum.values(flow=True)), np.sum(hden.values(flow=True))))
    hratio = hh.divideHists(hnum, hden)
    if not args.noSmoothing:
        hratio = hh.smoothTowardsOne(hratio)
    scale = np.sum(hden.values(flow=True)) / np.sum(hden.values(flow=True)*hratio.values(flow=True))
    logger.info('Adjustment after smoothing {0}'.format(scale))
    hratio = hh.scaleHist(hratio, scale)
    # hratio.values(flow=True)[...]    *= scale
    # hratio.variances(flow=True)[...] *= scale
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
    # hmirror.values(flow=True)[...]    *= mirrorscale
    # hmirror.variances(flow=True)[...] *= mirrorscale
    corrh[proc].values(flow=True)[...,0] = hones.values(flow=True)
    corrh[proc].values(flow=True)[...,1] = hratio.values(flow=True)
    corrh[proc].values(flow=True)[...,2] = hmirror.values(flow=True)

    logger.info("Make control plots")

    project = ["ewMll", "ewLogDeltaM"]

    def make_plot_2d(h, name, plot_error=False, cmin=None, cmax=None, flow=True):

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

        if plot_error:
            # plot relative errors instead
            h2d.values(flow=flow)[...] = np.sqrt(hh.relVariance(h2d.values(flow=flow), h2d.variances(flow=flow), fillOnes=True))

        fig, ax = plot_tools.figure(h2d, xlabel=xlabel, ylabel=ylabel, cms_label="Preliminary", automatic_scale=False, width_scale=1.2)

        if args.logPlots:
            cmin = min(h2d.values(flow=flow)[h2d.values(flow=flow)>0]) # smallest value that is not 0
            cmax = h2d.values(flow=flow).max()
            colormesh = ax.pcolormesh(xbins, ybins, h2d.values(flow=flow).T, norm=LogNorm(vmin=1, vmax=1000))
        else:
            cmin = 0 if cmin is None else cmin
            cmax = h2d.values(flow=flow).max() if cmax is None else cmax
            colormesh = ax.pcolormesh(xbins, ybins, h2d.values(flow=flow).T, vmin=cmin, vmax=cmax)

        cbar = fig.colorbar(colormesh, ax=ax)

        output_tools.make_plot_dir(*args.plotdir.rsplit("/", 1))
        plot_name = f"hist2d_{name}_{proc}"
        if args.noSmoothing and "_div_" in name:
            plot_name += "_noSmoothing"
        if args.logPlots:
            plot_name += "_log"
        plot_tools.save_pdf_and_png(args.plotdir, plot_name)
        plot_tools.write_index_and_log(args.plotdir, plot_name, args=args, analysis_meta_info=meta)

    def make_plot_1d(hists, name, axis, colors=["blue", "orange"], normalize=False, ymin=None, flow=False):
        
        if not isinstance(hists, list):
            hists = [hists]
            name = [name]

        nBinsTotal = hists[0].values(flow=flow).size
        nBins = hists[0].project(axis).values(flow=flow).size
        logger.info(f"Make plot {name} with axes {hists[0].axes.name} and {nBins}/{nBinsTotal} bins")
        # average over bins
        h1ds = [h.project(axis)*nBins/nBinsTotal for h in hists]

        if normalize:
            h1ds = [h/np.sum(h.values(flow=True)) for h in h1ds]

        ymax = max([max(h.values(flow=flow)) for h in h1ds])
        ymin = ymin if ymin is not None else min([min(h.values(flow=flow)) for h in h1ds])
        yrange = ymax - ymin

        fig, ax = plot_tools.figure(h1ds[0], xlabel=labels.get(axis,project[0]), ylabel="a.u.", cms_label="Preliminary", automatic_scale=False, width_scale=1.2,
            ylim=(ymin if ymin == 0 else ymin - yrange*0.2, ymax + yrange*0.2))

        outname = ""
        for i, h1d in enumerate(h1ds):
            x = h1d.axes.edges[0][:-1]
            y = h1d.values(flow=flow)
            err = np.sqrt(h1d.variances(flow=flow))

            if flow:
                # add extra bin with bin wdith of 1% of total width
                rangex = x[-1] - x[0]

                x = np.insert(x, 0, x[0] - rangex*0.02)
                x = np.append(x, x[-1] + rangex*0.02)

            ax.set_xlim((min(x),max(x)))

            ax.step(x, y, color=colors[i], label=name[i], where="post")
            ax.fill_between(x, y - err, y + err, alpha=0.3, color=colors[i], step="post")

            outname += f"_{name[i]}"

        plot_tools.addLegend(ax)

        output_tools.make_plot_dir(*args.plotdir.rsplit("/", 1))
        plot_name = f"hist_{outname}_{axis}_{proc}"
        if args.noSmoothing:
            plot_name += "_noSmoothing"
        plot_tools.save_pdf_and_png(args.plotdir, plot_name)
        plot_tools.write_index_and_log(args.plotdir, plot_name, args=args, analysis_meta_info=meta)

    if not args.noSmoothing:
        make_plot_2d(hnum, args.num)
        make_plot_2d(hden, args.den)

        make_plot_2d(hnum, f"{args.num}_err", plot_error=True)
        make_plot_2d(hden, f"{args.den}_err", plot_error=True)

    make_plot_2d(hratio, f"{args.num}_div_{args.den}", cmin=0, cmax=3)
    make_plot_2d(hratio, f"{args.num}_div_{args.den}_err", plot_error=True)

    if args.logPlots:
        # no log plots for 1d histgrams
        continue

    for ax in project:
        make_plot_1d(hratio, f"{args.num}_div_{args.den}", ax)

        if not args.noSmoothing:
            make_plot_1d([hnum, hden], [args.num, args.den], ax, normalize=True, ymin=0)

outname = args.num.replace('-', '') + 'ew'
with lz4.frame.open(f"{args.outpath}/{outname}CorrZ.pkl.lz4", "wb") as f:
    pickle.dump({
            'Z' : {
                f"{outname}_minnlo_ratio" : corrh['ZToMuMu'],
            },
            "meta_data" : output_tools.metaInfoDict(args=args),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

corrh['W'] = corrh['WplusToMuNu']+corrh['WminusToMuNu']
with lz4.frame.open(f"{args.outpath}/{outname}CorrW.pkl.lz4", "wb") as f:
    pickle.dump({
            'W' : {
                f"{outname}_minnlo_ratio" : corrh['W'],
            },
            "meta_data" : output_tools.metaInfoDict(args=args),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

# logger.info(corrh['ZToMuMu'])