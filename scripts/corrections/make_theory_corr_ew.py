import numpy as np
from wremnants import theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, logging
from utilities.io_tools import input_tools, output_tools
import hist
import argparse
import os
import h5py
import narf
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs="+", type=str, default=["w_z_gen_dists_scetlib_dyturboCorr_ewinput.hdf5"], help="File containing EW hists")
parser.add_argument("--nums", nargs="+", type=str, default=["horace-nlo"], help="Numerators")
parser.add_argument("--den", type=str, default="horace-lo-photos", help="Denominatos")
parser.add_argument("--normalize", action="store_true", default=False, help="Normalize distributions before computing ratio")
parser.add_argument("--noSmoothing", action="store_true", default=False, help="Disable smoothing of corrections")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--baseName", default="nominal_ew", type=str, help="histogram name")
parser.add_argument("--project", default=["ewMll", "ewLogDeltaM"], nargs="*", type=str, help="axes to project to")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")

args = parser.parse_args()

logger = logging.setup_logger("make_theory_corr_ew", 4 if args.debug else 3)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']
charge_dict = {'ZToMuMu': 0, 'WplusToMuNu': 1, 'WminusToMuNu': 0}

procs_dict = {
    "ZToMuMu": "ZmumuPostVFP",
    "WminusToMuNu": "WminusmunuPostVFP",
    "WplusToMuNu": "WplusmunuPostVFP",
}

project = args.project

# file created with `python WRemnants/scripts/histmakers/w_z_gen_dists.py --skipAngularCoeffs --filter horace -p ewinput`

res, meta, _ = input_tools.read_infile(args.input)

corrh = {}

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

        histo = res[proc_name]['output'][args.baseName].get()
        logger.info(f'Integrals for {name} {np.sum(histo.values(flow=True))}')

        if args.normalize:
            histo = hh.normalize(histo)
            logger.info(f'Integral for {name} after normalizing {np.sum(histo.values(flow=True))}')
        else:
            histo = hh.scaleHist(histo, res[proc_name]["dataset"]["xsec"]*10e6/res[proc_name]['weight_sum'], createNew=False)
            logger.info(f'Integral for {name} after scaling {np.sum(histo.values(flow=True))}')

        return histo
    
    nums = []
    hnums = []
    for num in args.nums:
        hnum = prepare_hist(num)
        if hnum is not None:
            nums.append(num)
            hnums.append(hnum)

    hden = prepare_hist(args.den)

    if hden is None:
        logger.warning(f"Denominator {args.den} for process {proc} not found! Continue with next process.")
        continue

    if len(hnums) < 1:
        logger.warning(f"No nominator found for {args.nums} for process {proc}! Continue with next process.")
        continue

    h2Dratios = []
    corrhs = {}

    def make_correction(h1, h2, name):
        hratio = hh.divideHists(h1, h2)
        if not args.noSmoothing and not len(hratio.axes)>1:
            # 2D smoothing
            hratio = hh.smoothTowardsOne(hratio)
            logger.info('Integrals after smoothing {0} {1}'.format(np.sum(h2.values(flow=True)), np.sum(h2.values(flow=True)*hratio.values(flow=True))))
            if args.normalize:
                scale = np.sum(h2.values(flow=True)) / np.sum(h2.values(flow=True)*hratio.values(flow=True))
                hratio = hh.scaleHist(hratio, scale)
                logger.info('Integrals after adjustment {0} {1}'.format(np.sum(h2.values(flow=True)), np.sum(h2.values(flow=True)*hratio.values(flow=True))))

        if len(hratio.axes)>1:
            h2Dratios.append(hratio)

        # Add charge axis
        if proc[0] == 'W':
            axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
        elif proc[0] == 'Z':
            axis_charge = hist.axis.Regular(1, -1., 1., underflow=False, overflow=False, name = "charge")
        hcharge = hist.Hist(*hratio.axes, axis_charge, storage=hratio._storage_type())
        hcharge.view(flow=True)[...,charge_dict[proc]] = hratio.view(flow=True)

        # Variations: 0=original MiNNLO, 1=Horace NLO, 2=mirrored (doubled)
        hones = hist.Hist(*hcharge.axes, storage=hist.storage.Double())
        hones.values(flow=True)[...,charge_dict[proc]] = np.ones(hratio.values(flow=True).shape)
        hmirror = hist.Hist(*hcharge.axes, storage=hist.storage.Double())
        hmirror.values(flow=True)[...] = 2*hcharge.values(flow=True) - hones.values(flow=True)
        if args.normalize:
            mirrorscale = np.sum(h2.values(flow=True)) / np.sum(h2.values(flow=True)*hmirror.values(flow=True)[...,0,charge_dict[proc]])
            logger.info(f'{proc} mirrorscale = {mirrorscale}')
            hmirror = hh.scaleHist(hmirror, mirrorscale)

        # Add syst axis
        if name not in corrh:
            corrh[name] = {}
        corrh[name][proc] = hist.Hist(*hcharge.axes, hist.axis.Regular(3, 0, 3, underflow=False, overflow=False, name="systIdx"), storage=hist.storage.Double())
        corrh[name][proc].values(flow=True)[...,0] = hones.values(flow=True)
        corrh[name][proc].values(flow=True)[...,1] = hcharge.values(flow=True)
        corrh[name][proc].values(flow=True)[...,2] = hmirror.values(flow=True)

        corrhs[name] = corrh[name]

    for num, hnum in zip(nums, hnums):
       
        make_correction(hnum, hden, f"{num}ew")

        for ax in project:
            hnum1D = hnum.project(ax)
            hden1D = hden.project(ax)

            make_correction(hnum1D, hden1D, f"{num}ew_{ax}")

for num, corrh in corrhs.items():
    outname = num.replace('-', '')
    if 'ZToMuMu' in corrh:
        outfile = f"{args.outpath}/{outname}"
        output_tools.write_theory_corr_hist(outfile, 'Z', {f"{outname}_minnlo_ratio" : corrh['ZToMuMu']}, args)

    if 'WplusToMuNu' in corrh and "WminusToMuNu" in corrh:
        corrh['W'] = corrh['WplusToMuNu']+corrh['WminusToMuNu']
        outfile = f"{args.outpath}/{outname}"
        output_tools.write_theory_corr_hist(outfile, 'Z', {f"{outname}_minnlo_ratio" : corrh['ZToMuMu']}, args)

