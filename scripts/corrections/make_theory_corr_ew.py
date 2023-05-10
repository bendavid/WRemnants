import numpy as np
import lz4.frame
import pickle
import logging
from wremnants import plot_tools, theory_corrections, theory_tools
from utilities import boostHistHelpers as hh
from utilities import common, input_tools, output_tools
import hist
import argparse
import os
import h5py
import narf

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default="w_z_gen_dists_scetlib_dyturboCorr_ewinput.hdf5", help="File containing EW hists")
parser.add_argument("--num", type=str, default="horace-nlo", help="Numerator")
parser.add_argument("--den", type=str, default="horace-lo-photos", help="Denominator")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']
charge_dict = {'ZToMuMu': 0, 'WplusToMuNu': 1, 'WminusToMuNu': 0}

# file created with `python WRemnants/scripts/histmakers/w_z_gen_dists.py --skipAngularCoeffs --filter horace -p ewinput`
f = h5py.File(args.input, 'r')
res = narf.ioutils.pickle_load_h5py(f["results"])

corrh = {}

for proc in procs:
    # Make 2D ratio
    print('------')
    hnum = res[f'{proc}_{args.num}']['output']['nominal_ew'].get()
    hden = res[f'{proc}_{args.den}']['output']['nominal_ew'].get()
    print('Integrals', np.sum(hnum.values(flow=True)), np.sum(hden.values(flow=True)))
    hnum = hh.normalize(hnum)
    hden = hh.normalize(hden)
    print('Integrals after normalizing', np.sum(hnum.values(flow=True)), np.sum(hden.values(flow=True)))
    hratio = hh.divideHists(hnum, hden)
    hratio = hh.smoothTowardsOne(hratio)
    scale = np.sum(hden.values(flow=True)) / np.sum(hden.values(flow=True)*hratio.values(flow=True))
    print('Adjustment after smoothing', scale)
    hratio.values(flow=True)[...]    *= scale
    hratio.variances(flow=True)[...] *= scale
    print('Integrals after adjustment', np.sum(hden.values(flow=True)), np.sum(hden.values(flow=True)*hratio.values(flow=True)))

    # Add dummy axis
    axis_dummy = hist.axis.Regular(1, -10., 10., underflow=False, overflow=False, name = "dummy")
    hdummy = hist.Hist(*hratio.axes, axis_dummy, storage=hist.storage.Double())
    hdummy.values(flow=True)[...,0] = hratio.values(flow=True)
    hratio = hdummy

    # Add charge axis
    if proc[0] == 'W':
        axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
    elif proc[0] == 'Z':
        axis_charge = hist.axis.Regular(1, -1., 1., underflow=False, overflow=False, name = "charge")
    hcharge = hist.Hist(*hratio.axes, axis_charge, storage=hist.storage.Double())
    hcharge.values(flow=True)[...,charge_dict[proc]] = hratio.values(flow=True)
    hratio = hcharge

    # Add syst axis
    corrh[proc] = hist.Hist(*hratio.axes, hist.axis.Regular(3, 0, 3, underflow=False, overflow=False, name="systIdx"), storage=hist.storage.Double())
    # Variations: 0=original MiNNLO, 1=Horace NLO, 2=mirrored
    hones = hist.Hist(*hratio.axes, storage=hist.storage.Double())
    hones.values(flow=True)[...,charge_dict[proc]] = np.ones(hdummy.values(flow=True).shape)
    hmirror = hist.Hist(*hratio.axes, storage=hist.storage.Double())
    hmirror.values(flow=True)[...] = 2*hratio.values(flow=True) - hones.values(flow=True)
    mirrorscale = np.sum(hden.values(flow=True)) / np.sum(hden.values(flow=True)*hmirror.values(flow=True)[...,0,charge_dict[proc]])
    # print(f'{proc} mirrorscale = {mirrorscale}')
    hmirror.values(flow=True)[...]    *= mirrorscale
    hmirror.variances(flow=True)[...] *= mirrorscale
    corrh[proc].values(flow=True)[...,0] = hones.values(flow=True)
    corrh[proc].values(flow=True)[...,1] = hratio.values(flow=True)
    corrh[proc].values(flow=True)[...,2] = hmirror.values(flow=True)

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

# print(corrh['ZToMuMu'])