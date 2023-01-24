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

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default="w_z_gen_dists_RawPFMET.pkl.lz4", help="File containing EW hists")
parser.add_argument("--num", type=str, default="horace-nlo", help="Numerator")
parser.add_argument("--den", type=str, default="horace-lo-photos", help="Denominator")
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

procs = ['ZToMuMu', 'WplusToMuNu', 'WminusToMuNu']
charge_dict = {'ZToMuMu': 0, 'WplusToMuNu': +1, 'WminusToMuNu': -1}

# file created with `python WRemnants/scripts/histmakers/w_z_gen_dists.py --skipAngularCoeffs --filter horace --ewHists`
with lz4.frame.open(args.input) as f:
    res = pickle.load(f)

corrh = {}

for proc in procs:
    # Make 2D ratio
    hnum = hh.normalize(res[f'{proc}_{args.num}']['output']['nominal_ew'])
    hden = hh.normalize(res[f'{proc}_{args.den}']['output']['nominal_ew'])
    hratio = hh.divideHists(hnum, hden)
    hratio = hh.smoothenTowardsOne(hratio)

    # Add dummy axis
    axis_dummy = hist.axis.Regular(1, -10., 10., underflow=False, overflow=False, name = "dummy")
    hdummy = hist.Hist(*hratio.axes, axis_dummy, storage=hist.storage.Double())
    hdummy.values()[...,0] = hratio.values()
    hratio = hdummy

    # Add charge axis
    if proc[0] == 'W':
        axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
    elif proc[0] == 'Z':
        axis_charge = hist.axis.Regular(1, -1., 1., underflow=False, overflow=False, name = "charge")
    hcharge = hist.Hist(*hratio.axes, axis_charge, storage=hist.storage.Double())
    hcharge.values()[...,charge_dict[proc]] = hratio.values()
    hratio = hcharge

    # Add syst axis
    corrh[proc] = hist.Hist(*hratio.axes, hist.axis.Regular(3, 0, 3, name="systIdx"), storage=hist.storage.Double())
    # Variations: 0=original MiNNLO, 1=Horace NLO, 2=mirrored
    horig = hist.Hist(*hratio.axes, storage=hist.storage.Double())
    horig.values()[...] = np.ones(hratio.shape)
    mirror = hh.mirrorHist(horig, hratio, cutoff=1e-5)
    corrh[proc].values()[...,0] = horig.values()
    corrh[proc].values()[...,1] = hratio.values()
    corrh[proc].values()[...,2] = mirror.values()

outname = args.num.replace('-', '') + 'ew'
with lz4.frame.open(f"{args.outpath}/{outname}CorrZ.pkl.lz4", "wb") as f:
    pickle.dump({
            'Z' : {
                f"{outname}_minnlo_ratio" : corrh['ZToMuMu'],
            },
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

with lz4.frame.open(f"{args.outpath}/{outname}CorrW.pkl.lz4", "wb") as f:
    pickle.dump({
            'W' : {
                f"{outname}_minnlo_ratio" : corrh['WplusToMuNu']+corrh['WminusToMuNu'],
            },
            "meta_data" : output_tools.metaInfoDict(),
        }, f, protocol = pickle.HIGHEST_PROTOCOL)

print(corrh['ZToMuMu'])