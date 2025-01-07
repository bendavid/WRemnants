#!/usr/bin/env python3

import argparse
import h5py
import numpy as np

from wremnants.combination_tools import mass_name_out, mass_reference_value, mass_scaling_factor, build_parm_map

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputFile")
parser.add_argument("-o", "--outputFile")

args = parser.parse_args()

with h5py.File(args.inputFile, "r") as f:
    parms = [parm.decode() for parm in f["parms"]]
    x = f["x"][...]
    cov = f["cov"][...]

#alphas not included for now
parm_map = build_parm_map(alphas=False)

# select the relevant subset of parameters
idxs = [parms.index(key) for (key, value) in parm_map.items() if key in parms]

parms = [parms[idx] for idx in idxs]
parms_out = [parm_map[parm] for parm in parms]
x = x[idxs]
cov = cov[idxs, :]
cov = cov[:, idxs]

mass_idx_out = parms_out.index(mass_name_out)

# scale parameter values and covariance matrix to MeV
x[mass_idx_out] *= mass_scaling_factor
cov[mass_idx_out, :] *= mass_scaling_factor
cov[:, mass_idx_out] *= mass_scaling_factor

errs = np.sqrt(np.diag(cov))

errsi = errs[:, None]
errsj = errs[None, :]
cor = cov/errsi/errsj

#FIXME make this more generic
if mass_idx_out != 0:
    raise ValueError(f"mass_idx_out is {mass_idx_out} but logic currently expects it to be zero")

cornp = cor[1:, 1:]

mass_nominal = mass_reference_value + x[mass_idx_out]
mass_uncertainty_total = errs[mass_idx_out]

impacts = cov[mass_idx_out, :]
impacts[mass_idx_out] = 0.

impacts_total_sq = np.sum(impacts**2)
mass_uncertainty_implicit = np.sqrt(mass_uncertainty_total**2 - impacts_total_sq)

# write to text output

if False:
    # unused test format
    with open(args.outputFile, 'w', encoding="utf-8") as f:
        f.write("# postfit values\n")
        for parm, val in zip(parms_out, x):
            f.write(f"{parm} {val}\n")

        f.write("\n")
        f.write("# postfit covariance\n")
        np.savetxt(f, cov)



out0 = f"{args.outputFile}.txt"
out1 = f"{args.outputFile}_NPs.txt"

# format used by ATLAS for now
with open(out0, 'w', encoding="utf-8") as f:
    f.write(f"Nominal {mass_idx_out} {mass_uncertainty_implicit}\n")
    for parm, val, impact in zip(parms_out, x, impacts):
        if parm != mass_name_out:
            f.write(f"{parm} {impact} {-impact}\n")

    f.write(f"Total {mass_uncertainty_total} {-mass_uncertainty_total}\n")
    f.write(f"ChiSquare {0.} {0.}\n")

with open(out1, 'w', encoding="utf-8") as f:
    for parm, val, err in zip(parms_out, x, errs):
        if parm != mass_name_out:
            f.write(f"{parm} {val} {err}\n")
    f.write("\n")
    f.write("CORRELATION MATRIX\n")
    np.savetxt(f, cornp)
