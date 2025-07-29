import h5py
import numpy as np

mass_name_out = "Nominal"
mass_reference_value = 80379.
mass_scaling_factor = 100.

def build_parm_map(alphas = True, asym = False):
    mass_name = "massShiftW100MeV"

    parm_map = {}
    parm_map[mass_name] = mass_name_out

    idxstart = 306000
    for i in range(100):
        # pdf nuisances are numbered starting from 1, and in addition the first one
        # is a spurious nuisance with alternate=nominal due to a (harmless) bug in the
        # treatment of the symmetric hessian sets so the offset of 2 is needed for the
        # parameter name here
        parm_map[f"pdf{i+2}NNPDF31"] = f"PDF_{idxstart + 1 + i}"

    idxstart = 331600
    for i in range(50):
        # pdf nuisances are numbered starting from 1, and in addition the first one
        # is a spurious nuisance with alternate=nominal due to a (harmless) bug in the
        # treatment of the symmetric hessian sets so the offset of 2 is needed for the
        # parameter name here
        parm_map[f"pdf{i+2}NNPDF40"] = f"PDF_{idxstart + 1 + i}"

    idxstart = 93300
    for i in range(40):
        # pdf nuisances are numbered starting from 1, and in addition the first one
        # is a spurious nuisance with alternate=nominal due to a (harmless) bug in the
        # treatment of the symmetric hessian sets so the offset of 2 is needed for the
        # parameter name here
        parm_map[f"pdf{i+2}PDF4LHC21"] = f"PDF_{idxstart + 1 + i}"

    idxstart = 14000
    for i in range(29):
        # pdf nuisance are numbered starting from 1
        # and are ordered up, down, up, down, ....
        parm_map[f"pdf{i+1}CT18SymAvg"] = f"PDF_{idxstart + 1 + 2*i}"

        if asym:
            parm_map[f"pdf{i+1}CT18SymDiff"] = f"PDF_{idxstart + 1 + 2*i}_diff"

    idxstart = 14100
    for i in range(29):
        # pdf nuisance are numbered starting from 1
        # and are ordered up, down, up, down, ....
        parm_map[f"pdf{i+1}CT18ZSymAvg"] = f"PDF_{idxstart + 1 + 2*i}"

        if asym:
            parm_map[f"pdf{i+1}CT18ZSymDiff"] = f"PDF_{idxstart + 1 + 2*i}_diff"

    idxstart = 27400
    for i in range(32):
        # pdf nuisance are numbered starting from 1
        # and are ordered up, down, up, down, ....
        parm_map[f"pdf{i+1}MSHT20SymAvg"] = f"PDF_{idxstart + 1 + 2*i}"

        if asym:
            parm_map[f"pdf{i+1}MSHT20SymDiff"] = f"PDF_{idxstart + 1 + 2*i}_diff"

    idxstart = 29100
    for i in range(52):
        # pdf nuisance are numbered starting from 1
        # and are ordered up, down, up, down, ....
        parm_map[f"pdf{i+1}MSHT20an3loSymAvg"] = f"PDF_{idxstart + 1 + 2*i}"

        if asym:
            parm_map[f"pdf{i+1}MSHT20an3loSymDiff"] = f"PDF_{idxstart + 1 + 2*i}_diff"

    if alphas:
        parm_map["pdfAlphaSSymAvg"] = "PDF_alphas"

        if asymm:
            parm_map["pdfAlphaSSymDiff"] = "PDF_alphas_diff"

    return parm_map


def makeTxtresults(inputFile, outputFile, asym=False):
    with h5py.File(inputFile, "r") as f:
        parms = [parm.decode() for parm in f["parms"]]
        x = f["x"][...]
        cov = f["cov"][...]

    #alphas not included for now
    parm_map = build_parm_map(alphas=False, asym=asym)

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
        with open(outputFile, 'w', encoding="utf-8") as f:
            f.write("# postfit values\n")
            for parm, val in zip(parms_out, x):
                f.write(f"{parm} {val}\n")

            f.write("\n")
            f.write("# postfit covariance\n")
            np.savetxt(f, cov)



    out0 = f"{outputFile}.txt"
    out1 = f"{outputFile}_NPs.txt"

    # format used by ATLAS for now
    with open(out0, 'w', encoding="utf-8") as f:
        f.write(f"Nominal {mass_nominal} {mass_uncertainty_implicit}\n")
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
