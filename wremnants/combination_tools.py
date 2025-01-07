
mass_name_out = "Nominal"
mass_reference_value = 80379.
mass_scaling_factor = 100.

def build_parm_map(alphas = True):
    mass_name = "massShiftW100MeV"

    parm_map = {}
    parm_map[mass_name] = mass_name_out

    for i in range(100):
        # pdf nuisances are numbered starting from 1, and in addition the first one
        # is a spurious nuisance with alternate=nominal due to a (harmless) bug in the
        # treatment of the symmetric hessian sets so the offset of 2 is needed for the
        # parameter name here
        parm_map[f"pdf{i+2}NNPDF31"] = f"PDF_{306000 + 1 + i}"

    if alphas:
        parm_map["pdfAlphaSSymAvg"] = "PDF_alphas"

    return parm_map
