common_groups = [
    "Total",
    "stat",
    "binByBinStat",
    "luminosity",
    "CMS_recoil",
    "CMS_background",
    "theory_ew",
]
nuisance_groupings = {
    "max": common_groups + [
        "massShift",
        "QCDscale", 
        "pdfMSHT20",
        "resum",
        "muon_eff_syst",
        "muon_eff_stat",
        "prefire",
        "muonCalibration",
    ],
    "min": common_groups + [
        "massShiftW", "massShiftZ",
        "QCDscalePtChargeMiNNLO", "QCDscaleZPtChargeMiNNLO", "QCDscaleWPtChargeMiNNLO", "QCDscaleZPtHelicityMiNNLO", "QCDscaleWPtHelicityMiNNLO", "QCDscaleZPtChargeHelicityMiNNLO", "QCDscaleWPtChargeHelicityMiNNLO",
        "pdfMSHT20NoAlphaS", "pdfMSHT20AlphaS",
        "resumTNP", "resumNonpert", "resumTransition", "resumScale",
        "muon_eff_stat_reco", "muon_eff_stat_trigger", "muon_eff_stat_iso", "muon_eff_stat_idip",
        "muon_eff_syst_reco", "muon_eff_syst_trigger", "muon_eff_syst_iso", "muon_eff_syst_idip",
        "muonPrefire", "ecalPrefire",
        "nonClosure", "resolutionCrctn",
    ]
}


text_dict = {
    "Zmumu": r"$\mathrm{Z}\rightarrow\mu\mu$",
    "ZToMuMu": r"$\mathrm{Z}\rightarrow\mu\mu$",
    "WplusToMuNu": r"$\mathrm{W}^+\rightarrow\mu\nu$",
    "WminusToMuNu": r"$\mathrm{W}^-\rightarrow\mu\nu$"
}

axis_labels = {
    "ewPTll": r"$\mathrm{Post\ FSR}\ p_\mathrm{T}^{\ell\ell}$",
    "ewMll": r"$\mathrm{Post\ FSR}\ m^{\ell\ell}$", 
    "ptVgen": r"$\mathrm{Pre\ FSR}\ p_\mathrm{T}^{\ell\ell}$",
    "massVgen": r"$\mathrm{Pre\ FSR}\ m^{\ell\ell}$", 
    "qT" : r"$\mathrm{Pre\ FSR}\ p_\mathrm{T}^{\ell\ell}$",
    "Q" : r"$\mathrm{Pre\ FSR}\ m^{\ell\ell}$", 
    "absY" : r"$\mathrm{Pre\ FSR}\ Y^{\ell\ell}$",
    "charge" : r"$\mathrm{Pre\ FSR\ charge}$", 
}

syst_labels = {
    "horacenloew" : {"0": "nominal", "1": "horace EW NLO/LO", "2": "horace EW NLO/LO doubled", },
    "virtual_ew" : {
        "0": r"NLOEW + HOEW, CMS, ($G_\mu, m_\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme",
        "1": r"NLOEW + HOEW, PS, ($G_\mu, m_\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme", 
        "2": r"NLOEW + HOEW, CMS, ($\alpha(m_\mathrm{Z}),m _\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme", }
}

syst_labels["virtual_ew_wlike"] = syst_labels["virtual_ew"]
