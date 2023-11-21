
process_colors = {
    "Data": "black",
    "Zmumu": "lightblue",
    "Zee": "lightblue",
    "Ztautau": "darkblue",
    "Wmunu": "darkred",
    "Wenu": "darkred",
    "Wtaunu": "orange",
    "DYlowMass": "deepskyblue",
    "Top": "green",
    "Diboson": "pink",
    "QCD": "grey",
    "Other": "grey",
    "Fake": "grey",
    "Fake_e": "grey",
    "Fake_mu": "grey",
}

process_labels = {
    "Data": "Data",
    "Zmumu": r"Z$\to\mu\mu$",
    "Zee": r"Z$\to ee$",
    "Ztautau": r"Z$\to\tau\tau$",
    "Wmunu":  r"W$^{\pm}\to\mu\nu$",
    "Wenu": r"W$^{\pm}\to e\nu$",
    "Wtaunu": r"W$^{\pm}\to\tau\nu$",
    "DYlowMass": r"Z$\to\mu\mu$, $10<m<50$ GeV",
    "Top": "Top",
    "Diboson": "Diboson",
    "QCD": "QCD MC",
    "Other": "Other",
    "Fake": "Nonprompt",
    "Fake_e": "Nonprompt (e)",
    "Fake_mu": r"Nonprompt (\mu)",
}

xlabels = {
    "pt" : r"p$_{T}^{\ell}$ (GeV)",
    "ptW" : r"p$_{T}^{\ell+p_{\mathrm{T}}^{miss}}$ (GeV)",
    "eta" : r"$\eta^{\ell}$",
    "ptll" : r"p$_{\mathrm{T}}^{\ell\ell}$ (GeV)",
    "yll" : r"y$^{\ell\ell}$",
    "mll" : r"m$_{\ell\ell}$ (GeV)",
    "ewMll" : r"m$^{\mathrm{EW}}_{\ell\ell}$ (GeV)",
    "costhetastarll" : r"$\cos{\phi^{\star}_{\ell\ell}}$",
    "phistarll" : r"$\phi^{\star}_{\ell\ell}$",
    "MET_pt" : r"p$_{\mathrm{T}}^{miss}$ (GeV)",
    "MET" : r"p$_{\mathrm{T}}^{miss}$ (GeV)",
    "met" : r"p$_{\mathrm{T}}^{miss}$ (GeV)",
    "mt" : r"m$_{T}^{\ell\nu}$ (GeV)",
    "etaPlus" : r"$\eta^{\ell(+)}$",
    "etaMinus" : r"$\eta^{\ell(-)}$",
    "ptPlus" : r"p$_{\mathrm{T}}^{\ell(+)}$ (GeV)",
    "ptMinus" : r"p$_{\mathrm{T}}^{\ell(-)}$ (GeV)",
    "etaSum":r"$\eta^{\ell(+)} + \eta^{\ell(-)}$",
    "etaDiff":r"$\eta^{\ell(+)} - \eta^{\ell(-)}$",
    "massVgen": "massVgen",
    "ewMll": "ewMll",
    "ewMlly": "ewMlly",
    "ewLogDeltaM": "ewLogDeltaM",
    # add 2d unrolled plots 
    "pt-eta" : r"(p$_{T}^{\ell}$, $\eta^{\ell}$) bin",
    "ptll-yll":r"p$_{\mathrm{T}}^{\ell\ell}$, y$^{\ell\ell}$ bin",
    "mll-yll":r"m$_{\ell\ell}$, y$^{\ell\ell}$ bin",
    "mll-ptll":r"m$_{\ell\ell}$, p$_{\mathrm{T}}^{\ell\ell}$ bin",
    "mll-etaPlus":r"m$_{\ell\ell}$, $\eta^{\ell(+)}$ bin",
    "mll-etaMinus":r"m$_{\ell\ell}$, $\eta^{\ell(-)}$ bin",
    "etaPlus-etaMinus":r"$\eta^{\ell(+)}$, $\eta^{\ell(-)}$ bin",
    "etaSum-etaDiff":r"$\eta^{\ell(+)} + \eta^{\ell(-)}$, $\eta^{\ell(+)} - \eta^{\ell(-)}$ bin",
    # add 3d unrolled plots 
    "mll-etaPlus-etaMinus":r"m$_{\ell\ell}$, $\eta^{\ell(+)}$, $\eta^{\ell(-)}$ bin",
    "mll-etaSum-etaDiff":r"m$_{\ell\ell}$, $\eta^{\ell(+)} + \eta^{\ell(-)}$, $\eta^{\ell(+)} - \eta^{\ell(-)}$ bin",
}

# uncertainties
common_groups = [
    "Total",
    "stat",
    "binByBinStat",
    "luminosity",
    "CMS_recoil",
    "CMS_background",
    "theory_ew",
    "normXsecW"
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
    "ptgen": r"$\mathrm{Pre\ FSR}\ p_\mathrm{T}^{\ell}$",
    "etagen": r"$\mathrm{Pre\ FSR}\ \eta^{\ell}$", 
    "ptVgen": r"$\mathrm{Pre\ FSR}\ p_\mathrm{T}^{\ell\ell}$",
    "massVgen": r"$\mathrm{Pre\ FSR}\ m^{\ell\ell}$", 
    "qT" : r"$\mathrm{Pre\ FSR}\ p_\mathrm{T}^{\ell\ell}$",
    "Q" : r"$\mathrm{Pre\ FSR}\ m^{\ell\ell}$", 
    "absY" : r"$\mathrm{Pre\ FSR}\ Y^{\ell\ell}$",
    "charge" : r"$\mathrm{Pre\ FSR\ charge}$", 
}

syst_labels = {
    "horacenloew" : {0: "nominal", 1: "horace EW NLO/LO", 2: "horace EW NLO/LO doubled", },
    "matrix_radish" : "MATRIX+RadISH",
    "virtual_ew" : {
        0: r"NLOEW + HOEW, CMS, ($G_\mu, m_\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme",
        1: r"NLOEW + HOEW, PS, ($G_\mu, m_\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme", 
        2: r"NLOEW + HOEW, CMS, ($\alpha(m_\mathrm{Z}),m _\mathrm{Z}, \mathrm{sin}^2\Theta_\mathrm{eff}$) scheme", }
}

syst_labels["virtual_ew_wlike"] = syst_labels["virtual_ew"]
