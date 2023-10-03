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
