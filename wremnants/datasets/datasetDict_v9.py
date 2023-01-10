import copy

BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
xsec_ZmmPostVFP = 2001.9
xsec_WpmunuPostVFP = 11765.9
xsec_WmmunuPostVFP = 8703.87
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)

dataDictV9 = {
    'dataPostVFP' : { 
                      'filepaths' : ["{BASE_PATH}/SingleMuon/NanoV9Run2016FDataPostVFP_{NANO_PROD_TAG}/*/*/*.root",
                                     "{BASE_PATH}/SingleMuon/NanoV9Run2016GDataPostVFP_{NANO_PROD_TAG}/*/*/*.root",
                                     "{BASE_PATH}/SingleMuon/NanoV9Run2016HDataPostVFP_{NANO_PROD_TAG}/*/*/*.root",
                                     ],
    },
    'ZmumuPostVFP' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'ZtautauPostVFP' : { 
                   'filepaths' : 
                   ["{BASE_PATH}/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                   # At least one tau->e or mu decay, so everything that's not all other decays
                   'xsec' : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO,
    },

    'WplusmunuPostVFP' : { 
                      'filepaths' : 
                      ["{BASE_PATH}/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                      'xsec' : xsec_WpmunuPostVFP
    },
    'WminusmunuPostVFP' : { 
                      'filepaths' : 
                      ["{BASE_PATH}/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                      'xsec' : xsec_WmmunuPostVFP
    },

    'WplustaunuPostVFP' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/WplusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                         'xsec' : BR_TAUToMU*xsec_WpmunuPostVFP,
    },
    
    'WminustaunuPostVFP' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/WminusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}/*/*/*.root"],
                         'xsec' : BR_TAUToMU*xsec_WmmunuPostVFP,
    },
    'TTLeptonicPostVFP' : { 
                          'filepaths' : 
                          ["{BASE_PATH}/BKGV9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*.root"],
                          'xsec' : 88.29,
                          'group' : "Top",
    },

    'TTSemileptonicPostVFP' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/BKGV9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*.root"],
                         'xsec' : 365.34,
                         'group' : "Top",
    },
    # TODO: these samples and cross sections are preliminary
    'SingleTschanLepDecaysPostVFP' : { 
                                          'filepaths' : 
                                          ["{BASE_PATH}/BKGV9/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/*.root"],
                                          'xsec' : 3.74,
                                          'group' : "Top",
    },
    'SingleTtWAntitopPostVFP' : { 
                                     'filepaths' : 
                                     ["{BASE_PATH}/BKGV9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                     'xsec' : 19.55,
                                     'group' : "Top",
    },
    'SingleTtchanAntitopPostVFP' : { 
                                        'filepaths' : 
                                        ["{BASE_PATH}/BKGV9/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                        'xsec' : 70.79,
                                        'group' : "Top",
    },
    'SingleTtchanTopPostVFP' : { 
                                    'filepaths' : 
                                    ["{BASE_PATH}/BKGV9/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                    'xsec' : 119.71,
                                    'group' : "Top",
    },   
    # TODO: should really use the specific decay channels
    'WWPostVFP' : { 
                    'filepaths' : 
                    ["{BASE_PATH}/BKGV9/WW_TuneCP5_13TeV-pythia8/*.root"],
                    'xsec' : 75.8,
                    'group' : "Diboson",
    },
    'WZPostVFP' : { 
                    'filepaths' : 
                    ["{BASE_PATH}/BKGV9/WZ_TuneCP5_13TeV-pythia8/*.root"],
                    'xsec' : 27.6,
                    'group' : "Diboson",
    },
    'ZZ2l2nuPostVFP' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/BKGV9/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/*.root"],
                         'xsec' : 0.564,
                         'group' : "Diboson",
    }
}
