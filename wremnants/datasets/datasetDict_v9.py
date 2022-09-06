import copy

BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
xsec_ZmmPostVFP = 2001.9
xsec_WpmunuPostVFP = 11765.9
xsec_WmmunuPostVFP = 8703.87
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)

dataDictV9 = {
    'dataPostVFP' : { 'name' :  "dataPostVFP",
                      'filepaths' : ["/scratch/shared/NanoAOD/TrackRefitv2/NanoV9DataPostVFP/Run2016F_220627_141813/*/*.root",
                                     "/scratch/shared/NanoAOD/TrackRefitv2/NanoV9DataPostVFP/Run2016G_220627_141950/*/*.root",
                                     "/scratch/shared/NanoAOD/TrackRefitv2/NanoV9DataPostVFP/Run2016H_220627_142357/*/*.root",]
    },
    'ZmmPostVFP' : { 'name' : "ZmumuPostVFP",
                   'filepaths' :
                    ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220729_152749/*/*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'ZttPostVFP' : { 'name' : "ZtautauPostVFP",
                   'filepaths' : 
                   ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220823_073001/000*/*.root"],
                   # At least one tau->e or mu decay, so everything that's not all other decays
                   'xsec' : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO,
    },

    'WpmunuPostVFP' : { 'name' : "WplusmunuPostVFP",
                      'filepaths' : 
                      ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220729_152523/*/*.root"],
                      'xsec' : xsec_WpmunuPostVFP
    },
    'WmmunuPostVFP' : { 'name' : "WminusmunuPostVFP",
                      'filepaths' : 
                      ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220729_152430/*/*.root"],
                      'xsec' : xsec_WmmunuPostVFP
    },

    'WptaunuPostVFP' : { 'name' : "WplustaunuPostVFP",
                         'filepaths' : 
                         ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/WplusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220823_121429/*/*.root"],
                         'xsec' : BR_TAUToMU*xsec_WpmunuPostVFP,
    },
    
    'WmtaunuPostVFP' : { 'name' : "WminustaunuPostVFP",
                         'filepaths' : 
                         ["/scratch/shared/NanoAOD/TrackRefitv2_GenPartPrecision/WminusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP/220823_073247/0000/*.root"],
                         'xsec' : BR_TAUToMU*xsec_WmmunuPostVFP,
    },
    'ttbarlnuPostVFP' : { 'name' : "TTLeptonicPostVFP",
                          'filepaths' : 
                          ["/scratch/shared/NanoAOD/BKGV9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*.root"],
                          'xsec' : 88.29,
                          'group' : "Top",
    },

    'ttbarlqPostVFP' : { 'name' : "TTSemileptonicPostVFP",
                         'filepaths' : 
                         ["/scratch/shared/NanoAOD/BKGV9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*.root"],
                         'xsec' : 365.34,
                         'group' : "Top",
    },
    # TODO: these samples and cross sections are preliminary
    'singleTop_schanLepDecaysPostVFP' : { 'name' : "SingleTschanLepDecaysPostVFP",
                                          'filepaths' : 
                                          ["/scratch/shared/NanoAOD/BKGV9/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/*.root"],
                                          'xsec' : 3.74,
                                          'group' : "Top",
    },
    'singleTop_tWAntitopPostVFP' : { 'name' : "SingleTtWAntitopPostVFP",
                                     'filepaths' : 
                                     ["/scratch/shared/NanoAOD/BKGV9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                     'xsec' : 19.55,
                                     'group' : "Top",
    },
    'singleTop_tchanAntitopPostVFP' : { 'name' : "SingleTtchanAntitopPostVFP",
                                        'filepaths' : 
                                        ["/scratch/shared/NanoAOD/BKGV9/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                        'xsec' : 70.79,
                                        'group' : "Top",
    },
    'singleTop_tchanTopPostVFP' : { 'name' : "SingleTtchanTopPostVFP",
                                    'filepaths' : 
                                    ["/scratch/shared/NanoAOD/BKGV9/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                    'xsec' : 119.71,
                                    'group' : "Top",
    },   
    # TODO: should really use the specific decay channels
    'wwPostVFP' : { 'name' : "WWPostVFP",
                    'filepaths' : 
                    ["/scratch/shared/NanoAOD/BKGV9/WW_TuneCP5_13TeV-pythia8/*.root"],
                    'xsec' : 75.8,
                    'group' : "Diboson",
    },
    'wzPostVFP' : { 'name' : "WZPostVFP",
                    'filepaths' : 
                    ["/scratch/shared/NanoAOD/BKGV9/WZ_TuneCP5_13TeV-pythia8/*.root"],
                    'xsec' : 27.6,
                    'group' : "Diboson",
    },
    'zz2l2nuPostVFP' : { 'name' : "ZZ2l2nuPostVFP",
                         'filepaths' : 
                         ["/scratch/shared/NanoAOD/BKGV9/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/*.root"],
                         'xsec' : 0.564,
                         'group' : "Diboson",
    }
}

###Pisa server
dataDictV9_pisa = copy.deepcopy(dataDictV9)

dataDictV9_pisa['dataPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/Run2016F_220627_141813/*/*.root", "/scratchnvme/wmass/NANOV9/postVFP/Run2016G_220627_141950//*/*.root", "/scratchnvme/wmass/NANOV9/postVFP/Run2016H_220627_142357/*/*.root"]

dataDictV9_pisa['ZmmPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"]
                   
dataDictV9_pisa['ZttPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/TrackRefitv1/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"]
                   
dataDictV9_pisa['WpmunuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"]

dataDictV9_pisa['WmmunuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"]

dataDictV9_pisa['WptaunuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/TrackRefitv1/WplusJetsToTauNu_TauToMu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"]
    
dataDictV9_pisa['WmtaunuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/TrackRefitv1/WminusJetsToTauNu_TauToMu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"]

dataDictV9_pisa['ttbarlnuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*.root"]

dataDictV9_pisa['ttbarlqPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*.root"]

dataDictV9_pisa['singleTop_schanLepDecaysPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/*.root"]

dataDictV9_pisa['singleTop_tWAntitopPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/*.root"]

dataDictV9_pisa['singleTop_tchanAntitopPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"]

dataDictV9_pisa['singleTop_tchanTopPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"]

dataDictV9_pisa['wwPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/WW_TuneCP5_13TeV-pythia8/*.root"]

dataDictV9_pisa['wzPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/WZ_TuneCP5_13TeV-pythia8/*.root"]

dataDictV9_pisa['zz2l2nuPostVFP']['filepaths']=["/scratchnvme/wmass/NANOV9/postVFP/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/*.root"]
