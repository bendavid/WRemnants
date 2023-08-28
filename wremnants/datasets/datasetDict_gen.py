from wremnants.datasets.datasetDict_v9 import xsec_ZmmPostVFP,xsec_WpmunuPostVFP,xsec_WmmunuPostVFP

horace_v1 = False
horace_v2 = False
horace_v3 = False

genDataDict = {
    'ZmumuMiNLO' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_TuneCP5_13TeV-powheg-NNLOPS-pythia8-photos/RunIISummer15wmLHEGS/221121_114507/000*/*.root"],
                   'xsec' : 1863.,
                   'group': "Zmumu",
    },
    'ZmumuNNLOPS' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_TuneCP5_13TeV-powheg-NNLOPS-pythia8-photos/RunIISummer15wmLHEGS/221121_114507/000*/*.root"],
                   'xsec' : 1863.,
                   'group': "Zmumu",
    },
    'ZToMuMu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_LO_TuneCP5_13TeV-horace-pythia8-photospp/*/*.root"],
                   'xsec' : xsec_ZmmPostVFP,
                   'group': "Zmumu",
    },
    'ZToMuMu_horace-lo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsJetsToMuMu_LO_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_ZmmPostVFP,
                   'group': "Zmumu",
    },
    'ZToMuMu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_NLOEW_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_ZmmPostVFP,
                   'group': "Zmumu",
    },
    'WplusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_LO_TuneCP5_13TeV-horace-pythia8-photospp/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WplusToMuNu_horace-lo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_LO_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WplusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_NLOEW_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WplusToMuNu_winhac-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_LO_TuneCP5_13TeV-winhac-pythia8-photospp/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WplusToMuNu_winhac-lo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_LO_TuneCP5_13TeV-winhac-pythia8/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WplusToMuNu_winhac-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusJetsToMuNu_NLOEW_TuneCP5_13TeV-winhac-pythia8/*/*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_LO_TuneCP5_13TeV-horace-pythia8-photospp/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_horace-lo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_LO_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_NLOEW_TuneCP5_13TeV-horace-pythia8/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_winhac-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_LO_TuneCP5_13TeV-winhac-pythia8-photospp/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_winhac-lo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_LO_TuneCP5_13TeV-winhac-pythia8/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
    'WminusToMuNu_winhac-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusJetsToMuNu_NLOEW_TuneCP5_13TeV-winhac-pythia8/*/*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
                    'group': "Wmunu",
    },
}

if horace_v1:
    genDataDict.update({
        'ZToMuMu_horace-v1-alpha-old-fsr-off-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-alpha-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-born-fsr-photos-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-born-fsr-photoslow-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photoslow-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-born-fsr-photosnopair-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photosnopair-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-born-fsr-pythia-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-pythia-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-exp-fsr-off-isr-off' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v1-exp-old-fsr-off-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v1/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
    })

if horace_v2:
    genDataDict.update({
        'ZToMuMu_horace-v2-born-fsr-photos_nopair-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v2/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v2-exp-new-fsr-off-isr-off' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v2/ZToMuMu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v2-exp-old-fsr-off-isr-pythia' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v2/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
    })

if horace_v3:
    genDataDict.update({
        'ZToMuMu_horace-v3-lo-photos' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v3-qed' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'ZToMuMu_horace-v3-nlo' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                'xsec' : xsec_ZmmPostVFP,
                'group': "Zmumu",
        },
        'WplusToMuNu_horace-v3-lo-photos' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                'xsec' : xsec_WpmunuPostVFP,
                'group': "Wmunu",
        },
        'WplusToMuNu_horace-v3-qed' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_WpmunuPostVFP,
                'group': "Wmunu",
        },
        'WplusToMuNu_horace-v3-nlo' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                'xsec' : xsec_WpmunuPostVFP,
                'group': "Wmunu",
        },
        'WminusToMuNu_horace-v3-lo-photos' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                'xsec' : xsec_WmmunuPostVFP,
                'group': "Wmunu",
        },
        'WminusToMuNu_horace-v3-qed' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                'xsec' : xsec_WmmunuPostVFP,
                'group': "Wmunu",
        },
        'WminusToMuNu_horace-v3-nlo' : { 
                'filepaths' :
                ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                'xsec' : xsec_WmmunuPostVFP,
                'group': "Wmunu",
        }
    })
