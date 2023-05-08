from wremnants.datasets.datasetDict_v9 import xsec_ZmmPostVFP,xsec_WpmunuPostVFP,xsec_WmmunuPostVFP

genDataDict = {
    'ZmumuMiNLO' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_TuneCP5_13TeV-powheg-NNLOPS-pythia8-photos/RunIISummer15wmLHEGS/221121_114507/000*/*.root"],
                   'xsec' : 1863.,
    },
    'ZmumuNNLOPS' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_TuneCP5_13TeV-powheg-NNLOPS-pythia8-photos/RunIISummer15wmLHEGS/221121_114507/000*/*.root"],
                   'xsec' : 1863.,
    },
    'ZToMuMu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'ZToMuMu_horace-qed' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'ZToMuMu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/ZToMuMu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'WplusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
    },
    'WplusToMuNu_horace-qed' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
    },
    'WplusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WplusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
    },
    'WminusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
    },
    'WminusToMuNu_horace-qed' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/job*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
    },
    'WminusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/Horace_v3/WminusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
    },
}
