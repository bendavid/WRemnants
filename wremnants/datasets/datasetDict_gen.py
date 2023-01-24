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
                    ["{BASE_PATH}/ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'ZToMuMu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/ZToMuMu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_ZmmPostVFP,
    },
    'WplusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
    },
    'WplusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WplusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_WpmunuPostVFP,
    },
    'WminusToMuNu_horace-lo-photos' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos_nopair-isr-pythia/job*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
    },
    'WminusToMuNu_horace-nlo' : { 
                   'filepaths' :
                    ["{BASE_PATH}/WminusToMuNu_TuneCP5_13TeV-horace-exp-new-fsr-off-isr-off/job*.root"],
                   'xsec' : xsec_WmmunuPostVFP,
    },
}
