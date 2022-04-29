import narf
import logging
import subprocess
import glob
import pathlib

lumicsv = f"{pathlib.Path(__file__).parent.parent}/data/bylsoutput.csv"
lumijson = f"{pathlib.Path(__file__).parent.parent}/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1, eos=False):
    filelist = []
    for path in paths:
        filelist.extend(glob.glob(path) if not eos else buildXrdFileList(path, "eoscms.cern.ch"))
    return filelist if maxFiles < 0 else filelist[:maxFiles]

def getDatasets(maxFiles=-1, filt=None, mode=None):
    dataPostVFP = narf.Dataset(name = "dataPostVFP",
        filepaths = makeFilelist(["/scratch/shared/NanoAOD/TrackRefitv1/SingleMuon/Run2016F_postVFP_220223_222034/*/*.root",
            "/scratch/shared/NanoAOD/TrackRefitv1/SingleMuon/Run2016G_220223_222128/*/*.root",
            "/scratch/shared/NanoAOD/TrackRefitv1/SingleMuon/Run2016H_220223_222223/*/*.root",
        ], maxFiles),
        is_data = True,
        lumi_csv = lumicsv,
        lumi_json = lumijson
    )

    ZmmPostVFP = narf.Dataset(name = "ZmumuPostVFP",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoAOD/TrackRefitv1/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 1976.1,
    )

    Zmm_bugfix = narf.Dataset(name = "Zmumu_bugfix",
        filepaths = makeFilelist(
            ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/DYJetsToMuMu_svn3900_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer15wmLHEGS/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 2001.9,
    )

    Zmm_bugfix_newprod = narf.Dataset(name = "Zmumu_newprod",
        filepaths = makeFilelist(
            ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/DYJetsToMuMu_svn3900_newprod_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 2001.9,
    )

    Zmm_bugfix_slc7 = narf.Dataset(name = "Zmumu_bugfix_slc7",
        filepaths = makeFilelist(
            ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/DYJetsToMuMu_svn3900_slc7_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer15wmLHEGS/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 1976.1,
    )

    Zmm_bugfix = narf.Dataset(name = "Zmumu_bugfix",
        filepaths = makeFilelist(
            ["/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/DYJetsToMuMu_svn3900_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer15wmLHEGS/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 2001.1,
    )

    BR_TAUToMU = 0.1739
    BR_TAUToE = 0.1782
    ZttPostVFP = narf.Dataset(name = "ZtautauPostVFP",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoAOD/TrackRefitv1/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
        is_data = False,
        # At least one tau->e or mu decay, so everything that's not all other decays
        xsec = ZmmPostVFP.xsec*(1.-(1. - BR_TAUToMU - BR_TAUToE)**2),
    )

    WpmunuPostVFP = narf.Dataset(name = "WplusmunuPostVFP",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoAOD/TrackRefitv1/WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 11572.19,
    )
    
    Wpmunu_bugfix = narf.Dataset(name = "Wplusmunu_bugfix",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoGen/WplusToMuNu_svn3900_slc7_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/*/*.root"], maxFiles),
        is_data = False,
        xsec = 11765.9,
    )

    Wpmunu_bugfix_reweight_h2 = narf.Dataset(name = "Wplusmunu_bugfix_reweight_h2",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoGen/WplusToMuNu_svn3900_slc7_ReweightToBugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/220413_121251/000*/*.root"], maxFiles),
        is_data = False,
        xsec = 11572.19,
    )

    WmmunuPostVFP = narf.Dataset(name = "WminusmunuPostVFP",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoAOD/TrackRefitv1/WminusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
            is_data = False,
            xsec = 8562.66,
    )

    Wmmunu_bugfix = narf.Dataset(name = "Wminusmunu_bugfix",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoGen/WminusToMuNu_svn3900_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/*/*.root", 
                "/scratch/shared/NanoGen/WminusToMuNu_svn3900_slc7_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/*/*.root"], maxFiles),
            is_data = False,
            xsec = 8703.87,
    )

    Wmmunu_bugfix_newprod = narf.Dataset(name = "Wminusmunu_bugfix_newprod",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoGen//WminusToMuNu_svn3900_newprod_BugFix_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/220421_232301/*/*.root"], maxFiles),
            is_data = False,
            xsec = 8703.87,
    )

    WptaunuPostVFP = narf.Dataset(name = "WplustaunuPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/NanoAOD/TrackRefitv1/WplusJetsToTauNu_TauToMu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
            is_data = False,
            xsec = BR_TAUToMU*WpmunuPostVFP.xsec,
    )

    WmtaunuPostVFP = narf.Dataset(name = "WminustaunuPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/NanoAOD/TrackRefitv1/WminusJetsToTauNu_TauToMu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
            is_data = False,
            xsec = BR_TAUToMU*WmmunuPostVFP.xsec,
    )

    ttbarlnuPostVFP = narf.Dataset(name = "TTLeptonicPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/TTbar_2l2nu_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 88.29,
    )

    ttbarlqPostVFP = narf.Dataset(name = "TTSemileptonicPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/TTbar_SemiLeptonic_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 365.64,
    )

    # TODO: should really use the specific decay channels
    wwPostVFP = narf.Dataset(name = "WWPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/WW_inclusive_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 75.8,
    )

    allPostVFP = [
        dataPostVFP, 
        WpmunuPostVFP,
        WmmunuPostVFP,
        WptaunuPostVFP,
        WmtaunuPostVFP,
        ZmmPostVFP,
        ZttPostVFP,
        ttbarlnuPostVFP,
        ttbarlqPostVFP,
        wwPostVFP
    ]

    allPostVFP_gen = allPostVFP[1:]
    allPostVFP_gen.extend([
        Zmm_bugfix_slc7,
        Zmm_bugfix_newprod,
        Wmmunu_bugfix,
        Wmmunu_bugfix_newprod,
        Wpmunu_bugfix,
        Wpmunu_bugfix_reweight_h2,
    ])

    samples = allPostVFP if mode != "gen" else allPostVFP_gen
    if filt:
        return list(filter(filt, samples))
    else:
        return samples

def buildXrdFileList(path, xrd):
    xrdpath = path[path.find('/store'):]
    logging.debug(f"Looking for path {xrdpath}")
    f = subprocess.check_output(['xrdfs', f'root://{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
    return filter(lambda x: "root" in x[-4:], f.split())
