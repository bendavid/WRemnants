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

def getDatasets(maxFiles=-1, filt=None):
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
    
    WmmunuPostVFP = narf.Dataset(name = "WminusmunuPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/NanoAOD/TrackRefitv1/WminusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MCPostVFPWeightFix/*/*/*.root"], maxFiles),
            is_data = False,
            xsec = 8562.66,
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
            group = "Top",
    )

    ttbarlqPostVFP = narf.Dataset(name = "TTSemileptonicPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/TTbar_SemiLeptonic_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 365.64,
            group = "Top",
    )

    # TODO: these samples and cross sections are preliminary
    singleTop_schanLepDecaysPostVFP = narf.Dataset(name = "SingleTschanLepDecaysPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/SingleTop_schan_lepDecays_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 3.74,
            group = "Top",
    )

    singleTop_tWAntitopPostVFP = narf.Dataset(name = "SingleTtWAntitopPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/SingleTop_tW_antitop_noFullyHadronic_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 19.55,
            group = "Top",
    )
    singleTop_tchanAntitopPostVFP = narf.Dataset(name = "SingleTtchanAntitopPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/SingleTop_tchan_antitop_inclusive_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 70.79,
            group = "Top",
    )
    singleTop_tchanTopPostVFP = narf.Dataset(name = "SingleTtchanTopPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/SingleTop_tchan_top_inclusive_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 119.71,
            group = "Top",
    )    

    # TODO: should really use the specific decay channels
    wwPostVFP = narf.Dataset(name = "WWPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/WW_inclusive_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 75.8,
            group = "Diboson",
    )

    wzPostVFP = narf.Dataset(name = "WZPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/WZ_inclusive_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 27.6,
            group = "Diboson",
    )

    zz2l2nuPostVFP = narf.Dataset(name = "ZZ2l2nuPostVFP",
        filepaths = makeFilelist(
                ["/scratch/shared/originalNANO/ZZ_2l2nu_postVFP/*.root"], maxFiles),
            is_data = False,
            xsec = 0.564,
            group = "Diboson",
    )

    allPostVFP = [dataPostVFP,
                  WpmunuPostVFP, WmmunuPostVFP, WptaunuPostVFP, WmtaunuPostVFP,
                  ZmmPostVFP, ZttPostVFP,
                  ttbarlnuPostVFP, ttbarlqPostVFP,
                  singleTop_schanLepDecaysPostVFP, singleTop_tWAntitopPostVFP, singleTop_tchanAntitopPostVFP, singleTop_tchanTopPostVFP,
                  wwPostVFP, wzPostVFP, zz2l2nuPostVFP]

    if filt:
        return list(filter(filt, allPostVFP))

    return allPostVFP

def buildXrdFileList(path, xrd):
    xrdpath = path[path.find('/store'):]
    logging.debug(f"Looking for path {xrdpath}")
    f = subprocess.check_output(['xrdfs', f'root://{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
    return filter(lambda x: "root" in x[-4:], f.split())
