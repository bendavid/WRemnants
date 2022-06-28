import narf
import logging
import subprocess
import glob
import pathlib

lumicsv = f"{pathlib.Path(__file__).parent.parent}/data/bylsoutput.csv"
lumijson = f"{pathlib.Path(__file__).parent.parent}/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1):
    filelist = []
    for path in paths:
        filelist.extend(glob.glob(path) if path[:4] != "/eos" else buildXrdFileList(path, "eoscms.cern.ch"))
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
            ["/scratch/shared/NanoAOD/TrackRefitv2/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 2001.9,
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
            ["/scratch/shared/NanoAOD/TrackRefitv2/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles),
        is_data = False,
        xsec = 11765.9,
    )
    
    WmmunuPostVFP = narf.Dataset(name = "WminusmunuPostVFP",
        filepaths = makeFilelist(
            ["/scratch/shared/NanoAOD/TrackRefitv2/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles),
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

    allPostVFP_gen = allPostVFP[1:]

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
