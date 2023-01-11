import narf
import os
import pathlib
import socket
import logging
import glob
from wremnants.datasets import datasets2016

logger = logging.getLogger("wremnants").getChild(__name__.split(".")[-1])

lumijson = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt"
lumicsv_mu = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/bylsoutput_HLT_HIMu17_Full.csv"
lumicsv_el = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/bylsoutput_HLT_HIEle20_Full.csv"

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
def getDatasets(maxFiles=-1, filt=None, flavor="",base_path=None):

    if not base_path:
        hostname = socket.gethostname()
        if hostname == "lxplus8s10.cern.ch":
            base_path = "/scratch/shared/lowPU/"
        elif "mit.edu" in hostname:
            base_path = "/data/submit/cms/store/wmass/lowPU/"

    logger.info(f"Loading samples from {base_path}.")

    BR_W_LEP = 3*0.1086 # PDG

    allProcs = [
    
        narf.Dataset(
            name="TTTo2L2Nu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=87.31483776,
            is_data=False,
        ),
        narf.Dataset(
            name="TTToSemiLeptonic",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=364.35,
            is_data=False,
        ),
        narf.Dataset(
            name="TTToHadronic",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=380.10,
            is_data=False,
        ),

        narf.Dataset(
            name="WWTo2L2Nu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=118.7*BR_W_LEP*BR_W_LEP,
            is_data=False,
        ),
        narf.Dataset(
            name="WZTo3LNu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=4.912,
            is_data=False,
        ),
        narf.Dataset(
            name="ZZ",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/ZZ_TuneCP5_13TeV-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=16.523,
            is_data=False,
        ),


        narf.Dataset(
            name="Zmumu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToMuMu_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=2025.74, # 1976.1
            is_data=False,
        ),
        narf.Dataset(
            name="Zee",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToEE_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=2025.74,
            is_data=False,
        ),
        
        narf.Dataset(
            name="Ztautau",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToTauTau_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=2025.74,
            is_data=False,
        ),
        
        
        narf.Dataset(
            name="WminusJetsToMuNu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=8677.3, # 8562.66
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToENu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToENu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=8677.3,
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToTauNu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToTauNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=8677.3,
            is_data=False,
        ),
        
        
        narf.Dataset(
            name="WplusJetsToMuNu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=11811.4, # 11572.19 
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToENu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToENu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=11572.19, # 
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToTauNu",
            filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToTauNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=11572.19,
            is_data=False,
        ),
    ]

    if flavor == "mu" or flavor == "mumu":
        allProcs.append(
        narf.Dataset(
                filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/SingleMuon/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                name="singlemuon",
                is_data=True,
                lumi_json = lumijson,
                lumi_csv = lumicsv_mu
            ),
        )
    if flavor == "e" or flavor == "ee":
        allProcs.append(
            narf.Dataset(
                filepaths=datasets2016.makeFilelist(["{BASE_PATH}/NanoAOD_v2/HighEGJet/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                name="singleelectron",
                is_data=True,
                lumi_json = lumijson,
                lumi_csv = lumicsv_el
            ),
        )
        
    if filt:
        allProcs = list(filter(filt, allProcs))

    for sample in allProcs:
        if not sample.filepaths:
            logger.warning(f"Failed to find any files for sample {sample.name}!")

    return allProcs
