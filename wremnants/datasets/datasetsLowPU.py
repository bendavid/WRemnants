import narf
import os
import pathlib
import socket
from wremnants.datasets import dataset_tools
from utilities import logging, common

logger = logging.child_logger(__name__)

lumijson = f"{common.data_dir}/lowPU/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt"
lumicsv_mu = f"{common.data_dir}/lowPU/bylsoutput_HLT_HIMu17_Full.csv"
lumicsv_el = f"{common.data_dir}/lowPU/bylsoutput_HLT_HIEle20_Full.csv"

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
def getDatasets(maxFiles=-1, filt=None, excl=None, flavor="",base_path=None):

    if not base_path:
        base_path = dataset_tools.getDataPath(lowpu=True)

    logger.info(f"Loading lowPU samples from {base_path}.")

    BR_W_LEP = 3*0.1086 # PDG

    allProcs = [
    
        narf.Dataset(
            name="TTTo2L2Nu",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=87.31483776,
            is_data=False,
            group="Top",
        ),
        narf.Dataset(
            name="TTToSemiLeptonic",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=364.35,
            is_data=False,
            group="Top",
        ),
        # narf.Dataset(
        #     name="TTToHadronic",
        #     filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
        #     xsec=380.10,
        #     is_data=False,
        #     group="Top",
        # ),

        narf.Dataset(
            name="WWTo2L2Nu",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=118.7*BR_W_LEP*BR_W_LEP,
            is_data=False,
            group="Diboson",
        ),
        narf.Dataset(
            name="WZTo3LNu",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=4.912,
            is_data=False,
            group="Diboson",
        ),
        narf.Dataset(
            name="ZZ",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/ZZ_TuneCP5_13TeV-pythia8/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=16.523,
            is_data=False,
            group="Diboson",
        ),
        narf.Dataset(
            name="Ztautau",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToTauTau_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=2025.74,
            is_data=False,
            group="Ztautau",
        ),
        narf.Dataset(
            name="WminusJetsToTauNu",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToTauNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=8677.3,
            is_data=False,
            group="Wtaunu",
        ),
        narf.Dataset(
            name="WplusJetsToTauNu",
            filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToTauNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
            xsec=11572.19,
            is_data=False,
            group="Wtaunu",
        ),
    ]

    if flavor == "" or flavor == "mu" or flavor == "mumu":
        allProcs += [
            narf.Dataset(
                name="WminusJetsToMuNu",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=8677.3, # 8562.66
                is_data=False,
                group="Wmunu",
            ),
            narf.Dataset(
                name="WplusJetsToMuNu",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=11811.4, # 11572.19 
                is_data=False,
                group="Wmunu",
            ),
            narf.Dataset(
                name="Zmumu",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToMuMu_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=2025.74, # 1976.1
                is_data=False,
                group="Zmumu",
            ),
            narf.Dataset(
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/SingleMuon/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                name="singlemuon",
                is_data=True,
                lumi_json = lumijson,
                lumi_csv = lumicsv_mu,
                group="Data",
            ),
        ]
    if flavor == "" or flavor == "e" or flavor == "ee":
        # electron only datasets
        allProcs += [
            narf.Dataset(
                name="WminusJetsToENu",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WminusJetsToENu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=8677.3,
                is_data=False,
                group="Wenu",
            ),
            narf.Dataset(
                name="WplusJetsToENu",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/WplusJetsToENu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=11572.19, # 
                is_data=False,
                group="Wenu",
            ),
            narf.Dataset(
                name="Zee",
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v3/DYJetsToEE_M-50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                xsec=2025.74,
                is_data=False,
                group="Zee",
            ),
            narf.Dataset(
                filepaths=dataset_tools.makeFilelist(["{BASE_PATH}/NanoAOD_v2/HighEGJet/*/*/*/*.root"], maxFiles, format_args=dict(BASE_PATH=base_path)),
                name="singleelectron",
                is_data=True,
                lumi_json = lumijson,
                lumi_csv = lumicsv_el,
                group="Data",
            ),
        ]

    allProcs = dataset_tools.filterProcs(filt, allProcs)
    allProcs = dataset_tools.excludeProcs(excl, allProcs)

    for sample in allProcs:
        if not sample.filepaths:
            logger.warning(f"Failed to find any files for dataset {sample.name}. Skipping!")

    return allProcs
