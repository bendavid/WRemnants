import narf
import os
import pathlib

lumijson = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt"
lumicsv_mu = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/bylsoutput_HLT_HIMu17_Full.csv"
lumicsv_el = f"{pathlib.Path(__file__).parent.parent}/data/lowPU/bylsoutput_HLT_HIEle20_Full.csv"

def findEOS(basedir, regex = ""):
    
    if ".root" in basedir: return basedir

    if regex != "":
    
        if basedir[-1] == "/": basedir = basedir[:-1]
        regex = basedir + "/" + regex

    rootFiles = []
    for root, directories, filenames in os.walk(basedir):
    
        for f in filenames:
       
            filePath = os.path.join(os.path.abspath(root), f)
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            if regex == "" or fnmatch.fnmatch(filePath, regex): rootFiles.append(filePath)
       
    return rootFiles  

# TODO: Allow filtering
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
def getDatasets(maxFiles=-1, filt=None):

    BR_W_LEP = 3*0.1086 # PDG

    allProcs = [
    
        narf.Dataset(
            name="TTTo2L2Nu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"),
            xsec=87.31483776,
            is_data=False,
        ),
        narf.Dataset(
            name="TTToSemiLeptonic",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"),
            xsec=364.35,
            is_data=False,
        ),
        narf.Dataset(
            name="TTToHadronic",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/TTToHadronic_TuneCP5_13TeV-powheg-pythia8"),
            xsec=380.10,
            is_data=False,
        ),
        
       
            

        
        narf.Dataset(
            name="WWTo2L2Nu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"),
            xsec=118.7*BR_W_LEP*BR_W_LEP,
            is_data=False,
        ),
        narf.Dataset(
            name="WZTo3LNu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8"),
            xsec=4.912,
            is_data=False,
        ),
        narf.Dataset(
            name="ZZ",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/ZZ_TuneCP5_13TeV-pythia8"),
            xsec=16.523,
            is_data=False,
        ),


        narf.Dataset(
            name="DYmumu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/DYJetsToMuMu_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/"),
            xsec=2025.74, # 1976.1
            is_data=False,
        ),
        narf.Dataset(
            name="DYee",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=2025.74,
            is_data=False,
        ),
        
        narf.Dataset(
            name="DYtautau",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/DYJetsToTauTau_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=2025.74,
            is_data=False,
        ),
        
        
        
        narf.Dataset(
            name="WminusJetsToMuNu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WminusJetsToMuNu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=8677.3, # 8562.66
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToENu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WminusJetsToENu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=8677.3,
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToTauNu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WminusJetsToTauNu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=8677.3,
            is_data=False,
        ),
        
        
        narf.Dataset(
            name="WplusJetsToMuNu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WplusJetsToMuNu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=11811.4, # 11572.19 
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToENu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WplusJetsToENu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=11572.19, # 
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToTauNu",
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/WplusJetsToTauNu_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"),
            xsec=11572.19,
            is_data=False,
        ),
        
        
        narf.Dataset(
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/SingleMuon/"),
            name="singlemuon",
            is_data=True,
            lumi_json = lumijson,
            lumi_csv = lumicsv_mu
        ),
        narf.Dataset(
            filepaths=findEOS("/scratch/shared/lowPU/NanoAOD_v2/HighEGJet/"),
            name="singleelectron",
            is_data=True,
            lumi_json = lumijson,
            lumi_csv = lumicsv_el
        ),
    ]

    if filt:
        return list(filter(filt, allProcs))

    return allProcs



def getDatasets_Z(maxFiles=-1, filt=None):
    allProcs = [
        narf.Dataset(
            name="TTTo2L2Nu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="TTToSemiLeptonic",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WWTo2L2Nu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WZTo3LNu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="DYmumu_MiNNLO",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToTauNu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToTauNu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WminusJetsToMuNu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            name="WplusJetsToMuNu",
            filepaths=[],
            xsec=1.,
            is_data=False,
        ),
        narf.Dataset(
            filepaths=[],
            name="singlemuon",
            xsec=1.,
            is_data=True,
        ),
    ]

    if filt:
        return list(filter(filt, allProcs))

    return allProcs
