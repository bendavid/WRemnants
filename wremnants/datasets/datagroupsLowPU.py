from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets import datasetsLowPU
from wremnants.datasets.datagroups import Datagroups
import hist

logger = logging.child_logger(__name__)

def make_datagroups_lowPU(input_file, combine=False, flavor="", excludeGroups=None, filterGroups=None):

    dg = Datagroups(input_file, combine, datasetsLowPU.getDatasets(flavor=flavor))

    dg.isW = True if flavor in ["mu", "e"] else False

    if dg.isW:
        sigOp = sel.signalHistWmass
    else:
        sigOp = None

    # individual procs
    dg.addGroup("TTTo2L2Nu",
        members=[dg.datasets["TTTo2L2Nu"]],
        label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
        color="#DE5A6A",
        selectOp = sigOp,
    )
    dg.addGroup("TTToSemiLeptonic",
        members=[dg.datasets["TTToSemiLeptonic"]],
        label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
        color="#DE5A6A",
        selectOp = sigOp,
    )
    dg.addGroup("TTToHadronic",
        members=[dg.datasets["TTToHadronic"]],
        label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
        color="#DE5A6A",
        selectOp = sigOp,
    )
        
        
    dg.addGroup("Zmumu",
        members=[dg.datasets["Zmumu"]],
        label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
        color="#F8CE68",
        selectOp = sigOp,
    )
    dg.addGroup("Zee",
        members=[dg.datasets["Zee"]],
        label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
        color="#F8CE68",
        selectOp = sigOp,
    )
    dg.addGroup("Ztautau",
        members = [dg.datasets[x] for x in ["Ztautau"]],
        label=r"DY #rightarrow #tau^{#plus}#tau^{#minus} (MiNNLO)",
        color="#64C0E8",
        selectOp = sigOp,
    )
        
    dg.addGroup("WplusJetsToMuNu",
        members = [dg.datasets[x] for x in ["WplusJetsToMuNu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WminusJetsToMuNu",
        members = [dg.datasets[x] for x in ["WminusJetsToMuNu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WplusJetsToTauNu",
        members = [dg.datasets[x] for x in ["WplusJetsToTauNu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WminusJetsToTauNu",
        members = [dg.datasets[x] for x in ["WminusJetsToTauNu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WplusJetsToENu",
        members = [dg.datasets[x] for x in ["WplusJetsToENu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WminusJetsToENu",
        members = [dg.datasets[x] for x in ["WminusJetsToENu"]],
        label="WJets",
        color="#64C0E8",
        selectOp = sigOp,
    )
        
        
        
    dg.addGroup("WZTo3LNu",
        members = [dg.datasets[x] for x in ["WZTo3LNu"]],
        label="WZTo3LNu",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("WWTo2L2Nu",
        members = [dg.datasets[x] for x in ["WWTo2L2Nu"]],
        label="WWTo2L2Nu",
        color="#64C0E8",
        selectOp = sigOp,
    )
    dg.addGroup("ZZ",
        members = [dg.datasets[x] for x in ["ZZ"]],
        label="ZZ",
        color="#64C0E8",
        selectOp = sigOp,
    )

        
        
        # grouped procs
    dg.addGroup("Top",
        members = [dg.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"]],
        label = "Top",
        color="#DE5A6A",
        selectOp = sigOp,
    )
        
    dg.addGroup("WJetsToMuNu",
        members = [dg.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
        label="W^{#pm} #rightarrow #mu^{#pm}#nu",
        color="#F8CE68",
        selectOp = sigOp,
    )
    dg.addGroup("WJetsToENu",
        members = [dg.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu"]],
        label="W^{#pm} #rightarrow e^{#pm}#nu",
        color="#F8CE68",
        selectOp = sigOp,
    )
    dg.addGroup("WJetsToTauNu",
        members = [dg.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
        label="W^{#pm} #rightarrow #tau^{#pm}#nu",
        color="#F8CE68",
        selectOp = sigOp,            
    )
    
    
    if flavor == "mumu": # Z->mumu

        dg.addGroup("EWK",
            members = [dg.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
            label="EWK (#tau^{#plus}#tau^{#minus}, VV)",
            color="#64C0E8",
            selectOp = sigOp,
            )
        dg.addGroup("EWK_noZZ",
            members = [dg.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
            label="EWK_noZZ (#tau#tau, W, VV)",
            color="#64C0E8",
            selectOp = sigOp,
            )
        dg.addGroup("Other",
            members = [x for x in dg.datasets.values() if not x.is_data and x.name not in ["Zmumu", "Ztautau"]],
            label="Other",
            color="#64C0E8",
            selectOp = sigOp,
        )
    if flavor == "ee": # Z->ee
            
        dg.addGroup("EWK",
            members = [dg.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
            label="EWK (Z #rightarrow #tau#tau, W, VV)",
            color="#64C0E8",
            selectOp = None,
            )
        dg.addGroup("EWK_noZZ",
            members = [dg.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
            label="EWK_noZZ (Z #rightarrow #tau#tau, W, VV)",
            color="#64C0E8",
            selectOp = None,
            )
        dg.addGroup("ZZ",
            members = [dg.datasets[x] for x in ["ZZ"]],
            label="ZZ",
            color="#64C0E8",
            selectOp = None,
        )
    

        
    if flavor == "mu": # W->mu
            
        dg.addGroup("VV",
            members = [dg.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
            label="VV",
            color="#64C0E8",
            selectOp = sigOp,
            )
        dg.addGroup("WJetsToTauNu",
            members = [dg.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
            label="EWK (DY, VV)",
            color="#64C0E8",
            selectOp = sigOp,
            )
        dg.addGroup("WJetsToMuNu",
            members = [dg.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
            label="W^{#pm} #rightarrow #mu^{#pm}#nu",
            color="#F8CE68",
            selectOp = sigOp,
            )
        dg.addGroup("DY",
            members = [dg.datasets[x] for x in ["Ztautau", "Zee", "Zmumu"]],
            label="DY",
            color="#64C0E8",
            selectOp = sigOp,
            )
        dg.addGroup("EWK",
            members = [dg.datasets[x] for x in ["Ztautau", "Zee", "Zmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
            label="EWK (#mu^{#plus}#mu^{#minus}, #tau^{#plus}#tau^{#minus}, VV, #tau^{#pm})",
            color="#64C0E8",
            selectOp = sigOp,
        )        

    dg.filterGroups(filterGroups)
    dg.excludeGroups(excludeGroups)

    if dg.wmass:
        dg.addGroup("Fake",
            members = [dg.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "Ztautau", "Zee", "Zmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
            label = "Nonprompt",
            scale = lambda x: 1. if x.is_data else -1,
            color="#A9A9A9",
            selectOp = sel.fakeHistABCD,
        )

    # data
    if flavor == "mu" or flavor == "mumu":  
        dg.addGroup("SingleMuon",
            members=[dg.datasets["singlemuon"]],
            label="Data",
            color="#000000",
            selectOp = sigOp,
        )
    if flavor == "e" or flavor == "ee":  
        dg.addGroup("SingleElectron",
            members=[dg.datasets["singleelectron"]],
            label="Data",
            color="#000000",
            selectOp = sigOp,
        )

    return dg
