from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets import datasetsLowPU
from wremnants.datasets.datagroups import Datagroups
import hist

logger = logging.child_logger(__name__)

def fakeHistABCD(h):
    if "mt" in [ax.name for ax in h.axes]:
        return h[{"passIso" : False, "passMT" : True}]*h[{"passIso" : True, "passMT" : False}].sum().value / h[{"passIso" : False, "passMT" : False}].sum().value

    return sel.fakeHistABCD(h) 

def signalHistSel(h, charge=None, genBin=None):
    if genBin != None:
        h = h[{"recoil_gen" : genBin}]

    sel = {"passIso" : True, "passMT": True}
    if charge in [-1, 1]:
        sel.update({"charge" : -1j if charge < 0 else 1j})
    for key in sel.copy().keys():
        if not key in [ax.name for ax in h.axes]: # remove ax slice if the ax does not exist
            del sel[key]
    return h[sel]

class DatagroupsLowPU(Datagroups):
    def __init__(self, infile, combine=False, flavor="", filterGroups=None, excludeGroups=None):
        super().__init__(infile, combine, datasetsLowPU.getDatasets(flavor=flavor))

        self.isW = True if flavor in ["mu", "e"] else False

        if self.isW:
            sigOp = signalHistSel
        else:
            sigOp = None

        # individual procs
        self.addGroup("TTTo2L2Nu",
            members=[self.datasets["TTTo2L2Nu"]],
            label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
            color="#DE5A6A",
            selectOp = sigOp,
        )
        self.addGroup("TTToSemiLeptonic",
            members=[self.datasets["TTToSemiLeptonic"]],
            label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
            color="#DE5A6A",
            selectOp = sigOp,
        )
        self.addGroup("TTToHadronic",
            members=[self.datasets["TTToHadronic"]],
            label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
            color="#DE5A6A",
            selectOp = sigOp,
        )
            
            
        self.addGroup("Zmumu",
            members=[self.datasets["Zmumu"]],
            label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
            color="#F8CE68",
            selectOp = sigOp,
        )
        self.addGroup("Zee",
            members=[self.datasets["Zee"]],
            label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
            color="#F8CE68",
            selectOp = sigOp,
        )
        self.addGroup("Ztautau",
            members = [self.datasets[x] for x in ["Ztautau"]],
            label=r"DY #rightarrow #tau^{#plus}#tau^{#minus} (MiNNLO)",
            color="#64C0E8",
            selectOp = sigOp,
        )
            
        self.addGroup("WplusJetsToMuNu",
            members = [self.datasets[x] for x in ["WplusJetsToMuNu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WminusJetsToMuNu",
            members = [self.datasets[x] for x in ["WminusJetsToMuNu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WplusJetsToTauNu",
            members = [self.datasets[x] for x in ["WplusJetsToTauNu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WminusJetsToTauNu",
            members = [self.datasets[x] for x in ["WminusJetsToTauNu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WplusJetsToENu",
            members = [self.datasets[x] for x in ["WplusJetsToENu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WminusJetsToENu",
            members = [self.datasets[x] for x in ["WminusJetsToENu"]],
            label="WJets",
            color="#64C0E8",
            selectOp = sigOp,
        )
            
            
            
        self.addGroup("WZTo3LNu",
            members = [self.datasets[x] for x in ["WZTo3LNu"]],
            label="WZTo3LNu",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("WWTo2L2Nu",
            members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
            label="WWTo2L2Nu",
            color="#64C0E8",
            selectOp = sigOp,
        )
        self.addGroup("ZZ",
            members = [self.datasets[x] for x in ["ZZ"]],
            label="ZZ",
            color="#64C0E8",
            selectOp = sigOp,
        )

            
            
            # grouped procs
        self.addGroup("Top",
            members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"]],
            label = "Top",
            color="#DE5A6A",
            selectOp = sigOp,
        )
            
        self.addGroup("WJetsToMuNu",
            members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
            label="W^{#pm} #rightarrow #mu^{#pm}#nu",
            color="#F8CE68",
            selectOp = sigOp,
        )
        self.addGroup("WJetsToENu",
            members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu"]],
            label="W^{#pm} #rightarrow e^{#pm}#nu",
            color="#F8CE68",
            selectOp = sigOp,
        )
        self.addGroup("WJetsToTauNu",
            members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
            label="W^{#pm} #rightarrow #tau^{#pm}#nu",
            color="#F8CE68",
            selectOp = sigOp,            
        )
        
        
        if flavor == "mumu": # Z->mumu

            self.addGroup("EWK",
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (#tau^{#plus}#tau^{#minus}, VV)",
                color="#64C0E8",
                selectOp = sigOp,
                )
            self.addGroup("EWK_noZZ",
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
                label="EWK_noZZ (#tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = sigOp,
                )
            self.addGroup("Other",
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["Zmumu", "Ztautau"]],
                label="Other",
                color="#64C0E8",
                selectOp = sigOp,
            )
        if flavor == "ee": # Z->ee
                
            self.addGroup("EWK",
                members = [self.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = None,
                )
            self.addGroup("EWK_noZZ",
                members = [self.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
                label="EWK_noZZ (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = None,
                )
            self.addGroup("ZZ",
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                selectOp = None,
            )
        

            
        if flavor == "mu": # W->mu
                
            self.addGroup("VV",
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="VV",
                color="#64C0E8",
                selectOp = sigOp,
                )
            self.addGroup("WJetsToTauNu",
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (DY, VV)",
                color="#64C0E8",
                selectOp = sigOp,
                )
            self.addGroup("WJetsToMuNu",
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="W^{#pm} #rightarrow #mu^{#pm}#nu",
                color="#F8CE68",
                selectOp = sigOp,
                )
            self.addGroup("DY",
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "Zmumu"]],
                label="DY",
                color="#64C0E8",
                selectOp = sigOp,
                )
            self.addGroup("EWK",
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "Zmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (#mu^{#plus}#mu^{#minus}, #tau^{#plus}#tau^{#minus}, VV, #tau^{#pm})",
                color="#64C0E8",
                selectOp = sigOp,
            )        

        self.filterGroups(filterGroups)
        self.excludeGroups(excludeGroups)

        if self.wmass:
            self.addGroup("Fake",
                members = [self.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "Ztautau", "Zee", "Zmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "Nonprompt",
                scale = lambda x: 1. if x.is_data else -1,
                color="#A9A9A9",
                selectOp = fakeHistABCD,
            )

        # data
        if flavor == "mu" or flavor == "mumu":  
            self.addGroup("SingleMuon",
                members=[self.datasets["singlemuon"]],
                label="Data",
                color="#000000",
                selectOp = sigOp,
            )
        if flavor == "e" or flavor == "ee":  
            self.addGroup("SingleElectron",
                members=[self.datasets["singleelectron"]],
                label="Data",
                color="#000000",
                selectOp = sigOp,
            )
    
        