from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets import datasetsLowPU
from wremnants.datasets.datagroups import datagroups
import lz4.frame
import pickle
import narf
import ROOT
import hist

logger = logging.child_logger(__name__)

class datagroupsLowPU(datagroups):
    isW = False
    def __init__(self, infile, combine=False, flavor="", excludeProcGroup=None, filterProcGroup=None):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets(flavor=flavor, filt=filterProcGroup, excl=excludeProcGroup)}
        super().__init__(infile, combine)
        #self.lumi = 0.199269742
        self.hists = {} # container storing temporary histograms
        self.isW = True if flavor in ["mu", "e"] else False
        self.groups = dict(
        
            # individual procs
            TTTo2L2Nu=dict(
                members=[self.datasets["TTTo2L2Nu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#DE5A6A",
                selectOp = self.signalHistSel,
            ),
            TTToSemiLeptonic=dict(
                members=[self.datasets["TTToSemiLeptonic"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#DE5A6A",
                selectOp = self.signalHistSel,
            ),
            TTToHadronic=dict(
                members=[self.datasets["TTToHadronic"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#DE5A6A",
                selectOp = self.signalHistSel,
            ),
            
            
            Zmumu=dict(
                members=[self.datasets["Zmumu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),
            Zee=dict(
                members=[self.datasets["Zee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),
            Ztautau=dict(
                members = [self.datasets[x] for x in ["Ztautau"]],
                label=r"DY #rightarrow #tau^{#plus}#tau^{#minus} (MiNNLO)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            
            WplusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WminusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WplusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WminusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WplusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WminusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            
            
            
            WZTo3LNu=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu"]],
                label="WZTo3LNu",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            WWTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
                label="WWTo2L2Nu",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),
            ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                selectOp = self.signalHistSel,
            ),

            
            
            # grouped procs
            Top=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"]],
                label = "Top",
                color="#DE5A6A",
                selectOp = self.signalHistSel,
            ),
            
            WJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="W^{#pm} #rightarrow #mu^{#pm}#nu",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),
            WJetsToENu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu"]],
                label="W^{#pm} #rightarrow e^{#pm}#nu",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),
            WJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="W^{#pm} #rightarrow #tau^{#pm}#nu",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),
            
        )
        
        
        if flavor == "mumu": # Z->mumu
            self.groups.update(

                EWK=dict(
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (#tau^{#plus}#tau^{#minus}, VV)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                EWK_noZZ=dict(
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
                label="EWK_noZZ (#tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                Other=dict(
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["Zmumu", "Ztautau"]],
                label="Other",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),

            )
        if flavor == "ee": # Z->ee
            self.groups.update(
                
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = None,
                ),
                EWK_noZZ=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
                label="EWK_noZZ (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = None,
                ),
                ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                selectOp = None,
                ),
            )
        

            
        if flavor == "mu": # W->mu
            self.groups.update(
                
                VV=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="VV",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                WJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (DY, VV)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                DY=dict(
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "Zmumu"]],
                label="DY",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "Zmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (#mu^{#plus}#mu^{#minus}, #tau^{#plus}#tau^{#minus}, VV, #tau^{#pm})",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                Fake=dict(
                    members = [self.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "Ztautau", "Zee", "Zmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
                    label = "Nonprompt",
                    scale = lambda x: 1. if x.is_data else -1,
                    color="#A9A9A9",
                    selectOp = self.fakeHistABCD,
                ),
            )        


        # data
        if flavor == "mu" or flavor == "mumu":  
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    selectOp = self.signalHistSel,
                ),
            )
        if flavor == "e" or flavor == "ee":  
            self.groups.update(
                SingleElectron=dict(
                    members=[self.datasets["singleelectron"]],
                    label="Data",
                    color="#000000",
                    selectOp = self.signalHistSel,
                ),
            )

        # this class doesn't implement exclusion of groups yet, and maybe it is not needed.
        # If if it is ever added, next line needs to be modified accordingly
        self.updateGroupNamesPostFilter(excludeGroup=[])
        #self.groupNamesPostFilter = list(x for x in self.groups.keys() if len(self.groups[x]["members"])) # and x not in excludeProcGroup)
        logger.debug(f"Filtered groups: {self.groupNamesPostFilter}")
            
    def signalHistSel(self, h, charge=None):
        s = hist.tag.Slicer()
        axes = [ax.name for ax in h.axes]
        if self.isW:
            sel = {"passIso" : True, "passMT": True}
            if charge in [-1, 1]:
                sel.update({"charge" : -1j if charge < 0 else 1j})
            for key in sel.copy().keys():
                if not key in axes:
                    del sel[key]
            return h[sel]
        else: return h
                
       
    def fakeHistABCD(self, h):
        axes = [ax.name for ax in h.axes]
        if "mt" in axes:
            s = hist.tag.Slicer()
            sf = h[{"passIso" : True, "passMT" : False}].sum().value / h[{"passIso" : False, "passMT" : False}].sum().value
            return h[{"passIso" : False, "passMT" : True}]*sf
        
        ret = hh.multiplyHists(
            hh.divideHists(h[{"passIso" : True, "passMT" : False}], 
                h[{"passIso" : False, "passMT" : False}],
                    cutoff=1
                ),
                    #where=h[{"passIso" : False, "passMT" : True}].values(flow=True)>1),
            h[{"passIso" : False, "passMT" : True}], 
        )
        return ret


    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]
        
        
    def histName(self, baseName, procName, syst):
        
        if baseName == "reco_mll" and (procName == "DYmumu" or procName == "DYee"): 
            baseName = "gen_reco_mll"

        # This is kind of hacky to deal with the different naming from combine
        if baseName != "x" and (syst == "" or syst == self.nominalName):
            return baseName
        if (baseName == "" or baseName == "x") and syst:
            return syst
        return "_".join([baseName,syst])
    
    # read single histogram (name, proc and syst)
    def readHist(self, baseName, proc, syst = "", scaleOp=None, forceNonzero=True, scaleToNewLumi=-1): 
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        #print(baseName, proc.name, histname, syst)
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        if isinstance(h, narf.ioutils.H5PickleProxy):
            h = h.get()
        #print(h)
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale
        

