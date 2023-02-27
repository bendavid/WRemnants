from utilities import boostHistHelpers as hh
from wremnants import histselections as sel
from wremnants.datasets import datasetsLowPU
from wremnants.datasets.datagroups import datagroups
import logging
import lz4.frame
import pickle
import narf
import ROOT
import hist

class datagroupsLowPU(datagroups):
    isW = False
    def __init__(self, infile, combine=False, flavor=""):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets(flavor=flavor)}
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
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                EWK_noZZ=dict(
                members = [self.datasets[x] for x in ["Ztautau", "Zee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
                label="EWK_noZZ (Z #rightarrow #tau#tau, W, VV)",
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
                label="EWK (#tau, VV, DY)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                Fake=dict(
                    members = [self.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "Ztautau", "Zee", "Zmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
                    label = "Nonprompt",
                    scale = lambda x: 1. if x.is_data else -1,
                    color="#A9A9A9", #ROOT.TColor.GetColor(222, 90, 106),  --> sel
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
            
    def signalHistSel(self, h, charge=None):
        s = hist.tag.Slicer()
        if self.isW:
            sel = {"passIso" : True, "passMT": True}
            if charge in [-1, 1]:
                sel.update({"charge" : -1j if charge < 0 else 1j})
            return h[sel]
        else: return h
                
       
    def fakeHistABCD(self, h):
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
        



class datagroupsLowPU_Z(datagroups):
    def __init__(self, infile, combine=False, flavor=""):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets()}
        super().__init__(infile, combine)
        #self.lumi = 0.199269742
        self.hists = {} # container storing temporary histograms
        self.groups = dict(
        
            # merged procs
            TTbar=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "TTbar",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                selectOp = None,
            ),
    
            # individual procs
            DYmumu=dict(
                members=[self.datasets["DYmumu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#F8CE68", #ROOT.TColor.GetColor(248, 206, 104),
                selectOp = None,
            ),
            DYee=dict(
                members=[self.datasets["DYee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color="#F8CE68",
                selectOp = None,
            ),
            
            WJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="WJets",
                color="#F8CE68",
                selectOp = None,
            ),
            
            
            WplusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            WminusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            WplusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            WminusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            WplusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            WminusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                selectOp = None,
            ),
            
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
                color="#64C0E8",
                selectOp = None,
            ),
            
            WZTo3LNu=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu"]],
                label="WZTo3LNu",
                color="#64C0E8",
                selectOp = None,
            ),
            WWTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
                label="WWTo2L2Nu",
                color="#64C0E8",
                selectOp = None,
            ),
            ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                selectOp = None,
            ),
            TTTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu"]],
                label = "TTTo2L2Nu",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                selectOp = None,
            ),
            TTToSemiLeptonic=dict(
                members = [self.datasets[x] for x in ["TTToSemiLeptonic"]],
                label = "TTToSemiLeptonic",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                selectOp = None,
            ),
            
        )
        
        if flavor == "mumu": # Z->mumu
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    selectOp = None,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                selectOp = None,
                ),
                EWK_noZZ=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu"]],
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
        if flavor == "ee": # Z->ee
            self.groups.update(
                SingleElectron=dict(
                    members=[self.datasets["singleelectron"]],
                    label="Data",
                    color="#000000",
                    selectOp = None,
                ),
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
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    selectOp = None,
                ),
                VV=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="VV",
                color="#64C0E8",
                selectOp = None,
                ),
                WJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (DY, VV)",
                color="#64C0E8",
                selectOp = None,
                ),
                DY=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu"]],
                label="DY",
                color="#64C0E8",
                selectOp = None,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK",
                color="#64C0E8",
                selectOp = None,
                ),
            )        

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
    def readHist(self, baseName, proc, syst = "", scaleOp=None, forceNonzero=True):
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        print(baseName, proc.name, histname, syst)
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        #print(h)
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale
        
        
        
        
        
class datagroupsLowPU_W(datagroups):
    def __init__(self, infile, combine=False, flavor=""):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets()}
        super().__init__(infile, combine)
        #self.lumi = 0.199269742
        self.hists = {} # container storing temporary histograms
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
                color="#F8CE68", #ROOT.TColor.GetColor(248, 206, 104),
                selectOp = self.signalHistSel,
            ),
            Zee=dict(
                members=[self.datasets["Zee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color="#F8CE68",
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
            
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
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

            
            
            # merged procs
            Top=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"]],
                label = "Top",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                selectOp = self.signalHistSel,
            ),
            
            WJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="W^{#pm} #rightarrow #mu^{#pm}#nu",
                color="#F8CE68",
                selectOp = self.signalHistSel,
            ),

        )
        
      
        if flavor == "mu": # W->mu
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    selectOp = self.signalHistSel,
                ),
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
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu"]],
                label="DY",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (#tau, VV, DY)",
                color="#64C0E8",
                selectOp = self.signalHistSel,
                ),
                Fake=dict(
                    members = [self.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "DYtautau", "DYee", "DYmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
                    label = "Nonprompt",
                    scale = lambda x: 1. if x.is_data else -1,
                    color="#A9A9A9", #ROOT.TColor.GetColor(222, 90, 106),  --> sel
                    selectOp = self.fakeHistABCD,
                ),
            )        

    def signalHistSel(self, h, charge=None):
        
        s = hist.tag.Slicer()
        # , "eta" : s[-2.4j:2.4j]
        ##sel = {"passIso" : 1, "mt":s[40j:10000j], "eta" : s[::hist.sum]}
        sel = {"passIso" : True, "passMT": True}
        if charge in [-1, 1]:
            sel.update({"charge" : -1j if charge < 0 else 1j})
        return h[sel]

    def fakeHistABCD(self, h):
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

    def fakeHistABCD_old(self, h):
        s = hist.tag.Slicer()
        #div = , , cutoff=1)
        #print(div)
        sf = h[{"passIso" : True, "mt":s[0j:40j], "eta" : s[::hist.sum]}].sum().value / h[{"passIso" : False, "mt":s[0j:40j], "eta" : s[::hist.sum]}].sum().value
        #sf = div.sum().value
        #print(sf)
        #sys.exit()
        #return 
        return h[{"passIso" : False, "mt":s[40j:10000j], "eta" : s[::hist.sum]}]*sf
        #hh.multiplyHists(hh.divideHists(h[{"passIso" : True, "mt":s[0j:40j]}], h[{"passIso" : False, "mt":s[0j:40j]}], cutoff=1), h[{"passIso" : False, "mt":s[40j:10000j]}])
        #return hh.multiplyHists(hh.divideHists(h[{"passIso" : True, "mt":s[0j:40j]}], h[{"passIso" : False, "mt":s[0j:40j]}], cutoff=1), h[{"passIso" : False, "mt":s[40j:10000j]}])


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
    def readHist(self, baseName, proc, syst = "", scaleOp=None, forceNonzero=True):
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
        
