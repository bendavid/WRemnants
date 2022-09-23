from wremnants import boostHistHelpers as hh
from wremnants import histselections as sel
from wremnants.datasets import datasetsLowPU
from wremnants.datasets.datagroups import datagroups
import logging
import lz4.frame
import pickle
import narf
import ROOT
import hist


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
                signalOp = None,
            ),
    
            # individual procs
            DYmumu=dict(
                members=[self.datasets["DYmumu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#F8CE68", #ROOT.TColor.GetColor(248, 206, 104),
                signalOp = None,
            ),
            DYee=dict(
                members=[self.datasets["DYee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color="#F8CE68",
                signalOp = None,
            ),
            
            WJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="WJets",
                color="#F8CE68",
                signalOp = None,
            ),
            
            
            WplusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            WminusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            WplusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            WminusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            WplusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            WminusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = None,
            ),
            
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
                color="#64C0E8",
                signalOp = None,
            ),
            
            WZTo3LNu=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu"]],
                label="WZTo3LNu",
                color="#64C0E8",
                signalOp = None,
            ),
            WWTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
                label="WWTo2L2Nu",
                color="#64C0E8",
                signalOp = None,
            ),
            ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                signalOp = None,
            ),
            TTTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu"]],
                label = "TTTo2L2Nu",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = None,
            ),
            TTToSemiLeptonic=dict(
                members = [self.datasets[x] for x in ["TTToSemiLeptonic"]],
                label = "TTToSemiLeptonic",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = None,
            ),
            
        )
        
        if flavor == "mumu": # Z->mumu
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    signalOp = None,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                signalOp = None,
                ),
            )
        if flavor == "ee": # Z->ee
            self.groups.update(
                SingleElectron=dict(
                    members=[self.datasets["singleelectron"]],
                    label="Data",
                    color="#000000",
                    signalOp = None,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYmumu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8",
                signalOp = None,
                ),
            )

        if flavor == "mu": # W->mu
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    signalOp = None,
                ),
                VV=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="VV",
                color="#64C0E8",
                signalOp = None,
                ),
                WJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (DY, VV)",
                color="#64C0E8",
                signalOp = None,
                ),
                DY=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu"]],
                label="DY",
                color="#64C0E8",
                signalOp = None,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK",
                color="#64C0E8",
                signalOp = None,
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
        
            # merged procs
            TTbar=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "TTbar",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = self.signalHistWmass,
            ),
    
            # individual procs
            DYmumu=dict(
                members=[self.datasets["DYmumu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color="#F8CE68", #ROOT.TColor.GetColor(248, 206, 104),
                signalOp = self.signalHistWmass,
            ),
            DYee=dict(
                members=[self.datasets["DYee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color="#F8CE68",
                signalOp = self.signalHistWmass,
            ),
            
            WJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="WJets",
                color="#F8CE68",
                signalOp = self.signalHistWmass,
            ),
            
            
            WplusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WminusJetsToMuNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToMuNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WplusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WminusJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToTauNu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WplusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WminusJetsToENu=dict(
                members = [self.datasets[x] for x in ["WminusJetsToENu"]],
                label="WJets",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            
            WZTo3LNu=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu"]],
                label="WZTo3LNu",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            WWTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
                label="WWTo2L2Nu",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
            ),
            TTTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu"]],
                label = "TTTo2L2Nu",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = self.signalHistWmass,
            ),
            TTToSemiLeptonic=dict(
                members = [self.datasets[x] for x in ["TTToSemiLeptonic"]],
                label = "TTToSemiLeptonic",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = self.signalHistWmass,
            ),
            
            

        )
        
      
        if flavor == "mu": # W->mu
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    signalOp = self.signalHistWmass,
                ),
                VV=dict(
                members = [self.datasets[x] for x in ["WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="VV",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
                ),
                WJetsToTauNu=dict(
                members = [self.datasets[x] for x in ["WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK (DY, VV)",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
                ),
                DY=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu"]],
                label="DY",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
                ),
                EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "DYee", "DYmumu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu"]],
                label="EWK",
                color="#64C0E8",
                signalOp = self.signalHistWmass,
                ),
                Fake=dict(
                    members = [self.datasets[x] for x in ["singlemuon", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu", "DYtautau", "DYee", "DYmumu", "TTTo2L2Nu", "TTToSemiLeptonic"]],
                    label = "Nonprompt",
                    scale = lambda x: 1. if x.is_data else -1,
                    color="#A9A9A9", #ROOT.TColor.GetColor(222, 90, 106),  --> sel
                    signalOp = self.fakeHistABCD,
                ),
            )        

    def signalHistWmass(self, h, charge=None):
        s = hist.tag.Slicer()
        # , "eta" : s[-2.4j:2.4j]
        sel = {"passIso" : 1, "mt":s[40j:10000j], "eta" : s[::hist.sum]}
        if charge in [-1, 1]:
            sel.update({"charge" : -1j if charge < 0 else 1j})
        return h[sel]

    def fakeHistABCD(self, h):
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
        