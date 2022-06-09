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
        self.lumi = 0.199269742
        self.hists = {} # container storing temporary histograms
        self.groups = dict(
        
            # merged procs
            TTbar=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "TTbar",
                color="#DE5A6A", #ROOT.TColor.GetColor(222, 90, 106),
                signalOp = None,
            ),
            EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (Z #rightarrow #tau#tau, W, VV)",
                color="#64C0E8", #ROOT.TColor.GetColor(100, 192, 232),
                signalOp = None,
            ),
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
            
            # individual procs
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
                color="#64C0E8",
                signalOp = None,
            ),
            WJets=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu"]],
                label="WJets",
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
            
        )
        
        # data
        if flavor == "mumu":
            self.groups.update(
                SingleMuon=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="#000000",
                    signalOp = None,
                ),
            )
        if flavor == "ee":
            self.groups.update(
                SingleElectron=dict(
                    members=[self.datasets["singleelectron"]],
                    label="Data",
                    color="#000000",
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
        #print(baseName, proc, histname, syst)
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
        