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

def signalOp(h):
    #print(h)
    return h

class datagroupsLowPU_Z(datagroups):
    def __init__(self, infile, combine=False, flavor=""):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets()}
        super().__init__(infile, combine)
        self.lumi = 0.199269742
        self.hists = {} # container storing temporary histograms
        self.groups = dict(
            TTbar=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "TTbar",
                color=ROOT.TColor.GetColor(222, 90, 106),
                signalOp = signalOp,
            ),
            EWK=dict(
                members = [self.datasets[x] for x in ["DYtautau", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu", "WZTo3LNu", "WWTo2L2Nu", "ZZ"]],
                label="EWK (W #rightarrow #tau, diboson)",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            DYmumu=dict(
                members=[self.datasets["DYmumu"]],
                label=r"DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO)",
                color=ROOT.TColor.GetColor(248, 206, 104),
                signalOp = signalOp,
            ),
            DYee=dict(
                members=[self.datasets["DYee"]],
                label=r"DY #rightarrow e^{#plus}e^{#minus} (MiNNLO)",
                color=ROOT.TColor.GetColor(248, 206, 104),
                signalOp = signalOp,
            ),
            SingleMuon=dict(
                members=[self.datasets["singlemuon"]],
                label="Data",
                color=ROOT.kBlack,
                signalOp = signalOp if flavor == "mumu" else None,
            ),
            SingleElectron=dict(
                members=[self.datasets["singleelectron"]],
                label="Data",
                color=ROOT.kBlack,
                signalOp = signalOp if flavor == "ee" else None,
            ),
            
            
            ## individual procs
            
            DYtautau=dict(
                members = [self.datasets[x] for x in ["DYtautau"]],
                label="DYtautau",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            WJets=dict(
                members = [self.datasets[x] for x in ["WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu", "WplusJetsToENu", "WminusJetsToENu"]],
                label="WJets",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            WZTo3LNu=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu"]],
                label="WZTo3LNu",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            WWTo2L2Nu=dict(
                members = [self.datasets[x] for x in ["WWTo2L2Nu"]],
                label="WWTo2L2Nu",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            ZZ=dict(
                members = [self.datasets[x] for x in ["ZZ"]],
                label="ZZ",
                color=ROOT.TColor.GetColor(100, 192, 232),
                signalOp = signalOp,
            ),
            
        )
        

    
        

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]
        
        
    def histName(self, baseName, procName, syst):
        
        if baseName == "reco_mll" and (procName == "DYmumu" or procName == "DYee"): 
            baseName = "gen_reco_mll"
        #print(baseName, procName)
        
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
        
    # read single proc (= sum of procs)
    def readProc(self, baseName, procName, syst = "", scaleOp=None, forceNonzero=True, axis=None):

        label = "%s_%s" % (baseName, procName)
        if syst != "": label += "_%s" % syst
        
        
        if not label in self.hists:
        
            group = self.groups[procName]
            hist = None
            for member in group["members"]:
                try:
                    h = self.readHist(baseName, member, syst, scaleOp=None, forceNonzero=forceNonzero)
                except ValueError as e:
                    if nominalIfMissing:
                        h = self.readHist(baseName, member, self.nominalName, scaleOp=None, forceNonzero=forceNonzero)
                    else:
                        logging.warning(str(e))
                        continue
                hist = h if not hist else hh.addHists(h, hist)
                
            self.hists[label] = hist

        if axis != None: hist = self.hists[label].project(axis)
        return hist

    def histName__(self, baseName, proc, syst):
        if proc in ["WplusJetsToMuNu", "WminusJetsToMuNu"] and "gen" not in baseName:
            baseName = baseName.replace("reco", "gen_reco")
        base = f"{baseName}_{proc}"
        return base if syst == "nominal" else f"{base}_{syst}_syst"

    def readHist__(self, baseName, proc, syst, scaleOp=None, forceNonzero=True):
        axisNames = None
        readname = self.histName(baseName, proc.name, syst)
        if "mt_reco_pf" in readname:
            axisNames = ["qTreco", "iso", "charge", "mt"]
        elif "mt_gen_reco_pf" in readname:
            axisNames = ["qTreco", "qTgen", "iso", "charge", "mt"]
        if syst != "nominal":
            axisNames.append("systAx")

        rthist = self.rtfile.Get(readname)
        if not rthist:
            raise ValueError(f"Histogram {readname} not found for process {proc.name}")
        h = narf.root_to_hist(rthist, axis_names=axisNames)
        # TODO: this is a hack. For this to be smoothly treated, qTgen should be the first axis
        if "qTgen" in axisNames:
            s = hist.tag.Slicer()
            h = h[{"qTgen" : s[::hist.sum]}]

        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale




class datagroupsLowPU_Z_old(datagroups):
    def __init__(self, infile, flavor="muon", combine=False):
        self.datasets = {x.name : x for x in datasetsLowPU.getDatasets_Z()}
        super().__init__(infile, combine)
        self.lumi = 0.199 # TODO More digits
        
        self.groups = dict(
            Top=dict(
                members = [self.datasets[x] for x in ["TTTo2L2Nu", "TTToSemiLeptonic"]],
                label = "Top",
                color="green",
                signalOp = sel.signalHistLowPileupZ,
            ),
            EWK=dict(
                members = [self.datasets[x] for x in ["WZTo3LNu", "WWTo2L2Nu", "WminusJetsToTauNu", "WplusJetsToTauNu", "WplusJetsToMuNu", "WminusJetsToMuNu"]],
                label="Diboson",
                color="pink",
                signalOp = sel.signalHistLowPileupZ,
            ),
            Zmumu=dict(
                members=[self.datasets["DYmumu_MiNNLO"]],
                label=r"Z$\to\mu\mu$",
                color="lightblue",
                signalOp = sel.signalHistLowPileupZ,
            ),
            
            
            Data=dict(
                members=[self.datasets["singlemuon"]],
                label="Data",
                color="black",
                signalOp = sel.signalHistLowPileupZ,
            ),
        )        
        if flavor == "muon":
        
            self.groups.update(
                Data=dict(
                    members=[self.datasets["singlemuon"]],
                    label="Data",
                    color="black",
                    signalOp = sel.signalHistLowPileupZ,
                )
            )
            
        if flavor == "electron":
        
            self.groups.update(
                Data=dict(
                    members=[self.datasets["singleelectron"]],
                    label="Data",
                    color="black",
                    signalOp = sel.signalHistLowPileupZ,
                )
            )
        
        

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi
        
    def histName(self, baseName, proc, syst):
        if proc in ["DYmumu_MiNNLO"] and "gen" not in baseName:
            baseName = baseName.replace("reco", "gen_reco") # m_reco --> m_gen_reco
        base = f"{baseName}_{proc}"
        
        if "recoilStatUnc" in syst: return f"{baseName}_{syst}_syst_{proc}"
        
        return base if syst == "nominal" else f"{base}_{syst}_syst"
        

    '''
    m_test                              : qTBinGen::qTBinReco_nom::m_ll::RecoilSystStatIdx              : 20,0,20::20,0,20::60,60,120::100,0,100   : proc=DYmumu_MiNNLO
m_recoilunc_stat_nom                : qTBinGen::qTBinReco_nom::m_ll::RecoilSystStatIdx              : 20,0,20::20,0,20::60,60,120::100,0,100   : proc=DYmumu_MiNNLO
m_recoilunc_stat_nom_para_data_m1   : qTBinGen::qTBinReco_para_data_m1::m_ll::RecoilSystStatIdx     : 20,0,20::20,0,20::60,60,120::100,0,100   : proc=DYmumu_MiNNLO

    '''
        
    def readHist(self, baseName, proc, syst, scaleOp=None, forceNonzero=True):
    
        axisNames = None
        readname = self.histName(baseName, proc.name, syst)
        
        if "m_reco" in readname:
            axisNames = ["qTreco", "m"]
        elif "m_gen_reco" in readname:
            axisNames = ["qTgen", "qTreco", "m"]
        elif "m_ll" in readname: # for inclusive analysis
            axisNames = ["m"]
        
        if "recoilStatUnc" in syst:
            axisNames = ["qTgen", "qTreco", "qTrecoPert", "m"] # , "RecoilSystStatIdx"

        if syst != "nominal":
            axisNames.append("systAx")

        rthist = self.rtfile.Get(readname)
        if not rthist:
            raise ValueError(f"Histogram {readname} not found for process {proc.name}")
        #print("Start conversion", rthist.GetNdimensions(), axisNames, baseName, proc.name, syst, "->", readname)
        h = narf.root_to_hist(rthist, axis_names=axisNames)
        #print("End conversion")
        # TODO: this is a hack. For this to be smoothly treated, qTgen should be the first axis
        #if "qTgen" in axisNames:
        #    s = hist.tag.Slicer()
        #    h = h[{"qTgen" : s[::hist.sum]}]

        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale

