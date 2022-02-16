from wremnants import boostHistHelpers as hh
from wremnants import histselections as sel
from wremnants.datasets import datasets2016
import logging
import lz4.frame
import pickle

class datagroups(object):
    def __init__(self, infile):
        with lz4.frame.open(infile) as f:
            self.results = pickle.load(f)
        if self.datasets:
            self.data = [x for x in self.datasets.values() if x.is_data]
        self.lumi = sum([self.results[x.name]["lumi"] for x in self.data])
        self.groups = {}
        self.nominalName = "nominal"

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def setHists(self, histname, procsToRead=None, label=None, nominalIfMissing=True, selectSignal=True, forceNonzero=True):
        if label == None:
            label = "hist"
        if not procsToRead:
            procsToRead = self.groups.keys()

        for procName, group in self.groups.items():
            # Safer to set this to null again to avoid mixing observables/systs
            group[label] = None
            if procName not in procsToRead:
                continue

            for member in group["members"]:
                try:
                    h = self.readHist(histname, member, group["scale"] if "scale" in group else None, forceNonzero)
                except ValueError as e:
                    if nominalIfMissing:
                        h = self.readHist(self.nominalName, member, group["scale"] if "scale" in group else None, forceNonzero)
                        pass
                    else:
                        logging.warning(str(e))
                        continue
                group[label] = h if not group[label] else hh.addHists(h, group[label])
            if selectSignal and group[label]:
                group[label] = group["signalOp"](group[label])

    def readHist(self, histname, proc, scaleOp=None, forceNonzero=True):
        output = self.results[proc.name]["output"]
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc}")
        h = output[histname]
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale


    def datagroupsForHist(self, histname, procsToRead=None, label="", dataHist="", selectSignal=True, forceNonzero=True):
        self.setHists(histname, procsToRead, label, dataHist, selectSignal, forceNonzero)
        return self.groups

    def resultsDict(self):
        return self.results

    def processes(self):
        return self.groups.keys()

class datagroups2016(datagroups):
    def __init__(self, infile):
        self.datasets = {x.name : x for x in datasets2016.getDatasets()}
        super().__init__(infile)
        self.groups =  {
            "Data" : dict(
                members = [self.datasets["dataPostVFP"]],
                color = "black",
                label = "Data",
                signalOp = sel.signalHistWmass,
            ),
            "Fake" : dict(
                members = list(self.datasets.values()),
                scale = lambda x: 1. if x.is_data else -1,
                label = "Nonprompt",
                color = "grey",
                signalOp = sel.fakeHistABCD,
            ),
            "Zmumu" : dict(
                members = [self.datasets["ZmumuPostVFP"]],
                label = r"Z$\to\mu\mu$",
                color = "lightblue",
                signalOp = sel.signalHistWmass,
            ),   
            "Wtau" : dict(
                members = [self.datasets["WminustaunuPostVFP"], self.datasets["WplustaunuPostVFP"]],
                label = r"W$^{\pm}\to\tau\nu$",
                color = "orange",
                signalOp = sel.signalHistWmass,
            ),
            "Wmunu" : dict(
                members = [self.datasets["WminusmunuPostVFP"], self.datasets["WplusmunuPostVFP"]],
                label = r"W$^{\pm}\to\mu\nu$",
                color = "darkred",
                signalOp = sel.signalHistWmass,
            ),
            "Ztt" : dict(
                members = [self.datasets["ZtautauPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
                signalOp = sel.signalHistWmass,
            ), 
            "Top" : dict(
                members = [self.datasets["TTSemileptonicPostVFP"], self.datasets["TTLeptonicPostVFP"]],
                label = "Top",
                color = "green",
                signalOp = sel.signalHistWmass,
            ), 
            "Diboson" : dict(
                members = [self.datasets["WWPostVFP"]],
                label = "Diboson",
                color = "pink",
                signalOp = sel.signalHistWmass,
            ), 
        } 
