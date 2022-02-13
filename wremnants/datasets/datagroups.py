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

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def setHists(self, histname, groups=None, selectSignal=True, forceNonzero=True):
        if not groups:
            groups = self.groups.values()
        for group in self.groups.values():
            # Safer to set this to null again to avoid mixing observables/systs
            group["hist"] = None
            if group not in groups:
                continue
            for member in group["members"]:
                output = self.results[member.name]["output"]
                if histname not in output:
                    logging.warning(f"Histogram {histname} not found for process {member.name}")
                    continue

                h = output[histname]
                if forceNonzero:
                    h = hh.clipNegativeVals(h)
                scale = self.processScaleFactor(member)
                if "scale" in group:
                    scale = scale*group["scale"](member)
                hscale = h*scale
                group["hist"] = hscale if not group["hist"] else hh.addHists(hscale, group["hist"])
            if selectSignal and group["hist"]:
                group["hist"] = group["signalOp"](group["hist"])

    def datagroupsForHist(self, histname, groups=None, selectSignal=True, forceNonzero=True):
        self.setHists(histname, groups, selectSignal, forceNonzero)
        return self.groups

    def resultsDict(self):
        return self.results

class datagroups2016(datagroups):
    def __init__(self, infile):
        self.datasets = {x.name : x for x in datasets2016.getDatasets()}
        super().__init__(infile)
        self.groups =  {
            "Data" : dict(
                members = [self.datasets["dataPostVFP"]],
                color = "black",
                label = "Data",
                hist = None,
                signalOp = sel.signalHistWmass,
            ),
            "Fake" : dict(
                members = list(self.datasets.values()),
                scale = lambda x: 1. if x.is_data else -1,
                label = "Nonprompt",
                color = "grey",
                hist = None,
                signalOp = sel.fakeHistABCD,
            ),
            "Zmumu" : dict(
                members = [self.datasets["ZmumuPostVFP"]],
                label = r"Z$\to\mu\mu$",
                color = "lightblue",
                hist = None,
                signalOp = sel.signalHistWmass,
            ),   
            "Wtau" : dict(
                members = [self.datasets["WminustaunuPostVFP"], self.datasets["WplustaunuPostVFP"]],
                label = r"W$^{\pm}\to\tau\nu$",
                color = "orange",
                hist = None,
                signalOp = sel.signalHistWmass,
            ),
            "Wmunu" : dict(
                members = [self.datasets["WminusmunuPostVFP"], self.datasets["WplusmunuPostVFP"]],
                label = r"W$^{\pm}\to\mu\nu$",
                color = "darkred",
                hist = None,
                signalOp = sel.signalHistWmass,
            ),
            "Ztt" : dict(
                members = [self.datasets["ZtautauPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
                hist = None,
                signalOp = sel.signalHistWmass,
            ), 
            "Top" : dict(
                members = [self.datasets["TTSemileptonicPostVFP"], self.datasets["TTLeptonicPostVFP"]],
                label = "Top",
                color = "green",
                hist = None,
                signalOp = sel.signalHistWmass,
            ), 
            "Diboson" : dict(
                members = [self.datasets["WWPostVFP"]],
                label = "Diboson",
                color = "pink",
                hist = None,
                signalOp = sel.signalHistWmass,
            ), 
        } 
