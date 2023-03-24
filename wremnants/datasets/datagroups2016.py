from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets.datagroups import Datagroups
from wremnants.datasets import datasets2016
import narf
import pandas as pd
import math


logger = logging.child_logger(__name__)

class Datagroups2016(Datagroups):
    def __init__(self, infile, combine=False, pseudodata_pdfset = None, applySelection=True,
                 excludeProcGroup=None, filterProcGroup=None
    ):
        super().__init__(infile, combine, datasets2016.getDatasets(filt=filterProcGroup, excl=excludeProcGroup))

        if self.wmass and applySelection:
            sigOp = sel.signalHistWmass
            fakeOp = sel.fakeHistABCD
        else:
            sigOp = None
            fakeOp = None

        self.groups =  {
            "Data" : dict(
                members = self.getSafeListFromDataset(["dataPostVFP"]),
                color = "black",
                label = "Data",
                selectOp = sigOp,
            ),
            "Zmumu" : dict(
                members = self.getSafeListFromDataset(["ZmumuPostVFP"]),
                label = r"Z$\to\mu\mu$",
                color = "lightblue",
                selectOp = sigOp,
            ),   
            "Ztautau" : dict(
                members = self.getSafeListFromDataset(["ZtautauPostVFP"]),
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
                selectOp = sigOp,
            ),            
        }
        if pseudodata_pdfset and combine:
            self.groups[f"pdf{pseudodata_pdfset.upper()}_sum"] = dict(
                label = f"pdf{pseudodata_pdfset.upper()}",
                color = "dimgray"
            )
        if self.wmass:
            self.groups.update({
                "Wmunu" : dict(
                    members = self.getSafeListFromDataset(["WminusmunuPostVFP", "WplusmunuPostVFP"]),
                    label = r"W$^{\pm}\to\mu\nu$",
                    color = "darkred",
                    selectOp = sigOp,
                ),
                }
            )
            # Reorder
            for k in ["Zmumu", "Ztautau"]:
                self.groups[k] = self.groups.pop(k)
            self.groups.update({
                "Wtaunu" : dict(
                    members = self.getSafeListFromDataset(["WminustaunuPostVFP", "WplustaunuPostVFP"]),
                    label = r"W$^{\pm}\to\tau\nu$",
                    color = "orange",
                    selectOp = sigOp,
                ),
                "Top" : dict(
                    members = list(filter(lambda y: y.group == "Top", self.datasets.values())),
                    label = "Top",
                    color = "green",
                    selectOp = sigOp,
                ), 
                "Diboson" : dict(
                    members = list(filter(lambda y: y.group == "Diboson", self.datasets.values())),
                    label = "Diboson",
                    color = "pink",
                    selectOp = sigOp,
                ), 
                "Fake" : dict(
                    members = list(filter(lambda y: y.group != "QCD", self.datasets.values())),
                    scale = lambda x: 1. if x.is_data else -1,
                    label = "Nonprompt",
                    color = "grey",
                    selectOp = fakeOp,
                ),
                "QCD" : dict(
                    members = list(filter(lambda y: y.group == "QCD", self.datasets.values())),
                    label = "QCD MC",
                    color = "grey",
                    selectOp = sigOp,
                ), 
           
            })
        else:
            self.groups["Other"] = dict(
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["ZmumuPostVFP", "ZtautauPostVFP"] and x.group != "QCD"],
                label = "Other",
                color = "grey",
            )
            
        # keep track of groups which have at least one process member after the filters
        # check also that the key is not in excludeProcGroup (e.g. Fake must be excluded here,
        # since it is created after filtering the input datasetDict)
        # if excludeProcGroup is None:
        #     excludeProcGroup = []

        # self.updateGroupNamesPostFilter(excludeGroup=excludeProcGroup)
        # #self.groupNamesPostFilter = list(x for x in self.groups.keys() if len(self.groups[x]["members"]) and x not in excludeProcGroup)
        # logger.debug(f"Filtered groups: {self.groupNamesPostFilter}")
        
    def make_yields_df(self, histName, procs, action):
        def sum_and_unc(h):
            return (h.sum().value, math.sqrt(h.sum().variance))
        df = pd.DataFrame([(k, *sum_and_unc(action(v[histName]))) for k,v in self.groups.items() if k in procs], 
                columns=["Process", "Yield", "Uncertainty"])
        return df

    def readHist(self, baseName, proc, syst, scaleOp=None, forceNonzero=True, scaleToNewLumi=-1):
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        logger.debug(f"Reading hist {histname} for proc {proc.name} and syst {syst}")
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        if isinstance(h, narf.ioutils.H5PickleProxy):
            h = h.get()
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        if scaleToNewLumi > 0:
            h = hh.scaleByLumi(h, scaleToNewLumi, createNew=True)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale
