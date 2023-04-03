from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets.datagroups import Datagroups
from wremnants.datasets import datasets2016
import narf
import math

logger = logging.child_logger(__name__)

class Datagroups2016(Datagroups):
    def __init__(self, infile, combine=False, pseudodata_pdfset = None, applySelection=True,
                 excludeGroups=None, filterGroups=None, splitWByCharge=False
    ):
        super().__init__(infile, combine, datasets2016.getDatasets())

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
                #label = r"Z$\to\mu\mu$ (N$^{3}LL+NNLO)$",
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
        
        self.filterGroups(filterGroups)
        self.excludeGroups(excludeGroups)

        if self.wmass:
            # add all processes to the fake contributions
            self.groups["Fake"] = dict(
                members = [member for sublist in [self.groups[x]["members"] for x in self.groups if x != "QCD"] for member in sublist],
                scale = lambda x: 1. if x.is_data else -1,
                label = "Nonprompt",
                color = "grey",
                selectOp = fakeOp,
            )
