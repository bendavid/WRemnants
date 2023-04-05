from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
from wremnants.datasets.datagroups import Datagroups
from wremnants.datasets.datagroup import Datagroup
from wremnants.datasets import datasets2016

logger = logging.child_logger(__name__)

class Datagroups2016(Datagroups):
    def __init__(self, infile, combine=False, pseudodata_pdfset = None, applySelection=True,
                 excludeGroups=None, filterGroups=None
    ):
        super().__init__(infile, combine, datasets2016.getDatasets())

        if self.wmass and applySelection:
            sigOp = sel.signalHistWmass
            fakeOp = sel.fakeHistABCD
        else:
            sigOp = None
            fakeOp = None

        self.addGroup("Data",
            members = self.getSafeListFromDataset(["dataPostVFP"]),
            color = "black",
            label = "Data",
            selectOp = sigOp,
        )
        self.addGroup("Zmumu",
            members = self.getSafeListFromDataset(["ZmumuPostVFP"]),
            #label = r"Z$\to\mu\mu$ (N$^{3}LL+NNLO)$",
            label = r"Z$\to\mu\mu$",
            color = "lightblue",
            selectOp = sigOp,
        ) 
        self.addGroup("Ztautau",
            members = self.getSafeListFromDataset(["ZtautauPostVFP"]),
            label = r"Z$\to\tau\tau$",
            color = "darkblue",
            selectOp = sigOp,
        )

        if pseudodata_pdfset and combine:
            self.addGroup(f"pdf{pseudodata_pdfset.upper()}_sum",
                label = f"pdf{pseudodata_pdfset.upper()}",
                color = "dimgray"
            )
        if self.wmass:
            self.addGroup("Wmunu",
                members = self.getSafeListFromDataset(["WminusmunuPostVFP", "WplusmunuPostVFP"]),
                label = r"W$^{\pm}\to\mu\nu$",
                color = "darkred",
                selectOp = sigOp,
            )
            self.addGroup("Wtaunu",
                members = self.getSafeListFromDataset(["WminustaunuPostVFP", "WplustaunuPostVFP"]),
                label = r"W$^{\pm}\to\tau\nu$",
                color = "orange",
                selectOp = sigOp,
            )
            self.addGroup("Top",
                members = list(filter(lambda y: y.group == "Top", self.datasets.values())),
                label = "Top",
                color = "green",
                selectOp = sigOp,
            )
            self.addGroup("Diboson",
                members = list(filter(lambda y: y.group == "Diboson", self.datasets.values())),
                label = "Diboson",
                color = "pink",
                selectOp = sigOp,
            )
            self.addGroup("QCD",
                members = list(filter(lambda y: y.group == "QCD", self.datasets.values())),
                label = "QCD MC",
                color = "grey",
                selectOp = sigOp,
            )
           
        else:
            self.addGroup("Other",
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["ZmumuPostVFP", "ZtautauPostVFP"] and x.group != "QCD"],
                label = "Other",
                color = "grey",
            )

        self.filterGroups(filterGroups)
        self.excludeGroups(excludeGroups)
        
        if self.wmass:
            # add all processes to the fake contributions after filtered and excluded groups
            self.addGroup("Fake",
                members = [member for sublist in [v.members for k, v in self.groups.items() if k != "QCD"] for member in sublist],
                scale = lambda x: 1. if x.is_data else -1,
                label = "Nonprompt",
                color = "grey",
                selectOp = fakeOp,
            )

        if self.wmass:
            self.gen_axes = ["etaGen", "ptGen"]
        elif self.wlike:
            self.gen_axes = ["qGen", "etaGen", "ptGen"]

