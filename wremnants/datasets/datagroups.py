from utilities import boostHistHelpers as hh
from wremnants import histselections as sel
from wremnants.datasets import datasets2016
import logging
import lz4.frame
import pickle
import narf
#import uproot
import ROOT
import re
import pandas as pd
import math
from utilities import common

logger = common.child_logger(__name__)

class datagroups(object):
    def __init__(self, infile, combine=False):
        self.combine = combine
        self.lumi = 1.
        if ".root" not in infile[-5:]:
            with lz4.frame.open(infile) as f:
                self.results = pickle.load(f)
            self.rtfile = None
        else:
            self.rtfile = ROOT.TFile.Open(infile)
            self.results = None

        self.lumi = None
        if self.datasets and self.results:
            self.data = [x for x in self.datasets.values() if x.is_data]
            if self.data:
                self.lumi = sum([self.results[x.name]["lumi"] for x in self.data if x.name in self.results])
        self.groups = {}

        if not self.lumi:
            self.lumi = 1
            
        self.nominalName = "nominal"

    def setNominalName(self, name):
        self.nominalName = name

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def getMetaInfo(self):
        return self.results["meta_info"]

    # for reading pickle files
    # as a reminder, the ND hists with tensor axes in the pickle files are organized as
    # pickle[procName]["output"][baseName] where
    ## procName are grouped into datagroups
    ## baseName takes values such as "nominal"
    def setHists(self, baseName, syst, procsToRead=None, label=None, nominalIfMissing=True, 
                 applySelection=True, forceNonzero=True, preOpMap=None, preOpArgs=None, scaleToNewLumi=-1):
        if not label:
            label = syst if syst else baseName
        if not procsToRead:
            procsToRead = self.groups.keys()

        foundExact = False
        for procName in procsToRead:
            logger.debug(f"Reading group {procName}")
            group = self.groups[procName]
            group[label] = None

            for member in group["members"]:
                logger.debug(f"Looking at group member {member.name}")
                scale = group["scale"] if "scale" in group else None
                try:
                    h = self.readHist(baseName, member, syst, scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    foundExact = True
                except ValueError as e:
                    if nominalIfMissing:
                        logger.info(f"{str(e)}. Using nominal hist {self.nominalName} instead")
                        h = self.readHist(self.nominalName, member, "", scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    else:
                        logger.warning(str(e))
                        continue

                if preOpMap and member.name in preOpMap:
                    logger.debug(f"Applying preOp to {member.name} after loading")
                    h = preOpMap[member.name](h, **preOpArgs)


                group[label] = h if not group[label] else hh.addHists(h, group[label])

            if not applySelection and "selectOp" in group and group["selectOp"]:
                logger.warning(f"Selection requested for process {procName} but applySelection=False, thus it will be ignored")
            if applySelection and group[label] and "selectOp" in group and group["selectOp"]:
                group[label] = group["selectOp"](group[label])
        # Avoid situation where the nominal is read for all processes for this syst
        if not foundExact:
            raise ValueError(f"Did not find systematic {syst} for any processes!")

    #TODO: Better organize to avoid duplicated code
    def setHistsCombine(self, baseName, syst, channel, procsToRead=None, excluded_procs=[], label=None):
        if type(excluded_procs) == str: excluded_procs = excluded_procs.split(",")
        #TODO Set axis names properly
        if baseName == "x":
            axisNames=["eta", "pt"]

        if not label:
            label = syst
        if not procsToRead:
            procsToRead = list(filter(lambda x: x not in excluded_procs, self.groups.keys()))

        for procName in procsToRead:
            group = self.groups[procName]
            group[label] = None
            if type(channel) == str: channel = channel.split(",")
            narf_hist = None
            for chn in channel:
                name = self.histNameCombine(procName, baseName, syst, chn)
                rthist = self.rtfile.Get(name)
                if not rthist:
                    raise RuntimeError(f"Failed to load hist {name} from file")
                if not narf_hist:
                    narf_hist = narf.root_to_hist(rthist, axis_names=axisNames)
                else:
                    narf_hist = hh.addHists(narf_hist, narf.root_to_hist(rthist, axis_names=axisNames))
            group[label] = narf_hist

    def histNameCombine(self, procName, baseName, syst, channel):
        name = f"{baseName}_{procName}"
        if syst != "nominal":
            name += "_"+syst
        if channel:
            name += "_"+channel
        if re.search("^pdf.*_sum", procName): # for pseudodata from alternative pdfset
            return("_".join([procName, channel])) 
        return name

    def loadHistsForDatagroups(self, baseName, syst, procsToRead=None, excluded_procs=None, channel="", label="", nominalIfMissing=True,
                               applySelection=True, forceNonzero=True, pseudodata=False, preOpMap={}, preOpArgs={}, scaleToNewLumi=-1):
        if self.rtfile and self.combine:
            self.setHistsCombine(baseName, syst, channel, procsToRead, excluded_procs, label)
        else:
            self.setHists(baseName, syst, procsToRead, label, nominalIfMissing, applySelection, forceNonzero, preOpMap, preOpArgs, scaleToNewLumi=scaleToNewLumi)

    def getDatagroups(self, excluded_procs=[]):
        if type(excluded_procs) == str:
            excluded_procs = list(excluded_procs)
        return dict(filter(lambda x: x[0] not in excluded_procs, self.groups.items()))

    def getNames(self, matches=[], exclude=False):
        if not matches:
            return list(x for x in self.groups.keys())
        else:
            if not exclude:
                return list(filter(lambda x: any([re.match(expr, x) for expr in matches]), self.groups.keys()))
            else:
                return list(filter(lambda x: x not in matches, self.groups.keys()))

    def getProcNames(self, to_expand=[], exclude_group=[]):
        procs = []
        if not to_expand:
            to_expand = self.groups.keys()
        for group_name in to_expand:
            if group_name not in exclude_group:
                for member in self.groups[group_name]["members"]:
                    procs.append(member.name)
        return procs

    def sortByYields(self, histName, nominalName="nominal"):
        def get_sum(h):
            return h.sum() if not hasattr(h.sum(), "value") else h.sum().value
        self.groups = dict(
            sorted(self.groups.items(), key=lambda x: get_sum(
                x[1][histName if histName in x[1] else nominalName])
                    if nominalName in x[1] or histName in x[1] else 0,
                reverse=True)
        )

    def getDatagroupsForHist(self, histName):
        filled = {}
        for k, v in self.groups.items():
            if histName in v:
                filled[k] = v
        return filled

    def resultsDict(self):
        return self.results

    def processes(self):
        return self.groups.keys()

    def addSummedProc(self, refname, name, label, color="red", exclude=["Data"], relabel=None, 
            procsToRead=None, reload=False, rename=None, action=None, preOpMap={}, preOpArgs={}):
        if reload:
            self.loadHistsForDatagroups(refname, syst=name, excluded_procs=exclude,
                procsToRead=procsToRead, preOpMap=preOpMap, preOpArgs=preOpArgs)

        if not rename:
            rename = name
        self.groups[rename] = dict(
            label=label,
            color=color,
            members=[],
        )
        tosum = []
        procs = procsToRead if procsToRead else self.groups.keys()
        for proc in filter(lambda x: x not in exclude+[rename], procs):
            h = self.groups[proc][name]
            if not h:
                raise ValueError(f"Failed to find hist for proc {proc}, histname {name}")
            if action:
                logger.debug(f"Applying action in addSummedProc! Before sum {h.sum()}")
                h = action(h)
                logger.debug(f"After action sum {h.sum()}")
            tosum.append(h)
        histname = refname if not relabel else relabel
        self.groups[rename][histname] = hh.sumHists(tosum)

    def copyWithAction(self, action, name, refproc, refname, label, color):
        self.groups[name] = dict(
            label=label,
            color=color,
            members=[],
        )
        self.groups[name][refname] = action(self.groups[refproc][refname])

    def setSelectOp(self, op, processes=None, exclude=False): 
        if processes == None:
            procs = self.getDatagroups()
        else:
            if exclude:
                procs = self.getDatagroups(excluded_procs=processes)
            else:
                procs = [processes] if isinstance(processes, str) else [p for p in processes]
        for proc in procs:
            if proc not in self.groups.keys():
                raise ValueError(f"In setSelectOp(): process {proc} not found")
            self.groups[proc]["selectOp"] = op
        
class datagroups2016(datagroups):
    def __init__(self, infile, combine=False, wlike=False, pseudodata_pdfset = None,
    ):
        self.datasets = {x.name : x for x in datasets2016.getDatasets()}
        super().__init__(infile, combine)
        if wlike:
            sigOp = None
            fakeOp = None
        else:
            sigOp = sel.signalHistWmass
            fakeOp = sel.fakeHistABCD
        ###
        self.hists = {} # container storing temporary histograms
        self.groups =  {
            "Data" : dict(
                members = [self.datasets["dataPostVFP"]],
                color = "black",
                label = "Data",
                selectOp = sigOp,
            ),
            "Zmumu" : dict(
                members = [self.datasets["ZmumuPostVFP"]],
                label = r"Z$\to\mu\mu$",
                color = "lightblue",
                selectOp = sigOp,
            ),   
            "Ztautau" : dict(
                members = [self.datasets["ZtautauPostVFP"]],
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
        if not wlike:
            self.groups.update({
                "Fake" : dict(
                    members = list(self.datasets.values()),
                    scale = lambda x: 1. if x.is_data else -1,
                    label = "Nonprompt",
                    color = "grey",
                    selectOp = fakeOp,
                ),
                "Wtau" : dict(
                    members = [self.datasets["WminustaunuPostVFP"], self.datasets["WplustaunuPostVFP"]],
                    label = r"W$^{\pm}\to\tau\nu$",
                    color = "orange",
                    selectOp = sigOp,
                ),
                "Wmunu" : dict(
                    members = [self.datasets["WminusmunuPostVFP"], self.datasets["WplusmunuPostVFP"]],
                    label = r"W$^{\pm}\to\mu\nu$",
                    color = "darkred",
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
            })
        else:
            self.groups["Other"] = dict(
                members = [x for x in self.datasets.values() if not x.is_data and x.name not in ["ZmumuPostVFP", "ZtautauPostVFP"]],
                label = "Other",
                color = "grey",
            )

    def make_yields_df(self, histName, procs, action):
        def sum_and_unc(h):
            return (h.sum().value, math.sqrt(h.sum().variance))
        df = pd.DataFrame([(k, *sum_and_unc(action(v[histName]))) for k,v in self.groups.items() if k in procs], 
                columns=["Process", "Yield", "Uncertainty"])
        return df

    def histName(self, baseName, procName, syst):
        # This is kind of hacky to deal with the different naming from combine
        if baseName != "x" and (syst == "" or syst == self.nominalName):
            return baseName
        if baseName in ["", "x", "nominal"] and syst:
            return syst
        if syst[:len(baseName)] == baseName:
            return syst
        return "_".join([baseName,syst])
    
    def readHist(self, baseName, proc, syst, scaleOp=None, forceNonzero=True, scaleToNewLumi=-1):
        print("PROC.NAME")
        print(proc.name)
        print("BASENAME")
        print(baseName)
        print("SYST")
        print(syst)
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        print(histname)
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        if forceNonzero:
            h = hh.clipNegativeVals(h)
        if scaleToNewLumi > 0:
            h = hh.scaleByLumi(h, scaleToNewLumi)
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale
