from utilities import boostHistHelpers as hh,common,logging
from wremnants import histselections as sel
from wremnants.datasets import datasets2016
import lz4.frame
import pickle
import hdf5plugin
import h5py
import narf
import hist
import ROOT
import re
import pandas as pd
import math
import os
import itertools

logger = logging.child_logger(__name__)

class datagroups(object):
    def __init__(self, infile, combine=False):
        self.combine = combine
        self.h5file = None
        self.rtfile = None
        if infile.endswith(".pkl.lz4"):
            with lz4.frame.open(infile) as f:
                self.results = pickle.load(f)
        elif infile.endswith(".hdf5"):
            self.h5file = h5py.File(infile, "r")
            self.results = narf.ioutils.pickle_load_h5py(self.h5file["results"])
        elif infile.endswith(".root"):
            self.rtfile = ROOT.TFile.Open(infile)
            self.results = None
        else:
            raise ValueError("Unsupported file type")

        self.wmass = os.path.basename(self.getScriptCommand().split()[0]).startswith("mw")
        self.wlike = os.path.basename(self.getScriptCommand().split()[0]).startswith("mz_wlike")

        self.lumi = None
        # FIXME: self.datasets is currently a data member of the inherited class, we should define it here as well
        # FIXME: if data is excluded the normalization will be lost
        if self.datasets and self.results:
            self.data = [x for x in self.datasets.values() if x.is_data]
            if self.data:
                self.lumi = sum([self.results[x.name]["lumi"] for x in self.data if x.name in self.results])
        self.groups = {}
        self.groupNamesPostFilter = [] # keep track of group names after any filter defined by the user

        self.unconstrainedProcesses = [] # processes that are treated as free floating (usually signals)

        if not self.lumi:
            logger.warning("")
            logger.warning("*"*30)
            logger.warning("No data sample selected: setting integrated luminosity to 1/fb")
            logger.warning("*"*30)
            logger.warning("")
            self.lumi = 1
            
        self.nominalName = "nominal"
        self.globalAction = None

    def __del__(self):
        if self.h5file:
            self.h5file.close()
        if self.rtfile:
            self.rtfile.Close()

    # To be used for applying a selection, rebinning, etc.
    def setGlobalAction(self, action):
        self.globalAction = action

    def setNominalName(self, name):
        self.nominalName = name

    def processScaleFactor(self, proc):
        if proc.is_data:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def getMetaInfo(self):
        if self.results:
            return self.results["meta_info"] if "meta_info" in self.results else self.results["meta_data"]
        raise NotImplementedError("Currently can't access meta data as dict for ROOT file")

    def getScriptCommand(self):
        if self.rtfile:
            return self.rtfile.Get("meta_info/command").GetTitle()
        else:
            meta_info = self.results["meta_info"] if "meta_info" in self.results else self.results["meta_data"]
            return meta_info["command"]

    def updateGroupNamesPostFilter(self, excludeGroup=[]):
        self.groupNamesPostFilter = list(x for x in self.groups.keys() if len(self.groups[x]["members"]) and x not in excludeGroup)
        logger.debug(f"Updating group names, new list is {self.groupNamesPostFilter}")
        
    # for reading pickle files
    # as a reminder, the ND hists with tensor axes in the pickle files are organized as
    # pickle[procName]["output"][baseName] where
    ## procName are grouped into datagroups
    ## baseName takes values such as "nominal"
    def setHists(self, baseName, syst, procsToRead=None, label=None, nominalIfMissing=True, 
                 applySelection=True, forceNonzero=True, preOpMap=None, preOpArgs=None, scaleToNewLumi=-1, 
                 excludeProcs=None, forceToNominal=[]):
        if not label:
            label = syst if syst else baseName
        logger.debug(f"In setHists(): procsToRead = {procsToRead}")

        if not procsToRead:
            if excludeProcs:
                procsToRead = list(filter(lambda x: x not in excludeProcs, self.groups.keys()))
            else:
                procsToRead = list(self.groups.keys())

        foundExact = False
        for procName in procsToRead:
            logger.debug(f"Reading group {procName}")
            group = self.groups[procName] if procName in self.groups else {}
            group[label] = None

            for member in group["members"]:
                logger.debug(f"Looking at group member {member.name}")
                scale = group["scale"] if "scale" in group else None
                read_syst = syst
                if member.name in forceToNominal:
                    read_syst = ""
                    logger.debug(f"Forcing group member {member.name} to read the nominal hist for syst {syst}")

                try:
                    h = self.readHist(baseName, member, read_syst, scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    foundExact = True
                except ValueError as e:
                    if nominalIfMissing:
                        logger.info(f"{str(e)}. Using nominal hist {self.nominalName} instead")
                        h = self.readHist(self.nominalName, member, "", scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    else:
                        logger.warning(str(e))
                        continue
                logger.debug(f"Hist axes are {h.axes.name}")

                if preOpMap and member.name in preOpMap:
                    logger.debug(f"Applying preOp to {member.name} after loading")
                    h = preOpMap[member.name](h, **preOpArgs)

                if self.globalAction:
                    h = self.globalAction(h)

                group[label] = h if not group[label] else hh.addHists(h, group[label])

            # Can use to apply common rebinning or selection on top of the usual one
            if "rebinOp" in group and group["rebinOp"]:
                group[label] = group["rebinOp"](group[label])

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
            if excludeProcs:
                procsToRead = list(filter(lambda x: x not in excludeProcs, self.groups.keys()))
            else:
                procsToRead = list(self.groups.keys())

        for procName in procsToRead:
            group = self.groups[procName] if procName in self.groups else {}
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

            if self.globalAction:
                narf_hist = self.globalAction(narf_hist)

            group[label] = narf_hist

    def histName(self, baseName, procName="", syst=""):
        return datagroups.histName(baseName, procName, syst, nominalName=self.nominalName)

    def histNameCombine(self, procName, baseName, syst, channel):
        return datagroups.histNameCombine(procName, baseName, syst, channel)

    def loadHistsForDatagroups(
        self, baseName, syst, procsToRead=None, excluded_procs=None, channel="", label="",
        nominalIfMissing=True, applySelection=True, forceNonzero=True, pseudodata=False,
        preOpMap={}, preOpArgs={}, scaleToNewLumi=-1, forceToNominal=[]
    ):
        logger.debug("Calling loadHistsForDatagroups()")
        logger.debug(f"the basename and syst is: {baseName}, {syst}")
        logger.debug(f"The procsToRead and excludedProcs are: {procsToRead}, {excluded_procs}")
        if self.rtfile and self.combine:
            self.setHistsCombine(baseName, syst, channel, procsToRead, excluded_procs, label)
        else:
            self.setHists(baseName, syst, procsToRead, label, nominalIfMissing, applySelection,
                          forceNonzero, preOpMap, preOpArgs,
                          scaleToNewLumi=scaleToNewLumi, 
                          excludeProcs=excluded_procs, forceToNominal=forceToNominal)

    def addGroup(self, keyname, dictToAdd, canReplaceKey=False):
        if canReplaceKey or keyname not in self.groups.keys():
            if keyname in self.groups.keys():
                logger.warning(f"Replacing {keyname} in groups")
            self.groups[keyname] = dictToAdd
            if keyname not in self.groupNamesPostFilter:
                self.groupNamesPostFilter.append(keyname)
            
    def deleteGroup(self, procs):
        if isinstance(procs, str):
            procs = [procs,]

        for p in procs:
            if p in self.groups.keys():
                del self.groups[p]
            if p in self.groupNamesPostFilter:
                self.groupNamesPostFilter.remove(p)
            
    def getDatagroups(self, excluded_procs=[], afterFilter=True):
        # usually excluded_procs is not needed if afterFilter=True, unless one has to filter further
        if type(excluded_procs) == str:
            excluded_procs = list(excluded_procs)
        filtDef = lambda x: x[0] not in excluded_procs
        if afterFilter:
            filtDef = lambda x: x[0] in self.groupNamesPostFilter
            if len(excluded_procs):
                filtDef = lambda x: x[0] in self.groupNamesPostFilter and x[0] not in excluded_procs
        return dict(filter(filtDef, self.groups.items()))

    # INFO: this method returns the list from the full set of defined groups, unless one filters further.
    # Instead, argument 'afterFilter' is used to return the names after the filter passed to the constructor
    def getNames(self, matches=[], exclude=False, afterFilter=True):
        listOfNames = list(x for x in self.groupNamesPostFilter) if afterFilter else list(x for x in self.groups.keys())
        if not matches:
            return listOfNames
        else:
            # matches uses regular expressions with search (and can be inverted when exclude is true),
            # thus a string will match if the process name contains that string anywhere inside it
            ##########
            # FIXME ? : allow for usage of simple 'string in name' syntax, with no regular expressions? Or exact names?
            #           Note that datasets2016.getDatasets currently accepts only exact names, so one should stay consistent
            ##########
            if exclude:
                return list(filter(lambda x: all([re.search(expr, x) is None for expr in matches]), listOfNames))
            else:
                return list(filter(lambda x: any([re.search(expr, x) for expr in matches]), listOfNames))
              
    def getProcNames(self, to_expand=[], exclude_group=[], afterFilter=True):
        procs = []
        if not to_expand:
            to_expand = self.groupNamesPostFilter if afterFilter else self.groups.keys()
        for group_name in to_expand:
            if afterFilter and group_name not in self.groupNamesPostFilter:
                continue
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

    def defineSignalBinsUnfolding(self, fitvar, base_process, add_overflow=False):
        # get gen bin names corresponding to fitvars
        genvar_dict = {
            "pt": "ptGen",
            "eta": "etaGen",
        }

        fitvars = fitvar.split("-")
        genvars = []
        gen_bins = []
        for fitvar in fitvars:
            if fitvar not in genvar_dict.keys():
                raise RuntimeError(f"No corresponding gen level definition for {fitvar} found!")
            genvar = genvar_dict[fitvar]
            genvars.append(genvar)
            gen_bin_edges = self.results[self.groups[base_process]["members"][0].name]["output"]["gen"].get().axes[genvar].edges

            if add_overflow:
                # to add separate overflow/underflow bins
                gen_bins.append([hist.underflow, *[i for i in range(len(gen_bin_edges)-1)], hist.overflow]) 
            else:
                gen_bins.append(range(len(gen_bin_edges)-1))

        for indices in itertools.product(*gen_bins):
            
            proc_genbin = dict(self.groups[base_process])
            proc_genbin['selectOp'] = lambda x, indices=indices, genvars=genvars: x[{var : i for var, i in zip(genvars, indices)}]

            proc_name = base_process
            for idx, var in zip(indices, fitvars):
                if idx == hist.underflow:
                    proc_name += f"_{var}U"
                elif idx == hist.overflow:
                    proc_name += f"_{var}O"
                else:
                    proc_name += f"_{var}{idx}"

            self.addGroup(proc_name, proc_genbin)
            self.unconstrainedProcesses.append(proc_name)

        # Remove inclusive signal
        self.deleteGroup(base_process)


    @staticmethod
    def histName(baseName, procName="", syst=""):
        # This is kind of hacky to deal with the different naming from combine
        if baseName != "x" and (syst == ""):
            return baseName
        if baseName in ["", "x"] and syst:
            return syst
        if syst[:len(baseName)] == baseName:
            return syst
        return "_".join([baseName,syst])
    
    @staticmethod
    def histNameCombine(procName, baseName, syst, channel):
        name = f"{baseName}_{procName}"
        if syst != "nominal":
            name += "_"+syst
        if channel:
            name += "_"+channel
        if re.search("^pdf.*_sum", procName): # for pseudodata from alternative pdfset
            return("_".join([procName, channel])) 
        return name

class datagroups2016(datagroups):
    def __init__(self, infile, combine=False, pseudodata_pdfset = None, applySelection=True,
                 excludeProcGroup=None, filterProcGroup=None
    ):
        self.datasets = {x.name : x for x in datasets2016.getDatasets(filt=filterProcGroup, excl=excludeProcGroup)}
        logger.debug(f"Getting these datasets: {self.datasets.keys()}")
        super().__init__(infile, combine)
        if self.wmass and applySelection:
            sigOp = sel.signalHistWmass
            fakeOp = sel.fakeHistABCD
        else:
            sigOp = None
            fakeOp = None

        ###
        self.hists = {} # container storing temporary histograms
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
        if excludeProcGroup is None:
            excludeProcGroup = []

        self.updateGroupNamesPostFilter(excludeGroup=excludeProcGroup)
        #self.groupNamesPostFilter = list(x for x in self.groups.keys() if len(self.groups[x]["members"]) and x not in excludeProcGroup)
        logger.debug(f"Filtered groups: {self.groupNamesPostFilter}")

    # TODO: move to base class
    def getSafeListFromDataset(self, procs):
        # return list of valid samples which belongs to the dataset or where not excluded elsewhere
        if isinstance(procs, str):
            return [self.datasets[procs]] if procs in self.datasets.keys() else []
        else:
            return list(self.datasets[x] for x in procs if x in self.datasets.keys())
        
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
