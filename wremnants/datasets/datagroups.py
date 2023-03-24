from utilities import boostHistHelpers as hh,common,logging
import lz4.frame
import pickle
import h5py
import narf
import ROOT
import re
import os
import itertools
import functools
import hist

logger = logging.child_logger(__name__)

class Datagroups(object):
    def __init__(self, infile, combine=False, datasets=None):
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

        self.lumi = 1

        if datasets:
            self.datasets = {x.name : x for x in datasets}
            logger.debug(f"Getting these datasets: {self.datasets.keys()}")

            if self.results:
                self.data = [x for x in self.datasets.values() if x.is_data]
                if self.data:
                    self.lumi = sum([self.results[x.name]["lumi"] for x in self.data if x.name in self.results])

        self.groups = {}
        self.nominalName = "nominal"
        self.globalAction = None
        self.unconstrainedProcesses = []

    def __del__(self):
        if self.h5file:
            self.h5file.close()
        if self.rtfile:
            self.rtfile.Close()

    def addGroup(self, group_name, dictToAdd, canReplaceKey=False):
        if canReplaceKey or group_name not in self.groups.keys():
            if group_name in self.groups.keys():
                logger.warning(f"Replacing {group_name} in groups")
            self.groups[group_name] = dictToAdd

    def deleteGroup(self, procs):
        if isinstance(procs, str):
            procs = [procs,]

        for p in procs:
            if p in self.groups.keys():
                del self.groups[p]

    def addGroupMember(self, group_name, member_name):
        # adds a process to the existing members of a given group
        if group_name not in self.groups:
            logger.warning(f"The group {group_name} is not defined in the datagroups object! Do nothing here.")
            return

        if self.datasets and member_name not in [d for d in self.datasets]:
            logger.warning(f"The member {member_name} can not be found in the current dataset, still add it to the list of members of the group and trust you.")

        self.groups[group_name]["members"] = [*self.groups[group_name]["members"], member_name]

    def deleteGroupMember(self, group_name, member_name):
        # deletes a process from the list of members of a given group
        if group_name not in self.groups:
            logger.warning(f"The group {group_name} is not defined in the datagroups object! Do nothing here.")
            return

        if member_name not in self.groups[group_name]["members"]:
            logger.warning(f"The member {member_name} can not be found in the group {group_name}! Do nothing here.")
            return

        self.groups[group_name]["members"] = [m for m in self.groups[group_name]["members"] if m != member_name]

    def selectGroups(self, selections):
        new_groupnames = []
        for selection in selections:
            new_groupnames += list(filter(lambda x, s=selection: x == s, self.groups.keys()))

        # remove duplicates selected by multiple filters
        return list(set(new_groupnames))

    def filterGroups(self, filters):
        if filters:
            if isinstance(filters, str):
                filters = [filters]

            if isinstance(filters, list):
                new_groupnames = self.selectGroups(filters)
            else:
                new_groupnames = list(filter(filters, self.groups.keys()))

            diff = list(set(new_groupnames) - set(self.groups.keys()))
            if diff:
                logger.debug(f"Filtered out following groups: {diff}")

            self.groups = {key: self.groups[key] for key in new_groupnames}

        if len(self.groups) == 0:
            logger.warning("Filtered groups but didn't find any match. Continue without any group.")

    def excludeGroups(self, excludes):
        if excludes:
            if isinstance(excludes, str):
                excludes = [excludes]

            if isinstance(excludes, list):
                # remove selected datasets
                new_groupnames = list(filter(lambda x: x not in self.selectGroups(excludes), self.groups))
            else:
                new_groupnames = list(filter(excludes, self.groups.keys()))

            diff = list(set(new_groupnames) - set(self.groups.keys()))
            if diff:
                logger.debug(f"Filtered out following groups: {diff}")

            self.groups = {key: self.groups[key] for key in new_groupnames}
        
        if len(self.groups) == 0:
            logger.warning("Excluded all groups. Continue without any group.")

    def getSafeListFromDataset(self, procs):
        # return list of valid samples which belongs to the dataset or where not excluded elsewhere
        if isinstance(procs, str):
            return [self.datasets[procs]] if procs in self.datasets.keys() else []
        else:
            return list(self.datasets[x] for x in procs if x in self.datasets.keys())
    
    def setGlobalAction(self, action):
        # To be used for applying a selection, rebinning, etc.
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
    def setHistsCombine(self, baseName, syst, channel, procsToRead=None, excludeProcs=[], label=None):
        if type(excludeProcs) == str: excludeProcs = excludeProcs.split(",")
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
        return Datagroups.histName(baseName, procName, syst, nominalName=self.nominalName)

    def histNameCombine(self, procName, baseName, syst, channel):
        return Datagroups.histNameCombine(procName, baseName, syst, channel)

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

    def getNames(self, matches=[], exclude=False):
        # This method returns the names from the defined groups, unless one selects further.
        listOfNames = list(x for x in self.groups.keys())
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

    def setSelectOp(self, op, processes=None): 
        if processes == None:
            procs = self.groups
        else:
            procs = [processes] if isinstance(processes, str) else processes

        for proc in procs:
            if proc not in self.groups.keys():
                raise ValueError(f"In setSelectOp(): process {proc} not found")
            self.groups[proc]["selectOp"] = op

    def defineSignalBinsUnfolding(self, fitvar, base_process):
        # get gen bin names corresponding to fitvars
        genvar_dict = {
            "pt": "ptGen",
            "eta": "etaGen"
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

            # # to add out of acceptance constributions for separate overflow/underflow bins
            # gen_bins.append([hist.underflow, *[i for i in range(len(gen_bin_edges)-1)], hist.overflow]) 

            gen_bins.append(range(len(gen_bin_edges)-1))

        for indices in itertools.product(*gen_bins):
            proc_genbin = dict(self.groups[base_process])
            if proc_genbin['selectOp'] is None:
                proc_genbin['selectOp'] = lambda x, indices=indices, genvars=genvars: x[{var : i for var, i in zip(genvars, indices)}]            
            else:
                # in case there is already a selection operation, e.g. for W, this one needs to be wrapped around the first
                proc_genbin['selectOp'] = lambda x, indices=indices, genvars=genvars, f0=proc_genbin['selectOp']: f0(x)[{var : i for var, i in zip(genvars, indices)}]

            proc_name = base_process
            for idx, var in zip(indices, fitvars):
                if idx == hist.underflow:
                    proc_name += f"_{var}U"
                elif idx == hist.overflow:
                    proc_name += f"_{var}O"
                else:
                    proc_name += f"_{var}{idx}"
            
            self.unconstrainedProcesses.append(proc_name)
            self.addGroup(proc_name, proc_genbin)
            
            if self.wmass:
                # Add fake distribution of gen level bin
                fake_genbin = dict(self.groups["Fake"])
                fake_genbin['selectOp'] = lambda x, indices=indices, genvars=genvars, f0=fake_genbin['selectOp']: f0(x)[{var : i for var, i in zip(genvars, indices)}]
                proc_name += "_Fake"
                self.addGroup(proc_name, fake_genbin)

            # self.addGroupMember("Fake", proc_name)

        # add one inclusive out of acceptance contribution and treat as background
        proc_genbin = dict(self.groups[base_process])

        conditions = [{var : hist.overflow} for var in genvars] + [{var : hist.underflow} for var in genvars]
        combined_condition = functools.reduce(lambda a, b: a | b, conditions)
        proc_genbin['selectOp'] = lambda x, f0=proc_genbin['selectOp']: f0(x)[combined_condition]

        self.addGroup(base_process, proc_genbin)

        # Remove inclusive signal
        for member in self.groups[base_process]["members"]:
            self.deleteGroupMember("Fake", member)
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
