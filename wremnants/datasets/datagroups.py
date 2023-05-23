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
import pandas as pd
import math

from wremnants.datasets.datagroup import Datagroup

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
                else:
                    logger.warning("No data process was selected, normalizing MC to to 1/fb")

        self.groups = {}
        self.nominalName = "nominal"
        self.globalAction = None
        self.gen_axes = None
        self.unconstrainedProcesses = []

    def __del__(self):
        if self.h5file:
            self.h5file.close()
        if self.rtfile:
            self.rtfile.Close()

    def addGroup(self, name, **kwargs):
        group = Datagroup(name, **kwargs)
        self.groups[name] = group
        
    def deleteGroups(self, names):
        for n in names:
            self.deleteGroup(n)

    def deleteGroup(self, name):
        if name in self.groups.keys():
            del self.groups[name]
        else:
            logger.warning(f"Try to delete group '{name}' but did not find this group.")

    def copyGroup(self, group_name, new_name, member_filter=None):
        self.groups[new_name] = self.groups[group_name].copy(new_name, member_filter)

    def selectGroups(self, selections):
        new_groupnames = []
        for selection in selections:
            new_groupnames += list(filter(lambda x, s=selection: x == s, self.groups.keys()))

        # remove duplicates selected by multiple filters
        return list(set(new_groupnames))

    def filterGroups(self, filters):
        if filters is None:
            return

        if isinstance(filters, str):
            filters = [filters]

        if isinstance(filters, list):
            new_groupnames = self.selectGroups(filters)
        else:
            new_groupnames = list(filter(filters, self.groups.keys()))

        diff = list(self.groups.keys() - set(new_groupnames))
        if diff:
            logger.info(f"Datagroups.filterGroups : filtered out following groups: {diff}")

        self.groups = {key: self.groups[key] for key in new_groupnames}

        if len(self.groups) == 0:
            logger.warning(f"Filtered groups using '{filters}' but didn't find any match. Continue without any group.")

    def excludeGroups(self, excludes):
        if excludes is None:
            return

        if isinstance(excludes, str):
            excludes = [excludes]

        if isinstance(excludes, list):
            # remove selected datasets
            new_groupnames = list(filter(lambda x: x not in self.selectGroups(excludes), self.groups))
        else:
            new_groupnames = list(filter(excludes, self.groups.keys()))

        diff = list(self.groups.keys() - set(new_groupnames))
        if diff:
            logger.info(f"Datagroups.excludeGroups: filtered out following groups: {diff}")

        self.groups = {key: self.groups[key] for key in new_groupnames}
        
        if len(self.groups) == 0:
            logger.warning(f"Excluded all groups using '{excludes}'. Continue without any group.")

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
        logger.info(f"In setHists(): for hist {syst} procsToRead = {procsToRead}")

        if not procsToRead:
            if excludeProcs:
                procsToRead = list(filter(lambda x: x not in excludeProcs, self.groups.keys()))
            else:
                procsToRead = list(self.groups.keys())

        foundExact = False

        # If fakes are present do them as last group, and when running on prompt group build the sum to be used for the fakes.
        # This makes the code faster and avoid possible bugs related to reading again the same processes
        # NOTE:
        # To speed up even more, one could directly use the per-group sum already computed for each group,
        # but this would need to assume that fakes effectively had all the single processes in each group as members
        # (usually it will be the case, but it is more difficult to handle in a fully general way and without bugs)
        histForFake = None # to store the data-MC sums used for the fakes, for each syst
        nameFake = "Fake" # TODO: actual name might/should be configurable
        if nameFake in procsToRead:
            procsToReadSort = [x for x in procsToRead if x != nameFake] + [nameFake]
            hasFake = True
            fakesMembers = [m.name for m in self.groups[nameFake].members]
            fakesMembersWithSyst = []
        else:
            hasFake = False
            procsToReadSort = [x for x in procsToRead]
        # Note: if 'hasFake' is kept as False (but Fake exists), the original behaviour for which Fake reads everything again is restored
            
        for procName in procsToReadSort:
            logger.debug(f"Reading group {procName}")
            if procName not in self.groups.keys():
                raise RuntimeError(f"Group {procName} not known. Defined groups are {list(self.groups.keys())}.")
            group = self.groups[procName]
            group.hists[label] = None

            for i, member in enumerate(group.members):
                if procName == nameFake and member.name in fakesMembersWithSyst:
                    # if we are here this process has been already use to build the fakes when running for other groups
                    continue
                logger.debug(f"Looking at group member {member.name}")
                scale = group.scale
                read_syst = syst
                if member.name in forceToNominal:
                    read_syst = ""
                    logger.debug(f"Forcing group member {member.name} to read the nominal hist for syst {syst}")

                try:
                    h = self.readHist(baseName, member, procName, read_syst, scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    foundExact = True
                except ValueError as e:
                    if nominalIfMissing:
                        logger.info(f"{str(e)}. Using nominal hist {self.nominalName} instead")
                        h = self.readHist(self.nominalName, member, procName, "", scaleOp=scale, forceNonzero=forceNonzero, scaleToNewLumi=scaleToNewLumi)
                    else:
                        logger.warning(str(e))
                        continue
                logger.debug(f"Hist axes are {h.axes.name}")

                ## TODO/FIXME: is the order of these operations important?
                ## Some operations are triggered based on the group and some might be specific for fakes,
                ## in which case they should be exported inside the "if" below, since they wouldn't be run here.
                ## Might have a function to contain all these operations to avoid code repetition,
                ## but not all are needed, so we should first assess whether the order is relevant before reshuffling
                if group.memberOp:
                    if group.memberOp[i] is not None:
                        logger.debug(f"Apply operation to member {i}: {member.name}/{procName}")
                        h = group.memberOp[i](h)
                    else:
                        logger.debug(f"No operation for member {i}: {member.name}/{procName}")

                if self.gen_axes:
                    # integrate over remaining gen axes 
                    projections = [a for a in h.axes.name if a not in self.gen_axes]
                    if len(projections) < len(h.axes.name):
                        h = h.project(*projections)

                if preOpMap and member.name in preOpMap:
                    logger.debug(f"Applying preOp to {member.name}/{procName} after loading")
                    h = preOpMap[member.name](h, **preOpArgs)
                    
                if self.globalAction:
                    h = self.globalAction(h)
                    logger.debug("Applying global action")

                hasPartialSumForFake = False
                if hasFake and procName != nameFake:
                    if member.name in fakesMembers:
                        if member.name not in fakesMembersWithSyst:
                            fakesMembersWithSyst.append(member.name)
                        hasPartialSumForFake = True
                        # apply the correct scale for fakes
                        scaleProcForFake = self.groups[nameFake].scale(member)
                        logger.debug(f"Summing hist {read_syst} for {member.name} to {nameFake} with scale = {scaleProcForFake}")
                        hProcForFake = scaleProcForFake * h.copy()
                        histForFake = hh.addHists(hProcForFake, histForFake, createNew=False) if histForFake else hProcForFake

                # The following must be done when the group is not Fake, or when the previous part for fakes was not done
                # For fake this essentially happens when the process doesn't have the syst, so that the nominal is used
                if procName != nameFake or (procName == nameFake and not hasPartialSumForFake):
                    if procName == nameFake:
                        logger.debug(f"Summing nominal hist instead of {syst} to {nameFake} for {member.name}")
                    else:
                        logger.debug(f"Summing {read_syst} to {procName} for {member.name}")
                    group.hists[label] = hh.addHists(h, group.hists[label], createNew=False) if group.hists[label] else h
                    logger.debug("Sum done")
                
            # now sum to fakes the partial sums which where not already done before
            # (group.hists[label] contains only the contribution from nominal histograms).
            # Then continue with the rest of the code as usual
            if hasFake and procName == nameFake:
                if histForFake is not None:
                    group.hists[label] = hh.addHists(histForFake, group.hists[label], createNew=False) if group.hists[label] else histForFake

            # Can use to apply common rebinning or selection on top of the usual one
            if group.rebinOp:
                group.hists[label] = group.rebinOp(group.hists[label])

            if group.selectOp:
                if not applySelection:
                    logger.warning(f"Selection requested for process {procName} but applySelection=False, thus it will be ignored")
                elif label in group.hists.keys():
                    group.hists[label] = group.selectOp(group.hists[label], **group.selectOpArgs)
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
            group.hists[label] = None
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

            group.hists[label] = narf_hist

    def loadHistsForDatagroups(
        self, baseName, syst, procsToRead=None, excluded_procs=None, channel="", label="",
        nominalIfMissing=True, applySelection=True, forceNonzero=True, pseudodata=False,
        preOpMap={}, preOpArgs={}, scaleToNewLumi=-1, forceToNominal=[]
    ):
        logger.debug("Calling loadHistsForDatagroups()")
        logger.debug(f"The basename and syst is: {baseName}, {syst}")
        logger.debug(f"The procsToRead and excludedProcs are: {procsToRead}, {excluded_procs}")
        if self.rtfile and self.combine:
            self.setHistsCombine(baseName, syst, channel, procsToRead, excluded_procs, label)
        else:
            self.setHists(baseName, syst, procsToRead, label, nominalIfMissing, applySelection,
                          forceNonzero, preOpMap, preOpArgs,
                          scaleToNewLumi=scaleToNewLumi, 
                          excludeProcs=excluded_procs, forceToNominal=forceToNominal)

    def getDatagroups(self):
        return self.groups

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
                for member in self.groups[group_name].members:
                    # protection against duplicates in the output list, they may arise from fakes
                    if member.name not in procs:
                        procs.append(member.name)
        return procs

    def sortByYields(self, histName, nominalName="nominal"):
        def get_sum(h):
            return h.sum() if not hasattr(h.sum(), "value") else h.sum().value
        self.groups = dict(
            sorted(self.groups.items(), key=lambda x: get_sum(
                x[1].hists[histName if histName in x[1].hists else nominalName])
                    if nominalName in x[1].hists or histName in x[1].hists else 0,
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
        self.addGroup(rename,
            label=label,
            color=color,
            members=[],
        )
        tosum = []
        procs = procsToRead if procsToRead else self.groups.keys()
        for proc in filter(lambda x: x not in exclude+[rename], procs):
            h = self.groups[proc].hists[name]
            if not h:
                raise ValueError(f"Failed to find hist for proc {proc}, histname {name}")
            if action:
                logger.debug(f"Applying action in addSummedProc! Before sum {h.sum()}")
                h = action(h)
                logger.debug(f"After action sum {h.sum()}")
            tosum.append(h)
        histname = refname if not relabel else relabel
        self.groups[rename].hists[histname] = hh.sumHists(tosum)

    def setSelectOp(self, op, processes=None): 
        if processes == None:
            procs = self.groups
        else:
            procs = [processes] if isinstance(processes, str) else processes

        for proc in procs:
            if proc not in self.groups.keys():
                raise ValueError(f"In setSelectOp(): process {proc} not found")
            self.groups[proc].selectOp = op

    def setGenAxes(self, gen_axes):
        if isinstance(gen_axes, str):
            self.gen_axes = [gen_axes,]
        else:
            self.gen_axes = gen_axes

    def defineSignalBinsUnfolding(self, group_name):
        if group_name not in self.groups.keys():
            raise RuntimeError(f"Base group {group_name} not found in groups {self.groups.keys()}!")

        nominal_hist = self.results[self.groups[group_name].members[0].name]["output"]["xnorm"].get()

        gen_bins = []
        for gen_axis in self.gen_axes:
            if gen_axis not in nominal_hist.axes.name:
                raise RuntimeError(f"Gen axis '{gen_axis}' not found in histogram axes '{nominal_hist.axes.name}'!")

            gen_bin_edges = nominal_hist.axes[gen_axis].edges
            gen_bins.append(range(len(gen_bin_edges)-1))

        base_members = self.groups[group_name].members[:]

        for indices in itertools.product(*gen_bins):

            proc_name = group_name
            for idx, var in zip(indices, self.gen_axes):
                proc_name += f"_{var}{idx}"

            self.copyGroup(group_name, proc_name)

            memberOp = lambda x, indices=indices, genvars=self.gen_axes: x[{var : i for var, i in zip(self.gen_axes, indices)}]
            self.groups[proc_name].memberOp = [memberOp for m in base_members]

            self.unconstrainedProcesses.append(proc_name)

        # add one out of acceptance group and treat as background
        ooa_name = group_name.split("_")[0]+"_bkg"
        if ooa_name not in self.groups.keys():
            self.copyGroup(group_name, ooa_name)
            self.groups[ooa_name].memberOp = []
            self.groups[ooa_name].members = []
        
        # list of possible slices for each axis
        slices = []
        for var in self.gen_axes:
            subslice = [{var : hist.tag.Slicer()[0:hist.overflow:hist.sum]}] 
            if nominal_hist.axes[var].traits.__dict__["underflow"]:
                subslice.append({var : hist.underflow})
            if nominal_hist.axes[var].traits.__dict__["overflow"]:
                subslice.append({var : hist.overflow})
            slices.append(subslice)

        # pick one slice for each axis from the list of possible slices
        for condition in [functools.reduce(lambda x, y: {**x, **y}, tup) for tup in itertools.product(*slices)]:
            # make sure that at least one axis takes the underflow/overflow, otherwise it would not be out of acceptance
            if not any([c in [hist.underflow, hist.overflow] for c in condition.values()]):
                continue
            logger.debug(f"Add members with condition {condition}")
            self.groups[ooa_name].addMembers(base_members, lambda x, c=condition: x[c])

        # Remove inclusive signal
        self.deleteGroup(group_name)

    def make_yields_df(self, histName, procs, action=lambda x: x, norm_proc=None):
        def sum_and_unc(h):
            if not hasattr(h.sum(), "value"):
                return (h.sum(), None)
            else:
                return (h.sum().value, math.sqrt(h.sum().variance))

        df = pd.DataFrame([(k, *sum_and_unc(action(v.hists[histName]))) for k,v in self.groups.items() if k in procs], 
                columns=["Process", "Yield", "Uncertainty"])

        if norm_proc and norm_proc in self.groups:
            hist = action(self.groups[norm_proc].hists[histName])
            denom = hist.sum() if not hasattr(hist.sum(), "value") else hist.sum().value
            df[f"Ratio to {norm_proc} (%)"] = df["Yield"]/denom*100
            
        return df

    def readHist(self, baseName, proc, group, syst, scaleOp=None, forceNonzero=True, scaleToNewLumi=-1):
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        logger.debug(f"Reading hist {histname} for proc/group {proc.name}/{group} and syst '{syst}'")
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")
        h = output[histname]
        if isinstance(h, narf.ioutils.H5PickleProxy):
            h = h.get()
        # Do a copy to detach the modified object from the original one, so the latter stays unmodified.
        # This is extremely important for fakes, since they reuse the other histograms to subtract from data
        # If fakes are not used, or one changes the code to avoid reading everything again for fakes, one could
        # actually modify directly the input histograms here, which might save some time compared to the copy.
        h = h.copy()
        if forceNonzero:
            h = hh.clipNegativeVals(h, createNew=False)
        if scaleToNewLumi > 0:
            h = hh.scaleHist(h, scaleToNewLumi, createNew=False)                        
        scale = self.processScaleFactor(proc)
        if scaleOp:
            scale = scale*scaleOp(proc)
        return h*scale

    def histName(self, baseName, procName="", syst=""):
        return Datagroups.histName(baseName, procName, syst, nominalName=self.nominalName)

    def histNameCombine(self, procName, baseName, syst, channel):
        return Datagroups.histNameCombine(procName, baseName, syst, channel)

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
