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
import numpy as np

from wremnants.datasets.datagroup import Datagroup
from wremnants.datasets.dataset_tools import getDatasets

logger = logging.child_logger(__name__)

class Datagroups(object):

    def __init__(self, infile, combine=False, filter_datasets=True, mode=None, **kwargs):
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
            raise ValueError(f"{infile} has unsupported file type")

        mode_map = {
            "w_z_gen_dists.py" : "vgen",
            "mz_dilepton.py" : "dilepton",
            "mz_wlike_with_mu_eta_pt.py" : "wlike",
            "mw_with_mu_eta_pt.py" : "wmass",
            "mw_lowPU.py" : "lowpu_w",
            "mz_lowPU.py" : "lowpu_z",
        }
        if mode == None:
            analysis_script = os.path.basename(self.getScriptCommand().split()[0])
            if analysis_script not in mode_map:
                raise ValueError(f"Unrecognized analysis script {analysis_script}! Expected one of {mode_map.keys()}")
            self.mode = mode_map[analysis_script]
        else:
            if mode not in mode_map.values():
                raise ValueError(f"Unrecognized mode '{mode}.' Must be one of {set(mode_map.values())}")
            self.mode = mode

        self.lumi = 1

        if self.results:
            try:
                args = self.getMetaInfo()["args"]
                self.flavor = args.get("flavor", None)
            except ValueError as e:
                logger.warning(e)
                self.flavor = None
        else:
            self.flavor = None

        self.groups = {}
        self.nominalName = "nominal"
        self.globalAction = None
        self.unconstrainedProcesses = []
        self.fakeName = "Fake"
        self.dataName = "Data"

        if "lowpu" in self.mode:
            from wremnants.datasets.datagroupsLowPU import make_datagroups_lowPU as make_datagroups
        else:
            from wremnants.datasets.datagroups2016 import make_datagroups_2016 as make_datagroups

        datasets = getDatasets(mode=self.mode)

        self.setDatasets(datasets, filter_datasets)
        self.setGenAxes()
                    
        make_datagroups(self, **kwargs)

    def __del__(self):
        if self.h5file:
            self.h5file.close()
        if self.rtfile:
            self.rtfile.Close()

    def setDatasets(self, datasets, filter_datasets=True):
        if self.results:
            # only keep datasets that are found in input file
            self.datasets = {x.name : x for x in datasets if not filter_datasets or x.name in self.results.keys() }
            
            # dictionary that maps dataset names to groups 
            dataset_to_group = {d_key: d.group for d_key, d in self.datasets.items()}

            for d_name, dataset in self.results.items():
                # if additional datasets are specified in results (for example aggregated groups or re-named datasets), get them
                if d_name in self.datasets.keys():
                    continue
                if d_name in ["meta_info",] or d_name.startswith("hist"):
                    continue
                
                g_name = d_name.replace("Bkg","") if d_name.startswith("Bkg") else d_name
                if g_name not in dataset_to_group.values():
                    g_name = dataset_to_group.get(g_name, g_name)
                
                logger.debug(f"Add dataset {d_name}")
                self.datasets[d_name] = narf.Dataset(**{
                    "name": d_name,
                    "group": g_name,
                    "filepaths": dataset["dataset"]["filepaths"],
                    "xsec": dataset["dataset"].get("xsec", None)
                    })

            self.data = [x for x in self.datasets.values() if x.is_data]
            if self.data:
                self.lumi = sum([self.results[x.name]["lumi"] for x in self.data if x.name in self.results])
                logger.info(f"Integrated luminosity from data: {self.lumi}/fb")
            else:
                logger.warning("No data process was selected, normalizing MC to 1/fb")

        else:
            self.datasets = {x.name : x for x in datasets}

        logger.debug(f"Getting these datasets: {self.datasets.keys()}")

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
        if self.globalAction is None:
            self.globalAction = action
        else:
            self.globalAction = lambda h, old_action=self.globalAction: action(old_action(h))

    def setNominalName(self, name):
        self.nominalName = name

    def processScaleFactor(self, proc):
        if proc.is_data or proc.xsec is None:
            return 1
        return self.lumi*1000*proc.xsec/self.results[proc.name]["weight_sum"]

    def getMetaInfo(self):
        if self.results:
            if "meta_info" not in self.results and "meta_data" not in self.results:
                raise ValueError("Did not find meta data in results file")
            return self.results["meta_info"] if "meta_info" in self.results else self.results["meta_data"]
        raise NotImplementedError("Currently can't access meta data as dict for ROOT file")

    def getScriptCommand(self):
        if self.rtfile:
            return self.rtfile.Get("meta_info/command").GetTitle()
        else:
            meta_info = self.getMetaInfo()
            return meta_info["command"]
        
    # for reading pickle files
    # as a reminder, the ND hists with tensor axes in the pickle files are organized as
    # pickle[procName]["output"][baseName] where
    ## procName are grouped into datagroups
    ## baseName takes values such as "nominal"
    def loadHistsForDatagroups(self, 
        baseName, syst, procsToRead=None, label=None, nominalIfMissing=True, 
        applySelection=True, fakerateIntegrationAxes=[], forceNonzero=True, preOpMap=None, preOpArgs=None, scaleToNewLumi=1, 
        excludeProcs=None, forceToNominal=[], sumFakesPartial=True
    ):
        logger.debug("Calling loadHistsForDatagroups()")
        logger.debug(f"The basename and syst is: {baseName}, {syst}")
        logger.debug(f"The procsToRead and excludedProcs are: {procsToRead}, {excludeProcs}")
        if not label:
            label = syst if syst else baseName
        # this line is annoying for the theory agnostic, too many processes for signal
        logger.debug(f"In loadHistsForDatagroups(): for hist {syst} procsToRead = {procsToRead}")

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
        if sumFakesPartial and self.fakeName in procsToRead:
            procsToReadSort = [x for x in procsToRead if x != self.fakeName] + [self.fakeName]
            hasFake = True
            fakesMembers = [m.name for m in self.groups[self.fakeName].members]
            fakesMembersWithSyst = []
            logger.debug(f"Has fake members: {fakesMembers}")
        else:
            hasFake = False
            procsToReadSort = [x for x in procsToRead]
        # Note: if 'hasFake' is kept as False (but Fake exists), the original behaviour for which Fake reads everything again is restored
        for procName in procsToReadSort:
            logger.debug(f"Reading group {procName}")

            if procName not in self.groups.keys():
                raise RuntimeError(f"Group {procName} not known. Defined groups are {list(self.groups.keys())}.")
            group = self.groups[procName]

            # Check if the histogram is already there and in case use it
            if Datagroups.histName(baseName, procName, syst) in group.hists:
                logger.debug(f"Existing histogram for proc {procName} base name {baseName} syst {syst} found.")
                group.hists[label] = group.hists[Datagroups.histName(baseName, procName, syst)]
                foundExact = True
                continue
            group.hists[label] = None

            for i, member in enumerate(group.members):
                if procName == self.fakeName and member.name in fakesMembersWithSyst:
                    # if we are here this process has been already used to build the fakes when running for other groups
                    continue
                logger.debug(f"Looking at group member {member.name}")
                read_syst = syst
                if member.name in forceToNominal:
                    read_syst = ""
                    logger.debug(f"Forcing group member {member.name} to read the nominal hist for syst {syst}")
                try:
                    h = self.readHist(baseName, member, procName, read_syst)
                    foundExact = True
                except ValueError as e:
                    if nominalIfMissing:
                        logger.info(f"{str(e)}. Using nominal hist {self.nominalName} instead")
                        h = self.readHist(self.nominalName, member, procName, "")
                    else:
                        logger.warning(str(e))
                        continue

                h_id = id(h)

                logger.debug(f"Hist axes are {h.axes.name}")

                if group.memberOp:
                    if group.memberOp[i] is not None:
                        logger.debug(f"Apply operation to member {i}: {member.name}/{procName}")
                        h = group.memberOp[i](h)
                    else:
                        logger.debug(f"No operation for member {i}: {member.name}/{procName}")

                if preOpMap and member.name in preOpMap:
                    logger.debug(f"Applying preOp to {member.name}/{procName} after loading")
                    h = preOpMap[member.name](h, **preOpArgs)

                if self.gen_axes != None:
                    # integrate over remaining gen axes 
                    logger.debug(f"Integrate over gen axes {self.gen_axes}")
                    projections = [a for a in h.axes.name if a not in self.gen_axes]
                    if len(projections) < len(h.axes.name):
                        h = h.project(*projections)
                    logger.debug(f"Integrated, Hist axes are {h.axes.name}")

                if h_id == id(h):
                    logger.debug(f"Make explicit copy")
                    h = h.copy()

                if self.globalAction:
                    logger.debug("Applying global action")
                    h = self.globalAction(h)

                if forceNonzero:
                    logger.debug("force non zero")
                    h = hh.clipNegativeVals(h, createNew=False)

                scale = self.processScaleFactor(member)
                scale *= scaleToNewLumi
                if group.scale:
                    scale *= group.scale(member)

                if not np.isclose(scale, 1, rtol=0, atol=1e-10):
                    logger.debug(f"Scale hist with {scale}")
                    h = hh.scaleHist(h, scale, createNew=False)

                hasPartialSumForFake = False
                if hasFake and procName != self.fakeName:
                    if member.name in fakesMembers:
                        logger.debug("Make partial sums for fakes")
                        if member.name not in fakesMembersWithSyst:
                            fakesMembersWithSyst.append(member.name)
                        hasPartialSumForFake = True
                        # apply the correct scale for fakes
                        scaleProcForFake = self.groups[self.fakeName].scale(member)
                        logger.debug(f"Summing hist {read_syst} for {member.name} to {self.fakeName} with scale = {scaleProcForFake}")
                        hProcForFake = scaleProcForFake * h
                        histForFake = hh.addHists(histForFake, hProcForFake, createNew=False) if histForFake else hProcForFake
                                
                # The following must be done when the group is not Fake, or when the previous part for fakes was not done
                # For fake this essentially happens when the process doesn't have the syst, so that the nominal is used
                if procName != self.fakeName or (procName == self.fakeName and not hasPartialSumForFake):
                    if procName == self.fakeName:
                        logger.debug(f"Summing nominal hist instead of {syst} to {self.fakeName} for {member.name}")
                    else:
                        logger.debug(f"Summing {read_syst} to {procName} for {member.name}")

                    group.hists[label] = hh.addHists(group.hists[label], h, createNew=False) if group.hists[label] else h

            if not nominalIfMissing and group.hists[label] is None:
                continue

            # now sum to fakes the partial sums which where not already done before
            # (group.hists[label] contains only the contribution from nominal histograms).
            # Then continue with the rest of the code as usual
            if hasFake and procName == self.fakeName:
                if histForFake is not None:
                    group.hists[label] = hh.addHists(group.hists[label], histForFake, createNew=False) if group.hists[label] else histForFake

            # Can use to apply common rebinning or selection on top of the usual one
            if group.rebinOp:
                logger.debug(f"Apply rebin operation for process {procName}")
                group.hists[label] = group.rebinOp(group.hists[label])

            if group.selectOp:
                if not applySelection:
                    logger.warning(f"Selection requested for process {procName} but applySelection=False, thus it will be ignored")
                elif label in group.hists.keys():
                    logger.debug(f"Apply selection for process {procName}")
                    if procName == self.fakeName and "fakerate_integration_axes" not in group.selectOpArgs and len(fakerateIntegrationAxes):
                        opArgs = {**group.selectOpArgs, "fakerate_integration_axes": fakerateIntegrationAxes}
                    else:
                        opArgs = group.selectOpArgs
                    if group.hists[label]:
                        group.hists[label] = group.selectOp(group.hists[label], **opArgs)

        # Avoid situation where the nominal is read for all processes for this syst
        if nominalIfMissing and not foundExact:
            raise ValueError(f"Did not find systematic {syst} for any processes!")

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
            if group_name not in self.groups:
                raise ValueError(f"Trying to expand unknown group {group_name}. Valid groups are {list(self.groups.keys())}")
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

    def addSummedProc(self, refname, name, label=None, color=None, exclude=["Data"], relabel=None, 
            procsToRead=None, reload=False, rename=None, action=None, preOpMap={}, preOpArgs={}, 
            fakerateIntegrationAxes=[], forceNonzero=True):
        if reload:
            self.loadHistsForDatagroups(refname, syst=name, excludeProcs=exclude,
                procsToRead=procsToRead, preOpMap=preOpMap, preOpArgs=preOpArgs, 
                fakerateIntegrationAxes=fakerateIntegrationAxes, forceNonzero=forceNonzero)

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

    def setGenAxes(self, gen_axes=None):
        if isinstance(gen_axes, str):
            gen_axes = [gen_axes]

        self.gen_axes = None

        if gen_axes != None:
            self.gen_axes = list(gen_axes)
        else:
            # infer gen axes from metadata
            try:
                args = self.getMetaInfo()["args"]
            except ValueError as e:
                logger.warning("No meta data found so no gen axes could be auto set")
                return
            if args.get("unfolding", False) is False and args.get("theoryAgnostic", False) is False:
                return

            if len(args.get("genVars", [])) > 0:
                self.gen_axes = args["genVars"]
            else:
                logger.warning(f"Unknown gen axes!")

        logger.debug(f"Gen axes are now {self.gen_axes}")

    def getGenBinIndices(self, h, axesToRead=[]):
        gen_bins = []
        for gen_axis in (axesToRead if len(axesToRead) else self.gen_axes):
            if gen_axis not in h.axes.name:
                raise RuntimeError(f"Gen axis '{gen_axis}' not found in histogram axes '{h.axes.name}'!")

            gen_bin_list = [i for i in range(h.axes[gen_axis].size)]
            if h.axes[gen_axis].traits.underflow:
                gen_bin_list.append(hist.underflow)
            if h.axes[gen_axis].traits.overflow:
                gen_bin_list.append(hist.overflow)
            gen_bins.append(gen_bin_list)
        return gen_bins

    def defineSignalBinsUnfolding(self, group_name, new_name=None, member_filter=None, histToReadAxes="xnorm", axesToRead=[]):
        if group_name not in self.groups.keys():
            raise RuntimeError(f"Base group {group_name} not found in groups {self.groups.keys()}!")

        base_members = self.groups[group_name].members[:]
        if member_filter is not None:
            base_members = [m for m in filter(lambda x, f=member_filter: f(x), base_members)]            

        if histToReadAxes not in self.results[base_members[0].name]["output"]:
            raise ValueError(f"Results for member {base_members[0].name} does not include xnorm. Found {self.results[base_members[0].name]['output'].keys()}")
        nominal_hist = self.results[base_members[0].name]["output"][histToReadAxes].get()

        gen_bin_indices = self.getGenBinIndices(nominal_hist, axesToRead=axesToRead)

        for indices in itertools.product(*gen_bin_indices):

            proc_name = group_name if new_name is None else new_name
            for idx, var in zip(indices, self.gen_axes):
                if idx == hist.underflow:
                    idx_str = "U"
                elif idx == hist.overflow:
                    idx_str = "O"
                else:
                    idx_str = str(idx)
                proc_name += f"_{var}{idx_str}"

            self.copyGroup(group_name, proc_name, member_filter=member_filter)

            memberOp = lambda x, indices=indices, genvars=self.gen_axes: x[{var : i for var, i in zip(genvars, indices)}]
            self.groups[proc_name].memberOp = [memberOp for m in base_members]

            self.unconstrainedProcesses.append(proc_name)

    def select_xnorm_groups(self, select_groups=None):
        # only keep members and groups where xnorm is defined
        logger.info("Select xnorm groups"+(f" {select_groups}" if select_groups else ""))
        if select_groups is not None:
            if isinstance(select_groups, str):
                select_groups = [select_groups]
            self.deleteGroups([g for g in self.groups.keys() if g not in select_groups])
        elif self.fakeName in self.groups:
            self.deleteGroup(self.fakeName)
        toDel_groups = []
        for g_name, group in self.groups.items():
            toDel_members = []
            for member in group.members:
                if member.name not in self.results.keys():
                    raise RuntimeError(f"The member {member.name} of group {g_name} was not found in the results!")
                if "xnorm" not in self.results[member.name]["output"].keys():
                    logger.debug(f"Member {member.name} has no xnorm and will be deleted")
                    toDel_members.append(member)
            if len(toDel_members) == len(group.members):
                logger.warning(f"All members of group {g_name} have no xnorm and the group will be deleted")
                toDel_groups.append(g_name)
            else:
                group.deleteMembers(toDel_members)
        self.deleteGroups(toDel_groups)

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

    def set_rebin_action(self, axes, ax_lim=[], ax_rebin=[], ax_absval=[]):
        if len(ax_lim) % 2 or len(ax_lim)/2 > len(axes) or len(ax_rebin) > len(axes):
            raise ValueError("Inconsistent rebin or axlim arguments. axlim must be at most two entries per axis, and rebin at most one")

        sel = {}
        for var,low,high,rebin in itertools.zip_longest(axes, ax_lim[::2], ax_lim[1::2], ax_rebin):
            s = hist.tag.Slicer()
            if low is not None and high is not None:
                logger.info(f"Restricting the axis '{var}' to range [{low}, {high}]")
                sel[var] = s[complex(0, low):complex(0, high):hist.rebin(rebin) if rebin else None]
            elif rebin:
                sel[var] = s[hist.rebin(rebin)]
            if rebin:
                logger.info(f"Rebinning the axis '{var}' by [{rebin}]")

        if len(sel) > 0:
            logger.info(f"Will apply the global selection {sel}")
            self.setGlobalAction(lambda h: h[sel])

        for i, (var, absval) in enumerate(itertools.zip_longest(axes, ax_absval)):
            if absval:
                logger.info(f"Taking the absolute value of axis '{var}'")
                self.setGlobalAction(lambda h, ax=var: hh.makeAbsHist(h, ax))
                axes[i] = f"abs{var}"

    def readHist(self, baseName, proc, group, syst):
        output = self.results[proc.name]["output"]
        histname = self.histName(baseName, proc.name, syst)
        logger.debug(f"Reading hist {histname} for proc/group {proc.name}/{group} and syst '{syst}'")
        if histname not in output:
            raise ValueError(f"Histogram {histname} not found for process {proc.name}")

        h = output[histname]
        if isinstance(h, narf.ioutils.H5PickleProxy):
            h = h.get()

        return h

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

