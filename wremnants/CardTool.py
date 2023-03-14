from collections import OrderedDict
from utilities import output_tools,boostHistHelpers as hh,common
import narf
import logging
import ROOT
import uproot
import time
import numpy as np
import collections.abc
import os
import itertools
import re
import hist

logger = common.child_logger(__name__)

def notImplemented(operation="Unknown"):
    raise NotImplementedError(f"Required operation '{operation}' is not implemented!")

class CardTool(object):
    def __init__(self, cardName="card.txt"):
        self.skipHist = False # don't produce/write histograms, file with them already exists
        self.outfile = None
        self.cardName = cardName
        self.systematics = {}
        self.lnNSystematics = {}
        self.procDict = {}
        self.predictedProcs = []
        self.fakeEstimate = None
        self.channels = ["plus", "minus"]
        self.cardContent = {}
        self.cardGroups = {}
        self.nominalTemplate = ""
        self.spacing = 28
        self.systTypeSpacing = 12
        self.procColumnsSpacing = 12
        self.fakeName = "Fake" # but better to set it explicitly
        self.dataName = "Data"
        self.nominalName = "nominal"
        self.datagroups = None
        self.pseudodata_datagroups = None
        self.unconstrainedProcesses = None
        self.noStatUncProcesses = []
        self.buildHistNameFunc = None
        self.histName = "x"
        self.nominalDim = None
        self.pseudoData = None
        self.pseudoDataIdx = None
        self.excludeSyst = None
        self.excludeProcGroups = None # exclude processes, by group name
        self.writeByCharge = True
        self.keepSyst = None # to override previous one with exceptions for special cases
        #self.loadArgs = {"operation" : "self.loadProcesses (reading hists from file)"}
        self.lumiScale = 1.
        self.project = None
        self.keepOtherChargeSyst = True
        self.chargeIdDict = {"minus" : {"val" : -1, "id" : "q0", "badId" : None},
                             "plus"  : {"val" : 1., "id" : "q1", "badId" : None},
                             "inclusive" : {"val" : "sum", "id" : "none", "badId" : "NOBADID"}, # for this channel there is no bad id, currently using random string to make sure it doesn't match
                             }

    def skipHistograms(self):
        self.skipHist = True
        if len(self.noStatUncProcesses):
            logger.info("Attention: histograms are not saved according to input options, thus statistical uncertainty won't be zeroed")
        
    def setSkipOtherChargeSyst(self):
        self.keepOtherChargeSyst = False
        self.chargeIdDict["plus"]["badId"] = "q0"
        self.chargeIdDict["minus"]["badId"] = "q1"

    def setProjectionAxes(self, project):
        self.project = project

    def setProcsNoStatUnc(self, procs, resetList=True):
        if self.skipHist:
            logger.info("Attention: trying to set statistical uncertainty to 0 for some processes, but histograms won't be saved according to input options")
        if resetList:
            self.noStatUncProcesses = []
        if isinstance(procs, str):
            self.noStatUncProcesses.append(procs)
        elif isinstance(procs, list):
            self.noStatUncProcesses.extend(procs)
        else:
            raise ValueError("In setNoStatUncForProcs(): expecting string or list argument")
    
    def setLumiScale(self, lumiScale):
        self.lumiScale = lumiScale
        
    def getProcsNoStatUnc(self):
        return self.noStatUncProcesses
        
    ## Functions to customize systs to be added in card, mainly for tests
    def setCustomSystForCard(self, exclude=None, keep=None):
        if exclude: self.excludeSyst = re.compile(exclude)
        if keep:    self.keepSyst    = re.compile(keep)
        
    def isExcludedNuisance(self, name):
        # note, re.match search for a match from the beginning, so if x="test" x.match("mytestPDF1") will NOT match 
        # might use re.search instead to be able to match from anywhere inside the name
        if self.excludeSyst != None and self.excludeSyst.match(name):
            if self.keepSyst != None and self.keepSyst.match(name):
                return False
            else:
                logger.info(f"   Excluding nuisance: {name}")
                return True
        else:
            return False
        
    # Function call to load hists for processes (e.g., read from a ROOT file)
    # Extra args will be passed to each call
    def setLoadDatagroups(self, datagroups, extraArgs={}):
        self.datagroups = datagroups
        self.loadArgs = extraArgs
    
    def setFakeName(self, name):
        self.fakeName = name

    def getFakeName(self):
        return self.fakeName

    def setPseudodata(self, pseudodata, idx = 0):
        self.pseudoData = pseudodata
        self.pseudoDataIdx = idx

    # Needs to be increased from default for long proc names
    def setSpacing(self, spacing):
        self.spacing = spacing

    def setDataName(self, name):
        self.dataName = name

    def setDatagroups(self, datagroups, resetGroups=False):
        self.datagroups = datagroups 
        if self.nominalName:
            self.datagroups.setNominalName(self.nominalName)
        # if processes are not set yet, do it now to skip one step
        # FIXME: should it allow one to reset the procDict by passing a flag, e.g. in case one loads a new datagroups?
        if resetGroups or not self.procDict:
            self.setProcesses(self.datagroups.getNames(afterFilter=True))
        
    def setPseudodataDatagroups(self, datagroups):
        self.pseudodata_datagroups = datagroups 
        if self.nominalName:
            self.pseudodata_datagroups.setNominalName(self.nominalName)
        
    def setChannels(self, channels):
        self.channels = channels
        
    def setWriteByCharge(self, writeByCharge):
        self.writeByCharge = writeByCharge

    def setNominalTemplate(self, template):
        if not os.path.isfile(template):
            raise IOError(f"Template file {template} is not a valid file")
        self.nominalTemplate = template

    def predictedProcesses(self):
        if self.predictedProcs:
            return self.predictedProcs
        return list(filter(lambda x: x != self.dataName, self.procDict.keys()))

    def setHistName(self, histName):
        self.histName = histName

    def setNominalName(self, histName):
        self.nominalName = histName
        if self.datagroups:
            self.datagroups.setNominalName(histName)

    # by default this returns True also for Fake since it has Data in the list of members
    # then self.isMC negates this one and thus will only include pure MC processes
    def isData(self, procName, onlyData=False):
        if onlyData:
            return all([x.is_data for x in self.datagroups.groups[procName]["members"]])
        else:
            return any([x.is_data for x in self.datagroups.groups[procName]["members"]])

    def isMC(self, procName):
        return not self.isData(procName)

    def addFakeEstimate(self, estimate):
        self.fakeEstimate = estimate

    def setProcesses(self, processes):
        self.procDict = {proc: {} for proc in processes}

    def getProcesses(self):
        return list(self.procDict.keys())

    def setExcludedProcs(self, proc_list):
        if isinstance(proc_list, list):
            self.excludeProcGroups = [x for x in proc_list]
        elif isinstance(proc_list, str):
            self.excludeProcGroups = [proc_list]
        elif proc_list is not None:
            logger.error(f"In CardTool.setExcludedProcs(): invalid argument with type {type(proc_list)}")
            quit()
            
    def filteredProcesses(self, filterExpr):
        return list(filter(filterExpr, self.procDict.keys()))

    def allMCProcesses(self):
        return self.filteredProcesses(lambda x: self.isMC(x))

    def addLnNSystematic(self, name, size, processes, group=None, groupFilter=None, decorrCharge=False):
        if not self.isExcludedNuisance(name):
            self.lnNSystematics.update({name : {"size" : size,
                                                "processes" : processes,
                                                "group" : group,
                                                "groupFilter" : groupFilter,
                                                "decorrCharge" : decorrCharge}
            })

    # action will be applied to the sum of all the individual samples contributing, actionMap should be used
    # to apply a separate action per process. this is needed for example for the scale uncertainty split
    # by pt or helicity
    # action takes place after mirroring
    # use doActionBeforeMirror to do something before it instead (so the mirroring will act on the modified histogram)
    # decorrelateByBin is to customize eta-pt decorrelation: pass dictionary with {axisName: [bin edges]}
    def addSystematic(self, name, systAxes, outNames=None, skipEntries=None, labelsByAxis=None, 
                      baseName="", mirror=False, scale=1, processes=None, group=None, noConstraint=False,
                      action=None, doActionBeforeMirror=False, actionArgs={}, actionMap={},
                      systNameReplace=[], groupFilter=None, passToFakes=False,
                      rename=None, splitGroup={}, decorrelateByCharge=False, decorrelateByBin={}):

        # Need to make an explicit copy of the array before appending
        procs_to_add = [x for x in (self.allMCProcesses() if processes is None else processes)]
        if passToFakes and self.getFakeName() not in procs_to_add and self.getFakeName() not in self.excludeProcGroups:
            procs_to_add.append(self.getFakeName())

        if action and actionMap:
            raise ValueError("Only one of action and actionMap args are allowed")

        # protection when the input list is empty because of filters but the systematic is built reading the nominal
        # since the nominal reads all filtered processes regardless whether a systematic is passed to them or not
        # this can happen when creating new systs by scaling of the nominal histogram
        if not len(procs_to_add):
            return 0

        if name == self.nominalName:
            logger.debug("Defining syst {rename} from nominal histogram")
            
        self.systematics.update({
            name if not rename else rename : {
                "outNames" : [] if not outNames else outNames,
                "baseName" : baseName,
                "processes" : procs_to_add,
                "systAxes" : systAxes,
                "labelsByAxis" : systAxes if not labelsByAxis else labelsByAxis,
                "group" : group,
                "groupFilter" : groupFilter,
                "splitGroup" : splitGroup if len(splitGroup) else {group : ".*"}, # dummy dictionary if splitGroup=None, to allow for uniform treatment
                "scale" : scale,
                "mirror" : mirror,
                "action" : action,
                "doActionBeforeMirror" : doActionBeforeMirror,
                "actionMap" : actionMap,
                "actionArgs" : actionArgs,
                "systNameReplace" : systNameReplace,
                "noConstraint" : noConstraint,
                "skipEntries" : [] if not skipEntries else skipEntries,
                "name" : name,
                "decorrCharge" : decorrelateByCharge,
                "decorrByBin": decorrelateByBin
            }
        })
        
    def setMirrorForSyst(self, syst, mirror=True):
        self.systematics[syst]["mirror"] = mirror

    def systLabelForAxis(self, axLabel, entry, axis):
        if type(axis) == hist.axis.StrCategory:
            if entry in axis:
                return entry
            else:
                raise ValueError(f"Did not find label {entry} in categorical axis {axis}")
        if axLabel == "mirror":
            return 'Down' if entry else 'Up' # first entry is the original, call it Up since it is usually defined by an actual scaling up of something (e.g. efficiencies)
        if axLabel == "downUpVar":
            return 'Up' if entry else 'Down'
        if "{i}" in axLabel:
            return axLabel.format(i=entry)
        return axLabel+str(entry)

    def excludeSystEntry(self, entry, skipEntries):
        for skipEntry in skipEntries:
            skip = False
            for e,match in zip(entry, skipEntry): 
                # Can use -1 to exclude all values of an axis
                if match == -1 or match == e or re.match(str(match), str(e)):
                    return True
        return False

    def expandSkipEntries(self, h, syst, skipEntries):
        updated_skip = []
        for skipEntry in skipEntries:
            nsyst_ax = len(self.systematics[syst]["systAxes"])+self.systematics[syst]["mirror"]
            if len(skipEntry) != nsyst_ax:
                raise ValueError(f"Error in syst {syst}. skipEntry must specify a value per axis. "
                        f"Found {nsyst_ax} axes ({self.systematics[syst]['systAxes']}) but {len(skipEntry)} "
                        "entries were given")
            # The lookup is handled by passing an imaginary number,
            # so detect these and then call the bin lookup on them
            # np.iscomplex returns false for 0.j, but still want to detect that
            to_lookup = np.array([isinstance(x, complex) for x in skipEntry])
            skip_arr = np.array(skipEntry)
            if to_lookup.any():
                bin_lookup = np.array([ax.index(x.imag) for x, ax in 
                    zip(skipEntry, h.axes[-nsyst_ax:]) if isinstance(x, complex)])
                skip_arr = skip_arr.real
                skip_arr[to_lookup] += bin_lookup
            updated_skip.append(skip_arr)

        return updated_skip

    def systHists(self, hvar, syst):
        if syst == self.nominalName:
            return {self.nominalName : hvar}

        systInfo = self.systematics[syst] 
        systAxes = systInfo["systAxes"]
        systAxesLabels = systAxes
        if "labelsByAxis" in systInfo:
            systAxesLabels = systInfo["labelsByAxis"]

        # Jan: moved above the mirror action, as this action can cause mirroring
        if systInfo["action"] and not systInfo["doActionBeforeMirror"]:
            hvar = systInfo["action"](hvar, **systInfo["actionArgs"])
        self.outfile.cd() # needed to restore the current directory in case the action opens a new root file
            
        axNames = systAxes[:]
        axLabels = systAxesLabels[:]
        if hvar.axes[-1].name == "mirror":
            axNames.append("mirror")
            axLabels.append("mirror")
        axes = [hvar.axes[ax] for ax in axNames]

        if not all([name in hvar.axes.name for name in axNames]):
            raise ValueError(f"Failed to find axis names {str(axNames)} in hist for syst {syst}. " \
                f"Axes in hist are {str(hvar.axes.name)}")

        # Converting to a list becasue otherwise if you print it for debugging you loose it
        entries = list(itertools.product(*[[x for x in ax] if type(ax) == hist.axis.StrCategory else range(ax.size) for ax in axes]))
        
        if len(systInfo["outNames"]) == 0:
            for entry in entries:
                skipEntries = None if "skipEntries" not in systInfo else self.expandSkipEntries(hvar, syst, systInfo["skipEntries"])
                if skipEntries and self.excludeSystEntry(entry, skipEntries):
                    systInfo["outNames"].append("")
                else:
                    name = systInfo["baseName"]
                    name += "".join([self.systLabelForAxis(al, entry[i], ax) for i,(al,ax) in enumerate(zip(axLabels,axes))])
                    if "systNameReplace" in systInfo and systInfo["systNameReplace"]:
                        for rep in systInfo["systNameReplace"]:
                            name = name.replace(*rep)
                    # Obviously there is a nicer way to do this...
                    if "Up" in name:
                        name = name.replace("Up", "")+"Up"
                    elif "Down" in name:
                        name = name.replace("Down", "")+"Down"
                    systInfo["outNames"].append(name)
            if not len(systInfo["outNames"]):
                raise RuntimeError(f"Did not find any valid variations for syst {syst}")

        variations = [hvar[{ax : binnum for ax,binnum in zip(axNames, entry)}] for entry in entries]
        if len(variations) != len(systInfo["outNames"]):
            logger.warning(f"The number of variations doesn't match the number of names for "
                f"syst {syst}. Found {len(systInfo['outNames'])} names and {len(variations)} variations.")

        return {name : var for name,var in zip(systInfo["outNames"], variations) if name}

    def variationName(self, proc, name):
        if name == self.nominalName:
            return f"{self.histName}_{proc}"
        else:
            return f"{self.histName}_{proc}_{name}"

    def checkSysts(self, hnom, var_map, proc, thresh=0.25):
        #if self.check_variations:
        var_names = set([name.replace("Up", "").replace("Down", "") for name in var_map.keys() if name])
        if len(var_names) != len(var_map.keys())/2:
            raise ValueError(f"Invalid syst names for process {proc}! Expected an up/down variation for each syst. "
                f"Found systs {var_names} and outNames {var_map.keys()}")
        for name in var_names:
            up = var_map[name+"Up"]
            down = var_map[name+"Down"]
            try:
                up_relsign = np.sign(up.values()-hnom.values())
            except ValueError as e:
                logger.error(f"Incompatible shapes between up and down for syst {name}")
                raise e
            down_relsign = np.sign(down.values()-hnom.values())
            vars_sameside = (up_relsign != 0) & (up_relsign == down_relsign)
            perc_sameside = np.count_nonzero(vars_sameside)/hnom.size 
            if perc_sameside > thresh:
                logger.warning(f"{perc_sameside:.0%} bins are one sided for syst {name} and process {proc}!")

    def writeForProcess(self, h, proc, syst):
        decorrelateByCharge = False
        decorrelateByBin = {}
        if syst != self.nominalName:
            systInfo = self.systematics[syst]
            if systInfo["doActionBeforeMirror"] and systInfo["action"]:
                h =systInfo["action"](h, **systInfo["actionArgs"])
                self.outfile.cd() # needed to restore the current directory in case the action opens a new root file
            if systInfo["mirror"]:
                hnom = self.procDict[proc][self.nominalName]
                h = hh.extendHistByMirror(h, hnom)
            if systInfo["decorrCharge"]:
                decorrelateByCharge = True
            if systInfo["decorrByBin"]:
                decorrelateByBin = systInfo["decorrByBin"]
        # Otherwise this is a processes not affected by the variation, don't write it out,
        # it's only needed for the fake subtraction
        logger.info(f"Writing systematic {syst} for process {proc}")
        var_map = self.systHists(h, syst) 
        # TODO: Make this optional
        if syst != self.nominalName:
            self.checkSysts(self.procDict[proc][self.nominalName], var_map, proc)
        setZeroStatUnc = False
        if proc in self.noStatUncProcesses:
            logger.info(f"Zeroing statistical uncertainty for process {proc}")
            setZeroStatUnc = True
        for name, var in var_map.items():
            if name != "":
                self.writeHist(var, self.variationName(proc, name), setZeroStatUnc=setZeroStatUnc,
                               decorrCharge=decorrelateByCharge, decorrByBin=decorrelateByBin)

    def addPseudodata(self, processes):
        datagroups = self.datagroups if not self.pseudodata_datagroups else self.pseudodata_datagroups
        datagroups.loadHistsForDatagroups(
            baseName=self.pseudoData, syst="", label=self.pseudoData,
            procsToRead=processes, excluded_procs=self.excludeProcGroups,
            scaleToNewLumi=self.lumiScale)
        procDict = datagroups.getDatagroups()
        logger.warning(f"Making pseudodata summing these proceses: {processes}")
        hists = [procDict[proc][self.pseudoData] for proc in processes]
        hdata = hh.sumHists(hists)
        # Kind of hacky, but in case the alt hist has uncertainties
        for systAxName in ["systIdx", "tensor_axis_0", "vars"]:
            if systAxName in [ax.name for ax in hdata.axes]:
                hdata = hdata[{systAxName : self.pseudoDataIdx }] 
        self.writeHist(hdata, self.pseudoData+"_sum")

    def writeForProcesses(self, syst, processes, label):
        for process in processes:
            hvar = self.procDict[process][label]
            if not hvar:
                raise RuntimeError(f"Failed to load hist for process {process}, systematic {syst}")
            self.writeForProcess(hvar, process, syst)
        if syst != self.nominalName:
            self.fillCardWithSyst(syst)

    def setOutfile(self, outfile):
        if type(outfile) == str:
            if self.skipHist:
                self.outfile = outfile # only store name, file will not be used and doesn't need to be opened
            else:
                self.outfile = ROOT.TFile(outfile, "recreate")
                self.outfile.cd()
        else:
            self.outfile = outfile
            self.outfile.cd()
            
    def writeOutput(self, args=None):
        self.datagroups.loadHistsForDatagroups(
            baseName=self.nominalName, syst=self.nominalName,
            procsToRead=self.procDict.keys(), excluded_procs=self.excludeProcGroups,
            label=self.nominalName, 
            scaleToNewLumi=self.lumiScale)
        self.procDict = self.datagroups.getDatagroups()
        self.writeForProcesses(self.nominalName, processes=self.procDict.keys(), label=self.nominalName)
        self.loadNominalCard()
        if self.pseudoData:
            self.addPseudodata([x for x in self.procDict.keys() if x != "Data"])

        self.writeLnNSystematics()
        for syst in self.systematics.keys():
            if self.isExcludedNuisance(syst): continue
            systMap = self.systematics[syst]
            systName = syst if not systMap["name"] else systMap["name"]
            processes = systMap["processes"]
            self.datagroups.loadHistsForDatagroups(
                self.nominalName, systName, label="syst",
                procsToRead=processes, excluded_procs=self.excludeProcGroups,
                forceNonzero=systName != "qcdScaleByHelicity",
                preOpMap=systMap["actionMap"], preOpArgs=systMap["actionArgs"],
                scaleToNewLumi=self.lumiScale)
            self.writeForProcesses(syst, label="syst", processes=processes)
            
        output_tools.writeMetaInfoToRootFile(self.outfile, exclude_diff='notebooks', args=args)
        if self.skipHist:
            logger.info("Histograms will not be written because 'skipHist' flag is set to True")
        self.writeCard()

        
    def writeCard(self):
        for chan in self.channels:
            with open(self.cardName.format(chan=chan), "w") as card:
                card.write(self.cardContent[chan])
                card.write("\n")
                card.write(self.cardGroups[chan])

    def addSystToGroup(self, groupName, chan, members, groupLabel="group"):
        group_expr = f"{groupName} {groupLabel} ="
        if group_expr in self.cardGroups[chan]:
            idx = self.cardGroups[chan].index(group_expr)+len(group_expr)
            self.cardGroups[chan] = self.cardGroups[chan][:idx] + " " + members + self.cardGroups[chan][idx:]
        else:
            self.cardGroups[chan] += f"\n{group_expr} {members}"                                              

    def writeLnNSystematics(self):
        nondata = self.predictedProcesses()
        # exit this function when a syst is applied to no process (can happen when some are excluded)
        for name,info in self.lnNSystematics.items():
            if all(x not in info["processes"] for x in nondata):
                logger.warning(f"Skipping syst {name}, procs to apply it to would be {info['processes']}, and predicted processes are {nondata}")
                return 0
            include = [(str(info["size"]) if x in info["processes"] else "-").ljust(self.procColumnsSpacing) for x in nondata]
            group = info["group"]
            groupFilter = info["groupFilter"]
            for chan in self.channels:
                if info["decorrCharge"]:
                    newname = f"{name}_{self.chargeIdDict[chan]['id']}"
                else:
                    newname = name
                if self.keepOtherChargeSyst or self.chargeIdDict[chan]["badId"] not in name:
                    self.cardContent[chan] += f'{name.ljust(self.spacing)} lnN{" "*(self.systTypeSpacing-2)} {"".join(include)}\n'
                    if group and not self.isExcludedNuisance(name) and len(list(filter(groupFilter, [name]))):
                        self.addSystToGroup(group, chan, name)

    def fillCardWithSyst(self, syst):
        systInfo = self.systematics[syst]
        scale = systInfo["scale"]
        procs = systInfo["processes"]
        group = systInfo["group"]
        groupFilter = systInfo["groupFilter"]
        label = "group" if not systInfo["noConstraint"] else "noiGroup"
        nondata = self.predictedProcesses()
        names = [x[:-2] if "Up" in x[-2:] else (x[:-4] if "Down" in x[-4:] else x) 
                    for x in filter(lambda x: x != "", systInfo["outNames"])]
        # exit this function when a syst is applied to no process (can happen when some are excluded)
        if all(x not in procs for x in nondata):
            return 0
        
        if type(scale) != dict:
            include = [(str(scale) if x in procs else "-").ljust(self.procColumnsSpacing) for x in nondata]

        splitGroupDict = systInfo["splitGroup"]
        shape = "shape" if not systInfo["noConstraint"] else "shapeNoConstraint"

        # Deduplicate while keeping order
        systNames = list(dict.fromkeys(names))
        if systInfo["decorrCharge"]:
            systNames = list([f"{x}_q{i}" for x in systNames for i in range(2)])

        systnamesPruned = [s for s in systNames if not self.isExcludedNuisance(s)]
        systNames = systnamesPruned[:]
        for systname in systNames:
            if type(scale) == dict:
                for reg in scale.keys():
                    if re.match(reg, systname):
                        thiscale = str(scale[reg])
                        include = [(thiscale if x in procs else "-").ljust(self.procColumnsSpacing) for x in nondata]
                        break # exit this inner loop when match is found, to save time
            shape = "shape" if not systInfo["noConstraint"] else "shapeNoConstraint"
            for chan in self.channels:
                # do not write systs which should only apply to other charge, to simplify card
                if self.keepOtherChargeSyst or self.chargeIdDict[chan]["badId"] not in systname:
                    self.cardContent[chan] += f"{systname.ljust(self.spacing)} {shape.ljust(self.systTypeSpacing)} {''.join(include)}\n"
        # unlike for LnN systs, here it is simpler to act on the list of these systs to form groups, rather than doing it syst by syst 
        if group:
            for chan in self.channels:
                if self.keepOtherChargeSyst:
                    systNamesForGroupPruned = systNames[:]
                else:
                    systNamesForGroupPruned = [s for s in systNames if self.chargeIdDict[chan]["badId"] not in s]
                systNamesForGroup = list(systNamesForGroupPruned if not groupFilter else filter(groupFilter, systNamesForGroupPruned))
                if len(systNamesForGroup):
                    for subgroup in splitGroupDict.keys():
                        matchre = re.compile(splitGroupDict[subgroup])
                        systNamesForSubgroup = list(filter(lambda x: matchre.match(x),systNamesForGroup))
                        if len(systNamesForSubgroup):
                            members = " ".join(systNamesForSubgroup)
                            self.addSystToGroup(subgroup, chan, members, groupLabel=label)

    def setUnconstrainedProcs(self, procs):
        self.unconstrainedProcesses = procs

    def processLabels(self):
        nondata = np.array(self.predictedProcesses())
        labels = np.arange(len(nondata))+1
        issig = np.isin(nondata, self.unconstrainedProcesses)
        labels[issig] = -np.arange(np.count_nonzero(issig))-1
        return labels

    def loadNominalCard(self):
        procs = self.predictedProcesses()
        nprocs = len(procs)
        for chan in self.channels:
            args = {
                "channel" :  chan,
                "channelPerProc" : chan.ljust(self.procColumnsSpacing)*nprocs,
                "processes" : "".join([x.ljust(self.procColumnsSpacing) for x in procs]),
                "labels" : "".join([str(x).ljust(self.procColumnsSpacing) for x in self.processLabels()]),
                # Could write out the proper normalizations pretty easily
                "rates" : "-1".ljust(self.procColumnsSpacing)*nprocs,
                "inputfile" : self.outfile if type(self.outfile) == str  else self.outfile.GetName(),
                "dataName" : self.dataName,
                "histName" : self.histName,
                "pseudodataHist" : self.pseudoData+"_sum" if self.pseudoData else f"{self.histName}_{self.dataName}"
            }
            self.cardContent[chan] = output_tools.readTemplate(self.nominalTemplate, args)
            self.cardGroups[chan] = ""
            
    def writeHistByCharge(self, h, name, decorrCharge=False):
        for charge in self.channels:
            newname = name
            if decorrCharge:
                # name should always have Up/Down here, since this acts on syst histograms, but let's allow for all cases
                upDown = "Up" if name.endswith("Up") else "Down" if name.endswith("Down") else ""
                newname = name.rstrip(upDown)
                newname = f"{newname}_{self.chargeIdDict[charge]['id']}{upDown}"
            if not self.keepOtherChargeSyst and self.chargeIdDict[charge]["badId"] in newname: continue
            q = self.chargeIdDict[charge]["val"]
            hout = narf.hist_to_root(h[{"charge" : h.axes["charge"].index(q) if q != "sum" else hist.sum}])
            hout.SetName(newname+f"_{charge}")
            hout.Write()

    def writeHistWithCharges(self, h, name):
        hout = narf.hist_to_root(h)
        hout.SetName(f"{name}_{self.channels[0]}")
        hout.Write()
    
    def writeHist(self, h, name, setZeroStatUnc=False, decorrCharge=False, decorrByBin={}):
        if self.skipHist:
            return
        
        if self.project:
            axes = self.project[:]
            if "charge" in h.axes.name:
                axes.append("charge")
            h = h.project(*axes)

        if not self.nominalDim:
            self.nominalDim = h.ndim
            if self.nominalDim-self.writeByCharge > 3:
                raise ValueError("Cannot write hists with > 3 dimensions as combinetf does not accept THn")
        
        if h.ndim != self.nominalDim:
            raise ValueError(f"Histogram {name} does not have the correct dimensions. Found {h.ndim}, expected {self.nominalDim}")

        if setZeroStatUnc:
            hist_no_error = h.copy()
            hist_no_error.variances(flow=True)[...] = 0.
            h = hist_no_error

        hists = {}
        if decorrByBin:
            s = hist.tag.Slicer()
            # TODO: if passing two axes this has to make all combinations, currently it assumes only one axis
            for ax in decorrByBin.keys():
                for ibin in range(len(decorrByBin[ax]) -1):
                    upDown = "Up" if name.endswith("Up") else "Down" if name.endswith("Down") else ""
                    newname = name.rstrip(upDown)
                    newname = f"{newname}_{ax}{ibin}_{upDown}"
                    # TODO: this takes a slice, but I actually want to use nominal everywhere except in this particular bin
                    # should I do alt-nomi and then sum this back?
                    hists[newname] = narf.hist_to_root(h[{ax : s[complex(0, decorrByBin[ax][ibin]):complex(0, decorrByBin[ax][ibin+1])]}])
        else:
            hists[name] = h    
                
        for hname,hist in hists.items():
            if self.writeByCharge:
                self.writeHistByCharge(hist, hname, decorrCharge=decorrCharge)
            else:
                self.writeHistWithCharges(hist, hname)
