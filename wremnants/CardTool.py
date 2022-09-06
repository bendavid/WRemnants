from collections import OrderedDict
from utilities import output_tools,boostHistHelpers as hh
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

logging.basicConfig(level=logging.INFO)

def notImplemented(operation="Unknown"):
    raise NotImplementedError(f"Required operation '{operation}' is not implemented!")

class CardTool(object):
    def __init__(self, cardName="card.txt"):
        self.skipHist = False # don't produce/write histograms, file with them already exists
        self.outfile = None
        self.cardName = cardName
        self.systematics = {}
        self.lnNSystematics = {}
        self.procDict = OrderedDict()
        self.predictedProcs = []
        self.fakeEstimate = None
        self.addMirrorForSyst = {}
        self.channels = ["plus", "minus"]
        self.cardContent = {}
        self.cardGroups = {}
        self.nominalTemplate = ""
        self.spacing = 28
        self.fakeName = "Fake" # but better to set it explicitly
        self.dataName = "Data"
        self.nominalName = "nominal"
        self.datagroups = None
        self.unconstrainedProcesses = None
        self.noStatUncProcesses = []
        self.buildHistNameFunc = None
        self.histName = "x"
        self.pseudoData = None
        self.excludeSyst = None
        self.writeByCharge = True
        self.keepSyst = None # to override previous one with exceptions for special cases
        #self.loadArgs = {"operation" : "self.loadProcesses (reading hists from file)"}
        self.keepOtherChargeSyst = True
        self.chargeIdDict = {"minus" : {"val" : -1, "id" : "q0", "badId" : None},
                             "plus"  : {"val" : 1., "id" : "q1", "badId" : None}
                             }

    def skipHistograms(self):
        self.skipHist = True
        if len(self.noStatUncProcesses):
            logging.info("Attention: histograms are not saved according to input options, thus statistical uncertainty won't be zeroed")
        
    def setSkipOtherChargeSyst(self):
        self.keepOtherChargeSyst = False
        self.chargeIdDict["plus"]["badId"] = "q0"
        self.chargeIdDict["minus"]["badId"] = "q1"

    def setProcsNoStatUnc(self, procs, resetList=True):
        if self.skipHist:
            logging.info("Attention: trying to set statistical uncertainty to 0 for some processes, but histograms won't be saved according to input options")
        if resetList:
            self.noStatUncProcesses = []
        if isinstance(procs, str):
            self.noStatUncProcesses.append(procs)
        elif isinstance(procs, list):
            self.noStatUncProcesses.extend(procs)
        else:
            raise ValueError("In setNoStatUncForProcs(): expecting string or list argument")

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
                logging.info(f"   Excluding nuisance: {name}")
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

    def setPseudodata(self, pseudodata):
        self.pseudoData = pseudodata

    # Needs to be increased from default for long proc names
    def setSpacing(self, spacing):
        self.spacing = spacing

    def setDataName(self, name):
        self.dataName = name

    def setDatagroups(self, datagroups):
        self.datagroups = datagroups 
        
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
            return self.predictedProces
        return list(filter(lambda x: x != self.dataName, self.procDict.keys()))

    def setHistName(self, histName):
        self.histName = histName

    def isData(self, procName):
        return any([x.is_data for x in self.datagroups.groups[procName]["members"]])

    def isMC(self, procName):
        return not self.isData(procName)

    def addFakeEstimate(self, estimate):
        self.fakeEstimate = estimate

    def setProcesses(self, processes):
        self.procDict = OrderedDict([ (proc, {}) for proc in processes])

    def filteredProcesses(self, filterExpr):
        return list(filter(filterExpr, self.datagroups.processes()))

    def allMCProcesses(self):
        return self.filteredProcesses(lambda x: self.isMC(x))

    def mirrorNames(self, baseName, size, offset=0):
        names = [""]*offset + [f"{baseName.format(i=i%size)}{'Up' if i % 2 else 'Down'}" for i in range(size*2)]
        return names

    def addLnNSystematic(self, name, size, processes, group=None, groupFilter=None):
        if not self.isExcludedNuisance(name):
            self.lnNSystematics.update({name : {"size" : size, "processes" : processes, "group" : group, "groupFilter" : groupFilter}})

    # action will be applied to the sum of all the individual samples contributing, actionMap should be used
    # to apply a separate action per process. this is needed for example for the scale uncertainty split
    # by pt or helicity
    def addSystematic(self, name, systAxes, outNames=None, skipEntries=None, labelsByAxis=None, 
                        baseName="", mirror=False, scale=1, processes=None, group=None, noConstraint=False,
                        action=None, actionArgs={}, actionMap={}, systNameReplace=[], groupFilter=None, passToFakes=False, 
                        rename=None, splitGroup={}):

        # Need to make an explicit copy of the array before appending
        procs_to_add = [x for x in (self.allMCProcesses() if not processes else processes)]
        if passToFakes and self.getFakeName() not in procs_to_add:
            procs_to_add.append(self.getFakeName())

        if action and actionMap:
            raise ValueError("Only one of action and actionMap args are allowed")

        self.systematics.update({
            name if not rename else rename : { "outNames" : [] if not outNames else outNames,
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
                     "actionMap" : actionMap,
                     "actionArgs" : actionArgs,
                     "systNameReplace" : systNameReplace,
                     "noConstraint" : noConstraint,
                     "skipEntries" : [] if not skipEntries else skipEntries,
                     "name" : name,
            }
        })

    def setMirrorForSyst(self, syst, mirror=True):
        self.systematics[syst]["mirror"] = mirror

    def systLabelForAxis(self, axLabel, entry):
        if axLabel == "mirror" or axLabel == "downUpVar":
            return 'Up' if entry else 'Down'
        if "{i}" in axLabel:
            return axLabel.format(i=entry)
        return axLabel+str(entry)

    def excludeSystEntry(self, entry, skipEntries):
        for skipEntry in skipEntries:
            # Can use -1 to exclude all values of an axis
            if all(y == -1 or x == y for x,y in zip(entry, skipEntry)):
                return True
        return False

    def systHists(self, hvar, syst):
        if syst == self.nominalName:
            return ([self.nominalName], [hvar])

        systInfo = self.systematics[syst] 
        systAxes = systInfo["systAxes"]
        systAxesLabels = systAxes
        if "labelsByAxis" in systInfo:
            systAxesLabels = systInfo["labelsByAxis"]

        # Jan: moved above the mirror action, as this action can cause mirroring
        if systInfo["action"]:
            hvar = systInfo["action"](hvar, **systInfo["actionArgs"])
            
        axNames = systAxes[:]
        axLabels = systAxesLabels[:]
        if hvar.axes[-1].name == "mirror":
            axNames.append("mirror")
            axLabels.append("mirror")

        if not all([name in [ax.name for ax in hvar.axes] for name in axNames]):
            raise ValueError(f"Failed to find axis names '{str(systAxes)} in hist. " \
                f"Axes in hist are {str([ax.name for ax in hvar.axes])}")
        entries = list(itertools.product(*[range(hvar.axes[ax].size) for ax in axNames]))
        
        if len(systInfo["outNames"]) == 0:
            for entry in entries:
                if "skipEntries" in systInfo and self.excludeSystEntry(entry, systInfo["skipEntries"]):
                    systInfo["outNames"].append("")
                else:
                    name = systInfo["baseName"]
                    name += "".join([self.systLabelForAxis(al, entry[i]) for i,al in enumerate(axLabels)])
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
                raise RuntimeError(f"All entries for syst {syst} were skipped!")

        variations = []
        for entry in entries:
            sel = {ax : binnum for ax,binnum in zip(axNames, entry)}
            variations.append(hvar[sel])
        return systInfo["outNames"], variations            

    def variationName(self, proc, name):
        if name == self.nominalName:
            return f"{self.histName}_{proc}"
        else:
            return f"{self.histName}_{proc}_{name}"

    # This logic used to be more complicated, leaving the function here for now even
    # though it's trivial
    def addMirror(self, h, proc, syst):
        return syst != self.nominalName and self.systematics[syst]["mirror"]

    def writeForProcess(self, h, proc, syst):
        if self.addMirror(h, proc, syst):
            hnom = self.procDict[proc][self.nominalName]
            h = hh.extendHistByMirror(h, hnom)
        # Otherwise this is a processes not affected by the variation, don't write it out,
        # it's only needed for the fake subtraction
        logging.info(f"Writing systematic {syst} for process {proc}")
        var_names, variations = self.systHists(h, syst) 
        if len(var_names) != len(variations):
            logging.warning("The number of variations doesn't match the number of names for "
                f"process {proc}, syst {syst}. Found {len(var_names)} names and {len(variations)} variations.")
        setZeroStatUnc = False
        if proc in self.noStatUncProcesses:
            logging.info(f"Zeroing statistical uncertainty for process {proc}")
            setZeroStatUnc = True
        for name, var in zip(var_names, variations):
            if name != "":
                self.writeHist(var, self.variationName(proc, name), setZeroStatUnc=setZeroStatUnc)

    def addPseudodata(self, processes):
        self.datagroups.loadHistsForDatagroups(
            baseName=self.pseudoData, syst="", label=self.pseudoData,
            procsToRead=processes)
        hists = [self.procDict[proc][self.pseudoData] for proc in processes]
        hdata = hh.sumHists(hists)
        # Kind of hacky, but in case the alt hist has uncertainties
        for systAxName in ["systIdx", "tensor_axis_0"]:
            if systAxName in [ax.name for ax in hdata.axes]:
                hdata = hdata[{systAxName : 0 }] 
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
        else:
            self.outfile = outfile

    def writeOutput(self):
        self.datagroups.loadHistsForDatagroups(
            baseName=self.histName, syst=self.nominalName, label=self.nominalName)
        self.procDict = self.datagroups.getDatagroups()
        self.writeForProcesses(self.nominalName, processes=self.procDict.keys(), label=self.nominalName)
        self.loadNominalCard()
        if self.pseudoData:
            self.addPseudodata(self.predictedProcesses())

        self.writeLnNSystematics()
        for syst in self.systematics.keys():
            if self.isExcludedNuisance(syst): continue
            systMap=self.systematics[syst]
            systName = syst if not systMap["name"] else systMap["name"]
            processes=systMap["processes"]
            self.datagroups.loadHistsForDatagroups(self.histName, systName, label="syst",
                    procsToRead=processes, forceNonzero=systName != "qcdScaleByHelicity",
                    preOpMap=systMap["actionMap"], preOpArgs=systMap["actionArgs"])
            self.writeForProcesses(syst, label="syst", processes=processes)    
        output_tools.writeMetaInfoToRootFile(self.outfile, exclude_diff='notebooks')
        if self.skipHist:
            logging.info("Histograms will not be written because 'skipHist' flag is set to True")
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
        for name,info in self.lnNSystematics.items():
            include = [(str(info["size"]) if x in info["processes"] else "-").ljust(self.spacing) for x in nondata]
            group = info["group"]
            groupFilter = info["groupFilter"]
            for chan in self.channels:
                if self.keepOtherChargeSyst or self.chargeIdDict[chan]["badId"] not in name:
                    self.cardContent[chan] += f'{name.ljust(self.spacing)}lnN{" "*(self.spacing-3)}{"".join(include)}\n'
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
        if type(scale) != dict:
            include = [(str(scale) if x in procs else "-").ljust(self.spacing) for x in nondata]

        splitGroupDict = systInfo["splitGroup"]
        shape = "shape" if not systInfo["noConstraint"] else "shapeNoConstraint"

        # Deduplicate while keeping order
        systNames = list(dict.fromkeys(names))
        systnamesPruned = [s for s in systNames if not self.isExcludedNuisance(s)]
        systNames = systnamesPruned[:]
        for systname in systNames:
            if type(scale) == dict:
                for reg in scale.keys():
                    if re.match(reg, systname):
                        thiscale = str(scale[reg])
                        include = [(thiscale if x in procs else "-").ljust(self.spacing) for x in nondata]
                        break # exit this inner loop when match is found, to save time
            shape = "shape" if not systInfo["noConstraint"] else "shapeNoConstraint"
            for chan in self.channels:
                # do not write systs which should only apply to other charge, to simplify card
                if self.keepOtherChargeSyst or self.chargeIdDict[chan]["badId"] not in systname:
                    self.cardContent[chan] += f"{systname.ljust(self.spacing)}{shape.ljust(self.spacing)}{''.join(include)}\n"
                # unlike for LnN systs, here it is simpler to act on the list of these systs to form groups, rather than doing it syst by syst 
                if group:
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
                "channelPerProc" : chan.ljust(self.spacing)*nprocs,
                "processes" : "".join([x.ljust(self.spacing) for x in procs]),
                "labels" : "".join([str(x).ljust(self.spacing) for x in self.processLabels()]),
                # Could write out the proper normalizations pretty easily
                "rates" : "-1".ljust(self.spacing)*nprocs,
                "inputfile" : self.outfile if type(self.outfile) == str  else self.outfile.GetName(),
                "dataName" : self.dataName,
                "histName" : self.histName,
                "pseudodataHist" : self.pseudoData+"_sum" if self.pseudoData else f"{self.histName}_{self.dataName}"
            }
            self.cardContent[chan] = output_tools.readTemplate(self.nominalTemplate, args)
            self.cardGroups[chan] = ""
            
    def writeHistByCharge(self, h, name):
        for charge in self.channels:
            if not self.keepOtherChargeSyst and self.chargeIdDict[charge]["badId"] in name: continue
            q = self.chargeIdDict[charge]["val"]
            hout = narf.hist_to_root(h[{"charge" : h.axes["charge"].index(q)}])
            hout.SetName(name+f"_{charge}")
            hout.Write()

    def writeHistWithCharges(self, h, name):
        hout = narf.hist_to_root(h)
        hout.SetName(f"{name}_{self.channels[0]}")
        hout.Write()
    
    def writeHist(self, h, name, setZeroStatUnc=False):
        if self.skipHist:
            return
        if setZeroStatUnc:
            hist_no_error = h.copy()
            hist_no_error.variances(flow=True)[...] = 0.
            h = hist_no_error
        if self.writeByCharge:    
            self.writeHistByCharge(h, name)
        else:
            self.writeHistWithCharges(h, name)
