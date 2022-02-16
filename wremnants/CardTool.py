from collections import OrderedDict
from . import boostHistHelpers as hh
from . import OutputTools
import narf
import logging
import ROOT
import time
import numpy as np
import collections.abc
import os
import itertools

logging.basicConfig(level=logging.INFO)

def notImplemented(operation="Unknown"):
    raise NotImplementedError(f"Required operation '{operation}' is not implemented!")

class CardTool(object):
    def __init__(self, cardName="card.txt", scale=1.):
        self.cardName = cardName
        self.systematics = {}
        self.lnNSystematics = {}
        self.procDict = OrderedDict()
        self.predictedProcs = []
        self.fakeEstimate = None
        self.addMirrorForSyst = {}
        self.writeSeparateCharges = False
        self.cardContent = ""
        self.cardGroups = ""
        self.nominalTemplate = ""
        self.spacing = 20
        self.scale = scale
        self.fakeName = "Fake"
        self.dataName = "data_obs"
        self.nominalName = "nominal"
        self.datagroups = None
        self.unconstrainedProcesses = None
        #self.loadArgs = {"operation" : "self.loadProcesses (reading hists from file)"}

    # Function call to load hists for processes (e.g., read from a ROOT file)
    # Extra args will be passed to each call
    def setLoadDatagroups(self, datagroups, extraArgs={}):
        self.datagroups = datagroups
        self.loadArgs = extraArgs
    
    def getFakeName(self):
        return self.fakeName

    def setScaleMC(self, scale):
        self.scale = scale

    # Needs to be increased from default for long proc names
    def setSpacing(self, spacing):
        self.spacing = spacing

    def setDataName(self, name):
        self.dataName = name

    def setDatagroups(self, datagroups):
        self.datagroups = datagroups 

    def setNominalTemplate(self, template):
        if not os.path.isfile(template):
            raise IOError(f"Template file {template} is not a valid file")
        self.nominalTemplate = template

    def predictedProcesses(self):
        if self.predictedProcs:
            return self.predictedProces
        return list(filter(lambda x: x != self.dataName, self.procDict.keys()))

    #def setInputFile(self, infile):
    #    if type(infile) == str:
    #        self.inputFile = ROOT.TFile(infile)
    def isData(self, procInfo):
        return all([x.is_data for x in procInfo["members"]])

    def isMC(self, procInfo):
        return any([x.is_data for x in procInfo["members"]])

    def addFakeEstimate(self, estimate):
        self.fakeEstimate = estimate

    def setProcesses(self, processes):
        self.procDict = OrderedDict([ (proc, {}) for proc in processes])

    def filteredProcesses(self, filterExpr):
        return list(filter(filterExpr, self.datagroups.processes()))

    def allMCProcesses(self):
        return self.filteredProcesses(lambda x: self.isMC(proc))

    def mirrorNames(self, baseName, size, offset=0):
        names = [""]*offset + [f"{baseName.format(i=i%size)}{'Up' if i % 2 else 'Down'}" for i in range(size*2)]
        return names

    def addLnNSystematic(self, name, size, processes):
        self.lnNSystematics.update({name : {"size" : size, "processes" : processes}})

    def addSystematic(self, name, outNames, systAxes, mirror=False, scale=1, processes=None, group=None, groupFilter=None):
        if not processes:
            processes = self.predictedProcesses()
        self.systematics.update({
            name : { "outNames" : outNames,
                "processes" : processes,
                "systAxes" : systAxes,
                "group" : group,
                "groupFilter" : groupFilter,
                "scale" : [scale]*len(outNames) if not isinstance(scale, collections.abc.Sequence) else scale,
                "mirror" : mirror
            }
        })

    def setMirrorForSyst(self, syst, mirror=True):
        self.systematics[syst]["mirror"] = mirror

    # TODO: Actually use the syst axes
    def systHists(self, hvar, syst):
        systAxes = self.systematics[syst]["systAxes"]
        axNames = systAxes[:]
        if hvar.axes[-1].name == "mirror":
            axNames.append("mirror")

        entries = itertools.product(*[range(hvar.axes[ax].size) for ax in axNames])
        variations = []
        for entry in entries:
            sel = {ax : binnum for ax,binnum in zip(axNames, entry)}
            variations.append(hvar[sel])
        return variations            

    def variationNames(self, proc, syst):
        if syst == self.nominalName:
            return [f"x_{proc}"]
        else:
            return [f"x_{proc}_{name}" if name != "" else name for name in self.systematics[syst]["outNames"]]

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
        var_names = self.variationNames(proc, syst) 
        variations = self.systHists(h, syst) if syst != self.nominalName else [h]
        if len(var_names) != len(variations):
            logging.warning("The number of variations doesn't match the number of names for "
                f"process {proc}, syst {syst}. Found {len(var_names)} names and {len(variations)} variations.")
        for name, var in zip(var_names, variations):
            if name != "":
                self.writeHist(var, name)

    def writeForProcesses(self, syst, processes, label):
        for process in processes:
            hvar = self.procDict[process][label]
            if not hvar:
                raise RuntimeError(f"Failed to load hist for process {process}, systematic {syst}")
            scale = self.scale if process != self.dataName else 1
            self.writeForProcess(scale*hvar, process, syst)
        if syst != self.nominalName:
            self.fillCardWithSyst(syst)

    def setOutfile(self, outfile):
        if type(outfile) == str:
            self.outfile = ROOT.TFile(outfile, "recreate")
        else:
            self.outfile = outfile

    def writeOutput(self):
        self.procDict = self.datagroups.datagroupsForHist(histname=self.nominalName, label=self.nominalName)
        self.writeForProcesses(self.nominalName, processes=self.procDict.keys(), label=self.nominalName)
        self.loadNominalCard()
        self.writeLnNSystematics()
        for syst in self.systematics.keys():
            processes=self.systematics[syst]["processes"]
            self.procDict = self.datagroups.datagroupsForHist(histname=syst, label="syst",
                dataHist=self.nominalName, procsToRead=processes)
            self.writeForProcesses(syst, label="syst", processes=processes)
        self.writeCard()

    def writeCard(self):
        with open(self.cardName, "w") as card:
            card.write(self.cardContent)
            card.write("\n")
            card.write(self.cardGroups)

    #def fillCardFromTemplate()

    def writeLnNSystematics(self):
        nondata = self.predictedProcesses()
        for name,info in self.lnNSystematics.items():
            include = [(str(info["size"]) if x in info["processes"] else "-").ljust(self.spacing) for x in nondata]
            self.cardContent += f'{name.ljust(self.spacing)}lnN{" "*(self.spacing-3)}{"".join(include)}\n'

    def fillCardWithSyst(self, syst):
        scale = self.systematics[syst]["scale"]
        procs = self.systematics[syst]["processes"]
        nondata = self.predictedProcesses()
        names = [x[:-2] if "Up" in x[-2:] else (x[:-4] if "Down" in x[-4:] else x) 
                    for x in filter(lambda x: x != "", self.systematics[syst]["outNames"])]
        include = [[(str(s) if x in procs else "-").ljust(self.spacing) for x in nondata] for s in scale]

        # Deduplicate
        systNames = list(dict.fromkeys(names))
        for systname, inc in zip(systNames, include):
            self.cardContent += f"{systname.ljust(self.spacing)}{'shape'.ljust(self.spacing)}{''.join(inc)}\n"

        group = self.systematics[syst]["group"]
        if group:
            # TODO: Make more general
            label = "group" if not "massShift" in systNames[0] else "noiGroup"
            filt = self.systematics[syst]["groupFilter"]
            members = " ".join(systNames if not filt else filter(filt, systNames))
            self.cardGroups += f"\n{group} {label} = {members}"

    def setUnconstrainedProcs(self, procs):
        self.unconstrainedProcesses = procs

    def processLabels(self):
        nondata = np.array(self.predictedProcesses())
        labels = np.arange(len(nondata))+1
        issig = np.isin(nondata, self.unconstrainedProcesses)
        labels[issig] = -np.arange(np.count_nonzero(issig))-1
        return labels

    def loadNominalCard(self):
        channel = "all"
        procs = self.predictedProcesses()
        nprocs = len(procs)
        args = {
            "channel" :  channel,
            "channelPerProc" : channel.ljust(self.spacing)*nprocs,
            "processes" : "".join([x.ljust(self.spacing) for x in procs]),
            "labels" : "".join([str(x).ljust(self.spacing) for x in self.processLabels()]),
            # Could write out the proper normalizations pretty easily
            "rates" : "-1".ljust(self.spacing)*nprocs,
            "inputfile" : self.outfile.GetName(),
        }
        self.cardContent = OutputTools.readTemplate(self.nominalTemplate, args)
        
    def writeHistByCharge(self, h, name, charges=["plus", "minus"]):
        for charge in charges:
            c = clabels.index(charge)
            hout = narf.hist_to_root(h[{"charge" : c}])
            hout.SetName(name+f"_{charge}")
            hout.Write()

    def writeHistWithCharges(self, h, name):
        hout = narf.hist_to_root(h)
        hout.SetName(name)
        hout.Write()
    
    def writeHist(self, h, name):
        if self.writeSeparateCharges:
            self.writeHistByCharge(h, name)
        else:
            self.writeHistWithCharges(h, name)
