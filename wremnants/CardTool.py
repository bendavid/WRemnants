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
        self.signalOp = lambda h: h
        self.fakeName = "Fake"
        self.dataName = "data_obs"
        self.loadProcesses = notImplemented
        self.unconstrainedProcesses = None
        self.loadArgs = {"operation" : "self.loadProcesses (reading hists from file)"}

    # Function call to load hists for processes (e.g., read from a ROOT file)
    # Extra args will be passed to each call
    def setLoadProcesses(self, loadProcesses, extraArgs={}):
        self.loadProcesses = loadProcesses
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

    def setNominalTemplate(self, template):
        if not os.path.isfile(template):
            raise IOError(f"Template file {template} is not a valid file")
        self.nominalTemplate = template

    def predictedProcesses(self):
        if self.predictedProcs:
            return self.predictedProces
        return list(filter(lambda x: x != self.dataName, self.procDict.keys()))

    def setInputFile(self, infile):
        if type(infile) == str:
            self.inputFile = ROOT.TFile(infile)

    def addFakeEstimate(self, estimate):
        self.fakeEstimate = estimate

    # Use to select signal region in ND hist
    def setSignalOperation(self, op):
        self.signalOp = op

    def setProcesses(self, processes):
        self.procDict = OrderedDict([ (proc, {}) for proc in processes])

    def addProcess(self, process, info={}):
        self.procDict.update({process : info})

    def filteredProcesses(self, filterExpr):
        return list(filter(filterExpr, self.procDict.keys()))

    def allMCProcesses(self):
        return self.filteredProcesses(lambda x: self.isMC(proc))

    def mirrorNames(self, baseName, size, offset=0):
        names = [""]*offset + [f"{baseName.format(i=i%size)}{'Up' if i < size else 'Down'}" for i in range(size*2)]
        return names

    def addLnNSystematic(self, name, size, processes):
        self.lnNSystematics.update({name : {"size" : size, "processes" : processes}})

    def addSystematic(self, name, outNames, mirror=False, scale=1, processes=None, group=None, groupFilter=None):
        if not processes:
            processes = self.predictedProcesses()
        self.systematics.update({
            name : { "outNames" : outNames,
                "processes" : processes,
                "group" : group,
                "groupFilter" : groupFilter,
                "scale" : [scale]*len(outNames) if not isinstance(scale, collections.abc.Sequence) else scale,
                "mirror" : mirror
            }
        })

    def addToFakeHist(self, hvar, syst):
        label = "syst" if syst != "nominal" else "nominal"
        self.procDict[self.fakeName][label] = hh.addHists(self.procDict[self.fakeName][label], -1*hvar)

    def setMirrorForSyst(self, syst, mirror=True):
        self.systematics[syst]["mirror"] = mirror

    def systHists(self, hvar):
        if hvar.axes[-1].name != "systIdx":
            raise ValueError(f"Hist has no axis systIdx. Hist is {hvar}")
        variations = [hvar[{"systIdx" : i}] for i in range(hvar.axes["systIdx"].size)]
        return variations            

    def variationNames(self, proc, syst):
        if syst == "nominal":
            return [f"x_{proc}"]
        else:
            return [f"x_{proc}_{name}" for name in self.systematics[syst]["outNames"]]

    def addMirror(self, h, proc, syst):
        isFake = proc == self.fakeName
        mirrorSyst = syst != "nominal" and self.systematics[syst]["mirror"]
        isSyst = h.axes[-1].name == "systIdx"
        return not isFake and mirrorSyst and isSyst

    def writeForProcess(self, h, proc, syst):
        isVar = h.axes[-1].name == "systIdx"
        if self.addMirror(h, proc, syst):
            hnom = self.procDict[proc]["nominal"]
            h = hh.extendHistByMirror(h, hnom)
        if self.fakeEstimate and self.isMC(proc):
            self.addToFakeHist(h, syst)
        # Otherwise this is a processes not affected by the variation, don't write it out,
        # it's only needed for the fake subtraction
        if isVar or syst == "nominal":
            logging.info(f"Writing systematic {syst} for process {proc}")
            hwrite = self.signalOp(h) if proc != self.fakeName else self.fakeEstimate(h)
            var_names = self.variationNames(proc, syst) 
            variations = self.systHists(hwrite) if isVar else [hwrite]
            if len(var_names) != len(variations):
                logging.warning("The number of variations doesn't match the number of names for "
                    f"syst {syst}. Found {len(var_names)} names and {len(variations)} variations.")
            for name, var in zip(var_names, variations):
                if name != "":
                    self.writeHist(var, name)

    def isMC(self, proc):
        return proc not in [self.dataName, self.fakeName]

    def writeForAllProcesses(self, syst):
        if self.fakeEstimate:
            label = "syst" if syst != "nominal" else "nominal"
            dataHist = self.procDict[self.dataName]["nominal"]
            if not dataHist:
                raise RuntimeError(f"Failed to read data hist {self.dataName} from file")
            self.addProcess(self.fakeName, {label : dataHist.copy()})
        for process, info in self.procDict.items():
            if syst == "nominal" and not syst in info:
                raise RuntimeError(f"Failed to load nominal hist for process {process}")
            hvar = info["syst"] if ("syst" in info and info["syst"]) else info["nominal"]
            if not hvar:
                raise RuntimeError(f"Failed to load hist for process {process}, systematic {syst}")
            scale = self.scale if self.isMC(process) else 1
            self.writeForProcess(scale*hvar, process, syst)
        if syst != "nominal":
            self.fillCardWithSyst(syst)

    def setOutfile(self, outfile):
        if type(outfile) == str:
            self.outfile = ROOT.TFile(outfile, "recreate")
        else:
            self.outfile = outfile

    def writeOutput(self):
        self.loadProcesses(self.procDict, self.inputFile, syst="nominal", **self.loadArgs)
        self.writeForAllProcesses("nominal")
        self.loadNominalCard()
        self.writeLnNSystematics()
        for syst in self.systematics.keys():
            self.loadProcesses(self.procDict, self.inputFile, syst=syst, label="syst", 
                processes=self.systematics[syst]["processes"], **self.loadArgs)
            self.writeForAllProcesses(syst)
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
            label = "group" if not "massShift" in systNames[0] else "noiGroup"
            filt = self.systematics[syst]["groupFilter"]
            members = " ".join(systNames if not filt else filter(filt, systNames))
            self.cardGroups += f"\n{group} {label} = {members}"

    def setUnconstrainedProcs(self, procs):
        self.unconstrainedProcesses = procs
        print(procs)

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
