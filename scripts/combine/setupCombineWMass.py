#!/usr/bin/env python3
from Utilities import CardTool,datasets,theoryTools,selection
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/kelong/CombineStudies")
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

templateDir = "Templates/WMass"
cardTool = CardTool.CardTool(f"{args.outfolder}/Wmass.txt")
cardTool.setInputFile(args.inputFile)
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(f"{args.outfolder}/WMassCombineInput.root")
cardTool.setProcesses(datasets.processDict().keys())
cardTool.addFakeEstimate(selection.fakeHistABCD)
cardTool.setSignalOperation(selection.signalHistWmass)

cardTool.addSystematic("pdfNNPDF31", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x),
    outNames=theoryTools.pdfNames(cardTool, "NNPDF31"),
    mirror=True,
    group="pdfNNPDF31",
)
for name,num in zip(["effSystTnP", "effStatTnP",], [1, 624*2]):
    cardTool.addSystematic(name, 
        outNames=cardTool.mirrorNames(f"{name}"+("{i}" if num > 1 else ""), num),
        mirror=True,
        group=name,
    )
#if args.qcdByHelicity:
#    #cardTool.addSystematic("minnloQCDUncByHelicity", 
#    #    processes=cardTool.filteredProcesses(lambda x: "W" in x[0]),
#    #    outNames=qcdByHelicityLabels(),
#    #    group="QCDscaleByHelicity",
#    #)
#    #cardTool.addSystematic("qcdScale", 
#    #    processes=cardTool.filteredProcesses(lambda x: "Z" in x[0]),
#    #    outNames=theoryTools.qcdScaleNames(),
#    #    group="QCDscale",
#    #)
#else:
#    cardTool.addSystematic("qcdScale", 
#        processes=cardTool.filteredProcesses(lambda x: "W" in x[0] or "DY" in x[:2]),
#        outNames=theoryTools.qcdScaleNames(),
#        group="QCDscale",
#    )
cardTool.addSystematic("muonScaleFromMassWeights10MeV", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0]),
    outNames=cardTool.mirrorNames("dummyMuonScale5MeV", 1),
    group="massShift",
    scale=0.5,
)
cardTool.addSystematic("massWeight", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0]),
    outNames=theoryTools.massWeightNames(["massShift100MeV"]),
    group="massShift",
    groupFilter=lambda x: x == "massShift100MeV",
)
cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05)
cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
# This needs to be handled by shifting the norm before subtracting from the fakes
# cardTool.addSystematic("lumi", outNames=["", "lumiDown", "lumiUp"], group="luminosity")
# TODO: Allow to be appended to previous group
cardTool.writeOutput()

