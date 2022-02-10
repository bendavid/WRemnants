#!/usr/bin/env python3
from wremnants import CardTool,theoryTools
from wremnants import combineDatasets as data
from wremnants import histselections as sel
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

templateDir = "Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPileupW.txt")
cardTool.setInputFile(args.inputFile)
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(f"{args.outfolder}/LowPileupWCombineInput.root")
cardTool.setProcesses(data.processDictLowPileup(
    signalPtSplit=[0, 9, 20, 37, 60, 89, r"$\infty$"], removeUnsplit=True).keys()
)
# Best to keep this a multiple of 4
cardTool.setSpacing(28)
cardTool.setSignalOperation(sel.signalHistLowPileupW)
cardTool.setLoadProcesses(data.fillLowPileupProcDictFromRoot,
    extraArgs={"makeVariable" : 
        {"qTreco" : [0.0, 5.0, 9.0, 15.0, 20.0, 25.0, 37.0, 45.0, 60.0, 89.0, 150.],
        "iso" : [0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.6, 0.65, 100]}
    }
)
cardTool.addFakeEstimate(sel.fakeHistIsoRegion)
cardTool.setUnconstrainedProcs([cardTool.getFakeName()]+cardTool.filteredProcesses(lambda x: "ToMuNu" in x))
cardTool.setScaleMC(0.199)

cardTool.addSystematic("PDF", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x),
    outNames=theoryTools.pdfNames(cardTool, "NNPDF31", skipFirst=True),
    mirror=True,
    group="pdfNNPDF31",
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
cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05)
cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
# This needs to be handled by shifting the norm before subtracting from the fakes
# cardTool.addSystematic("lumi", outNames=["", "lumiDown", "lumiUp"], group="luminosity")
# TODO: Allow to be appended to previous group
cardTool.writeOutput()


