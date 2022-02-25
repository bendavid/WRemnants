#!/usr/bin/env python3
from wremnants import CardTool,theory_tools
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU
from wremnants import histselections as sel
import argparse
import os
import pathlib

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/kelong/CombineStudies")
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

def histName(proc, syst, variable="mt_reco_pf"):
    if proc == "Wmunu":
        variable = "mt_gen_reco_pf"
    base = f"{variable}_{proc}"
    return base if syst == "nominal" else f"{base}_{syst}_syst"

datagroups = datagroupsLowPU(args.inputFile)
templateDir = "Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPileupW.txt")

templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPileupW_{{chan}}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/LowPileupWCombineInput.root"))
cardTool.setDatagroups(datagroups)
cardTool.setBuildHistNameFunction(histName)
cardTool.setUnconstrainedProcs([cardTool.getFakeName(), "Wmunu"])

#cardTool.addSystematic("PDF", 
#    processes=cardTool.filteredProcesses(lambda x: "W" in x),
#    outNames=theoryTools.pdfNames(cardTool, "NNPDF31", skipFirst=True),
#    mirror=True,
#    group="pdfNNPDF31",
#)
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


