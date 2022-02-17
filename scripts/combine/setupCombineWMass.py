#!/usr/bin/env python3
from wremnants import CardTool,theory_tools
from wremnants import histselections as sel
from wremnants.datasets.datagroups import datagroups2016
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

datagroups = datagroups2016(args.inputFile)
templateDir = f"{scriptdir}/Templates/WMass"
cardTool = CardTool.CardTool(f"{args.outfolder}/Wmass_{{chan}}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/WMassCombineInput.root"))
cardTool.setDatagroups(datagroups)

#TODO: Change the mirrorNames function so it gives the right order for multiple axes
cardTool.addSystematic("pdfNNPDF31", 
    processes=cardTool.filteredProcesses(lambda x: x[0] == "W" or x == "Fake"),
    mirror=True,
    group="pdfNNPDF31",
    systAxes=["tensor_axis_0"],
    labelsByAxis=["pdf{i}NNPDF31"],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    skipEntries=[(0, 0), (0,1)],
)
for name,num in zip(["effSystIsoTnP", "effStatTnP",], [2, 624*4]):
    axes = ["idiptrig-iso"] if num == 2 else ["SF eta", "SF pt", "SF charge", "idiptrig-iso"]
    axlabels = ["IDIPTrig"] if num == 2 else ["eta", "pt", "q", "Trig"]
    cardTool.addSystematic(name, 
        mirror=True,
        group=name,
        systAxes=axes,
        labelsByAxis=axlabels,
        baseName=name+"_",
        processes=cardTool.filteredProcesses(lambda x: x != "Data"),
    )
#if args.qcdByHelicity:
#    #cardTool.addSystematic("minnloQCDUncByHelicity", 
#    #    processes=cardTool.filteredProcesses(lambda x: "W" in x[0]),
#    #    outNames=qcdByHelicityLabels(),
#    #    group="QCDscaleByHelicity",
#    #)
#    #cardTool.addSystematic("qcdScale", 
#    #    processes=cardTool.filteredProcesses(lambda x: "Z" in x[0]),
#    #    outNames=theory_tools.qcdScaleNames(),
#    #    group="QCDscale",
#    #)
#else:
#    cardTool.addSystematic("qcdScale", 
#        processes=cardTool.filteredProcesses(lambda x: "W" in x[0] or "DY" in x[:2]),
#        outNames=theory_tools.qcdScaleNames(),
#        group="QCDscale",
#    )
cardTool.addSystematic("muonScaleSyst", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0] and "mu" in x),
    group="muonScale",
    baseName="CMS_scale_m_",
    systAxes=["downUpVar", "scaleEtaSlice"],
    labelsByAxis=["downUpVar", "ieta"],
)
cardTool.addSystematic("massWeight", 
    # TODO: Add the mass weights to the tau samples
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0] and "mu" in x),
    outNames=theory_tools.massWeightNames(["massShift100MeV"]),
    group="massShift",
    groupFilter=lambda x: x == "massShift100MeV",
    mirror=False,
    #TODO: Name this
    systAxes=["tensor_axis_0"],
)
cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05)
cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
# This needs to be handled by shifting the norm before subtracting from the fakes
# cardTool.addSystematic("lumi", outNames=["", "lumiDown", "lumiUp"], group="luminosity")
# TODO: Allow to be appended to previous group
cardTool.writeOutput()

