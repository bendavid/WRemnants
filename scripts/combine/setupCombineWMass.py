#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools
from wremnants import histselections as sel
from wremnants.datasets.datagroups import datagroups2016
import argparse
import os
import pathlib

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/kelong/CombineStudies")
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--qcdScale", choices=["byHelicityPt", "byPt", "integrated"], default="byHelicityPt", 
        help="Decorrelation for QCDscale (additionally always by charge)")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

datagroups = datagroups2016(args.inputFile)
templateDir = f"{scriptdir}/Templates/WMass"
cardTool = CardTool.CardTool(f"{args.outfolder}/Wmass_{{chan}}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/WMassCombineInput.root"))
cardTool.setDatagroups(datagroups)
cardTool.setSpacing(36)

cardTool.addSystematic("pdfNNPDF31", 
    processes=cardTool.filteredProcesses(lambda x: x[0] == "W" or x == "Fake"),
    mirror=True,
    group="pdfNNPDF31",
    systAxes=["tensor_axis_0"],
    labelsByAxis=["pdf{i}NNPDF31"],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    # -1 means all possible values of the mirror axis
    skipEntries=[(0, -1)],
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

scaleSystAxes = ["chargeVgen", "muRfact", "muFfact"] 
scaleLabelsByAxis = ["q", "muR", "muF"]
scaleGroupName = "QCDscale"
scaleActionArgs = {"sum_ptV" : True, "sum_helicity" : True}
scaleSkipEntries = [(-1, 1, 1), (-1, 0, 2), (-1, 2, 0)]
# TODO: reuse some code here
if "Helicity" in args.qcdScale:
    scaleActionArgs.pop("sum_helicity")
    scaleGroupName += "ByHelicity"
    scaleSystAxes.insert(0, "helicity")
    scaleLabelsByAxis.insert(0, "Coeff")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries]
if "Pt" in args.qcdScale:
    scaleActionArgs.pop("sum_ptV")
    scaleGroupName += "ByPtV"
    scaleSystAxes.insert(0, "ptVgen")
    scaleLabelsByAxis.insert(0, "genPtV")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries]

cardTool.addSystematic("qcdScaleByHelicity", 
    action=syst_tools.scale_helicity_hist_to_variations,
    actionArgs=scaleActionArgs,
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0]),
    group=scaleGroupName,
    systAxes=scaleSystAxes,
    labelsByAxis=scaleLabelsByAxis,
    # Exclude all combinations where muR = muF = 1 (nominal) or where
    # they are extreme values (ratio = 4 or 1/4)
    skipEntries=scaleSkipEntries,
    # This is hacky but it's the best idea I have for now...
    systNameReplace=[("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
        ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp")],
    baseName="QCDscale_",
    )

cardTool.addSystematic("muonScaleSyst", 
    processes=cardTool.filteredProcesses(lambda x: "W" in x[0] and "mu" in x),
    group="muonScale",
    baseName="CMS_scale_m_",
    systAxes=["downUpVar", "scaleEtaSlice"],
    labelsByAxis=["downUpVar", "ieta"],
)
cardTool.addSystematic("muonL1PrefireSyst", 
    processes=cardTool.filteredProcesses(lambda x: x != "Data"),
    group="muonPrefire",
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
    noConstraint=True,
    systAxes=["tensor_axis_0"],
)
cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05)
cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
# This needs to be handled by shifting the norm before subtracting from the fakes
# cardTool.addSystematic("lumi", outNames=["", "lumiDown", "lumiUp"], group="luminosity")
# TODO: Allow to be appended to previous group
cardTool.writeOutput()

