#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z
from wremnants import histselections as sel
import argparse
import os
import pathlib

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/jaeyserm/CombineStudies")
parser.add_argument("-i", "--inputFile", type=str, default="")
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
parser.add_argument("--qcdScale", choices=["byHelicityPt", "byPt", "integrated", "byHelicity"], default="byHelicityPt", 
        help="Decorrelation for QCDscale (additionally always by charge)")
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None, required=True)
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

if args.inputFile == "": args.inputFile = "mz_lowPU_%s.pkl.lz4" % args.flavor

datagroups = datagroupsLowPU_Z(args.inputFile, flavor=args.flavor)
if args.flavor == "mumu":
    dataName = "SingleMuon"
    groups = ["TTbar", "EWK", "SingleMuon", "DYmumu"]
if args.flavor == "ee":
    dataName = "SingleElectron"
    groups = ["TTbar", "EWK", "SingleElectron", "DYee"]

recoBins = [0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 10000]
genBins = [0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 10000]
unconstrainedProcs = [] # all gen bins
unconstrainedProc = "DYmumu" if args.flavor == "mumu" else "DYee"
doInclusive = False
if not doInclusive:
    proc_base = dict(datagroups.groups[unconstrainedProc]) 
    for i in range(len(genBins)-1): # add gen bin processes
    
        proc_name = "DY_genBin%d" % (i+1)
        proc_genbin = dict(proc_base)
        proc_genbin['signalOp'] = lambda x, i=i: x[{"recoil_gen" : i}] # index correct? Start from 1?
        datagroups.groups[proc_name] = proc_genbin
        unconstrainedProcs.append(proc_name)

    # remove Zmumu processes
    del datagroups.groups[unconstrainedProc]
    
    # remove non-used procs (hack)
    toDel = []
    for group in datagroups.groups: 
        if not group in groups+unconstrainedProcs: toDel.append(group)
    for group in toDel: del datagroups.groups[group]
 


templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPU_Z{args.flavor}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/LowPU_Z{args.flavor}.root"))
cardTool.setDatagroups(datagroups)
cardTool.setHistName("reco_mll") # for signal gen_reco_mll
cardTool.setChannels(["all"])
cardTool.setDataName(dataName)
cardTool.setUnconstrainedProcs(unconstrainedProcs)


pdfName = theory_tools.pdfMap["nnpdf31"]["name"]
cardTool.addSystematic(pdfName, 
    processes=unconstrainedProcs,
    mirror=True,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    labelsByAxis=[pdfName.replace("pdf", "pdf{i}")],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    # -1 means all possible values of the mirror axis
    skipEntries=[(0, -1)],
)

cardTool.addSystematic(f"alphaS002{pdfName}", 
    processes=unconstrainedProcs,
    mirror=False,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    outNames=[pdfName+"AlphaSUp", pdfName+"AlphaSDown"],
    scale=0.75,
)



scaleSystAxes = ["chargeVgen", "muRfact", "muFfact"]  ##  "ptVgen", "chargeVgen", "helicityWeight_tensor"
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
    processes=unconstrainedProcs,
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


#cardTool.addSystematic("PDF", 
#    processes=cardTool.filteredProcesses(lambda x: x[0] == "W" or x == "Fake"),
#    mirror=True,
#    group="pdfNNPDF31",
#    systAxes=["systAx"],
#    labelsByAxis=["pdf{i}NNPDF31"],
#    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
#    skipEntries=[(0, 0), (0,1)],
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
#cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05)
#cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
#cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
# This needs to be handled by shifting the norm before subtracting from the fakes
# cardTool.addSystematic("lumi", outNames=["", "lumiDown", "lumiUp"], group="luminosity")
# TODO: Allow to be appended to previous group

cardTool.addLnNSystematic("CMS_lumi", processes=cardTool.allMCProcesses(), size=1.02)

cardTool.writeOutput()


