#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools,combine_helpers
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU
from utilities import common_lowPU as common
import argparse
import os
import pathlib
import hist
import numpy as np

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="combineResults/lowPU/differentialZ/")
parser.add_argument("-i", "--inputFile", type=str, default="")
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None, required=True)
parser.add_argument("--fitType", choices=["differential", "wmass", "wlike", "inclusive"], default="differential", 
        help="Fit type, defines POI and fit observable (recoil or mT)")
parser.add_argument("--xsec", dest="xsec", action='store_true', 
        help="Write card for masked xsec normalization")
parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")
parser.add_argument("--statOnly", dest="statOnly", action='store_true', 
        help="Stat-only cards")
parser.add_argument("--lumiScale", dest="lumiScale", help="Luminosity scale", type=float, default=1.0)
parser.add_argument("--doStatOnly", action="store_true", default=False, 
        help="Set up fit to get stat-only uncertainty (currently combinetf with -S 0 doesn't work)")
parser.add_argument("--qcdScale", choices=["byHelicityPtAndByPt", "byHelicityPt", "byHelicityCharge", "byPt", "byCharge", "integrated",], default="integrated", 
        help="Decorrelation for QCDscale (additionally always by charge). With 'byHelicityPtAndByPt' two independent histograms are stored, split and not split by helicities (for tests)")
parser.add_argument("--scetlibUnc", action='store_true', help="Include SCETlib uncertainties")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

if args.inputFile == "": args.inputFile = "lowPU_%s_%s.pkl.lz4" % (args.flavor, args.met)


datagroups = datagroupsLowPU(args.inputFile, flavor=args.flavor)
unconstrainedProcs = [] # POI processes
constrainedProcs = []   # constrained signal procs
bkgProcs = ["Other", "Ztautau"] if args.fitType == "wmass" else ["Top", "EWK"] # background procs
dataProc = "SingleMuon" if args.flavor == "mumu" else "SingleElectron"


histName = "reco_mll"
s = hist.tag.Slicer()
if args.fitType == "differential":
    proc_base = dict(datagroups.groups["Zmumu" if args.flavor == "mumu" else "Zee"]) # signal process (Zmumu or Zee)
    for i in range(len(common.axis_recoil_gen_ptZ)): # add gen bin processes to the dict
        proc_name = "Zmumu_genBin%d" % (i+1)
        proc_genbin = dict(proc_base)
        proc_genbin['signalOp'] = lambda x, i=i: x[{"recoil_gen" : i}]
        datagroups.groups[proc_name] = proc_genbin
        if args.fitType == "differential": unconstrainedProcs.append(proc_name)
elif args.fitType == "wmass": 
    proc_name = "Zmumu" if args.flavor == "mumu" else "Zee"
    constrainedProcs.append(proc_name)
    datagroups.groups[proc_name]['signalOp'] = lambda x: x[{"recoil_gen" : s[::hist.sum]}] if "recoil_gen" in x.axes.name else x
    constrainedProcs.append(proc_name)
elif args.fitType == "wlike":
    histName = "mT_corr_rec"
    constrainedProcs.append("Zmumu" if args.flavor == "mumu" else "Zee") # need sum over gen bins
elif args.fitType == "inclusive":
    unconstrainedProcs.append("Zmumu" if args.flavor == "mumu" else "Zee") # need sum over gen bins
   
suffix = ""
if args.xsec:
    suffix = "_xsec"
    bkgProcs = [] # for xsec norm card, remove all bkg procs but keep the data
  
# hack: remove non-used procs/groups, as there can be more procs/groups defined than defined above
# need to remove as cardTool takes all procs in the datagroups
toDel = []
for group in datagroups.groups: 
    if not group in constrainedProcs+unconstrainedProcs+bkgProcs+[dataProc]: toDel.append(group)
for group in toDel: del datagroups.groups[group]    

templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/lowPU_Z{args.flavor}_{args.met}_{args.fitType}{suffix}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/lowPU_Z{args.flavor}_{args.met}_{args.fitType}{suffix}.root"))
cardTool.setDatagroups(datagroups)
cardTool.setHistName(histName) 
cardTool.setNominalName(histName) 
cardTool.setChannels([f"{args.flavor}{suffix}"])
cardTool.setDataName(dataProc)
cardTool.setProcsNoStatUnc(procs=[])
cardTool.setSpacing(36)
cardTool.setWriteByCharge(False)
cardTool.setUnconstrainedProcs(unconstrainedProcs)
cardTool.setLumiScale(args.lumiScale)

Zmumu_procs = cardTool.filteredProcesses(lambda x: "Zmumu" in x)
Zmumu_procsIncTau = Zmumu_procs + ["Ztautau"]

if args.fitType == "wlike" or args.fitType == "wmass":

    cardTool.addSystematic("massWeight", 
        processes=constrainedProcs,
        outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=False),
        group="massShift",
        groupFilter=lambda x: x == "massShift100MeV",
        mirror=False,
        #TODO: Name this
        noConstraint=True,
        systAxes=["tensor_axis_0"],
    )

if args.doStatOnly:
    cardTool.addLnNSystematic("dummy", processes=cardTool.allMCProcesses(), size=1.001, group="dummy")
    cardTool.writeOutput()
    print("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
    quit()


pdfName = theory_tools.pdfMap["nnpdf31"]["name"]
pdfAction = {x : lambda h: h[{"recoil_gen" : s[::hist.sum]}] for x in Zmumu_procs if "gen" not in x},
cardTool.addSystematic(pdfName, 
    processes=Zmumu_procs,
    mirror=True,
    group=pdfName,
    actionMap=pdfAction,
    systAxes=["tensor_axis_0"],
    labelsByAxis=[pdfName.replace("pdf", "pdf{i}")],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    # -1 means all possible values of the mirror axis
    skipEntries=[(0, -1)],
)
cardTool.addSystematic(f"alphaS002{pdfName}", 
    processes=Zmumu_procs,
    mirror=False,
    actionMap=pdfAction,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    outNames=[pdfName+"AlphaSUp", pdfName+"AlphaSDown"],
    scale=0.75,
)

combine_helpers.add_scale_uncertainty(cardTool, args.qcdScale, constrainedProcs+unconstrainedProcs, to_fakes=False, use_hel_hist=True, scetlib=args.scetlibUnc)

if not args.xsec:

    cardTool.addSystematic("prefireCorr",
        processes=cardTool.allMCProcesses(),
        mirror = False,
        group="CMS_prefire17",
        baseName="CMS_prefire17",
        systAxes = ["downUpVar"],
        labelsByAxis = ["downUpVar"],
    )
    
    recoil_vars = ["target_para", "target_perp", "source_para", "source_perp", "target_para_bkg", "target_perp_bkg"]
    recoil_grps = ["recoil_stat", "recoil_stat", "recoil_stat", "recoil_stat", "recoil_syst", "recoil_syst"]
    for i, tag in enumerate(recoil_vars):
        cardTool.addSystematic("recoilSyst_%s" % tag,
            processes=Zmumu_procs,
            mirror = True,
            group = recoil_grps[i],
            systAxes = ["tensor_axis_0"],
            labelsByAxis = ["recoil_%s_{i}" % tag],
            #action=scale_recoil_hist_to_variations,
            #scale = 1./math.sqrt(args.lumiScale) if not "bkg" in tag else 1.0,
        )


    for lepEff in ["lepSF_HLT_DATA_stat", "lepSF_HLT_DATA_syst", "lepSF_HLT_MC_stat", "lepSF_HLT_MC_syst", "lepSF_ISO_stat", "lepSF_ISO_DATA_syst", "lepSF_ISO_MC_syst", "lepSF_IDIP_stat", "lepSF_IDIP_DATA_syst", "lepSF_IDIP_MC_syst"]:
        
        cardTool.addSystematic(lepEff,
            processes=cardTool.allMCProcesses(),
            mirror = True,
            group="CMS_lepton_eff",
            baseName=lepEff,
            systAxes = ["tensor_axis_0"],
            labelsByAxis = [""],
        )

    if args.fitType == "wmass":
        cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15, group="CMS_background")
    else:
        cardTool.addLnNSystematic("CMS_Top", processes=["TTbar"], size=1.06, group="CMS_background")
        cardTool.addLnNSystematic("CMS_VV", processes=["EWK"], size=1.16, group="CMS_background")
    
    cardTool.addLnNSystematic("CMS_lumi_lowPU", processes=cardTool.allMCProcesses(), size=1.02, group="CMS_lumi_lowPU")



cardTool.writeOutput()
