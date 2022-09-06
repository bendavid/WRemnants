#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools
from wremnants import histselections as sel
from wremnants.datasets.datagroups import datagroups2016
from utilities import common
import argparse
import os
import pathlib
import logging
import hist
import copy

logging.basicConfig(level=logging.INFO)

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/kelong/CombineStudies", help="Main output folder, with the root file storing all histograms and datacards for single charge")
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--qcdScale", choices=["byHelicityPtAndByPt", "byHelicityPt", "byHelicityCharge", "byPt", "byCharge", "integrated",], default="byHelicityPt", 
        help="Decorrelation for QCDscale (additionally always by charge). With 'byHelicityPtAndByPt' two independent histograms are stored, split and not split by helicities (for tests)")
parser.add_argument("--rebinPtV", type=float, default=0, nargs='*', help="Rebin axis with gen boson pt by this value (default does nothing)")
parser.add_argument("--wlike", action='store_true', help="Run W-like analysis of mZ")
parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
parser.add_argument("--pdf", type=str, default="nnpdf31", choices=theory_tools.pdfMap.keys(), help="PDF to use")
parser.add_argument("-b", "--fitObs", type=str, default="nominal", help="Observable to fit") # TODO: what does it do?
parser.add_argument("-p", "--pseudoData", type=str, help="Hist to use as pseudodata")
parser.add_argument("-x",  "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
parser.add_argument("-k",  "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
parser.add_argument("--qcdProcessName", dest="qcdProcessName" , type=str, default="Fake",   help="Name for QCD process")
parser.add_argument("--noStatUncFakes", dest="noStatUncFakes" , action="store_true",   help="Set bin error for QCD background templates to 0, to check MC stat uncertainties for signal only")
parser.add_argument("--noQCDscaleFakes", dest="noQCDscaleFakes" , action="store_true",   help="Do not apply QCd scale uncertainties on fakes, mainly for debugging")
parser.add_argument("--skipOtherChargeSyst", dest="skipOtherChargeSyst" , action="store_true",   help="Skip saving histograms and writing nuisance in datacard for systs defined for a given charge but applied on the channel with the other charge")
parser.add_argument("--skipSignalSystOnFakes", dest="skipSignalSystOnFakes" , action="store_true", help="Do not propagate signal uncertainties on fakes, mainly for checks.")
parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
parser.add_argument("--correlateEffStatIsoByCharge", action='store_true', help="Correlate isolation efficiency uncertanties between the two charges (by default they are decorrelated)")
parser.add_argument("--doStatOnly", action="store_true", default=False, help="Set up fit to get stat-only uncertainty (currently combinetf with -S 0 doesn't work)")
parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.makedirs(args.outfolder)

if args.noHist and args.noStatUncFakes:
    logging.warning("Option --noHist would override --noStatUncFakes. Please select only one of them")
    quit()
    
datagroups = datagroups2016(args.inputFile, wlike=args.wlike)
templateDir = f"{scriptdir}/Templates/WMass"
name = "WMass" if not args.wlike else "ZMassWLike"
cardTool = CardTool.CardTool(f"{args.outfolder}/{name}_{{chan}}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
if args.noHist:
    cardTool.skipHistograms()
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/{name}CombineInput.root"))
cardTool.setDatagroups(datagroups)
cardTool.setFakeName(args.qcdProcessName)
cardTool.setSpacing(36)
if args.noStatUncFakes:
    cardTool.setProcsNoStatUnc(procs=args.qcdProcessName, resetList=False)
cardTool.setCustomSystForCard(args.excludeNuisances, args.keepNuisances)
if args.skipOtherChargeSyst:
    cardTool.setSkipOtherChargeSyst()
if args.pseudoData:
    cardTool.setPseudodata(args.pseudoData)

passSystToFakes = not args.wlike and not args.skipSignalSystOnFakes
    
single_v_samples = cardTool.filteredProcesses(lambda x: x[0] in ["W", "Z"])
single_v_nonsig_samples = cardTool.filteredProcesses(lambda x: x[0] == ("W" if args.wlike else "Z"))
single_vmu_samples = list(filter(lambda x: "mu" in x, single_v_samples))
signal_samples = list(filter(lambda x: x[0] == ("Z" if args.wlike else "W"), single_vmu_samples))
signal_samples_inctau = list(filter(lambda x: x[0] == ("Z" if args.wlike else "W"), single_v_samples))

logging.info(f"All MC processes {cardTool.allMCProcesses()}")
logging.info(f"Single V samples: {single_v_samples}")
logging.info(f"Single V no signal samples: {single_v_nonsig_samples}")
logging.info(f"Signal samples: {signal_samples}")

pdfInfo = theory_tools.pdf_info_map(signal_samples[0], args.pdf)
pdfName = pdfInfo["name"]

# keep mass weights here as first systematic, in case one wants to run stat-uncertainty only with --doStatOnly
cardTool.addSystematic("massWeight", 
    # TODO: Add the mass weights to the tau samples ## FIXME: isn't it done?
    processes=signal_samples_inctau,
    outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=args.wlike),
    group="massShift",
    groupFilter=lambda x: x == "massShift100MeV",
    mirror=False,
    #TODO: Name this
    noConstraint=True,
    systAxes=["tensor_axis_0"],
    passToFakes=passSystToFakes,
)

if args.doStatOnly:
    # print a card with only mass weights and a dummy syst
    cardTool.addLnNSystematic("dummy", processes=["Other"] if args.wlike else ["Top", "Diboson"], size=1.001, group="dummy")
    cardTool.writeOutput()
    print("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
    quit()
    
if args.wlike:
    # TOCHECK: no fakes here, most likely
    cardTool.addLnNSystematic("luminosity", processes=cardTool.allMCProcesses(), size=1.012, group="luminosity")
else:
    cardTool.addSystematic("luminosity",
                           processes=cardTool.allMCProcesses(),
                           outNames=["lumiDown", "lumiUp"],
                           group="luminosity",
                           systAxes=["downUpVar"],
                           labelsByAxis=["downUpVar"],
                           passToFakes=passSystToFakes)

if pdfInfo["combine"] == "symHessian":
    cardTool.addSystematic(pdfName, 
        processes=single_v_samples,
        mirror=True,
        group=pdfName,
        systAxes=["tensor_axis_0"],
        labelsByAxis=[pdfName.replace("pdf", "pdf{i}")],
        # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
        # -1 means all possible values of the mirror axis
        skipEntries=[(0, -1)],
        passToFakes=passSystToFakes,
    )
else:
    cardTool.addSystematic(pdfName, 
        processes=single_v_samples,
        mirror=False,
        group=pdfName,
        systAxes=["tensor_axis_0"],
        outNames=theory_tools.pdfNamesAsymHessian(pdfInfo["entries"]),
        passToFakes=passSystToFakes,
        scale=pdfInfo["scale"] if "scale" in pdfInfo else 1,
    )

cardTool.addSystematic(f"alphaS002{pdfName}", 
    processes=single_v_samples,
    mirror=False,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    outNames=[pdfName+"AlphaSUp", pdfName+"AlphaSDown"],
    scale=0.75, # TODO: this depends on the set, should be provided in theory_tools.py
    passToFakes=passSystToFakes,
)
if not args.noEfficiencyUnc:
    for name,num in zip(["effSystTnP", "effStatTnP",], [2, 624*4]):
        ## TODO: this merged implementation for the effstat makes it very cumbersome to do things differently for iso and trigidip!!
        ## the problem is that I would need custom actions inside based on actual nuisance names, which needs to be commanded from outside, and this is not straightforward
        axes = ["idiptrig-iso"] if num == 2 else ["SF eta", "SF pt", "SF charge", "idiptrig-iso"]
        axlabels = ["Trig"] if num == 2 else ["eta", "pt", "q", "Trig"]  # WARNING: Trig0/Trig1 actually stands for trigger/isolation, the axis name was intended to indicate the order "axis_idiptrig_iso"
        cardTool.addSystematic(name, 
            mirror=True,
            group="muon_eff_syst" if "Syst" in name else "muon_eff_stat", # TODO: for now better checking them separately
            systAxes=axes,
            labelsByAxis=axlabels,
            baseName=name+"_",
            processes=cardTool.allMCProcesses(),
            passToFakes=passSystToFakes,
            systNameReplace=[("Trig0", "IDIPTrig"), ("q0Trig1", "Iso"), ("q1Trig1", "Iso")] if args.correlateEffStatIsoByCharge else [("Trig0", "IDIPTrig"), ("Trig1", "Iso")], # replace with better names
            scale=1.0 if "Syst" in name  else {".*effStatTnP.*Iso" : "1.414", ".*effStatTnP.*IDIPTrig" : "1.0"} if not args.correlateEffStatIsoByCharge else 1.0 # only for iso, scale up by sqrt(2) when decorrelating between charges and efficiencies were derived inclusively
        )

inclusiveScale = args.qcdScale == "integrated"
helicity = "Helicity" in args.qcdScale

scale_hist = "qcdScale"
scaleSystAxes = ["muRfact", "muFfact"] 
scaleLabelsByAxis = ["muR", "muF"]
scaleGroupName = "QCDscale"
# Exclude all combinations where muR = muF = 1 (nominal) or where
# they are extreme values (ratio = 4 or 1/4)
scaleSkipEntries = [(1, 1), (0, 2), (2, 0)]
# This is hacky but it's the best idea I have for now...
systNameReplaceVec = [("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
                      ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp"),
                      ("genQ0", "genVminus"), ("genQ1", "genVplus")]
scaleActionArgs = {"sum_axis" : ["ptVgen"]}

scale_action_map = {proc : syst_tools.scale_helicity_hist_to_variations for proc in common.vprocs}

cardTool.addSystematic("qcdScale",
    rename="qcdScaleNonSignal",
    action=syst_tools.scale_helicity_hist_to_variations,
    actionArgs=scaleActionArgs,
    processes=single_v_nonsig_samples,
    group=f"qcdScale{'W' if args.wlike else 'Z'}",
    systAxes=copy.deepcopy(scaleSystAxes),
    labelsByAxis=scaleLabelsByAxis[:],
    skipEntries=scaleSkipEntries[:],
    systNameReplace=systNameReplaceVec[:],
    baseName=f"QCDscale{'W' if args.wlike else 'Z'}_",
    passToFakes=False if args.noQCDscaleFakes else passSystToFakes,
)

if args.qcdScale == "byCharge":
    scale_action = syst_tools.scale_helicity_hist_to_variations
    scaleActionArgs = {}
    scaleGroupName += "ByChargeV"
    scaleSystAxes.insert(0, "chargeVgen")
    scaleLabelsByAxis.insert(0, "genQ")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before
    
if "Pt" in args.qcdScale:
    scaleActionArgs = {"rebinPtV" : args.rebinPtV}
    scaleGroupName += "ByPtV"
    scaleSystAxes.insert(0, "ptVgen")
    scaleLabelsByAxis.insert(0, "PtBin")
    scaleSystAxes.insert(0, "chargeVgen")
    scaleLabelsByAxis.insert(0, "genQ")
    scaleSkipEntries = [(-1, -1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before

    if args.qcdScale == "byHelicityPtAndByPt":
        # note: the following uses a different histogram compared to the one for helicity splitting
        # when integrating on the coefficients for the helicity split version they should give the same results modulo statistical fluctuations
        print("Option --qcdScale byHelicityPtAndByPt was chosen: doing additional qcdScale histogram for tests")
        print(scaleActionArgs)
        print(scaleLabelsByAxis)
        cardTool.addSystematic("qcdScale",
                               actionMap=scale_action_map,
                               actionArgs=copy.deepcopy(scaleActionArgs), # to avoid possible undesired updates below
                               processes=signal_samples_inctau,
                               group=scaleGroupName,
                               systAxes=scaleSystAxes[:],
                               labelsByAxis=scaleLabelsByAxis[:],
                               skipEntries=scaleSkipEntries[:],
                               systNameReplace=systNameReplaceVec,
                               baseName="QCDscaleByPt_",
                               passToFakes=False if args.noQCDscaleFakes else passSystToFakes,
        )

if helicity:
    scale_hist = "qcdScaleByHelicity"
    scaleActionArgs = {"rebinPtV" : args.rebinPtV}
    scaleSystAxes.insert(0, "helicity")
    scaleLabelsByAxis.insert(0, "Coeff")
    if args.qcdScale == "byHelicityCharge":
        # mainly for tests
        scaleActionArgs.update({"sum_axis" : ["ptVgen"]})
        scaleGroupName += "ByHelicityCharge"
        scaleSystAxes.insert(0, "chargeVgen")
        scaleLabelsByAxis.insert(0, "genQ")
        scaleSkipEntries = [(-1, -1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before
    else:
        scaleGroupName += "ByHelicity"
        scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before

print("Inclusive scale", inclusiveScale)
print(scaleActionArgs if not inclusiveScale else None)
print(scaleLabelsByAxis)

cardTool.addSystematic(scale_hist,
    actionMap=scale_action_map,
    actionArgs=scaleActionArgs,
    processes=signal_samples_inctau,
    group=scaleGroupName,
    # splitGroup={f"{scaleGroupName}_coeff{i}" : f".*Coeff{i}" for i in range(9)}, # key is the new group name to make it unique, value is the pattern to filter nuisances
    systAxes=scaleSystAxes,
    labelsByAxis=scaleLabelsByAxis,
    # Exclude all combinations where muR = muF = 1 (nominal) or where
    # they are extreme values (ratio = 4 or 1/4)
    skipEntries=scaleSkipEntries,
    systNameReplace=systNameReplaceVec,
    baseName="QCDscale_",
    passToFakes=False if args.noQCDscaleFakes else passSystToFakes,
)

cardTool.addSystematic("muonScaleSyst", 
    processes=single_vmu_samples,
    group="muonScale",
    baseName="CMS_scale_m_",
    systAxes=["downUpVar", "scaleEtaSlice"],
    labelsByAxis=["downUpVar", "ieta"],
    passToFakes=passSystToFakes,
    scale = args.scaleMuonCorr,
)
cardTool.addSystematic("muonL1PrefireSyst", 
    processes=cardTool.allMCProcesses(),
    group="muonPrefire",
    baseName="CMS_prefire_syst_m",
    systAxes=["downUpVar"],
    labelsByAxis=["downUpVar"],
    passToFakes=passSystToFakes,
)
# TODO: Allow to be appended to previous group ## FIXME: doesn't it do it already?
cardTool.addSystematic("muonL1PrefireStat", 
    processes=cardTool.allMCProcesses(),
    group="muonPrefire",
    baseName="CMS_prefire_stat_m",
    systAxes=["downUpVar", "etaPhiRegion"],
    labelsByAxis=["downUpVar", "etaPhiReg"],
    passToFakes=passSystToFakes,
)

cardTool.addSystematic("ecalL1Prefire", 
    processes=cardTool.allMCProcesses(),
    group="ecalPrefire",
    baseName="CMS_prefire_ecal",
    systAxes=["downUpVar"],
    labelsByAxis=["downUpVar"],
    passToFakes=passSystToFakes,
)

if not args.wlike:
    cardTool.addLnNSystematic("CMS_Fakes", processes=[args.qcdProcessName], size=1.05)
    cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
    cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
else:
    cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15)

cardTool.writeOutput()

