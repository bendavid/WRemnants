#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools,combine_helpers
from wremnants import histselections as sel
from wremnants.datasets.datagroups import datagroups2016
from utilities import common
import argparse
import os
import pathlib
import logging
import hist
import copy
import math

scriptdir = f"{pathlib.Path(__file__).parent}"

def make_parser(parser=None):
    if not parser:
        parser = common.common_parser_combine()
    parser.add_argument("--wlike", action='store_true', help="Run W-like analysis of mZ")
    parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
    parser.add_argument("-p", "--pseudoData", type=str, help="Hist to use as pseudodata")
    parser.add_argument("-x",  "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
    parser.add_argument("-k",  "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
    parser.add_argument("--skipOtherChargeSyst", dest="skipOtherChargeSyst" , action="store_true",   help="Skip saving histograms and writing nuisance in datacard for systs defined for a given charge but applied on the channel with the other charge")
    parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
    parser.add_argument("--decorrelateEffStatIsoByCharge", dest="decorrelateEffStatIsoByCharge", action='store_true', help="Don't correlate isolation efficiency uncertanties between the two charges (by default they are correlated). Obsolete option, one should rather use charge dependent efficiencies directly when they exist")
    parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing)")
    parser.add_argument("--effStatLumiScale", type=float, default=None, help="Rescale equivalent luminosity for efficiency stat uncertainty by this value (e.g. 10 means ten times more data from tag and probe)")
    parser.add_argument("--binnedScaleFactors", action='store_true', help="Use binned scale factors (different helpers and nuisances)")
    return parser

def main(args):
    logging.basicConfig()
    base_logger = logging.getLogger("wremnants")
    base_logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
    logger = base_logger.getChild("setupCombineWMass")

    outfolder = "/".join([args.baseDir, args.outfolder])
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    if args.noHist and args.noStatUncFakes:
        raise ValueError("Option --noHist would override --noStatUncFakes. Please select only one of them")

    wlike = args.wlike
    datagroups = datagroups2016(args.inputFile, wlike=wlike)

    templateDir = f"{scriptdir}/Templates/WMass"
    name = "WMass" if not wlike else "ZMassWLike"
    cardTool = CardTool.CardTool(f"{outfolder}/{name}_{{chan}}.txt")
    cardTool.setNominalTemplate(f"{templateDir}/main.txt")
    if args.noHist:
        cardTool.skipHistograms()
    cardTool.setOutfile(os.path.abspath(f"{outfolder}/{name}CombineInput.root"))
    cardTool.setDatagroups(datagroups)
    cardTool.setFakeName(args.qcdProcessName)
    cardTool.setSpacing(52)
    if args.noStatUncFakes:
        cardTool.setProcsNoStatUnc(procs=args.qcdProcessName, resetList=False)
    cardTool.setCustomSystForCard(args.excludeNuisances, args.keepNuisances)
    if args.skipOtherChargeSyst:
        cardTool.setSkipOtherChargeSyst()
    if args.pseudoData:
        cardTool.setPseudodata(args.pseudoData)
    if args.lumiScale:
        cardTool.setLumiScale(args.lumiScale)
        
    passSystToFakes = not (wlike or args.skipSignalSystOnFakes)
        
    single_v_samples = cardTool.filteredProcesses(lambda x: x[0] in ["W", "Z"])
    single_v_nonsig_samples = cardTool.filteredProcesses(lambda x: x[0] == ("W" if wlike else "Z"))
    single_vmu_samples = list(filter(lambda x: "mu" in x, single_v_samples))
    signal_samples = list(filter(lambda x: x[0] == ("Z" if wlike else "W"), single_vmu_samples))
    signal_samples_inctau = list(filter(lambda x: x[0] == ("Z" if wlike else "W"), single_v_samples))

    logger.info(f"All MC processes {cardTool.allMCProcesses()}")
    logger.info(f"Single V samples: {single_v_samples}")
    logger.info(f"Single V no signal samples: {single_v_nonsig_samples}")
    logger.info(f"Signal samples: {signal_samples}")

    if not wlike and "wlike" in args.inputFile:
        logger.error("You appear to be running with a Wlike input file without the wlike flag! This will probably fail!")
    elif "wlike" not in args.inputFile and args.wlike:
        logger.error("You appear to be running with on a non-Wlike input file with the wlike flag! This will probably fail!")
        
    pdfInfo = theory_tools.pdf_info_map("ZmumuPostVFP", args.pdf)
    pdfName = pdfInfo["name"]

    # keep mass weights here as first systematic, in case one wants to run stat-uncertainty only with --doStatOnly
    cardTool.addSystematic("massWeight", 
        processes=signal_samples_inctau,
        outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=wlike),
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
        cardTool.addLnNSystematic("dummy", processes=["Other"] if wlike else ["Top", "Diboson"], size=1.001, group="dummy")
        cardTool.writeOutput()
        logger.info("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
        quit()
        
    if wlike:
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
        chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
        effTypesNoIso = ["reco", "tracking", "idip", "trigger"]
        effStatTypes = [x for x in effTypesNoIso]
        if args.binnedScaleFactors:
            effStatTypes.extend(["iso"])
        else:
            effStatTypes.extend(["iso_effData", "iso_effMC"])
        allEffTnP = [f"effStatTnP_sf_{eff}" for eff in effStatTypes] + ["effSystTnP"]
        for name in allEffTnP:
            if "Syst" in name:
                axes = ["reco-tracking-idip-trigger-iso"]
                axlabels = ["WPSYST"]
                nameReplace = [("WPSYST0", "reco"), ("WPSYST1", "tracking"), ("WPSYST2", "idip"), ("WPSYST3", "trigger"), ("WPSYST4", "iso"), ("effSystTnP", "effSyst")]
                scale = 1.0
                mirror = True
                groupName = "muon_eff_syst"
                splitGroupDict = {f"{groupName}_{x}" : f".*effSyst.*{x}" for x in list(effTypesNoIso + ["iso"])}
                splitGroupDict[groupName] = ".*effSyst.*" # add also the group with everything
            else:
                nameReplace = [] if any(x in name for x in chargeDependentSteps) else [("q0", ""), ("q1", "")]  # this part correlates nuisances between charges
                if args.binnedScaleFactors:
                    axes = ["SF eta", "nPtBins", "SF charge"]
                    axlabels = ["eta", "pt", "q"]
                    mirror = True
                    nameReplace = nameReplace + [("effStatTnP_sf_", "effStatBinned_")]
                else:
                    axes = ["SF eta", "nPtEigenBins", "SF charge", "downUpVar"]
                    axlabels = ["eta", "pt", "q", "downUpVar"]
                    mirror = False
                    nameReplace = nameReplace + [("effStatTnP_sf_", "effStatSmooth_")]
                scale = 1.0
                groupName = "muon_eff_stat"
                splitGroupDict = {f"{groupName}_{x}" : f".*effStat.*{x}" for x in effStatTypes}
                splitGroupDict[groupName] = ".*effStat.*" # add also the group with everything
            if args.effStatLumiScale and "Syst" not in name:
                scale /= math.sqrt(args.effStatLumiScale)

            cardTool.addSystematic(name, 
                mirror=mirror,
                group=groupName,
                systAxes=axes,
                labelsByAxis=axlabels,
                baseName=name+"_",
                processes=cardTool.allMCProcesses(),
                passToFakes=passSystToFakes,
                systNameReplace=nameReplace,
                scale=scale,
                splitGroup=splitGroupDict
            )

    to_fakes = not (wlike or args.noQCDscaleFakes)
    combine_helpers.add_scale_uncertainty(cardTool, args.qcdScale, signal_samples_inctau, to_fakes, pdf=args.pdf, scetlib=args.scetlibUnc)
    # for Z background in W mass case (W background for Wlike is essentially 0, useless to apply QCD scales there)
    if not wlike:
        combine_helpers.add_scale_uncertainty(cardTool, "integrated", single_v_nonsig_samples, False, pdf=args.pdf, name_append="Z", scetlib=args.scetlibUnc)

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
    cardTool.addSystematic("muonL1PrefireStat", 
        processes=cardTool.allMCProcesses(),
        group="muonPrefire",
        baseName="CMS_prefire_stat_m_",
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

    if not wlike:
        cardTool.addLnNSystematic("CMS_Fakes", processes=[args.qcdProcessName], size=1.05, group="MultijetBkg")
        cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
        cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)

        # FIXME: it doesn't really make sense to mirror this one since the systematic goes only in one direction
        cardTool.addSystematic(f"qcdJetPt45", 
                               processes=["Fake"],
                               mirror=True,
                               group="MultijetBkg",
                               systAxes=[],
                               outNames=["qcdJetPt45Down", "qcdJetPt45Up"],
                               passToFakes=passSystToFakes,
        )
    else:
        cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15)

    cardTool.writeOutput()

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    main(args)
