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

scriptdir = f"{pathlib.Path(__file__).parent}"

def make_parser(parser=None):
    if not parser:
        parser = common.common_parser_combine()
    parser.add_argument("--wlike", action='store_true', help="Run W-like analysis of mZ")
    parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
    parser.add_argument("-p", "--pseudoData", type=str, help="Hist to use as pseudodata")
    parser.add_argument("-x",  "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
    parser.add_argument("-k",  "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
    parser.add_argument("--qcdProcessName", dest="qcdProcessName" , type=str, default="Fake",   help="Name for QCD process")
    parser.add_argument("--noStatUncFakes", dest="noStatUncFakes" , action="store_true",   help="Set bin error for QCD background templates to 0, to check MC stat uncertainties for signal only")
    parser.add_argument("--skipOtherChargeSyst", dest="skipOtherChargeSyst" , action="store_true",   help="Skip saving histograms and writing nuisance in datacard for systs defined for a given charge but applied on the channel with the other charge")
    parser.add_argument("--skipSignalSystOnFakes", dest="skipSignalSystOnFakes" , action="store_true", help="Do not propagate signal uncertainties on fakes, mainly for checks.")
    parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
    parser.add_argument("--correlateEffStatIsoByCharge", action='store_true', help="Correlate isolation efficiency uncertanties between the two charges (by default they are decorrelated)")
    parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing")
    parser.add_argument("--muonScaleVariation", choices=["smearing_weights", "massweights"], default="smearing_weights", help="the method with whicht the distributions for the muon scale variations is derived")
    return parser

def main(args):
    logger = common.setup_base_logger('setupCombineWMass', args.debug)

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
        # TODO: Add the mass weights to the tau samples ## FIXME: isn't it done?
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
        for name,num in zip(["effSystTnP", "effStatTnP", "effStatTnP_tracking", "effStatTnP_reco"], [4, 624*4, 384*2, 144*2]):
            ## TODO: this merged implementation for the effstat makes it very cumbersome to do things differently for iso and trigidip!!
            ## the problem is that I would need custom actions inside based on actual nuisance names, which needs to be commanded from outside, and this is not straightforward
            if "Syst" in name:
                axes = ["reco-tracking-idiptrig-iso"]
                axlabels = ["WPSYST"]
                nameReplace = [("WPSYST0", "Reco"), ("WPSYST1", "Tracking"), ("WPSYST2", "IDIPTrig"), ("WPSYST3", "Iso")]
                scale = 1.0
            elif "tracking" in name:
                axes = ["SF eta", "SF pt", "SF charge"]
                axlabels = ["eta", "pt", "q"]
                nameReplace = [("effStatTnP_tracking", "effStatTnP"), ("q0", "Tracking"), ("q1", "Tracking")] # this serves two purposes: it correlates nuisances between charges and add a sensible labels to nuisances
                scale = 1.0
            elif "reco" in name:
                axes = ["SF eta", "SF pt", "SF charge"]
                axlabels = ["eta", "pt", "q"]
                nameReplace = [("effStatTnP_reco", "effStatTnP"), ("q0", "Reco"), ("q1", "Reco")] # this serves two purposes: it correlates nuisances between charges and add a sensible labels to nuisances
                scale = 1.0
            else:
                axes = ["SF eta", "SF pt", "SF charge", "idiptrig-iso"]
                axlabels = ["eta", "pt", "q", "Trig"]
                nameReplace = [("Trig0", "IDIPTrig"), ("q0Trig1", "Iso"), ("q1Trig1", "Iso")] if args.correlateEffStatIsoByCharge else [("Trig0", "IDIPTrig"), ("Trig1", "Iso")] # replace with better names
                scale = 1.0 if args.correlateEffStatIsoByCharge else {".*effStatTnP.*Iso" : "1.414", ".*effStatTnP.*IDIPTrig" : "1.0"} # only for iso, scale up by sqrt(2) when decorrelating between charges and efficiencies were derived inclusively
            cardTool.addSystematic(name, 
                mirror=True,
                group="muon_eff_syst" if "Syst" in name else "muon_eff_stat", # TODO: for now better checking them separately
                systAxes=axes,
                labelsByAxis=axlabels,
                baseName=name+"_",
                processes=cardTool.allMCProcesses(),
                passToFakes=passSystToFakes,
                systNameReplace=nameReplace,
                scale=scale
            )

    to_fakes = not (wlike or args.noQCDscaleFakes)
    combine_helpers.add_scale_uncertainty(cardTool, args.qcdScale, signal_samples_inctau, to_fakes, pdf=args.pdf, scetlib=args.scetlibUnc)
    # for Z background in W mass case (W background for Wlike is essentially 0, useless to apply QCD scales there)
    if not wlike:
        combine_helpers.add_scale_uncertainty(cardTool, "integrated", single_v_nonsig_samples, False, pdf=args.pdf, name_append="Z", scetlib=args.scetlibUnc)

    msv_config_dict = {
        "smearing_weights":{
            "hist_name": "muonScaleSyst_responseWeights",
            "syst_axes": ["downUpVar"],
            "syst_axes_labels": ["downUpVar"]
        },
        "massweights":{
            "hist_name": "muonScaleSyst",
            "syst_axes": ["downUpVar", "scaleEtaSlice"],
            "syst_axes_labels": ["downUpVar", "ieta"]
        }
    }
    msv_config = msv_config_dict[args.muonScaleVariation]

    cardTool.addSystematic(msv_config['hist_name'], 
        processes=single_vmu_samples,
        group="muonScale",
        baseName="CMS_scale_m_",
        systAxes=msv_config['syst_axes'],
        labelsByAxis=msv_config['syst_axes_labels'],
        passToFakes=passSystToFakes,
        scale = args.scaleMuonCorr
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
