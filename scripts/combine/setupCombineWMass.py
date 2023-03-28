#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools,combine_helpers
from wremnants import histselections as sel
from wremnants.datasets.datagroups2016 import Datagroups2016
from utilities import common, logging
import argparse
import os
import pathlib
import hist
import copy
import math

import pdb

scriptdir = f"{pathlib.Path(__file__).parent}"
data_dir = f"{pathlib.Path(__file__).parent}/../../wremnants/data/"

def make_parser(parser=None):
    if not parser:
        parser = common.common_parser_combine()
    parser.add_argument("--fitvar", help="Variable to fit", default="eta-pt")
    parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
    parser.add_argument("--ewUnc", action='store_true', help="Include EW uncertainty")
    parser.add_argument("--pseudoData", type=str, help="Hist to use as pseudodata")
    parser.add_argument("--pseudoDataIdx", type=int, default=0, help="Variation index to use as pseudodata")
    parser.add_argument("--pseudoDataFile", type=str, help="Input file for pseudodata (if it should be read from a different file", default=None)
    parser.add_argument("-x",  "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
    parser.add_argument("-k",  "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
    parser.add_argument("--skipOtherChargeSyst", action="store_true", help="Skip saving histograms and writing nuisance in datacard for systs defined for a given charge but applied on the channel with the other charge")
    parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
    parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing)")
    parser.add_argument("--effStatLumiScale", type=float, default=None, help="Rescale equivalent luminosity for efficiency stat uncertainty by this value (e.g. 10 means ten times more data from tag and probe)")
    parser.add_argument("--binnedScaleFactors", action='store_true', help="Use binned scale factors (different helpers and nuisances)")
    parser.add_argument("--directIsoSFsmoothing", action='store_true', help="If isolation SF were smoothed directly instead of being derived from smooth efficiencies")
    parser.add_argument("--xlim", type=float, nargs=2, default=None, help="Restrict x axis to this range")
    parser.add_argument("--constrainMass", action='store_true', help="Constrain mass parameter in the fit (e.g. for ptll fit)")
    parser.add_argument("--unfold", action='store_true', help="Prepare datacard for unfolding")
    parser.add_argument("--fitXsec", action='store_true', help="Fit signal inclusive cross section")
    
    return parser

def main(args):

    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    # NOTE: args.filterProcGroups and args.excludeProcGroups should in principle not be used together
    #       (because filtering is equivalent to exclude something), however the exclusion is also meant to skip
    #       processes which are defined in the original process dictionary but are not supposed to be (always) run on
    if args.addQCDMC or "QCD" in args.filterProcGroups:
        logger.warning("Adding QCD MC to list of processes for the fit setup")
    else:
        if "QCD" not in args.excludeProcGroups:
            logger.warning("Automatic removal of QCD MC from list of processes. Use --filterProcGroups 'QCD' or --addQCDMC to keep it")
            args.excludeProcGroups.append("QCD")
    filterGroup = args.filterProcGroups if args.filterProcGroups else None
    excludeGroup = args.excludeProcGroups if args.excludeProcGroups else None
    logger.debug(f"Filtering these groups of processes: {args.filterProcGroups}")
    logger.debug(f"Excluding these groups of processes: {args.excludeProcGroups}")
    
    datagroups = Datagroups2016(args.inputFile, excludeGroups=excludeGroup, filterGroups=filterGroup, splitWByCharge=args.unfold)
    
    if args.xlim:
        if len(args.fitvar.split("-")) > 1:
            raise ValueError("Restricting the x axis not supported for 2D hist")
        s = hist.tag.Slicer()
        datagroups.setGlobalAction(lambda h: h[{args.fitvar : s[complex(0, args.xlim[0]):complex(0, args.xlim[1])]}])

    wmass = datagroups.wmass
    wlike = datagroups.wlike

    if wmass:
        name = "WMass"
    elif wlike:
        name = "ZMassWLike"
    else:
        name = "ZMassDilepton"

    tag = name+"_"+args.fitvar.replace("-","_")
    if args.doStatOnly:
        tag += "_statOnly"
    if args.postfix:
        tag += f"_{args.postfix}"

    outfolder = f"{args.outfolder}/{tag}/"
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    templateDir = f"{scriptdir}/Templates/WMass"

    if args.unfold and args.fitXsec:
        raise ValueError("Options --unfolding and --fitXsec are incompatible. Please choose one or the other")
    elif args.fitXsec:
        datagroups.unconstrainedProcesses.append("Wmunu" if wmass else "Zmumu")
    elif args.unfold:
        if not args.constrainMass:
            logger.warning("Unfolding is specified but the mass is treated free floating, to constrain the mass add '--constrainMass'")

        base_procs = ["Wmunu_q0", "Wmunu_q1"] if wmass else ["Zmumu"]

        for base_proc in base_procs:
            datagroups.defineSignalBinsUnfolding(base_proc)

    if args.noHist and args.noStatUncFakes:
        raise ValueError("Option --noHist would override --noStatUncFakes. Please select only one of them")

    # Start to create the CardTool object, customizing everything
    cardTool = CardTool.CardTool(f"{outfolder}/{name}_{{chan}}.txt")
    cardTool.setDatagroups(datagroups)
    logger.debug(f"Making datacards with these processes: {cardTool.getProcesses()}")
    cardTool.setNominalTemplate(f"{templateDir}/main.txt")
    if args.sumChannels:
        cardTool.setChannels(["inclusive"])
    cardTool.setProjectionAxes(args.fitvar.split("-"))
    if args.noHist:
        cardTool.skipHistograms()
    cardTool.setOutfile(os.path.abspath(f"{outfolder}/{name}CombineInput.root"))
    cardTool.setFakeName(args.qcdProcessName)
    cardTool.setSpacing(52)
    if args.noStatUncFakes:
        cardTool.setProcsNoStatUnc(procs=args.qcdProcessName, resetList=False)
    cardTool.setCustomSystForCard(args.excludeNuisances, args.keepNuisances)
    if args.skipOtherChargeSyst:
        cardTool.setSkipOtherChargeSyst()
    if args.pseudoData:
        cardTool.setPseudodata(args.pseudoData, args.pseudoDataIdx)
        if args.pseudoDataFile:
            cardTool.setPseudodataDatagroups(Datagroups2016(args.pseudoDataFile,
                                                excludeProcGroup=excludeGroup,
                                                filterProcGroup=filterGroup))
    if args.lumiScale:
        cardTool.setLumiScale(args.lumiScale)

    logger.info(f"cardTool.allMCProcesses(): {cardTool.allMCProcesses()}")
        
    passSystToFakes = wmass and not args.skipSignalSystOnFakes and args.qcdProcessName not in excludeGroup

    single_v_samples = cardTool.filteredProcesses(lambda x: x[0] in ["W", "Z"])
    single_v_nonsig_samples = cardTool.filteredProcesses(lambda x: x[0] == ("Z" if wmass else "W"))
    single_vmu_samples = list(filter(lambda x: "mu" in x, single_v_samples))
    signal_samples = list(filter(lambda x: x[0] == ("W" if wmass else "Z"), single_vmu_samples))
    signal_samples_inctau = list(filter(lambda x: x[0] == ("W" if wmass else "Z"), single_v_samples))

    allMCprocesses_noQCDMC = [x for x in cardTool.allMCProcesses() if x != "QCD"]
    
    logger.info(f"All MC processes {allMCprocesses_noQCDMC}")
    logger.info(f"Single V samples: {single_v_samples}")
    logger.info(f"Single V no signal samples: {single_v_nonsig_samples}")
    logger.info(f"Signal samples: {signal_samples}")

    if not args.constrainMass:
        # keep mass weights here as first systematic, in case one wants to run stat-uncertainty only with --doStatOnly
        cardTool.addSystematic("massWeight", 
                               processes=signal_samples_inctau,
                               group="massShift",
                               groupFilter=lambda x: x == "massShift100MeV",
                               skipEntries=[(f"^massShift{i}MeV.*",) for i in range(0, 100, 10)]+[("^massShift2p1MeV.*",)],
                               mirror=False,
                               #TODO: Name this
                               noConstraint=True,
                               systAxes=["massShift"],
                               passToFakes=passSystToFakes,
    )
    
    if args.doStatOnly:
        # print a card with only mass weights and a dummy syst
        cardTool.addLnNSystematic("dummy", processes=["Top", "Diboson"] if wmass else ["Other"], size=1.001, group="dummy")
        cardTool.writeOutput(args=args)
        logger.info("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
        quit()
    
    if args.constrainMass:
        # add an uncertainty on the mass, e.g. for ptll fits
        cardTool.addSystematic("massWeight", 
            processes=signal_samples_inctau,
            group="massShift",
            groupFilter=lambda x: x == "massShift2p1MeV",
            skipEntries=[(f"^massShift{i}MeV.*",) for i in range(0, 110, 10)],
            mirror=False,
            systAxes=["massShift"],
            passToFakes=passSystToFakes,
        )

    if wmass:
        cardTool.addSystematic("luminosity",
                               processes=allMCprocesses_noQCDMC,
                               outNames=["lumiDown", "lumiUp"],
                               group="luminosity",
                               systAxes=["downUpVar"],
                               labelsByAxis=["downUpVar"],
                               passToFakes=passSystToFakes)
    else:
        # TOCHECK: no fakes here, most likely
        cardTool.addLnNSystematic("luminosity", processes=allMCprocesses_noQCDMC, size=1.012, group="luminosity")

    if args.ewUnc:
        cardTool.addSystematic(f"horacenloewCorr_unc", 
            processes=single_v_samples,
            mirror=True,
            group="theory_ew",
            systAxes=["systIdx"],
            labelsByAxis=["horacenloewCorr_unc"],
            skipEntries=[(0, -1), (2, -1)],
            passToFakes=passSystToFakes,
        )

    if not args.noEfficiencyUnc:
        chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
        effTypesNoIso = ["reco", "tracking", "idip", "trigger"]
        effStatTypes = [x for x in effTypesNoIso]
        if args.binnedScaleFactors or args.directIsoSFsmoothing:
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
                nameReplace = []# if any(x in name for x in chargeDependentSteps) else [("q0", ""), ("q1", "")]  # this part correlates nuisances between charges
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
            silentCheckOtherCharge = False
            if wmass or (wlike and "trigger" in name):
                # For wmass the histograms associated to the wrong reco charge are present but equal to nominal one.
                # To avoid annoying warnings of this expected behaviour we silent the related check in cardtool.
                # This also has to happen for effStat trigger in the wlike case
                # Note: should be the case also for trigger effSyst, but since it is packed with the other steps we can't
                # silent it only for that from here (in any case we can live with a single unnecessary warning)
                silentCheckOtherCharge = True
                
            cardTool.addSystematic(
                name, 
                mirror=mirror,
                group=groupName,
                systAxes=axes,
                labelsByAxis=axlabels,
                baseName=name+"_",
                processes=allMCprocesses_noQCDMC,
                passToFakes=passSystToFakes,
                systNameReplace=nameReplace,
                scale=scale,
                splitGroup=splitGroupDict,
                silentCheckOtherCharge=silentCheckOtherCharge
            )

    to_fakes = passSystToFakes and not args.noQCDscaleFakes
    combine_helpers.add_pdf_uncertainty(cardTool, single_v_samples, passSystToFakes)
    combine_helpers.add_scale_uncertainty(cardTool, args.minnloScaleUnc, signal_samples_inctau, to_fakes, resum=args.resumUnc)
    # for Z background in W mass case (W background for Wlike is essentially 0, useless to apply QCD scales there)
    if wmass:
        combine_helpers.add_scale_uncertainty(cardTool, "integrated", single_v_nonsig_samples, False, name_append="Z", resum=args.resumUnc)

    msv_config_dict = {
        "smearingWeights":{
            "hist_name": "muonScaleSyst_responseWeights",
            "syst_axes": ["unc", "downUpVar"],
            "syst_axes_labels": ["unc", "downUpVar"]
        },
        "massWeights":{
            "hist_name": "muonScaleSyst",
            "syst_axes": ["downUpVar", "scaleEtaSlice"],
            "syst_axes_labels": ["downUpVar", "ieta"]
        },
        "manualShift":{
            "hist_name": "muonScaleSyst_manualShift",
            "syst_axes": ["downUpVar"],
            "syst_axes_labels": ["downUpVar"]
        }
    }

    # FIXME: remove this once msv from smearing weights is implemented for the Z
    msv_config = msv_config_dict[args.muonScaleVariation] if wmass else msv_config_dict["massWeights"]

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
        processes=allMCprocesses_noQCDMC,
        group="muonPrefire",
        baseName="CMS_prefire_syst_m",
        systAxes=["downUpVar"],
        labelsByAxis=["downUpVar"],
        passToFakes=passSystToFakes,
    )
    cardTool.addSystematic("muonL1PrefireStat", 
        processes=allMCprocesses_noQCDMC,
        group="muonPrefire",
        baseName="CMS_prefire_stat_m_",
        systAxes=["downUpVar", "etaPhiRegion"],
        labelsByAxis=["downUpVar", "etaPhiReg"],
        passToFakes=passSystToFakes,
    )
    cardTool.addSystematic("ecalL1Prefire", 
        processes=allMCprocesses_noQCDMC,
        group="ecalPrefire",
        baseName="CMS_prefire_ecal",
        systAxes=["downUpVar"],
        labelsByAxis=["downUpVar"],
        passToFakes=passSystToFakes,
    )
    if wmass:
        #cardTool.addLnNSystematic("CMS_Fakes", processes=[args.qcdProcessName], size=1.05, group="MultijetBkg")
        cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
        cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)

        ## FIXME 1: with the jet cut removed this syst is probably no longer needed, but one could still consider
        ## it to cover for how much the fake estimate changes when modifying the composition of the QCD region
        ## FIXME 2: it doesn't really make sense to mirror this one since the systematic goes only in one direction
        # cardTool.addSystematic(f"qcdJetPt30", 
        #                        processes=["Fake"],
        #                        mirror=True,
        #                        group="MultijetBkg",
        #                        systAxes=[],
        #                        outNames=["qcdJetPt30Down", "qcdJetPt30Up"],
        #                        passToFakes=passSystToFakes,
        # )
        #
        if "Fake" not in excludeGroup:
            cardTool.addSystematic(f"nominal", # this is the histogram to read
                                   systAxes=[],
                                   processes=["Fake"],
                                   mirror=True,
                                   group="MultijetBkg",
                                   outNames=["mtCorrFakesDown", "mtCorrFakesUp"],
                                   decorrelateByCharge=True,
                                   rename="mtCorrFakes", # this is the name used to identify the syst in the list of systs
                                   action=sel.applyCorrection,
                                   doActionBeforeMirror=True,
                                   actionArgs={"scale": 1.0,
                                               "corrFile" : f"{data_dir}/fakesWmass/fakerateFactorMtBasedCorrection_vsEtaPt.root",
                                               "corrHist": "etaPtCharge_mtCorrection",
                                               "offsetCorr": 1.0,
                                               "createNew": True}
                               # add action to multiply by correction
        )

    else:
        cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15)

    cardTool.writeOutput(args=args)
    logger.info(f"Output stored in {outfolder}")
    
if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    main(args)
    
