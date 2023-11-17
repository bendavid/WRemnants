#!/usr/bin/env python3
from wremnants import CardTool,combine_helpers,combine_theory_helper, HDF5Writer
from wremnants.syst_tools import massWeightNames
from wremnants.datasets.datagroups import Datagroups
from utilities import common, logging, boostHistHelpers as hh
from utilities.io_tools import input_tools
import argparse
import hist
import math, copy
import h5py
import narf.ioutils

def make_parser(parser=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfolder", type=str, default=".", help="Output folder with the root file storing all histograms and datacards for single charge (subfolder WMass or ZMassWLike is created automatically inside)")
    parser.add_argument("-i", "--inputFile", nargs="+", type=str)
    parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4],
                        help="Set verbosity level with logging, the larger the more verbose")
    parser.add_argument("--noColorLogger", action="store_true", help="Do not use logging with colors")
    parser.add_argument("--hdf5", action="store_true", help="Write out datacard in hdf5")
    parser.add_argument("--sparse", action="store_true", help="Write out datacard in sparse mode (only for when using hdf5)")
    parser.add_argument("--excludeProcGroups", type=str, nargs="*", help="Don't run over processes belonging to these groups (only accepts exact group names)", default=["QCD"])
    parser.add_argument("--filterProcGroups", type=str, nargs="*", help="Only run over processes belonging to these groups", default=[])
    parser.add_argument("-x", "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
    parser.add_argument("-k", "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
    parser.add_argument("--absolutePathInCard", action="store_true", help="In the datacard, set Absolute path for the root file where shapes are stored")
    parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
    parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing)")
    parser.add_argument("--qcdProcessName" , type=str, default="Fake", help="Name for QCD process (must be consistent with what is used in datagroups2016.py")
    # setting on the fit behaviour
    parser.add_argument("--realData", action='store_true', help="Store real data in datacards")
    parser.add_argument("--fitvar", nargs="+", help="Variable to fit", default=["eta-pt-charge"])
    parser.add_argument("--rebin", type=int, nargs='*', default=[], help="Rebin axis by this value (default, 1, does nothing)")
    parser.add_argument("--absval", type=int, nargs='*', default=[], help="Take absolute value of axis if 1 (default, 0, does nothing)")
    parser.add_argument("--axlim", type=float, default=[], nargs='*', help="Restrict axis to this range (assumes pairs of values by axis, with trailing axes optional)")
    parser.add_argument("--lumiScale", type=float, default=1.0, help="Rescale equivalent luminosity by this value (e.g. 10 means ten times more data and MC)")
    parser.add_argument("--sumChannels", action='store_true', help="Only use one channel")
    parser.add_argument("--fitXsec", action='store_true', help="Fit signal inclusive cross section")
    parser.add_argument("--fitMassDiff", type=str, default=None, choices=["charge", "eta-sign", "eta-range"], help="Fit an additional POI for the difference in the boson mass")
    parser.add_argument("--fitresult", type=str, default=None ,help="Use data and covariance matrix from fitresult (for making a theory fit)")
    parser.add_argument("--noMCStat", action='store_true', help="Do not include MC stat uncertainty in covariance for theory fit (only when using --fitresult)")
    parser.add_argument("--fakerateAxes", nargs="+", help="Axes for the fakerate binning", default=["eta","pt","charge"])
    parser.add_argument("--ABCD", action="store_true", help="Produce datacard for simultaneous fit of ABCD regions")
    # settings on the nuisances itself
    parser.add_argument("--doStatOnly", action="store_true", default=False, help="Set up fit to get stat-only uncertainty (currently combinetf with -S 0 doesn't work)")
    parser.add_argument("--minnloScaleUnc", choices=["byHelicityPt", "byHelicityPtCharge", "byHelicityCharge", "byPtCharge", "byPt", "byCharge", "integrated", "none"], default="byHelicityPt",
            help="Decorrelation for QCDscale")
    parser.add_argument("--resumUnc", default="tnp", type=str, choices=["scale", "tnp", "none"], help="Include SCETlib uncertainties")
    parser.add_argument("--npUnc", default="Delta_Lambda", type=str, choices=combine_theory_helper.TheoryHelper.valid_np_models, help="Nonperturbative uncertainty model")
    parser.add_argument("--tnpMagnitude", default=1, type=float, help="Variation size for the TNP")
    parser.add_argument("--scaleTNP", default=5, type=float, help="Scale the TNP uncertainties by this factor")
    parser.add_argument("--scalePdf", default=1, type=float, help="Scale the PDF hessian uncertainties by this factor")
    parser.add_argument("--pdfUncFromCorr", action='store_true', help="Take PDF uncertainty from correction hist (Requires having run that correction)")
    parser.add_argument("--ewUnc", type=str, nargs="*", default=["default"], choices=["default","horacenloew", "winhacnloew", "virtual_ew", "virtual_ew_wlike"], help="Include EW uncertainty")
    parser.add_argument("--widthUnc", action='store_true', help="Include uncertainty on W and Z width")
    parser.add_argument("--skipSignalSystOnFakes" , action="store_true", help="Do not propagate signal uncertainties on fakes, mainly for checks.")
    parser.add_argument("--noQCDscaleFakes", action="store_true",   help="Do not apply QCd scale uncertainties on fakes, mainly for debugging")
    parser.add_argument("--addQCDMC", action="store_true", help="Include QCD MC when making datacards (otherwise by default it will always be excluded)")
    parser.add_argument("--muonScaleVariation", choices=["smearingWeights", "massWeights", "manualShift"], default="smearingWeights", help="the method with which the muon scale variation histograms are derived")
    parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
    parser.add_argument("--correlatedNonClosureNuisances", action='store_true', help="get systematics from histograms for the Z non-closure nuisances without decorrelation in eta and pt")
    parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
    parser.add_argument("--effStatLumiScale", type=float, default=None, help="Rescale equivalent luminosity for efficiency stat uncertainty by this value (e.g. 10 means ten times more data from tag and probe)")
    parser.add_argument("--binnedScaleFactors", action='store_true', help="Use binned scale factors (different helpers and nuisances)")
    parser.add_argument("--isoEfficiencySmoothing", action='store_true', help="If isolation SF was derived from smooth efficiencies instead of direct smoothing")
    # pseudodata
    parser.add_argument("--pseudoData", type=str, nargs="+", help="Histograms to use as pseudodata")
    parser.add_argument("--pseudoDataAxes", type=str, nargs="+", default=[None], help="Variation axes to use as pseudodata for each of the histograms")
    parser.add_argument("--pseudoDataIdxs", type=str, nargs="+", default=[None], help="Variation indices to use as pseudodata for each of the histograms")
    parser.add_argument("--pseudoDataFile", type=str, help="Input file for pseudodata (if it should be read from a different file)", default=None)
    parser.add_argument("--pseudoDataProcsRegexp", type=str, default=".*", help="Regular expression for processes taken from pseudodata file (all other processes are automatically got from the nominal file). Data is excluded automatically as usual")
    # unfolding/differential/theory agnostic
    parser.add_argument("--unfolding", action='store_true', help="Prepare datacard for unfolding")
    parser.add_argument("--genAxes", type=str, default=None, nargs="+", help="Specify which gen axis should be used in unfolding, if 'None', use all (inferred from metadata).")
    parser.add_argument("--theoryAgnostic", action='store_true', help="Prepare datacard for theory agnostic analysis, similar to unfolding but different axis and possibly other differences")
    parser.add_argument("--poiAsNoi", action='store_true', help="Experimental option only with --theoryAgnostic or --unfolding, to treat POIs ad NOIs, with a single signal histogram")
    parser.add_argument("--priorNormXsec", type=float, default=1, help="Prior for shape uncertainties on cross sections for theory agnostic or unfolding analysis with POIs as NOIs (1 means 100\%). If negative, it will use shapeNoConstraint in the fit")
    parser.add_argument("--scaleNormXsecHistYields", type=float, default=None, help="Scale yields of histogram with cross sections variations for theory agnostic analysis with POIs as NOIs. Can be used together with --priorNormXsec")
    parser.add_argument("--addTauToSignal", action='store_true', help="Events from the same process but from tau final states are added to the signal")
    # utility options to deal with charge when relevant, mainly for theory agnostic but also unfolding
    parser.add_argument("--recoCharge", type=str, default=["plus", "minus"], nargs="+", choices=["plus", "minus"], help="Specify reco charge to use, default uses both. This is a workaround for unfolding/theory-agnostic fit when running a single reco charge, as gen bins with opposite gen charge have to be filtered out")
    parser.add_argument("--forceRecoChargeAsGen", action="store_true", help="Force gen charge to match reco charge in CardTool, this only works when the reco charge is used to define the channel")
    # TODO: some options that should exist only for a specific case, can implement a subparser to substitute --unfolding and --theoryAgnostic
    # some options are actually in common between them, so an intermediate subparser might be used
    ##parsers = parser.add_subparsers(dest='analysisFitSetup')
    ##theoryAgnosticParser = parsers.add_parser("unfolding", help="Activate unfolding analysis")
    ##theoryAgnosticParser = parsers.add_parser("theoryAgnostic", help="Activate theory agnostic analysis")
    ##theoryAgnosticParser.add_argument("--noPDFandQCDtheorySystOnSignal", action='store_true', help="Removes PDF and theory uncertainties on signal processes with norm uncertainties when using --poiAsNoi")
    #
    # WIP parser.add_argument("--noNormNuisanceOOA", action='store_true', help="Remove normalization uncertainty on out-of-acceptance template bins. Only implemented with --poiAsNoi")
    parser.add_argument("--noPDFandQCDtheorySystOnSignal", action='store_true', help="Removes PDF and theory uncertainties on signal processes with norm uncertainties when using --poiAsNoi")

    return parser

def setup(args, inputFile, fitvar, xnorm=False):

    # NOTE: args.filterProcGroups and args.excludeProcGroups should in principle not be used together
    #       (because filtering is equivalent to exclude something), however the exclusion is also meant to skip
    #       processes which are defined in the original process dictionary but are not supposed to be (always) run on
    if args.addQCDMC or "QCD" in args.filterProcGroups:
        logger.warning("Adding QCD MC to list of processes for the fit setup")
    elif "QCD" not in args.excludeProcGroups:
        logger.warning("Automatic removal of QCD MC from list of processes. Use --filterProcGroups 'QCD' or --addQCDMC to keep it")
        args.excludeProcGroups.append("QCD")
    filterGroup = args.filterProcGroups if args.filterProcGroups else None
    excludeGroup = args.excludeProcGroups if args.excludeProcGroups else None
    if args.ABCD and (excludeGroup is None or "Fake" not in excludeGroup):
        excludeGroup.append("Fake")
    logger.debug(f"Filtering these groups of processes: {args.filterProcGroups}")
    logger.debug(f"Excluding these groups of processes: {args.excludeProcGroups}")

    datagroups = Datagroups(inputFile, excludeGroups=excludeGroup, filterGroups=filterGroup, applySelection= not xnorm and not args.ABCD)

    if not xnorm and (args.axlim or args.rebin or args.absval):
        datagroups.set_rebin_action(fitvar, args.axlim, args.rebin, args.absval)

    wmass = datagroups.mode in ["wmass", "lowpu_w"]
    wlike = datagroups.mode == "wlike"
    lowPU = "lowpu" in datagroups.mode
    # Detect lowpu dilepton
    dilepton = "dilepton" in datagroups.mode or any(x in ["ptll", "mll"] for x in fitvar)

    constrainMass = (dilepton and not "mll" in fitvar) or args.fitXsec 

    if wmass:
        base_group = "Wenu" if datagroups.flavor == "e" else "Wmunu"
    else:
        base_group = "Zee" if datagroups.flavor == "ee" else "Zmumu"

    if args.addTauToSignal:
        # add tau signal processes to signal group
        datagroups.groups[base_group].addMembers(datagroups.groups[base_group.replace("mu","tau")].members)
        datagroups.deleteGroup(base_group.replace("mu","tau"))

    if args.fitXsec:
        datagroups.unconstrainedProcesses.append(base_group)

    if xnorm:
        datagroups.select_xnorm_groups(base_group)
        datagroups.globalAction = None # reset global action in case of rebinning or such
        if not args.unfolding:
            # creating the xnorm model (e.g. for the theory fit)
            if wmass:
                # add gen charge as additional axis
                datagroups.groups[base_group].add_member_axis("qGen", datagroups.results, 
                    member_filters={-1: lambda x: x.name.startswith("Wminus"), 1: lambda x: x.name.startswith("Wplus")}, 
                    hist_filter=lambda x: x.startswith("xnorm"))
                xnorm_axes = ["qGen", *datagroups.gen_axes]
            else:
                xnorm_axes = datagroups.gen_axes[:]
            datagroups.setGenAxes([a for a in xnorm_axes if a not in fitvar])
    
    if args.poiAsNoi:
        constrainMass = False if args.theoryAgnostic else True
        poi_axes = datagroups.gen_axes if args.genAxes is None else args.genAxes
        # remove specified gen axes from set of gen axes in datagroups so that those are integrated over
        datagroups.setGenAxes([a for a in datagroups.gen_axes if a not in poi_axes])

        # FIXME: temporary customization of signal and out-of-acceptance process names for theory agnostic with POI as NOI
        # There might be a better way to do it more homogeneously with the rest.
        if args.theoryAgnostic:
            # Important: don't set the gen axes with datagroups.setGenAxes(args.genAxes) when doing poiAsNoi (to be checked, it is currently done few lines above) 
            constrainMass = False
            hasSeparateOutOfAcceptanceSignal = False
            # check if the out-of-acceptance signal process exists as an independent process
            if any(m.name.startswith("Bkg") for m in datagroups.groups[base_group].members):
                hasSeparateOutOfAcceptanceSignal = True
                if wmass:
                    # out of acceptance contribution
                    datagroups.copyGroup(base_group, f"BkgWmunu", member_filter=lambda x: x.name.startswith("Bkg"))
                    datagroups.groups[base_group].deleteMembers([m for m in datagroups.groups[base_group].members if m.name.startswith("Bkg")])
                else:
                    # out of acceptance contribution
                    datagroups.copyGroup(base_group, f"BkgZmumu", member_filter=lambda x: x.name.startswith("Bkg"))
                    datagroups.groups[base_group].deleteMembers([m for m in datagroups.groups[base_group].members if m.name.startswith("Bkg")])
            # FIXME: at some point we should decide what name to use
            if any(x in args.excludeProcGroups for x in ["BkgWmunu", "outAccWmunu"]) and hasSeparateOutOfAcceptanceSignal:
                datagroups.deleteGroup("BkgWmunu") # remove out of acceptance signal
    elif args.unfolding or args.theoryAgnostic:
        constrainMass = False if args.theoryAgnostic else True
        datagroups.setGenAxes(args.genAxes)

        if wmass:
            # gen level bins, split by charge
            if "minus" in args.recoCharge:
                datagroups.defineSignalBinsUnfolding(base_group, f"W_qGen0", member_filter=lambda x: x.name.startswith("Wminus"))
            if "plus" in args.recoCharge:
                datagroups.defineSignalBinsUnfolding(base_group, f"W_qGen1", member_filter=lambda x: x.name.startswith("Wplus"))
        else:
            datagroups.defineSignalBinsUnfolding(base_group, "Z", member_filter=lambda x: x.name.startswith(base_group))
        
        # out of acceptance contribution
        to_del = [m for m in datagroups.groups[base_group].members if not m.name.startswith("Bkg")]
        if len(datagroups.groups[base_group].members) == len(to_del):
            datagroups.deleteGroup(base_group)
        else:
            datagroups.groups[base_group].deleteMembers(to_del)    


    # Start to create the CardTool object, customizing everything
    cardTool = CardTool.CardTool(xnorm=xnorm, ABCD=wmass and args.ABCD and not xnorm, real_data=args.realData)
    cardTool.setDatagroups(datagroups)
    if args.qcdProcessName:
        cardTool.setFakeName(args.qcdProcessName)
    logger.debug(f"Making datacards with these processes: {cardTool.getProcesses()}")
    if args.absolutePathInCard:
        cardTool.setAbsolutePathShapeInCard()
    cardTool.setFitAxes(fitvar)
    cardTool.setFakerateAxes(args.fakerateAxes)
    if wmass and args.ABCD:
        # In case of ABCD we need to have different fake processes for e and mu to have uncorrelated uncertainties
        cardTool.setFakeName(datagroups.fakeName + (datagroups.flavor if datagroups.flavor else "")) 
        cardTool.unroll=True
    
    if args.sumChannels or xnorm or dilepton or (wmass and args.ABCD) or "charge" not in fitvar:
        cardTool.setWriteByCharge(False)
    else:
        cardTool.setChannels(args.recoCharge)
        if args.forceRecoChargeAsGen:
            cardTool.setExcludeProcessForChannel("plus", ".*qGen0")
            cardTool.setExcludeProcessForChannel("minus", ".*qGen1")
    
    if xnorm:
        histName = "xnorm"
        cardTool.setHistName(histName)
        cardTool.setNominalName(histName)
    else:
        cardTool.setHistName(args.baseName)
        cardTool.setNominalName(args.baseName)
        
    # define sumGroups for integrated cross section
    if args.unfolding and not args.poiAsNoi:
        # TODO: make this less hardcoded to filter the charge (if the charge is not present this will duplicate things)
        if wmass:
            if "plus" in args.recoCharge:
                cardTool.addPOISumGroups(genCharge="qGen1")
            if "minus" in args.recoCharge:
                cardTool.addPOISumGroups(genCharge="qGen0")
        else:
            cardTool.addPOISumGroups()

    if args.noHist:
        cardTool.skipHistograms()
    cardTool.setSpacing(28)
    cardTool.setCustomSystForCard(args.excludeNuisances, args.keepNuisances)
    if args.pseudoData:
        cardTool.setPseudodata(args.pseudoData, args.pseudoDataAxes, args.pseudoDataIdxs, args.pseudoDataProcsRegexp)
        if args.pseudoDataFile:
            # FIXME: should make sure to apply the same customizations as for the nominal datagroups so far
            pseudodataGroups = Datagroups(args.pseudoDataFile, excludeGroups=excludeGroup, filterGroups=filterGroup, applySelection= not xnorm and not args.ABCD)
            cardTool.setPseudodataDatagroups(pseudodataGroups)
    cardTool.setLumiScale(args.lumiScale)

    if not args.theoryAgnostic:
        logger.info(f"cardTool.allMCProcesses(): {cardTool.allMCProcesses()}")
        
    passSystToFakes = wmass and not args.skipSignalSystOnFakes and args.qcdProcessName not in excludeGroup and (filterGroup == None or args.qcdProcessName in filterGroup) and not xnorm

    # TODO: move to a common place if it is  useful
    def assertSample(name, startsWith=["W", "Z"], excludeMatch=[]):
        return any(name.startswith(init) for init in startsWith) and all(excl not in name for excl in excludeMatch)

    dibosonMatch = ["WW", "WZ", "ZZ"] 
    WMatch = ["W", "BkgW"] # TODO: the name of out-of-acceptance might be changed at some point, maybe to WmunuOutAcc, so W will match it as well (and can exclude it using "OutAcc" if needed)
    ZMatch = ["Z", "BkgZ"]
    signalMatch = WMatch if wmass else ZMatch

    cardTool.addProcessGroup("single_v_samples", lambda x: assertSample(x, startsWith=[*WMatch, *ZMatch], excludeMatch=dibosonMatch))
    if wmass:
        cardTool.addProcessGroup("w_samples", lambda x: assertSample(x, startsWith=WMatch, excludeMatch=dibosonMatch))
        cardTool.addProcessGroup("wtau_samples", lambda x: assertSample(x, startsWith=["Wtaunu"]))
        if not xnorm:
            cardTool.addProcessGroup("single_v_nonsig_samples", lambda x: assertSample(x, startsWith=ZMatch, excludeMatch=dibosonMatch))
    cardTool.addProcessGroup("single_vmu_samples",    lambda x: assertSample(x, startsWith=[*WMatch, *ZMatch], excludeMatch=[*dibosonMatch, "tau"]))
    cardTool.addProcessGroup("signal_samples",        lambda x: assertSample(x, startsWith=signalMatch,        excludeMatch=[*dibosonMatch, "tau"]))
    cardTool.addProcessGroup("signal_samples_inctau", lambda x: assertSample(x, startsWith=signalMatch,        excludeMatch=[*dibosonMatch]))
    cardTool.addProcessGroup("MCnoQCD", lambda x: x not in ["QCD", "Data"])
    # FIXME/FOLLOWUP: the following groups may actually not exclude the OOA when it is not defined as an independent process with specific name
    cardTool.addProcessGroup("signal_samples_noOutAcc",        lambda x: assertSample(x, startsWith=["W" if wmass else "Z"], excludeMatch=[*dibosonMatch, "tau"]))
    cardTool.addProcessGroup("signal_samples_inctau_noOutAcc", lambda x: assertSample(x, startsWith=["W" if wmass else "Z"], excludeMatch=[*dibosonMatch]))

    if not (args.theoryAgnostic or args.unfolding) :
        logger.info(f"All MC processes {cardTool.procGroups['MCnoQCD']}")
        logger.info(f"Single V samples: {cardTool.procGroups['single_v_samples']}")
        if wmass and not xnorm:
            logger.info(f"Single V no signal samples: {cardTool.procGroups['single_v_nonsig_samples']}")
        logger.info(f"Signal samples: {cardTool.procGroups['signal_samples']}")

    signal_samples_forMass = ["signal_samples_inctau"]
    if args.theoryAgnostic and not args.poiAsNoi:
        logger.error("Temporarily not using mass weights for Wtaunu. Please update when possible")
        signal_samples_forMass = ["signal_samples"]

    label = 'W' if wmass else 'Z'
    if not (args.doStatOnly and constrainMass):
        cardTool.addSystematic(f"massWeight{label}",
                            processes=signal_samples_forMass,
                            group=f"massShift",
                            noi=not constrainMass,
                            skipEntries=massWeightNames(proc=label, exclude=100),
                            mirror=False,
                            noConstraint=not constrainMass,
                            systAxes=["massShift"],
                            passToFakes=passSystToFakes,
        )

        if args.fitMassDiff:
            suffix = "".join([a.capitalize() for a in args.fitMassDiff.split("-")])
            mass_diff_args = dict(
                name=f"massWeight{label}",
                processes=signal_samples_forMass,
                rename=f"massDiff{suffix}{label}",
                group=f"massDiff{label}",
                systNameReplace=[("Shift",f"Diff{suffix}")],
                skipEntries=massWeightNames(proc=label, exclude=50),
                noi=not constrainMass,
                noConstraint=not constrainMass,
                mirror=False,
                systAxes=["massShift"],
                passToFakes=passSystToFakes,
            )
            if args.fitMassDiff == "charge":
                cardTool.addSystematic(**mass_diff_args,
                                    # # on gen level based on the sample, only possible for mW
                                    # actionMap={m.name: (lambda h, swap=swap_bins: swap(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown")) 
                                    #     for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members if "minus" in m.name},
                                    # on reco level based on reco charge
                                    actionMap={m.name: (lambda h: 
                                            hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "charge", 0)) 
                                        for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
                )
            elif args.fitMassDiff == "eta-sign":
                cardTool.addSystematic(**mass_diff_args, 
                                    actionMap={m.name: (lambda h: 
                                            hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "eta", hist.tag.Slicer()[0:complex(0,0):]))
                                        for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
                )
            elif args.fitMassDiff == "eta-range":
                cardTool.addSystematic(**mass_diff_args, 
                                    actionMap={m.name: (lambda h: 
                                            hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "eta", hist.tag.Slicer()[complex(0,-0.9):complex(0,0.9):]))
                                        for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
                )

    # this appears within doStatOnly because technically these nuisances should be part of it
    if args.poiAsNoi:
        noi_args = dict(
            group=f"normXsec{label}",
            scale=1 if args.priorNormXsec < 0 else args.priorNormXsec, # histogram represents an (args.priorNormXsec*100)% prior
            mirror=False,
            passToFakes=passSystToFakes,
        )
        if args.theoryAgnostic:
            # open file with theory bands
            with h5py.File("theoryband_variations.hdf5", "r") as ff:
                scale_hists = narf.ioutils.pickle_load_h5py(ff["theorybands"])
            print(scale_hists)
            nuisanceBaseName = f"norm{label}CHANNEL_"
            # First do in acceptance bins, then OOA later (for OOA we need to group bins into macro regions)
            cardTool.addSystematic("yieldsTheoryAgnostic",
                                   **noi_args,
                                   systAxes=poi_axes, 
                                   processes=["signal_samples"],
                                   baseName=nuisanceBaseName,
                                   noConstraint=True if args.priorNormXsec < 0 else False,
                                   #customizeNuisanceAttributes={".*AngCoeff4" : {"scale" : 1, "shapeType": "shapeNoConstraint"}},
                                   labelsByAxis=["PtVBin", "YVBin", "AngCoeff"],
                                   systAxesFlow=[], # only bins in acceptance in this call
                                   skipEntries=[{"helicitySig" : [6,7,8]}], # removing last three indices out of 9 (0,1,...,7,8) corresponding to A5,6,7
                                   actionMap={
                                    #    m.name: (lambda h, scale_hist=scale_hists[m.name]: print(h,scale_hist)) for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members},
                                       m.name: (lambda h, scale_hist=scale_hists[m.name]: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}], hh.multiplyHists(h, scale_hist,allowBroadcast=False,flow=False), scale2=args.scaleNormXsecHistYields)) for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members},
                                   )
            # now OOA
            nuisanceBaseNameOOA = f"{nuisanceBaseName}OOA_"
            # TODO: implement a loop to generalize it
            #
            # ptV OOA, yV in acceptance, integrate helicities 
            cardTool.addSystematic("yieldsTheoryAgnostic",
                                   rename=f"yieldsTheoryAgnostic_OOA_ptV",
                                   **noi_args,
                                   processes=["signal_samples"],
                                   baseName=nuisanceBaseNameOOA,
                                   noConstraint=True if args.priorNormXsec < 0 else False,
                                   systAxes=["ptVgenSig"],
                                   labelsByAxis=["PtVBin"],
                                   systAxesFlow=["ptVgenSig"], # this can activate nuisances on overflow bins, mainly just ptV and yV since the helicity axis has no overflow bins
                                   actionMap={
                                       m.name: (lambda h: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}],
                                                                      h[{"ptVgenSig": hist.tag.Slicer()[hist.overflow:],
                                                                         "absYVgenSig": hist.tag.Slicer()[0:h.axes["absYVgenSig"].size:hist.sum],
                                                                         "helicitySig": hist.tag.Slicer()[::hist.sum]}],
                                                                      scale2=args.scaleNormXsecHistYields)
                                                ) for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members
                                   },
                                   )
            # ptV in acceptance, yV OOA, integrate helicities
            cardTool.addSystematic("yieldsTheoryAgnostic",
                                   rename=f"yieldsTheoryAgnostic_OOA_yV",
                                   **noi_args,
                                   processes=["signal_samples"],
                                   baseName=nuisanceBaseNameOOA,
                                   noConstraint=True if args.priorNormXsec < 0 else False,
                                   systAxes=["absYVgenSig"],
                                   labelsByAxis=["YVBin"],
                                   systAxesFlow=["absYVgenSig"], # this can activate nuisances on overflow bins, mainly just ptV and yV since the helicity axis has no overflow bins
                                   actionMap={
                                       m.name: (lambda h: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}],
                                                                      h[{"ptVgenSig": hist.tag.Slicer()[0:h.axes["ptVgenSig"].size:hist.sum],
                                                                         "absYVgenSig": hist.tag.Slicer()[hist.overflow:],
                                                                         "helicitySig": hist.tag.Slicer()[::hist.sum]}],
                                                                      scale2=args.scaleNormXsecHistYields)
                                                ) for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members
                                   },
                                   )
            # ptV OOA and yV OOA, integrate helicities
            cardTool.addSystematic("yieldsTheoryAgnostic",
                                   rename=f"yieldsTheoryAgnostic_OOA_ptVyV",
                                   **noi_args,
                                   processes=["signal_samples"],
                                   baseName=nuisanceBaseNameOOA,
                                   noConstraint=True if args.priorNormXsec < 0 else False,
                                   systAxes=["ptVgenSig", "absYVgenSig"],
                                   labelsByAxis=["PtVBin", "YVBin"],
                                   systAxesFlow=["ptVgenSig", "absYVgenSig"], # this can activate nuisances on overflow bins, mainly just ptV and yV since the helicity axis has no overflow bins
                                   actionMap={
                                       m.name: (lambda h: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}],
                                                                      h[{"ptVgenSig": hist.tag.Slicer()[hist.overflow:],
                                                                         "absYVgenSig": hist.tag.Slicer()[hist.overflow:],
                                                                         "helicitySig": hist.tag.Slicer()[::hist.sum]}],
                                                                      scale2=args.scaleNormXsecHistYields)
                                                ) for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members
                                   },
                                   )

        elif args.unfolding:
            noi_args.update(dict(
                name=f"yieldsUnfolding",
                systAxes=poi_axes,
                processes=["signal_samples"],
                noConstraint=True,
                noi=True,
                systAxesFlow=[a for a in poi_axes if a in ["ptGen"]], # use underflow/overflow bins for ptGen
                labelsByAxis=poi_axes,
            ))
            if wmass:
                # add two sets of systematics, one for each charge
                cardTool.addSystematic(**noi_args,
                    rename=f"noiWminus",
                    baseName=f"W_qGen0",
                    actionMap={
                        m.name: (lambda h: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}], h, scale2=args.scaleNormXsecHistYields))
                            if "minus" in m.name else (lambda h: h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}])
                        for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members},
                )
                cardTool.addSystematic(**noi_args,
                    rename=f"noiWplus",
                    baseName=f"W_qGen1",
                    actionMap={
                        m.name: (lambda h: hh.addHists(h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}], h, scale2=args.scaleNormXsecHistYields))
                            if "plus" in m.name else (lambda h: h[{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}])
                        for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members},
                )
            else:
                cardTool.addSystematic(**noi_args,
                    baseName=f"{label}_",
                    actionMap={
                        m.name: (lambda h: hh.addHists(
                            h[{**{ax: hist.tag.Slicer()[::hist.sum] for ax in poi_axes}, "acceptance": hist.tag.Slicer()[::hist.sum]}], h[{"acceptance":True}], scale2=args.scaleNormXsecHistYields))
                        for g in cardTool.procGroups["signal_samples"] for m in cardTool.datagroups.groups[g].members},
                )

    if args.doStatOnly:
        # print a card with only mass weights
        logger.info("Using option --doStatOnly: the card was created with only mass nuisance parameter")
        return cardTool

    if wmass and not xnorm:
        cardTool.addSystematic(f"massWeightZ",
                                processes=['single_v_nonsig_samples'],
                                group=f"massShift",
                                skipEntries=massWeightNames(proc="Z", exclude=2.1),
                                mirror=False,
                                noConstraint=False,
                                systAxes=["massShift"],
                                passToFakes=passSystToFakes,
        )

    if args.widthUnc:
        widthSkipZ = [("widthZ2p49333GeV",), ("widthZ2p49493GeV",), ("widthZ2p4952GeV",)] 
        widthSkipW = [("widthW2p09053GeV",), ("widthW2p09173GeV",), ("widthW2p085GeV",)]
        if wmass and not xnorm:
            cardTool.addSystematic(f"widthWeightZ",
                                    processes=["single_v_nonsig_samples"],
                                    group=f"widthZ",
                                    skipEntries=widthSkipZ[:],
                                    mirror=False,
                                    systAxes=["width"],
                                    passToFakes=passSystToFakes,
            )
        cardTool.addSystematic(f"widthWeight{label}",
                                processes=["signal_samples_inctau"],
                                skipEntries=widthSkipZ[:] if label=="Z" else widthSkipW[:],
                                group=f"width{label}",
                                mirror=False,
                                #TODO: Name this
                                systAxes=["width"],
                                passToFakes=passSystToFakes,
        )


    ewUncs = args.ewUnc
    if "default" in ewUncs:
        # set default EW uncertainty depending on the analysis type
        if wlike:
            ewUnc = "virtual_ew_wlike"
        elif dilepton:
            ewUnc = "virtual_ew"
        else:
            ewUnc = "horacenloew"
        ewUncs = [ewUnc if u=="default" else u for u in ewUncs]

    for ewUnc in ewUncs:
        if datagroups.flavor == "e":
            logger.warning("EW uncertainties are not implemented for electrons, proceed w/o EW uncertainty")
            continue
        if ewUnc=="winhacnloew" and not wmass:
            logger.warning("Winhac is not implemented for any other process than W, proceed w/o winhac EW uncertainty")
            continue
        elif ewUnc.startswith("virtual_ew") and wmass:
            logger.warning("Virtual EW corrections are not implemented for any other process than Z, uncertainty only applied to this background")

        cardTool.addSystematic(f"{ewUnc}Corr", 
            processes=['w_samples'] if ewUnc=="winhacnloew" else ['single_v_samples'],
            mirror=True,
            group="theory_ew",
            systAxes=["systIdx"],
            labelsByAxis=[f"{ewUnc}Corr"],
            scale=2,
            skipEntries=[(1, -1), (2, -1)] if ewUnc.startswith("virtual_ew") else [(0, -1), (2, -1)],
            passToFakes=passSystToFakes,
        )

    to_fakes = passSystToFakes and not args.noQCDscaleFakes and not xnorm
    
    theory_helper = combine_theory_helper.TheoryHelper(cardTool, hasNonsigSamples=(wmass and not xnorm))
    theory_helper.configure(resumUnc=args.resumUnc, 
        propagate_to_fakes=to_fakes,
        np_model=args.npUnc,
        tnp_magnitude=args.tnpMagnitude,
        tnp_scale = args.scaleTNP,
        mirror_tnp=True,
        pdf_from_corr=args.pdfUncFromCorr,
        scale_pdf_unc=args.scalePdf,
        minnlo_unc=args.minnloScaleUnc,
    )

    theorySystSamples = ["signal_samples_inctau", "single_v_nonsig_samples"]
    if xnorm:
        theorySystSamples = ["signal_samples"]
    if args.noPDFandQCDtheorySystOnSignal:
        theorySystSamples = ["wtau_samples", "single_v_nonsig_samples"]

    theory_helper.add_all_theory_unc(theorySystSamples, skipFromSignal=args.noPDFandQCDtheorySystOnSignal)

    if xnorm or datagroups.mode == "vgen":
        return cardTool

    # Below: experimental uncertainties

    if wmass:
        #cardTool.addLnNSystematic("CMS_Fakes", processes=[args.qcdProcessName], size=1.05, group="MultijetBkg")
        cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06, group="CMS_background")
        cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16, group="CMS_background")
        cardTool.addSystematic("luminosity",
                                processes=['MCnoQCD'],
                                outNames=["lumiDown", "lumiUp"],
                                group="luminosity",
                                systAxes=["downUpVar"],
                                labelsByAxis=["downUpVar"],
                                passToFakes=passSystToFakes)
    else:
        cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15, group="CMS_background")
        cardTool.addLnNSystematic("luminosity", processes=['MCnoQCD'], size=1.017 if lowPU else 1.012, group="luminosity")

    if not args.noEfficiencyUnc:

        if not lowPU:

            chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
            effTypesNoIso = ["reco", "tracking", "idip", "trigger"]
            effStatTypes = [x for x in effTypesNoIso]
            if args.binnedScaleFactors or not args.isoEfficiencySmoothing:
                effStatTypes.extend(["iso"])
            else:
                effStatTypes.extend(["iso_effData", "iso_effMC"])
            allEffTnP = [f"effStatTnP_sf_{eff}" for eff in effStatTypes] + ["effSystTnP"]
            for name in allEffTnP:
                if "Syst" in name:
                    axes = ["reco-tracking-idip-trigger-iso", "n_syst_variations"]
                    axlabels = ["WPSYST", "_etaDecorr"]
                    nameReplace = [("WPSYST0", "reco"), ("WPSYST1", "tracking"), ("WPSYST2", "idip"), ("WPSYST3", "trigger"), ("WPSYST4", "iso"), ("effSystTnP", "effSyst"), ("etaDecorr0", "fullyCorr") ]
                    scale = 1.0
                    mirror = True
                    mirrorDownVarEqualToNomi=False
                    groupName = "muon_eff_syst"
                    splitGroupDict = {f"{groupName}_{x}" : f".*effSyst.*{x}" for x in list(effTypesNoIso + ["iso"])}
                    # decorrDictEff = {                        
                    #     "x" : {
                    #         "label" : "eta",
                    #         "edges": [round(-2.4+i*0.1,1) for i in range(49)]
                    #     }
                    # }

                else:
                    nameReplace = [] if any(x in name for x in chargeDependentSteps) else [("q0", "qall")] # for iso change the tag id with another sensible label
                    mirror = True
                    mirrorDownVarEqualToNomi=False
                    if args.binnedScaleFactors:
                        axes = ["SF eta", "nPtBins", "SF charge"]
                    else:
                        axes = ["SF eta", "nPtEigenBins", "SF charge"]
                    axlabels = ["eta", "pt", "q"]
                    nameReplace = nameReplace + [("effStatTnP_sf_", "effStat_")]           
                    scale = 1.0
                    groupName = "muon_eff_stat"
                    splitGroupDict = {f"{groupName}_{x}" : f".*effStat.*{x}" for x in effStatTypes}
                if args.effStatLumiScale and "Syst" not in name:
                    scale /= math.sqrt(args.effStatLumiScale)

                cardTool.addSystematic(
                    name, 
                    mirror=mirror,
                    mirrorDownVarEqualToNomi=mirrorDownVarEqualToNomi,
                    group=groupName,
                    systAxes=axes,
                    labelsByAxis=axlabels,
                    baseName=name+"_",
                    processes=['MCnoQCD'],
                    passToFakes=passSystToFakes,
                    systNameReplace=nameReplace,
                    scale=scale,
                    splitGroup=splitGroupDict,
                    decorrelateByBin = {}
                )
                # if "Syst" in name and decorrDictEff != {}:
                #     # add fully correlated version again
                #     cardTool.addSystematic(
                #         name,
                #         rename=f"{name}_EtaDecorr",
                #         mirror=mirror,
                #         mirrorDownVarEqualToNomi=mirrorDownVarEqualToNomi,
                #         group=groupName,
                #         systAxes=axes,
                #         labelsByAxis=axlabels,
                #         baseName=name+"_",
                #         processes=['MCnoQCD'],
                #         passToFakes=passSystToFakes,
                #         systNameReplace=nameReplace,
                #         scale=scale,
                #         splitGroup=splitGroupDict,
                #         decorrelateByBin = decorrDictEff
                #     )
        else:
            if datagroups.flavor in ["mu", "mumu"]:
                lepEffs = ["muSF_HLT_DATA_stat", "muSF_HLT_DATA_syst", "muSF_HLT_MC_stat", "muSF_HLT_MC_syst", "muSF_ISO_stat", "muSF_ISO_DATA_syst", "muSF_ISO_MC_syst", "muSF_IDIP_stat", "muSF_IDIP_DATA_syst", "muSF_IDIP_MC_syst"]
            else:
                lepEffs = [] # ["elSF_HLT_syst", "elSF_IDISO_stat"]

            for lepEff in lepEffs:
                cardTool.addSystematic(lepEff,
                    processes=cardTool.allMCProcesses(),
                    mirror = True,
                    group="CMS_lepton_eff",
                    baseName=lepEff,
                    systAxes = ["tensor_axis_0"],
                    labelsByAxis = [""], 
                )

    if (wmass or wlike) and input_tools.args_from_metadata(cardTool, "recoilUnc"):
        combine_helpers.add_recoil_uncertainty(cardTool, ["signal_samples"],
            passSystToFakes=passSystToFakes,
            flavor=datagroups.flavor if datagroups.flavor else "mu",
            pu_type="lowPU" if lowPU else "highPU")

    if lowPU:
        if datagroups.flavor in ["e", "ee"]:
            # disable, prefiring for muons currently broken? (fit fails)
            cardTool.addSystematic("prefireCorr",
                processes=cardTool.allMCProcesses(),
                mirror = False,
                group="CMS_prefire17",
                baseName="CMS_prefire17",
                systAxes = ["downUpVar"],
                labelsByAxis = ["downUpVar"], 
            )

        return cardTool

    # Below: all that is highPU specific

    # msv_config_dict = {
    #     "smearingWeights":{
    #         "hist_name": "muonScaleSyst_responseWeights",
    #         "syst_axes": ["unc", "downUpVar"],
    #         "syst_axes_labels": ["unc", "downUpVar"]
    #     },
    #     "massWeights":{
    #         "hist_name": "muonScaleSyst",
    #         "syst_axes": ["downUpVar", "scaleEtaSlice"],
    #         "syst_axes_labels": ["downUpVar", "ieta"]
    #     },
    #     "manualShift":{
    #         "hist_name": "muonScaleSyst_manualShift",
    #         "syst_axes": ["downUpVar"],
    #         "syst_axes_labels": ["downUpVar"]
    #     }
    # }

    # msv_config = msv_config_dict[args.muonScaleVariation]

    # cardTool.addSystematic(msv_config['hist_name'], 
    #     processes=['single_v_samples' if wmass else 'single_vmu_samples'],
    #     group="muonCalibration",
    #     baseName="CMS_scale_m_",
    #     systAxes=msv_config['syst_axes'],
    #     labelsByAxis=msv_config['syst_axes_labels'],
    #     passToFakes=passSystToFakes,
    #     scale = args.scaleMuonCorr,
    # )
    cardTool.addSystematic("muonL1PrefireSyst", 
        processes=['MCnoQCD'],
        group="muonPrefire",
        splitGroup = {f"prefire" : f".*"},
        baseName="CMS_prefire_syst_m",
        systAxes=["downUpVar"],
        labelsByAxis=["downUpVar"],
        passToFakes=passSystToFakes,
    )
    cardTool.addSystematic("muonL1PrefireStat", 
        processes=['MCnoQCD'],
        group="muonPrefire",
        splitGroup = {f"prefire" : f".*"},
        baseName="CMS_prefire_stat_m_",
        systAxes=["downUpVar", "etaPhiRegion"],
        labelsByAxis=["downUpVar", "etaPhiReg"],
        passToFakes=passSystToFakes,
    )
    cardTool.addSystematic("ecalL1Prefire", 
        processes=['MCnoQCD'],
        group="ecalPrefire",
        splitGroup = {f"prefire" : f".*"},
        baseName="CMS_prefire_ecal",
        systAxes=["downUpVar"],
        labelsByAxis=["downUpVar"],
        passToFakes=passSystToFakes,
    )

    non_closure_scheme = input_tools.args_from_metadata(cardTool, "nonClosureScheme")
    correlated_non_closure = input_tools.args_from_metadata(cardTool, "correlatedNonClosureNP")
    if non_closure_scheme in ["A-M-separated", "A-only"]:
        cardTool.addSystematic("Z_non_closure_parametrized_A", 
            processes=['single_v_samples'],
            group="nonClosure",
            splitGroup={f"muonCalibration" : f".*"},
            baseName="Z_nonClosure_parametrized_A_",
            systAxes=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            labelsByAxis=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            passToFakes=passSystToFakes
        )
    if non_closure_scheme in ["A-M-separated", "M-only", "binned-plus-M"]:
        cardTool.addSystematic("Z_non_closure_parametrized_M", 
            processes=['single_v_samples'],
            group="nonClosure",
            splitGroup={f"muonCalibration" : f".*"},
            baseName="Z_nonClosure_parametrized_M_",
            systAxes=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            labelsByAxis=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            passToFakes=passSystToFakes
        )            
    if non_closure_scheme == "A-M-combined":
        cardTool.addSystematic("Z_non_closure_parametrized", 
            processes=['single_v_samples'],
            group="nonClosure",
            splitGroup={f"muonCalibration" : f".*"},
            baseName="Z_nonClosure_parametrized_",
            systAxes=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            labelsByAxis=["unc", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            passToFakes=passSystToFakes
        )
    if non_closure_scheme in ["binned", "binned-plus-M"]:
        cardTool.addSystematic("Z_non_closure_binned", 
            processes=['single_v_samples'],
            group="nonClosure",
            splitGroup={f"muonCalibration" : f".*"},
            baseName="Z_nonClosure_binned_",
            systAxes=["unc_ieta", "unc_ipt", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            labelsByAxis=["unc_ieta", "unc_ipt", "downUpVar"] if not correlated_non_closure else ["downUpVar"],
            passToFakes=passSystToFakes
        )
    if not input_tools.args_from_metadata(cardTool, "noSmearing"):
        cardTool.addSystematic("muonResolutionSyst_responseWeights", 
            mirror = True,
            processes=['single_v_samples'],
            group="resolutionCrctn",
            splitGroup={f"muonCalibration" : f".*"},
            baseName="Resolution_correction_",
            systAxes=["smearing_variation"],
            passToFakes=passSystToFakes
        )
       
    
    # Previously we had a QCD uncertainty for the mt dependence on the fakes, see: https://github.com/WMass/WRemnants/blob/f757c2c8137a720403b64d4c83b5463a2b27e80f/scripts/combine/setupCombineWMass.py#L359

    return cardTool

def analysis_label(card_tool):
    analysis_name_map = {
        "wmass" : "WMass",
        "vgen" : "ZGen" if card_tool.getProcesses()[0][0] == "Z" else "WGen",
        "wlike" : "ZMassWLike", 
        "dilepton" : "ZMassDilepton",
        "lowpu_w" : "WMass_lowPU",
        "lowpu_z" : "ZMass_lowPU",
    }

    if card_tool.datagroups.mode not in analysis_name_map:
        raise ValueError(f"Invalid datagroups mode {datagroups.mode}")

    return analysis_name_map[card_tool.datagroups.mode]

def outputFolderName(outfolder, card_tool, doStatOnly, postfix):
    to_join = [analysis_label(card_tool)]+card_tool.fit_axes

    if doStatOnly:
        to_join.append("statOnly")
    if card_tool.datagroups.flavor:
        to_join.append(card_tool.datagroups.flavor)
    if postfix is not None:
        to_join.append(postfix)

    return f"{outfolder}/{'_'.join(to_join)}/"

def main(args, xnorm=False):
    fitvar = args.fitvar[0].split("-") if not xnorm else ["count"]
    cardTool = setup(args, args.inputFile[0], fitvar, xnorm)
    cardTool.setOutput(outputFolderName(args.outfolder, cardTool, args.doStatOnly, args.postfix), analysis_label(cardTool))
    cardTool.writeOutput(args=args, forceNonzero=not args.unfolding, check_systs=not args.unfolding)
    return

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    if args.poiAsNoi and args.theoryAgnostic == args.unfolding:
        raise ValueError("Option --poiAsNoi requires either --theoryAgnostic or --unfolding but not both")    
    if args.noHist and args.noStatUncFakes:
        raise ValueError("Option --noHist would override --noStatUncFakes. Please select only one of them")
    if args.unfolding and args.fitXsec:
        raise ValueError("Options --unfolding and --fitXsec are incompatible. Please choose one or the other")

    if args.theoryAgnostic:
        if args.genAxes is None:
            args.genAxes = ["ptVgenSig", "absYVgenSig", "helicitySig"]
            logger.warning(f"Automatically setting '--genAxes {' '.join(args.genAxes)}' for theory agnostic analysis")
            if args.poiAsNoi:
                logger.warning("This is only needed to properly get the systematic axes")
                
        if not args.poiAsNoi:
            # The following is temporary, just to avoid passing the option explicitly
            logger.warning("For now setting --theoryAgnostic activates --doStatOnly")
            args.doStatOnly = True
    
    if args.hdf5: 
        writer = HDF5Writer.HDF5Writer()

        # loop over all files
        outnames = []
        for i, ifile in enumerate(args.inputFile):
            fitvar = args.fitvar[i].split("-")
            cardTool = setup(args, ifile, fitvar, xnorm=args.fitresult is not None)
            outnames.append( (outputFolderName(args.outfolder, cardTool, args.doStatOnly, args.postfix), analysis_label(cardTool)) )

            writer.add_channel(cardTool)
            if args.unfolding and not args.poiAsNoi:
                cardTool = setup(args, ifile, ["count"], xnorm=True)
                writer.add_channel(cardTool)
        if args.fitresult:
            writer.set_fitresult(args.fitresult, mc_stat=not args.noMCStat)

        if len(outnames) == 1:
            outfile, outfolder = outnames[0]
        else:
            outfile, outfolder = f"{args.outfolder}/Combination{'_statOnly' if args.doStatOnly else ''}_{args.postfix}/", "Combination"
        logger.info(f"Writing HDF5 output to {outfile}")
        writer.write(args, outfile, outfolder)
    else:
        if len(args.inputFile) > 1:
            raise IOError(f"Multiple input files only supported within --hdf5 mode")

        main(args)
        if args.unfolding and not args.poiAsNoi:
            logger.warning("Now running with xnorm = True")
            # in case of unfolding and hdf5, the xnorm histograms are directly written into the hdf5
            main(args, xnorm=True)

    logging.summary()
