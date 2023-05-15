from utilities import boostHistHelpers as hh, common, output_tools, logging, differential

parser,initargs = common.common_parser(True)

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation, muon_calibration, muon_selections, unfolding_tools
import hist
import lz4.frame
import math
import time
import os

parser = common.set_parser_default(parser, "pt", [34, 26, 60])

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
    
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              nanoVersion="v8" if args.v8 else "v9", base_path=args.dataPath)

if args.validateByMassWeights:
    raise NotImplementedError("Validation of muon scale hists. by massWeights are not implemented!")

era = args.era

# custom template binning
template_neta = int(args.eta[0])
template_mineta = args.eta[1]
template_maxeta = args.eta[2]
logger.info(f"Eta binning: {template_neta} bins from {template_mineta} to {template_maxeta}")
template_npt = int(args.pt[0])
template_minpt = args.pt[1]
template_maxpt = args.pt[2]
logger.info(f"Pt binning: {template_npt} bins from {template_minpt} to {template_maxpt}")

# standard regular axes
axis_eta = hist.axis.Regular(template_neta, template_mineta, template_maxeta, name = "eta", overflow=not args.excludeFlow, underflow=not args.excludeFlow)
axis_pt = hist.axis.Regular(template_npt, template_minpt, template_maxpt, name = "pt", overflow=not args.excludeFlow, underflow=not args.excludeFlow)

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

nominal_axes = [axis_eta, axis_pt, axis_charge]

unfolding_axes, unfolding_cols = differential.get_pt_eta_charge_axes(args.genBins, template_minpt, template_maxpt, template_maxeta)

# axes for mT measurement
axis_mt = hist.axis.Regular(200, 0., 200., name = "mt",underflow=False, overflow=True)
axis_eta_mT = hist.axis.Variable([-2.4, 2.4], name = "eta")

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = axis_pt.edges[-1],
                                                                                                                                     is_w_like = True) 
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = axis_pt.edges[-1],
                                                                                                                                     is_w_like = True, directIsoSFsmoothing=args.directIsoSFsmoothing)
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)

mc_jpsi_crctn_helper, data_jpsi_crctn_helper = muon_calibration.make_jpsi_crctn_helpers(args)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)

smearing_helper = muon_calibration.make_muon_smearing_helpers() if args.smearing else None

bias_helper = muon_calibration.make_muon_bias_helpers(args) if args.biasCalibration else None

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theoryCorr)

# recoil initialization
if not args.noRecoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", args, flavor="mumu")


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    unfold = args.unfolding and dataset.name == "ZmumuPostVFP"
    apply_theory_corr = args.theoryCorr and dataset.name in corr_helpers

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    if unfold:
        df = unfolding_tools.define_gen_level(df, args.genLevel, dataset.name, mode="wlike")
        unfolding_tools.add_xnorm_histograms(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols)

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, nMuons=2, use_trackerMuons=args.trackerMuons, use_isolation=True)

    df = muon_selections.define_trigger_muons(df)

    df = muon_selections.select_z_candidate(df, args.pt[1], args.pt[2])

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "trigMuons")
    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "nonTrigMuons")

    df = muon_selections.apply_triggermatching_muon(df, dataset, "trigMuons_eta0", "trigMuons_phi0")

    if not dataset.is_data:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0",
                                                                                      "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_charge0"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        df = df.Define("exp_weight", "weight_pu*weight_fullMuonSF_withTrackingReco*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom")
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")

    results.append(df.HistoBoost("weight", [hist.axis.Regular(100, -2, 2)], ["nominal_weight"], storage=hist.storage.Double()))

    if not args.noRecoil:
        df = df.Define("yZ", "ll_mom4.Rapidity()")
        lep_cols = ["Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]"]
        trg_cols = ["trigMuons_pt0", "trigMuons_phi0", "nonTrigMuons_pt0", "nonTrigMuons_phi0"]
        df = recoilHelper.recoil_Z(df, results, dataset, common.zprocs_recoil, lep_cols, trg_cols) # produces corrected MET as MET_corr_rec_pt/phi
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    #TODO improve this to include muon mass?
    ###########
    # utility plots of transverse mass, with or without recoil corrections
    ###########
    met_vars = ("MET_pt", "MET_phi")
    df = df.Define("transverseMass_uncorr", f"wrem::mt_wlike_nano(trigMuons_pt0, trigMuons_phi0, nonTrigMuons_pt0, nonTrigMuons_phi0, {', '.join(met_vars)})")
    results.append(df.HistoBoost("transverseMass_uncorr", [axis_mt], ["transverseMass_uncorr", "nominal_weight"]))
    ###########
    met_vars = ("MET_corr_rec_pt", "MET_corr_rec_phi")
    df = df.Define("transverseMass", f"wrem::mt_wlike_nano(trigMuons_pt0, trigMuons_phi0, nonTrigMuons_pt0, nonTrigMuons_phi0, {', '.join(met_vars)})")
    results.append(df.HistoBoost("transverseMass", [axis_mt], ["transverseMass", "nominal_weight"]))
    results.append(df.HistoBoost("MET", [hist.axis.Regular(20, 0, 100, name="MET")], ["MET_corr_rec_pt", "nominal_weight"]))
    ###########
    
    df = df.Filter("transverseMass >= 45.") # 40 for Wmass, thus be 45 here (roughly half the boson mass)
    
    nominal_cols = ["trigMuons_eta0", "trigMuons_pt0", "trigMuons_charge0"]

    if unfold:
        axes = [*nominal_axes, *unfolding_axes] 
        cols = [*nominal_cols, *unfolding_cols]
    else:
        axes = nominal_axes
        cols = nominal_cols

    nominal = df.HistoBoost("nominal", axes, [*cols, "nominal_weight"])
    results.append(nominal)

    if not args.noRecoil:
        df = recoilHelper.add_recoil_unc_Z(df, results, dataset, cols, axes, "nominal")

    if not dataset.is_data and not args.onlyMainHistograms:

        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, axes, cols, is_w_like=True)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, axes, cols)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isW or isZ:

            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, axes, cols, for_wmass=False)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                syst_tools.add_muonscale_hist(
                    results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, axes, cols,
                    muon_eta="trigMuons_eta0")

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args, update_name=not args.forceDefaultName)
