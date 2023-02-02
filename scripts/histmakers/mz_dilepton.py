from utilities import boostHistHelpers as hh,common,output_tools

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation, muon_calibration, muon_selections
import hist
import lz4.frame
import logging
import math
import time
import os

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--muonCorr", type=str, default="massfit", choices=["lbl", "trackfit_only_mctruth", "none", "massfit", "massfit_lbl", "trackfit_only"], 
    help="Type of correction to apply to the muons")
parser.add_argument("--muScaleMag", type=float, default=1e-4, help="Magnitude of dummy muon scale uncertainty")
parser.add_argument("--muScaleBins", type=int, default=1, help="Number of bins for muon scale uncertainty")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--csvars_hist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--uncertainty-axes", type=str, nargs="*", default=["mll",], help="")
parser.add_argument("--finePtBinning", action='store_true', help="Use fine binning for ptll")
parser.add_argument("--bias-calibration", action='store_true', help="Adjust central value by calibration bias hist")
parser.add_argument("--smearing", action='store_true', help="Smear pT such that resolution matches data")

parser = common.set_parser_default(parser, "pt", [44,26.,70.])
parser = common.set_parser_default(parser, "eta", [6,-2.4,2.4])

args = parser.parse_args()

if args.noColorLogger:
    logger = common.setup_base_logger(os.path.basename(__file__), args.debug)
else:
    logger = common.setup_color_logger(os.path.basename(__file__), args.verbose)

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9", base_path=args.data_path)

era = args.era


# available axes for dilepton validation plots
axes = {
    "mll": hist.axis.Regular(60, 60., 120., name = "mll"),
    "yll": hist.axis.Regular(25, -2.5, 2.5, name = "yll"),
    "ptll": hist.axis.Variable(common.ptV_binning if not args.finePtBinning else range(60), name = "ptll"),
    "etaPlus": hist.axis.Regular(int(args.eta[0]), args.eta[1], args.eta[2], name = "etaPlus"),
    "etaMinus": hist.axis.Regular(int(args.eta[0]), args.eta[1], args.eta[2], name = "etaMinus"),
    "ptPlus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptPlus"),
    "ptMinus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptMinus"),
    "cosThetaStarll": hist.axis.Regular(20, -1., 1., name = "cosThetaStarll"),
    "phiStarll": hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phiStarll"),
    "charge": hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge") # categorical axes in python bindings always have an overflow bin, so use a regular
}

for a in args.uncertainty_axes:
    if a not in axes.keys():
        logger.error(f" {a} is not a known axes! Supported axes choices are {list(axes.keys())}")

nominal_cols = args.uncertainty_axes

if args.csvars_hist:
    nominal_cols += ["cosThetaStarll", "phiStarll"]

nominal_cols.append("charge")

nominal_axes = [axes[a] for a in nominal_cols] 

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False
)

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     is_w_like = True) 
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     is_w_like = True)
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)

mc_jpsi_crctn_helper, data_jpsi_crctn_helper = muon_validation.make_jpsi_crctn_helpers(args.muonCorr)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers()

smearing_helper = muon_calibration.make_muon_smearing_helpers() if args.smearing else None

bias_helper = muon_calibration.make_muon_bias_helpers(args.muonCorr, args.smearing) if args.bias_calibration else None

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mumu", met=args.met)


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args.muonCorr, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, nMuons=2, use_trackerMuons=args.trackerMuons, use_isolation=True)

    df = muon_selections.define_trigger_muons(df)

    df = muon_selections.select_z_candidate(df, args.pt[1], args.pt[2])

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "trigMuons")
    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "nonTrigMuons")

    df = muon_selections.apply_triggermatching_muon(df, dataset, "trigMuons_eta0", "trigMuons_phi0")

    df = df.Define("ptll", "ll_mom4.pt()")
    df = df.Define("yll", "ll_mom4.Rapidity()")
    df = df.Define("absYll", "std::fabs(yll)")
    df = df.Define("csSineCosThetaPhill", "trigMuons_charge0 == -1 ? wrem::csSineCosThetaPhi(trigMuons_mom4, nonTrigMuons_mom4) : wrem::csSineCosThetaPhi(nonTrigMuons_mom4, trigMuons_mom4)")
    
    # "renaming" to write out corresponding axis
    df = df.Alias("charge", "trigMuons_charge0")
    df = df.Define("etaPlus", "trigMuons_charge0 == -1 ? nonTrigMuons_eta0 : trigMuons_eta0") 
    df = df.Define("etaMinus", "trigMuons_charge0 == 1 ? nonTrigMuons_eta0 : trigMuons_eta0") 
    df = df.Define("ptPlus", "trigMuons_charge0 == -1 ? nonTrigMuons_pt0 : trigMuons_pt0") 
    df = df.Define("ptMinus", "trigMuons_charge0 == 1 ? nonTrigMuons_pt0 : trigMuons_pt0") 

    df = df.Define("cosThetaStarll", "csSineCosThetaPhill.costheta")
    df = df.Define("phiStarll", "std::atan2(csSineCosThetaPhill.sinphi, csSineCosThetaPhill.cosphi)")

    if not dataset.is_data:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0",
                                                                                      "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_charge0"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        weight_expr = "weight*weight_pu*weight_fullMuonSF_withTrackingReco*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"

        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)
        if isW or isZ:
            df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)
            df = theory_tools.define_scale_tensor(df)
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")

    results.append(df.HistoBoost("weight", [hist.axis.Regular(100, -2, 2)], ["nominal_weight"]))
    results.append(df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"]))

    if not dataset.is_data and not args.onlyMainHistograms:


        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, nominal_axes, nominal_cols, is_w_like=True)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, nominal_axes, nominal_cols)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isW or isZ:

            if args.theory_corr and dataset.name in corr_helpers:
                results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                    corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only))

            scale_axes = [*nominal_axes, axis_ptVgen, axis_chargeVgen]
            scale_cols = [*nominal_cols, "ptVgen", "chargeVgen"]
            syst_tools.add_qcdScale_hist(results, df, scale_axes, scale_cols)
            syst_tools.add_pdf_hists(results, df, dataset.name, nominal_axes, nominal_cols, args.pdfs)


            df = syst_tools.define_mass_weights(df)
            if isZ:
                syst_tools.add_massweights_hist(results, df, nominal_axes, nominal_cols)
                # there is no W backgrounds for the Wlike, make QCD scale histograms only for Z
                # should probably remove the charge here, because the Z only has a single charge and the pt distribution does not depend on which charged lepton is selected
                if not args.skipHelicity:
                    # TODO: Should have consistent order here with the scetlib correction function
                    syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, scale_axes, scale_cols)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                syst_tools.add_muonscale_hist(
                    results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, nominal_axes, nominal_cols,
                    muon_eta="trigMuons_eta0")


    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, "mz_dilepton.pkl.lz4", args)
