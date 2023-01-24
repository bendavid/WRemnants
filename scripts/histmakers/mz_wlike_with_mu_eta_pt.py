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
import pdb

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--muonCorr", type=str, default="trackfit_only", choices=["lbl", "trackfit_only_mctruth", "none", "massfit", "massfit_lbl", "trackfit_only"], 
    help="Type of correction to apply to the muons")
parser.add_argument("--muScaleMag", type=float, default=1e-4, help="Magnitude of dummy muon scale uncertainty")
parser.add_argument("--muScaleBins", type=int, default=1, help="Number of bins for muon scale uncertainty")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--csvars_hist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--uncertainty-hist", type=str, choices=["dilepton", "nominal"], default="nominal", help="Histogram to store uncertainties for")
parser.add_argument("--finePtBinning", action='store_true', help="Use fine binning for ptll")
parser.add_argument("--dileptonIntegrateAxes", type=str, nargs="*", choices=["ptll", "mll", "yll",], 
    default=[], help="Collapse axes in dilepton hist (to avoid overly bloated output)")
args = parser.parse_args()

logging.basicConfig(level=logging.INFO)

f = next((x for x in parser._actions if x.dest == "pt"), None)
if f:
    f.default = [34,26.,60.]

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9", base_path=args.data_path)

era = args.era

muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False
)

# custom template binning
template_neta = int(args.eta[0])
template_mineta = args.eta[1]
template_maxeta = args.eta[2]
logging.info(f"Eta binning: {template_neta} bins from {template_mineta} to {template_maxeta}")
template_npt = int(args.pt[0])
template_minpt = args.pt[1]
template_maxpt = args.pt[2]
logging.info(f"Pt binning: {template_npt} bins from {template_minpt} to {template_maxpt}")

# standard regular axes
axis_eta = hist.axis.Regular(template_neta, template_mineta, template_maxeta, name = "eta")
axis_pt = hist.axis.Regular(template_npt, template_minpt, template_maxpt, name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

nominal_axes = [axis_eta, axis_pt, axis_charge]

# axes for mT measurement
axis_mt = hist.axis.Regular(200, 0., 200., name = "mt",underflow=False, overflow=True)
axis_eta_mT = hist.axis.Variable([-2.4, 2.4], name = "eta")

# extra axes for dilepton validation plots
axis_mll = hist.axis.Regular(60, 60., 120., name = "mll")
axis_yll = hist.axis.Regular(25, -2.5, 2.5, name = "yll")

axis_ptll = hist.axis.Variable(common.ptV_binning if not args.finePtBinning else range(60),
    name = "ptll"
)

if "mll" in args.dileptonIntegrateAxes:
    axis_mll = hist.axis.Variable([0, math.inf], name=axis_mll.name, flow=False)
if "yll" in args.dileptonIntegrateAxes:
    axis_yll = hist.axis.Variable([-math.inf, math.inf], name=axis_yll.name, flow=False)
if "ptll" in args.dileptonIntegrateAxes:
    axis_ptll = hist.axis.Variable([0, math.inf], name=axis_ptll.name, flow=False)

axis_costhetastarll = hist.axis.Regular(20, -1., 1., name = "costhetastarll")
axis_phistarll = hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phistarll")

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logging.info("Using binned scale factors and uncertainties")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = axis_pt.edges[-1],
                                                                                                                                     is_w_like = True) 
else:
    logging.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = axis_pt.edges[-1],
                                                                                                                                     is_w_like = True)
logging.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)

mass_fit = "massfit" in args.muonCorr
if mass_fit:
    mc_corrfile = "calibrationJMC_smeared_v718_nominalLBL.root" if "lbl" in args.muonCorr else "calibrationJMC_smeared_v718_nominal.root"
    data_corrfile = "calibrationJDATA_smeared_v718_LBL.root" if "lbl" in args.muonCorr else "calibrationJDATA_smeared_v718.root"
    jpsi_crctn_MC_helper = muon_validation.make_jpsi_crctn_helper(filepath=f"{common.data_dir}/calibration/{mc_corrfile}")
    jpsi_crctn_data_helper = muon_validation.make_jpsi_crctn_helper(filepath=f"{common.data_dir}/calibration/{data_corrfile}")
    # FIXME fix uncertainty helpers
    #jpsi_crctn_unc_MC_helper = wremnants.make_jpsi_crctn_unc_helper(filepath=f"{common.data_dir}/calibration/calibrationJMC_smeared_v718_nominal.root")
    #jpsi_crctn_unc_data_helper = wremnants.make_jpsi_crctn_unc_helper(filepath=f"{common.data_dir}/calibration/calibrationJDATA_smeared_v718.root")

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = wremnants.make_muon_calibration_helpers()

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mumu", met=args.met)


def build_graph(df, dataset):
    logging.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")

    if dataset.is_data:
        cvh_helper = data_calibration_helper
        jpsi_helper = jpsi_crctn_data_helper if mass_fit else None
    else:
        cvh_helper = mc_calibration_helper
        jpsi_helper = jpsi_crctn_MC_helper if mass_fit else None

    df = muon_calibration.define_corrected_muons(df, cvh_helper, args.muonCorr, dataset)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, nMuons=2, use_trackerMuons=args.trackerMuons, use_isolation=True)

    df = muon_calibration.define_trigger_muons(df, jpsi_helper, args.muonCorr)

    df = df.Filter("Sum(trigMuons) == 1 && Sum(nonTrigMuons) == 1")
    df = df.Filter("nonTrigMuons_pt0 > 26.") # since this is our "neutrino", should there also be an upper cut (and eta cuts)?

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "trigMuons")
    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "nonTrigMuons")

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)
    df = muon_selections.select_trigger_muon(df, dataset, "trigMuons_eta0", "trigMuons_phi0")

    df = df.Define("TrigMuon_mom4", "ROOT::Math::PtEtaPhiMVector(trigMuons_pt0, trigMuons_eta0, trigMuons_phi0, wrem::muon_mass)")
    df = df.Define("NonTrigMuon_mom4", "ROOT::Math::PtEtaPhiMVector(nonTrigMuons_pt0, nonTrigMuons_eta0, nonTrigMuons_phi0, wrem::muon_mass)")
    df = df.Define("Z_mom4", "ROOT::Math::PxPyPzEVector(TrigMuon_mom4)+ROOT::Math::PxPyPzEVector(NonTrigMuon_mom4)")
    df = df.Define("ptZ", "Z_mom4.pt()")
    df = df.Define("massZ", "Z_mom4.mass()")
    df = df.Define("yZ", "Z_mom4.Rapidity()")
    df = df.Define("absYZ", "std::fabs(yZ)")
    df = df.Define("csSineCosThetaPhiZ", "trigMuons_charge0 == -1 ? wrem::csSineCosThetaPhi(TrigMuon_mom4, NonTrigMuon_mom4) : wrem::csSineCosThetaPhi(NonTrigMuon_mom4, TrigMuon_mom4)")

    df = df.Define("cosThetaStarZ", "csSineCosThetaPhiZ.costheta")
    df = df.Define("phiStarZ", "std::atan2(csSineCosThetaPhiZ.sinphi, csSineCosThetaPhiZ.cosphi)")

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

    if isW or isZ:
        df = syst_tools.define_mass_weights(df, isW)

    # Move the mass cut to apply to the dilepton plot
    df = df.Filter("massZ >= 60. && massZ < 120.")

    df_dilepton = df
    df_dilepton = df_dilepton.Filter("trigMuons_pt0 > 26.") # should there also be an upper cut (and eta cuts)?

    dilepton_axes = [axis_mll, axis_yll, axis_ptll, axis_costhetastarll, axis_phistarll]
    dilepton_cols = ["massZ", "yZ", "ptZ", "cosThetaStarZ", "phiStarZ"]
    if not args.csvars_hist:
        dilepton_axes = dilepton_axes[:-2]
        dilepton_cols = dilepton_cols[:-2]
    dilepton_cols.append("trigMuons_charge0")
    dilepton_axes.append(axis_charge)
    dilepton = df_dilepton.HistoBoost("dilepton", dilepton_axes, [*dilepton_cols, "nominal_weight"])
    results.append(dilepton)
    
    if not args.no_recoil:
        df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
        df = recoilHelper.setup_recoil_Z(df, results, dataset)
        df = recoilHelper.auxHists(df, results)
        df = recoilHelper.apply_recoil_Z(df, results, dataset, ["ZmumuPostVFP"])  # produces corrected MET as MET_corr_rec_pt/phi
        #if isZ: df = recoilHelper.recoil_Z_unc_lowPU(df, results, "", "", axis_mt, axis_mll)

        df = df.Define("transverseMass_uncorr", f"wrem::mt_wlike_nano(trigMuons_pt0, trigMuons_phi0, nonTrigMuons_pt0, nonTrigMuons_phi0, MET_corr_xy_pt, MET_corr_xy_phi)")
        results.append(df.HistoBoost("transverseMass_uncorr", [axis_mt], ["transverseMass_uncorr", "nominal_weight"]))
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    df = df.Define("transverseMass", f"wrem::mt_wlike_nano(trigMuons_pt0, trigMuons_phi0, nonTrigMuons_pt0, nonTrigMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
    results.append(df.HistoBoost("transverseMass", [axis_mt], ["transverseMass", "nominal_weight"]))
    
    df = df.Filter("transverseMass >= 45.") # 40 for Wmass, thus be 45 here (roughly half the boson mass)
    
    nominal_cols = ["trigMuons_eta0", "trigMuons_pt0", "trigMuons_charge0"]

    nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
    results.append(nominal)

    if not dataset.is_data and not args.onlyMainHistograms:
        if args.uncertainty_hist == "nominal":
            unc_df = df
            unc_cols = nominal_cols
            unc_axes = nominal_axes
        else:
            unc_df = df_dilepton
            unc_cols = dilepton_cols
            unc_axes = dilepton_axes

        unc_df = syst_tools.add_muon_efficiency_unc_hists(results, unc_df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, unc_axes, unc_cols, args.uncertainty_hist, is_w_like=True)
        unc_df = syst_tools.add_L1Prefire_unc_hists(results, unc_df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, unc_axes, unc_cols, args.uncertainty_hist)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isW or isZ:

            if args.theory_corr and dataset.name in corr_helpers:
                results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                    corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only,
                    with_uncertainties=args.uncertainty_hist == "nominal")
                )
                results.extend(theory_tools.make_theory_corr_hists(df_dilepton, "dilepton", dilepton_axes, dilepton_cols, 
                    corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only,
                    with_uncertainties=args.uncertainty_hist == "dilepton")
                )

            scale_axes = [*unc_axes, axis_ptVgen, axis_chargeVgen]
            scale_cols = [*unc_cols, "ptVgen", "chargeVgen"]
            syst_tools.add_scale_hist(results, unc_df, scale_axes, scale_cols, args.uncertainty_hist)
            syst_tools.add_pdf_hists(results, unc_df, dataset.name, unc_axes, unc_cols, args.pdfs, args.uncertainty_hist)

            if isZ:
                syst_tools.add_massweights_hist(results, unc_df, scale_axes, scale_cols, args.uncertainty_hist)
                # there is no W backgrounds for the Wlike, make QCD scale histograms only for Z
                # should probably remove the charge here, because the Z only has a single charge and the pt distribution does not depend on which charged lepton is selected
                if not args.skipHelicity:
                    # TODO: Should have consistent order here with the scetlib correction function
                    syst_tools.add_qcdScaleByHelicityUnc_hist(results, unc_df, qcdScaleByHelicity_helper, scale_axes, scale_cols, args.uncertainty_hist)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                syst_tools.add_muonscale_hist(results, unc_df, args.muonCorrEtaBins, args.muonCorrMag, isW, unc_axes, unc_cols, args.uncertainty_hist)


        df = unc_df

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, "mz_wlike_with_mu_eta_pt.pkl.lz4", args)
