from utilities import boostHistHelpers as hh,common,output_tools

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation
import hist
import lz4.frame
import logging
import math
import time

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--muonCorr", type=str, default="trackfit_only", choices=["lbl", "none", "trackfit_only"], help="Type of correction to apply to the muons")
parser.add_argument("--muScaleMag", type=float, default=1e-4, help="Magnitude of dummy muon scale uncertainty")
parser.add_argument("--muScaleBins", type=int, default=1, help="Number of bins for muon scale uncertainty")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--csvars_hist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--uncertainty-hist", type=str, choices=["dilepton", "nominal"], default="nominal", help="Histogram to store uncertainties for")
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

axis_ptll = hist.axis.Variable(
    [0, 2, 3, 4, 4.75, 5.5, 6.5, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100, 150], name = "ptll"
)

axis_costhetastarll = hist.axis.Regular(20, -1., 1., name = "costhetastarll")
axis_phistarll = hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phistarll")

# extra axes which can be used to label tensor_axes

muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile, era = era, max_pt = axis_pt.edges[-1], is_w_like = True)
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

if muon_corr == "mass_fit":
    jpsi_crctn_data_helper = wremnants.make_jpsi_crctn_helper(filepath=)
    jpsi_crctn_unc_data_helper = wremnants.make_jpsi_crctn_unc_helper(filepath=)
    jpsi_crctn_MC_helper = wremnants.make_jpsi_crctn_helper(filepath = args.jpsiCrctnMCInput)
    jpsi_crctn_unc_MC_helper = wremnants.make_jpsi_crctn_unc_helper(filepath = args.jpsiCrctnMCInput)
mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = wremnants.make_muon_calibration_helpers()

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mumu", met=args.met)


#def add_plots_with_systematics

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

    calibration_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    df = wremnants.define_corrected_muons(df, calibration_helper, args.muonCorr, dataset)

    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")
    df = df.Filter("Sum(vetoMuons) == 2")
    
    if args.trackerMuons:
        if dataset.group in ["Top", "Diboson"]:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_highPurity")
        else:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_innerTrackOriginalAlgo != 13 && Muon_innerTrackOriginalAlgo != 14 && Muon_highPurity")
    else:
        df = df.Define("Muon_category", "Muon_isGlobal")

    df = df.Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_category && Muon_pfRelIso04_all < 0.15")
    df = df.Filter("Sum(goodMuons) == 2")

    # mu- for even event numbers, mu+ for odd event numbers
    df = df.Define("TrigMuon_charge", "event % 2 == 0 ? -1 : 1")
    df = df.Define("NonTrigMuon_charge", "-TrigMuon_charge")

    df = df.Define("trigMuons", "goodMuons && Muon_correctedCharge == TrigMuon_charge")
    df = df.Define("nonTrigMuons", "goodMuons && Muon_correctedCharge == NonTrigMuon_charge")

    df = df.Filter("Sum(trigMuons) == 1 && Sum(nonTrigMuons) == 1")

    df = df.Filter("NonTrigMuon_pt > 26.")

    df = df.Define("TrigMuon_isStandalone", "Muon_isStandalone[trigMuons][0]")
    df = df.Define("NonTrigMuon_isStandalone", "Muon_isStandalone[nonTrigMuons][0]")
    # central NanoAOD for backgrounds do not have standalone variables, therefore cannot use cut on them when running backgrounds (for now we neglect this bias)
    # when using tracker muons the presence of standalone muons is not guaranteed, so again cannot use or cut on those (also in data or Z, regardless whether branches exist)
    # it was verified that tracking scale factors can still be used as a function of inner track coordinates, despite having been measured versus standalone ones
    # (this is true also in the case of global muons, no significant difference is observed applying tracking SF using inner or outer tracks)
    if dataset.group in ["Top", "Diboson"]:
        df = df.Alias("TrigMuon_SApt",  "TrigMuon_pt")
        df = df.Alias("TrigMuon_SAeta", "TrigMuon_eta")
        df = df.Alias("TrigMuon_SAphi", "TrigMuon_phi")
        df = df.Alias("NonTrigMuon_SApt",  "NonTrigMuon_pt")
        df = df.Alias("NonTrigMuon_SAeta", "NonTrigMuon_eta")
        df = df.Alias("NonTrigMuon_SAphi", "NonTrigMuon_phi")
    elif args.trackerMuons:
        # try to use standalone variables when possible
        df = df.Define("TrigMuon_SApt",  "TrigMuon_isStandalone ? Muon_standalonePt[trigMuons][0] : TrigMuon_pt")
        df = df.Define("TrigMuon_SAeta", "TrigMuon_isStandalone ? Muon_standaloneEta[trigMuons][0] : TrigMuon_eta")
        df = df.Define("TrigMuon_SAphi", "TrigMuon_isStandalone ? Muon_standalonePhi[trigMuons][0] : TrigMuon_phi")
        df = df.Define("NonTrigMuon_SApt",  "NonTrigMuon_isStandalone ? Muon_standalonePt[nonTrigMuons][0] : NonTrigMuon_pt")
        df = df.Define("NonTrigMuon_SAeta", "NonTrigMuon_isStandalone ? Muon_standaloneEta[nonTrigMuons][0] : NonTrigMuon_eta")
        df = df.Define("NonTrigMuon_SAphi", "NonTrigMuon_isStandalone ? Muon_standalonePhi[nonTrigMuons][0] : NonTrigMuon_phi")
    else:
        df = df.Define("TrigMuon_SApt",  "Muon_standalonePt[trigMuons][0]")
        df = df.Define("TrigMuon_SAeta", "Muon_standaloneEta[trigMuons][0]")
        df = df.Define("TrigMuon_SAphi", "Muon_standalonePhi[trigMuons][0]")
        df = df.Define("NonTrigMuon_SApt",  "Muon_standalonePt[nonTrigMuons][0]")
        df = df.Define("NonTrigMuon_SAeta", "Muon_standaloneEta[nonTrigMuons][0]")
        df = df.Define("NonTrigMuon_SAphi", "Muon_standalonePhi[nonTrigMuons][0]")

    # these next cuts are mainly needed for consistency with the reco efficiency measurement for the case with global muons
    # note, when SA does not exist this cut is still fine because of how we define these variables
    df = df.Filter("TrigMuon_SApt > 15.0 && wrem::deltaR2(TrigMuon_SAeta, TrigMuon_SAphi, TrigMuon_eta, TrigMuon_phi) < 0.09")
    df = df.Filter("NonTrigMuon_SApt > 15.0 && wrem::deltaR2(NonTrigMuon_SAeta, NonTrigMuon_SAphi, NonTrigMuon_eta, NonTrigMuon_phi) < 0.09")
    if common.muonEfficiency_standaloneNumberOfValidHits > 0 and not args.trackerMuons and not dataset.group in ["Top", "Diboson"]:
        nHitsSA = common.muonEfficiency_standaloneNumberOfValidHits
        df = df.Filter(f"Muon_standaloneNumberOfValidHits[trigMuons][0] >= {nHitsSA} && Muon_standaloneNumberOfValidHits[nonTrigMuons][0] >= {nHitsSA}")
    
    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")

    df = df.Filter("Sum(vetoElectrons) == 0")

    if dataset.group in ["Top", "Diboson"]:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    else:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Filter("wrem::hasTriggerMatch(TrigMuon_eta,TrigMuon_phi,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")


    df = df.Define("TrigMuon_mom4", "ROOT::Math::PtEtaPhiMVector(TrigMuon_pt, TrigMuon_eta, TrigMuon_phi, wrem::muon_mass)")
    df = df.Define("NonTrigMuon_mom4", "ROOT::Math::PtEtaPhiMVector(NonTrigMuon_pt, NonTrigMuon_eta, NonTrigMuon_phi, wrem::muon_mass)")
    df = df.Define("Z_mom4", "ROOT::Math::PxPyPzEVector(TrigMuon_mom4)+ROOT::Math::PxPyPzEVector(NonTrigMuon_mom4)")
    df = df.Define("ptZ", "Z_mom4.pt()")
    df = df.Define("massZ", "Z_mom4.mass()")
    df = df.Define("yZ", "Z_mom4.Rapidity()")
    df = df.Define("absYZ", "std::fabs(yZ)")
    #if dataset.is_data or (isZ or isW):
    #    df = df.Define("TrigMuon_cvh_mom4",
    #        (
    #            "ROOT::Math::PtEtaPhiMVector("
    #            "TrigMuon_cvh_pt, TrigMuon_cvh_eta, TrigMuon_cvh_phi, wrem::muon_mass)"
    #        )
    #    )
    #    df = df.Define("NonTrigMuon_cvh_mom4",
    #        (
    #            "ROOT::Math::PtEtaPhiMVector("
    #            "NonTrigMuon_cvh_pt, NonTrigMuon_cvh_eta, NonTrigMuon_cvh_phi, wrem::muon_mass)"
    #        )
    #    )
    #    df = df.Define("Z_cvh_mom4", "ROOT::Math::PxPyPzEVector(TrigMuon_cvh_mom4)+ROOT::Math::PxPyPzEVector(NonTrigMuon_cvh_mom4)")
    #    df = df.Define("massZ_cvh", "Z_cvh_mom4.mass()")
    #    if (dataset.is_data and args.dataCrctn == 'jpsi_crctd') or ((isZ or isW) and args.MCCrctn == 'jpsi_crctd'):
    #        wremnants.define_jpsi_crctd_z_mass(df)
    #        wremnants.define_jpsi_crctd_unc_z_mass(df)
    df = df.Define("csSineCosThetaPhiZ", "TrigMuon_charge == -1 ? wrem::csSineCosThetaPhi(TrigMuon_mom4, NonTrigMuon_mom4) : wrem::csSineCosThetaPhi(NonTrigMuon_mom4, TrigMuon_mom4)")

    df = df.Define("cosThetaStarZ", "csSineCosThetaPhiZ.costheta")
    df = df.Define("phiStarZ", "std::atan2(csSineCosThetaPhiZ.sinphi, csSineCosThetaPhiZ.cosphi)")

    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers

    if not dataset.is_data:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["TrigMuon_pt", "TrigMuon_eta", "TrigMuon_SApt", "TrigMuon_SAeta", "TrigMuon_charge",
                                                                                      "NonTrigMuon_pt", "NonTrigMuon_eta", "NonTrigMuon_SApt", "NonTrigMuon_SAeta", "NonTrigMuon_charge"])
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

    # dilepton plots go here, before mass or transverse mass cuts
    df = df.Filter("massZ >= 60. && massZ < 120.")

    df_dilepton = df
    df_dilepton = df_dilepton.Filter("TrigMuon_pt > 26.")

    dilepton_axes = [axis_mll, axis_yll, axis_ptll, axis_costhetastarll, axis_phistarll]
    #if (dataset.is_data and args.dataCrctn):
    #    dilepton_cols = [f"massZ_{args.dataCrctn}", "yZ", "ptZ", "cosThetaStarZ", "phiStarZ"]
    #    if args.dataCrctn == 'jpsi_crctd':
    #        dilepton_unc_cols = [
    #            f"massZ_{args.dataCrctn}", "yZ", "ptZ", "cosThetaStarZ", "phiStarZ", "massZ_jpsi_crctd_unc"]
    #elif ((isZ or isW) and args.MCCrctn):
    #    dilepton_cols = [f"massZ_{args.MCCrctn}", "yZ", "ptZ", "cosThetaStarZ", "phiStarZ"]
    #    if args.MCCrctn == 'jpsi_crctd':
    #        dilepton_unc_cols = [f"massZ_jpsi_crctd_unc"]
    #else:
    dilepton_cols = ["massZ", "yZ", "ptZ", "cosThetaStarZ", "phiStarZ"]
    if not args.csvars_hist:
        dilepton_axes = dilepton_axes[:-2]
        dilepton_cols = dilepton_cols[:-2]
    dilepton_cols.append("TrigMuon_charge")
    dilepton_axes.append(axis_charge)
    
    dilepton = df_dilepton.HistoBoost("dilepton", dilepton_axes, [*dilepton_cols, "nominal_weight"])
    results.append(dilepton)

    if not args.no_recoil:
        df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
        df = recoilHelper.setup_recoil_Z(df, results, dataset)
        df = recoilHelper.auxHists(df, results)
        df = recoilHelper.apply_recoil_Z(df, results, dataset, ["ZmumuPostVFP"])  # produces corrected MET as MET_corr_rec_pt/phi
        #if isZ: df = recoilHelper.recoil_Z_unc_lowPU(df, results, "", "", axis_mt, axis_mll)
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    if apply_theory_corr:
        results.extend(theory_tools.make_theory_corr_hists(df_dilepton, "dilepton", dilepton_axes, dilepton_cols, 
            corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only,
            with_uncertainties=True)
        )

    dilepton = df_dilepton.HistoBoost("dilepton", dilepton_axes, [*dilepton_cols, "nominal_weight"])
    results.append(dilepton)

    #TODO improve this to include muon mass?
    met_vars = ("MET_pt", "MET_phi")
    df = df.Define("transverseMass_uncorr", f"wrem::mt_wlike_nano(TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi, {', '.join(met_vars)})")
    results.append(df.HistoBoost("transverseMass_uncorr", [axis_mt], ["transverseMass_uncorr", "nominal_weight"]))
    met_vars = (x.replace("xy", "rec") for x in met_vars)

    df = df.Define("transverseMass", f"wrem::mt_wlike_nano(TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi, {', '.join(met_vars)})")
    results.append(df.HistoBoost("transverseMass", [axis_mt], ["transverseMass", "nominal_weight"]))
    
    df = df.Filter("transverseMass >= 45.") # 40 for Wmass, thus be 45 here (roughly half the boson mass)
    
    nominal_cols = ["TrigMuon_eta", "TrigMuon_pt", "TrigMuon_charge"]

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

        unc_df = syst_tools.add_muon_efficiency_unc_hists(results, unc_df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, unc_axes, unc_cols)
        unc_df = syst_tools.add_muon_prefire_unc_hists(results, unc_df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, unc_axes, unc_cols)

        df = unc_df
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isW or isZ:

            if apply_theory_corr:
                results.extend(theory_tools.make_theory_corr_hists(df, "nominal", axes=nominal_axes, cols=nominal_cols, 
                    helpers=corr_helpers[dataset.name], generators=args.theory_corr, modify_central_weight=not args.theory_corr_alt_only,
                    with_uncertainties=True)
                )

            scale_axes = unc_axes if args.uncertainty_hist == "dilepton" else [*unc_axes, axis_ptVgen, axis_chargeVgen]
            scale_cols = unc_cols if args.uncertainty_hist == "dilepton" else [*nominal_cols, "ptVgen", "chargeVgen"]
            syst_tools.add_scale_hist(results, unc_df, scale_axes, scale_cols)
            syst_tools.add_pdf_hists(results, unc_df, dataset.name, nominal_axes, nominal_cols, args.pdfs)

            if isZ:
                # there is no W backgrounds for the Wlike, make QCD scale histograms only for Z
                # should probably remove the charge here, because the Z only has a single charge and the pt distribution does not depend on which charged lepton is selected
                if not args.skipHelicity:
                    # TODO: Should have consistent order here with the scetlib correction function
                    df = df.Define("helicityWeight_tensor", qcdScaleByHelicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
                    qcdScaleByHelicityUnc = df.HistoBoost("qcdScaleByHelicity", [*nominal_axes, axis_ptVgen, axis_chargeVgen], [*nominal_cols, "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=qcdScaleByHelicity_helper.tensor_axes)
                    results.append(qcdScaleByHelicityUnc)
                syst_tools.add_massweights_hist(results, df_dilepton, "dilepton", dilepton_axes, dilepton_cols)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                # TODO: Move to syst_tools
                netabins = args.muonCorrEtaBins
                nweights = 21
                mag = args.muonCorrMag
                df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, TrigMuon_eta, {mag}, {str(isW).lower()})")
                scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
                dummyMuonScaleSyst = df.HistoBoost("muonScaleSyst", nominal_axes, [*nominal_cols, f"muonScaleDummy{netabins}Bins"],
                    tensor_axes=[common.down_up_axis, scale_etabins_axis])

                results.append(dummyMuonScaleSyst)

            df = df.Define("Muon_cvhMomCov", "wrem::splitNestedRVec(Muon_cvhMomCov_Vals, Muon_cvhMomCov_Counts)")

            '''
            df = df.Define("muonScaleSyst_responseWeights_tensor", calibration_uncertainty_helper,
                           ["Muon_correctedPt",
                            "Muon_correctedEta",
                            "Muon_correctedPhi",
                            "Muon_correctedCharge",
                            "Muon_genPartIdx",
                            "Muon_cvhMomCov",
                            "vetoMuonsPre",
                            "GenPart_pt",
                            "GenPart_eta",
                            "GenPart_phi",
                            "GenPart_pdgId",
                            "GenPart_statusFlags",
                            "nominal_weight"])

            dummyMuonScaleSyst_responseWeights = df.HistoBoost("muonScaleSyst_responseWeights", nominal_axes, [*nominal_cols, "muonScaleSyst_responseWeights_tensor"], tensor_axes = calibration_uncertainty_helper.tensor_axes)
            results.append(dummyMuonScaleSyst_responseWeights)
            '''

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, "mz_wlike_with_mu_eta_pt.pkl.lz4", args)
