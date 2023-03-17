import argparse
from utilities import output_tools, common, rdf_tools, logging

parser,initargs = common.common_parser(True)

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_calibration, muon_selections, muon_validation
import hist
import lz4.frame
import math
import time
from utilities import boostHistHelpers as hh
import pathlib
import os

data_dir = f"{pathlib.Path(__file__).parent}/../../wremnants/data/"
parser.add_argument("--noScaleFactors", action="store_true", help="Don't use scale factors for efficiency (legacy option for tests)")
parser.add_argument("--lumiUncertainty", type=float, help="Uncertainty for luminosity in excess to 1 (e.g. 1.012 means 1.2\%)", default=1.012)
parser.add_argument("--vqtTest", action="store_true", help="Test of isolation SFs dependence on V q_T projection (at the moment just for the W)")
parser.add_argument("--sfFileVqtTest", type=str, help="File with muon scale factors as a function of V q_T projection", default=f"{data_dir}/testMuonSF/fits_2.root")
parser.add_argument("--vqtTestIntegrated", action="store_true", help="Test of isolation SFs dependence on V q_T projection, integrated (would be the same as default SF, but pt-eta binning is different)")
parser.add_argument("--vqtTestReal", action="store_true", help="Test of isolation SFs dependence on V q_T projection, using 3D SFs directly (instead of the Vqt fits)")
parser.add_argument("--vqtTestIncludeTrigger", action="store_true", help="Test of isolation SFs dependence on V q_T projection. Including trigger")
args = parser.parse_args()

if args.vqtTestIntegrated:
    sfFileVqtTest = f"{data_dir}/testMuonSF/IsolationEfficienciesCoarseBinning.root"
else:
    sfFileVqtTest = args.sfFileVqtTest

logger = logging.setup_logger(__file__, args.verbose, args.color_logger)

datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              nanoVersion="v8" if args.v8 else "v9", base_path=args.data_path)

era = args.era

# custom template binning
template_neta = int(args.eta[0])
template_mineta = args.eta[1]
template_maxeta = args.eta[2]
print(f"Eta binning: {template_neta} bins from {template_mineta} to {template_maxeta}")
template_npt = int(args.pt[0])
template_minpt = args.pt[1]
template_maxpt = args.pt[2]
print(f"Pt binning: {template_npt} bins from {template_minpt} to {template_maxpt}")

# standard regular axes
axis_eta = hist.axis.Regular(template_neta, template_mineta, template_maxeta, name = "eta")
axis_pt = hist.axis.Regular(template_npt, template_minpt, template_maxpt, name = "pt")

axis_charge = common.axis_charge
axis_passIso = common.axis_passIso
axis_passMT = common.axis_passMT

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

# axes for study of fakes
axis_mt_fakes = hist.axis.Regular(120, 0., 120., name = "mt", underflow=False, overflow=True)
axis_iso_fakes = hist.axis.Regular(60, 0., 0.6, name = "PFrelIso04", underflow=False, overflow=True)
axis_hasjet_fakes = hist.axis.Boolean(name = "hasJets") # only need case with 0 jets or > 0 for now
mTStudyForFakes_axes = [axis_eta, axis_pt, axis_charge, axis_mt_fakes, axis_passIso, axis_hasjet_fakes]

axis_met = hist.axis.Regular(200, 0., 200., name = "met", underflow=False, overflow=True)

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper()
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False
)

if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile, era = era, max_pt = axis_pt.edges[-1], usePseudoSmoothing=True) 

    if args.vqtTest:
        if args.vqtTestReal:
            includeTrigger = False
            if args.vqtTestIncludeTrigger:
                includeTrigger = True
            muon_efficiency_helper_vqt, dummy_helper1, dummy_helper2 = wremnants.make_muon_efficiency_helpers_binned_vqt_real(filename = args.sfFile,
                                                                                                                              era = era,
                                                                                                                              max_pt = axis_pt.edges[-1],
                                                                                                                              includeTrigger = includeTrigger)
        else:
            if not args.vqtTestIntegrated:
                muon_efficiency_helper_vqt, dummy_helper1, dummy_helper2 = wremnants.make_muon_efficiency_helpers_binned_vqt(filename = args.sfFile, filenamevqt = sfFileVqtTest,
                                                                                                                             era = era,
                                                                                                                             max_pt = axis_pt.edges[-1]) 
            else:
                includeTrigger = False
                if args.vqtTestIncludeTrigger:
                    includeTrigger = True
                muon_efficiency_helper_vqt, muon_efficiency_helper_vqt_syst, dummy_helper2 = wremnants.make_muon_efficiency_helpers_binned_vqt_integrated(filename = args.sfFile, filenamevqt = sfFileVqtTest,
                                                                                                                                                          era = era,
                                                                                                                                                          max_pt = axis_pt.edges[-1],
                                                                                                                                                          includeTrigger = includeTrigger) 
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile, era = era, max_pt = axis_pt.edges[-1])
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

mc_jpsi_crctn_helper, data_jpsi_crctn_helper, jpsi_crctn_MC_unc_helper, jpsi_crctn_data_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, make_uncertainty_helper=True)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)

smearing_helper = muon_calibration.make_muon_smearing_helpers() if args.smearing else None

bias_helper = muon_calibration.make_muon_bias_helpers(args) if args.bias_calibration else None

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", args, flavor="mu")

def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"
    require_prompt = "tau" not in dataset.name # for muon GEN-matching

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")
    #df = df.Filter("event % 2 == 1") # test with odd/even events

    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=1)
    df = muon_selections.select_good_muons(df, nMuons=1, use_trackerMuons=args.trackerMuons, use_isolation=False)

    # the corrected RECO muon kinematics, which is intended to be used as the nominal
    df = muon_calibration.define_corrected_reco_muon_kinematics(df)

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "goodMuons")
 
    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)
    df = muon_selections.apply_triggermatching_muon(df, dataset, "goodMuons_eta0", "goodMuons_phi0")

    # gen match to bare muons to select only prompt muons from top processes
    if isTop:
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 13")
        df = df.Filter("wrem::hasMatchDR2(goodMuons_eta0,goodMuons_phi0,GenPart_eta[postFSRmuons],GenPart_phi[postFSRmuons],0.09)")

    if isW or isZ:
        df = muon_calibration.define_cvh_reco_muon_kinematics(df)
        #FIXME: make the smearing weights work without filtering on taus
        if args.smearingWeights:
            reco_sel = "vetoMuonsPre"
            df = muon_calibration.define_genFiltered_recoMuonSel(df, reco_sel, require_prompt)
            reco_sel_GF = muon_calibration.getColName_genFiltered_recoMuonSel(reco_sel, require_prompt)
            df = muon_calibration.define_covMatFiltered_recoMuonSel(df, reco_sel_GF)
            df = muon_calibration.define_matched_gen_muons_covMat(df, reco_sel_GF)
            df = muon_calibration.define_matched_gen_muons_kinematics(df, reco_sel_GF)
            df = muon_calibration.calculate_matched_gen_muon_kinematics(df, reco_sel_GF)
            df = muon_calibration.define_matched_genSmeared_muon_kinematics(df, reco_sel_GF)

            reco_sel = "goodMuons"
            df = muon_calibration.define_matched_gen_muons_kinematics(df, reco_sel)
            df = muon_calibration.calculate_matched_gen_muon_kinematics(df, reco_sel)
            df = muon_calibration.define_matched_gen_muons_covMat(df, reco_sel)
            df = muon_calibration.define_matched_genSmeared_muon_kinematics(df, reco_sel)

            for var in ['Pt', 'Eta', 'Charge']:
                df = df.Define(f"{reco_sel}_{var.lower()}0_gen", f"{reco_sel}_gen{var.capitalize()}[0]")
                df = df.Define(f"{reco_sel_GF}_{var.lower()}0_gen", f"{reco_sel_GF}_gen{var.capitalize()}[0]")

                df = df.Define(f"{reco_sel}_{var.lower()}0_gen_smeared", f"{reco_sel}_genSmeared{var.capitalize()}[0]")
                df = df.Define(f"{reco_sel_GF}_{var.lower()}0_gen_smeared", f"{reco_sel_GF}_genSmeared{var.capitalize()}[0]")

#            reco_sel = "vetoMuonsPre"
#            df = muon_calibration.define_matched_gen_muons_kinematics(df, reco_sel)
#            df = muon_calibration.calculate_matched_gen_muon_kinematics(df, reco_sel)
#            df = muon_calibration.define_matched_gen_muons_covMat(df, reco_sel)
#            df = muon_calibration.define_matched_genSmeared_muon_kinematics(df, reco_sel)

        if args.validationHists:
            for reco_type in ['crctd', 'cvh', 'uncrct', 'gen_smeared']:
                df = muon_calibration.define_reco_over_gen_cols(df, reco_type)

    df = df.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")

    # Jet collection actually has a pt threshold of 15 GeV in MiniAOD 
    df = df.Define("goodCleanJetsNoPt", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && abs(Jet_eta) < 2.4 && wrem::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_correctedEta[vetoMuons],Muon_correctedPhi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    ########################################################################
    # define event weights here since they are needed below for some helpers
    if dataset.is_data:
        df = df.Define("nominal_weight", "1.0")            
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        weight_expr = "weight*weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        if not args.noScaleFactors:
            if not (args.vqtTest and isW):
                df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
            else:
                df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 13")
                df = df.Define("postFSRnus", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 14")
                df = df.Define("postFSRnusIdx", "wrem::postFSRLeptonsIdx(postFSRnus)")
                df = df.Define("goodMuons_zqtproj0","wrem::zqtproj0(goodMuons_pt0, goodMuons_eta0, goodMuons_phi0, GenPart_pt, GenPart_eta, GenPart_phi, postFSRnusIdx)")
                if args.vqtTestIntegrated:
                    df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper_vqt, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
                else:
                    df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper_vqt, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "goodMuons_zqtproj0"])
            weight_expr += "*weight_fullMuonSF_withTrackingReco"
        if args.vertex_weight:
            weight_expr += "*weight_vtx"
        
        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)
        df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)
    ########################################################################
    
    if not args.no_recoil:
        lep_cols = ["goodMuons_pt0", "goodMuons_phi0", "goodMuons_charge0", "Muon_pt[goodMuons][0]"]
        df = recoilHelper.recoil_W(df, results, dataset, common.vprocs, lep_cols) # produces corrected MET as MET_corr_rec_pt/phi  vprocs_lowpu wprocs_recoil_lowpu
        df = recoilHelper.recoil_W_unc(df, results, dataset, common.vprocs)
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
    df = df.Define("hasCleanJet", "Sum(goodCleanJetsNoPt && Jet_pt > 30) >= 1")

    # couple of histograms specific for tests with fakes
    mTStudyForFakes = df.HistoBoost("mTStudyForFakes", mTStudyForFakes_axes, ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "transverseMass", "passIso", "hasCleanJet", "nominal_weight"])
    results.append(mTStudyForFakes)
    mtIsoJetCharge = df.HistoBoost("mtIsoJetCharge", [axis_mt_fakes, axis_iso_fakes, axis_hasjet_fakes, axis_charge], ["transverseMass", "goodMuons_pfRelIso04_all0", "hasCleanJet", "goodMuons_charge0", "nominal_weight"])
    results.append(mtIsoJetCharge)
    
    df = df.Define("passMT", "transverseMass >= 40.0")
    # no longer cut on jet at low mT, it biases the fakes estimate

    # utility plot, mt and met together in 2D, to plot them later
    mtAndMET = df.HistoBoost("mtAndMET", [axis_mt_fakes, axis_met, axis_charge, axis_passIso, axis_passMT], ["transverseMass", "MET_corr_rec_pt", "goodMuons_charge0", "passIso", "passMT", "nominal_weight"])
    results.append(mtAndMET)
    ##
    
    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if dataset.is_data:
        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)

        if not args.onlyMainHistograms:
            syst_tools.add_QCDbkg_jetPt_hist(results, df, nominal_axes, nominal_cols, jet_pt=30)

    else:  
        results.append(df.HistoBoost("nominal_weight", [hist.axis.Regular(200, -4, 4)], ["nominal_weight"]))

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
            )
        if args.smearingWeights and (isW or isZ): 
            nominal_cols_gen, nominal_cols_gen_smeared = muon_calibration.make_alt_reco_and_gen_hists(df, results, nominal_axes, reco_sel_GF)
            if args.validationHists: 
                wremnants.make_reco_over_gen_hists(df, results)

    if not dataset.is_data and not args.onlyMainHistograms:
        
        syst_tools.add_QCDbkg_jetPt_hist(results, df, nominal_axes, nominal_cols, jet_pt=30)

        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, nominal_axes, nominal_cols)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, nominal_axes, nominal_cols)

        if args.vqtTest:
            if args.vqtTestIntegrated:
                df = df.Define("effSystTnP_weight_vqt", muon_efficiency_helper_vqt_syst, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "nominal_weight"])
                effSystTnP_vqt = df.HistoBoost("effSystTnP_vqt", nominal_axes, [*nominal_cols, "effSystTnP_weight_vqt"], tensor_axes = muon_efficiency_helper_vqt_syst.tensor_axes)
                results.append(effSystTnP_vqt)
        
        # luminosity, done here as shape variation despite being a flat scaling so to facilitate propagating to fakes afterwards
        df = df.Define("luminosityScaling", f"wrem::constantScaling(nominal_weight, {args.lumiUncertainty})")
        luminosity = df.HistoBoost("nominal_luminosity", nominal_axes, [*nominal_cols, "luminosityScaling"], tensor_axes = [common.down_up_axis])
        results.append(luminosity)
                
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        
        if isW or isZ:

            df = theory_tools.define_scale_tensor(df)

            scale_axes = [*nominal_axes, axis_ptVgen, axis_chargeVgen]
            scale_cols = [*nominal_cols, "ptVgen", "chargeVgen"]
            syst_tools.add_qcdScale_hist(results, df, scale_axes, scale_cols)

            if isW and not args.skipHelicity:
                # TODO: Should have consistent order here with the scetlib correction function                    
                syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, scale_axes, scale_cols)

            syst_tools.add_pdf_hists(results, df, dataset.name, nominal_axes, nominal_cols, args.pdfs)

            df = syst_tools.define_mass_weights(df, dataset.name)
            if isW:
                syst_tools.add_massweights_hist(results, df, nominal_axes, nominal_cols, proc=dataset.name)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                df = syst_tools.add_muonscale_hist(results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, nominal_axes, nominal_cols)

                if args.smearingWeights:
                    df = syst_tools.add_muonscale_smeared_hist(results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, nominal_axes, nominal_cols_gen_smeared)

                # TODO: Move to syst_tools
                netabins = args.muonCorrEtaBins
                nweights = 23 if isZ else 21
                mag = args.muonCorrMag

                df = df.Define("unity", "1.0")
                df = df.Define(
                    f"muonScaleDummy{netabins}Bins_PerSe",
                    f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(unity, massWeight_tensor, goodMuons_eta0, {mag}, {str(isW).lower()})"
                )
                df = df.Define("massweights_down", f"muonScaleDummy{netabins}Bins_PerSe(0,0)")
                df = df.Define("massweights_up", f"muonScaleDummy{netabins}Bins_PerSe(1,0)")

                axis_mass_weight = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "axis_mass_weight")
                dummyMuonScaleSystPerSeDown = df.HistoBoost(
                    "nominal_muonScaleSystPerSeDown",
                    [axis_mass_weight],
                    ["massweights_down"]
                )
                dummyMuonScaleSystPerSeUp = df.HistoBoost(
                    "nominal_muonScaleSystPerSeUp",
                    [axis_mass_weight],
                    ["massweights_up"]
                )
                results.append(dummyMuonScaleSystPerSeDown)
                results.append(dummyMuonScaleSystPerSeUp)
                if args.smearingWeights:
                    df = df.DefinePerSample("bool_true", "true")
                    df = df.DefinePerSample("bool_false", "false")
                    if args.muonCorrData == "massfit" or "massfit_lbl":
                        jpsi_unc_helper = jpsi_crctn_data_unc_helper
                        df = df.Define("muonScaleSyst_responseWeights_tensor_gensmear", jpsi_unc_helper,
                            [
                            f"{reco_sel_GF}_genQop",
                            f"{reco_sel_GF}_genPhi",
                            f"{reco_sel_GF}_genEta",
                            f"{reco_sel_GF}_genSmearedQop",
                            f"{reco_sel_GF}_genSmearedPhi",
                            f"{reco_sel_GF}_genSmearedEta",
                            f"{reco_sel_GF}_genSmearedCharge",
                            f"{reco_sel_GF}_genSmearedPt",
                            f"{reco_sel_GF}_covMat",
                            "nominal_weight",
                            "bool_false"
                            ]
                        )
                        dummyMuonScaleSyst_responseWeights = df.HistoBoost(
                            "muonScaleSyst_responseWeights_gensmear", nominal_axes,
                            [*nominal_cols_gen_smeared, "muonScaleSyst_responseWeights_tensor_gensmear"],
                            tensor_axes = jpsi_unc_helper.tensor_axes
                        )
                        results.append(dummyMuonScaleSyst_responseWeights)
                    else:
                        df = df.Define("muonScaleSyst_responseWeights_tensor_gensmear", calibration_uncertainty_helper,
                            [
                            "goodMuons_qop0_gen_smeared",
                            "goodMuons_pt0_gen_smeared_a_la_qop",
                            "goodMuons_eta0_gen_smeared",
                            "goodMuons_phi0_gen_smeared",
                            "goodMuons_charge0_gen_smeared",
                            "covMat_goodGenMuons0",
                            "goodMuons_qop0_gen",
                            "goodMuons_pt0_gen",
                            "goodMuons_eta0_gen",
                            "goodMuons_phi0_gen",
                            "goodMuons_charge0_gen",
                            "nominal_weight"
                            ]
                        )
                        dummyMuonScaleSyst_responseWeights = df.HistoBoost("muonScaleSyst_responseWeights_gensmear", nominal_axes, [*nominal_cols_gen_smeared, "muonScaleSyst_responseWeights_tensor_gensmear"], tensor_axes = calibration_uncertainty_helper.tensor_axes)
                        results.append(dummyMuonScaleSyst_responseWeights)
            if args.validationHists and args.smearingWeights:
                df = wremnants.define_cols_for_smearing_weights(df, calibration_uncertainty_helper)
                wremnants.make_hists_for_smearing_weights(df, nominal_axes, nominal_cols, results)

            if args.smearingWeights:
                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_mil", "goodMuons_pt0_gen_smeared * 1.001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_mil", "goodMuons_pt0_gen_smeared / 1.001")
                df = df.Define("goodMuons_pt0_scaleUp_tenthmil", "goodMuons_pt0 * 1.0001")
                df = df.Define("goodMuons_pt0_scaleDn_tenthmil", "goodMuons_pt0 / 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_tenthmil", "goodMuons_pt0_gen_smeared * 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_tenthmil", "goodMuons_pt0_gen_smeared / 1.0001")
                muonScaleVariationUpMil = df.HistoBoost(
                    "nominal_muonScaleVariationUpMil", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationDnMil = df.HistoBoost(
                    "nominal_muonScaleVariationDnMil", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleDn_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationUpTenthmil = df.HistoBoost(
                    "nominal_muonScaleVariationUpTenthmil", 
                    nominal_axes,
                    [nominal_cols[0], "goodMuons_pt0_scaleUp_tenthmil", *nominal_cols[2:], "nominal_weight"]
                )
                muonScaleVariationDnTenthmil = df.HistoBoost(
                    "nominal_muonScaleVariationDnTenthmil", 
                    nominal_axes,
                    [nominal_cols[0], "goodMuons_pt0_scaleDn_tenthmil", *nominal_cols[2:], "nominal_weight"]
                )
                muonScaleVariationUpTenthmil_gen_smear = df.HistoBoost(
                    "nominal_muonScaleVariationUpTenthmil_gen_smear", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_tenthmil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationDnTenthmil_gen_smear = df.HistoBoost(
                    "nominal_muonScaleVariationDnTenthmil_gen_smear", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleDn_tenthmil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                results.append(muonScaleVariationUpMil)
                results.append(muonScaleVariationDnMil)
                results.append(muonScaleVariationUpTenthmil)
                results.append(muonScaleVariationDnTenthmil)
                results.append(muonScaleVariationUpTenthmil_gen_smear)
                results.append(muonScaleVariationDnTenthmil_gen_smear)

            df = df.Define("Muon_cvhMomCov", "wrem::splitNestedRVec(Muon_cvhMomCov_Vals, Muon_cvhMomCov_Counts)")

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
if args.smearingWeights:
    muon_calibration.transport_smearing_weights_to_reco(resultdict)
    muon_calibration.muon_scale_variation_from_manual_shift(resultdict)

output_tools.write_analysis_output(resultdict, "mw_with_mu_eta_pt.hdf5", args)
