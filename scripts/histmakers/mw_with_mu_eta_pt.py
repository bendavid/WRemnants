import argparse
from utilities import output_tools, common, rdf_tools, logging, differential

parser,initargs = common.common_parser(True)

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_calibration, muon_selections, muon_validation, unfolding_tools
import hist
import lz4.frame
import math
import time
from utilities import boostHistHelpers as hh
import pathlib
import os
import numpy as np

data_dir = f"{pathlib.Path(__file__).parent}/../../wremnants/data/"
parser.add_argument("--noScaleFactors", action="store_true", help="Don't use scale factors for efficiency (legacy option for tests)")
parser.add_argument("--lumiUncertainty", type=float, help="Uncertainty for luminosity in excess to 1 (e.g. 1.012 means 1.2\%)", default=1.012)
parser.add_argument("--vqtTest", action="store_true", help="Test of isolation SFs dependence on V q_T projection (at the moment just for the W)")
parser.add_argument("--sfFileVqtTest", type=str, help="File with muon scale factors as a function of V q_T projection", default=f"{data_dir}/testMuonSF/fits_2.root")
parser.add_argument("--vqtTestIntegrated", action="store_true", help="Test of isolation SFs dependence on V q_T projection, integrated (would be the same as default SF, but pt-eta binning is different)")
parser.add_argument("--vqtTestReal", action="store_true", help="Test of isolation SFs dependence on V q_T projection, using 3D SFs directly (instead of the Vqt fits)")
parser.add_argument("--vqtTestIncludeTrigger", action="store_true", help="Test of isolation SFs dependence on V q_T projection. Including trigger")
parser.add_argument("--noGenMatchMC", action='store_true', help="Don't use gen match filter for prompt muons with MC samples (note: QCD MC never has it anyway)")
parser.add_argument("--dphiMuonMetCut", type=float, help="Threshold to cut |deltaPhi| > thr*np.pi between muon and met", default=0.25)
args = parser.parse_args()

if args.vqtTestIntegrated:
    sfFileVqtTest = f"{data_dir}/testMuonSF/IsolationEfficienciesCoarseBinning.root"
else:
    sfFileVqtTest = args.sfFileVqtTest

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              nanoVersion="v8" if args.v8 else "v9", base_path=args.dataPath)

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
axis_eta = hist.axis.Regular(template_neta, template_mineta, template_maxeta, name = "eta", overflow=not args.excludeFlow, underflow=not args.excludeFlow)
axis_pt = hist.axis.Regular(template_npt, template_minpt, template_maxpt, name = "pt", overflow=not args.excludeFlow, underflow=not args.excludeFlow)

axis_charge = common.axis_charge
axis_passIso = common.axis_passIso
axis_passMT = common.axis_passMT

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

unfolding_axes, unfolding_cols = differential.get_pt_eta_axes(args.genBins, template_minpt, template_maxpt, template_maxeta)

# axes for study of fakes
axis_mt_fakes = hist.axis.Regular(120, 0., 120., name = "mt", underflow=False, overflow=True)
axis_iso_fakes = hist.axis.Regular(60, 0., 0.6, name = "PFrelIso04", underflow=False, overflow=True)
axis_dphi_fakes = hist.axis.Regular(16, 0., np.pi, name = "DphiMuonMet", underflow=False, overflow=False)
axis_hasjet_fakes = hist.axis.Boolean(name = "hasJets") # only need case with 0 jets or > 0 for now
mTStudyForFakes_axes = [axis_eta, axis_pt, axis_charge, axis_mt_fakes, axis_passIso, axis_hasjet_fakes, axis_dphi_fakes]

axis_met = hist.axis.Regular(200, 0., 200., name = "met", underflow=False, overflow=True)

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper()

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
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile, era = era, max_pt = axis_pt.edges[-1], directIsoSFsmoothing=args.directIsoSFsmoothing)
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

mc_jpsi_crctn_helper, data_jpsi_crctn_helper, jpsi_crctn_MC_unc_helper, jpsi_crctn_data_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, make_uncertainty_helper=True)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)

smearing_helper = muon_calibration.make_muon_smearing_helpers() if args.smearing else None

bias_helper = muon_calibration.make_muon_bias_helpers(args) if args.biasCalibration else None

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theoryCorr)

if args.nonClosureScheme == "binned":
    z_non_closure_binned_helper = muon_calibration.make_Z_non_closure_binned_helper(
        correlate = args.correlatedNonClosureNP
    )
else:
    z_non_closure_parametrized_helper = muon_calibration.make_Z_non_closure_parametrized_helper(
        correlate = args.correlatedNonClosureNP
    )


# recoil initialization
if not args.noRecoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", args, flavor="mu")

def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"
    isQCDMC = dataset.group == "QCD"
    require_prompt = "tau" not in dataset.name # for muon GEN-matching

    unfold = args.unfolding and dataset.name in ["WplusmunuPostVFP", "WminusmunuPostVFP"]

    # disable auxiliary histograms when unfolding to reduce memory consumptions
    auxiliary_histograms = not args.unfolding

    apply_theory_corr = args.theoryCorr and dataset.name in corr_helpers

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    if unfold:
        df = unfolding_tools.define_gen_level(df, args.genLevel, dataset.name, mode="wmass")
        unfolding_tools.add_xnorm_histograms(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols)

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")
    #df = df.Filter("event % 2 == 1") # test with odd/even events

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=1)
    df = muon_selections.select_good_muons(df, nMuons=1, use_trackerMuons=args.trackerMuons, use_isolation=False)

    # the corrected RECO muon kinematics, which is intended to be used as the nominal
    df = muon_calibration.define_corrected_reco_muon_kinematics(df)

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "goodMuons")
 
    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)
    df = muon_selections.apply_triggermatching_muon(df, dataset, "goodMuons_eta0", "goodMuons_phi0")

    # gen match to bare muons to select only prompt muons from MC processes, but also including tau decays
    # status flags in NanoAOD: https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2016ULpostVFP/doc_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1.html
    if not dataset.is_data and not isQCDMC and not args.noGenMatchMC:
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (5<<1)) && abs(GenPart_pdgId) == 13")
        df = df.Filter("wrem::hasMatchDR2(goodMuons_eta0,goodMuons_phi0,GenPart_eta[postFSRmuons],GenPart_phi[postFSRmuons],0.09)")

    if isW or isZ:
        df = muon_calibration.define_cvh_reco_muon_kinematics(df)
        #FIXME: make the smearing weights work without filtering on taus
        if args.muonScaleVariation == 'smearingWeights':
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

            for var in ['Pt', 'Eta', 'Charge', 'Qop']:
                df = df.Define(f"{reco_sel}_{var.lower()}0_gen", f"{reco_sel}_gen{var.capitalize()}[0]")
                df = df.Define(f"{reco_sel_GF}_{var.lower()}0_gen", f"{reco_sel_GF}_gen{var.capitalize()}[0]")

                df = df.Define(f"{reco_sel}_{var.lower()}0_gen_smeared", f"{reco_sel}_genSmeared{var.capitalize()}[0]")
                df = df.Define(f"{reco_sel_GF}_{var.lower()}0_gen_smeared", f"{reco_sel_GF}_genSmeared{var.capitalize()}[0]")
            df = df.Define(f"{reco_sel_GF}_covMat0", f"{reco_sel_GF}_covMat[0]")

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
        df = df.DefinePerSample("nominal_weight", "1.0")            
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        weight_expr = "weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        if not args.noScaleFactors:
            if not (args.vqtTest and isW):
                df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
            else:
                df = df.Define("postFSRnus", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 14")
                df = df.Define("postFSRnusIdx", "wrem::postFSRLeptonsIdx(postFSRnus)")
                df = df.Define("goodMuons_zqtproj0","wrem::zqtproj0(goodMuons_pt0, goodMuons_eta0, goodMuons_phi0, GenPart_pt, GenPart_eta, GenPart_phi, postFSRnusIdx)")
                if args.vqtTestIntegrated:
                    df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper_vqt, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
                else:
                    df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper_vqt, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "goodMuons_zqtproj0"])
            weight_expr += "*weight_fullMuonSF_withTrackingReco"
        if not args.noVertexWeight:
            weight_expr += "*weight_vtx"
        
        df = df.Define("exp_weight", weight_expr)
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)
    ########################################################################
    
    if not args.noRecoil:
        lep_cols = ["goodMuons_pt0", "goodMuons_phi0", "goodMuons_charge0", "Muon_pt[goodMuons][0]"]
        df = recoilHelper.recoil_W(df, results, dataset, common.vprocs, lep_cols) # produces corrected MET as MET_corr_rec_pt/phi  vprocs_lowpu wprocs_recoil_lowpu
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
    df = df.Define("hasCleanJet", "Sum(goodCleanJetsNoPt && Jet_pt > 30) >= 1")

    df = df.Define("deltaPhiMuonMet", "std::abs(wrem::deltaPhi(goodMuons_phi0,MET_corr_rec_phi))")

    if auxiliary_histograms: 
        # couple of histograms specific for tests with fakes
        mTStudyForFakes = df.HistoBoost("mTStudyForFakes", mTStudyForFakes_axes, ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "transverseMass", "passIso", "hasCleanJet", "deltaPhiMuonMet", "nominal_weight"])
        results.append(mTStudyForFakes)

    dphiMuonMetCut = args.dphiMuonMetCut * np.pi
    # add filter of deltaPhi(muon,met) before other histograms (but before the previous histogram for test with fakes)
    df = df.Filter(f"deltaPhiMuonMet > {dphiMuonMetCut}") # pi/4 was found to be a good threshold for signal with mT > 40 GeV

    if auxiliary_histograms:
        mtIsoJetCharge = df.HistoBoost("mtIsoJetCharge", [axis_mt_fakes, axis_iso_fakes, axis_hasjet_fakes, axis_charge], ["transverseMass", "goodMuons_pfRelIso04_all0", "hasCleanJet", "goodMuons_charge0", "nominal_weight"])
        results.append(mtIsoJetCharge)
    
    df = df.Define("passMT", "transverseMass >= 40.0")

    if auxiliary_histograms:
        # utility plot, mt and met, to plot them later
        results.append(df.HistoBoost("MET", [axis_met, axis_charge, axis_passIso, axis_passMT], ["MET_corr_rec_pt", "goodMuons_charge0", "passIso", "passMT", "nominal_weight"]))
        results.append(df.HistoBoost("transverseMass", [axis_mt_fakes, axis_charge, axis_passIso, axis_passMT], ["transverseMass", "goodMuons_charge0", "passIso", "passMT", "nominal_weight"]))
    
    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if unfold:
        axes = [*nominal_axes, *unfolding_axes] 
        cols = [*nominal_cols, *unfolding_cols]
    else:
        axes = nominal_axes
        cols = nominal_cols

    if not args.onlyMainHistograms:
        syst_tools.add_QCDbkg_jetPt_hist(results, df, axes, cols, jet_pt=30)

    if dataset.is_data:
        nominal = df.HistoBoost("nominal", axes, cols)
        results.append(nominal)

    else:  
        nominal = df.HistoBoost("nominal", axes, [*cols, "nominal_weight"])
        results.append(nominal)

        results.append(df.HistoBoost("nominal_weight", [hist.axis.Regular(200, -4, 4)], ["nominal_weight"], storage=hist.storage.Double()))

        if not args.noRecoil:
            df = recoilHelper.add_recoil_unc_W(df, results, dataset, cols, axes, "nominal")

        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal", axes, cols, 
                corr_helpers[dataset.name], args.theoryCorr, modify_central_weight=not args.theoryCorrAltOnly, isW = isW)
            )
        if args.muonScaleVariation == 'smearingWeights' and (isW or isZ): 
            nominal_cols_gen, nominal_cols_gen_smeared = muon_calibration.make_alt_reco_and_gen_hists(df, results, axes, cols, reco_sel_GF)
            if args.validationHists: 
                wremnants.make_reco_over_gen_hists(df, results)

    if not dataset.is_data and not args.onlyMainHistograms:
        
        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, axes, cols)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, axes, cols)

        if args.vqtTest:
            if args.vqtTestIntegrated:
                df = df.Define("effSystTnP_weight_vqt", muon_efficiency_helper_vqt_syst, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "nominal_weight"])
                effSystTnP_vqt = df.HistoBoost("effSystTnP_vqt", axes, [*cols, "effSystTnP_weight_vqt"], tensor_axes = muon_efficiency_helper_vqt_syst.tensor_axes, storage=hist.storage.Double())
                results.append(effSystTnP_vqt)
        
        # luminosity, done here as shape variation despite being a flat scaling so to facilitate propagating to fakes afterwards
        df = df.Define("luminosityScaling", f"wrem::constantScaling(nominal_weight, {args.lumiUncertainty})")
        luminosity = df.HistoBoost("nominal_luminosity", axes, [*cols, "luminosityScaling"], tensor_axes = [common.down_up_axis], storage=hist.storage.Double())
        results.append(luminosity)
                
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        
        if isW or isZ:

            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, axes, cols, for_wmass=True)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                df = syst_tools.add_muonscale_hist(results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, axes, cols)

                if args.muonScaleVariation == 'smearingWeights':
                    df = syst_tools.add_muonscale_smeared_hist(results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, axes, nominal_cols_gen_smeared)

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
                    ["massweights_down"], 
                    storage=hist.storage.Double()
                )
                dummyMuonScaleSystPerSeUp = df.HistoBoost(
                    "nominal_muonScaleSystPerSeUp",
                    [axis_mass_weight],
                    ["massweights_up"], 
                    storage=hist.storage.Double()
                )
                results.append(dummyMuonScaleSystPerSeDown)
                results.append(dummyMuonScaleSystPerSeUp)

                if (
                    (args.muonCorrData == "massfit" or "massfit_lbl") and 
                    (args.muonScaleVariation == 'smearingWeights')
                ):
                    # muon scale variation from stats. uncertainty on the jpsi massfit
                    df = df.DefinePerSample("bool_true", "true")
                    df = df.DefinePerSample("bool_false", "false")
                    if args.validateByMassWeights:
                        jpsi_unc_helper = muon_validation.make_jpsi_crctn_unc_helper_massweights(
                            "wremnants/data/calibration/calibrationJDATA_rewtgr_3dmap_LBL.root",
                            nweights,
                            scale = 3.04
                        )
                        df = df.Define("muonScaleSyst_responseWeights_tensor_gensmear", jpsi_unc_helper,
                            [
                                f"{reco_sel_GF}_eta0_gen_smeared",
                                f"{reco_sel_GF}_charge0_gen_smeared",
                                f"{reco_sel_GF}_pt0_gen_smeared",
                                "massWeight_tensor",
                                "nominal_weight",
                                f"bool_{str(isW).lower()}"
                            ]
                        )
                    else:
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
                        "muonScaleSyst_responseWeights_gensmear", axes,
                        [*nominal_cols_gen_smeared, "muonScaleSyst_responseWeights_tensor_gensmear"],
                        tensor_axes = jpsi_unc_helper.tensor_axes, storage=hist.storage.Double()
                    )
                    results.append(dummyMuonScaleSyst_responseWeights)

                    # for the Z non-closure nuisances
                    if args.nonClosureScheme == "A-M-separated":
                        df = df.DefinePerSample("AFlag", "0x01")
                        df = df.DefinePerSample("MFlag", "0x04")

                        df = df.Define("Z_non_closure_parametrized_A", z_non_closure_parametrized_helper,
                            [
                                f"{reco_sel_GF}_qop0_gen",
                                f"{reco_sel_GF}_eta0_gen",
                                f"{reco_sel_GF}_qop0_gen_smeared",
                                f"{reco_sel_GF}_eta0_gen_smeared",
                                f"{reco_sel_GF}_charge0_gen_smeared",
                                f"{reco_sel_GF}_pt0_gen_smeared",
                                f"{reco_sel_GF}_covMat0",
                                "nominal_weight",
                                "AFlag"
                            ]
                        )
                        hist_Z_non_closure_parametrized_A = df.HistoBoost(
                            "Z_non_closure_parametrized_A_gensmear",
                            nominal_axes,
                            [*nominal_cols_gen_smeared, "Z_non_closure_parametrized_A"],
                            tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double()
                        )
                        results.append(hist_Z_non_closure_parametrized_A)

                        df = df.Define("Z_non_closure_parametrized_M", z_non_closure_parametrized_helper,
                            [
                                f"{reco_sel_GF}_qop0_gen",
                                f"{reco_sel_GF}_eta0_gen",
                                f"{reco_sel_GF}_qop0_gen_smeared",
                                f"{reco_sel_GF}_eta0_gen_smeared",
                                f"{reco_sel_GF}_charge0_gen_smeared",
                                f"{reco_sel_GF}_pt0_gen_smeared",
                                f"{reco_sel_GF}_covMat0",
                                "nominal_weight",
                                "MFlag"
                            ]
                        )
                        hist_Z_non_closure_parametrized_M = df.HistoBoost(
                            "Z_non_closure_parametrized_M_gensmear",
                            nominal_axes,
                            [*nominal_cols_gen_smeared, "Z_non_closure_parametrized_M"],
                            tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double()
                        )
                        results.append(hist_Z_non_closure_parametrized_M)
                    elif args.nonClosureScheme == "A-M-combined":
                        df = df.DefinePerSample("AMFlag", "0x01 | 0x04")
                        df = df.Define("Z_non_closure_parametrized", z_non_closure_parametrized_helper,
                            [
                                f"{reco_sel_GF}_qop0_gen",
                                f"{reco_sel_GF}_eta0_gen",
                                f"{reco_sel_GF}_qop0_gen_smeared",
                                f"{reco_sel_GF}_eta0_gen_smeared",
                                f"{reco_sel_GF}_charge0_gen_smeared",
                                f"{reco_sel_GF}_pt0_gen_smeared",
                                f"{reco_sel_GF}_covMat0",
                                "nominal_weight",
                                "AMFlag"
                            ])
                        hist_Z_non_closure_parametrized = df.HistoBoost(
                            "Z_non_closure_parametrized_gensmear",
                            nominal_axes,
                            [*nominal_cols_gen_smeared, "Z_non_closure_parametrized"],
                            tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double()
                        )
                        results.append(hist_Z_non_closure_parametrized)
                    elif args.nonClosureScheme == "binned":
                        df = df.Define("Z_non_closure_binned", z_non_closure_binned_helper,
                            [
                                f"{reco_sel_GF}_qop0_gen",
                                f"{reco_sel_GF}_pt0_gen",
                                f"{reco_sel_GF}_eta0_gen",
                                f"{reco_sel_GF}_charge0_gen",
                                f"{reco_sel_GF}_qop0_gen_smeared",
                                f"{reco_sel_GF}_pt0_gen_smeared",
                                f"{reco_sel_GF}_eta0_gen_smeared",
                                f"{reco_sel_GF}_charge0_gen_smeared",
                                f"{reco_sel_GF}_covMat0",
                                "nominal_weight"
                            ]
                        )
                        hist_Z_non_closure_binned = df.HistoBoost(
                            "Z_non_closure_binned_gensmear",
                            nominal_axes,
                            [*nominal_cols_gen_smeared, "Z_non_closure_binned"],
                            tensor_axes = z_non_closure_binned_helper.tensor_axes,
                            storage=hist.storage.Double()
                        )
                        results.append(hist_Z_non_closure_binned)

            if args.muonScaleVariation == 'smearingWeights':
                if args.validationHists:
                    df = wremnants.define_cols_for_smearing_weights(df, calibration_uncertainty_helper)
                    wremnants.make_hists_for_smearing_weights(df, axes, cols, results)

                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_mil", "goodMuons_pt0_gen_smeared * 1.001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_mil", "goodMuons_pt0_gen_smeared / 1.001")
                df = df.Define("goodMuons_pt0_scaleUp_tenthmil", "goodMuons_pt0 * 1.0001")
                df = df.Define("goodMuons_pt0_scaleDn_tenthmil", "goodMuons_pt0 / 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_tenthmil", "goodMuons_pt0_gen_smeared * 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_tenthmil", "goodMuons_pt0_gen_smeared / 1.0001")
                muonScaleVariationUpMil = df.HistoBoost(
                    "nominal_muonScaleVariationUpMil", 
                    axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
                )
                muonScaleVariationDnMil = df.HistoBoost(
                    "nominal_muonScaleVariationDnMil", 
                    axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleDn_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
                )
                muonScaleVariationUpTenthmil = df.HistoBoost(
                    "nominal_muonScaleVariationUpTenthmil", 
                    axes,
                    [cols[0], "goodMuons_pt0_scaleUp_tenthmil", *cols[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
                )
                muonScaleVariationDnTenthmil = df.HistoBoost(
                    "nominal_muonScaleVariationDnTenthmil", 
                    axes,
                    [cols[0], "goodMuons_pt0_scaleDn_tenthmil", *cols[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
                )
                muonScaleVariationUpTenthmil_gen_smear = df.HistoBoost(
                    "nominal_muonScaleVariationUpTenthmil_gen_smear", 
                    axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_tenthmil", *nominal_cols_gen_smeared[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
                )
                muonScaleVariationDnTenthmil_gen_smear = df.HistoBoost(
                    "nominal_muonScaleVariationDnTenthmil_gen_smear", 
                    axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleDn_tenthmil", *nominal_cols_gen_smeared[2:], "nominal_weight"], 
                    storage=hist.storage.Double()
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
if not args.onlyMainHistograms and args.muonScaleVariation == 'smearingWeights':
    muon_calibration.transport_smearing_weights_to_reco(resultdict, nonClosureScheme = args.nonClosureScheme)
    muon_calibration.muon_scale_variation_from_manual_shift(resultdict)

output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args, update_name=not args.forceDefaultName)
