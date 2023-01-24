import argparse
from utilities import output_tools, common, rdf_tools

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_calibration, muon_selections
import hist
import lz4.frame
import logging
import math
import time
from utilities import boostHistHelpers as hh
import pathlib

data_dir = f"{pathlib.Path(__file__).parent}/../../wremnants/data/"

logging.basicConfig(level=logging.INFO)

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--muonCorr", type=str, default="lbl", choices=["lbl", "none", "trackfit_only"], help="Type of correction to apply to the muons")
parser.add_argument("--noScaleFactors", action="store_true", help="Don't use scale factors for efficiency")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--lumiUncertainty", type=float, help="Uncertainty for luminosity in excess to 1 (e.g. 1.012 means 1.2\%)", default=1.012)
parser.add_argument(
    "--dataCrctn", type = str, choices = [None, 'cvh', 'jpsi_crctd'],
    default=None, help = "alternative correction to be applied to data"
)
parser.add_argument(
    "--MCCrctn", type = str, choices = [None, 'cvh', 'cvhbs', 'truth', 'jpsi_crctd'], 
    default = None, help = "alternative correction to be applied to MC"
)
parser.add_argument("--jpsiCrctnDataInput", type = str, default = None, help = "path to the root file for jpsi corrections on the data")
parser.add_argument("--jpsiCrctnMCInput", type = str, default = None, help = "path to the root file for jpsi corrections on the MC")
args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts]) 
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9", base_path=args.data_path)
    
era = args.era

muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper()
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False
)

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

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

# TODO: get from common
axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")
nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

# axes for recoil/MET
axis_mt = hist.axis.Variable([0] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)
axis_MET_pt = hist.axis.Regular(300, 0, 300, name = "MET_pt", underflow=False)
axis_recoil_magn = hist.axis.Regular(300, 0, 300, name = "recoil_magn", underflow=False)


# axes for study of fakes
axis_mt_fakes = hist.axis.Regular(60, 0., 120., name = "mt", underflow=False, overflow=True)
axis_njet_fakes = hist.axis.Regular(2, -0.5, 1.5, name = "Numbr of jets", underflow=False, overflow=False) # only need case with 0 jets or > 0

if args.binnedScaleFactors:
    logging.info("Using binned scale factors and uncertainties")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = axis_pt.edges[-1], usePseudoSmoothing=True) 
else:
    logging.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era, max_pt =
                                                                                                                                     axis_pt.edges[-1])
logging.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = wremnants.make_muon_calibration_helpers()

if args.dataCrctn == 'jpsi_crctd':
    jpsi_crctn_data_helper = muon_validation.make_jpsi_crctn_helper(args.jpsiCrctnDataInput)
    jpsi_crctn_unc_data_helper = muon_validation.make_jpsi_crctn_unc_helper(args.jpsiCrctnDataInput)
if args.MCCrctn == 'jpsi_crctd':
    jpsi_crctn_MC_helper = muon_validation.make_jpsi_crctn_helper(args.jpsiCrctnMCInput)
    jpsi_crctn_unc_MC_helper = muon_validation.make_jpsi_crctn_unc_helper(args.jpsiCrctnMCInput)

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mu", met=args.met)

# FIXME: Currently breaks the taus
smearing_weights = False

def build_graph(df, dataset):
    logging.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")
    #df = df.Filter("event % 2 == 1") # test with odd/even events

    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers

    calibration_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    df = muon_calibration.define_corrected_muons(df, calibration_helper, args.muonCorr, dataset)

    df = muon_selections.select_veto_muons(df, nMuons=1)
    df = muon_selections.select_good_muons(df, nMuons=1, use_trackerMuons=args.trackerMuons, use_isolation=False)

    # the corrected RECO muon kinematics, which is intended to be used as the nominal
    df = muon_calibration.define_corrected_reco_muon_kinematics(df)

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "goodMuons")
 
    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)
    df = muon_selections.select_trigger_muon(df, dataset, "goodMuons_eta0", "goodMuons_phi0")

    # gen match to bare muons to select only prompt muons from top processes
    if isTop:
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 13")
        df = df.Filter("wrem::hasMatchDR2(goodMuons_eta0,goodMuons_phi0,GenPart_eta[postFSRmuons],GenPart_phi[postFSRmuons],0.09)")

    if isW or isZ:
        df = muon_calibration.define_cvh_reco_muon_kinematics(df)
        #FIXME: make the smearing weights work without filtering on taus
        if smearing_weights:
            df = muon_calibration.define_uncrct_reco_muon_kinematics(df)
            df = muon_calibration.get_good_gen_muons_idx_in_GenPart(df, reco_subset = "goodMuons")
            df = muon_calibration.define_good_gen_muon_kinematics(df)
            df = muon_calibration.calculate_good_gen_muon_kinematics(df)
            df = muon_calibration.define_gen_smeared_muon_kinematics(df)

        if args.validationHists:
            for reco_type in ['crctd', 'cvh', 'uncrct', 'gen_smeared']:
                df = muon_calibration.define_reco_over_gen_cols(df, reco_type)

    df = df.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")

    df = df.Define("goodCleanJets", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && Jet_pt > 30 && abs(Jet_eta) < 2.4 && wrem::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_correctedEta[vetoMuons],Muon_correctedPhi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")
    df = df.Define("goodCleanJetsPt45", "goodCleanJets && Jet_pt > 45")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    if not args.no_recoil:
        if dataset.is_data:
            df = recoilHelper.setup_MET(df, results, dataset, "goodMuons_pt0", "goodMuons_phi0", "Muon_pt[goodMuons][0]")
            df = df.Alias("MET_corr_rec_pt", "MET_corr_xy_pt")
            df = df.Alias("MET_corr_rec_phi", "MET_corr_xy_phi")
        else:
            df = recoilHelper.setup_MET(df, results, dataset, "goodMuons_pt0", "goodMuons_phi0", "Muon_pt[goodMuons][0]")
            df = recoilHelper.setup_recoil_gen(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"])
            df = recoilHelper.apply_recoil_W(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"]) # produces corrected MET as MET_corr_rec_pt/phi
    else:
        df = df.Alias("MET_corr_rec_pt", "MET_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_phi")

    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
    df = df.Define("passMT", "transverseMass >= 40.0")
    df = df.Filter("passMT || Sum(goodCleanJets)>=1")


    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if dataset.is_data:
        df = df.Define("nominal_weight", "1.0")

        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)

        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, nominal_cols)
        results.append(qcdJetPt45)
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        weight_expr = "weight*weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        if not args.noScaleFactors:
            df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
            weight_expr += "*weight_fullMuonSF_withTrackingReco"
        if args.vertex_weight:
            weight_expr += "*weight_vtx"
        
        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)
        
        results.append(df.HistoBoost("nominal_weight", [hist.axis.Regular(200, -4, 4)], ["nominal_weight"]))

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
            )
        if smearing_weights and (isW or isZ): 
            nominal_cols_cvh, nominal_cols_uncrct, nominal_cols_gen, nominal_cols_gen_smeared = muon_calibration.make_alt_reco_and_gen_hists(df, results, nominal_axes)
            if args.validationHists: 
                wremnants.make_reco_over_gen_hists(df, results)

        
    results.append(df.HistoBoost("MET_pt", [axis_MET_pt, axis_charge, axis_passMT, axis_passIso], ["MET_corr_rec_pt", "goodMuons_charge0", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("transverseMass", [axis_mt, axis_charge, axis_passMT, axis_passIso], ["transverseMass", "goodMuons_charge0", "passMT", "passIso", "nominal_weight"]))   
    if isW or isZ:
        results.extend(theory_tools.make_theory_corr_hists(df, "transverseMass", [axis_mt, axis_charge, axis_passMT, axis_passIso], ["transverseMass", "goodMuons_charge0", "passMT", "passIso"],
            corr_helpers[dataset.name], args.theory_corr, modify_central_weight=False)
        )

    if not dataset.is_data and not args.onlyMainHistograms:
        
        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(qcdJetPt45)

        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, nominal_axes, nominal_cols)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, nominal_axes, nominal_cols)
        
        # luminosity, done here as shape variation despite being a flat scaling so to facilitate propagating to fakes afterwards
        df = df.Define("luminosityScaling", f"wrem::constantScaling(nominal_weight, {args.lumiUncertainty})")
        luminosity = df.HistoBoost("luminosity", nominal_axes, [*nominal_cols, "luminosityScaling"], tensor_axes = [common.down_up_axis])
        results.append(luminosity)
                
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        
        if isW or isZ:

            df = theory_tools.define_scale_tensor(df)
            df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)

            scale_axes = [*nominal_axes, axis_ptVgen, axis_chargeVgen]
            scale_cols = [*nominal_cols, "ptVgen", "chargeVgen"]
            syst_tools.add_scale_hist(results, df, scale_axes, scale_cols)

            if isW and not args.skipHelicity:
                # TODO: Should have consistent order here with the scetlib correction function                    
                syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, scale_axes, scale_cols)

            syst_tools.add_pdf_hists(results, df, dataset.name, nominal_axes, nominal_cols, args.pdfs)

            df = syst_tools.define_mass_weights(df, isW)
            if isW:
                syst_tools.add_massweights_hist(results, df, nominal_axes, nominal_cols)


            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                # TODO: Move to syst_tools
                netabins = args.muonCorrEtaBins
                nweights = 21
                mag = args.muonCorrMag
                df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, goodMuons_eta0, {mag}, {str(isW).lower()})")
                scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
                if smearing_weights:
                    dummyMuonScaleSyst_gen_smear = df.HistoBoost(
                        "muonScaleSyst_gen_smear",
                        nominal_axes,
                        [*nominal_cols_gen_smeared, f"muonScaleDummy{netabins}Bins"],
                        tensor_axes=[common.down_up_axis, scale_etabins_axis]
                    )
                dummyMuonScaleSyst = df.HistoBoost(
                    "muonScaleSyst",
                    nominal_axes,
                    [*nominal_cols, f"muonScaleDummy{netabins}Bins"],
                    tensor_axes=[common.down_up_axis, scale_etabins_axis]
                )
                if smearing_weights:
                    results.append(dummyMuonScaleSyst_gen_smear)
                results.append(dummyMuonScaleSyst)
                df = df.Define("unity", "1.0")
                df = df.Define(
                    f"muonScaleDummy{netabins}Bins_PerSe",
                    f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(unity, massWeight_tensor, goodMuons_eta0, {mag}, {str(isW).lower()})"
                )
                df = df.Define("massweights_down", f"muonScaleDummy{netabins}Bins_PerSe(0,0)")
                df = df.Define("massweights_up", f"muonScaleDummy{netabins}Bins_PerSe(1,0)")

                axis_mass_weight = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "axis_mass_weight")
                dummyMuonScaleSystPerSeDown = df.HistoBoost(
                    "muonScaleSystPerSeDown",
                    [axis_mass_weight],
                    ["massweights_down"]
                )
                dummyMuonScaleSystPerSeUp = df.HistoBoost(
                    "muonScaleSystPerSeUp",
                    [axis_mass_weight],
                    ["massweights_up"]
                )
                results.append(dummyMuonScaleSystPerSeDown)
                results.append(dummyMuonScaleSystPerSeUp)
                if smearing_weights:
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
            if args.validationHists and smearing_weights:
                df = wremnants.define_cols_for_smearing_weights(df, calibration_uncertainty_helper)
                wremnants.make_hists_for_smearing_weights(df, nominal_axes, nominal_cols, results)

            if smearing_weights:
                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_mil", "goodMuons_pt0_gen_smeared_a_la_qop * 1.001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_mil", "goodMuons_pt0_gen_smeared_a_la_qop / 1.001")
                df = df.Define("goodMuons_pt0_scaleUp_tenthmil", "goodMuons_pt0 * 1.0001")
                df = df.Define("goodMuons_pt0_scaleDn_tenthmil", "goodMuons_pt0 / 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleUp_tenthmil", "goodMuons_pt0_gen_smeared_a_la_qop * 1.0001")
                df = df.Define("goodMuons_pt0_gen_smeared_scaleDn_tenthmil", "goodMuons_pt0_gen_smeared_a_la_qop / 1.0001")
                muonScaleVariationUpMil = df.HistoBoost(
                    "muonScaleVariationUpMil", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationDnMil = df.HistoBoost(
                    "muonScaleVariationDnMil", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleDn_mil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationUpTenthmil = df.HistoBoost(
                    "muonScaleVariationUpTenthmil", 
                    nominal_axes,
                    [nominal_cols[0], "goodMuons_pt0_scaleUp_tenthmil", *nominal_cols[2:], "nominal_weight"]
                )
                muonScaleVariationDnTenthmil = df.HistoBoost(
                    "muonScaleVariationDnTenthmil", 
                    nominal_axes,
                    [nominal_cols[0], "goodMuons_pt0_scaleDn_tenthmil", *nominal_cols[2:], "nominal_weight"]
                )
                muonScaleVariationUpTenthmil_gen_smear = df.HistoBoost(
                    "muonScaleVariationUpTenthmil_gen_smear", 
                    nominal_axes,
                    [nominal_cols_gen_smeared[0], "goodMuons_pt0_gen_smeared_scaleUp_tenthmil", *nominal_cols_gen_smeared[2:], "nominal_weight"]
                )
                muonScaleVariationDnTenthmil_gen_smear = df.HistoBoost(
                    "muonScaleVariationDnTenthmil_gen_smear", 
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
if smearing_weights:
    wremnants.transport_smearing_weights_to_reco(resultdict)
    wremnants.muon_scale_variation_from_manual_shift(resultdict)

output_tools.write_analysis_output(resultdict, "mw_with_mu_eta_pt.pkl.lz4", args)
