import argparse
from utilities import output_tools, common, rdf_tools

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist
import lz4.frame
import logging
import math
import time

logging.basicConfig(level=logging.INFO)

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--noMuonCorr", action="store_true", help="Don't use corrected pt-eta-phi-charge")
parser.add_argument("--noScaleFactors", action="store_true", help="Don't use scale factors for efficiency")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--lumiUncertainty", type=float, help="Uncertainty for luminosity in excess to 1 (e.g. 1.012 means 1.2\%)", default=1.012)
args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts]) 
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9")

if not args.no_recoil:
    logging.warning("Recoil correction for high PU is not yet supported! Setting false")
    args.no_recoil = True

era = args.era
noMuonCorr = args.noMuonCorr

muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)
#qcdScaleByHelicity_Zhelper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
qcdScaleByHelicity_Whelper = wremnants.makeQCDScaleByHelicityHelper()

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

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")


nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

axis_chargeVgen = qcdScaleByHelicity_Whelper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False
)

# axes for mT measurement
axis_mt = hist.axis.Variable([0] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)
axis_eta_mT = hist.axis.Variable([-2.4, 2.4], name = "eta")

muon_efficiency_helper, muon_efficiency_helper_stat, muon_efficiency_helper_stat_tracking, muon_efficiency_helper_stat_reco, muon_efficiency_helper_syst = wremnants.make_muon_efficiency_helpers(era = era, max_pt = axis_pt.edges[-1])

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

calibration_helper, calibration_uncertainty_helper = wremnants.make_muon_calibration_helpers()

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    import ROOT
    ROOT.gInterpreter.Declare('#include "lowpu_recoil.h"')
    ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mu", met="RawPFMET") # DeepMETReso RawPFMET

def build_graph(df, dataset):
    print("build graph", dataset.name)
    results = []

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")

    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"
    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers
    if noMuonCorr:
        df = df.Alias("Muon_correctedPt", "Muon_pt")
        df = df.Alias("Muon_correctedEta", "Muon_eta")
        df = df.Alias("Muon_correctedPhi", "Muon_phi")
        df = df.Alias("Muon_correctedCharge", "Muon_charge")
    else:
        if dataset.is_data:
            #TODO corrections not available for data yet
            df = df.Alias("Muon_correctedPt", "Muon_cvhbsPt")
            df = df.Alias("Muon_correctedEta", "Muon_cvhbsEta")
            df = df.Alias("Muon_correctedPhi", "Muon_cvhbsPhi")
            df = df.Alias("Muon_correctedCharge", "Muon_cvhbsCharge")
        elif isW or isZ:
            df = df.Define("Muon_cvhbsMomCov", "wrem::splitNestedRVec(Muon_cvhbsMomCov_Vals, Muon_cvhbsMomCov_Counts)")
            df = wremnants.define_corrected_muons(df, calibration_helper)
        else:
            # no track refit available for background monte carlo samples and this is "good enough"
            df = df.Alias("Muon_correctedPt", "Muon_pt")
            df = df.Alias("Muon_correctedEta", "Muon_eta")
            df = df.Alias("Muon_correctedPhi", "Muon_phi")
            df = df.Alias("Muon_correctedCharge", "Muon_charge")
        
    #standalone quantities, currently only in data and W/Z samples
            
    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")
    df = df.Filter("Sum(vetoMuons) == 1")
    df = df.Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_isGlobal")
    df = df.Filter("Sum(goodMuons) == 1")

    # the corrected RECO muon kinematics, which is intended to be used as the nominal
    df = wremnants.define_corrected_reco_muon_kinematics(df)
    df = df.Define("goodMuons_abseta0", "abs(goodMuons_eta0)")
    df = df.Define("goodMuons_phi0", "Muon_correctedPhi[goodMuons][0]")
    df = df.Define("goodMuons_charge0", "Muon_correctedCharge[goodMuons][0]")

    if dataset.group in ["Top", "Diboson"]:
        df = df.Alias("goodMuons_SApt0",  "goodMuons_pt0")
        df = df.Alias("goodMuons_SAeta0", "goodMuons_eta0")
        df = df.Alias("goodMuons_SAphi0", "goodMuons_phi0")
    else:
        df = df.Define("goodMuons_SApt0",  "Muon_standalonePt[goodMuons][0]")
        df = df.Define("goodMuons_SAeta0", "Muon_standaloneEta[goodMuons][0]")
        df = df.Define("goodMuons_SAphi0", "Muon_standalonePhi[goodMuons][0]")
    df = df.Filter("goodMuons_SApt0 > 15.0 && wrem::deltaR2(goodMuons_SAeta0, goodMuons_SAphi0, goodMuons_eta0, goodMuons_phi0) < 0.09")
    
    # the cvhbs RECO muon kinematics, on top of which the corrections are applied to get nominal
    df = df.Define("goodMuons_pt0_cvhbs", "Muon_cvhbsPt[goodMuons][0]")
    df = df.Define("goodMuons_eta0_cvhbs", "Muon_cvhbsEta[goodMuons][0]")
    df = df.Define("goodMuons_phi0_cvhbs", "Muon_cvhbsPhi[goodMuons][0]")
    df = df.Define("goodMuons_charge0_cvhbs", "Muon_cvhbsCharge[goodMuons][0]")

    # the uncorrected RECO muon kinematics
    df = df.Define("goodMuons_pt0_uncrct", "Muon_pt[goodMuons][0]")
    df = df.Define("goodMuons_eta0_uncrct", "Muon_eta[goodMuons][0]")
    df = df.Define("goodMuons_phi0_uncrct", "Muon_phi[goodMuons][0]")
    df = df.Define("goodMuons_charge0_uncrct", "Muon_charge[goodMuons][0]")

    # define GEN muon kinematics
    #if dataset.is_data:
    if True:
        df = wremnants.define_gen_muons_from_reco_match(df, reco_subset = "goodMuons")
    
    if isW or isZ:
        # the cvhbs RECO muon kinematics, on top of which the corrections are applied to get nominal
        df = wremnants.define_cvhbs_reco_muon_kinematics(df)
        # the uncorrected RECO muon kinematics
        df = wremnants.define_uncrct_reco_muon_kinematics(df)

        # define GEN muon kinematics
        df = wremnants.get_good_gen_muons_idx_in_GenPart(df, reco_subset = "goodMuons")
        df = wremnants.define_good_gen_muons_kinematics(df)
        df = wremnants.define_good_gen_muon_kinematics(df)
        # calculate GEN quantities that are not in the nano
        df = wremnants.calculate_good_gen_muon_kinematics(df)

        df = df.Define("covMat_goodGenMuons0",
            ("wrem::getCovMatForGoodMuons0("
            "    Muon_cvhbsMomCov_Vals, Muon_cvhbsMomCov_Counts," 
            "    goodMuons, goodMuonsByGenTruth"
            ")")
        )
        df = df.Define("goodMuons_pt0_gen_smeared", 
            (
            "wrem::smearGenPt(covMat_goodGenMuons0, goodMuons_charge0_gen, "
            "goodMuons_pt0_gen, goodMuons_theta0_gen)"
            )
        )
        df = df.Define("goodMuons_qop0_gen_smeared", 
            "wrem::smearGenQop(covMat_goodGenMuons0, goodMuons_qop0_gen)"
        )
        df = df.Define("goodMuons_pt0_gen_smeared_a_la_qop", 
            "goodMuons_charge0_gen * std::sin(goodMuons_theta0_gen) / goodMuons_qop0_gen_smeared"
        )
        df = df.Define("goodMuons_qop0_gen_smeared_a_la_pt", 
            "goodMuons_charge0_gen * std::sin(goodMuons_theta0_gen) / goodMuons_pt0_gen_smeared"
        )
        df = df.Filter("covMat_goodGenMuons0[0] > 0 && covMat_goodGenMuons0[0] < 1")
        df = df.Define("goodMuons_eta0_gen_smeared", "goodMuons_eta0_gen")
        df = df.Define("goodMuons_phi0_gen_smeared", "goodMuons_phi0_gen")
        df = df.Define("goodMuons_charge0_gen_smeared", "goodMuons_charge0_gen")

        if args.validationHists:
            # corrected RECO / GEN
            df = df.Define("goodMuons_pt0_crctd_over_gen", "goodMuons_pt0/goodMuons_pt0_gen")
            df = df.Define("goodMuons_eta0_crctd_over_gen", "goodMuons_eta0/goodMuons_eta0_gen")
            df = df.Define("goodMuons_phi0_crctd_over_gen", "goodMuons_phi0/goodMuons_phi0_gen")
            # cvhbs RECO / GEN
            df = df.Define("goodMuons_pt0_cvhbs_over_gen", "goodMuons_pt0_cvhbs/goodMuons_pt0_gen")
            df = df.Define("goodMuons_eta0_cvhbs_over_gen", "goodMuons_eta0_cvhbs/goodMuons_eta0_gen")
            df = df.Define("goodMuons_phi0_cvhbs_over_gen", "goodMuons_phi0_cvhbs/goodMuons_phi0_gen")
            # uncorrected RECO / GEN
            df = df.Define("goodMuons_pt0_uncrct_over_gen", "goodMuons_pt0_uncrct/goodMuons_pt0_gen")
            df = df.Define("goodMuons_eta0_uncrct_over_gen", "goodMuons_eta0_uncrct/goodMuons_eta0_gen")
            df = df.Define("goodMuons_phi0_uncrct_over_gen", "goodMuons_phi0_uncrct/goodMuons_phi0_gen")
            # smeared GEN / GEN
            df = df.Define("goodMuons_pt0_gen_smeared_over_gen", "goodMuons_pt0_gen_smeared/goodMuons_pt0_gen")
            df = df.Define("goodMuons_eta0_gen_smeared_over_gen", "goodMuons_eta0_gen_smeared/goodMuons_eta0_gen")
            df = df.Define("goodMuons_phi0_gen_smeared_over_gen", "goodMuons_phi0_gen_smeared/goodMuons_phi0_gen")
            df = df.Define("goodMuons_qop0_gen_smeared_over_gen", "goodMuons_qop0_gen_smeared/goodMuons_qop0_gen")

    df = df.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")
    #TODO improve this to include muon mass?
    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")

    df = df.Filter("Sum(vetoElectrons) == 0")

    df = df.Define("goodCleanJets", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && Jet_pt > 30 && abs(Jet_eta) < 2.4 && wrem::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_correctedEta[vetoMuons],Muon_correctedPhi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")
    df = df.Define("goodCleanJetsPt45", "goodCleanJets && Jet_pt > 45")

    df = df.Define("passMT", "transverseMass >= 40.0")
    df = df.Filter("passMT || Sum(goodCleanJets)>=1")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    if dataset.group in ["Top", "Diboson"]:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    else:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Filter("wrem::hasTriggerMatch(goodMuons_eta0,goodMuons_phi0,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")
    
    # gen match to bare muons to select only prompt muons from top processes
    if isTop:
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 13")
        df = df.Filter("wrem::hasMatchDR2(goodMuons_eta0,goodMuons_phi0,GenPart_eta[postFSRmuons],GenPart_phi[postFSRmuons],0.09)")

    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]
    if dataset.is_data:
        df = df.Define("nominal_weight", "1.0")
        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)
        
        if not args.no_recoil:
            df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
            df = df.Define("mT_corr_rec", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_xy_pt, MET_corr_xy_phi)")
            results.append(df.HistoBoost("mT_corr_rec", [axis_mt, axis_eta_mT, axis_charge, axis_passIso], ["mT_corr_rec", "goodMuons_abseta0", "goodMuons_charge0", "passIso"]))

        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, nominal_cols)
        results.append(qcdJetPt45)
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId"])

        weight_expr = "weight*weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        if not args.noScaleFactors:
            weight_expr += "*weight_fullMuonSF_withTrackingReco"
        if args.vertex_weight:
            weight_expr += "*weight_vtx"
            
        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
            )
    if isW or isZ:
        nominal_cols_cvhbs = [
            "goodMuons_eta0_cvhbs", "goodMuons_pt0_cvhbs", "goodMuons_charge0_cvhbs",
            "passIso", "passMT"
        ]
        nominal_cols_uncrct = [
            "goodMuons_eta0_uncrct", "goodMuons_pt0_uncrct", "goodMuons_charge0_uncrct", 
            "passIso", "passMT"
        ]
        nominal_cols_gen = [
            "goodMuons_eta0_gen", "goodMuons_pt0_gen", "goodMuons_charge0_gen", 
            "passIso", "passMT"
        ]
        nominal_cols_gen_smeared = [
            "goodMuons_eta0_gen_smeared", "goodMuons_pt0_gen_smeared_a_la_qop", "goodMuons_charge0_gen_smeared",
            "passIso", "passMT"
        ]

        nominal_cvhbs =       df.HistoBoost("nominal_cvhbs", nominal_axes, [*nominal_cols_cvhbs, "nominal_weight"])
        nominal_uncrct =      df.HistoBoost("nominal_uncrct", nominal_axes, [*nominal_cols_uncrct, "nominal_weight"])
        nominal_gen =         df.HistoBoost("nominal_gen", nominal_axes, [*nominal_cols_gen, "nominal_weight"])
        nominal_gen_smeared = df.HistoBoost("nominal_gen_smeared", nominal_axes, [*nominal_cols_gen_smeared, "nominal_weight"])

        results.append(nominal_cvhbs)
        results.append(nominal_uncrct)
        results.append(nominal_gen)
        results.append(nominal_gen_smeared)

        if args.validationHists:
            nominal_cols_crctd_over_gen = [
                "goodMuons_pt0_crctd_over_gen"
            ]
            nominal_cols_cvhbs_over_gen = [
                "goodMuons_pt0_cvhbs_over_gen"
            ]
            nominal_cols_uncrct_over_gen = [
                "goodMuons_pt0_uncrct_over_gen"
            ]
            nominal_cols_gen_smeared_over_gen = [
                "goodMuons_pt0_gen_smeared_over_gen",
                "goodMuons_qop0_gen_smeared_over_gen"
            ]
            axis_pt_reco_over_gen = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "reco_pt_over_gen")
            axis_qop_reco_over_gen = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "reco_qop_over_gen")
            crctd_over_gen =  df.HistoBoost("crctd_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_crctd_over_gen, "nominal_weight"])
            cvhbs_over_gen =  df.HistoBoost("cvhbs_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_cvhbs_over_gen, "nominal_weight"])
            uncrct_over_gen = df.HistoBoost("uncrct_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_uncrct_over_gen, "nominal_weight"])
            gen_smeared_over_gen = df.HistoBoost("gen_smeared_over_gen", [axis_pt_reco_over_gen, axis_qop_reco_over_gen], [*nominal_cols_gen_smeared_over_gen, "nominal_weight"])
    
            results.append(crctd_over_gen)
            results.append(cvhbs_over_gen)
            results.append(uncrct_over_gen)
            results.append(gen_smeared_over_gen)

        if not args.no_recoil:
            df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
            df = recoilHelper.setup_gen(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"])
            df = recoilHelper.apply_W(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"]) # produces corrected MET as MET_corr_rec_pt/phi

            #df = df.Define("mT_corr_rec", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
            #results.append(df.HistoBoost("mT_corr_rec", [axis_mt, axis_eta_mT, axis_charge, axis_passIso], ["mT_corr_rec", "goodMuons_abseta0", "goodMuons_charge0", "passIso", "nominal_weight"]))
            #if dataset.name in ["WplusmunuPostVFP", "WminusmunuPostVFP"]: df = recoilHelper.recoil_W_unc_lowPU(df, results, axis_charge, axis_mt, axis_eta, axis_passIso)

        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(qcdJetPt45)

        
        df = df.Define("effStatTnP_tensor", muon_efficiency_helper_stat, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0", "passIso", "nominal_weight"])
        effStatTnP = df.HistoBoost("effStatTnP", nominal_axes, [*nominal_cols, "effStatTnP_tensor"], tensor_axes = muon_efficiency_helper_stat.tensor_axes)
        results.append(effStatTnP)

        # temporary solution, so ugly and awkward, but does what I want ...
        df = df.DefinePerSample("zero", "0")
        df = df.DefinePerSample("unity", "1")
        df = df.Define("effStatTnP_tracking_tensor", muon_efficiency_helper_stat_tracking, ["goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "unity", "nominal_weight"])
        effStatTnP_tracking = df.HistoBoost("effStatTnP_tracking", nominal_axes, [*nominal_cols, "effStatTnP_tracking_tensor"], tensor_axes = muon_efficiency_helper_stat_tracking.tensor_axes)
        results.append(effStatTnP_tracking)
        df = df.Define("effStatTnP_reco_tensor", muon_efficiency_helper_stat_reco, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0", "zero", "nominal_weight"])
        effStatTnP_reco = df.HistoBoost("effStatTnP_reco", nominal_axes, [*nominal_cols, "effStatTnP_reco_tensor"], tensor_axes = muon_efficiency_helper_stat_reco.tensor_axes)
        results.append(effStatTnP_reco)
        
        df = df.Define("effSystTnP_weight", muon_efficiency_helper_syst, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "nominal_weight"])
        effSystTnP = df.HistoBoost("effSystTnP", nominal_axes, [*nominal_cols, "effSystTnP_weight"], tensor_axes = muon_efficiency_helper_syst.tensor_axes)
        results.append(effSystTnP)

        df = df.Define("muonL1PrefireStat_tensor", muon_prefiring_helper_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId", "nominal_weight"])
        muonL1PrefireStat = df.HistoBoost("muonL1PrefireStat", nominal_axes, [*nominal_cols, "muonL1PrefireStat_tensor"], tensor_axes = muon_prefiring_helper_stat.tensor_axes)
        results.append(muonL1PrefireStat)

        df = df.Define("muonL1PrefireSyst_tensor", muon_prefiring_helper_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId", "nominal_weight"])
        muonL1PrefireSyst = df.HistoBoost("muonL1PrefireSyst", nominal_axes, [*nominal_cols, "muonL1PrefireSyst_tensor"], tensor_axes = [common.down_up_axis])
        results.append(muonL1PrefireSyst)

        df = df.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")
        ecalL1Prefire = df.HistoBoost("ecalL1Prefire", nominal_axes, [*nominal_cols, "ecalL1Prefire_tensor"], tensor_axes = [common.down_up_axis])
        results.append(ecalL1Prefire)
        
        # luminosity, done here as shape variation despite being a flat scaling so to facilitate propagating to fakes afterwards
        df = df.Define("luminosityScaling", f"wrem::constantScaling(nominal_weight, {args.lumiUncertainty})")
        luminosity = df.HistoBoost("luminosity", nominal_axes, [*nominal_cols, "luminosityScaling"], tensor_axes = [common.down_up_axis])
        results.append(luminosity)
                
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        
        if isW or isZ:

            df = theory_tools.define_scale_tensor(df)
            results.append(theory_tools.make_scale_hist(df, [*nominal_axes, axis_ptVgen, axis_chargeVgen], [*nominal_cols, "ptVgen", "chargeVgen"]))
            if isW and not args.skipHelicity:
                helicity_helper = qcdScaleByHelicity_Whelper
                # TODO: Should have consistent order here with the scetlib correction function
                df = df.Define("helicityWeight_tensor", helicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
                qcdScaleByHelicityUnc = df.HistoBoost("qcdScaleByHelicity", [*nominal_axes, axis_ptVgen, axis_chargeVgen], [*nominal_cols, "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=helicity_helper.tensor_axes)
                results.append(qcdScaleByHelicityUnc)

            df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)
            results.extend(theory_tools.make_pdf_hists(df, dataset.name, nominal_axes, nominal_cols, args.pdfs))

            masswargs = (nominal_axes, nominal_cols) if isW else (None, None)
            df, masswhist = syst_tools.define_mass_weights(df, isW, *masswargs)
            if masswhist:
                results.append(masswhist)
        
            df = df.Define("unity", "1.0")
            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                # TODO: Move to syst_tools
                netabins = args.muonCorrEtaBins
                nweights = 21
                mag = args.muonCorrMag
                df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, goodMuons_eta0, {mag}, {str(isW).lower()})")
                scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
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
                results.append(dummyMuonScaleSyst_gen_smear)
                results.append(dummyMuonScaleSyst)
            
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

            # TODO: put this type of histograms into a separate histmaker for just validation
            df = df.Define("muonScaleSyst_smearingWeightsPerSe_tensor", calibration_uncertainty_helper,
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
                "unity"
                ]
            )

            df = df.Define("smearing_weights_down", "muonScaleSyst_smearingWeightsPerSe_tensor(0,0)")
            df = df.Define("smearing_weights_up", "muonScaleSyst_smearingWeightsPerSe_tensor(0,1)")
    
            axis_smearing_weight = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "smearing_weight")
            smearing_weights_down = df.HistoBoost("smearing_weights_down", [*nominal_axes, axis_smearing_weight], [*nominal_cols, "smearing_weights_down"])
            smearing_weights_up = df.HistoBoost("smearing_weights_up", [*nominal_axes, axis_smearing_weight], [*nominal_cols, "smearing_weights_up"])
            results.append(smearing_weights_down)
            results.append(smearing_weights_up)
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

            df = df.Define("nominal_weights_perse", "nominal_weight")
            nominal_weights_hist = df.HistoBoost(
                "nominal_weights_hist",
                [*nominal_axes, axis_smearing_weight],
                [*nominal_cols, "nominal_weights_perse"]
            )
            results.append(nominal_weights_hist)

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
wremnants.transport_smearing_weights_to_reco(resultdict)

output_tools.write_analysis_output(resultdict, "mw_with_mu_eta_pt.pkl.lz4", args)
