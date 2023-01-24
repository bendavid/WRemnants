import ROOT
import pathlib
import hist
import narf
from utilities import rdf_tools
from utilities import common
from utilities import boostHistHelpers as hh
from . import muon_validation

ROOT.gInterpreter.Declare('#include "muon_calibration.h"')
ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_calibration_helpers(mc_filename=data_dir+"/calibration/correctionResults_v718_idealgeom_gensim.root", 
        data_filename=data_dir+"/calibration/correctionResults_v718_recjpsidata.root", 
        era = None):

    mc_helper = ROOT.wrem.CVHCorrector(mc_filename)
    data_helper = ROOT.wrem.CVHCorrector(data_filename)

    uncertainty_hist = get_dummy_uncertainties()
    uncertainty_hist_cpp = narf.hist_to_pyroot_boost(uncertainty_hist, tensor_rank = 2)
    # min gen pt = 9 GeV to avoid threshold effects
    # max weight = 10 to protect against outliers
    uncertainty_helper = ROOT.wrem.calibration_uncertainty_helper[type(uncertainty_hist_cpp)](ROOT.std.move(uncertainty_hist_cpp), 9., 10.)

    down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
    uncertainty_helper.tensor_axes = (uncertainty_hist.axes["calvar"], down_up_axis)

    return mc_helper, data_helper, uncertainty_helper

def make_muon_bias_helpers(filename=data_dir+"/calibration/correctionResults_v718_idealgeom_gensim.root"):
    # this helper adds a correction for the bias from nonclosure to the muon pT

    helper = ROOT.wrem.BiasCorrector(filename)

    uncertainty_hist = get_dummy_uncertainties()
    uncertainty_hist_cpp = narf.hist_to_pyroot_boost(uncertainty_hist, tensor_rank = 2)
    # min gen pt = 9 GeV to avoid threshold effects
    # max weight = 10 to protect against outliers
    uncertainty_helper = ROOT.wrem.calibration_uncertainty_helper[type(uncertainty_hist_cpp)](ROOT.std.move(uncertainty_hist_cpp), 9., 10.)

    down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
    uncertainty_helper.tensor_axes = (uncertainty_hist.axes["calvar"], down_up_axis)

    return helper, uncertainty_helper

def get_dummy_uncertainties():
    axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
    axis_calvar = hist.axis.Integer(0, 1, underflow=False, overflow=False, name = "calvar")
    axis_calparm = hist.axis.Integer(0, 4, underflow=False, overflow=False, name = "calparm")

    h = hist.Hist(axis_eta, axis_calvar, axis_calparm)

    # forward b-field-like
    h.values()[..., 0, 0] = 1e-4

    # set underflow and overflow to match boundaries
    h.values(flow=True)[0, ...] = h.values(flow=True)[1, ...]
    h.values(flow=True)[-1, ...] = h.values(flow=True)[-2, ...]

    return h

def define_corrected_muons_wmass(df, helper, corr_type, dataset):
    if not (dataset.is_data or dataset.name in common.vprocs):
        corr_type = "none" 

    if corr_type == "none":
        df = df.Alias("Muon_correctedPt", "Muon_pt")
        df = df.Alias("Muon_correctedEta", "Muon_eta")
        df = df.Alias("Muon_correctedPhi", "Muon_phi")
        df = df.Alias("Muon_correctedCharge", "Muon_charge")
    elif "trackfit_only" in corr_type:
        # TODO: Is this a reasonable configuration?
        fit = "cvhideal" if corr_type == "trackfit_only_mctruth" and not dataset.is_data else "cvh"
        print("Trackfit", fit)
        df = df.Define("Muon_correctedPt", f"Muon_{fit}Pt")
        df = df.Define("Muon_correctedEta", f"Muon_{fit}Eta")
        df = df.Define("Muon_correctedPhi", f"Muon_{fit}Phi")
        df = df.Define("Muon_correctedCharge", f"Muon_{fit}Charge")
    elif corr_type == "lbl":
        corr_branch = "cvh" if dataset.is_data else "cvhideal"

        # split the nested vectors
        df = df.Define("Muon_cvhmergedGlobalIdxs", "wrem::splitNestedRVec(Muon_cvhmergedGlobalIdxs_Vals, Muon_cvhmergedGlobalIdxs_Counts)")
        df = df.Define(f"Muon_{corr_branch}JacRef", f"wrem::splitNestedRVec(Muon_{corr_branch}JacRef_Vals, Muon_{corr_branch}JacRef_Counts)")

        df = df.Define("Muon_correctedMom4Charge", helper, [f"Muon_{corr_branch}Pt", f"Muon_{corr_branch}Eta", f"Muon_{corr_branch}Phi", f"Muon_{corr_branch}Charge", "Muon_cvhmergedGlobalIdxs", f"Muon_{corr_branch}JacRef"])

        # split into individual vectors
        df = df.Define("Muon_correctedPt", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Pt(); } ); return res;")
        df = df.Define("Muon_correctedEta", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Eta(); } ); return res;")
        df = df.Define("Muon_correctedPhi", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Phi(); } ); return res;")
        df = df.Define("Muon_correctedCharge", "ROOT::VecOps::RVec<int> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.second; }); return res;")
    elif "massfit" in corr_type:
        raise ValueError(f"mass_fit not defined for W analysis")
    else:
        raise ValueError(f"Invalid correction type choice {corr_type}")
    return df

def define_corrected_muons_wlike(df, cvh_helper, jpsi_helper, corr_type, dataset, trackerMuons):
    if not (dataset.is_data or dataset.name in common.vprocs):
        corr_type = "none" 

    if corr_type == "none":
        df = df.Alias("Muon_correctedPt", "Muon_pt")
        df = df.Alias("Muon_correctedEta", "Muon_eta")
        df = df.Alias("Muon_correctedPhi", "Muon_phi")
        df = df.Alias("Muon_correctedCharge", "Muon_charge")
    elif "lbl" in corr_type:
        corr_branch = "cvh" if dataset.is_data else "cvhideal"

        # split the nested vectors
        df = df.Define("Muon_cvhmergedGlobalIdxs", "wrem::splitNestedRVec(Muon_cvhmergedGlobalIdxs_Vals, Muon_cvhmergedGlobalIdxs_Counts)")
        df = df.Define(f"Muon_{corr_branch}JacRef", f"wrem::splitNestedRVec(Muon_{corr_branch}JacRef_Vals, Muon_{corr_branch}JacRef_Counts)")

        df = df.Define("Muon_correctedMom4Charge", cvh_helper, [f"Muon_{corr_branch}Pt", f"Muon_{corr_branch}Eta", f"Muon_{corr_branch}Phi", f"Muon_{corr_branch}Charge", "Muon_cvhmergedGlobalIdxs", f"Muon_{corr_branch}JacRef"])

        # split into individual vectors
        df = df.Define("Muon_correctedPt", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Pt(); } ); return res;")
        df = df.Define("Muon_correctedEta", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Eta(); } ); return res;")
        df = df.Define("Muon_correctedPhi", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Phi(); } ); return res;")
        df = df.Define("Muon_correctedCharge", "ROOT::VecOps::RVec<int> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.second; }); return res;")
    else:
        fit = "cvhideal" if corr_type == "trackfit_only_mctruth" and not dataset.is_data else "cvh"
        df = df.Define("Muon_correctedPt", f"Muon_{fit}Pt")
        df = df.Define("Muon_correctedEta", f"Muon_{fit}Eta")
        df = df.Define("Muon_correctedPhi", f"Muon_{fit}Phi")
        df = df.Define("Muon_correctedCharge", f"Muon_{fit}Charge")
    
    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")

    if trackerMuons:
        if dataset.group in ["Top", "Diboson"]:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_highPurity")
        else:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_innerTrackOriginalAlgo != 13 && Muon_innerTrackOriginalAlgo != 14 && Muon_highPurity")
    else:
        df = df.Define("Muon_category", "Muon_isGlobal")

    df = df.Filter("Sum(vetoMuons) == 2")
    df = df.Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_category && Muon_pfRelIso04_all < 0.15")
    df = df.Filter("Sum(goodMuons) == 2")

    # mu- for even event numbers, mu+ for odd event numbers
    df = df.Define("TrigMuon_charge", "event % 2 == 0 ? -1 : 1")
    df = df.Define("NonTrigMuon_charge", "-TrigMuon_charge")

    df = df.Define("trigMuons", "goodMuons && Muon_correctedCharge == TrigMuon_charge")
    df = df.Define("nonTrigMuons", "goodMuons && Muon_correctedCharge == NonTrigMuon_charge")

    df = df.Define("TrigMuon_preCorr_pt", "Muon_correctedPt[trigMuons][0]")
    df = df.Define("TrigMuon_eta", "Muon_correctedEta[trigMuons][0]")
    df = df.Define("TrigMuon_phi", "Muon_correctedPhi[trigMuons][0]")

    df = df.Define("NonTrigMuon_preCorr_pt", "Muon_correctedPt[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_eta", "Muon_correctedEta[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_phi", "Muon_correctedPhi[nonTrigMuons][0]")

    # For corr_type == "massfit_lbl" use the truth-assisted lbl for MC
    if corr_type == "massfit":
        df = muon_validation.define_jpsi_crctd_muons_pt(df, jpsi_helper)
        df = df.Alias("TrigMuon_pt", "TrigMuon_jpsi_crctd_pt")
        df = df.Alias("NonTrigMuon_pt", "NonTrigMuon_jpsi_crctd_pt")
    else:
        df = df.Alias("TrigMuon_pt", "TrigMuon_preCorr_pt")
        df = df.Alias("NonTrigMuon_pt", "NonTrigMuon_preCorr_pt")

    return df

def get_good_gen_muons_idx_in_GenPart(df, reco_subset = "goodMuons"):
    df = df.Define("goodMuons_idx", f"Muon_genPartIdx[{reco_subset}]")
    df = df.Define("goodMuonsByGenTruth",
        (
        "ROOT::VecOps::RVec<bool> res(goodMuons_idx.size());"
        "for (int i = 0; i < goodMuons_idx.size(); i++) {"
        "    res[i] = ("
        "        (GenPart_statusFlags[goodMuons_idx[i]] & 0x01) &&"
        "        (abs(GenPart_pdgId[goodMuons_idx[i]]) == 13) &&"
        "        (GenPart_status[goodMuons_idx[i]] == 1)"
        "    );"
        "}"
        "return res;"
        )
    )
    df = df.Define("goodGenMuons_idx","goodMuons_idx[goodMuonsByGenTruth]")
    df = df.Filter("goodGenMuons_idx.size() > 0")
    return df

def define_good_gen_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        if var == "charge": var = "pdgId"
        col_name = f"goodMuons_{var.lower()}0_gen"
        selection = f"ROOT::VecOps::Take(GenPart_{var}, goodGenMuons_idx)[0]"
        df = df.Define(col_name, selection)
        if var == "pdgId":
            df = df.Define("goodMuons_charge0_gen", "goodMuons_pdgid0_gen > 0 ? -1 : 1")
    return df

def calculate_good_gen_muon_kinematics(df):
    df = df.Define("goodMuons_theta0_gen",
        "2. * std::atan(std::exp(-double(goodMuons_eta0_gen)))")
    df = df.Define("goodMuons_p0_gen", "double(goodMuons_pt0_gen) / std::sin(goodMuons_theta0_gen)")
    df = df.Define("goodMuons_qop0_gen", "double(goodMuons_charge0_gen) / goodMuons_p0_gen")
    return df

def define_gen_smeared_muon_kinematics(df):
    df = df.Define("covMat_goodGenMuons0",
        ("wrem::getCovMatForGoodMuons0("
        "    Muon_cvhMomCov_Vals, Muon_cvhMomCov_Counts," 
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
    return df

def define_corrected_reco_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        df = df.Define(
            f"goodMuons_{var.lower()}0",
            f"Muon_corrected{var.capitalize()}[goodMuons][0]"
        )
    return df

def define_cvh_reco_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        df = df.Define(
            f"goodMuons_{var.lower()}0_cvh",
            f"Muon_cvh{var.capitalize()}[goodMuons][0]"
        )
    return df

def define_uncrct_reco_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        df = df.Define(
            f"goodMuons_{var.lower()}0_uncrct",
            f"Muon_{var.lower()}[goodMuons][0]"
        )
    return df

def transport_smearing_weights_to_reco(
    resultdict, procs = ['WplusmunuPostVFP', 'WminusmunuPostVFP', 'ZmumuPostVFP'],
    proj_axes = ['eta', 'pt', 'charge', 'passIso', 'passMT']
):
    for proc in procs:
        proc_hists = resultdict[proc]['output']
        nominal_gen_smear = (
            proc_hists['nominal_gen_smeared'].project(*proj_axes))
        msv_sw_gen_smear = [
            (proc_hists['muonScaleSyst_responseWeights_gensmear'][...,0,0].project(*proj_axes)),
            (proc_hists['muonScaleSyst_responseWeights_gensmear'][...,0,1].project(*proj_axes))
        ]
        sw_per_bin_gen_smear = [hh.divideHists(x, nominal_gen_smear) for x in msv_sw_gen_smear]
        nominal_reco = proc_hists['nominal'].project(*proj_axes)
        msv_sw_reco = [hh.multiplyHists(nominal_reco, x) for x in sw_per_bin_gen_smear]
        proc_hists['muonScaleSyst_responseWeights'] = hh.combineUpDownVarHists(*msv_sw_reco)

def muon_scale_variation_from_manual_shift(
    resultdict, procs = ['WplusmunuPostVFP', 'WminusmunuPostVFP', 'ZmumuPostVFP'],
):
    for proc in procs:
        proc_hists = resultdict[proc]['output']
        manual_shift_hists = [proc_hists['muonScaleVariationDnTenthmil'], proc_hists['muonScaleVariationUpTenthmil']]
        proc_hists['muonScaleSyst_manualShift'] = hh.combineUpDownVarHists(*manual_shift_hists)

def make_alt_reco_and_gen_hists(df, results, nominal_axes):
    nominal_cols_cvh = [
        "goodMuons_eta0_cvh", "goodMuons_pt0_cvh", "goodMuons_charge0_cvh",
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

    nominal_cvh =       df.HistoBoost("nominal_cvh", nominal_axes, [*nominal_cols_cvh, "nominal_weight"])
    nominal_uncrct =      df.HistoBoost("nominal_uncrct", nominal_axes, [*nominal_cols_uncrct, "nominal_weight"])
    nominal_gen =         df.HistoBoost("nominal_gen", nominal_axes, [*nominal_cols_gen, "nominal_weight"])
    nominal_gen_smeared = df.HistoBoost("nominal_gen_smeared", nominal_axes, [*nominal_cols_gen_smeared, "nominal_weight"])

    results.append(nominal_cvh)
    results.append(nominal_uncrct)
    results.append(nominal_gen)
    results.append(nominal_gen_smeared)
    return [nominal_cols_cvh, nominal_cols_uncrct, nominal_cols_gen, nominal_cols_gen_smeared]

####################
## FOR VALIDATION ##
####################

def define_reco_over_gen_cols(df, reco_type, kinematic_vars = ['pt', 'eta']):
    kinematic_vars = common.string_to_list(kinematic_vars)
    for var in kinematic_vars:
        reco_col = f"goodMuons_{var.lower()}0" if reco_type == 'crctd' \
                   else f"goodMuons_{var.lower()}0_{reco_type}"
        df = df.Define(
            f"goodMuons_{var.lower()}0_{reco_type}_over_gen",
            f"{reco_col}/goodMuons_{var.lower()}0_gen"
        )
    return df

def make_reco_over_gen_hists(df, results):
    nominal_cols_crctd_over_gen = [
        "goodMuons_pt0_crctd_over_gen"
    ]
    nominal_cols_cvh_over_gen = [
        "goodMuons_pt0_cvh_over_gen"
    ]
    nominal_cols_uncrct_over_gen = [
        "goodMuons_pt0_uncrct_over_gen"
    ]
    nominal_cols_gen_smeared_over_gen = [
        "goodMuons_pt0_gen_smeared_over_gen",
    ]
    axis_pt_reco_over_gen = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "reco_pt_over_gen")
    axis_qop_reco_over_gen = hist.axis.Regular(1000, 0.9, 1.1, underflow=True, overflow=True, name = "reco_qop_over_gen")
    crctd_over_gen =  df.HistoBoost("crctd_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_crctd_over_gen, "nominal_weight"])
    cvh_over_gen =  df.HistoBoost("cvh_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_cvh_over_gen, "nominal_weight"])
    uncrct_over_gen = df.HistoBoost("uncrct_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_uncrct_over_gen, "nominal_weight"])
    gen_smeared_over_gen = df.HistoBoost("gen_smeared_over_gen", [axis_pt_reco_over_gen], [*nominal_cols_gen_smeared_over_gen, "nominal_weight"])
    
    results.append(crctd_over_gen)
    results.append(cvh_over_gen)
    results.append(uncrct_over_gen)
    results.append(gen_smeared_over_gen)

def define_cols_for_smearing_weights(df, helper_func):
    df = df.Define("unity_as_double", "1.0")
    df = df.Define("muonScaleSyst_smearingWeightsPerSe_tensor", helper_func,
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
        "unity_as_double"
        ]
    )
    df = df.Define("smearing_weights_down", "muonScaleSyst_smearingWeightsPerSe_tensor(0,0)")
    df = df.Define("smearing_weights_up", "muonScaleSyst_smearingWeightsPerSe_tensor(0,1)")
    return df

def make_hists_for_smearing_weights(df, nominal_axes, nominal_cols, results):
    axis_smearing_weight = hist.axis.Regular(1000, 0.99, 1.01, underflow=True, overflow=True, name = "smearing_weight")
    smearing_weights_down = df.HistoBoost("smearing_weights_down", [*nominal_axes, axis_smearing_weight], [*nominal_cols, "smearing_weights_down"])
    smearing_weights_up = df.HistoBoost("smearing_weights_up", [*nominal_axes, axis_smearing_weight], [*nominal_cols, "smearing_weights_up"])
    results.append(smearing_weights_down)
    results.append(smearing_weights_up)
