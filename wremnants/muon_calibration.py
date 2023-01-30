import ROOT
import pathlib
import hist
import narf
from utilities import rdf_tools
from utilities import common
from utilities import boostHistHelpers as hh
from . import muon_validation
import uproot

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

def make_muon_bias_helpers(corr_type, smearing=True):
    # this helper adds a correction for the bias from nonclosure to the muon pT
    if (smearing and corr_type !="massfit") or corr_type not in ["massfit", "massfit_lbl"]:
        print("The bias correction is not available for this configuration, procees without bias correction helper!")
        return None
    elif smearing:
        filename = "closureZ_smeared"
    elif corr_type == "massfit":
        filename = "closureZ"
    elif corr_type == "massfit_lbl":
        filename = "closureZ_LBL"

    h2d = uproot.open(f"wremnants/data/closure/{filename}.root:closure").to_hist()
    # Drop the uncertainty because the Weight storage type doesn't play nice with ROOT
    h2d_nounc = hist.Hist(*h2d.axes, data=h2d.values(flow=True))
    # Set overflow to closest values
    h2d_nounc[hist.underflow,:][...] = h2d_nounc[0,:].view(flow=True)
    h2d_nounc[:,hist.underflow][...] = h2d_nounc[:,0].view(flow=True)
    h2d_nounc[hist.overflow,:][...] = h2d_nounc[-1,:].view(flow=True)
    h2d_nounc[:,hist.overflow][...] = h2d_nounc[:,-1].view(flow=True)

    h2d_cpp = narf.hist_to_pyroot_boost(h2d_nounc, tensor_rank=0)
    helper = ROOT.wrem.BiasCalibrationHelper[type(h2d_cpp).__cpp_name__](ROOT.std.move(h2d_cpp))

    return helper

def make_muon_smearing_helpers():
    # this helper smears muon pT to match the resolution in data
    h2d = uproot.open(f"wremnants/data/calibration/smearing.root:smearing").to_hist()
    # Drop the uncertainty because the Weight storage type doesn't play nice with ROOT
    h2d_nounc = hist.Hist(*h2d.axes, data=h2d.values(flow=True))
    # Set overflow to closest values
    h2d_nounc[hist.underflow,:][...] = h2d_nounc[0,:].view(flow=True)
    h2d_nounc[:,hist.underflow][...] = h2d_nounc[:,0].view(flow=True)
    h2d_nounc[hist.overflow,:][...] = h2d_nounc[-1,:].view(flow=True)
    h2d_nounc[:,hist.overflow][...] = h2d_nounc[:,-1].view(flow=True)

    h2d_cpp = narf.hist_to_pyroot_boost(h2d_nounc, tensor_rank=0)
    helper = ROOT.wrem.SmearingHelper[type(h2d_cpp).__cpp_name__](ROOT.std.move(h2d_cpp))

    return helper

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

def define_lblcorr_muons(df, cvh_helper, dataset):
    corr_branch = "cvh" if dataset.is_data else "cvhideal"

    # split the nested vectors
    df = df.Define("Muon_cvhmergedGlobalIdxs", "wrem::splitNestedRVec(Muon_cvhmergedGlobalIdxs_Vals, Muon_cvhmergedGlobalIdxs_Counts)")
    df = df.Define(f"Muon_{corr_branch}JacRef", f"wrem::splitNestedRVec(Muon_{corr_branch}JacRef_Vals, Muon_{corr_branch}JacRef_Counts)")

    df = df.Define("Muon_lblMom4Charge", cvh_helper, [f"Muon_{corr_branch}Pt", f"Muon_{corr_branch}Eta", f"Muon_{corr_branch}Phi", f"Muon_{corr_branch}Charge", "Muon_cvhmergedGlobalIdxs", f"Muon_{corr_branch}JacRef"])

    # split into individual vectors
    df = df.Define("Muon_lblPt", "ROOT::VecOps::RVec<float> res(Muon_lblMom4Charge.size()); std::transform(Muon_lblMom4Charge.begin(), Muon_lblMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Pt(); } ); return res;")
    df = df.Define("Muon_lblEta", "ROOT::VecOps::RVec<float> res(Muon_lblMom4Charge.size()); std::transform(Muon_lblMom4Charge.begin(), Muon_lblMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Eta(); } ); return res;")
    df = df.Define("Muon_lblPhi", "ROOT::VecOps::RVec<float> res(Muon_lblMom4Charge.size()); std::transform(Muon_lblMom4Charge.begin(), Muon_lblMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Phi(); } ); return res;")
    df = df.Define("Muon_lblCharge", "ROOT::VecOps::RVec<int> res(Muon_lblMom4Charge.size()); std::transform(Muon_lblMom4Charge.begin(), Muon_lblMom4Charge.end(), res.begin(), [](const auto &x) { return x.second; }); return res;")
    return df


def define_corrected_muons(df, cvh_helper, jpsi_helper, corr_type, dataset, smearing_helper=None, bias_helper=None):
    if not (dataset.is_data or dataset.name in common.vprocs):
        corr_type = "none" 

    muon = "Muon"
    if "lbl" in corr_type:
        df = define_lblcorr_muons(df, cvh_helper, dataset)
        muon = "Muon_lbl"
    elif corr_type != "none":
        fit = "cvhideal" if corr_type == "trackfit_only_mctruth" and not dataset.is_data else "cvh"
        muon = f"Muon_{fit}"

    muon_pt = muon

    if "massfit" in corr_type:
        df = df.Define("Muon_jpsiCorrectedPt", jpsi_helper, [muon_var_name(muon, var) for var in ["pt", "eta", "charge"]])
        muon_pt = "Muon_jpsiCorrected"

    if smearing_helper and not dataset.is_data:
        df = df.Define("Muon_smearedPt", smearing_helper, [muon_var_name(muon_pt, "pt"), muon_var_name(muon, "eta")])
        muon_pt = "Muon_smeared"

    if bias_helper and not dataset.is_data:
        df = df.Define("Muon_biasedPt", bias_helper, [muon_var_name(muon_pt, "pt"), muon_var_name(muon, "eta")])
        muon_pt = "Muon_biased"

    for var in ["pt", "eta", "phi", "charge"]:
        mu_name = muon_pt if var == "pt" else muon
        df = df.Alias(muon_var_name("Muon_corrected", var), muon_var_name(mu_name, var))

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

def muon_var_name(mu_type, var):
        return mu_type+(f"_{var}" if mu_type == "Muon" else var.capitalize())

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

def define_corrected_reco_muon_kinematics(df, muons="goodMuons", kinematic_vars = ["pt", "eta", "phi", "charge"], index=0):
    for var in kinematic_vars:
        df = df.Define(
            f"{muons}_{var.lower()}{index}",
            f"Muon_corrected{var.capitalize()}[{muons}][{index}]"
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
