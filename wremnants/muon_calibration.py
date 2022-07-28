import ROOT
import pathlib
import hist
import narf
from utilities import rdf_tools
from utilities import boostHistHelpers as hh

ROOT.gInterpreter.Declare('#include "muon_calibration.h"')
ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_calibration_helpers(filename = data_dir + "/calibration/correctionResults_v207_gen.root", era = None):

    helper = ROOT.wrem.CVHCorrector(filename)

    uncertainty_hist = get_dummy_uncertainties()
    uncertainty_hist_cpp = narf.hist_to_pyroot_boost(uncertainty_hist, tensor_rank = 2)
    # min gen pt = 9 GeV to avoid threshold effects
    # max weight = 10 to protect against outliers
    uncertainty_helper = ROOT.wrem.calibration_uncertainty_helper[type(uncertainty_hist_cpp)](ROOT.std.move(uncertainty_hist_cpp), 9., 10.)

    down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
    uncertainty_helper.tensor_axes = (uncertainty_hist.axes["calvar"], down_up_axis)

    return helper, uncertainty_helper

'''
def get_dummy_uncertainties():
    axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
    axis_calvar = hist.axis.Integer(0, 6, underflow=False, overflow=False, name = "calvar")
    axis_calparm = hist.axis.Integer(0, 4, underflow=False, overflow=False, name = "calparm")

    h = hist.Hist(axis_eta, axis_calvar, axis_calparm)

    # forward b-field-like
    h.values()[:axis_eta.index(-1.2), 0, 0] = 1e-4

    # forward alignment-like
    h.values()[:axis_eta.index(-1.2), 1, 2] = 1e-4/40.

    # forward resolution
    h.values()[:axis_eta.index(-1.2), 2, 3] = 0.05

    # central b-field-like
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 3, 0] = 1e-4

    # central alignment-like
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 4, 2] = 1e-4/40.

    # central resolution
    h.values()[axis_eta.index(-1.2):axis_eta.index(0.), 5, 3] = 0.05

    #mirror to positive eta
    h.values()[axis_eta.index(0.):, ...] = h.values()[:axis_eta.index(0.), ...]

    # set underflow and overflow to match boundaries
    h.values(flow=True)[0, ...] = h.values(flow=True)[1, ...]
    h.values(flow=True)[-1, ...] = h.values(flow=True)[-2, ...]

    return h
'''

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

def define_corrected_muons(df, helper):
    # split the nested vectors
    df = df.Define("Muon_cvhmergedGlobalIdxs", "wrem::splitNestedRVec(Muon_cvhmergedGlobalIdxs_Vals, Muon_cvhmergedGlobalIdxs_Counts)")
    df = df.Define("Muon_cvhbsJacRef", "wrem::splitNestedRVec(Muon_cvhbsJacRef_Vals, Muon_cvhbsJacRef_Counts)")

    df = df.Define("Muon_correctedMom4Charge", helper, ["Muon_cvhbsPt", "Muon_cvhbsEta", "Muon_cvhbsPhi", "Muon_cvhbsCharge", "Muon_cvhmergedGlobalIdxs", "Muon_cvhbsJacRef"])

    # split into individual vectors
    df = df.Define("Muon_correctedPt", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Pt(); } ); return res;")
    df = df.Define("Muon_correctedEta", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Eta(); } ); return res;")
    df = df.Define("Muon_correctedPhi", "ROOT::VecOps::RVec<float> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.first.Phi(); } ); return res;")
    df = df.Define("Muon_correctedCharge", "ROOT::VecOps::RVec<int> res(Muon_correctedMom4Charge.size()); std::transform(Muon_correctedMom4Charge.begin(), Muon_correctedMom4Charge.end(), res.begin(), [](const auto &x) { return x.second; }); return res;")
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

def define_corrected_reco_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        df = df.Define(
            f"goodMuons_{var.lower()}0",
            f"Muon_corrected{var.capitalize()}[goodMuons][0]"
        )
    return df

def define_cvhbs_reco_muon_kinematics(df, kinematic_vars = ["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        df = df.Define(
            f"goodMuons_{var.lower()}0_cvhbs",
            f"Muon_cvhbs{var.capitalize()}[goodMuons][0]"
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
    proj_axes = ['pt', 'eta']
):
    for proc in procs:
        proc_hists = resultdict[proc]['output']
        if 'plus' in proc.lower(): charge = 1
        elif 'minus' in proc.lower(): charge = -1
        nominal_gen_smear = (
            proc_hists['nominal_gen_smeared'][{'charge':hist.loc(charge)}].project(*proj_axes))
        msv_sw_gen_smear = [
            (proc_hists['muonScaleSyst_responseWeights_gensmear']
                       [{'charge':hist.loc(charge)}][...,0,0].project(*proj_axes)),
            (proc_hists['muonScaleSyst_responseWeights_gensmear']
                       [{'charge':hist.loc(charge)}][...,0,1].project(*proj_axes))
        ]
        sw_per_bin_gen_smear = [hh.divideHists(x, nominal_gen_smear) for x in msv_sw_gen_smear]
        nominal_reco = proc_hists['nominal'][{'charge':hist.loc(charge)}].project(*proj_axes)
        msv_sw_reco = [hh.multiplyHists(nominal_reco, x) for x in sw_per_bin_gen_smear]
        proc_hists['muonScaleSyst_responseWeights'] = hh.combineUpDownVarHists(*msv_sw_reco)


