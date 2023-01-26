import ROOT
import hist
import narf
import numpy as np
import uproot
from functools import reduce
from utilities import common

ROOT.gInterpreter.Declare('#include "muon_validation.h"')

def make_jpsi_crctn_helper(filepath):
    f = uproot.open(filepath)

    # TODO: convert variable axis to regular if the bin width is uniform
    A, e, M = [x.to_hist() for x in [f['A'], f['e'], f['M']]]

    # TODO: make this into a function in utilities/boostHistHelpers
    if (A.axes != e.axes) or (e.axes != M.axes):
        raise RuntimeError("A, e, M histograms have different axes!")
    else:
        axis_param = hist.axis.Regular(3, 0, 3, underflow = False, overflow = False, name = 'param')
        axis_cenval = hist.axis.Regular(1, 0, 1, name = 'cenval')
        hist_comb = hist.Hist(*A.axes, axis_param, axis_cenval, storage = hist.storage.Double())
        hist_comb.view()[...,0] = np.stack([x.values() for x in [A, e, M]], axis = -1)

    hist_comb_cpp = narf.hist_to_pyroot_boost(hist_comb, tensor_rank = 2)
    jpsi_crctn_helper = ROOT.wrem.JpsiCorrectionsRVecHelper[type(hist_comb_cpp).__cpp_name__](
        ROOT.std.move(hist_comb_cpp)
    )
    return jpsi_crctn_helper

def make_jpsi_crctn_unc_helper(filepath, n_scale_params = 3, n_tot_params = 6, n_eta_bins = 48):
    f = uproot.open(filepath)
    cov = f['covariance_matrix'].to_hist()
    cov_scale_params = get_jpsi_scale_param_cov_mat(cov, n_scale_params, n_tot_params, n_eta_bins)

    #
    w,v = np.linalg.eigh(cov_scale_params)    
    var_mat = np.sqrt(w) * v
    axis_eta = hist.axis.Regular(n_eta_bins, 0, 1, name = 'eta')
    axis_scale_params = hist.axis.Regular(n_scale_params, 0, 1, name = 'scale_params')
    axis_scale_params_unc = hist.axis.Regular(n_eta_bins * n_scale_params, 0, 1, name = 'unc')
    hist_scale_params_unc = hist.Hist(axis_eta, axis_scale_params, axis_scale_params_unc)
    for i in range(n_eta_bins):
        lb, ub = i * n_scale_params, (i + 1) * n_scale_params
        hist_scale_params_unc.view()[i,...] = var_mat[lb:ub][:]
    hist_scale_params_unc_cpp = narf.hist_to_pyroot_boost(hist_scale_params_unc, tensor_rank = 2)
    jpsi_crctn_unc_helper = ROOT.wrem.JpsiCorrectionsHelper[type(hist_scale_params_unc_cpp).__cpp_name__](
        ROOT.std.move(hist_scale_params_unc_cpp)
    )
    jpsi_crctn_unc_helper.tensor_axes = (hist_scale_params_unc.axes['unc'], common.down_up_axis)
    return jpsi_crctn_unc_helper

# returns the cov mat of only scale parameters in eta bins, in the form of a 2D numpy array
# there are 3 scale params (A, e, M) + 3 resolution params for each eta bin in the jpsi calib file
def get_jpsi_scale_param_cov_mat(cov, n_scale_params = 3, n_tot_params = 6, n_eta_bins = 48):
    cov_dim = len(cov.axes[0].edges) - 1
    if cov_dim != n_tot_params * n_eta_bins:
        raise ValueError(
            f"dimension of the covariance matrix {cov_dim} doesn't match the input number of "
            f"total parameters {n_tot_params} times the number of eta bins {n_eta_bins}"
        )
    idx_first_param = np.arange(0, cov_dim, n_tot_params)
    idx_scale_params = np.sort(reduce(
        np.append, [(idx_first_param + i) for i in range(n_scale_params)]
    ))
    cov_scale_params = np.empty([n_scale_params * n_eta_bins, n_scale_params * n_eta_bins])
    for irow_scale_params, irow_all_params in enumerate(idx_scale_params):
        cov_scale_params[irow_scale_params, :] = np.take(
            cov.values()[irow_all_params], idx_scale_params
        )
    return cov_scale_params

# "muon" is for mw; "muons" is for wlike, for which we select one of the trig/nonTrig muons

def define_cvh_muon_kinematics(df):
    df = df.Define("goodMuons_cvh_pt0", "Muon_cvhPt[gooodMuons][0]")
    df = df.Define("goodMuons_cvh_eta0", "Muon_cvhEta[goodMuons][0]")
    df = df.Define("goodMuons_cvh_phi0", "Muon_cvhPhi[goodMuons][0]")
    return df

def define_cvh_muons_kinematics(df):
    df = df.Define("trigMuons_cvh_pt0", "Muon_correctedPt[trigMuons][0]")
    df = df.Define("trigMuons_cvh_eta0", "Muon_correctedEta[trigMuons][0]")
    df = df.Define("trigMuons_cvh_phi0", "Muon_correctedPhi[trigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_pt0", "Muon_correctedPt[nonTrigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_eta0", "Muon_correctedEta[nonTrigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_phi0", "Muon_correctedPhi[nonTrigMuons][0]")
    return df

def define_jpsi_crctd_muons_pt_unc(df, helper):
    df = df.Define("trigMuons_jpsi_crctd_pt_unc", helper,
        [
            "trigMuons_cvh_eta",
            "trigMuons_cvh_pt",
            "trigMuons_charge",
            "trigMuons_jpsi_crctd_pt"
        ]
    )
    df = df.Define("nonTrigMuons_jpsi_crctd_pt_unc", helper,
        [
            "nonTrigMuons_cvh_eta",
            "nonTrigMuons_cvh_pt",
            "nonTrigMuons_charge",
            "nonTrigMuons_jpsi_crctd_pt"
        ]
    )
    return df

def define_jpsi_crctd_z_mass(df):
    df = df.Define("trigMuons_jpsi_crctd_mom4",
        (
            "ROOT::Math::PtEtaPhiMVector("
            "trigMuons_jpsi_crctd_pt, trigMuons_cvh_eta, trigMuons_cvh_phi, wrem::muon_mass)"
        )
    )
    df = df.Define("nonTrigMuons_jpsi_crctd_mom4",
        (
            "ROOT::Math::PtEtaPhiMVector("
            "nonTrigMuons_jpsi_crctd_pt, nonTrigMuons_cvh_eta, nonTrigMuons_cvh_phi, wrem::muon_mass)"
        )
    )
    df = df.Define("Z_jpsi_crctd_mom4", "ROOT::Math::PxPyPzEVector(trigMuons_jpsi_crctd_mom4)+ROOT::Math::PxPyPzEVector(nonTrigMuons_jpsi_crctd_mom4)")
    df = df.Define("massZ_jpsi_crctd", "Z_jpsi_crctd_mom4.mass()")
    return df

def define_jpsi_crctd_unc_z_mass(df):
    df = df.Define("trigMuons_jpsi_crctd_mom4_unc",
        (
            "ROOT::VecOps::RVec<double> res(trigMuons_jpsi_crctd_pt_unc.size());"
            "for (int i = 0; i < trigMuons_jpsi_crctd_pt_unc.size(); i++) {"
            "    res[i] = ("
            "       ROOT::Math::PtEtaPhiMVector("
            "           trigMuons_jpsi_crctd_pt_unc[i],"
            "           trigMuons_cvh_eta,"
            "           trigMuons_cvh_phi,"
            "           wrem::muon_mass"
            "       )"
            "    );"
            "}"
            "return res;"
        )
    )
    df = df.Define("nonTrigMuons_jpsi_crctd_mom4_unc",
        (
            "ROOT::VecOps::RVec<double> res(nonTrigMuons_jpsi_crctd_pt_unc.size());"
            "for (int i = 0; i < nonTrigMuons_jpsi_crctd_pt_unc.size(); i++) {"
            "    res[i] = ("
            "        ROOT::Math::PtEtaPhiMVector("
            "            nonTrigMuons_jpsi_crctd_pt_unc[i],"
            "            nonTrigMuons_cvh_eta," 
            "            nonTrigMuons_cvh_phi,"
            "            wrem::muon_mass"
            "        )"
            "    );"
            "}"
            "return res;"
        )
    )
    df = df.Define("Z_jpsi_crctd_mom4_unc", 
        (
            "ROOT::VecOps::RVec<double> res(trigMuons_jpsi_crctd_mom4_unc.size());"
            "for (int i = 0; i < trigMuons_jpsi_crctd_mom4_unc.size(); i++) {"
            "    res[i] = ("
            "        ROOT::Math::PxPyPzEVector(trigMuons_jpsi_crctd_mom4_unc[i]) +"
            "        ROOT::Math::PxPyPzEVector(nonTrigMuons_jpsi_crctd_mom4_unc[i])"
            "    )"
            "}"
            "return res"
        )
    )
    df = df.Define("massZ_jpsi_crctd_unc", 
        (
            "ROOT::VecOps::RVec<double> res(Z_jpsi_crctd_mom4_unc.size());"
            "for (int i = 0; i < Z_jpsi_crctd_mom4_unc.size(); i++) {"
            "    res[i] = Z_jpsi_crctd_mom4_unc[i].mass()"
            "}"
            "return res"
        )
    )
    return df
