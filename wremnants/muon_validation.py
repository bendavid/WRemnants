import ROOT
import hist
import narf
import numpy as np
import uproot
from functools import reduce
from utilities import common, logging
from wremnants.muon_calibration import get_jpsi_scale_param_cov_mat

narf.clingutils.Declare('#include "muon_validation.h"')

logger = logging.child_logger(__name__)

def make_jpsi_crctn_unc_helper_massweights(
    filepath, n_massweights,
    n_scale_params = 3, n_tot_params = 4, n_eta_bins = 48, scale = 1.0
):
    f = uproot.open(filepath)
    cov = f['covariance_matrix'].to_hist()
    cov_scale_params = get_jpsi_scale_param_cov_mat(cov, n_scale_params, n_tot_params, n_eta_bins, scale)

    w,v = np.linalg.eigh(cov_scale_params)    
    var_mat = np.sqrt(w) * v
    axis_eta = hist.axis.Regular(n_eta_bins, -2.4, 2.4, name = 'eta')
    axis_scale_params = hist.axis.Regular(n_scale_params, 0, 1, name = 'scale_params')
    axis_scale_params_unc = hist.axis.Regular(
        n_eta_bins * n_scale_params, 0, 1,
        underflow = False, overflow = False,  name = 'unc'
    )
    hist_scale_params_unc = hist.Hist(axis_eta, axis_scale_params, axis_scale_params_unc)
    for i in range(n_eta_bins):
        lb, ub = i * n_scale_params, (i + 1) * n_scale_params
        hist_scale_params_unc.view()[i,...] = var_mat[lb:ub][:]
    hist_scale_params_unc_cpp = narf.hist_to_pyroot_boost(hist_scale_params_unc, tensor_rank = 2)
    jpsi_crctn_unc_helper = ROOT.wrem.JpsiCorrectionsUncHelper_massWeights[type(hist_scale_params_unc_cpp).__cpp_name__, n_massweights](
        ROOT.std.move(hist_scale_params_unc_cpp)
    )
    jpsi_crctn_unc_helper.tensor_axes = (hist_scale_params_unc.axes['unc'], common.down_up_axis)
    return jpsi_crctn_unc_helper

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
