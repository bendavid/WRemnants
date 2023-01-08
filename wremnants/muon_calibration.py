import ROOT
import pathlib
import hist
import narf
from utilities import common

ROOT.gInterpreter.Declare('#include "muon_calibration.h"')

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

def define_corrected_muons(df, helper, corr_type, dataset):
    if not (dataset.is_data or dataset.name in common.vprocs):
        corr_type = "none" 

    if corr_type == "none":
        df = df.Alias("Muon_correctedPt", "Muon_pt")
        df = df.Alias("Muon_correctedEta", "Muon_eta")
        df = df.Alias("Muon_correctedPhi", "Muon_phi")
        df = df.Alias("Muon_correctedCharge", "Muon_charge")
    elif corr_type == "trackfit_only":
        df = df.Define("Muon_correctedPt", "Muon_cvhPt")
        df = df.Define("Muon_correctedEta", "Muon_cvhEta")
        df = df.Define("Muon_correctedPhi", "Muon_cvhPhi")
        df = df.Define("Muon_correctedCharge", "Muon_cvhCharge")
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
    else:
        raise ValueError(f"Invalid correction type choice {corr_type}")

    return df
