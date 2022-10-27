import ROOT
import hist
import narf
import numpy as np
import uproot

ROOT.gInterpreter.Declare('#include "muon_validation.h"')

def make_jpsi_crctn_helper(filepath = "/scratch/shared/MuonCalibration/calibrationJDATA_aftersm.root"):
    f = uproot.open(filepath)

    # TODO: convert variable axis to regular if the bin width is uniform
    A, e, M = [x.to_hist() for x in [f['A'], f['e'], f['M']]]

    # TODO: make this into a function in utilities/boostHistHelpers
    if (A.axes != e.axes) or (e.axes != M.axes):
        raise RuntimeError("A, e, M histograms have different axes!")
    else:
        axis_param = hist.axis.Regular(3, 0, 3, underflow = False, overflow = False, name = 'param')
        hist_comb = hist.Hist(*A.axes, axis_param, storage = A._storage_type())
        hist_comb.view(flow=True)[...] = np.stack([x.view(flow=True) for x in [A, e, M]], axis = -1)

    hist_comb_cpp = narf.hist_to_pyroot_boost(hist_comb, tensor_rank = 1)
    jpsi_crctn_helper = ROOT.wrem.JpsiCorrectionsHelper[type(hist_comb_cpp).__cpp_name__](ROOT.std.move(hist_comb_cpp))
    print("jpsi_helper is", jpsi_crctn_helper)
    return jpsi_crctn_helper

def define_cvh_muons_kinematics(df):
    df = df.Define("TrigMuon_cvh_pt", "Muon_cvhPt[trigMuons][0]")
    df = df.Define("TrigMuon_cvh_eta", "Muon_cvhEta[trigMuons][0]")
    df = df.Define("TrigMuon_cvh_phi", "Muon_cvhPhi[trigMuons][0]")
    df = df.Define("NonTrigMuon_cvh_pt", "Muon_cvhPt[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_cvh_eta", "Muon_cvhEta[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_cvh_phi", "Muon_cvhPhi[nonTrigMuons][0]")
    return df

def define_jpsi_crctd_muons_pt(df, helper):
    df = df.Define("TrigMuon_jpsi_crctd_pt", helper,
        [
            "TrigMuon_cvh_eta",
            "TrigMuon_cvh_pt",
            "TrigMuon_charge"
        ]
    )
    df = df.Define("NonTrigMuon_jpsi_crctd_pt", helper,
        [
            "NonTrigMuon_cvh_eta",
            "NonTrigMuon_cvh_pt",
            "NonTrigMuon_charge"
        ]
    )
    return df
