import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=None)
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None)
args = parser.parse_args()

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
if not args.nThreads:
    ROOT.ROOT.EnableImplicitMT()
elif args.nThreads != 1:
    ROOT.ROOT.EnableImplicitMT(args.nThreads)

print(ROOT.GetThreadPoolSize())
#ROOT.DisableImplicitMT() # need to disable it first in order to correctly set the number of threads
#ROOT.EnableImplicitMT(1)
print(ROOT.GetThreadPoolSize())

import pickle
import gzip

import narf
import wremnants
from wremnants import theory_tools
import hist
import lz4.frame
import logging
import math
import sys

import decimal
def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)

###################################
flavor = args.flavor # mumu, ee

if flavor == "mumu":

    procs = ["DYmumu", "singlemuon", "DYee", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]
    #procs = ["DYmumu"]
    
elif flavor == "ee":

    procs = ["DYee", "singleelectron", "DYmumu", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]
    #procs = ["DYee"]
    
else: sys.exit("Flavor must be ee or mumu")

filt = lambda x,filts=procs: any([f in x.name for f in filts]) 
datasets = wremnants.datasetsLowPU.getDatasets(maxFiles=args.maxFiles, filt=filt)


# load lowPU specific libs
#ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")
ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')
ROOT.gInterpreter.Declare('#include "lowpu_efficiencies.h"')
ROOT.gInterpreter.Declare('#include "lowpu_prefire.h"')
ROOT.gInterpreter.Declare('#include "lowpu_rochester.h"')
ROOT.gInterpreter.Declare('#include "lowpu_recoil.h"')

# load recoil hists
ROOT.wrem.recoil_init()


# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_mll = hist.axis.Regular(60, 60., 120., name = "mll")
axis_yll = hist.axis.Regular(50, -2.5, 2.5, name = "yll")
axis_ptl = hist.axis.Regular(100, 0., 200., name = "ptl")
axis_etal = hist.axis.Regular(50, -2.5, 2.5, name = "etal")

# recoil/MET axes
axis_MET_pt = hist.axis.Regular(300, 0, 300, name = "MET_pt")
axis_MET_phi = hist.axis.Regular(50, 0, 4, name = "MET_phi")
axis_recoil_magn = hist.axis.Regular(300, 0, 300, name = "recoil_magn")
axis_recoil_para = hist.axis.Regular(1000, -500, 500, name = "recoil_para")
axis_recoil_perp = hist.axis.Regular(1000, -500, 500, name = "recoil_perp")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge


nominal_axes = [axis_eta, axis_pt, axis_charge]
axis_lin = hist.axis.Regular(5, 0, 5, name = "lin")


axis_ptll = hist.axis.Variable([0, 2, 3, 4, 4.75, 5.5, 6.5, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100], name = "ptll")
axis_costhetastarll = hist.axis.Regular(20, -1., 1., name = "costhetastarll")
axis_phistarll = hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phistarll")


qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
axis_ptVgen = qcdScaleByHelicity_helper.hist.axes["ptVgen"]
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]


# unfolding axes
axis_recoil_reco = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 10000], name = "recoil_reco", underflow=False)
axis_recoil_gen = hist.axis.Variable([0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 10000], name = "recoil_gen", underflow=False)

axis_recoil_reco_pert = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 10000], name = "recoil_reco_pert", underflow=False)

# axes for final cards/fitting
reco_mll_axes = [axis_recoil_reco, axis_mll]
gen_reco_mll_axes = [axis_recoil_gen, axis_recoil_reco, axis_mll]

# extra axes which can be used to label tensor_axes
down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")


qTbins = list(range(0, 30, 1)) + list(range(30, 60, 2)) + list(range(60, 100, 5)) + list(range(100, 210, 10)) + [10000]
qTbins = list(drange(0, 30, 0.5)) + list(range(30, 60, 2)) + list(range(60, 100, 5)) + list(range(100, 210, 10)) + [10000]
axis_qT = hist.axis.Variable(qTbins, name = "qT")


# set the recoil correction bins in the analyzer
qTbins_vec = ROOT.std.vector["float"]()
for v in qTbins: qTbins_vec.push_back(v)
setattr(ROOT.wrem, "qTbins", qTbins_vec)


#axis_recoil_stat_unc = hist.axis.Integer(0, len(qTbins)*11, underflow=False, overflow=False, name = "recoil_stat_unc_var")
#axis_recoil_stat_unc = hist.axis.Regular(len(qTbins)*11, 0, len(qTbins)*11, underflow=False, overflow=False, name = "recoil_stat_unc_var")
axis_recoil_stat_unc = hist.axis.Regular(len(qTbins), 0, len(qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")

def build_graph(df, dataset):

    print("build graph")
    results = []

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")
    
    if flavor == "mumu":
    
        if not dataset.is_data: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterMC(Muon_pt, Muon_eta, Muon_phi, Muon_charge, Muon_genPartIdx, GenPart_pt, Muon_nTrackerLayers)")
            df = df.Filter("HLT_Mu17")
            
        else: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterData(Muon_pt, Muon_eta, Muon_phi, Muon_charge)")
            df = df.Filter("HLT_HIMu17")
            
        
        df = df.Define("vetoMuons", "Muon_pt_corr > 10 && Muon_looseId && abs(Muon_eta) < 2.4 && abs(Muon_dxybs) < 0.05")
        df = df.Filter("Sum(vetoMuons) == 2")
        
        df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4")
        df = df.Filter("Sum(vetoElectrons) == 0")
        
        df = df.Define("goodLeptons", "vetoMuons && Muon_pt_corr > 25 && Muon_mediumId && Muon_pfRelIso04_all < 0.15")
        df = df.Define("goodLeptonsPlus", "goodLeptons && Muon_charge > 0")
        df = df.Define("goodLeptonsMinus", "goodLeptons && Muon_charge < 0")
        df = df.Filter("Sum(goodLeptons) == 2")
        df = df.Filter("(Muon_charge[goodLeptons][0] + Muon_charge[goodLeptons][1]) == 0")
        
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidateLowPU(TrigObj_id, TrigObj_pt, TrigObj_l1pt, TrigObj_l2pt, TrigObj_filterBits)")
        df = df.Define("trigMatch", "wrem::hasTriggerMatchLowPU(Muon_eta[goodLeptons], Muon_phi[goodLeptons], TrigObj_eta[goodTrigObjs], TrigObj_phi[goodTrigObjs])")
        df = df.Define("nonTrigMatch", "wrem::inverse(trigMatch)")
        df = df.Filter("Sum(trigMatch) > 0")

        df = df.Define("Lep_pt_uncorr", "Muon_pt[goodLeptons]")
        df = df.Define("Lep_pt", "Muon_pt_corr[goodLeptons]")
        df = df.Define("Lep_eta", "Muon_eta[goodLeptons]")
        df = df.Define("Lep_phi", "Muon_phi[goodLeptons]")
        df = df.Define("Lep_charge", "Muon_charge[goodLeptons]")
        df = df.Define("Lep_mass", "Muon_mass[goodLeptons]")
        
        if not dataset.is_data:
            df = df.Define("lepSF_ISO", "wrem::lepSF(Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons], 1)")
            df = df.Define("lepSF_IDIP", "wrem::lepSF(Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons], 2)") # largest effect
            df = df.Define("lepSF_HLT", "wrem::lepSF_HLT_q(Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons], 13)")
            df = df.Define("prefireCorr", "wrem::prefireCorr(0, Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_phi[goodLeptons])")
            df = df.Define("SFMC", "lepSF_IDIP*lepSF_ISO*lepSF_HLT*prefireCorr")
        
        else: df = df.Define("SFMC", "1.0")
    
    else:

        if not dataset.is_data: 
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(0, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 0)")
            df = df.Filter("HLT_Ele20_WPLoose_Gsf")
            
        else: 
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(1, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 0)")
            df = df.Filter("HLT_HIEle20_WPLoose_Gsf")
            
            
        # by default the E/gamma smearings are applied, but are they propagated to the MET? 
        #df = df.Define("Electron_pt_uncorr", "wrem::Egamma_undoCorrection(Electron_pt, Electron_eta, Electron_ecalCorr)")    
        #df = df.Alias("Electron_pt_corr", "Electron_pt")
        
        df = df.Define("vetoElectrons", "Electron_pt_corr > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4")
        df = df.Filter("Sum(vetoElectrons)==2")
        
        df = df.Define("vetoMuons", "Muon_pt > 10 && Muon_looseId && abs(Muon_eta) < 2.4 && abs(Muon_dxybs) < 0.05 && abs(Muon_dz)< 0.2")
        df = df.Filter("Sum(vetoMuons) == 0")
        
        df = df.Define("goodLeptons", "vetoElectrons && Electron_pt_corr > 25 && Electron_cutBased >= 3 && !(abs(Electron_eta) > 1.4442 && abs(Electron_eta) < 1.566)")
        df = df.Define("goodLeptonsPlus", "goodLeptons && Electron_charge > 0")
        df = df.Define("goodLeptonsMinus", "goodLeptons && Electron_charge < 0")
        df = df.Filter("Sum(goodLeptons) == 2")
        
        df = df.Filter("(Electron_charge[goodLeptons][0] + Electron_charge[goodLeptons][1]) == 0")
        df = df.Define("goodTrigObjs", "wrem::goodElectronTriggerCandidateLowPU(TrigObj_id, TrigObj_pt, TrigObj_l1pt, TrigObj_l2pt, TrigObj_filterBits)")
        df = df.Define("trigMatch", "wrem::hasTriggerMatchLowPU(Electron_eta[goodLeptons], Electron_phi[goodLeptons], TrigObj_eta[goodTrigObjs], TrigObj_phi[goodTrigObjs])")
        df = df.Define("nonTrigMatch", "wrem::inverse(trigMatch)")
        df = df.Filter("Sum(trigMatch) > 0")
        
        df = df.Define("Lep_pt_uncorr", "Electron_pt[goodLeptons]")
        df = df.Define("Lep_pt", "Electron_pt_corr[goodLeptons]")
        df = df.Define("Lep_eta", "Electron_eta[goodLeptons]")
        df = df.Define("Lep_phi", "Electron_phi[goodLeptons]")
        df = df.Define("Lep_charge", "Electron_charge[goodLeptons]")
        df = df.Define("Lep_mass", "Electron_mass[goodLeptons]")
        
        if not dataset.is_data:
            df = df.Define("lepSF_IDISO", "wrem::lepSF(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 3)")
            df = df.Define("lepSF_HLT", "wrem::lepSF_HLT_q(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 11)")
            df = df.Define("prefireCorr", "wrem::prefireCorr(0, Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_phi[goodLeptons])")
            df = df.Define("SFMC", "lepSF_IDISO*lepSF_HLT*prefireCorr")
        
        else: df = df.Define("SFMC", "1.0")    
    
   

    df = df.Filter("Lep_pt[0] > 25")
    df = df.Filter("Lep_pt[1] > 25")
        
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")

    df = df.Define("Lep1_mom4", "ROOT::Math::PtEtaPhiMVector(Lep_pt[0], Lep_eta[0], Lep_phi[0], Lep_mass[0])")
    df = df.Define("Lep2_mom4", "ROOT::Math::PtEtaPhiMVector(Lep_pt[1], Lep_eta[1], Lep_phi[1], Lep_mass[0])")
    df = df.Define("Z_mom4", "ROOT::Math::PxPyPzEVector(Lep1_mom4) + ROOT::Math::PxPyPzEVector(Lep2_mom4)")
    df = df.Define("ptZ", "Z_mom4.pt()")
    df = df.Define("massZ", "Z_mom4.mass()")
    df = df.Define("yZ", "Z_mom4.Rapidity()")
    df = df.Define("absYZ", "std::fabs(yZ)")
    df = df.Filter("massZ > 60 && massZ < 120")
    
    if not dataset.is_data:
        df = df.Define("nominal_weight", "weight*SFMC")
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")


    results.append(df.HistoBoost("massZ", [axis_mll], ["massZ", "nominal_weight"]))
    results.append(df.HistoBoost("yZ", [axis_yll], ["yZ", "nominal_weight"]))
    results.append(df.HistoBoost("ptZ", [axis_ptll], ["ptZ", "nominal_weight"]))
    
    
    # lepton kinematics: leading/subleading, triggered/nontriggered
    df = df.Define("Lep_pt_leading", "Lep_pt[0] > Lep_pt[1] ? Lep_pt[0] : Lep_pt[1]")
    df = df.Define("Lep_pt_subleading", "Lep_pt[0] > Lep_pt[1] ? Lep_pt[1] : Lep_pt[0]")
    df = df.Define("Lep_pt_plus", "Lep_charge[0] > 0 ? Lep_pt[0] : Lep_pt[1]")
    df = df.Define("Lep_pt_minus", "Lep_charge[0] > 0? Lep_pt[1] : Lep_pt[0]")
    df = df.Define("Lep_pt_trg", "Lep_pt[trigMatch]")    
    df = df.Define("Lep_pt_nontrg", "Lep_pt[nonTrigMatch]")
    
    df = df.Define("Lep_eta_leading", "Lep_pt[0] > Lep_pt[1] ? Lep_eta[0] : Lep_eta[1]")
    df = df.Define("Lep_eta_subleading", "Lep_pt[0] > Lep_pt[1] ? Lep_eta[1] : Lep_eta[0]")
    df = df.Define("Lep_eta_plus", "Lep_charge[0] > 0 ? Lep_eta[0] : Lep_eta[1]")
    df = df.Define("Lep_eta_minus", "Lep_charge[0] > 0? Lep_eta[1] : Lep_eta[0]")
    df = df.Define("Lep_eta_trg", "Lep_eta[trigMatch]")
    df = df.Define("Lep_eta_nontrg", "Lep_eta[nonTrigMatch]")
    
    results.append(df.HistoBoost("lep_pt", [axis_ptl], ["Lep_pt", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_leading", [axis_ptl], ["Lep_pt_leading", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_subleading", [axis_ptl], ["Lep_pt_subleading", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_plus", [axis_ptl], ["Lep_pt_plus", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_minus", [axis_ptl], ["Lep_pt_minus", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_trg", [axis_ptl], ["Lep_pt_trg", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt_nontrg", [axis_ptl], ["Lep_pt_nontrg", "nominal_weight"]))

    results.append(df.HistoBoost("lep_eta", [axis_etal], ["Lep_eta", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_leading", [axis_etal], ["Lep_eta_leading", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_subleading", [axis_etal], ["Lep_eta_subleading", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_plus", [axis_etal], ["Lep_eta_plus", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_minus", [axis_etal], ["Lep_eta_minus", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_trg", [axis_etal], ["Lep_eta_trg", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta_nontrg", [axis_etal], ["Lep_eta_nontrg", "nominal_weight"]))

    df = df.Define("noTrigMatch", "Sum(trigMatch)")
    results.append(df.HistoBoost("noTrigMatch", [axis_lin], ["noTrigMatch", "nominal_weight"]))


    # Recoil calibrations
    df = df.Define("Lep1_mom2", "ROOT::Math::Polar2DVectorD(Lep_pt[0], Lep_phi[0])")
    df = df.Define("Lep2_mom2", "ROOT::Math::Polar2DVectorD(Lep_pt[1], Lep_phi[1])")
    df = df.Define("Z_mom2", "Lep1_mom2 + Lep2_mom2") # 2D vector sum of both leptons
    df = df.Define("qT", "Z_mom2.R()")
    
    df = df.Alias("MET_uncorr_pt", "DeepMETResolutionTune_pt")
    df = df.Alias("MET_uncorr_phi", "DeepMETResolutionTune_phi")

    df = df.Define("MET_corr_lep", "wrem::METLeptonCorrection(MET_uncorr_pt, MET_uncorr_phi, Lep_pt_uncorr, Lep_pt, Lep_phi)")
    df = df.Define("MET_corr_lep_pt", "MET_corr_lep[0]")
    df = df.Define("MET_corr_lep_phi", "MET_corr_lep[1]")
    
    df = df.Define("recoil_uncorr", "wrem::recoilComponents(MET_corr_lep_pt, MET_corr_lep_phi, qT, Z_mom2.Phi())")
    df = df.Define("recoil_uncorr_magn", "recoil_uncorr[0]")
    df = df.Define("recoil_uncorr_para", "recoil_uncorr[1]")
    df = df.Define("recoil_uncorr_perp", "recoil_uncorr[2]")
    
    
    results.append(df.HistoBoost("MET_uncorr_pt", [axis_MET_pt], ["MET_uncorr_pt", "nominal_weight"]))
    results.append(df.HistoBoost("MET_uncorr_phi", [axis_MET_phi], ["MET_uncorr_phi", "nominal_weight"]))
    results.append(df.HistoBoost("MET_corr_lep_pt", [axis_MET_pt], ["MET_corr_lep_pt", "nominal_weight"]))
    results.append(df.HistoBoost("MET_corr_lep_phi", [axis_MET_phi], ["MET_corr_lep_phi", "nominal_weight"]))
    
    results.append(df.HistoBoost("recoil_uncorr_magn", [axis_recoil_magn], ["recoil_uncorr_magn", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_uncorr_para", [axis_recoil_para], ["recoil_uncorr_para", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_uncorr_perp", [axis_recoil_perp], ["recoil_uncorr_perp", "nominal_weight"]))
    
    results.append(df.HistoBoost("recoil_uncorr_magn_qt", [axis_qT, axis_recoil_magn], ["qT", "recoil_uncorr_magn", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_uncorr_para_qt", [axis_qT, axis_recoil_para], ["qT", "recoil_uncorr_para", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_uncorr_perp_qt", [axis_qT, axis_recoil_perp], ["qT", "recoil_uncorr_perp", "nominal_weight"]))

    results.append(df.HistoBoost("qT", [axis_qT], ["qT", "nominal_weight"]))
    
    # apply recoil correction (only for DY samples)
    if dataset.name == "DYmumu" or dataset.name == "DYee":
    
        df = df.Define("qTbin", "wrem::getqTbin(qT)")
        df = df.Define("recoil_corr", "wrem::recoilCorrectionBinned(recoil_uncorr_para, recoil_uncorr_perp, qTbin)")
        #df = df.Define("recoil_corr", "wrem::recoilCorrectionParametric(recoil_uncorr_para, recoil_uncorr_perp, qT)")
        df = df.Define("recoil_corr_magn", "recoil_corr[0]")
        df = df.Define("recoil_corr_para", "recoil_corr[1]")
        df = df.Define("recoil_corr_para_qT", "recoil_corr[1] + qT")
        df = df.Define("recoil_corr_perp", "recoil_corr[2]")
        
        df = df.Define("MET_corr_rec", "wrem::METCorrection(MET_corr_lep_pt, MET_corr_lep_phi, recoil_corr_para, recoil_corr_perp, qT, Z_mom2.Phi()) ")
        df = df.Define("MET_corr_rec_pt", "MET_corr_rec[0]")
        df = df.Define("MET_corr_rec_phi", "MET_corr_rec[1]")

    else:
    
        df = df.Alias("recoil_corr_magn", "recoil_uncorr_magn")
        df = df.Alias("recoil_corr_para", "recoil_uncorr_para")
        df = df.Define("recoil_corr_para_qT", "recoil_uncorr_para + qT")
        df = df.Alias("recoil_corr_perp", "recoil_uncorr_perp")
        df = df.Alias("MET_corr_rec_pt", "MET_corr_lep_pt")
        df = df.Alias("MET_corr_rec_phi", "MET_corr_lep_phi")
       

    results.append(df.HistoBoost("recoil_corr_magn", [axis_recoil_magn], ["recoil_corr_magn", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_corr_para", [axis_recoil_para], ["recoil_corr_para", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_corr_para_qT", [axis_recoil_perp], ["recoil_corr_para_qT", "nominal_weight"]))
    results.append(df.HistoBoost("recoil_corr_perp", [axis_recoil_perp], ["recoil_corr_perp", "nominal_weight"]))
    results.append(df.HistoBoost("MET_corr_rec_pt", [axis_MET_pt], ["MET_corr_rec_pt", "nominal_weight"]))

    reco_mll_cols = ["recoil_corr_magn", "massZ"]
    results.append(df.HistoBoost("reco_mll", reco_mll_axes, [*reco_mll_cols, "nominal_weight"]))
    
    # 
    
    
    if dataset.name == "DYmumu" or dataset.name == "DYee":
    
        df = wremnants.define_prefsr_vars(df)
    
        gen_reco_mll_cols = ["ptVgen", "recoil_corr_magn", "massZ"]
        results.append(df.HistoBoost("gen_reco_mll", gen_reco_mll_axes, [*gen_reco_mll_cols, "nominal_weight"]))
        
        
        #nweights = 21 if isW else 23
        nweights = 23
        df = df.Define("massWeight_tensor", f"auto res = wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight); res = nominal_weight*res; return res;")
        df = df.Define("massWeight_tensor_wnom", f"auto res = massWeight_tensor; res = nominal_weight*res; return res;")
        massWeight = df.HistoBoost("massWeight", gen_reco_mll_axes, [*nominal_cols, "massWeight_tensor_wnom"])
        results.append(massWeight)
    
        # pdfs
        results.extend(theory_tools.define_and_make_pdf_hists(df, gen_reco_mll_axes, gen_reco_mll_cols, hname="gen_reco_mll"))
        
        #df = theory_tools.define_scale_tensor(df)
        #results.append(theory_tools.make_scale_hist(df, gen_reco_mll_axes, gen_reco_mll_cols, hname="gen_reco_mll"))
        
        df = theory_tools.define_scale_tensor(df)
        df = df.Define("helicityWeight_tensor", qcdScaleByHelicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
        qcdScaleByHelicityUnc = df.HistoBoost("gen_reco_mll_qcdScaleByHelicity", gen_reco_mll_axes+[axis_ptVgen, axis_chargeVgen], [*gen_reco_mll_cols, "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=qcdScaleByHelicity_helper.tensor_axes)
        results.append(qcdScaleByHelicityUnc)


        # recoil stat uncertainty
        df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(qTbins)))
        recoil_vars = [(1,2), (1,3), (1,4),   (2,2), (2,3),   (3,2), (3,3), (3,4),    (4,2), (4,3)]
        for k in recoil_vars:
            df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_StatUnc(recoil_uncorr_para, recoil_uncorr_perp, qTbin, %d, %d)" % (k[0], k[1]))
            results.append(df.HistoBoost("gen_reco_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_magn", "recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "nominal_weight"]))
        
        
        # lepton efficiencies
            
        '''
        lepSF_HLT_DATA_stat_syst    : lepSF_HLT_var_mu(1, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_HLT_DATA_syst_syst    : lepSF_HLT_var_mu(2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_HLT_MC_stat_syst      : lepSF_HLT_var_mu(-1, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_HLT_MC_syst_syst      : lepSF_HLT_var_mu(-2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly

        lepSF_ISO_stat_syst         : lepSF_ISO_var_mu(1, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_ISO_DATA_syst_syst    : lepSF_ISO_var_mu(2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_ISO_MC_syst_syst      : lepSF_ISO_var_mu(-2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly

        lepSF_IDIP_stat_syst        : lepSF_IDIP_var_mu(1, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_IDIP_DATA_syst_syst   : lepSF_IDIP_var_mu(2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        lepSF_IDIP_MC_syst_syst     : lepSF_IDIP_var_mu(-2, Muon_pt_corr[goodLeptons], Muon_eta[goodLeptons], Muon_charge[goodLeptons]): mcOnly
        '''
                
                
        
        # prefire
        df = df.Define("prefireCorr_syst", "wrem::prefireCorr_syst(Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Lep_pt, Lep_eta, Lep_phi)")
        df = df.Define("prefireCorr_syst_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; auto w = nominal_weight*prefireCorr_syst; std::copy(std::begin(w), std::end(w), res.data()); return res;")
        results.append(df.HistoBoost("gen_reco_mll_prefireCorr", [*gen_reco_mll_axes], [*gen_reco_mll_cols, "prefireCorr_syst_tensor"], tensor_axes = [down_up_axis]))

        
    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

fname = "mz_lowPU_%s.pkl.lz4" % flavor

print("writing output")
with lz4.frame.open(fname, "wb") as f:
    pickle.dump(resultdict, f, protocol = pickle.HIGHEST_PROTOCOL)
