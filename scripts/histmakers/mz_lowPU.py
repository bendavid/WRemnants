import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=None)
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None)
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
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
from wremnants import theory_tools,output_tools
import hist
import lz4.frame
import logging
import math
import sys

import scripts.lowPU.config as lowPUcfg


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




# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_mll = hist.axis.Regular(60, 60., 120., underflow=False, overflow=False, name = "mll")
axis_yll = hist.axis.Regular(50, -2.5, 2.5, name = "yll")
axis_ptl = hist.axis.Regular(100, 0., 200., name = "ptl")
axis_etal = hist.axis.Regular(50, -2.5, 2.5, name = "etal")


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
axis_recoil_reco = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 150], name = "recoil_reco", underflow=False, overflow=True)
axis_recoil_gen = hist.axis.Variable([0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 150], name = "recoil_gen", underflow=False, overflow=True)


# axes for final cards/fitting
reco_mll_axes = [axis_recoil_reco, axis_mll]
gen_reco_mll_axes = [axis_recoil_gen, axis_recoil_reco, axis_mll]
axis_mt = hist.axis.Regular(200, 0., 200., name = "mt")

# extra axes which can be used to label tensor_axes
down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")



# recoil initialization
from wremnants import recoil_tools
recoilHelper = recoil_tools.Recoil("lowPU")



def build_graph(df, dataset):

    print("build graph")
    results = []

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")
    
    if flavor == "mumu":
    
        if not dataset.is_data: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterMC(Muon_pt, Muon_eta, Muon_phi, Muon_charge, Muon_genPartIdx, GenPart_pt, Muon_nTrackerLayers)")
            #df = df.Alias("Muon_pt_corr", "Muon_pt")
            df = df.Filter("HLT_Mu17")
            
        else: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterData(Muon_pt, Muon_eta, Muon_phi, Muon_charge)")
            #df = df.Alias("Muon_pt_corr", "Muon_pt")
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
        if dataset.name == "DYmumu" or dataset.name == "DYee":
            df = df.Define("nominal_pdf_cen", theory_tools.pdf_central_weight(dataset.name, "nnpdf31"))
            df = df.Define("nominal_weight", "weight*SFMC*nominal_pdf_cen")
        else: 
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
    #df = recoilHelper.recoil_setup_Z(df, results, "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi", "Lep_pt", "Lep_phi", "Lep_pt_uncorr")
    df = recoilHelper.recoil_setup_Z(df, results, "MET_pt", "MET_phi", "Lep_pt", "Lep_phi", "Lep_pt_uncorr")
    df = recoilHelper.recoil_apply_Z(df, results, dataset.name, ["DYee", "DYmumu"]) # produces corrected MET as MET_corr_rec_pt/phi
    
   
    reco_mll_cols = ["recoil_corr_magn", "massZ"]
    results.append(df.HistoBoost("reco_mll", reco_mll_axes, [*reco_mll_cols, "nominal_weight"]))
    
    #results.extend(theory_tools.define_and_make_pdf_hists(df, reco_mll_axes, reco_mll_cols, dataset.name, hname="reco_mll"))
    

    
    # W-like
    df = df.Define("TrigMuon_charge", "event % 2 == 0 ? -1 : 1")
    df = df.Define("NonTrigMuon_charge", "-TrigMuon_charge")
    df = df.Define("trigMuons", "Lep_charge == TrigMuon_charge")
    df = df.Define("nonTrigMuons", "Lep_charge == NonTrigMuon_charge")
    df = df.Define("TrigMuon_pt", "Lep_pt[trigMuons][0]")
    df = df.Define("TrigMuon_eta", "Lep_eta[trigMuons][0]")
    df = df.Define("TrigMuon_phi", "Lep_phi[trigMuons][0]")

    df = df.Define("NonTrigMuon_pt", "Lep_pt[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_eta", "Lep_eta[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_phi", "Lep_phi[nonTrigMuons][0]")
    df = df.Define("transverseMass", "wrem::mt_wlike_nano(TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi, MET_corr_rec_pt, MET_corr_rec_phi)")
    #df = df.Filter("transverseMass >= 40.")
    
    results.append(df.HistoBoost("mt", [axis_mt], ["transverseMass", "nominal_weight"]))
    
    
    if dataset.is_data: return results, weightsum
    
    if dataset.name == "DYmumu" or dataset.name == "DYee":
    
        df = wremnants.define_prefsr_vars(df)
    
        gen_reco_mll_cols = ["ptVgen", "recoil_corr_magn", "massZ"]
        results.append(df.HistoBoost("gen_reco_mll", gen_reco_mll_axes, [*gen_reco_mll_cols, "nominal_weight"]))
        

        # pdfs
        results.extend(theory_tools.define_and_make_pdf_hists(df, gen_reco_mll_axes, gen_reco_mll_cols, dataset.name, hname="gen_reco_mll"))
        results.extend(theory_tools.define_and_make_pdf_hists(df, [axis_mt], ["transverseMass"], dataset.name, hname="mt"))

        # QCD scale
        df = theory_tools.define_scale_tensor(df)
        df = df.Define("helicityWeight_tensor", qcdScaleByHelicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
        qcdScaleByHelicityUnc = df.HistoBoost("gen_reco_mll_qcdScaleByHelicity", gen_reco_mll_axes+[axis_ptVgen, axis_chargeVgen], [*gen_reco_mll_cols, "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=qcdScaleByHelicity_helper.tensor_axes)
        results.append(qcdScaleByHelicityUnc)
        qcdScaleByHelicityUnc = df.HistoBoost("mt_qcdScaleByHelicity", [axis_mt]+[axis_ptVgen, axis_chargeVgen], ["transverseMass", "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=qcdScaleByHelicity_helper.tensor_axes)
        results.append(qcdScaleByHelicityUnc)


        # recoil stat uncertainty (for recoil, MET and mT)
        df = recoilHelper.recoil_Z_statUnc_lowPU(df, results, axis_recoil_gen, axis_recoil_reco, axis_mt, axis_mll)

        # lepton efficiencies
        if dataset.name == "DYmumu":
        
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_stat", 120, "wrem::lepSF_HLT_var_mu(1, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll",  gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_syst", 120, "wrem::lepSF_HLT_var_mu(2, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll",  gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_stat",   120, "wrem::lepSF_HLT_var_mu(-1, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_syst",   120, "wrem::lepSF_HLT_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_stat",      36,  "wrem::lepSF_ISO_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_DATA_syst", 36,  "wrem::lepSF_ISO_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_MC_syst",   36,  "wrem::lepSF_ISO_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_stat",     36,  "wrem::lepSF_IDIP_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_DATA_syst",36,  "wrem::lepSF_IDIP_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_MC_syst",  36,  "wrem::lepSF_IDIP_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "gen_reco_mll", gen_reco_mll_axes, gen_reco_mll_cols)
   
        
        # prefire
        df = df.Define("prefireCorr_syst", "wrem::prefireCorr_syst(Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Lep_pt, Lep_eta, Lep_phi)")
        df = df.Define("prefireCorr_syst_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; auto w = nominal_weight*prefireCorr_syst; std::copy(std::begin(w), std::end(w), res.data()); return res;")
        results.append(df.HistoBoost("gen_reco_mll_prefireCorr", [*gen_reco_mll_axes], [*gen_reco_mll_cols, "prefireCorr_syst_tensor"], tensor_axes = [down_up_axis]))

        
        # Breit-Wigner mass weights
        nweights = 21
        df = df.Define("MEParamWeight", "wrem::breitWignerWeights(massVgen, 0)")
        df = df.Define("massWeight_tensor", f"auto res = wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight); res = nominal_weight*res; return res;")

        results.append(df.HistoBoost("gen_reco_mll_massWeight", gen_reco_mll_axes, [*gen_reco_mll_cols, "massWeight_tensor"]))
        results.append(df.HistoBoost("mt_massWeight", [axis_mt], ["transverseMass", "massWeight_tensor"]))
        
    else:
    
        # prefire
        df = df.Define("prefireCorr_syst", "wrem::prefireCorr_syst(Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Lep_pt, Lep_eta, Lep_phi)")
        df = df.Define("prefireCorr_syst_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; auto w = nominal_weight*prefireCorr_syst; std::copy(std::begin(w), std::end(w), res.data()); return res;")
        results.append(df.HistoBoost("reco_mll_prefireCorr", reco_mll_axes, [*reco_mll_cols, "prefireCorr_syst_tensor"], tensor_axes = [down_up_axis]))

        # lepton efficiencies
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_stat", 120, "wrem::lepSF_HLT_var_mu(1, Lep_pt, Lep_eta, Lep_charge)", "reco_mll",  reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_syst", 120, "wrem::lepSF_HLT_var_mu(2, Lep_pt, Lep_eta, Lep_charge)", "reco_mll",  reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_stat",   120, "wrem::lepSF_HLT_var_mu(-1, Lep_pt, Lep_eta, Lep_charge)", "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_syst",   120, "wrem::lepSF_HLT_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_stat",      36,  "wrem::lepSF_ISO_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_DATA_syst", 36,  "wrem::lepSF_ISO_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_MC_syst",   36,  "wrem::lepSF_ISO_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_stat",     36,  "wrem::lepSF_IDIP_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_DATA_syst",36,  "wrem::lepSF_IDIP_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mll", reco_mll_axes, reco_mll_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_MC_syst",  36,  "wrem::lepSF_IDIP_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mll", reco_mll_axes, reco_mll_cols)
   
                
        
    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, fname, args.postfix)
