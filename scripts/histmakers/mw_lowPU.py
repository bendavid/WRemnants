import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=None)
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None)
parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")
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

import scripts.lowPU.config as lowPUcfg


###################################
flavor = args.flavor # mu, e
met = args.met # mumu, ee

if flavor == "mu": 
 
    procs = ["DYmumu", "singlemuon", "DYee", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]
    #procs = ["WplusJetsToMuNu", "WminusJetsToMuNu"]
    
elif flavor == "e":

    procs = ["DYee", "singleelectron", "DYmumu", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]
    #procs = ["DYee"]
    
else: sys.exit("Flavor must be e or mu")

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
#axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_passMT = hist.axis.Boolean(name = "passMT")
axis_mll = hist.axis.Regular(60, 60., 120., underflow=False, overflow=False, name = "mll")
axis_yll = hist.axis.Regular(50, -2.5, 2.5, name = "yll")
axis_ptl = hist.axis.Regular(100, 0., 200., name = "ptl")
axis_etal = hist.axis.Regular(50, -2.5, 2.5, name = "etal")

axis_iso = hist.axis.Regular(100, 0, 5, underflow=False, overflow=True, name = "iso")

axis_passIso = hist.axis.Boolean(name = "passIso")

axis_MET_pt = hist.axis.Regular(300, 0, 300, name = "MET_pt", underflow=False)

axis_recoil_magn = hist.axis.Regular(300, 0, 300, name = "recoil_magn", underflow=False)

#axis_eta = hist.axis.Variable([0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4], name = "eta")
axis_eta = hist.axis.Variable([-2.4, 2.4], name = "eta")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge


nominal_axes = [axis_eta, axis_pt, axis_charge]
axis_lin = hist.axis.Regular(5, 0, 5, name = "lin")


axis_ptll = hist.axis.Variable([0, 2, 3, 4, 4.75, 5.5, 6.5, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100], name = "ptll")
axis_costhetastarll = hist.axis.Regular(20, -1., 1., name = "costhetastarll")
axis_phistarll = hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phistarll")


qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = False)
axis_ptVgen = qcdScaleByHelicity_helper.hist.axes["ptVgen"]
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]


# unfolding axes
axis_recoil_reco = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 150], name = "recoil_reco", underflow=False, overflow=True)
axis_recoil_gen = hist.axis.Variable([0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 150], name = "recoil_gen", underflow=False, overflow=True)


# axes for final cards/fitting
reco_mll_axes = [axis_recoil_reco, axis_mll]
gen_reco_mll_axes = [axis_recoil_gen, axis_recoil_reco, axis_mll]
#axis_mt = hist.axis.Regular(200, 0., 200., name = "mt", underflow=False)
#axis_mt = hist.axis.Variable([0, 10, 15, 20, 25, 30, 35,] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)
axis_mt = hist.axis.Variable([0] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)

#axis_mt = hist.axis.Variable(list(range(40, 100, 2)), name = "mt",underflow=False, overflow=False)

# extra axes which can be used to label tensor_axes
down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
axis_cutFlow = hist.axis.Regular(1, 0, 1, name = "cutFlow")


# recoil initialization
from wremnants import recoil_tools
recoilHelper = recoil_tools.Recoil("lowPU", flavor, met)


def build_graph(df, dataset):

    print("build graph")
    results = []

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")
    weightsum = df.SumAndCount("weight")
    
    if flavor == "mu":
    
        if not dataset.is_data: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterMC(Muon_pt, Muon_eta, Muon_phi, Muon_charge, Muon_genPartIdx, GenPart_pt, Muon_nTrackerLayers)")
            #df = df.Alias("Muon_pt_corr", "Muon_pt")
            df = df.Filter("HLT_Mu17")
            
        else: 
        
            df = df.Define("Muon_pt_corr", "wrem::applyRochesterData(Muon_pt, Muon_eta, Muon_phi, Muon_charge)")
            #df = df.Alias("Muon_pt_corr", "Muon_pt")
            df = df.Filter("HLT_HIMu17")
        
        
        df = df.Define("vetoMuons", "Muon_pt_corr > 10 && Muon_looseId && abs(Muon_eta) < 2.4 && abs(Muon_dxybs) < 0.05")
        df = df.Filter("Sum(vetoMuons) == 1")
        
        df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4")
        df = df.Filter("Sum(vetoElectrons) == 0")
        
        df = df.Define("goodLeptons", "vetoMuons && Muon_pt_corr > 25 && Muon_mediumId") # ISO requirement comes later   && Muon_pfRelIso04_all < 0.15
        df = df.Define("goodLeptonsPlus", "goodLeptons && Muon_charge > 0")
        df = df.Define("goodLeptonsMinus", "goodLeptons && Muon_charge < 0")
        df = df.Filter("Sum(goodLeptons) == 1")
        
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidateLowPU(TrigObj_id, TrigObj_pt, TrigObj_l1pt, TrigObj_l2pt, TrigObj_filterBits)")
        df = df.Define("trigMatch", "wrem::hasTriggerMatchLowPU(Muon_eta[goodLeptons], Muon_phi[goodLeptons], TrigObj_eta[goodTrigObjs], TrigObj_phi[goodTrigObjs])")
        df = df.Define("nonTrigMatch", "wrem::inverse(trigMatch)")
        df = df.Filter("Sum(trigMatch) > 0")

        df = df.Define("Lep_pt_uncorr", "Muon_pt[goodLeptons][0]")
        df = df.Define("Lep_pt", "Muon_pt_corr[goodLeptons][0]")
        df = df.Define("Lep_eta", "Muon_eta[goodLeptons][0]")
        df = df.Define("Lep_abs_eta", "abs(Lep_eta)")
        df = df.Define("Lep_phi", "Muon_phi[goodLeptons][0]")
        df = df.Define("Lep_charge", "Muon_charge[goodLeptons][0]")
        df = df.Define("Lep_mass", "Muon_mass[goodLeptons][0]")
        df = df.Define("Lep_iso", "Muon_pfRelIso04_all[goodLeptons][0]")
        
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
        df = df.Define("Lep_charge", "Electron_charge[goodLeptons][0]")
        df = df.Define("Lep_mass", "Electron_mass[goodLeptons]")
        
        if not dataset.is_data:
            df = df.Define("lepSF_IDISO", "wrem::lepSF(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 3)")
            df = df.Define("lepSF_HLT", "wrem::lepSF_HLT_q(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 11)")
            df = df.Define("prefireCorr", "wrem::prefireCorr(0, Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_phi[goodLeptons])")
            df = df.Define("SFMC", "lepSF_IDISO*lepSF_HLT*prefireCorr")
        
        else: df = df.Define("SFMC", "1.0")    
    
   

    df = df.Filter("Lep_pt > 25")
    
    
        
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")
    
    df = df.Define("Lep1_mom4", "ROOT::Math::PtEtaPhiMVector(Lep_pt, Lep_eta, Lep_phi, Lep_mass)")
    #df = df.Define("Z_mom4", "ROOT::Math::PxPyPzEVector(Lep1_mom4) + ROOT::Math::PxPyPzEVector(Lep2_mom4)")
    #df = df.Define("ptZ", "Z_mom4.pt()")
    #df = df.Define("massZ", "Z_mom4.mass()")
    #df = df.Define("yZ", "Z_mom4.Rapidity()")
    #df = df.Define("absYZ", "std::fabs(yZ)")
    #df = df.Filter("massZ > 60 && massZ < 120")

    if not dataset.is_data: 
        weight_expr = "weight*SFMC"
        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")

    # isolation requirement
    #df = df.Filter("Lep_iso < 0.15")
    df = df.Define("passIso", "Lep_iso < 0.15")
    
    
    #nominal_cols = ["mT_corr_rec", "Lep_abs_eta", "Lep_charge", "passIso"]
    #nominal_axes = [axis_mt, axis_eta, axis_charge, axis_passIso]

    
    results.append(df.HistoBoost("lep_pt", [axis_ptl], ["Lep_pt", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta", [axis_etal], ["Lep_eta", "nominal_weight"]))
    results.append(df.HistoBoost("lep_iso", [axis_iso], ["Lep_iso", "nominal_weight"]))
    
    df = df.Define("noTrigMatch", "Sum(trigMatch)")
    results.append(df.HistoBoost("noTrigMatch", [axis_lin], ["noTrigMatch", "nominal_weight"]))

    if no is_data:

    # Recoil calibrations
    df = recoilHelper.setup_MET(df, results, dataset, "Lep_pt", "Lep_phi", "Lep_pt_uncorr")
    df = recoilHelper.setup_gen(df, results, dataset, ["WplusJetsToMuNu", "WminusJetsToMuNu"])
    df = recoilHelper.apply_recoil_W(df, results, dataset, ["WplusJetsToMuNu", "WminusJetsToMuNu"]) # produces corrected MET as MET_corr_rec_pt/phi


    df = df.Define("mT_corr_rec", "wrem::mt_2(Lep_pt, Lep_phi, MET_corr_rec_pt, MET_corr_rec_phi)")
    df = df.Define("passMT", "mT_corr_rec > 40")
    

    results.append(df.HistoBoost("MET_corr_rec_pt", [axis_MET_pt, axis_charge, axis_passMT, axis_passIso], ["MET_corr_rec_pt", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("mT_corr_rec", [axis_mt, axis_charge, axis_passMT, axis_passIso], ["mT_corr_rec", "Lep_charge", "passMT", "passIso", "nominal_weight"]))      
    results.append(df.HistoBoost("recoil_corr_rec_magn", [axis_recoil_magn, axis_charge, axis_passMT, axis_passIso], ["recoil_corr_rec_magn", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
        
    

    if dataset.is_data: return results, weightsum


    # QCD scales same for Tau and Mu
    if dataset.name == "WplusJetsToMuNu" or dataset.name == "WminusJetsToMuNu" or dataset.name == "WplusJetsToTauNu" or dataset.name == "WminusJetsToTauNu":
    
        # QCD scale
        df = theory_tools.define_scale_tensor(df)
        df = df.Define("helicityWeight_tensor", qcdScaleByHelicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
        qcdScaleByHelicityUnc = df.HistoBoost("mT_corr_rec_qcdScaleByHelicity", [axis_mt, axis_charge, axis_passMT, axis_passIso, axis_ptVgen, axis_chargeVgen], ["mT_corr_rec", "Lep_charge", "passMT", "passIso", "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=qcdScaleByHelicity_helper.tensor_axes)
        results.append(qcdScaleByHelicityUnc)
        
    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers
    if apply_theory_corr:
        results.extend(theory_tools.make_theory_corr_hists(df, "mll_reco", axes=gen_reco_mll_axes, cols=gen_reco_mll_cols, 
            helpers=corr_helpers[dataset.name], generators=args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
        )

    ###return results, weightsum    
    if dataset.name == "WplusJetsToMuNu" or dataset.name == "WminusJetsToMuNu":
    
        # recoil uncertainties
        df = recoilHelper.recoil_W_unc_lowPU(df, results, axis_charge, axis_mt, axis_recoil_magn, axis_eta, axis_passMT, axis_passIso)
        
        # pdfs
        results.extend(theory_tools.define_and_make_pdf_hists(df, [axis_mt, axis_charge, axis_passMT, axis_passIso], ["mT_corr_rec", "Lep_charge", "passMT", "passIso"], dataset.name, hname="mT_corr_rec"))        


        '''

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
        '''
        
        # Breit-Wigner mass weights
        nweights = 21
        df = df.Define("MEParamWeight", "wrem::breitWignerWeights(massVgen, 1)")
        df = df.Define("massWeight_tensor", f"auto res = wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight); res = nominal_weight*res; return res;")
        df = df.Define("massWeight_tensor_unscaled", f"auto res = wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight); res = res; return res;")
        #results.append(df.HistoBoost("gen_reco_mll_massWeight", gen_reco_mll_axes, [*gen_reco_mll_cols, "massWeight_tensor"]))
        #results.append(df.HistoBoost("mT_uncorr_massWeight", [axis_mt, axis_eta, axis_charge], ["mT_uncorr", "Lep_abs_eta", "Lep_charge", "massWeight_tensor"]))
        ###results.append(df.HistoBoost("mT_corr_xy_massWeight", [axis_mt, axis_eta, axis_charge, axis_passIso], ["mT_corr_xy", "Lep_abs_eta", "Lep_charge", "passIso", "massWeight_tensor"]))
        results.append(df.HistoBoost("mT_corr_rec_massWeight", [axis_mt, axis_charge, axis_passMT, axis_passIso], ["mT_corr_rec", "Lep_charge", "passMT", "passIso", "massWeight_tensor"]))



        # Muon momentum scale
        netabins = 1
        nweights = 21
        mag = 1.e-4
        df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor_unscaled, Lep_abs_eta, {mag})")
        scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
        dummyMuonScaleSyst = df.HistoBoost("mT_corr_rec_muonScaleSyst", [axis_mt, axis_charge, axis_passMT, axis_passIso], ["mT_corr_rec", "Lep_charge", "passMT", "passIso", f"muonScaleDummy{netabins}Bins"], tensor_axes=[down_up_axis, scale_etabins_axis])
        results.append(dummyMuonScaleSyst)



        '''
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
   
    '''            
        
    return results, weightsum


def build_graph_cutFlow(df, dataset):

    print("build graph")
    results = []

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")
    df = df.Define("cutFlow", "1.0")

    weightsum = df.SumAndCount("weight")
    results.append(df.HistoBoost("cutflow_0", [axis_cutFlow], ["cutFlow", "weight"]))
    
    if flavor == "mu":
    
        # at least 1 muon in the collection
        df = df.Filter("Sum(Muon_pt) > 0")
        results.append(df.HistoBoost("cutflow_1", [axis_cutFlow], ["cutFlow", "weight"]))
    
        # loose muon
        df = df.Define("vetoMuons_1", "Muon_pt > 10 && Muon_looseId && abs(Muon_dxybs) < 0.05")
        df = df.Filter("Sum(vetoMuons_1) >= 0")
        results.append(df.HistoBoost("cutflow_2", [axis_cutFlow], ["cutFlow", "weight"]))
    
        # eta acceptance
        df = df.Define("vetoMuons_2", "vetoMuons_1 && abs(Muon_eta) < 2.4")
        df = df.Filter("Sum(vetoMuons_2) > 0")
        results.append(df.HistoBoost("cutflow_3", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # ID
        df = df.Define("vetoMuons_3", "vetoMuons_2 && Muon_mediumId")
        df = df.Filter("Sum(vetoMuons_3) > 0")
        results.append(df.HistoBoost("cutflow_4", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # ISO
        df = df.Define("vetoMuons_4", "vetoMuons_3 && Muon_pfRelIso04_all < 0.15")
        df = df.Filter("Sum(vetoMuons_4) > 0")
        results.append(df.HistoBoost("cutflow_5", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # pT
        df = df.Define("goodLeptons", "vetoMuons_4 && Muon_pt > 25")
        df = df.Filter("Sum(goodLeptons) > 0")
        results.append(df.HistoBoost("cutflow_6", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # 1 muon
        df = df.Filter("Sum(goodLeptons) == 1")
        results.append(df.HistoBoost("cutflow_7", [axis_cutFlow], ["cutFlow", "weight"]))
        
        
        # electron veto
        df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4")
        df = df.Filter("Sum(vetoElectrons) == 0")
        results.append(df.HistoBoost("cutflow_8", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # trigger
        if not dataset.is_data:  df = df.Filter("HLT_Mu17") 
        else: df = df.Filter("HLT_HIMu17")      
        results.append(df.HistoBoost("cutflow_9", [axis_cutFlow], ["cutFlow", "weight"]))
        
        # trigger match
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidateLowPU(TrigObj_id, TrigObj_pt, TrigObj_l1pt, TrigObj_l2pt, TrigObj_filterBits)")
        df = df.Define("trigMatch", "wrem::hasTriggerMatchLowPU(Muon_eta[goodLeptons], Muon_phi[goodLeptons], TrigObj_eta[goodTrigObjs], TrigObj_phi[goodTrigObjs])")
        df = df.Define("nonTrigMatch", "wrem::inverse(trigMatch)")
        df = df.Filter("Sum(trigMatch) > 0")
        results.append(df.HistoBoost("cutflow_10", [axis_cutFlow], ["cutFlow", "weight"]))
        

        df = df.Define("Lep_pt", "Muon_pt[goodLeptons][0]")
        df = df.Define("Lep_eta", "Muon_eta[goodLeptons][0]")
        df = df.Define("Lep_phi", "Muon_phi[goodLeptons][0]")
        df = df.Define("Lep_charge", "Muon_charge[goodLeptons][0]")
        df = df.Define("Lep_mass", "Muon_mass[goodLeptons][0]")
        df = df.Define("Lep_iso", "Muon_pfRelIso04_all[goodLeptons][0]")
        
          
          
        # MET filters       
        df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")
        results.append(df.HistoBoost("cutflow_11", [axis_cutFlow], ["cutFlow", "weight"]))
    

        # mT cut
        if met == "DeepMETReso": df = df.Define("mT_uncorr", "wrem::mt_2(Lep_pt, Lep_phi, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
        if met == "RawPFMET": df = df.Define("mT_uncorr", "wrem::mt_2(Lep_pt, Lep_phi, RawMET_pt, RawMET_phi)")
        df = df.Filter("mT_uncorr > 40")
        results.append(df.HistoBoost("cutflow_12", [axis_cutFlow], ["cutFlow", "weight"]))


    else:

        if not dataset.is_data: 
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(0, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 0)")
            df = df.Filter("HLT_Ele20_WPLoose_Gsf")
            
        else: 
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(1, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 0)")
            df = df.Filter("HLT_HIEle20_WPLoose_Gsf")
            
        results.append(df.HistoBoost("cutflow_1", [axis_cutFlow], ["cutFlow", "weight"]))    
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
        df = df.Define("Lep_charge", "Electron_charge[goodLeptons][0]")
        df = df.Define("Lep_mass", "Electron_mass[goodLeptons]")
        
        if not dataset.is_data:
            df = df.Define("lepSF_IDISO", "wrem::lepSF(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 3)")
            df = df.Define("lepSF_HLT", "wrem::lepSF_HLT_q(Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_charge[goodLeptons], 11)")
            df = df.Define("prefireCorr", "wrem::prefireCorr(0, Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Electron_pt_corr[goodLeptons], Electron_eta[goodLeptons], Electron_phi[goodLeptons])")
            df = df.Define("SFMC", "lepSF_IDISO*lepSF_HLT*prefireCorr")
        
        else: df = df.Define("SFMC", "1.0")    
    
   


    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph)
fname = "lowPU_%s_%s.pkl.lz4" % (flavor, met)
output_tools.write_analysis_output(resultdict, fname, args)
