import argparse
from utilities import output_tools
from utilities import common, logging

parser,initargs = common.common_parser()
parser.add_argument("--flavor", type=str, choices=["e", "mu"], help="Flavor (e or mu)", default="mu")
args = parser.parse_args()


import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import math
import hist
import ROOT
import scripts.lowPU.config as lowPUcfg


###################################
flavor = args.flavor # mu, e
sigProcs = ["WminusJetsToMuNu", "WplusJetsToMuNu"] if flavor == "mu" else ["WminusJetsToENu", "WplusJetsToENu"]

corr_helpers = theory_corrections.load_corr_helpers(common.wprocs_lowpu, args.theoryCorr)

datasets = wremnants.datasetsLowPU.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              flavor=flavor)

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

for d in datasets: logger.info(f"Dataset {d.name}")

# load lowPU specific libs
#ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")
ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')
ROOT.gInterpreter.Declare('#include "lowpu_efficiencies.h"')
ROOT.gInterpreter.Declare('#include "lowpu_prefire.h"')
ROOT.gInterpreter.Declare('#include "lowpu_rochester.h"')
ROOT.gInterpreter.Declare('#include "lowpu_recoil.h"')




# standard regular axes
axis_eta = hist.axis.Regular(24, -2.4, 2.4, name = "eta", underflow=False, overflow=False)
axis_pt = hist.axis.Regular(75, 25., 100., name = "pt", underflow=False)
axis_phi = hist.axis.Regular(50, -4, 4, name = "phi")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_iso = hist.axis.Regular(50, 0, 1, underflow=False, overflow=True, name = "iso")

axis_passMT = hist.axis.Boolean(name = "passMT")
axis_passIso = hist.axis.Boolean(name = "passIso")
axis_lin = hist.axis.Regular(5, 0, 5, name = "lin")



qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = False)
axis_ptVgen = qcdScaleByHelicity_helper.hist.axes["ptVgen"]
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]



# axes for final cards/fitting
axis_mT = hist.axis.Variable([0] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)
axis_mT = hist.axis.Regular(100, 0, 200, name = "mt", underflow=False)
reco_mT_axes = [common.axis_recoil_reco_ptW_lowpu, common.axis_mt_lowpu, axis_charge, axis_passMT, axis_passIso]
gen_reco_mT_axes = [common.axis_recoil_gen_ptW_lowpu, common.axis_recoil_reco_ptW_lowpu, common.axis_mt_lowpu, axis_charge, axis_passMT, axis_passIso]
axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)


# extra axes which can be used to label tensor_axes
down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
axis_cutFlow = hist.axis.Regular(1, 0, 1, name = "cutFlow")


# recoil initialization
from wremnants import recoil_tools
recoilHelper = recoil_tools.Recoil("lowPU", args, flavor)


def build_graph(df, dataset):

    print("build graph")
    results = []

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")
    weightsum = df.SumAndCount("weight")
    
    
    # normalization xsecs (propagate pdfs/qcdscales)
    if dataset.name in sigProcs:

        axes_xnorm = [axis_xnorm, common.axis_recoil_gen_ptZ_lowpu, axis_charge]
        cols_xnorm = ["xnorm", "ptVgen", "chargeVgen"]
        
        df_xnorm = df
        df_xnorm = df_xnorm.DefinePerSample("exp_weight", "1.0")
        df_xnorm = theory_tools.define_theory_weights_and_corrs(df_xnorm, dataset.name, corr_helpers, args)
        df_xnorm = df_xnorm.Define("xnorm", "0.5")
        results.append(df_xnorm.HistoBoost("xnorm", axes_xnorm, [*cols_xnorm, "nominal_weight"]))

        df_xnorm = theory_tools.define_scale_tensor(df_xnorm)        
        syst_tools.add_pdf_hists(results, df_xnorm, dataset.name, axes_xnorm, cols_xnorm, args.pdfs, "xnorm")
        syst_tools.add_qcdScale_hist(results, df_xnorm, [*axes_xnorm, axis_ptVgen, axis_chargeVgen], [*cols_xnorm, "ptVgen", "chargeVgen"], "xnorm")
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df_xnorm, qcdScaleByHelicity_helper, [*axes_xnorm, axis_ptVgen, axis_chargeVgen], [*cols_xnorm, "ptVgen", "chargeVgen"], base_name="xnorm")
     
        
    
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
        df = df.Filter("Sum(goodLeptons) == 1")
        
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
    df = df.Define("passIso", "Lep_iso < 0.15")


    if not dataset.is_data: 
        df = df.Define("exp_weight", "SFMC")
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")

    df = df.Define("noTrigMatch", "Sum(trigMatch)")
    results.append(df.HistoBoost("noTrigMatch", [axis_lin], ["noTrigMatch", "nominal_weight"]))


    # Recoil calibrations
    lep_cols = ["Lep_pt", "Lep_phi", "Lep_charge", "Lep_pt_uncorr"]
    df = recoilHelper.recoil_W(df, results, dataset, common.vprocs_lowpu, lep_cols) # produces corrected MET as MET_corr_rec_pt/phi  vprocs_lowpu wprocs_recoil_lowpu
   
    df = df.Alias("mT", "mT_corr_rec")
    df = df.Define("passMT", "mT > 40")
    
    results.append(df.HistoBoost("lep_pt", [axis_pt, axis_charge, axis_passMT, axis_passIso], ["Lep_pt", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("lep_eta", [axis_eta, axis_charge, axis_passMT, axis_passIso], ["Lep_eta", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("lep_phi", [axis_phi, axis_charge, axis_passMT, axis_passIso], ["Lep_phi", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pt", [axis_pt, axis_charge, axis_passMT, axis_passIso], ["Lep_pt", "Lep_charge", "passMT", "passIso", "nominal_weight"]))
    results.append(df.HistoBoost("lep_iso", [axis_iso], ["Lep_iso", "nominal_weight"]))
   
    results.append(df.HistoBoost("qcd_space", [axis_pt, axis_eta, axis_iso, axis_charge, axis_mT], ["Lep_pt", "Lep_eta", "Lep_iso", "Lep_charge", "mT", "nominal_weight"]))

    
    gen_reco_mT_cols = ["ptVgen", "recoil_corr_rec_magn", "mT", "Lep_charge", "passMT", "passIso"]
    reco_mT_cols = ["recoil_corr_rec_magn", "mT", "Lep_charge", "passMT", "passIso"]
    
    
    
    
    if dataset.name in common.vprocs_lowpu:
    
        # pdfs
        if dataset.name in sigProcs:
            syst_tools.add_pdf_hists(results, df, dataset.name, gen_reco_mT_axes, gen_reco_mT_cols, args.pdfs, "reco_mT")
        else:
            syst_tools.add_pdf_hists(results, df, dataset.name, reco_mT_axes, reco_mT_cols, args.pdfs, "reco_mT")
        syst_tools.add_pdf_hists(results, df, dataset.name, [axis_mT], ["mT"], args.pdfs, "mT")

        # QCD scale
        df = theory_tools.define_scale_tensor(df)
        syst_tools.add_qcdScale_hist(results, df, [*gen_reco_mT_axes, axis_ptVgen, axis_chargeVgen], [*gen_reco_mT_cols, "ptVgen", "chargeVgen"], "reco_mT") 
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, [*gen_reco_mT_axes, axis_ptVgen, axis_chargeVgen], [*gen_reco_mT_cols, "ptVgen", "chargeVgen"], base_name="reco_mT")
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, [axis_mT, axis_ptVgen, axis_chargeVgen], ["mT", "ptVgen", "chargeVgen"], base_name="mT")


    if dataset.name in sigProcs:
    
        results.append(df.HistoBoost("reco_mT", gen_reco_mT_axes, [*gen_reco_mT_cols, "nominal_weight"]))
        df = recoilHelper.add_recoil_unc_W(df, results, dataset, gen_reco_mT_cols, gen_reco_mT_axes, "reco_mT")

        
        # mass weights (Breit-Wigner and nominal)
        #df = syst_tools.define_mass_weights(df, dataset.name)
        #syst_tools.add_massweights_hist(results, df, gen_reco_mT_axes, [*gen_reco_mT_cols], base_name="reco_mT_massWeight", proc=dataset.name)
        #syst_tools.add_massweights_hist(results, df, [axis_mT], ["mT"], base_name="mT_massWeight", proc=dataset.name)


        # Muon momentum scale
        '''
        netabins = 1
        nweights = 21
        mag = 1.e-4
        df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor_unscaled, Lep_abs_eta, {mag})")
        scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
        dummyMuonScaleSyst = df.HistoBoost("mT_corr_rec_muonScaleSyst", [axis_mT, axis_charge, axis_passMT, axis_passIso], ["mT", "Lep_charge", "passMT", "passIso", f"muonScaleDummy{netabins}Bins"], tensor_axes=[down_up_axis, scale_etabins_axis])
        results.append(dummyMuonScaleSyst)
        '''



    else:
    
        results.append(df.HistoBoost("reco_mT", reco_mT_axes, [*reco_mT_cols, "nominal_weight"]))
        if dataset.is_data: return results, weightsum
        
        df = recoilHelper.add_recoil_unc_W(df, results, dataset, reco_mT_cols, reco_mT_axes, "reco_mT")

    
        
     
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
fname = "lowPU_%s.hdf5" % flavor
output_tools.write_analysis_output(resultdict, fname, args)
