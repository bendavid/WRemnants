import argparse
from utilities import output_tools
from utilities import common, logging

parser,initargs = common.common_parser()
parser.add_argument("--flavor", type=str, choices=["ee", "mumu"], help="Flavor (ee or mumu)", default="mumu")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist
import scripts.lowPU.config as lowPUcfg


corr_helpers = theory_corrections.load_corr_helpers(common.zprocs_lowpu, args.theoryCorr)

###################################
flavor = args.flavor # mumu, ee
sigProcs = ["Zmumu"] if flavor == "mumu" else ["Zee"]

datasets = wremnants.datasetsLowPU.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              flavor=flavor)


for d in datasets: logger.info(f"Dataset {d.name}")


# load lowPU specific libs
#ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")
narf.clingutils.Declare('#include "lowpu_utils.h"')
narf.clingutils.Declare('#include "lowpu_efficiencies.h"')
narf.clingutils.Declare('#include "lowpu_prefire.h"')
narf.clingutils.Declare('#include "lowpu_rochester.h"')
narf.clingutils.Declare('#include "lowpu_recoil.h"')


# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
#axis_mll = hist.axis.Regular(60, 60., 120., underflow=False, overflow=False, name = "mll")
axis_yll = hist.axis.Regular(50, -2.5, 2.5, name = "yll")
axis_ptll = hist.axis.Regular(300, 0, 300,  name = "ptll")

axis_ptl = hist.axis.Regular(100, 0., 200., name = "ptl")
axis_etal = hist.axis.Regular(50, -2.5, 2.5, name = "etal")


bins_mll = [60, 65, 70, 72, 74, 76, 78] + list(range(80, 100, 1)) + [100, 102, 104, 106, 108, 110, 115, 120]
axis_mll = hist.axis.Variable(bins_mll, name = "mll") 
axis_lin = hist.axis.Regular(5, 0, 5, name = "lin")



qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
axis_ptVgen = qcdScaleByHelicity_helper.hist.axes["ptVgen"]
axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]




# axes for final cards/fitting
reco_mll_axes = [common.axis_recoil_reco_ptZ_lowpu, common.axis_mt_lowpu, axis_charge]
gen_reco_mll_axes = [common.axis_recoil_gen_ptZ_lowpu, common.axis_recoil_reco_ptZ_lowpu, common.axis_mt_lowpu, axis_charge]
axis_mt = hist.axis.Regular(200, 0., 200., name = "mt", underflow=False)
axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)



# recoil initialization
from wremnants import recoil_tools
recoilHelper = recoil_tools.Recoil("lowPU", args, flavor)



def build_graph(df, dataset):

    print("build graph")
    results = []

    isW = dataset.name in common.wprocs_lowpu
    isZ = dataset.name in common.zprocs_lowpu

    if dataset.is_data: df = df.DefinePerSample("weight", "1.0")
    else: df = df.Define("weight", "std::copysign(1.0, genWeight)")
  
    weightsum = df.SumAndCount("weight")
    df = df.Define("TrigMuon_charge", "event % 2 == 0 ? -1 : 1") # wlike charge
    
    # normalization xsecs (propagate pdfs/qcdscales)
    apply_theory_corr = args.theoryCorr and dataset.name in corr_helpers

    if dataset.name in sigProcs:

        axes_xnorm = [axis_xnorm, common.axis_recoil_gen_ptZ_lowpu, axis_charge]
        cols_xnorm = ["xnorm", "ptVgen", "TrigMuon_charge"]
        
        df_xnorm = df
        df_xnorm = df_xnorm.DefinePerSample("exp_weight", "1.0")
        df_xnorm = theory_tools.define_theory_weights_and_corrs(df_xnorm, dataset.name, corr_helpers, args)
        df_xnorm = df_xnorm.Define("xnorm", "0.5")
        results.append(df_xnorm.HistoBoost("xnorm", axes_xnorm, [*cols_xnorm, "nominal_weight"]))

        df_xnorm = theory_tools.define_scale_tensor(df_xnorm)        
        syst_tools.add_pdf_hists(results, df_xnorm, dataset.name, axes_xnorm, cols_xnorm, args.pdfs, "xnorm")
        syst_tools.add_qcdScale_hist(results, df_xnorm, [*axes_xnorm, axis_ptVgen, axis_chargeVgen], [*cols_xnorm, "ptVgen", "chargeVgen"], "xnorm")
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df_xnorm, qcdScaleByHelicity_helper, [*axes_xnorm, axis_ptVgen, axis_chargeVgen], [*cols_xnorm, "ptVgen", "chargeVgen"], base_name="xnorm")
     
        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df_xnorm, "xnorm", axes=axes_xnorm, cols=cols_xnorm, 
                helpers=corr_helpers[dataset.name], generators=args.theoryCorr, modify_central_weight=not args.theoryCorrAltOnly, isW=isW)
            )
  
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
    
        # undo the scale/smearing corrections, needed to correct RawMET
        df = df.Define("Lep_pt_uncorr", "wrem::Egamma_undoCorrection(Electron_pt, Electron_eta, Electron_ecalCorr)")

        if not dataset.is_data: 
            
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(0, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 1)")
            df = df.Filter("HLT_Ele20_WPLoose_Gsf")
            
        else: 
            df = df.Define("Electron_pt_corr", "wrem::applyEGammaScaleSmearingUnc(1, Electron_pt, Electron_eta, Electron_dEscaleUp, Electron_dEscaleDown, Electron_dEsigmaUp, Electron_dEsigmaDown, 1)")
            df = df.Filter("HLT_HIEle20_WPLoose_Gsf")
            
            
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
        
        #df = df.Define("Lep_pt_uncorr", "Electron_pt[goodLeptons]")
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
        df = df.Define("exp_weight", "SFMC")
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)
    else:
        df = df.DefinePerSample("nominal_weight", "1.0")

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
    
    
    

    
    # W-like
    #df = df.Define("TrigMuon_charge", "event % 2 == 0 ? -1 : 1")
    df = df.Define("NonTrigMuon_charge", "-TrigMuon_charge")
    df = df.Define("trigMuons", "Lep_charge == TrigMuon_charge")
    df = df.Define("nonTrigMuons", "Lep_charge == NonTrigMuon_charge")
    df = df.Define("TrigMuon_pt", "Lep_pt[trigMuons][0]")
    df = df.Define("TrigMuon_eta", "Lep_eta[trigMuons][0]")
    df = df.Define("TrigMuon_phi", "Lep_phi[trigMuons][0]")

    df = df.Define("NonTrigMuon_pt", "Lep_pt[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_eta", "Lep_eta[nonTrigMuons][0]")
    df = df.Define("NonTrigMuon_phi", "Lep_phi[nonTrigMuons][0]")

    lep_cols = ["Lep_pt", "Lep_phi", "Lep_pt_uncorr"]
    trg_cols = ["TrigMuon_pt", "TrigMuon_phi", "NonTrigMuon_pt", "NonTrigMuon_phi"]
    df = recoilHelper.recoil_Z(df, results, dataset, common.zprocs_recoil_lowpu, lep_cols, trg_cols) # produces corrected MET as MET_corr_rec_pt/phi
    
    results.append(df.HistoBoost("mZ", [axis_mll], ["massZ", "nominal_weight"]))
    results.append(df.HistoBoost("yZ", [axis_yll], ["yZ", "nominal_weight"]))
    results.append(df.HistoBoost("ptZ", [axis_ptll], ["ptZ", "nominal_weight"]))
    
    df = df.Define("mT_wlike", f"wrem::mt_wlike_nano(TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi, MET_corr_rec_pt, MET_corr_rec_phi)")
    results.append(df.HistoBoost("mT_wlike", [axis_mt], ["mT_wlike", "nominal_weight"]))
    
    results.append(df.HistoBoost("lep_pT", [axis_pt], ["TrigMuon_pt", "nominal_weight"]))
    results.append(df.HistoBoost("lep_pT_qTrw", [axis_pt], ["TrigMuon_pt", "nominal_weight_qTrw"]))
    
    gen_reco_mt_cols = ["ptVgen", "recoil_corr_rec_magn", "mT_wlike", "TrigMuon_charge"]
    reco_mt_cols = ["recoil_corr_rec_magn", "mT_wlike", "TrigMuon_charge"]

    
    if dataset.name in common.zprocs_lowpu:
    
        if dataset.name in sigProcs:
            syst_tools.add_pdf_hists(results, df, dataset.name, gen_reco_mll_axes, gen_reco_mt_cols, args.pdfs, "reco_mT")
        else:
            syst_tools.add_pdf_hists(results, df, dataset.name, reco_mll_axes, reco_mt_cols, args.pdfs, "reco_mT")
        syst_tools.add_pdf_hists(results, df, dataset.name, [axis_mt], ["mT_corr_rec"], args.pdfs, "mt")

        # QCD scale
        df = theory_tools.define_scale_tensor(df)
        syst_tools.add_qcdScale_hist(results, df, [*gen_reco_mll_axes, axis_ptVgen, axis_chargeVgen], [*gen_reco_mt_cols, "ptVgen", "chargeVgen"], "reco_mT") 
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, [*gen_reco_mll_axes, axis_ptVgen, axis_chargeVgen], [*gen_reco_mt_cols, "ptVgen", "chargeVgen"], base_name="reco_mT")
        syst_tools.add_qcdScaleByHelicityUnc_hist(results, df, qcdScaleByHelicity_helper, [axis_mt, axis_ptVgen, axis_chargeVgen], ["mT_corr_rec", "ptVgen", "chargeVgen"], base_name="mt")
    
    # TODO: Should this also be added for the mT hist?

    if apply_theory_corr:
        results.extend(theory_tools.make_theory_corr_hists(df, "reco_mT", axes=gen_reco_mll_axes, cols=gen_reco_mt_cols, 
            helpers=corr_helpers[dataset.name], generators=args.theoryCorr, modify_central_weight=not args.theoryCorrAltOnly, isW=isW)
        )

    if dataset.name in sigProcs:
    
        results.append(df.HistoBoost("reco_mT", gen_reco_mll_axes, [*gen_reco_mt_cols, "nominal_weight"]))
        df = recoilHelper.add_recoil_unc_Z(df, results, dataset, gen_reco_mt_cols, gen_reco_mll_axes, "reco_mT")

        # lepton efficiencies
        if dataset.name == "Zmumu":
        
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_stat", 120, "wrem::lepSF_HLT_var_mu(1, Lep_pt, Lep_eta, Lep_charge)", "reco_mT",  gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_syst", 120, "wrem::lepSF_HLT_var_mu(2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT",  gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_stat",   120, "wrem::lepSF_HLT_var_mu(-1, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_syst",   120, "wrem::lepSF_HLT_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_stat",      36,  "wrem::lepSF_ISO_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_DATA_syst", 36,  "wrem::lepSF_ISO_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_MC_syst",   36,  "wrem::lepSF_ISO_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_stat",     36,  "wrem::lepSF_IDIP_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_DATA_syst",36,  "wrem::lepSF_IDIP_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
            df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_MC_syst",  36,  "wrem::lepSF_IDIP_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", gen_reco_mll_axes, gen_reco_mt_cols)
   
        
        # prefire
        df = df.Define("prefireCorr_syst", "wrem::prefireCorr_syst(Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Lep_pt, Lep_eta, Lep_phi)")
        df = df.Define("prefireCorr_syst_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; auto w = nominal_weight*prefireCorr_syst; std::copy(std::begin(w), std::end(w), res.data()); return res;")
        results.append(df.HistoBoost("reco_mT_prefireCorr", [*gen_reco_mll_axes], [*gen_reco_mt_cols, "prefireCorr_syst_tensor"], tensor_axes = [common.down_up_axis]))

        
        # mass weights (Breit-Wigner and nominal)
        df = syst_tools.define_mass_weights(df, dataset.name)
        syst_tools.add_massweights_hist(results, df, [common.axis_mll_lowpu], ["massZ"], base_name="mll_massWeight", proc=dataset.name)
        syst_tools.add_massweights_hist(results, df, gen_reco_mll_axes, [*gen_reco_mt_cols], base_name="reco_mll_massWeight", proc=dataset.name)
        syst_tools.add_massweights_hist(results, df, [axis_mt], ["mT_corr_rec"], base_name="mT_corr_rec_massWeight", proc=dataset.name)

        
    else:
        
        
        results.append(df.HistoBoost("reco_mT", reco_mll_axes, [*reco_mt_cols, "nominal_weight"]))
        if dataset.is_data: return results, weightsum
    
        # prefire
        df = df.Define("prefireCorr_syst", "wrem::prefireCorr_syst(Jet_pt, Jet_eta, Jet_phi, Jet_muEF, Jet_neEmEF, Jet_chEmEF, Photon_pt, Photon_eta, Photon_phi, Lep_pt, Lep_eta, Lep_phi)")
        df = df.Define("prefireCorr_syst_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; auto w = nominal_weight*prefireCorr_syst; std::copy(std::begin(w), std::end(w), res.data()); return res;")
        results.append(df.HistoBoost("reco_mT_prefireCorr", reco_mll_axes, [*reco_mt_cols, "prefireCorr_syst_tensor"], tensor_axes = [common.down_up_axis]))

        # lepton efficiencies
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_stat", 120, "wrem::lepSF_HLT_var_mu(1, Lep_pt, Lep_eta, Lep_charge)", "reco_mT",  reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_DATA_syst", 120, "wrem::lepSF_HLT_var_mu(2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT",  reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_stat",   120, "wrem::lepSF_HLT_var_mu(-1, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_HLT_MC_syst",   120, "wrem::lepSF_HLT_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_stat",      36,  "wrem::lepSF_ISO_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_DATA_syst", 36,  "wrem::lepSF_ISO_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_ISO_MC_syst",   36,  "wrem::lepSF_ISO_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_stat",     36,  "wrem::lepSF_IDIP_var_mu(1, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_DATA_syst",36,  "wrem::lepSF_IDIP_var_mu(2, Lep_pt, Lep_eta, Lep_charge)",  "reco_mT", reco_mll_axes, reco_mt_cols)
        df = lowPUcfg.lepSF_systs(df, results, "lepSF_IDIP_MC_syst",  36,  "wrem::lepSF_IDIP_var_mu(-2, Lep_pt, Lep_eta, Lep_charge)", "reco_mT", reco_mll_axes, reco_mt_cols)
   
                
        
    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
fname = "lowPU_%s.hdf5" % flavor
output_tools.write_analysis_output(resultdict, fname, args, update_name=not args.forceDefaultName)
