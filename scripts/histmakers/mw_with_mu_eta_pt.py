import argparse
from utilities import output_tools,common

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist
import lz4.frame
import logging
import math
import time
import pathlib

data_dir = f"{pathlib.Path(__file__).parent}/../../wremnants/data/"

logging.basicConfig(level=logging.INFO)

parser.add_argument("-e", "--era", type=str, choices=["2016PreVFP","2016PostVFP"], help="Data set to process", default="2016PostVFP")
parser.add_argument("--noMuonCorr", action="store_true", help="Don't use corrected pt-eta-phi-charge")
parser.add_argument("--noScaleFactors", action="store_true", help="Don't use scale factors for efficiency")
parser.add_argument("--muonCorrMag", default=1.e-4, type=float, help="Magnitude of dummy muon momentum calibration uncertainty")
parser.add_argument("--muonCorrEtaBins", default=1, type=int, help="Number of eta bins for dummy muon momentum calibration uncertainty")
parser.add_argument("--lumiUncertainty", type=float, help="Uncertainty for luminosity in excess to 1 (e.g. 1.012 means 1.2\%)", default=1.012)
args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts]) 
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9")

if not args.no_recoil:
    logging.warning("Recoil correction for high PU is not yet supported! Setting false")
    args.no_recoil = True

era = args.era
noMuonCorr = args.noMuonCorr

muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)
#qcdScaleByHelicity_Zhelper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)
qcdScaleByHelicity_Whelper = wremnants.makeQCDScaleByHelicityHelper()

# custom template binning
template_neta = int(args.eta[0])
template_mineta = args.eta[1]
template_maxeta = args.eta[2]
print(f"Eta binning: {template_neta} bins from {template_mineta} to {template_maxeta}")
template_npt = int(args.pt[0])
template_minpt = args.pt[1]
template_maxpt = args.pt[2]
print(f"Pt binning: {template_npt} bins from {template_minpt} to {template_maxpt}")

# standard regular axes
axis_eta = hist.axis.Regular(template_neta, template_mineta, template_maxeta, name = "eta")
axis_pt = hist.axis.Regular(template_npt, template_minpt, template_maxpt, name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

# TODO: get from common
axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")
nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

axis_chargeVgen = qcdScaleByHelicity_Whelper.hist.axes["chargeVgen"]
axis_ptVgen = hist.axis.Variable(
    common.ptV_10quantiles_binning, 
    name = "ptVgen", underflow=False,
)

# axes for mT measurement
axis_mt = hist.axis.Variable([0] + list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200], name = "mt",underflow=False, overflow=True)
axis_eta_mT = hist.axis.Variable([-2.4, 2.4], name = "eta")

# axes for study of fakes
axis_mt_fakes = hist.axis.Regular(60, 0., 120., name = "mt", underflow=False, overflow=True)
axis_njet_fakes = hist.axis.Regular(2, -0.5, 1.5, name = "Numbr of jets", underflow=False, overflow=False) # only need case with 0 jets or > 0

# TODO: add a way to select binned or smoothed SF from command option (need to call different helpers)
muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile, era = era, max_pt = axis_pt.edges[-1])
print(args.sfFile)

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

calibration_helper, calibration_uncertainty_helper = wremnants.make_muon_calibration_helpers()

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

# recoil initialization
if not args.no_recoil:
    from wremnants import recoil_tools
    import ROOT
    ROOT.gInterpreter.Declare('#include "lowpu_recoil.h"')
    ROOT.gInterpreter.Declare('#include "lowpu_utils.h"')
    recoilHelper = recoil_tools.Recoil("highPU", flavor="mu", met="RawPFMET") # DeepMETReso RawPFMET


def build_graph(df, dataset):
    print("build graph", dataset.name)
    results = []

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")
    df = df.Filter("event % 2 == 1") # test with odd/even events

    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"
    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers
    if noMuonCorr:
        df = df.Alias("Muon_correctedPt", "Muon_pt")
        df = df.Alias("Muon_correctedEta", "Muon_eta")
        df = df.Alias("Muon_correctedPhi", "Muon_phi")
        df = df.Alias("Muon_correctedCharge", "Muon_charge")
    else:
        if dataset.is_data:
            #TODO corrections not available for data yet
            df = df.Alias("Muon_correctedPt", "Muon_cvhbsPt")
            df = df.Alias("Muon_correctedEta", "Muon_cvhbsEta")
            df = df.Alias("Muon_correctedPhi", "Muon_cvhbsPhi")
            df = df.Alias("Muon_correctedCharge", "Muon_cvhbsCharge")
        elif isW or isZ:
            df = wremnants.define_corrected_muons(df, calibration_helper)
        else:
            # no track refit available for background monte carlo samples and this is "good enough"
            df = df.Alias("Muon_correctedPt", "Muon_pt")
            df = df.Alias("Muon_correctedEta", "Muon_eta")
            df = df.Alias("Muon_correctedPhi", "Muon_phi")
            df = df.Alias("Muon_correctedCharge", "Muon_charge")
                    
    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")
    df = df.Filter("Sum(vetoMuons) == 1")
    if args.trackerMuons:
        if dataset.group in ["Top", "Diboson"]:
            df = df.Define("Muon_category", "Muon_isTracker")
        else:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_innerTrackOriginalAlgo != 13 && Muon_innerTrackOriginalAlgo != 14")
    else:
        df = df.Define("Muon_category", "Muon_isGlobal")
    df = df.Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_category")
    df = df.Filter("Sum(goodMuons) == 1")

    df = df.Define("goodMuons_pt0", "Muon_correctedPt[goodMuons][0]")
    df = df.Define("goodMuons_eta0", "Muon_correctedEta[goodMuons][0]")
    df = df.Define("goodMuons_abseta0", "abs(goodMuons_eta0)")
    df = df.Define("goodMuons_phi0", "Muon_correctedPhi[goodMuons][0]")
    df = df.Define("goodMuons_charge0", "Muon_correctedCharge[goodMuons][0]")

    df = df.Define("goodMuons_isStandalone", "Muon_isStandalone[goodMuons][0]")

    #standalone quantities, currently only in data and W/Z samples
    if dataset.group in ["Top", "Diboson"]:
        df = df.Alias("goodMuons_SApt0",  "goodMuons_pt0")
        df = df.Alias("goodMuons_SAeta0", "goodMuons_eta0")
        df = df.Alias("goodMuons_SAphi0", "goodMuons_phi0")
    elif args.trackerMuons:
        # try to use standalone variables when possible
        df = df.Define("goodMuons_SApt0",  "goodMuons_isStandalone ? Muon_standalonePt[goodMuons][0] : goodMuons_pt0")
        df = df.Define("goodMuons_SAeta0", "goodMuons_isStandalone ? Muon_standaloneEta[goodMuons][0] : goodMuons_eta0")
        df = df.Define("goodMuons_SAphi0", "goodMuons_isStandalone ? Muon_standalonePhi[goodMuons][0] : goodMuons_phi0")
    else:
        df = df.Define("goodMuons_SApt0",  "Muon_standalonePt[goodMuons][0]")
        df = df.Define("goodMuons_SAeta0", "Muon_standaloneEta[goodMuons][0]")
        df = df.Define("goodMuons_SAphi0", "Muon_standalonePhi[goodMuons][0]")

    # the next cut is mainly needed for consistency with the reco efficiency measurement for the case with global muons
    # note, when SA does not exist this cut is still fine because of how we define these variables
    df = df.Filter("goodMuons_SApt0 > 15.0 && wrem::deltaR2(goodMuons_SAeta0, goodMuons_SAphi0, goodMuons_eta0, goodMuons_phi0) < 0.09")
    
    df = df.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")

    #TODO improve this to include muon mass?
    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")

    df = df.Filter("Sum(vetoElectrons) == 0")

    df = df.Define("goodCleanJets", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && Jet_pt > 30 && abs(Jet_eta) < 2.4 && wrem::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_correctedEta[vetoMuons],Muon_correctedPhi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")
    df = df.Define("goodCleanJetsPt45", "goodCleanJets && Jet_pt > 45")

    df = df.Define("passMT", "transverseMass >= 40.0")
    df = df.Filter("passMT || Sum(goodCleanJets)>=1")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    if dataset.group in ["Top", "Diboson"]:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    else:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Filter("wrem::hasTriggerMatch(goodMuons_eta0,goodMuons_phi0,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")

    # gen match to bare muons to select only prompt muons from top processes
    if isTop:
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && abs(GenPart_pdgId) == 13")
        df = df.Filter("wrem::hasMatchDR2(goodMuons_eta0,goodMuons_phi0,GenPart_eta[postFSRmuons],GenPart_phi[postFSRmuons],0.09)")
    
    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if dataset.is_data:
        df = df.Define("nominal_weight", "1.0")
        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)
        
        if not args.no_recoil:
            df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
            df = df.Define("mT_corr_rec", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_xy_pt, MET_corr_xy_phi)")
            results.append(df.HistoBoost("mT_corr_rec", [axis_mt, axis_eta_mT, axis_charge, axis_passIso], ["mT_corr_rec", "goodMuons_abseta0", "goodMuons_charge0", "passIso"]))

        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, nominal_cols)
        results.append(qcdJetPt45)
        
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId"])

        weight_expr = "weight*weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        if not args.noScaleFactors:
            df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"])
            weight_expr += "*weight_fullMuonSF_withTrackingReco"
        if args.vertex_weight:
            weight_expr += "*weight_vtx"
            
        df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)

        if apply_theory_corr:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal", nominal_axes, nominal_cols, 
                corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
            )

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)
        
        if not args.no_recoil:
            df = recoilHelper.setup_MET(df, results, dataset, "Muon_pt[goodMuons]", "Muon_phi[goodMuons]", "Muon_pt[goodMuons]")
            df = recoilHelper.setup_gen(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"])
            df = recoilHelper.apply_W(df, results, dataset, ["WplusmunuPostVFP", "WminusmunuPostVFP"]) # produces corrected MET as MET_corr_rec_pt/phi

            #df = df.Define("mT_corr_rec", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")
            #results.append(df.HistoBoost("mT_corr_rec", [axis_mt, axis_eta_mT, axis_charge, axis_passIso], ["mT_corr_rec", "goodMuons_abseta0", "goodMuons_charge0", "passIso", "nominal_weight"]))
            #if dataset.name in ["WplusmunuPostVFP", "WminusmunuPostVFP"]: df = recoilHelper.recoil_W_unc_lowPU(df, results, axis_charge, axis_mt, axis_eta, axis_passIso)


    if not dataset.is_data and not args.onlyMainHistograms:
        
        dQCDbkGVar = df.Filter("passMT || Sum(goodCleanJetsPt45)>=1")
        qcdJetPt45 = dQCDbkGVar.HistoBoost("qcdJetPt45", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(qcdJetPt45)

        for key,helper in muon_efficiency_helper_stat.items():
            if "iso" in key:
                df = df.Define(f"effStatTnP_{key}_tensor", helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0", "passIso", "nominal_weight"])
            else:
                df = df.Define(f"effStatTnP_{key}_tensor", helper, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0", "nominal_weight"])
            effStatTnP = df.HistoBoost(f"effStatTnP_{key}", nominal_axes, [*nominal_cols, f"effStatTnP_{key}_tensor"], tensor_axes = helper.tensor_axes)
            results.append(effStatTnP)
        
        df = df.Define("effSystTnP_weight", muon_efficiency_helper_syst, ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso", "nominal_weight"])
        effSystTnP = df.HistoBoost("effSystTnP", nominal_axes, [*nominal_cols, "effSystTnP_weight"], tensor_axes = muon_efficiency_helper_syst.tensor_axes)
        results.append(effSystTnP)

        df = df.Define("muonL1PrefireStat_tensor", muon_prefiring_helper_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId", "nominal_weight"])
        muonL1PrefireStat = df.HistoBoost("muonL1PrefireStat", nominal_axes, [*nominal_cols, "muonL1PrefireStat_tensor"], tensor_axes = muon_prefiring_helper_stat.tensor_axes)
        results.append(muonL1PrefireStat)

        df = df.Define("muonL1PrefireSyst_tensor", muon_prefiring_helper_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_looseId", "nominal_weight"])
        muonL1PrefireSyst = df.HistoBoost("muonL1PrefireSyst", nominal_axes, [*nominal_cols, "muonL1PrefireSyst_tensor"], tensor_axes = [common.down_up_axis])
        results.append(muonL1PrefireSyst)

        df = df.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")
        ecalL1Prefire = df.HistoBoost("ecalL1Prefire", nominal_axes, [*nominal_cols, "ecalL1Prefire_tensor"], tensor_axes = [common.down_up_axis])
        results.append(ecalL1Prefire)
        
        # luminosity, done here as shape variation despite being a flat scaling so to facilitate propagating to fakes afterwards
        df = df.Define("luminosityScaling", f"wrem::constantScaling(nominal_weight, {args.lumiUncertainty})")
        luminosity = df.HistoBoost("luminosity", nominal_axes, [*nominal_cols, "luminosityScaling"], tensor_axes = [common.down_up_axis])
        results.append(luminosity)
                
        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isW or isZ:

            df = theory_tools.define_scale_tensor(df)
            results.append(theory_tools.make_scale_hist(df, [*nominal_axes, axis_ptVgen, axis_chargeVgen], [*nominal_cols, "ptVgen", "chargeVgen"]))
            if isW and not args.skipHelicity:
                helicity_helper = qcdScaleByHelicity_Whelper
                # TODO: Should have consistent order here with the scetlib correction function
                df = df.Define("helicityWeight_tensor", helicity_helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
                qcdScaleByHelicityUnc = df.HistoBoost("qcdScaleByHelicity", [*nominal_axes, axis_ptVgen, axis_chargeVgen], [*nominal_cols, "ptVgen", "chargeVgen", "helicityWeight_tensor"], tensor_axes=helicity_helper.tensor_axes)
                results.append(qcdScaleByHelicityUnc)

            df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)
            results.extend(theory_tools.make_pdf_hists(df, dataset.name, nominal_axes, nominal_cols, args.pdfs))

            masswargs = (nominal_axes, nominal_cols) if isW else (None, None)
            df, masswhist = syst_tools.define_mass_weights(df, isW, *masswargs)
            if masswhist:
                results.append(masswhist)

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                # TODO: Move to syst_tools
                netabins = args.muonCorrEtaBins
                nweights = 21
                mag = args.muonCorrMag
                df = df.Define(f"muonScaleDummy{netabins}Bins", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, goodMuons_eta0, {mag}, {str(isW).lower()})")
                scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
                dummyMuonScaleSyst = df.HistoBoost("muonScaleSyst", nominal_axes, [*nominal_cols, f"muonScaleDummy{netabins}Bins"],
                                                   tensor_axes=[common.down_up_axis, scale_etabins_axis])

                results.append(dummyMuonScaleSyst)

            df = df.Define("Muon_cvhbsMomCov", "wrem::splitNestedRVec(Muon_cvhbsMomCov_Vals, Muon_cvhbsMomCov_Counts)")

            df = df.Define("muonScaleSyst_responseWeights_tensor", calibration_uncertainty_helper,
                           ["Muon_correctedPt",
                            "Muon_correctedEta",
                            "Muon_correctedPhi",
                            "Muon_correctedCharge",
                            "Muon_genPartIdx",
                            "Muon_cvhbsMomCov",
                            "vetoMuonsPre",
                            "GenPart_pt",
                            "GenPart_eta",
                            "GenPart_phi",
                            "GenPart_pdgId",
                            "GenPart_statusFlags",
                            "nominal_weight"])

            dummyMuonScaleSyst_responseWeights = df.HistoBoost("muonScaleSyst_responseWeights", nominal_axes, [*nominal_cols, "muonScaleSyst_responseWeights_tensor"], tensor_axes = calibration_uncertainty_helper.tensor_axes)
            results.append(dummyMuonScaleSyst_responseWeights)



    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
output_tools.write_analysis_output(resultdict, "mw_with_mu_eta_pt.pkl.lz4", args)
