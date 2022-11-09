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
parser.add_argument("--nano", type=str, help="NanoAOD version to run on", default="tnp")
args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts]) 
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion=args.nano)

if not args.no_recoil:
    logging.warning("Recoil correction for high PU is not yet supported! Setting false")
    args.no_recoil = True

era = args.era
noMuonCorr = args.noMuonCorr

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

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

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

    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isTop = dataset.group == "Top"
    apply_theory_corr = args.theory_corr and dataset.name in corr_helpers
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

    df = df.Define("passMT", "1")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    if dataset.group in ["Top", "Diboson"]:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    else:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Filter("wrem::hasTriggerMatch(goodMuons_eta0,goodMuons_phi0,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")

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
        
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])

        weight_expr = "weight*weight_pu"
        if args.vertex_weight:
            weight_expr += "*weight_vtx"
            
        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
output_tools.write_analysis_output(resultdict, "mw_with_mu_eta_pt.pkl.lz4", args)
