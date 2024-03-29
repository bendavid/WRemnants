import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
args = parser.parse_args()

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
if args.nThreads is not None:
    if args.nThreads > 1:
        ROOT.ROOT.EnableImplicitMT(args.nThreads)
else:
    ROOT.ROOT.EnableImplicitMT()


import pickle
import gzip

import narf
import wremnants
import hist
import lz4.frame
import numba


ROOT.wrem.initializeScaleFactors(wremnants.data_dir, wremnants.data_dir + "/testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root")

datasets = wremnants.datasets2016.allDatasets(istest=False)

era = "GToH"

muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

wprocs = ["WplusmunuPostVFP", "WminusmunuPostVFP"]
zprocs = ["ZmumuPostVFP"]

# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

# extra axes which can be used to label tensor_axes

down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "down-up variation")

down_nom_up_axis = hist.axis.Regular(3, -1.5, 1.5, underflow=False, overflow=False, name = "down-nominal-up variation")


def build_graph(df, dataset):
    results = []

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    df = df.Define("vetoMuons", "Muon_pt > 10 && Muon_looseId && abs(Muon_eta) < 2.4 && abs(Muon_dxybs) < 0.05")

    df = df.Filter("Sum(vetoMuons) == 1")

    df = df.Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_isGlobal")

    df = df.Filter("Sum(goodMuons) == 1")

    df = df.Define("goodMuons_pt0", "Muon_pt[goodMuons][0]")
    df = df.Define("goodMuons_eta0", "Muon_eta[goodMuons][0]")
    df = df.Define("goodMuons_phi0", "Muon_phi[goodMuons][0]")
    df = df.Define("goodMuons_charge0", "Muon_charge[goodMuons][0]")

    df = df.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")

    df = df.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_pt, MET_phi)")

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")

    df = df.Define("goodCleanJets", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && Jet_pt > 30 && abs(Jet_eta) < 2.4 && wrem::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_eta[vetoMuons],Muon_phi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")

    df = df.Define("passMT", "transverseMass >= 40.0")

    df = df.Filter("passMT || Sum(goodCleanJets)>=1")

    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if dataset.is_data:
        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)
    else:
        df = df.DefinePerSample("eraVFP", f"wrem::{era}")

        df = df.Define("weight_pu", "wrem::puw_2016UL_era(Pileup_nTrueInt,eraVFP)")
        df = df.Define("weight_fullMuonSF", "wrem::_get_fullMuonSF(goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)")
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_eta", "Muon_pt", "Muon_phi", "Muon_looseId"])
        df = df.Define("weight_tnpTrackingRecoSF", "wrem::_get_tnpTrackingRecoSF(goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP)")

        df = df.Define("nominal_weight", "weight*weight_pu*weight_fullMuonSF*weight_newMuonPrefiringSF*weight_tnpTrackingRecoSF")

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

        # slice 101 elements starting from 0 and clip values at += 10.0
        # extra assignment is to force the correct return type
        df = df.Define("pdfWeights_tensor", "auto res = wrem::clip_tensor(wrem::vec_to_tensor_t<double, 101>(LHEPdfWeight), 10.); res = nominal_weight*res; return res;")

        pdfNNPDF31 = df.HistoBoost("pdfNNPDF31", nominal_axes, [*nominal_cols, "pdfWeights_tensor"])
        results.append(pdfNNPDF31)

        # slice 2 elements starting from 101
        # extra assignment is to force the correct return type
        df = df.Define("pdfWeightsAS_tensor", "auto res = wrem::vec_to_tensor_t<double, 2>(LHEPdfWeight, 101); res = nominal_weight*res; return res;")

        alphaS002NNPDF31 = df.HistoBoost("alphaS002NNPDF31", nominal_axes, [*nominal_cols, "pdfWeightsAS_tensor"])
        results.append(alphaS002NNPDF31)

        #TODO convert _get_fullMuonSFvariation_splitIso to produce a tensor natively

        # extra assignment is to force the correct return type
        df = df.Define("effStatTnP_tensor", "auto res = wrem::vec_to_tensor_t<double, 1248>(wrem::_get_fullMuonSFvariation_splitIso(624, goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)); res = (nominal_weight/weight_fullMuonSF)*res; return res;")

        effStatTnP = df.HistoBoost("effStatTnP", nominal_axes, [*nominal_cols, "effStatTnP_tensor"])
        results.append(effStatTnP)

        df = df.Define("effSystIsoTnP_weight", "(nominal_weight/weight_fullMuonSF)*_get_fullMuonSF_dataAltSig_splitIso(1,goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)")

        effSystIsoTnP = df.HistoBoost("effSystIsoTnP", nominal_axes, [*nominal_cols, "effSystIsoTnP_weight"])
        results.append(effSystIsoTnP)

        df = df.Define("effSystTrigAndIdipTnP_weight", "(nominal_weight/weight_fullMuonSF)*_get_fullMuonSF_dataAltSig_splitIso(0,goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)")

        effSystTrigAndIdipTnP = df.HistoBoost("effSystTrigAndIdipTnP", nominal_axes, [*nominal_cols, "effSystTrigAndIdipTnP_weight"])
        results.append(effSystIsoTnP)

        #FIXME skipping EffTrackingRecoTnP_ since it's not consistently defined yet

        df = df.Define("muonL1PrefireStat_tensor", muon_prefiring_helper_stat, ["Muon_eta", "Muon_pt", "Muon_phi", "Muon_looseId", "nominal_weight"])

        muonL1PrefireStat = df.HistoBoost("muonL1PrefireStat", nominal_axes, [*nominal_cols, "muonL1PrefireStat_tensor"], tensor_axes = ["muon prefiring eta-phi regions", down_up_axis])
        results.append(muonL1PrefireStat)

        df = df.Define("muonL1PrefireSyst_tensor", muon_prefiring_helper_syst, ["Muon_eta", "Muon_pt", "Muon_phi", "Muon_looseId", "nominal_weight"])

        muonL1PrefireSyst = df.HistoBoost("muonL1PrefireSyst", nominal_axes, [*nominal_cols, "muonL1PrefireSyst_tensor"], tensor_axes = [down_up_axis])
        results.append(muonL1PrefireSyst)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples
        if dataset.name in wprocs:
            # extra assignment is to force the correct return type
            df = df.Define("massWeight_tensor", "auto res = wrem::vec_to_tensor_t<double, 21>(MEParamWeight); res = nominal_weight*res; return res;")

            massWeight = df.HistoBoost("massWeight", nominal_axes, [*nominal_cols, "massWeight_tensor"])
            results.append(massWeight)

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

fname = "mw_with_mu_eta_pt.pkl.lz4"

print("writing output")
#with gzip.open(fname, "wb") as f:
with lz4.frame.open(fname, "wb") as f:
    pickle.dump(resultdict, f)
