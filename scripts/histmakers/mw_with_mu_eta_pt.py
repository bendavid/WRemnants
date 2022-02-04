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


#TODO add the right arguments here
ROOT.wremnants.initializeScaleFactors(wremnants.data_dir, wremnants.data_dir + "/testMuonSF/scaleFactorProduct_28Oct2021_nodz_dxybs_genMatchDR01.root") #other args

datasets = wremnants.datasets2016.allDatasets(istest=False)


# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]


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

    df = df.Define("transverseMass", "wremnants::mt_2(goodMuons_pt0, goodMuons_phi0, MET_pt, MET_phi)")

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")

    df = df.Define("goodCleanJets", "Jet_jetId >= 6 && (Jet_pt > 50 || Jet_puId >= 4) && Jet_pt > 30 && abs(Jet_eta) < 2.4 && wremnants::cleanJetsFromLeptons(Jet_eta,Jet_phi,Muon_eta[vetoMuons],Muon_phi[vetoMuons],Electron_eta[vetoElectrons],Electron_phi[vetoElectrons])")

    df = df.Define("passMT", "transverseMass >= 40.0")

    df = df.Filter("passMT || Sum(goodCleanJets)>=1")

    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    if dataset.is_data:
        nominal = df.HistoBoost("nominal", nominal_axes, nominal_cols)
        results.append(nominal)
    else:
        #TODO what's the right tag here?
        df = df.DefinePerSample("eraVFP", "wremnants::GToH")

        df = df.Define("weight_pu", "wremnants::puw_2016UL_era(Pileup_nTrueInt,eraVFP)")
        df = df.Define("weight_fullMuonSF", "wremnants::_get_fullMuonSF(goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)")
        df = df.Define("weight_newMuonPrefiringSF", "wremnants::_get_newMuonPrefiringSF(Muon_eta,Muon_pt,Muon_phi,Muon_looseId,eraVFP)")
        df = df.Define("weight_tnpRecoSF", "wremnants::_get_tnpRecoSF(goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,0, wremnants::reco)")
        df = df.Define("weight_tnpTrackingSF", "_get_tnpTrackingSF(goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP)")

        df = df.Define("nominal_weight", "weight*weight_pu*weight_fullMuonSF*weight_newMuonPrefiringSF*weight_tnpRecoSF*weight_tnpTrackingSF")

        nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "nominal_weight"])
        results.append(nominal)

        # slice 101 elements starting from 0 and clip values at += 10.0
        df = df.Define("pdfWeights_tensor", "wremnants::clip_tensor(wremnants::vec_to_tensor_t<float, 101>(LHEPdfWeight), 10.)")

        pdfNNPDF31 = df.HistoBoost("pdfNNPDF31", nominal_axes, [*nominal_cols, "pdfWeights_tensor"])
        results.append(pdfNNPDF31)

        # slice 2 elements starting from 101
        df = df.Define("pdfWeightsAS_tensor", "wremnants::vec_to_tensor_t<double, 2>(LHEPdfWeight, 101)")

        alphaS002NNPDF31 = df.HistoBoost("alphaS002NNPDF31", nominal_axes, [*nominal_cols, "pdfWeightsAS_tensor"])
        results.append(alphaS002NNPDF31)

        #TODO check that I actually have a consistent version of functions.h and functionsWMass.h
        #TODO convert _get_fullMuonSFvariation_splitIso to produce a tensor natively
        df = df.Define("effStatTnP_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<1248>> res = (nominal_weight/weight_fullMuonSF)*wremnants::vec_to_tensor_t<double, 1248>(wremnants::_get_fullMuonSFvariation_splitIso(624, goodMuons_pt0 ,goodMuons_eta0,goodMuons_charge0,-1,-1,eraVFP,passIso)); return res;")

        effStatTnP = df.HistoBoost("effStatTnP", nominal_axes, [*nominal_cols, "effStatTnP_tensor"])
        results.append(effStatTnP)








    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

fname = "mw_with_mu_eta_pt.pkl.lz4"

print("writing output")
#with gzip.open(fname, "wb") as f:
with lz4.frame.open(fname, "wb") as f:
    pickle.dump(resultdict, f)
