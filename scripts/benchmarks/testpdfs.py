import time
time0 = time.time()

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=["WplusmunuPostVFP"])
parser.add_argument("--useBoost", type=bool, help="use boost histograms", default=False)
parser.add_argument("--useTensor", type=bool, help="use tensors", default=False)
parser.add_argument("--nHists", type=int, help="number of hist copies", default=1)
parser.add_argument("--useSingles", type=bool, help="use tensors", default=False)


args = parser.parse_args()

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
if not args.nThreads:
    ROOT.ROOT.EnableImplicitMT()
elif args.nThreads != 1:
    ROOT.ROOT.EnableImplicitMT(args.nThreads)

#ROOT.gErrorIgnoreLevel = ROOT.kPrint
ROOT.ROOT.Experimental.RLogManager.Get().SetVerbosity(ROOT.ROOT.Experimental.ELogLevel.kInfo)

import pickle
import gzip

import narf
import wremnants
import hist
import lz4.frame
import logging

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None)

for dataset in datasets:
  print(dataset.name)


wprocs = ["WplusmunuPostVFP", "WminusmunuPostVFP", "WminustaunuPostVFP", "WplustaunuPostVFP"]
zprocs = ["ZmumuPostVFP", "ZtautauPostVFP"]

# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]
ptV_axis = hist.axis.Variable([0.0, 2.9, 4.7, 6.7, 9.0, 11.8, 15.3, 20.1, 27.2, 40.2, 13000.0], name="genPtV")

axis_pdfsyst_idx = hist.axis.Integer(0, 103, name = "pdfsyst_idx")

#assert(0)

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

    df = df.Define("transverseMass", "std::sqrt(2.*goodMuons_pt0*MET_pt*(1-std::cos(goodMuons_phi0-MET_phi)))");

    df = df.Define("passMT", "transverseMass >= 40.0")
    df = df.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

    nominal_cols = ["goodMuons_eta0", "goodMuons_pt0", "goodMuons_charge0", "passIso", "passMT"]

    events_pass = df.Count()
    results.append(events_pass)

    #nominal = df.HistoBoost("nominal", nominal_axes, [*nominal_cols, "weight"])
    #results.append(nominal)

    if args.useTensor:
        df = df.Define("pdfWeights_tensor", "auto res = wrem::vec_to_tensor_t<double, 103>(LHEPdfWeight); res = weight*res; return res;")
    else:
        if args.useSingles:
            for ipdf in range(103):
                df = df.Define(f"pdfWeights_{ipdf}", f"weight*LHEPdfWeight[{ipdf}];")
        else:
            df = df.DefinePerSample("pdfsyst_idx", "std::array<int, 103> res;  std::iota(res.begin(), res.end(), 0); return res;")

            df = df.Define("pdfWeights_rvec", "weight*LHEPdfWeight")

    for ihist in range(args.nHists):
        if args.useTensor:
            if not args.useBoost:
                raise ValueError("Invalid options")
            pdfNNPDF31 = df.HistoBoost(f"pdfNNPDF31_{ihist}", nominal_axes, [*nominal_cols, "pdfWeights_tensor"])
            results.append(pdfNNPDF31)
        else:
            if args.useBoost:
                #pdfNNPDF31 = df.HistoBoost(f"pdfNNPDF31_{ihist}", nominal_axes + [axis_pdfsyst_idx], nominal_cols + [ "pdfsyst_idx", "pdfWeights_rvec"])
                pdfNNPDF31 = df.HistoBoost(f"pdfNNPDF31_{ihist}", [axis_pdfsyst_idx] + nominal_axes , [ "pdfsyst_idx" ] + nominal_cols + [ "pdfWeights_rvec"])
                results.append(pdfNNPDF31)
            else:
                print("ROOT path")
                if args.useSingles:
                    for ipdf in range(103):
                        pdfNNPDF31 = df.HistoND((f"pdfNNPDF31_{ipdf}_{ihist}", "", 5, [48, 29, 2, 2, 2], [-2.4, 26.,-2.,-0.5, -0.5], [2.4, 55., 2., 1.5, 1.5]), nominal_cols + [ f"pdfWeights_{ipdf}"])
                        results.append(pdfNNPDF31)
                else:
                    #pdfNNPDF31 = df.HistoND((f"pdfNNPDF31_{ihist}", "", 6, [48, 29, 2, 2, 2, 103], [-2.4, 26.,-2.,-0.5, -0.5, -0.5], [2.4, 55., 2., 1.5, 1.5, 102.5]), nominal_cols + [ "pdfsyst_idx", "pdfWeights_rvec"])
                    pdfNNPDF31 = df.HistoND((f"pdfNNPDF31_{ihist}", "", 6, [103, 48, 29, 2, 2, 2], [-0.5, -2.4, 26.,-2.,-0.5, -0.5], [102.5, 2.4, 55., 2., 1.5, 1.5]), ["pdfsyst_idx"] + nominal_cols + ["pdfWeights_rvec"])
                    results.append(pdfNNPDF31)



    return results, weightsum

time_to_graph_setup = time.time()

resultdict = narf.build_and_run(datasets, build_graph)

print(resultdict)

print("init time:", time_to_graph_setup-time0)

