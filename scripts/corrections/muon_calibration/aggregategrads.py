import os

import ROOT

import narf.matrix_utils
import wremnants.production
from wremnants.production.module_corrections import (
    book_grad_helper,
    book_hess_helper_sparse,
)

# os.environ["XRD_PARALLELEVTLOOP"] = "24"
# os.environ["ROOT_TTREE_CACHE_PREFETCH"] = "1"


# ROOT.gEnv.SetValue("TFile.AsyncPrefetching", 1)

# ROOT.gInterpreter.ProcessLine(".O3")
# ROOT.ROOT.EnableImplicitMT(64)
ROOT.ROOT.EnableImplicitMT()


import h5py
import hdf5plugin
import numpy as np

# from utils import lumitools
import narf
import narf.lumitools

# assert(0)

# hlt_paths = ['HLT_Dimuon20_Jpsi','HLT_DoubleMu4_JpsiTrk_Displaced','HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing','HLT_Mu7p5_Track2_Jpsi','HLT_Mu7p5_Track3p5_Jpsi','HLT_Dimuon0_Jpsi_Muon','HLT_Dimuon0er16_Jpsi_NoVertexing','HLT_Dimuon10_Jpsi_Barrel','HLT_Dimuon16_Jpsi','HLT_DoubleMu4_3_Jpsi_Displaced','HLT_Mu7p5_Track7_Jpsi']

hlt_paths = [
    "HLT_Dimuon20_Jpsi",
    "HLT_DoubleMu4_JpsiTrk_Displaced",
    "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
    "HLT_Mu7p5_Track2_Jpsi",
    "HLT_Mu7p5_Track3p5_Jpsi",
    "HLT_Dimuon0er16_Jpsi_NoVertexing",
    "HLT_Dimuon10_Jpsi_Barrel",
    "HLT_Dimuon16_Jpsi",
    "HLT_DoubleMu4_3_Jpsi_Displaced",
    "HLT_Mu7p5_Track7_Jpsi",
]


chainjpsi = ROOT.TChain("tree")
chainjpsi.Add(
    "/scratch/submit/cms/wmass/muoncal2/MuonGunUL2016_v719_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_grads/230112_043259/0000/*.root"
)

wremdir = os.environ["WREM_BASE"]

jsonhelper = narf.lumitools.make_jsonhelper(
    f"{wremdir}/wremnants-data/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
)


# print("first file:",
# chainjpsi.GetListOfFiles()[0])

filenameinfo = chainjpsi.GetListOfFiles()[0].GetTitle()
finfo = ROOT.TFile.Open(filenameinfo)
runtree = finfo.Get("runtree")
nparms = int(runtree.GetEntries())

dj = ROOT.ROOT.RDataFrame(chainjpsi)

cols = dj.GetColumnNames()
print(cols)
# quit()


ROOT.RDF.Experimental.AddProgressBar(dj)


# print(dj.Sum("edmvalref").GetValue())
# print(dj.Max("gradmax").GetValue())

# assert(0)

dj = dj.Filter(jsonhelper, ["run", "lumi"], "jsonfilter")

# dj = dj.Filter(" || ".join(hlt_paths))


# dj = dj.Filter("Mupluscons_pt > 1.1 && Muminuscons_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.8 && Jpsi_mass<3.4");

# dj = dj.Filter("Mupluscons_pt > 1.1 && Muminuscons_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Mupluscons_pt > 1.5 && Muminuscons_pt > 1.5 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonLoose && Muminus_muonLoose");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonLoose && Muminus_muonLoose");

# dj = dj.Filter("Muplus_pt > 4.0 && Muminus_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && (Muplus_muonIsGlobal || Muplus_muonIsTracker) && (Muminus_muonIsGlobal || Muminus_muonIsTracker)");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonLoose && Muminus_muonLoose");

# dj = dj.Filter("Mupluscons_pt > 1.1 && (Mupluscons_pt > 4.0 || abs(Mupluscons_eta) > 1.2) && Muminuscons_pt > 1.1 && (Muminuscons_pt > 4.0 || abs(Muminuscons_eta) > 1.2) && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.8 && Jpsi_mass<3.3 && Muplus_muonLoose && Muminus_muonLoose");

# dj = dj.Filter("Muplus_pt > 1.1 && (Muplus_pt > 4.0 || abs(Muplus_eta) > 1.2) && Muminus_pt > 1.1 && (Muminus_pt > 4.0 || abs(Muminus_eta) > 1.2) && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.8 && Jpsi_mass<3.3 && Muplus_muonLoose && Muminus_muonLoose")


# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1. && Muminus_pt > 1. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0");

# dj = dj.Filter("Muplus_pt > 1.1 && (Muplus_pt > 4.0 || abs(Muplus_eta) > 1.2) && Muminus_pt > 1.1 && (Muminus_pt > 4.0 || abs(Muminus_eta) > 1.2) && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 5 && Muplus_nvalidpixel>1 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Jpsi_pt>5.5 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 10 && Muplus_nvalidpixel>1 && Muminus_nvalid > 10 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && (Muplus_muonIsGlobal || Muplus_muonIsTracker || abs(Muplus_eta)>1.2) && (Muminus_muonIsGlobal || Muminus_muonIsTracker || abs(Muminus_eta)>1.2)");

# dj = dj.Filter("Muplus_pt > 1. && Muminus_pt > 1. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1. && Muminus_pt > 1. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");

# dj = dj.Filter("Muplus_pt > 1. && Muminus_pt > 1. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");
# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplusgen_pt > 1.1 && Muminusgen_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");
# dj = dj.Filter("Muplus_pt > 3.5 && Muminus_pt > 3.5 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");


# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");


# dj = dj.Filter("Mupluscons_pt > 1.1 && Muminuscons_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Muplus_pt > 5. && Muminus_pt > 5. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1.3 && Muminus_pt > 1.3 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9 && Muplus_nvalid > 8 && Muplus_nvalidpixel>0 && Muminus_nvalid > 8 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplusgen_pt > 0.9 && Muminusgen_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0");

# dj = dj.Filter("Muplusgen_pt > 0.9 && Muminusgen_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0");

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0");

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass > 2.85 && Jpsi_mass < 3.35");

# dj = dj.Filter("Mupluscons_pt > 5.0 && Muminuscons_pt > 5.0 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass > 2.85 && Jpsi_mass < 3.35");

# dj = dj.Filter("Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass > 2.9 && Jpsi_mass < 3.3");

dj = dj.Filter("Jpsi_pt > 8.2")

dj = dj.Filter("Muplus_pt > 4.0 && Muminus_pt > 4.0")

dj = dj.Filter(
    "Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass > 2.92 && Jpsi_mass < 3.28"
)


# dj = dj.Filter("Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsigen_mass > 2.9 && Jpsigen_mass < 3.3");


# dj = dj.Filter("Muplusgen_pt > 5.0 && Muminusgen_pt > 5.0")
# dj = dj.Filter("Mupluscons_pt > 5.0 && Muminuscons_pt > 5.0")
# dj = dj.Filter("Muplus_pt > 5.0 && Muminus_pt > 5.0")
# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9")


# dj = dj.Filter("Muplus_pt/Jpsi_mass > (4.0/2.9) && Muminus_pt/Jpsi_mass > (4.0/2.9)")
# dj = dj.Filter("Mupluscons_pt/Jpsicons_mass > (4.0/2.9) && Muminuscons_pt/Jpsicons_mass > (4.0/2.9)")

# dj = dj.Filter("Muplus_pt > 4.0 && Muminus_pt > 4.0")
# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0")
# dj = dj.Filter("Muplusgen_pt > 4.0 && Muminusgen_pt > 4.0")
# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9")
# dj = dj.Filter("Muplus_pt > 4.86 && Muminus_pt > 4.86")
# dj = dj.Filter("Muplus_pt < 13.0 && Muminus_pt < 13.0")
# dj = dj.Filter("std::fabs(Muplus_eta) < 2.4 && std::fabs(Muminus_eta) < 2.4")

##dj = dj.Filter("Mupluscons_pt > 5.0 && Muminuscons_pt > 5.0")
# dj = dj.Filter("Muplus_pt/Jpsi_mass > 4.5/2.85 && Muminus_pt/Jpsi_mass > 4.5/2.85")
# dj = dj.Filter("Mupluscons_pt/Jpsi_mass > 5.0/2.85 && Muminuscons_pt/Jpsi_mass > 5.0/2.85")

# dj = dj.Filter("Mupluscons_pt > 1.1 && Muminuscons_pt > 1.1")

# dj = dj.Filter("(std::fabs(Mupluscons_eta) < 1.2 && Mupluscons_pt > 4.5) || (std::fabs(Mupluscons_eta) >= 1.2 && Mupluscons_pt > (-2.5*std::fabs(Mupluscons_eta)/1.2 + 7.0))")
# dj = dj.Filter("(std::fabs(Muminuscons_eta) < 1.2 && Muminuscons_pt > 4.5) || (std::fabs(Muminuscons_eta) >= 1.2 && Muminuscons_pt > (-2.5*std::fabs(Muminuscons_eta)/1.2 + 7.0))")


# dj = dj.Filter("dmassconvval < 6e-4")

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("abs(Muplus_eta)>1.2 || Muplus_pt>5.0")
# dj = dj.Filter("abs(Muminus_eta)>1.2 || Muminus_pt>5.0")


# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("abs(Mupluscons_eta)>1.2 || Mupluscons_pt>5.0")
# dj = dj.Filter("abs(Muminuscons_eta)>1.2 || Muminuscons_pt>5.0")

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 10 && Muplus_nvalidpixel>1 && Muminus_nvalid > 10 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 0.9 && Muminus_pt > 0.9 && Muplus_nvalid > 8 && Muplus_nvalidpixel>1 && Muminus_nvalid > 8 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 5.0 && Muminus_pt > 5.0 && Muplus_nvalid > 8 && Muplus_nvalidpixel>1 && Muminus_nvalid > 8 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Mupluscons_pt > 5.0 && Muminuscons_pt > 5.0 && Muplus_nvalid > 8 && Muplus_nvalidpixel>1 && Muminus_nvalid > 8 && Muminus_nvalidpixel > 1 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");
# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && abs(Mupluscons_eta)<2.4 && abs(Muminuscons_eta)<2.4 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");
# dj = dj.Filter("Mupluscons_pt > 0.9 && Muminuscons_pt > 0.9 && abs(Mupluscons_eta)<2.4 && abs(Muminuscons_eta)<2.4 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal");
# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && Muplus_muonIsGlobal && Muminus_muonIsGlobal && std::fabs(Mupluscons_eta) < 1. && std::fabs(Muminuscons_eta) < 1.");


# dj = dj.Filter("Mupluscons_pt > 4.0 && Muminuscons_pt > 4.0 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");

# dj = dj.Filter("Muplus_pt > 1.1 && Muminus_pt > 1.1 && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3 && (Muplus_muonIsGlobal || fabs(Muplus_eta) > 1.2) && (Muminus_muonIsGlobal || fabs(Muplus_eta) > 1.2)");


# dj = dj.Filter("std::fabs(Mupluscons_eta) < 2.4 && std::fabs(Muminuscons_eta) < 2.4")
# dj = dj.Filter("std::fabs(Jpsi_z) < 10.")
# dj = dj.Filter("maxelement(gradv) < 1e5")

# dj = dj.Filter("Muplus_pt > 1.5 && Muminus_pt > 1.5 && Muplus_pt < 23. && Muminus_pt < 23. && Muplus_nvalid > 3 && Muplus_nvalidpixel>0 && Muminus_nvalid > 3 && Muminus_nvalidpixel > 0 && Jpsi_mass>2.9 && Jpsi_mass<3.3");


# dj = dj.Filter("Muplus_muonLoose && Muminus_muonLoose")

# dj = dj.Filter("dmassconvval >= 0. && dmassconvval/3.09692 < 1e-3")

# dj = dj.Filter("Muplus_highpurity && Muminus_highpurity")

dj = dj.Filter("edmvalref < 1e-5")
# dj = dj.Filter("edmvalref < 1e-5 && edmvalref >= 0.")
# dj = dj.Filter("chisqval/ndof < 3.")
# dj = dj.Filter("chisqval/ndof < 10.")
# dj = dj.Filter("chisqval/ndof < 3.")


# dj = dj.Filter("Jpsigen_mass > 3.0968")
# dj = dj.Filter("Jpsigen_mass>3.0960 && Jpsigen_mass<3.0978")


# dj = dj.Filter("Jpsigen_mass > 3.096 && Jpsigen_mass>0. && Jpsigen_mass<3.0978")

# dj = dj.Filter("gradmax < 1e5")

# dj = dj.Filter("valid(gradv) && valid(hesspackedv)");
# dj = dj.Filter("valid(gradv) && valid(hesspackedv) && abs(deltachisqval) < 1e-2");
# dj = dj.Filter("valid(gradv) && valid(hesspackedv) && edmval < 1e-5");
# dj = dj.Filter("valid(gradv) && valid(hesspackedv) && edmval_cons0 < 1e-5");


# dj = dj.Define("massweightval", "Numba::massweight(Muplus_pt, Muplus_eta, Muplus_phi, Muminus_pt, Muminus_eta, Muminus_phi, Jpsi_mass, run)")


dj = dj.Define("massweightval", "1.0")

dj = dj.Filter("massweightval > 0.")


djcount = dj.Count()
maxgradient = dj.Max("gradmax")


print("nparms", nparms)
# quit()

grad_res = book_grad_helper(
    dj, nparms, ["gradv", "globalidxv", "massweightval"]
)
hess_res = book_hess_helper_sparse(
    dj, nparms, ["hesspackedv", "globalidxv", "massweightval"]
)


# gradval = grad.GetResult()
# hessval = hess.GetResult()

# print(gradval[0])

print("djcount", djcount.GetValue())
print("maxgradient", maxgradient.GetValue())


# chunksize = 32

# fout = h5py.File("combinedgrads.hdf5", "w", rdcc_nbytes = nparms*8*chunksize*4, rdcc_nslots = nparms//chunksize*10)

grad = grad_res.GetValue()  # numpy 1-D array, length nparms
hess = hess_res.GetValue()  # scipy CSR (nparms, nparms), symmetric

import wums.ioutils

with h5py.File("combinedgrads.hdf5", "w") as fout:
    wums.ioutils.pickle_dump_h5py("grad", grad, fout)
    wums.ioutils.pickle_dump_h5py("hess", hess, fout)
