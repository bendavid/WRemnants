import os

# import wums.fitutils
import hist
import numpy as np
import ROOT

import narf
import narf.lumitools
import wremnants
import wremnants.production.datasets.dataset_tools
import wremnants.production.muon_calibration
import wremnants.production.pileup
import wremnants.production.vertex
import wremnants.utilities

# import matplotlib.pyplot as plt
from wremnants.production.histmaker_tools import write_analysis_output
from wremnants.utilities import common, parsing

# os.environ["XRD_NETWORKSTACK"] = "IPv4"
# os.environ["XRD_PARALLELEVTLOOP"] = "24"


analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

args = parser.parse_args()


ROOT.ROOT.EnableImplicitMT()
# ROOT.ROOT.EnableImplicitMT(128)


print("done declarations")

# mode = "z"
mode = "jpsi"


# hlt_paths_weights = ['HLT_Dimuon20_Jpsi','HLT_DoubleMu4_JpsiTrk_Displaced','HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing','HLT_Mu7p5_Track2_Jpsi','HLT_Mu7p5_Track3p5_Jpsi','HLT_Dimuon0_Jpsi_Muon','HLT_Dimuon0er16_Jpsi_NoVertexing','HLT_Dimuon10_Jpsi_Barrel','HLT_Dimuon16_Jpsi','HLT_DoubleMu4_3_Jpsi_Displaced','HLT_Mu7p5_Track7_Jpsi']


hlt_paths = ["HLT_Dimuon20_Jpsi"]

# mmin = 2.97
# mmax = 3.23

netabins = 24
nphibins = 16

axis_eta = hist.axis.Regular(netabins, -2.4, 2.4, name="eta")
axis_phi = hist.axis.Regular(nphibins, 0.0, 2.0 * np.pi, circular=True, name="phi")


mmin = 2.92
mmax = 3.28
muptmin = 4.0

axis_mass = hist.axis.Regular(20, mmin, mmax, name="mass")


print("definining chains")


def chain_to_dataset(chain, name):
    filepaths = [str(f.GetTitle()) for f in chain.GetListOfFiles()]
    # filepaths = filespaths[:2]
    return narf.Dataset(name, filepaths=filepaths)


eoscms = "eoscms.cern.ch"
# eoscms = "128.142.52.33"

# chainideal = ROOT.TChain("tree")
# chainideal.Add("/scratch/submit/cms/wmass/muoncalproc/MuonGunUL2016_v725_RecJpsiPythiaPhotosPt0to8WithExt1_quality_novtx_noconstraint_idealgeom_corgensim/jpsicor.root")
# chainideal.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_idealgeom/230214_153350/0000/globalcor_0_1.root")
# chainideal.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_idealgeom/230214_153350/0000/*.root")
# chainideal.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint_idealgeom/230214_153512/0000/*.root")
# chainideal.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint_idealgeom/230214_153512/0001/*.root")

pathideal = f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_idealgeom/230214_153350/"
fideal = wremnants.production.datasets.dataset_tools.buildFileListXrd(
    pathideal, num_clients=64
)
# print(fideal)

# quit()


# chainnom = ROOT.TChain("tree")
# chainnom.Add("/scratch/submit/cms/wmass/muoncal2/MuonGunUL2016_v719_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_grads/230112_043259/0000/*.root")
# chainnom.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint/230214_151859/0000/globalcor_0_1.root")
# chainnom.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint/230214_151859/0000/*.root")
# chainnom.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint/230214_152107/0000/*.root")
# chainnom.Add(f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint/230214_152107/0001/*.root")

pathnom = f"root://{eoscms}//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint/230214_151859/"
fnom = wremnants.production.datasets.dataset_tools.buildFileListXrd(
    pathnom, num_clients=64
)

print(fideal)
print(fnom)


# filenameinfo = chainnom.GetListOfFiles()[0].GetTitle()
filenameinfo = fideal[0]
finfo = ROOT.TFile.Open(filenameinfo)
finfo.ls()
runtree = finfo.Get("runtree")
print(runtree)
nparms = int(runtree.GetEntries())
finfo.Close()

axis_corparms = hist.axis.Integer(
    0, nparms, underflow=False, overflow=False, name="corparms"
)

# chain.Add("/scratch/submit/cms/wmass/muoncal2/MuonGunUL2016_v719_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_grads/230112_043259/0000/globalcor_0_112.root")

# dmassconvval/mass 3e-5 - 6e-4

# dataset_ideal = chain_to_dataset(chainideal, "jpsi_ideal")
# dataset_nom = chain_to_dataset(chainnom, "jpsi_nom")

dataset_ideal = narf.Dataset("jpsi_ideal", fideal)
dataset_nom = narf.Dataset("jpsi_ideal", fnom)

datasets = [dataset_ideal, dataset_nom]

# print(dataset_ideal.filepaths)
# print(dataset_nom.filepaths)

# quit()

era = "2016PostVFP"
# data_dir = wremnants.utilities.common.data_dir
# jsonhelper = narf.lumitools.make_jsonhelper(f"{data_dir}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")

pileup_helper = wremnants.production.pileup.make_pileup_helper(era=era)
vertex_helper = wremnants.production.vertex.make_vertex_helper(era=era)


print("defining helpers")

# corfilenom = "/eos/cms/store/cmst3/group/wmass/muoncal_corrections/correctionResults_v718_recjpsi.root"


# corfilenom = "/eos/cms/store/cmst3/group/wmass/muoncal_corrections/correctionResults_v722_recjpsiextnommc.root"
# corfilenom = "/eos/cms/store/cmst3/group/wmass/muoncal_corrections/correctionResults_v721_recjpsidata.root"

# corfilenom = "correctionResults_v722_recjpsiextnommc.root"
# corfilenom = "correctionResults_v721_recjpsidata.root"
#
helper = wremnants.production.muon_calibration.make_muon_calibration_helper_single()
# helper_nom = wremnants.muon_calibration.make_muon_calibration_helper_single(corfilenom)

# axis_eta_corr = hist.axis.Regular(netabins, -2.4, 2.4, name="eta_corr")
# axis_corr = hist.axis.StrCategory(["A", "e", "M"], overflow=False)
# hcorr = hist.Hist(axis_eta_corr, axis_corr)
# hcorr_cpp = narf.histutils.hist_to_pyroot_boost(hcorr, tensor_rank=1)
# helper_corr = ROOT.wrem.JpsiCorrectionsHelper[type(hcorr_cpp)](ROOT.std.move(hcorr_cpp))

print("done defining helpers")

# nparmslocal = 5
# ncorparmslocal = 3


def build_graph_base(df, dataset, max_events=-1):

    if max_events > 0:
        df = df.Filter(f"rdfentry_ < {max_events}")

    df = df.DefinePerSample("weight", "1.0")

    weightsum = df.SumAndCount("weight")

    if dataset.is_data:
        df = df.DefinePerSample("nominal_weight", "weight")
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["Jpsigen_z", "Pileup_nTrueInt"])
        df = df.Define("nominal_weight", "weight*weight_pu*weight_vtx")

    is_ideal = dataset.name == "jpsi_ideal"

    if is_ideal:
        print("ideal: applying corrections")
        df = wremnants.production.muon_calibration.define_lbl_corrections_jpsi_calibration_ntuples(
            df, helper
        )
    else:
        print("nominal: no corrections")
        df = wremnants.production.muon_calibration.define_passthrough_corrections_jpsi_calibration_ntuples(
            df
        )

    df = df.Filter("std::fabs(Mupluscor_eta) < 2.4 && std::fabs(Muminuscor_eta) < 2.4")

    df = df.Filter(f"Mupluscor_pt > {muptmin}")
    df = df.Filter(f"Muminuscor_pt > {muptmin}")
    df = df.Filter("Jpsicor_pt > 8.2")

    df = df.Filter(f"Jpsicor_mass > {mmin} && Jpsicor_mass < {mmax}")

    return df, weightsum


doquantiles = False

if doquantiles:
    chainideal = ROOT.TChain("tree")
    dfideal = ROOT.ROOT.RDataFrame(chainideal)
    ROOT.ROOT.RDF.Experimental.AddProgressBar(dfideal)

    dfideal, weightsum = build_graph_base(dfideal, dataset_ideal, max_events=int(10e6))

    axis_pt_quant = hist.axis.Regular(
        5, 0.0, 1.0, underflow=False, overflow=False, name="pt_quant"
    )

    count = dfideal.Count()
    # TODO the dynamic quantile binning can work in principle, but to be compatible with continuous shifts
    # would need to be modified to return continous quantiles instead of just integers
    # use simple quantiles for now instead
    # quantile_hists = narf.histutils.build_quantile_hists(dfideal,
    #                                               cols = ["Muplus_eta", "Muplus_phi", "Muplus_pt"],
    #                                               condaxes = [axis_eta, axis_phi],
    #                                               quantaxes = [axis_pt_quant]
    # )
    quantile_hists = narf.histutils.build_quantile_hists(
        dfideal, cols=["Mupluscor_pt"], condaxes=[], quantaxes=[axis_pt_quant]
    )

    print("neventspost", count.GetValue())
    print(quantile_hists)
    ptquants = quantile_hists[0].values()
    print(ptquants)
    ptquants = [muptmin, *ptquants[:-1], np.inf]
    print(ptquants)
else:
    # [  4.58014965   5.2434926    6.22155809   8.08470631 451.72012329]
    ptquants = [muptmin, 4.58014965, 5.2434926, 6.22155809, 8.08470631, np.inf]


axis_pt = hist.axis.Variable(ptquants, underflow=True, overflow=False, name="pt")
print(axis_pt)

nominal_axes = [axis_eta, axis_phi, axis_pt, axis_mass]
nominal_cols_plus = ["Mupluscor_eta", "Mupluscor_phi", "Mupluscor_pt", "Jpsicor_mass"]
nominal_cols_minus = [
    "Muminuscor_eta",
    "Muminuscor_phi",
    "Muminuscor_pt",
    "Jpsicor_mass",
]


def build_graph(df, dataset):
    # df = df.DefinePerSample("weighttmp", "1.0")
    # weightsum = df.SumAndCount("weighttmp")
    # df = build_graph_base(df, dataset, max_events=int(10e6))

    df, weightsum = build_graph_base(df, dataset)

    results = []

    df = df.DefinePerSample("idxstest", "std::array{0, 1, 2, 3, 4, 5}")

    # df, quantile_axes_plus, quantile_cols_plus = narf.histutils.define_quantile_ints(df, ["Muplus_eta", "Muplus_phi", "Muplus_pt"], quantile_hists)
    # df, quantile_axes_minus, quantile_cols_minus = narf.histutils.define_quantile_ints(df, ["Muminus_eta", "Muminus_phi", "Muminus_pt"], quantile_hists)

    # hmuplus = df.HistoBoost("hmuplus",
    #                         axes = quantile_axes_plus + [axis_mass],
    #                         cols = quantile_cols_plus + ["Jpsi_mass"] + ["nominal_weight"])
    # results.append(hmuplus)

    # hmuminus = df.HistoBoost("hmuminus",
    #                         axes = quantile_axes_minus + [axis_mass],
    #                         cols = quantile_cols_minus + ["Jpsi_mass"] + ["nominal_weight"])
    # results.append(hmuminus)

    hmuplus = df.HistoBoost(
        "hmuplus", axes=nominal_axes, cols=nominal_cols_plus + ["nominal_weight"]
    )
    results.append(hmuplus)

    hmuminus = df.HistoBoost(
        "hmuminus", axes=nominal_axes, cols=nominal_cols_minus + ["nominal_weight"]
    )
    results.append(hmuminus)

    is_ideal = dataset.name == "jpsi_ideal"

    if is_ideal:
        hmuplus_corparms = df.HistoBoost(
            "hmuplus_corparms",
            axes=nominal_axes + [axis_corparms],
            cols=nominal_cols_plus + ["idxstest"] + ["nominal_weight"],
            storage=hist.storage.Double(),
        )
        results.append(hmuplus_corparms)

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree="tree")

print(resultdict)

fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
write_analysis_output(resultdict, fout, args)
