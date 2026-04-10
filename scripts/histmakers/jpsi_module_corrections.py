import os

# os.environ["XRD_NETWORKSTACK"] = "IPv4"
# os.environ["XRD_PARALLELEVTLOOP"] = "24"


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
from wremnants.production.histmaker_tools import write_analysis_output
from wremnants.utilities import common, parsing

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

args = parser.parse_args()


ROOT.ROOT.EnableImplicitMT()
# ROOT.ROOT.EnableImplicitMT(64)


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
# muptmin = 4.0
muptmin = 6.2

# Continuous quantile-transformed pt and mass axes: the quantile helpers
# built below map the original pt / mass values to CDF-style values in
# [0, 1], so the nominal histogram axes are Regular over [0, 1].
npt_quant = 5
nmass_quant = 20
axis_pt_quant = hist.axis.Regular(
    npt_quant, 0.0, 1.0, underflow=False, overflow=False, name="pt_quant"
)
axis_mass_quant = hist.axis.Regular(
    nmass_quant, 0.0, 1.0, underflow=False, overflow=False, name="mass_quant"
)


def paths_to_filenames(paths):
    filenames = []
    for path in paths:
        filenames += wremnants.production.datasets.dataset_tools.buildFileListXrd(
            path,
            #    num_clients=64,
        )
    return filenames


pathsideal = []
pathsideal.append(
    "root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_idealgeom/230214_153350/"
)
# pathsideal.append("root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint_idealgeom/230214_153512/")
fideal = paths_to_filenames(pathsideal)


pathsnom = []
pathsnom.append(
    "root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint/230214_151859/"
)
# pathsnom.append("root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint/230214_152107/")
fnom = paths_to_filenames(pathsnom)


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

dataset_ideal = narf.Dataset("jpsi_ideal", fideal)
dataset_nom = narf.Dataset("jpsi_nom", fnom)

datasets = [dataset_ideal, dataset_nom]


era = "2016PostVFP"

pileup_helper = wremnants.production.pileup.make_pileup_helper(era=era)
vertex_helper = wremnants.production.vertex.make_vertex_helper(era=era)


print("defining helpers")
#

helper = wremnants.production.muon_calibration.make_muon_calibration_helper_single()
helper_var = ROOT.wrem.CVHCorrectorUncertainty[5]()


print("done defining helpers")


def build_graph_base(df, dataset, max_events=-1):

    if max_events > 0:
        df = df.Filter(f"rdfentry_ < {max_events}")

    df = df.DefinePerSample("weight", "1.0")

    weightsum = df.SumAndCount("weight")

    df = df.Filter("||".join(hlt_paths))

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

    # df = df.Filter(f"Mupluscor_pt > {muptmin}")
    # df = df.Filter(f"Muminuscor_pt > {muptmin}")
    df = df.Filter("max(Mupluscor_pt, Muminuscor_pt) > 13.2")
    df = df.Filter(f"min(Mupluscor_pt, Muminuscor_pt) > {muptmin}")
    df = df.Filter("Jpsicor_pt > 8.2")

    df = df.Filter(f"Jpsicor_mass > {mmin} && Jpsicor_mass < {mmax}")

    return df, weightsum


# Pre-pass over a subset of the ideal sample to build continuous quantile
# helpers for pt and mass. The helpers map raw pt / mass values to CDF-style
# values in [0, 1] via linear interpolation between the learned quantile
# edges, and are then applied (inside build_graph) to both the nominal muon
# columns and the shifted-variation RVecs.
dfideal_quant = ROOT.ROOT.RDataFrame("tree", fideal[:100])
ROOT.ROOT.RDF.Experimental.AddProgressBar(dfideal_quant)
dfideal_quant, _ = build_graph_base(dfideal_quant, dataset_ideal)

quantile_hists_pt, pt_centers_hist, pt_volume_hist = (
    narf.histutils.build_quantile_hists(
        dfideal_quant,
        cols=["Mupluscor_pt"],
        condaxes=[],
        quantaxes=[axis_pt_quant],
        continuous=True,
    )
)

quantile_hists_mass, mass_centers_hist, mass_volume_hist = (
    narf.histutils.build_quantile_hists(
        dfideal_quant,
        cols=["Jpsicor_mass"],
        condaxes=[],
        quantaxes=[axis_mass_quant],
        continuous=True,
    )
)

nominal_axes = [axis_eta, axis_phi, axis_pt_quant, axis_mass_quant]
nominal_cols_plus = [
    "Mupluscor_eta",
    "Mupluscor_phi",
    "Mupluscor_pt_quant",
    "Jpsicor_mass_quant",
]
nominal_cols_minus = [
    "Muminuscor_eta",
    "Muminuscor_phi",
    "Muminuscor_pt_quant",
    "Jpsicor_mass_quant",
]


def build_graph(df, dataset):
    # df = df.DefinePerSample("weighttmp", "1.0")
    # weightsum = df.SumAndCount("weighttmp")
    # df = build_graph_base(df, dataset, max_events=int(10e6))

    df, weightsum = build_graph_base(df, dataset)

    # Apply the continuous quantile transform to the nominal pt / mass columns
    # used as axis inputs.
    df, _, _ = narf.histutils.define_quantile_ints(
        df, cols=["Mupluscor_pt"], quantile_hists=quantile_hists_pt
    )
    df, _, _ = narf.histutils.define_quantile_ints(
        df, cols=["Muminuscor_pt"], quantile_hists=quantile_hists_pt
    )
    df, _, _ = narf.histutils.define_quantile_ints(
        df, cols=["Jpsicor_mass"], quantile_hists=quantile_hists_mass
    )

    results = []

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
        df = narf.rdfutils.flexible_define(
            df,
            "mom4_var_plus",
            helper_var,
            [
                "Mupluscor_pt",
                "Mupluscor_eta",
                "Mupluscor_phi",
                "Muplus_charge",
                "globalidxv",
                "Muplus_jacRef",
            ],
        )

        df = narf.rdfutils.flexible_define(
            df,
            "mom4_var_minus",
            helper_var,
            [
                "Muminuscor_pt",
                "Muminuscor_eta",
                "Muminuscor_phi",
                "Muminus_charge",
                "globalidxv",
                "Muminus_jacRef",
            ],
        )

        df = df.Define("mom4var_jpsi", "mom4_var_plus + mom4_var_minus")

        df = df.Define(
            "Muplus_pt_var",
            "return ROOT::VecOps::Map(mom4_var_plus, [](auto const &x){ return x.Pt(); })",
        )
        df = df.Define(
            "Muplus_eta_var",
            "return ROOT::VecOps::Map(mom4_var_plus, [](auto const &x){ return x.Eta(); })",
        )
        df = df.Define(
            "Muplus_phi_var",
            "return ROOT::VecOps::Map(mom4_var_plus, [](auto const &x){ return x.Phi(); })",
        )

        df = df.Define(
            "Muminus_pt_var",
            "return ROOT::VecOps::Map(mom4_var_minus, [](auto const &x){ return x.Pt(); })",
        )
        df = df.Define(
            "Muminus_eta_var",
            "return ROOT::VecOps::Map(mom4_var_minus, [](auto const &x){ return x.Eta(); })",
        )
        df = df.Define(
            "Muminus_phi_var",
            "return ROOT::VecOps::Map(mom4_var_minus, [](auto const &x){ return x.Phi(); })",
        )

        df = df.Define(
            "Jpsi_mass_var",
            "return ROOT::VecOps::Map(mom4var_jpsi, [](auto const &x){ return x.M(); })",
        )

        # Apply the continuous quantile transform element-wise to the shifted
        # pt / mass variation RVecs (MapWrapper broadcasts over containers).
        df, _, _ = narf.histutils.define_quantile_ints(
            df, cols=["Muplus_pt_var"], quantile_hists=quantile_hists_pt
        )
        df, _, _ = narf.histutils.define_quantile_ints(
            df, cols=["Muminus_pt_var"], quantile_hists=quantile_hists_pt
        )
        df, _, _ = narf.histutils.define_quantile_ints(
            df, cols=["Jpsi_mass_var"], quantile_hists=quantile_hists_mass
        )

        df = narf.histutils.shifted_smeared_hist_weight(
            df,
            "Muplus_shift_weight",
            axes=nominal_axes,
            original_cols=nominal_cols_plus,
            shifted_cols=[
                "Muplus_eta_var",
                "Muplus_phi_var",
                "Muplus_pt_var_quant",
                "Jpsi_mass_var_quant",
            ],
            nominal_weight_col="nominal_weight",
        )

        df = df.Define("Muplus_shift_weight_delta", "Muplus_shift_weight - nominal_weight")

        hmuplus_corparms = df.HistoBoost(
            "hmuplus_corparms",
            axes=nominal_axes + [axis_corparms],
            cols=nominal_cols_plus + ["globalidxv", "Muplus_shift_weight_delta"],
            # storage=hist.storage.Double(),
            storage = narf.histutils.SparseStorage(0.02),
        )
        results.append(hmuplus_corparms)

        df = narf.histutils.shifted_smeared_hist_weight(
            df,
            "Muminus_shift_weight",
            axes=nominal_axes,
            original_cols=nominal_cols_minus,
            shifted_cols=[
                "Muminus_eta_var",
                "Muminus_phi_var",
                "Muminus_pt_var_quant",
                "Jpsi_mass_var_quant",
            ],
            nominal_weight_col="nominal_weight",
        )

        df = df.Define("Muminus_shift_weight_delta", "Muminus_shift_weight - nominal_weight")

        hmuminus_corparms = df.HistoBoost(
            "hmuminus_corparms",
            axes=nominal_axes + [axis_corparms],
            cols=nominal_cols_minus + ["globalidxv", "Muminus_shift_weight_delta"],
            # storage=hist.storage.Double(),
            storage = narf.histutils.SparseStorage(0.02),
        )
        results.append(hmuminus_corparms)

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree="tree")

# Also persist the quantile-transform bin centers and volumes so downstream
# analysis can map the [0, 1] quantile axes back to the original pt / mass.
resultdict["quantile_pt_centers"] = pt_centers_hist
resultdict["quantile_pt_volume"] = pt_volume_hist
resultdict["quantile_mass_centers"] = mass_centers_hist
resultdict["quantile_mass_volume"] = mass_volume_hist

fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
write_analysis_output(resultdict, fout, args)
