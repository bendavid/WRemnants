"""Produce a flat snapshot of the ideal sample for normalizing-flow training.

Reads the J/psi ideal ntuple(s), applies the LBL muon calibration to get
post-correction Mu{plus,minus}cor_{pt,eta,phi}, and writes a slim TTree
containing only the per-muon gen/reco columns needed by
``train_muon_response_flow.py``. Everything else (LBL corrections, RDF,
narf, wremnants, ROOT) lives in this one step; the training script then
reads the snapshot via plain uproot in a minimal Python environment.

Branches written (one row per J/psi event):
    Mupluscor_{pt,eta,phi}, Muplusgen_{pt,eta,phi}
    Muminuscor_{pt,eta,phi}, Muminusgen_{pt,eta,phi}

Multi-path inputs are accepted: each ``--input-paths`` entry may be a
single .root file or a top-level directory / xrootd URL to be walked
recursively via
``wremnants.production.datasets.dataset_tools.buildFileList``.
"""

import argparse
import os
import sys
import time
from typing import List

import ROOT

import narf
import wremnants
import wremnants.production.muon_calibration
import wremnants.production.pileup
import wremnants.production.vertex
from wremnants.production.datasets import dataset_tools


DEFAULT_INPUT_PATHS = [
    "root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/"
    "JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_"
    "RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint_idealgeom/"
    "230214_153350/",
    "root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/"
    "JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_"
    "RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint_idealgeom/"
    "230214_153512/",
]


OUTPUT_BRANCHES = [
    "Mupluscor_pt",
    "Mupluscor_eta",
    "Mupluscor_phi",
    "Muminuscor_pt",
    "Muminuscor_eta",
    "Muminuscor_phi",
    "Muplusgen_pt",
    "Muplusgen_eta",
    "Muplusgen_phi",
    "Muminusgen_pt",
    "Muminusgen_eta",
    "Muminusgen_phi",
    "nominal_weight",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-paths",
        nargs="+",
        default=DEFAULT_INPUT_PATHS,
        help="Top-level xrootd/posix paths or explicit .root files. "
        "Each entry is resolved via "
        "wremnants.production.datasets.dataset_tools.buildFileList. "
        "Defaults to the two ideal-geometry J/psi MC samples "
        "used by jpsi_module_corrections.py.",
    )
    p.add_argument(
        "--tree", default="tree", help="TTree name inside the input files"
    )
    p.add_argument(
        "--output",
        default="flow_training_snapshot.root",
        help="Output ROOT file path.",
    )
    p.add_argument(
        "--output-tree",
        default="tree",
        help="TTree name inside the output snapshot.",
    )
    p.add_argument(
        "--max-files",
        type=int,
        default=-1,
        help="Cap on number of input .root files after resolution. "
        "-1 uses all.",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Number of threads for ImplicitMT. 0 leaves ROOT to choose "
        "(typically all available cores).",
    )
    p.add_argument(
        "--pt-min",
        type=float,
        default=2.0,
        help="Reject events where either gen pt is below this (GeV). "
        "Kept generous; the training script can tighten further.",
    )
    p.add_argument(
        "--eta-max",
        type=float,
        default=2.5,
        help="|gen eta| cut, applied to both muons. Reco eta is "
        "kept unrestricted so the training script sees the full "
        "response distribution at any reco eta given the gen "
        "kinematics.",
    )
    p.add_argument(
        "--era",
        default="2016PostVFP",
        help="Era label for the pileup and vertex helpers.",
    )
    p.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        default=True,
        help="Disable the RDF progress bar.",
    )
    return p.parse_args()


def resolve_input_paths(paths: List[str]) -> List[str]:
    out: List[str] = []
    for p in paths:
        if p.lower().endswith(".root"):
            out.append(p)
            continue
        found = dataset_tools.buildFileList(p)
        if not found:
            print(
                f"warning: no .root files found under {p}", file=sys.stderr
            )
        out.extend(found)
    return out


def main():
    args = parse_args()

    print("resolving input paths")
    files = resolve_input_paths(args.input_paths)
    if not files:
        print("error: no input files resolved", file=sys.stderr)
        return 1
    if args.max_files > 0 and len(files) > args.max_files:
        print(f"  capping from {len(files)} to {args.max_files} files")
        files = files[: args.max_files]
    print(f"  {len(files)} file(s) to read")

    if args.threads == 0:
        ROOT.ROOT.EnableImplicitMT()
    elif args.threads > 1:
        ROOT.ROOT.EnableImplicitMT(args.threads)

    print("building RDF and applying LBL corrections")
    t0 = time.time()
    df = ROOT.ROOT.RDataFrame(args.tree, files)
    if args.progress:
        ROOT.ROOT.RDF.Experimental.AddProgressBar(df)

    helper = (
        wremnants.production.muon_calibration.make_muon_calibration_helper_single()
    )
    df = (
        wremnants.production.muon_calibration
        .define_lbl_corrections_jpsi_calibration_ntuples(df, helper)
    )

    pileup_helper = wremnants.production.pileup.make_pileup_helper(era=args.era)
    vertex_helper = wremnants.production.vertex.make_vertex_helper(era=args.era)
    df = df.DefinePerSample("weight", "1.0")
    df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
    df = df.Define("weight_vtx", vertex_helper, ["Jpsigen_z", "Pileup_nTrueInt"])
    df = df.Define("nominal_weight", "weight*weight_pu*weight_vtx")

    df = df.Filter(
        f"std::fabs(Muplusgen_eta) < {args.eta_max} && "
        f"std::fabs(Muminusgen_eta) < {args.eta_max}",
        "eta_acceptance",
    )
    df = df.Filter(
        f"Muplusgen_pt > {args.pt_min} && Muminusgen_pt > {args.pt_min} && "
        f"Mupluscor_pt > 0. && Muminuscor_pt > 0.",
        "pt_positivity_and_gen_minimum",
    )

    out = args.output
    out_dir = os.path.dirname(os.path.abspath(out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"snapshotting to {out}")
    snapshot_options = ROOT.RDF.RSnapshotOptions()
    snapshot_options.fCompressionAlgorithm = (
        ROOT.RCompressionSetting.EAlgorithm.kZSTD
    )
    snapshot_options.fCompressionLevel = 5

    cols_vec = ROOT.std.vector("string")()
    for c in OUTPUT_BRANCHES:
        cols_vec.push_back(c)

    df.Snapshot(args.output_tree, out, cols_vec, snapshot_options)

    dt = time.time() - t0
    # Report the final counts from the filter cutflow.
    print("done.")
    print(f"elapsed: {dt:.1f}s")
    print(f"output file size: {os.path.getsize(out) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
