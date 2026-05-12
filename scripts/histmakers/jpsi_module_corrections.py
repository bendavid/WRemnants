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

parser.add_argument(
    "--computeQuantiles",
    action="store_true",
    help="Compute quantile edges from a pre-pass over the ideal sample and "
    "save them to a file. When omitted the pre-computed quantile file is "
    "loaded instead.",
)
parser.add_argument(
    "--nomCorrections",
    default=None,
    type=str,
    help="Path to a correctionResults ROOT file to apply to jpsi_nom "
    "instead of the passthrough (no-correction) default. This allows "
    "closure tests where fitted corrections are applied back to the "
    "nominal sample.",
)
parser.add_argument(
    "--injectA",
    default=0.0,
    type=float,
    help="Inject a known A-type curvature shift (dimensionless) into the "
    "jpsi_nom sample, applied uniformly across all eta bins before the "
    "kinematic filters. The parameterized fit should recover this value "
    "in its A output. Convention: k_new = k*(1+A). Example: 1e-4 for a "
    "100 ppm curvature bias.",
)
parser.add_argument(
    "--injectE",
    default=0.0,
    type=float,
    help="Inject a known e-type curvature shift (units of GeV) into the "
    "jpsi_nom sample, applied uniformly across all eta bins. "
    "Convention: k_new += -e*k^2. Example: 1e-3.",
)
parser.add_argument(
    "--injectM",
    default=0.0,
    type=float,
    help="Inject a known M-type charge-dependent curvature shift (units "
    "of GeV^-1) into the jpsi_nom sample, applied uniformly across all "
    "eta bins. Convention: k_new += charge*M. Example: 1e-5.",
)
parser.add_argument(
    "--skipCorparms",
    action="store_true",
    default=False,
    help="Skip building the per-corparm SparseHist variations "
    "(hmuplus_corparms / hmuminus_corparms) on the ideal sample. Useful "
    "for parameterized-only fits (A/e/M scale shifts), which are faster "
    "to build and require less memory.",
)
parser.add_argument(
    "--injectIdeal",
    action="store_true",
    default=False,
    help="Also apply the --injectA/E/M shift to the ideal sample "
    "(applied as an overlay on top of the LBL-corrected Mupluscor_pt "
    "columns; Jpsi kinematics recomputed from shifted mom4s). "
    "Diagnostic only: lets us compare (injected ideal) against "
    "(nominal ideal + stored hmu*_scale at A=+1sigma) to check the "
    "injection path and the helper_scale_shift path produce consistent "
    "4-vectors.",
)
parser.add_argument(
    "--rawAxes",
    action="store_true",
    default=False,
    help="Use raw Regular(pt) and Regular(mass) axes instead of the "
    "chained quantile-transformed (pt_quant, mass_quant) axes. "
    "Diagnostic only: lets us check whether the quantile-transform "
    "layer is the source of the closure-test bias.",
)
parser.add_argument(
    "--skipMaxPtCut",
    action="store_true",
    default=False,
    help="Skip the max(Mupluscor_pt, Muminuscor_pt) > 13.2 GeV filter. "
    "Diagnostic only: isolates filter-boundary effects in closure tests.",
)
parser.add_argument(
    "--unchainedQuantiles",
    action="store_true",
    default=False,
    help="Build two INDEPENDENT quantile transforms (pt_quant conditioned "
    "on eta only; mass_quant conditioned on eta only), instead of the "
    "default chained version where mass_quant is conditioned on "
    "(eta, pt_quant). Diagnostic only.",
)
parser.add_argument(
    "--explicitShift",
    action="store_true",
    default=False,
    help="Fill hmu{plus,minus}_scale by explicitly filling each "
    "(scale_eta, param) slot with the shifted (pt, mass) coordinates "
    "and the nominal weight, instead of via the dipole-approximation "
    "shift_weight tensor. Produces absolute shifted histograms in both "
    "paths (drop the delta subtraction). Diagnostic for whether the "
    "dipole approximation in shifted_smeared_hist_weight is biasing "
    "the closure test. Parameterized variations only; corparms path "
    "continues to use the shift_weight approach.",
)

args = parser.parse_args()

inject_any = (
    args.injectA != 0.0 or args.injectE != 0.0 or args.injectM != 0.0
)
if inject_any and args.nomCorrections is not None:
    raise RuntimeError(
        "--injectA/--injectE/--injectM cannot be combined with "
        "--nomCorrections (the two would stack in non-obvious ways). "
        "Use them exclusively."
    )


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
pathsideal.append("root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint_idealgeom/230214_153512/")
fideal = paths_to_filenames(pathsideal)


pathsnom = []
pathsnom.append(
    "root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInf_quality_novtx_noconstraint/230214_151859/"
)
pathsnom.append("root://eoscms.cern.ch//store/group/phys_smp/ec/bendavid/muoncal/JPsiToMuMu_Pt8toInf-pythia8/MuonGunUL2016_v722_RecJpsiPythiaPhotosPt8toInfExt1_quality_novtx_noconstraint/230214_152107/")
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

helper_nom = None
if args.nomCorrections is not None:
    print(f"Loading nominal corrections from {args.nomCorrections}")
    helper_nom = wremnants.production.muon_calibration.make_muon_calibration_helper_single(
        filename=args.nomCorrections
    )

axis_eta_scale = hist.axis.Regular(
    netabins, -2.4, 2.4, name="scale_eta", underflow=False, overflow=False
)
helper_scale_shift = wremnants.production.muon_calibration.make_parameterized_scale_shift_helper(
    eta_axis=axis_eta_scale,
)

if args.explicitShift:
    ROOT.gInterpreter.Declare(
        """
        #ifndef WREM_JPSI_EXPLICIT_SHIFT_HPP
        #define WREM_JPSI_EXPLICIT_SHIFT_HPP
        #include <ROOT/RVec.hxx>
        #include <array>
        #include <cstddef>
        #include <string>

        namespace wrem_jpsi_explicit_shift {
          constexpr std::size_t NETA = 24;
          constexpr std::size_t NPARAM = 3;

          struct FillCols {
            ROOT::VecOps::RVec<double> pt;
            ROOT::VecOps::RVec<double> mass;
            ROOT::VecOps::RVec<double> sc_eta;
            ROOT::VecOps::RVec<std::string> param;
          };

          template <typename PtShiftT, typename MassShiftT>
          FillCols make_fill_cols(
              const PtShiftT &pt_shift,
              const MassShiftT &mass_shift,
              double sc_eta_low,
              double sc_eta_high) {
            constexpr std::size_t N = NETA * NPARAM;
            const std::array<std::string, NPARAM> param_names = {"A", "e", "M"};
            const double width =
                (sc_eta_high - sc_eta_low) / static_cast<double>(NETA);

            FillCols cols;
            cols.pt.resize(N);
            cols.mass.resize(N);
            cols.sc_eta.resize(N);
            cols.param.resize(N);

            for (std::size_t i = 0; i < NETA; ++i) {
              const double center =
                  sc_eta_low + (static_cast<double>(i) + 0.5) * width;
              for (std::size_t j = 0; j < NPARAM; ++j) {
                const std::size_t idx = i * NPARAM + j;
                cols.pt[idx] = pt_shift(i, j);
                cols.mass[idx] = mass_shift(i, j);
                cols.sc_eta[idx] = center;
                cols.param[idx] = param_names[j];
              }
            }
            return cols;
          }
        }
        #endif
        """
    )

print("done defining helpers")


def _apply_injected_scale_shift_overlay(df, A, E, M):
    """Apply injection on top of existing Mupluscor_pt/Muminuscor_pt
    columns (i.e. after LBL corrections have already run). Used by
    --injectIdeal for the construction-consistency diagnostic.

    Redefines the corrected pt columns and the downstream
    Mupluscor_mom4 / Muminuscor_mom4 / Jpsi* columns so the
    histogram-fill path sees the shifted values.
    """
    shift_expr = (
        "[](double pt, int charge) {"
        f"  const double A = {A:.17g};"
        f"  const double E = {E:.17g};"
        f"  const double M = {M:.17g};"
        "  const double k = 1.0 / pt;"
        "  const double k_new = k * (1.0 + A) - E * k * k"
        "    + static_cast<double>(charge) * M;"
        "  return 1.0 / k_new;"
        "}"
    )
    df = df.Redefine(
        "Mupluscor_pt",
        f"return ({shift_expr})(Mupluscor_pt, Muplus_charge);",
    )
    df = df.Redefine(
        "Muminuscor_pt",
        f"return ({shift_expr})(Muminuscor_pt, Muminus_charge);",
    )
    df = df.Redefine(
        "Mupluscor_mom4",
        "ROOT::Math::PtEtaPhiMVector(Mupluscor_pt, Mupluscor_eta, "
        "Mupluscor_phi, wrem::muon_mass)",
    )
    df = df.Redefine(
        "Muminuscor_mom4",
        "ROOT::Math::PtEtaPhiMVector(Muminuscor_pt, Muminuscor_eta, "
        "Muminuscor_phi, wrem::muon_mass)",
    )
    df = df.Redefine(
        "Jpsicor_mom4",
        "ROOT::Math::PxPyPzEVector(Mupluscor_mom4) "
        "+ ROOT::Math::PxPyPzEVector(Muminuscor_mom4)",
    )
    df = df.Redefine("Jpsicor_pt", "Jpsicor_mom4.Pt()")
    df = df.Redefine("Jpsicor_eta", "Jpsicor_mom4.Eta()")
    df = df.Redefine("Jpsicor_phi", "Jpsicor_mom4.Phi()")
    df = df.Redefine("Jpsicor_mass", "Jpsicor_mom4.M()")
    return df


def _apply_injected_scale_shift_nom(df, A, E, M):
    """Passthrough-with-injection path for jpsi_nom.

    Applies a fixed A/E/M curvature shift (uniform across eta bins) to
    the raw muon pt's and rebuilds the corrected muon and Jpsi 4-vectors
    from the shifted values. The shift convention matches
    ``ParameterizedScaleShiftHelper``:
        k_new = k*(1 + A) - E*k^2 + charge*M,   pt_new = 1 / k_new
    where ``k = 1/pt``. Applied before the kinematic filters so cuts see
    the shifted pt's. The parameterized fit run against these histograms
    should recover A/E/M as the per-eta-bin scale coefficients.
    """
    shift_expr = (
        "[](double pt, int charge) {"
        f"  const double A = {A:.17g};"
        f"  const double E = {E:.17g};"
        f"  const double M = {M:.17g};"
        "  const double k = 1.0 / pt;"
        "  const double k_new = k * (1.0 + A) - E * k * k"
        "    + static_cast<double>(charge) * M;"
        "  return 1.0 / k_new;"
        "}"
    )

    df = df.DefinePerSample("Muplus_charge", "1")
    df = df.DefinePerSample("Muminus_charge", "-1")

    # Explicit ``return`` is required here: RDF otherwise auto-prepends
    # ``return`` to the Define expression, but it skips that when it
    # detects a ``return`` token anywhere in the string (including
    # inside the lambda body), which would leave the IIFE result
    # unreturned and type the column as void.
    df = df.Define(
        "Mupluscor_pt",
        f"return ({shift_expr})(Muplus_pt, Muplus_charge);",
    )
    df = df.Alias("Mupluscor_eta", "Muplus_eta")
    df = df.Alias("Mupluscor_phi", "Muplus_phi")

    df = df.Define(
        "Muminuscor_pt",
        f"return ({shift_expr})(Muminus_pt, Muminus_charge);",
    )
    df = df.Alias("Muminuscor_eta", "Muminus_eta")
    df = df.Alias("Muminuscor_phi", "Muminus_phi")

    df = df.Define(
        "Mupluscor_mom4",
        "ROOT::Math::PtEtaPhiMVector(Mupluscor_pt, Mupluscor_eta, "
        "Mupluscor_phi, wrem::muon_mass)",
    )
    df = df.Define(
        "Muminuscor_mom4",
        "ROOT::Math::PtEtaPhiMVector(Muminuscor_pt, Muminuscor_eta, "
        "Muminuscor_phi, wrem::muon_mass)",
    )
    df = df.Define(
        "Jpsicor_mom4",
        "ROOT::Math::PxPyPzEVector(Mupluscor_mom4) "
        "+ ROOT::Math::PxPyPzEVector(Muminuscor_mom4)",
    )
    df = df.Define("Jpsicor_pt", "Jpsicor_mom4.Pt()")
    df = df.Define("Jpsicor_eta", "Jpsicor_mom4.Eta()")
    df = df.Define("Jpsicor_phi", "Jpsicor_mom4.Phi()")
    df = df.Define("Jpsicor_mass", "Jpsicor_mom4.M()")

    return df


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
        df = df.Filter("event%2 == 0")
    else:
        df = df.Filter("event%2 == 1")

    if is_ideal:
        print("ideal: applying corrections")
        df = wremnants.production.muon_calibration.define_lbl_corrections_jpsi_calibration_ntuples(
            df, helper
        )
        if inject_any and args.injectIdeal:
            print(
                f"ideal: also overlaying injection A={args.injectA:.3e}, "
                f"e={args.injectE:.3e}, M={args.injectM:.3e} (diagnostic)"
            )
            df = _apply_injected_scale_shift_overlay(
                df, args.injectA, args.injectE, args.injectM
            )
    elif helper_nom is not None:
        print("nominal: applying corrections from --nomCorrections")
        df = wremnants.production.muon_calibration.define_lbl_corrections_jpsi_calibration_ntuples(
            df, helper_nom
        )
    elif inject_any:
        print(
            f"nominal: injecting scale shift A={args.injectA:.3e}, "
            f"e={args.injectE:.3e}, M={args.injectM:.3e}"
        )
        df = _apply_injected_scale_shift_nom(
            df, args.injectA, args.injectE, args.injectM
        )
    else:
        print("nominal: no corrections")
        df = wremnants.production.muon_calibration.define_passthrough_corrections_jpsi_calibration_ntuples(
            df
        )

    df = df.Filter("std::fabs(Mupluscor_eta) < 2.4 && std::fabs(Muminuscor_eta) < 2.4")

    # df = df.Filter(f"Mupluscor_pt > {muptmin}")
    # df = df.Filter(f"Muminuscor_pt > {muptmin}")
    if not args.skipMaxPtCut:
        df = df.Filter("max(Mupluscor_pt, Muminuscor_pt) > 13.2")
    df = df.Filter(f"min(Mupluscor_pt, Muminuscor_pt) > {muptmin}")
    df = df.Filter("Jpsicor_pt > 8.2")

    df = df.Filter(f"Jpsicor_mass > {mmin} && Jpsicor_mass < {mmax}")

    return df, weightsum


import h5py
import wums.ioutils

quantile_file = "jpsi_module_corrections_quantiles.hdf5"

if args.rawAxes:
    # Raw linear axes for the Test A diagnostic. Use a Variable pt axis
    # with an explicit bin boundary at 13.2 GeV so the max(pt+,pt-)>13.2
    # filter aligns with a bin edge — cleanly separating acceptance-
    # boundary events into the lowest pt bin.
    pt_edges = [muptmin, 13.2]
    hi = muptmin + 50.0
    n_above = npt_quant - 1  # remaining bins above 13.2
    step = (hi - 13.2) / n_above
    pt_edges += [13.2 + (i + 1) * step for i in range(n_above)]
    axis_pt_raw = hist.axis.Variable(
        pt_edges, underflow=False, overflow=False, name="pt",
    )
    axis_mass_raw = hist.axis.Regular(
        2*nmass_quant, mmin, mmax,
        underflow=False, overflow=False, name="mass",
    )
    print(f"rawAxes pt edges: {pt_edges}")
    # Dummy empty quantile_hists / helpers so the rest of the script
    # doesn't fail when we skip the quantile computation.
    quantile_hists = None
    centers_hist = None
    volume_hist = None
elif args.computeQuantiles:
    # Pre-pass over a subset of the ideal sample to build continuous quantile
    # helpers for pt and mass, chained conditional on eta. A finely-binned
    # histogram (eta, pt_fine, mass_fine) is filled in a single multi-threaded
    # event loop, then the chained quantile edges are extracted offline from
    # cumulative sums — no sorting or data materialization needed.
    dfideal_quant = ROOT.ROOT.RDataFrame("tree", fideal)
    ROOT.ROOT.RDF.Experimental.AddProgressBar(dfideal_quant)
    dfideal_quant, _ = build_graph_base(dfideal_quant, dataset_ideal)

    # Fine axes covering the filtered data range.
    axis_pt_fine = hist.axis.Regular(
        500, muptmin, muptmin + 50., name="pt_fine"
    )
    axis_mass_fine = hist.axis.Regular(
        500, mmin, mmax, name="mass_fine"
    )

    fine_hist = dfideal_quant.HistoBoost(
        "fine_quant",
        [axis_eta, axis_pt_fine, axis_mass_fine],
        ["Mupluscor_eta", "Mupluscor_pt", "Jpsicor_mass", "nominal_weight"],
        storage=hist.storage.Double(),
    )
    fine_hist = fine_hist.GetValue()

    if args.unchainedQuantiles:
        # Build two independent quantile helpers, each conditioned on eta
        # only (not chained via pt_quant_int).
        fine_hist_pt = fine_hist.project("eta", "pt_fine")
        fine_hist_mass = fine_hist.project("eta", "mass_fine")
        qhists_pt, centers_pt, volume_pt = (
            narf.histutils.build_quantile_hists_from_fine(
                fine_hist_pt, condaxes=[axis_eta],
                quantaxes=[axis_pt_quant], continuous=True,
            )
        )
        qhists_mass, centers_mass, volume_mass = (
            narf.histutils.build_quantile_hists_from_fine(
                fine_hist_mass, condaxes=[axis_eta],
                quantaxes=[axis_mass_quant], continuous=True,
            )
        )
        # Stored as two independent 1-quantile helpers.
        quantile_hists = [qhists_pt[0], qhists_mass[0]]
        # Use pt volume/centers for the volume hist (diagnostic only;
        # not used downstream for this unchained-quantiles test).
        centers_hist = centers_pt
        volume_hist = volume_pt
    else:
        quantile_hists, centers_hist, volume_hist = (
            narf.histutils.build_quantile_hists_from_fine(
                fine_hist,
                condaxes=[axis_eta],
                quantaxes=[axis_pt_quant, axis_mass_quant],
                continuous=True,
            )
        )

    quantile_data = {
        "quantile_hists": wums.ioutils.H5PickleProxy(quantile_hists),
        "centers_hist": wums.ioutils.H5PickleProxy(centers_hist),
        "volume_hist": wums.ioutils.H5PickleProxy(volume_hist),
    }

    with h5py.File(quantile_file, "w") as fq:
        wums.ioutils.pickle_dump_h5py("quantile_data", quantile_data, fq)
    print(f"Quantile histograms saved to {quantile_file}")
else:
    with h5py.File(quantile_file, "r") as fq:
        quantile_data = wums.ioutils.pickle_load_h5py(fq["quantile_data"])

        quantile_hists = quantile_data["quantile_hists"].get()
        centers_hist = quantile_data["centers_hist"].get()
        volume_hist = quantile_data["volume_hist"].get()

if args.rawAxes:
    nominal_axes = [axis_eta, axis_phi, axis_pt_raw, axis_mass_raw]
    nominal_cols_plus = [
        "Mupluscor_eta",
        "Mupluscor_phi",
        "Mupluscor_pt",
        "Jpsicor_mass",
    ]
    nominal_cols_minus = [
        "Muminuscor_eta",
        "Muminuscor_phi",
        "Muminuscor_pt",
        "Jpsicor_mass",
    ]
else:
    nominal_axes = [axis_eta, axis_phi, axis_pt_quant, axis_mass_quant]
    nominal_cols_plus = [
        "Mupluscor_eta",
        "Mupluscor_phi",
        "Mupluscor_pt_plus_quant",
        "Jpsicor_mass_plus_quant",
    ]
    nominal_cols_minus = [
        "Muminuscor_eta",
        "Muminuscor_phi",
        "Muminuscor_pt_minus_quant",
        "Jpsicor_mass_minus_quant",
    ]


def build_graph(df, dataset):
    # df = df.DefinePerSample("weighttmp", "1.0")
    # weightsum = df.SumAndCount("weighttmp")
    # df = build_graph_base(df, dataset, max_events=int(10e6))

    df, weightsum = build_graph_base(df, dataset)

    if not args.rawAxes:
        # Apply the continuous quantile transform to the nominal
        # pt / mass columns. With default (chained) quantiles, a single
        # call threads pt_quant through as conditioning for mass_quant.
        # With --unchainedQuantiles, each quantile variable is
        # conditioned only on eta — call define_quantiles once per
        # variable so no chaining happens.
        if args.unchainedQuantiles:
            df, _, _ = narf.histutils.define_quantiles(
                df, cols=["Mupluscor_eta", "Mupluscor_pt"],
                quantile_hists=[quantile_hists[0]], label="plus",
            )
            df, _, _ = narf.histutils.define_quantiles(
                df, cols=["Mupluscor_eta", "Jpsicor_mass"],
                quantile_hists=[quantile_hists[1]], label="plus",
            )
            df, _, _ = narf.histutils.define_quantiles(
                df, cols=["Muminuscor_eta", "Muminuscor_pt"],
                quantile_hists=[quantile_hists[0]], label="minus",
            )
            df, _, _ = narf.histutils.define_quantiles(
                df, cols=["Muminuscor_eta", "Jpsicor_mass"],
                quantile_hists=[quantile_hists[1]], label="minus",
            )
        else:
            df, _, _ = narf.histutils.define_quantiles(
                df,
                cols=["Mupluscor_eta", "Mupluscor_pt", "Jpsicor_mass"],
                quantile_hists=quantile_hists,
                label="plus",
            )
            df, _, _ = narf.histutils.define_quantiles(
                df,
                cols=["Muminuscor_eta", "Muminuscor_pt", "Jpsicor_mass"],
                quantile_hists=quantile_hists,
                label="minus",
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
        if not args.skipCorparms:
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

            # Apply the chained continuous quantile transform element-wise to
            # the shifted pt / mass variation RVecs (MapWrapper broadcasts over
            # containers).
            df, _, _ = narf.histutils.define_quantiles(
                df,
                cols=["Muplus_eta_var", "Muplus_pt_var", "Jpsi_mass_var"],
                quantile_hists=quantile_hists,
                label="plus",
            )
            df, _, _ = narf.histutils.define_quantiles(
                df,
                cols=["Muminus_eta_var", "Muminus_pt_var", "Jpsi_mass_var"],
                quantile_hists=quantile_hists,
                label="minus",
            )

            df = narf.histutils.shifted_smeared_hist_weight(
                df,
                "Muplus_shift_weight",
                axes=nominal_axes,
                original_cols=nominal_cols_plus,
                shifted_cols=[
                    "Muplus_eta_var",
                    "Muplus_phi_var",
                    "Muplus_pt_var_plus_quant",
                    "Jpsi_mass_var_plus_quant",
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
                metadata = {"dparm": helper_var.dparm()},
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
                    "Muminus_pt_var_minus_quant",
                    "Jpsi_mass_var_minus_quant",
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
                metadata = {"dparm": helper_var.dparm()},
            )
            results.append(hmuminus_corparms)

        # --- A, e, M scale shift variations ---

        # Shifted 4-vector tensors (NEtaBins, 3) for each muon.
        df = narf.rdfutils.flexible_define(
            df,
            "mom4_scale_shift_plus",
            helper_scale_shift,
            ["Mupluscor_pt", "Mupluscor_eta", "Mupluscor_phi", "Muplus_charge"],
        )
        df = narf.rdfutils.flexible_define(
            df,
            "mom4_scale_shift_minus",
            helper_scale_shift,
            ["Muminuscor_pt", "Muminuscor_eta", "Muminuscor_phi", "Muminus_charge"],
        )

        # Shifted J/psi mass tensor: add both muons' shifted 4-vector
        # tensors element-wise.  For each (eta_bin, param) variation, the
        # helper fills non-matching eta bins with the nominal 4-vector, so
        # the sum correctly shifts both muons when both are in the varied
        # eta bin, or just one when only one is.
        df = df.Define(
            "Jpsi_mass_scale_shift",
            "return narf::eval_if_tensor("
            "  (mom4_scale_shift_plus + mom4_scale_shift_minus)"
            "  .unaryExpr([](const ROOT::Math::PxPyPzEVector &v)"
            "  { return v.M(); }))",
        )

        # Shifted pt tensors.
        df = df.Define(
            "Muplus_pt_scale_shift",
            "return narf::eval_if_tensor(mom4_scale_shift_plus.unaryExpr("
            "  [](const ROOT::Math::PxPyPzEVector &v){ return v.Pt(); }))",
        )
        df = df.Define(
            "Muminus_pt_scale_shift",
            "return narf::eval_if_tensor(mom4_scale_shift_minus.unaryExpr("
            "  [](const ROOT::Math::PxPyPzEVector &v){ return v.Pt(); }))",
        )

        if not args.rawAxes:
            # Quantile-transform the shifted pt and mass tensors (the
            # TensorMapWrapper inside the quantile helpers broadcasts the
            # scalar conditioning variables over tensor elements).
            if args.unchainedQuantiles:
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Mupluscor_eta", "Muplus_pt_scale_shift"],
                    quantile_hists=[quantile_hists[0]],
                    label="scale_plus",
                )
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Mupluscor_eta", "Jpsi_mass_scale_shift"],
                    quantile_hists=[quantile_hists[1]],
                    label="scale_plus",
                )
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Muminuscor_eta", "Muminus_pt_scale_shift"],
                    quantile_hists=[quantile_hists[0]],
                    label="scale_minus",
                )
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Muminuscor_eta", "Jpsi_mass_scale_shift"],
                    quantile_hists=[quantile_hists[1]],
                    label="scale_minus",
                )
            else:
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Mupluscor_eta", "Muplus_pt_scale_shift", "Jpsi_mass_scale_shift"],
                    quantile_hists=quantile_hists,
                    label="scale_plus",
                )
                df, _, _ = narf.histutils.define_quantiles(
                    df,
                    cols=["Muminuscor_eta", "Muminus_pt_scale_shift", "Jpsi_mass_scale_shift"],
                    quantile_hists=quantile_hists,
                    label="scale_minus",
                )

        if args.rawAxes:
            scale_shifted_cols_plus = [
                "Mupluscor_eta",
                "Mupluscor_phi",
                "Muplus_pt_scale_shift",
                "Jpsi_mass_scale_shift",
            ]
            scale_shifted_cols_minus = [
                "Muminuscor_eta",
                "Muminuscor_phi",
                "Muminus_pt_scale_shift",
                "Jpsi_mass_scale_shift",
            ]
        else:
            scale_shifted_cols_plus = [
                "Mupluscor_eta",
                "Mupluscor_phi",
                "Muplus_pt_scale_shift_scale_plus_quant",
                "Jpsi_mass_scale_shift_scale_plus_quant",
            ]
            scale_shifted_cols_minus = [
                "Muminuscor_eta",
                "Muminuscor_phi",
                "Muminus_pt_scale_shift_scale_minus_quant",
                "Jpsi_mass_scale_shift_scale_minus_quant",
            ]

        scale_shift_axes = helper_scale_shift.tensor_axes
        scale_eta_ax = scale_shift_axes[0]
        sc_eta_low = float(scale_eta_ax.edges[0])
        sc_eta_high = float(scale_eta_ax.edges[-1])

        if args.explicitShift:
            # Scalar-fill each (scale_eta, param) slot with the shifted pt
            # and mass (and nominal eta/phi/weight). The histogram stores
            # the *absolute* shifted template; downstream fit must set
            # as_difference=False.
            for label, shift_cols in [
                ("plus", scale_shifted_cols_plus),
                ("minus", scale_shifted_cols_minus),
            ]:
                cols_col = f"Mu{label}_explicit_fill_cols"
                df = df.Define(
                    cols_col,
                    "wrem_jpsi_explicit_shift::make_fill_cols("
                    f"{shift_cols[2]}, {shift_cols[3]}, "
                    f"{sc_eta_low}, {sc_eta_high})",
                )
                for field in ("pt", "mass", "sc_eta", "param"):
                    df = df.Define(
                        f"Mu{label}_explicit_{field}",
                        f"{cols_col}.{field}",
                    )

            hmuplus_scale = df.HistoBoost(
                "hmuplus_scale",
                axes=nominal_axes + list(scale_shift_axes),
                cols=[
                    nominal_cols_plus[0],
                    nominal_cols_plus[1],
                    "Muplus_explicit_pt",
                    "Muplus_explicit_mass",
                    "Muplus_explicit_sc_eta",
                    "Muplus_explicit_param",
                    "nominal_weight",
                ],
            )
            results.append(hmuplus_scale)

            hmuminus_scale = df.HistoBoost(
                "hmuminus_scale",
                axes=nominal_axes + list(scale_shift_axes),
                cols=[
                    nominal_cols_minus[0],
                    nominal_cols_minus[1],
                    "Muminus_explicit_pt",
                    "Muminus_explicit_mass",
                    "Muminus_explicit_sc_eta",
                    "Muminus_explicit_param",
                    "nominal_weight",
                ],
            )
            results.append(hmuminus_scale)
        else:
            # Shift weights (tensor-valued): the HistShiftHelper computes
            # element-wise weight corrections for each (eta_bin, param)
            # variation. The histogram stores the *absolute* shifted
            # template (nominal_weight times the correction tensor);
            # downstream fit must set as_difference=False.
            df = narf.histutils.shifted_smeared_hist_weight(
                df,
                "Muplus_scale_shift_weight",
                axes=nominal_axes,
                original_cols=nominal_cols_plus,
                shifted_cols=scale_shifted_cols_plus,
                nominal_weight_col="nominal_weight",
            )

            hmuplus_scale = df.HistoBoost(
                "hmuplus_scale",
                axes=nominal_axes,
                cols=nominal_cols_plus + ["Muplus_scale_shift_weight"],
                tensor_axes=scale_shift_axes,
            )
            results.append(hmuplus_scale)

            df = narf.histutils.shifted_smeared_hist_weight(
                df,
                "Muminus_scale_shift_weight",
                axes=nominal_axes,
                original_cols=nominal_cols_minus,
                shifted_cols=scale_shifted_cols_minus,
                nominal_weight_col="nominal_weight",
            )

            hmuminus_scale = df.HistoBoost(
                "hmuminus_scale",
                axes=nominal_axes,
                cols=nominal_cols_minus + ["Muminus_scale_shift_weight"],
                tensor_axes=scale_shift_axes,
            )
            results.append(hmuminus_scale)

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree="tree")

# Also persist the quantile-transform bin centers and volumes so downstream
# analysis can map the [0, 1] quantile axes back to the original pt / mass.
# Skipped under --rawAxes (no quantile transform was done).
if not args.rawAxes:
    resultdict["quantile_info"] = {
        "centers": wums.ioutils.H5PickleProxy(centers_hist),
        "volume": wums.ioutils.H5PickleProxy(volume_hist),
    }

fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
write_analysis_output(resultdict, fout, args)
