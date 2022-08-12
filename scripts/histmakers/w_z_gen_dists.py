import argparse
import pickle
import gzip
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
initargs,_ = parser.parse_known_args()

ROOT.gInterpreter.ProcessLine(".O3")
if not initargs.nThreads:
    ROOT.ROOT.EnableImplicitMT()
elif initargs.nThreads != 1:
    ROOT.ROOT.EnableImplicitMT(initargs.nThreads)
import narf
import wremnants
from wremnants import theory_tools,syst_tools
import hist
import lz4.frame
import logging
import math

logging.basicConfig(level=logging.INFO)

parser.add_argument("--pdfs", type=str, nargs="*", default=["nnpdf31"], choices=theory_tools.pdfMapExtended.keys(), help="PDF sets to produce error hists for")
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=["Wplus", "Wminus", "Zmumu", "Ztautau"])
parser.add_argument("--skipAngularCoeffs", action='store_true', help="Skip the conversion of helicity moments to angular coeff fractions")
parser.add_argument("--singleLeptonHists", action='store_true', help="Also store single lepton kinematics")
parser.add_argument("--v8", action='store_true', help="Use NanoAODv8. Default is v9")
args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])

datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, mode="gen")

print('Use v8?', args.v8)
if args.v8:
    datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, nanoVersion = "v8")


axis_massWgen = hist.axis.Variable([0., math.inf], name="massVgen")

axis_massZgen = hist.axis.Regular(12, 60., 120., name="massVgen")

axis_absYVgen = hist.axis.Variable(
    [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 10], name = "absYVgen"
)
axis_ptVgen = hist.axis.Variable(
#    [0, 2, 3, 4, 4.75, 5.5, 6.5, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100], name = "ptVgen"
    range(0,121), name = "ptVgen"
)


axis_chargeWgen = hist.axis.Regular(
    2, -2, 2, name="chargeVgen", underflow=False, overflow=False
)

axis_chargeZgen = hist.axis.Integer(
    0, 1, name="chargeVgen", underflow=False, overflow=False
)

axis_l_eta_gen = hist.axis.Regular(48, -2.4, 2.4, name = "prefsr_lepton_eta_gen")
axis_l_pt_gen = hist.axis.Regular(29, 26., 55., name = "prefsr_lepton_pt_gen")
axes_l_gen = [axis_l_eta_gen, axis_l_pt_gen]

wprocs_bugged = [
    "WplusmunuPostVFP", 
    "WminusmunuPostVFP",
    "WminustaunuPostVFP",
    "WplustaunuPostVFP"
]
wprocs_bugfix = [
    "Wplusmunu_bugfix", 
    "Wminusmunu_bugfix_newprod",
    "Wplusmunu_bugfix",
    "Wplusmunu_bugfix_reweight_h2",
]
wprocs = [*wprocs_bugged, *wprocs_bugfix]
wprocs_bugged_to_check = [ # mu only for making comparison plots between bugged and bugfix samples
    "WplusmunuPostVFP", 
    "WminusmunuPostVFP",
]
zprocs_bugged = ["ZmumuPostVFP", "ZtautauPostVFP"]
zprocs_bugfix = ["ZmumuPostVFP_bugfix", "ZmumuPostVFP_bugfix_slc7"]
zprocs = [*zprocs_bugged, *zprocs_bugfix]
zprocs_bugged_to_check = ["ZmumuPostVFP"]

def build_graph(df, dataset):
    print("build graph")
    print(dataset.name)
    results = []
    
    if dataset.is_data:
        raise RuntimeError("Running GEN analysis over data is not supported")

    df = df.Define("nominal_pdf_cen", theory_tools.pdf_central_weight(dataset.name, args.pdfs[0]))
    weight_expr = "std::copysign(1.0, genWeight)"
    weight_expr = f"{weight_expr}*nominal_pdf_cen"

    if "reweight_h2" in dataset.name:
        weight_expr = f"{weight_expr}*H2BugFixWeight[0]"

    df = df.Define("nominal_weight", weight_expr)
    weightsum = df.SumAndCount("nominal_weight")

    df = theory_tools.define_scale_tensor(df)

    if dataset.name in zprocs:
        nominal_axes = [axis_massZgen, axis_absYVgen, axis_ptVgen, axis_chargeZgen]
    else:
        nominal_axes = [axis_massWgen, axis_absYVgen, axis_ptVgen, axis_chargeWgen]

    nominal_cols = ["massVgen", "absYVgen", "ptVgen", "chargeVgen"]
    df = wremnants.define_prefsr_vars(df)

    if args.singleLeptonHists:
        if dataset.name == 'WplusmunuPostVFP':
            df = df.Define('ptPrefsrMuon', 'genlanti.pt()')
            df = df.Define('etaPrefsrMuon', 'genlanti.eta()')
        elif dataset.name == 'WminusmunuPostVFP' or 'ZmumuPostVFP' in dataset.name:
            df = df.Define('ptPrefsrMuon', 'genl.pt()')
            df = df.Define('etaPrefsrMuon', 'genl.eta()')
        if dataset.name in ['WplusmunuPostVFP', 'WminusmunuPostVFP'] or 'ZmumuPostVFP' in dataset.name:
            nominal_cols = [*nominal_cols, 'etaPrefsrMuon', 'ptPrefsrMuon']
            nominal_axes = [*nominal_axes, axis_l_eta_gen, axis_l_pt_gen]

    nominal_gen = df.HistoBoost("nominal_gen", nominal_axes, [*nominal_cols, "nominal_weight"])
    results.append(nominal_gen)

    results.append(theory_tools.make_scale_hist(df, nominal_axes, nominal_cols))

    for pdf in args.pdfs:
        results.extend(theory_tools.define_and_make_pdf_hists(df, nominal_axes, nominal_cols, dataset.name, pdfset=pdf))

    if dataset.name != "Zmumu_bugfix":
        isW = dataset.name in wprocs

        df, masswhist = syst_tools.define_mass_weights(df, isW, nominal_axes, nominal_cols)
        if masswhist:
            results.append(masswhist)

    df = df.Define("helicity_moments_scale_tensor", "wrem::makeHelicityMomentScaleTensor(csSineCosThetaPhi, scaleWeights_tensor, nominal_weight)")
    helicity_moments_scale = df.HistoBoost("helicity_moments_scale", nominal_axes, [*nominal_cols, "helicity_moments_scale_tensor"], tensor_axes = [wremnants.axis_helicity, *wremnants.scale_tensor_axes])
    results.append(helicity_moments_scale)

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

fname = "w_z_gen_dists.pkl.lz4"

print("writing output")
with lz4.frame.open(fname, "wb") as f:
    pickle.dump(resultdict, f, protocol = pickle.HIGHEST_PROTOCOL)

print("computing angular coefficients")

z_moments_bugged = None
z_moments_bugfix = None
w_moments_bugged = None
w_moments_bugfix = None

if args.skipAngularCoeffs:
    exit(0)

# TODO: This should be cleaned up
for key, val in resultdict.items():
    moments = val["output"]["helicity_moments_scale"]
    if key in zprocs_bugged_to_check:
        if z_moments_bugged is None:
            z_moments_bugged = moments
        else:
            # gen level kinematics, stack tau to mu channel to increase stats
            z_moments_bugged += moments
    elif key in wprocs_bugged_to_check:
        if w_moments_bugged is None:
            w_moments_bugged = moments
        else:
            w_moments_bugged += moments
    elif key in zprocs_bugfix:
        if z_moments_bugfix is None:
            z_moments_bugfix = moments
        else:
            z_moments_bugfix += moments
    elif key in wprocs_bugfix:
        if w_moments_bugfix is None:
            w_moments_bugfix = moments
        else:
            w_moments_bugfix += moments


if z_moments_bugged:
    z_coeffs_bugged = wremnants.moments_to_angular_coeffs(z_moments_bugged)
    with lz4.frame.open("z_coeffs_bugged.pkl.lz4", "wb") as f:
        pickle.dump(z_coeffs_bugged, f, protocol = pickle.HIGHEST_PROTOCOL)
if w_moments_bugged:
    w_coeffs_bugged = wremnants.moments_to_angular_coeffs(w_moments_bugged)
    with lz4.frame.open("w_coeffs_bugged.pkl.lz4", "wb") as f:
        pickle.dump(w_coeffs_bugged, f, protocol = pickle.HIGHEST_PROTOCOL)
if w_moments_bugfix:
    w_coeffs_bugfix = wremnants.moments_to_angular_coeffs(w_moments_bugfix)
    with lz4.frame.open("w_coeffs_bugfix.pkl.lz4", "wb") as f:
        pickle.dump(w_coeffs_bugfix, f, protocol = pickle.HIGHEST_PROTOCOL)
if z_moments_bugfix:
    z_coeffs_bugfix = wremnants.moments_to_angular_coeffs(z_moments_bugfix)
    with lz4.frame.open("z_coeffs_bugfix.pkl.lz4", "wb") as f:
        pickle.dump(z_coeffs_bugfix, f, protocol = pickle.HIGHEST_PROTOCOL)
