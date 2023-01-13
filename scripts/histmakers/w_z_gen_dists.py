from utilities import boostHistHelpers as hh, common, output_tools

parser,initargs = common.common_parser()

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist
import math


parser.add_argument("--skipAngularCoeffs", action='store_true', help="Skip the conversion of helicity moments to angular coeff fractions")
parser.add_argument("--singleLeptonHists", action='store_true', help="Also store single lepton kinematics")
parser.add_argument("--ewHists", action='store_true', help="Also store histograms for EW reweighting. Use with --filter horace")

f = next((x for x in parser._actions if x.dest == "filterProcs"), None)
if f:
    f.default = common.vprocs

args = parser.parse_args()

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None, 
    nanoVersion="v8" if args.v8 else "v9", base_path=args.data_path, mode='gen')

axis_massWgen = hist.axis.Variable([5., 13000.], name="massVgen", underflow=True, overflow=False)

#axis_massZgen = hist.axis.Regular(12, 60., 120., name="massVgen")
axis_massZgen = hist.axis.Regular(10, 60., 120., name="massVgen")

axis_absYVgen = hist.axis.Variable(
    [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 5], 
    name = "absYVgen", underflow=False
)

axis_ygen = hist.axis.Regular(200, -5., 5., name="y")

axis_ptVgen = hist.axis.Variable(
    list(range(0,151))+[160., 190.0, 220.0, 250.0, 300.0, 400.0, 500.0, 800.0, 1500.0], 
    name = "ptVgen", underflow=False,
)

axis_chargeWgen = hist.axis.Regular(
    2, -2, 2, name="chargeVgen", underflow=False, overflow=False
)

axis_chargeZgen = hist.axis.Integer(
    0, 1, name="chargeVgen", underflow=False, overflow=False
)

axis_l_eta_gen = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_l_pt_gen = hist.axis.Regular(29, 26., 55., name = "pt")

corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theory_corr)

def build_graph(df, dataset):
    print("build graph")
    print(dataset.name)
    results = []
    
    if dataset.is_data:
        raise RuntimeError("Running GEN analysis over data is not supported")

    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs

    weight_expr = "std::copysign(1.0, genWeight)"
    df = df.Define("weight", weight_expr)
    # This sum should happen before any change of the weight
    weightsum = df.SumAndCount("weight")

    if "reweight_h2" in dataset.name:
        weight_expr = f"{weight_expr}*H2BugFixWeight[0]"
    elif "NNLOPS" in dataset.name:
        weight_expr = f"{weight_expr}*LHEScaleWeightAltSet1[4]"

    df = theory_tools.define_weights_and_corrs(df, weight_expr, dataset.name, corr_helpers, args)

    if isZ:
        nominal_axes = [axis_massZgen, axis_ygen, axis_ptVgen, axis_chargeZgen]
        lep_axes = [axis_l_eta_gen, axis_l_pt_gen, axis_chargeZgen]
    else:
        nominal_axes = [axis_massWgen, axis_ygen, axis_ptVgen, axis_chargeWgen]
        lep_axes = [axis_l_eta_gen, axis_l_pt_gen, axis_chargeWgen]

    nominal_cols = ["massVgen", "yVgen", "ptVgen", "chargeVgen"]
    lep_cols = ["etaPrefsrLep", "ptPrefsrLep", "chargeVgen"]

    if args.singleLeptonHists and isW or isZ:
        if isW:
            df = df.Define('ptPrefsrLep', 'genlanti.pt()')
            df = df.Define('etaPrefsrLep', 'genlanti.eta()')
        else:
            df = df.Define('ptPrefsrLep', 'genl.pt()')
            df = df.Define('etaPrefsrLep', 'genl.eta()')
        results.append(df.HistoBoost("nominal_genlep", lep_axes, [*lep_cols, "nominal_weight"]))

    if args.ewHists and (isW or isZ):
        df = theory_tools.define_ew_vars(df)
        if isZ:
            massBins = theory_tools.make_ew_binning(mass = 91.1535, width = 2.4932, initialStep=0.010)
        else:
            massBins = theory_tools.make_ew_binning(mass = 80.3815, width = 2.0904, initialStep=0.010)
        ew_cols = ['ewMll', 'ewLogDeltaM']
        axis_ewMll = hist.axis.Variable(massBins, name = "ewMll")
        axis_ewLogDeltaM = hist.axis.Regular(100, -5, 5, name = "ewLogDeltaM")
        ew_axes = [axis_ewMll, axis_ewLogDeltaM]
        results.append(df.HistoBoost("nominal_ew", ew_axes, [*ew_cols, "nominal_weight"]))

    nominal_gen = df.HistoBoost("nominal_gen", nominal_axes, [*nominal_cols, "nominal_weight"])

    results.append(nominal_gen)

    if 'horace' in dataset.name:
        return results, weightsum

    if "LHEScaleWeight" in df.GetColumnNames():
        df = theory_tools.define_scale_tensor(df)
        syst_tools.add_scale_hist(results, df, nominal_axes, nominal_cols, "nominal_gen")
        df = df.Define("helicity_moments_scale_tensor", "wrem::makeHelicityMomentScaleTensor(csSineCosThetaPhi, scaleWeights_tensor, nominal_weight)")
        helicity_moments_scale = df.HistoBoost("helicity_moments_scale", nominal_axes, [*nominal_cols, "helicity_moments_scale_tensor"], tensor_axes = [wremnants.axis_helicity, *wremnants.scale_tensor_axes])
        results.append(helicity_moments_scale)

    if "LHEPdfWeight" in df.GetColumnNames():
        df = theory_tools.define_pdf_columns(df, dataset.name, args.pdfs, args.altPdfOnlyCentral)
        syst_tools.add_pdf_hists(results, df, dataset.name, nominal_axes, nominal_cols, args.pdfs, "nominal_gen")

    if args.theory_corr and dataset.name in corr_helpers:
        results.extend(theory_tools.make_theory_corr_hists(df, "nominal_gen", nominal_axes, nominal_cols,
            corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
        )
        if args.singleLeptonHists:
            results.extend(theory_tools.make_theory_corr_hists(df, "nominal_genlep", lep_axes, lep_cols, 
                corr_helpers[dataset.name], args.theory_corr, modify_central_weight=not args.theory_corr_alt_only)
            )

    if "MEParamWeight" in df.GetColumnNames():
        df = syst_tools.define_mass_weights(df, isW)
        syst_tools.add_massweights_hist(results, df, nominal_axes, nominal_cols, "nominal_gen")

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

output_tools.write_analysis_output(resultdict, "w_z_gen_dists.pkl.lz4", args)

print("computing angular coefficients")

z_moments = None
w_moments = None

if not args.skipAngularCoeffs:
    for dataset in datasets:
        name = dataset.name
        moments = resultdict[name]["output"]["helicity_moments_scale"]
        if name in common.zprocs:
            if z_moments is None:
                z_moments = moments
            else:
                z_moments += moments
        elif name in common.wprocs:
            if w_moments is None:
                w_moments = moments
            else:
                w_moments += moments

if z_moments and w_moments:
    # REMINDER: common.ptV_binning is not the one using 10% quantiles, and the quantiles are not a subset of this binning, but apparently it doesn't matter
    yax_name = axis_absYVgen.name
    if "y" in z_moments.axes.name and "y" in w_moments.axes.name:
        z_moments = hh.makeAbsHist(z_moments, "y")
        w_moments = hh.makeAbsHist(w_moments, "y")
        yax_name = "absy"

    z_moments = hh.rebinHist(z_moments, axis_ptVgen.name, common.ptV_binning)
    z_moments = hh.rebinHist(z_moments, axis_massZgen.name, axis_massZgen.edges[::2])
    z_moments = hh.rebinHist(z_moments, yax_name, axis_absYVgen.edges[:-1])
    w_moments = hh.rebinHist(w_moments, axis_ptVgen.name, common.ptV_binning)
    w_moments = hh.rebinHist(w_moments, yax_name, axis_absYVgen.edges[:-1])

    coeffs = {"Z" : wremnants.moments_to_angular_coeffs(z_moments) if z_moments else None,
            "W" : wremnants.moments_to_angular_coeffs(w_moments) if w_moments else None,
    }

    output_tools.write_analysis_output(coeffs, "w_z_coeffs.pkl.lz4", args)
