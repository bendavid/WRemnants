from wremnants.utilities import binning, common, parsing, samples
from wums import logging

parser, initargs = parsing.common_parser("w_mass")

import os

import hist
import ROOT

import narf
from wremnants.production import (
    muon_calibration,
    muon_selections,
    pileup,
    vertex,
)
from wremnants.production.datasets.dataset_tools import getDatasets
from wremnants.production.histmaker_tools import write_analysis_output

parser.add_argument(
    "--testHelpers", action="store_true", help="Test the smearing weights helper"
)

args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    extended="msht20an3lo" not in args.pdfs,
    nanoVersion="v9",
    base_path=args.dataPath,
)

era = args.era


axis_genPt = hist.axis.Regular(45, 9.0, 81.0, name="genPt")
axis_genEta = hist.axis.Regular(50, -2.5, 2.5, name="genEta")
axis_genCharge = hist.axis.Regular(
    2, -2.0, 2.0, underflow=False, overflow=False, name="genCharge"
)
# axis_qopr = hist.axis.Regular(1001, 0.0, 2.0, name="qopr")
axis_qopr = hist.axis.Regular(100, 0.0, 100.0, name="qopr")

axis_eta = hist.axis.Regular(args.eta[0], args.eta[1], args.eta[2], name="eta")
axis_charge = binning.axis_charge
axis_nvalidpixel = hist.axis.Integer(0, 10, name="nvalidpixel")

# response_axes = [axis_genPt, axis_genEta, axis_genCharge, axis_qopr]
response_axes = [axis_genPt, axis_genCharge, axis_qopr]

axis_globalparms = hist.axis.Integer(0, int(109e3), name="globalparms")

pileup_helper = pileup.make_pileup_helper(era=era)
vertex_helper = vertex.make_vertex_helper(era=era)

calib_filepaths = common.calib_filepaths
(
    mc_jpsi_crctn_helper,
    data_jpsi_crctn_helper,
    jpsi_crctn_MC_unc_helper,
    jpsi_crctn_data_unc_helper,
) = muon_calibration.make_jpsi_crctn_helpers(
    args, calib_filepaths, make_uncertainty_helper=True
)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = (
    muon_calibration.make_muon_calibration_helpers(args)
)

smearing_helper, smearing_uncertainty_helper = (
    (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()
)
bias_helper = (
    muon_calibration.make_muon_bias_helpers(args) if args.biasCalibration else None
)

sigmarel = 5e-3
scalerel = 5e-4
nreps = 100

smearing_helper_simple_multi = ROOT.wrem.SmearingHelperSimpleMulti[nreps](sigmarel)

if args.testHelpers:
    response_helper = ROOT.wrem.SplinesDifferentialWeightsHelper(
        f"{common.data_dir}/calibration/muon_response.tflite"
    )

    smearing_helper_simple = ROOT.wrem.SmearingHelperSimple(
        sigmarel, ROOT.ROOT.GetThreadPoolSize()
    )
    smearing_helper_simple_weights = ROOT.wrem.SmearingHelperSimpleWeight(sigmarel)
    smearing_helper_simple_gaussian_weights = ROOT.wrem.SmearingHelperSimpleGaussianWeight(sigmarel)
    smearing_helper_simple_transform = ROOT.wrem.SmearingHelperSimpleTransform(sigmarel)
    scale_helper_simple_weights = ROOT.wrem.ScaleHelperSimpleWeight(scalerel)
    scale_helper_simple_moment_weights = ROOT.wrem.ScaleHelperSimpleMomentWeight(scalerel)
    scale_helper_simple_gaussian_weights = ROOT.wrem.ScaleHelperSimpleGaussianWeight(scalerel)
    gaussian_weight_helper = ROOT.wrem.GaussianWeightHelper()
    module_weight_helper = ROOT.wrem.ModuleWeightHelper()


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in samples.wprocs
    isZ = dataset.name in samples.zprocs
    isTop = dataset.group == "Top"
    isQCDMC = dataset.group == "QCD"

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")
    df = df.Define("isEvenEvent", "event % 2 == 0")

    if dataset.is_data:
        df = df.DefinePerSample("nominal_weight", "1.0")
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])

        weight_expr = "weight*weight_pu*weight_vtx"
        df = df.Define("nominal_weight", weight_expr)

    df = muon_calibration.define_corrected_muons(
        df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper
    )

    df = muon_selections.select_veto_muons(df, nMuons=-1)
    df = muon_selections.select_good_muons(
        df, ptLow=0.0, ptHigh=1e6, nMuons=-1, datasetGroup=None
    )

    df = df.Define(
        "goodTrigObjs",
        f"wrem::goodMuonTriggerCandidate<wrem::Era::Era_{era}>(TrigObj_id,TrigObj_filterBits)",
    )
    df = df.Define(
        "trigMuons",
        "wrem::hasTriggerMatch(Muon_correctedEta,Muon_correctedPhi,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])",
    )

    if not dataset.is_data:
        df = df.Define(
            "Muon_cvhidealMomCov",
            "wrem::splitCovariance<3>(Muon_cvhidealMomCov_Vals, Muon_cvhidealMomCov_Counts)",
        )


        df = df.Define(
            "genMatchedMuons", "Muon_genPartFlav == 1 || Muon_genPartFlav==15"
        )
        df = df.Define("selMuons", "vetoMuonsPre && genMatchedMuons")
        df = df.Define("selGoodMuons", "goodMuons && selMuons")
        df = df.Define("selGoodTrigMuons", "selGoodMuons && trigMuons")

        df = df.Define("selMuons_correctedPt", "Muon_correctedPt[selMuons]")
        df = df.Define("selMuons_correctedEta", "Muon_correctedEta[selMuons]")
        df = df.Define("selMuons_correctedPhi", "Muon_correctedPhi[selMuons]")
        df = df.Define("selMuons_correctedCharge", "Muon_correctedCharge[selMuons]")
        df = df.Define(
            "selMuons_cvhidealNValidPixelHits", "Muon_cvhidealNValidPixelHits[selMuons]"
        )
        df = df.Define("selMuons_genPartIdx", "Muon_genPartIdx[selMuons]")
        df = df.Define("selMuons_genPt", "Take(GenPart_pt, selMuons_genPartIdx)")
        df = df.Define("selMuons_genEta", "Take(GenPart_eta, selMuons_genPartIdx)")
        df = df.Define("selMuons_genPhi", "Take(GenPart_phi, selMuons_genPartIdx)")
        df = df.Define("selMuons_genPdgId", "Take(GenPart_pdgId, selMuons_genPartIdx)")
        df = df.Define(
            "selMuons_genCharge",
            "-1*(selMuons_genPdgId > 0) + 1*(selMuons_genPdgId < 0)",
        )
        df = df.Define(
            "selMuons_qop",
            "selMuons_correctedCharge*1.0/(selMuons_correctedPt*cosh(selMuons_correctedEta))",
        )
        df = df.Define(
            "selMuons_genQop",
            "selMuons_genCharge*1.0/(selMuons_genPt*cosh(selMuons_genEta))",
        )
        df = df.Define("selMuons_qopr", "selMuons_qop/selMuons_genQop")

        df = df.Define("selMuons_cvhidealMomCov", "Muon_cvhidealMomCov[selMuons]")

        df = df.Define("selGoodMuons_correctedPt", "Muon_correctedPt[selGoodMuons]")
        df = df.Define("selGoodMuons_correctedEta", "Muon_correctedEta[selGoodMuons]")
        df = df.Define("selGoodMuons_correctedPhi", "Muon_correctedPhi[selGoodMuons]")
        df = df.Define(
            "selGoodMuons_correctedCharge", "Muon_correctedCharge[selGoodMuons]"
        )
        df = df.Define(
            "selGoodMuons_cvhidealNValidPixelHits",
            "Muon_cvhidealNValidPixelHits[selGoodMuons]",
        )
        df = df.Define("selGoodMuons_genPartIdx", "Muon_genPartIdx[selGoodMuons]")
        df = df.Define(
            "selGoodMuons_genPt", "Take(GenPart_pt, selGoodMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodMuons_genEta", "Take(GenPart_eta, selGoodMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodMuons_genPhi", "Take(GenPart_phi, selGoodMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodMuons_genPdgId", "Take(GenPart_pdgId, selGoodMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodMuons_genCharge",
            "-1*(selGoodMuons_genPdgId > 0) + 1*(selGoodMuons_genPdgId < 0)",
        )
        df = df.Define(
            "selGoodMuons_qop",
            "selGoodMuons_correctedCharge*1.0/(selGoodMuons_correctedPt*cosh(selGoodMuons_correctedEta))",
        )
        df = df.Define(
            "selGoodMuons_genQop",
            "selGoodMuons_genCharge*1.0/(selGoodMuons_genPt*cosh(selGoodMuons_genEta))",
        )
        df = df.Define("selGoodMuons_qopr", "selGoodMuons_qop/selGoodMuons_genQop")

        df = df.Define(
            "selGoodTrigMuons_correctedPt", "Muon_correctedPt[selGoodTrigMuons]"
        )
        df = df.Define(
            "selGoodTrigMuons_correctedEta", "Muon_correctedEta[selGoodTrigMuons]"
        )
        df = df.Define(
            "selGoodTrigMuons_correctedPhi", "Muon_correctedPhi[selGoodTrigMuons]"
        )
        df = df.Define(
            "selGoodTrigMuons_correctedCharge", "Muon_correctedCharge[selGoodTrigMuons]"
        )
        df = df.Define(
            "selGoodTrigMuons_cvhidealNValidPixelHits",
            "Muon_cvhidealNValidPixelHits[selGoodTrigMuons]",
        )
        df = df.Define(
            "selGoodTrigMuons_genPartIdx", "Muon_genPartIdx[selGoodTrigMuons]"
        )
        df = df.Define(
            "selGoodTrigMuons_genPt", "Take(GenPart_pt, selGoodTrigMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodTrigMuons_genEta", "Take(GenPart_eta, selGoodTrigMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodTrigMuons_genPhi", "Take(GenPart_phi, selGoodTrigMuons_genPartIdx)"
        )
        df = df.Define(
            "selGoodTrigMuons_genPdgId",
            "Take(GenPart_pdgId, selGoodTrigMuons_genPartIdx)",
        )
        df = df.Define(
            "selGoodTrigMuons_genCharge",
            "-1*(selGoodTrigMuons_genPdgId > 0) + 1*(selGoodTrigMuons_genPdgId < 0)",
        )
        df = df.Define(
            "selGoodTrigMuons_qop",
            "selGoodTrigMuons_correctedCharge*1.0/(selGoodTrigMuons_correctedPt*cosh(selGoodTrigMuons_correctedEta))",
        )
        df = df.Define(
            "selGoodTrigMuons_genQop",
            "selGoodTrigMuons_genCharge*1.0/(selGoodTrigMuons_genPt*cosh(selGoodTrigMuons_genEta))",
        )
        df = df.Define(
            "selGoodTrigMuons_qopr", "selGoodTrigMuons_qop/selGoodTrigMuons_genQop"
        )

    if isW or isZ:
        response_cols = [
            "selMuons_genPt",
            # "selMuons_genEta",
            "selMuons_genCharge",
            # "selMuons_qopr",
            "selMuons_correctedPt",
        ]
        hist_qopr = df.HistoBoost(
            "hist_qopr", response_axes, [*response_cols, "nominal_weight"]
        )
        results.append(hist_qopr)

        df = df.Define("selMuons_shiftedqopr", f"(1. + {scalerel})*selMuons_qopr")
        df = df.Define("selMuons_shiftedpt", f"(1. - {scalerel}*selMuons_correctedEta*selMuons_correctedEta)*selMuons_correctedPt")

        response_cols_shifted = [
            "selMuons_genPt",
            # "selMuons_genEta",
            "selMuons_genCharge",
            # "selMuons_shiftedqopr",
            "selMuons_shiftedpt",
        ]
        hist_qopr_shifted = df.HistoBoost(
            "hist_qopr_shifted",
            response_axes,
            [*response_cols_shifted, "nominal_weight"],
        )
        hist_qopr_shifted._hist.metadata = {"scalerel": scalerel}
        results.append(hist_qopr_shifted)

        df = df.Define(
            "selMuonsMulti_smearedmqop",
            smearing_helper_simple_multi,
            [
                "run",
                "luminosityBlock",
                "event",
                "selMuons_correctedPt",
                "selMuons_correctedEta",
                "selMuons_correctedCharge",
            ],
        )

        df = df.Define(
            "selMuonsMulti_genPt", f"wrem::replicate_rvec(selMuons_genPt, {nreps})"
        )
        df = df.Define(
            "selMuonsMulti_genEta", f"wrem::replicate_rvec(selMuons_genEta, {nreps})"
        )
        df = df.Define(
            "selMuonsMulti_genCharge",
            f"wrem::replicate_rvec(selMuons_genCharge, {nreps})",
        )
        df = df.Define(
            "selMuonsMulti_genQop", f"wrem::replicate_rvec(selMuons_genQop, {nreps})"
        )
        df = df.Define(
            "selMuonsMulti_smearedqopr",
            "selMuonsMulti_smearedmqop/selMuonsMulti_genQop",
        )

        response_cols_smeared_multi = [
            "selMuonsMulti_genPt",
            "selMuonsMulti_genEta",
            "selMuonsMulti_genCharge",
            "selMuonsMulti_smearedqopr",
        ]
        hist_qopr_smearedmulti = df.HistoBoost(
            "hist_qopr_smearedmulti",
            response_axes,
            [*response_cols_smeared_multi, "nominal_weight"],
        )
        hist_qopr_smearedmulti._hist.metadata = {"sigmarel": sigmarel}
        results.append(hist_qopr_smearedmulti)

        hist_nvalidpixel_gen = df.HistoBoost(
            "hist_nvalidpixel_gen",
            [axis_genEta, axis_genPt, axis_genCharge, axis_nvalidpixel],
            [
                "selMuons_genEta",
                "selMuons_genPt",
                "selMuons_genCharge",
                "selMuons_cvhidealNValidPixelHits",
                "nominal_weight",
            ],
        )
        results.append(hist_nvalidpixel_gen)

        hist_nvalidpixel_nontrig = df.HistoBoost(
            "hist_nvalidpixel_nontrig",
            [axis_eta, axis_genPt, axis_charge, axis_nvalidpixel],
            [
                "selGoodMuons_correctedEta",
                "selGoodMuons_genPt",
                "selGoodMuons_correctedCharge",
                "selGoodMuons_cvhidealNValidPixelHits",
                "nominal_weight",
            ],
        )
        results.append(hist_nvalidpixel_nontrig)

        hist_nvalidpixel_trig = df.HistoBoost(
            "hist_nvalidpixel_trig",
            [axis_eta, axis_genPt, axis_charge, axis_nvalidpixel],
            [
                "selGoodTrigMuons_correctedEta",
                "selGoodTrigMuons_genPt",
                "selGoodTrigMuons_correctedCharge",
                "selGoodTrigMuons_cvhidealNValidPixelHits",
                "nominal_weight",
            ],
        )
        results.append(hist_nvalidpixel_trig)

        if args.testHelpers:

            df = df.Define(
                "selMuons_response_weight",
                response_helper,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_genPt",
                    "selMuons_genEta",
                    "selMuons_genCharge",
                ],
            )

            df = df.Define(
                "weight_smear",
                smearing_helper_simple_weights,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_response_weight",
                    "nominal_weight",
                ],
            )

            df = df.Define(
                "weight_smear_gaussian",
                smearing_helper_simple_gaussian_weights,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_genPt",
                    "selMuons_genEta",
                    "selMuons_genCharge",
                    "selMuons_cvhidealMomCov",
                    "nominal_weight",
                ],
            )

            df = df.DefineSlot(
                "selMuons_smearedPt",
                smearing_helper_simple,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                ],
            )

            df = df.Define(
                "selMuons_smearedqop",
                "selMuons_correctedCharge*1.0/(selMuons_smearedPt*cosh(selMuons_correctedEta))",
            )
            df = df.Define(
                "selMuons_smearedqopr", "selMuons_smearedqop/selMuons_genQop"
            )

            df = df.Define(
                "selMuons_transformedPt",
                smearing_helper_simple_transform,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_response_weight",
                ],
            )

            df = df.Define(
                "selMuons_transformedqop",
                "selMuons_correctedCharge*1.0/(selMuons_transformedPt*cosh(selMuons_correctedEta))",
            )
            df = df.Define(
                "selMuons_transformedqopr", "selMuons_transformedqop/selMuons_genQop"
            )

            df = df.Define(
                "weight_scale",
                scale_helper_simple_weights,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_response_weight",
                    "nominal_weight",
                ],
            )

            df = df.Define(
                "weight_scale_moment",
                scale_helper_simple_moment_weights,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_genPt",
                    "selMuons_genEta",
                    "selMuons_genCharge",
                    "nominal_weight",
                ],
            )

            df = df.Define(
                "weight_scale_gaussian",
                scale_helper_simple_gaussian_weights,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedCharge",
                    "selMuons_genPt",
                    "selMuons_genEta",
                    "selMuons_genCharge",
                    "selMuons_cvhidealMomCov",
                    "nominal_weight",
                ],
            )

            df = df.Define(
                "weight_grads_gaussian",
                gaussian_weight_helper,
                [
                    "selMuons_correctedPt",
                    "selMuons_correctedEta",
                    "selMuons_correctedPhi",
                    "selMuons_correctedCharge",
                    "selMuons_genPt",
                    "selMuons_genEta",
                    "selMuons_genPhi",
                    "selMuons_genCharge",
                    "selMuons_cvhidealMomCov",
                ],
            )

            df = df.Define(
                "module_weights",
                module_weight_helper,
                [
                    "Muon_cvhmergedGlobalIdxs",
                    "Muon_cvhidealJacRef",
                    "weight_grads_gaussian",
                    "nominal_weight",
                ],
            )

            df = df.Define(
                "module_weight_idxs",
                "module_weights.first"
            )

            df = df.Define(
                "module_weight_vals",
                "module_weights.second"
            )

            response_cols_smeared = [
                "selMuons_genPt",
                "selMuons_genEta",
                "selMuons_genCharge",
                "selMuons_smearedqopr",
            ]
            hist_qopr_smeared = df.HistoBoost(
                "hist_qopr_smeared",
                response_axes,
                [*response_cols_smeared, "nominal_weight"],
            )
            results.append(hist_qopr_smeared)

            response_cols_transformed = [
                "selMuons_genPt",
                "selMuons_genEta",
                "selMuons_genCharge",
                "selMuons_transformedqopr",
            ]
            hist_qopr_transformed = df.HistoBoost(
                "hist_qopr_transformed",
                response_axes,
                [*response_cols_transformed, "nominal_weight"],
            )
            results.append(hist_qopr_transformed)

            hist_qopr_smeared_weight = df.HistoBoost(
                "hist_qopr_smeared_weight",
                response_axes,
                [*response_cols, "weight_smear"],
            )
            results.append(hist_qopr_smeared_weight)

            hist_qopr_smeared_weight_gaussian = df.HistoBoost(
                "hist_qopr_smeared_weight_gaussian",
                response_axes,
                [*response_cols, "weight_smear_gaussian"],
            )
            results.append(hist_qopr_smeared_weight_gaussian)

            hist_qopr_scaled_weight = df.HistoBoost(
                "hist_qopr_scaled_weight",
                response_axes,
                [*response_cols, "weight_scale"],
            )
            results.append(hist_qopr_scaled_weight)

            # hist_qopr_scaled_weight_moment = df.HistoBoost(
            #     "hist_qopr_scaled_weight_moment",
            #     response_axes,
            #     [*response_cols, "weight_scale_moment"],
            # )
            hist_qopr_scaled_weight_moment = narf.histutils.shifted_hist(
                df,
                "hist_qopr_scaled_weight_moment",
                response_axes,
                response_cols,
                response_cols_shifted,
                "nominal_weight"
            )
            results.append(hist_qopr_scaled_weight_moment)

            hist_qopr_scaled_weight_gaussian = df.HistoBoost(
                "hist_qopr_scaled_weight_gaussian",
                response_axes,
                [*response_cols, "weight_scale_gaussian"],
            )
            results.append(hist_qopr_scaled_weight_gaussian)

            dfsingle = df.Filter("Sum(selGoodMuons) > 0")
            dfsingle = dfsingle.Define("selMuons0_genPt", "selMuons_genPt[0]")
            dfsingle = dfsingle.Define("selMuons0_genEta", "selMuons_genEta[0]")
            dfsingle = dfsingle.Define("selMuons0_genCharge", "selMuons_genCharge[0]")
            dfsingle = dfsingle.Define("selMuons0_qopr", "selMuons_qopr[0]")

            response_cols_single = [
                "selMuons0_genPt",
                "selMuons0_genEta",
                "selMuons0_genCharge",
                "selMuons0_qopr",
            ]

            hist_qopr_global_parms = dfsingle.HistoBoost(
                "hist_qopr_global_parms",
                [response_axes[-1], axis_globalparms],
                [response_cols_single[-1], "module_weight_idxs", "module_weight_vals"], storage=hist.storage.Double(),
            )
            results.append(hist_qopr_global_parms)

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)
