import ROOT
import hist
import numpy as np
import copy
from utilities import boostHistHelpers as hh,common
from wremnants import theory_corrections
import logging
from scipy import ndimage

ROOT.gInterpreter.Declare('#include "theoryTools.h"')

# integer axis for -1 through 7
axis_helicity = hist.axis.Integer(
    -1, 8, name="helicity", overflow=False, underflow=False
)

# this puts the bin centers at 0.5, 1.0, 2.0
axis_muRfact = hist.axis.Variable(
    [0.25, 0.75, 1.25, 2.75], name="muRfact", underflow=False, overflow=False
)
axis_muFfact = hist.axis.Variable(
    [0.25, 0.75, 1.25, 2.75], name="muFfact", underflow=False, overflow=False
)

scale_tensor_axes = (axis_muRfact, axis_muFfact)

pdfMap = {
    "nnpdf31" : {
        "name" : "pdfNNPDF31",
        "branch" : "LHEPdfWeight",
        "combine" : "symHessian",
        "entries" : 101,
        "alphas" : ["LHEPdfWeight[101]", "LHEPdfWeight[102]"],
    },
    "ct18" : {
        # This has CT18 + CT18Z in it :-/
        "name" : "pdfCT18",
        "branch" : "LHEPdfWeightAltSet18",
        "combine" : "asymHessian",
        "entries" : 59,
        "alphas" : ["LHEPdfWeightAltSet18[59]", "LHEPdfWeightAltSet18[60]"],
    "alphaRange" : "002",
        "scale" : 1/1.645 # Convert from 90% CL to 68%
    },
    "mmht" : {
        "name" : "pdfMMHT",
        "branch" : "LHEPdfWeightAltSet19",
        "combine" : "asymHessian",
        "entries" : 51,
        "alphas" : ["LHEPdfWeightAltSet20[1]", "LHEPdfWeightAltSet20[2]"],
    "alphaRange" : "001",
    },
    "nnpdf30" : {
    "name" : "pdfNNPDF30",
    "branch" : "LHEPdfWeightAltSet13",
        "combine" : "symHessian",
    "entries" : 101,
    "alphas" : ["LHEPdfWeightAltSet15[0]", "LHEPdfWeightAltSet16[0]"],
    "alphaRange" : "001",
    },
}

pdfMapExtended = copy.deepcopy(pdfMap)
pdfMapExtended["ct18"]["branch"] = "LHEPdfWeightAltSet11"
pdfMapExtended["ct18"]["alphas"] = ["LHEPdfWeightAltSet11[59]", "LHEPdfWeightAltSet11[62]"]
pdfMapExtended["mmht"]["branch"] = "LHEPdfWeightAltSet13"
pdfMapExtended["mmht"]["alphas"] = ["LHEPdfWeightAltSet14[1]", "LHEPdfWeightAltSet14[2]"]
pdfMapExtended.update({
    "nnpdf40" : {
        "name" : "pdfNNPDF40",
        "branch" : "LHEPdfWeightAltSet3",
        "combine" : "symHessian",
        "entries" : 53,
        "alphas" : ["LHEPdfWeightAltSet3[51]", "LHEPdfWeightAltSet3[52]"],
        "alphasRange" : "002", # TODO: IS that true?
    },
    "pdf4lhc21" : {
        "name" : "pdfPDF4LHC21",
        "branch" : "LHEPdfWeightAltSet10",
        "combine" : "symHessian",
        "entries" : 41,
        "alphas" : ["LHEPdfWeightAltSet10[41]", "LHEPdfWeightAltSet10[42]"],
        "alphasRange" : "002", # TODO: IS that true?
    },
    "msht20" : {
        "name" : "pdfMSHT20",
        "branch" : "LHEPdfWeightAltSet12",
        "combine" : "symHessian",
        "entries" : 51,
        "alphas" : ["LHEPdfWeight[51]", "LHEPdfWeightAltSet12[52]"],
        "alphasRange" : "002", # TODO: IS that true?
    },
    "atlasWZj20" : {
        "name" : "pdfATLASWZJ20",
        "branch" : "LHEPdfWeightAltSet19",
        "combine" : "symHessian",
        "entries" : 33,
        "alphas" : ["LHEPdfWeight[41]", "LHEPdfWeight[42]"],
        "alphasRange" : "002", # TODO: IS that true?
    },


})

only_central_pdf_datasets = [
    "Wplusmunu_bugfix",
    "Wminusmunu_bugfix",
    "Zmumu_bugfix",
    "Zmumu_bugfix_slc7",
]

extended_pdf_datasets = common.vprocs

def define_prefsr_vars(df):
    df = df.Define("prefsrLeps", "wrem::prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother)")
    df = df.Define("genl", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[prefsrLeps[0]], GenPart_eta[prefsrLeps[0]], GenPart_phi[prefsrLeps[0]], GenPart_mass[prefsrLeps[0]])")
    df = df.Define("genlanti", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[prefsrLeps[1]], GenPart_eta[prefsrLeps[1]], GenPart_phi[prefsrLeps[1]], GenPart_mass[prefsrLeps[1]])")
    df = df.Define("genV", "ROOT::Math::PxPyPzEVector(genl)+ROOT::Math::PxPyPzEVector(genlanti)")
    df = df.Define("ptVgen", "genV.pt()")
    df = df.Define("massVgen", "genV.mass()")
    df = df.Define("yVgen", "genV.Rapidity()")
    df = df.Define("phiVgen", "genV.Phi()")
    df = df.Define("absYVgen", "std::fabs(yVgen)")
    df = df.Define("chargeVgen", "GenPart_pdgId[prefsrLeps[0]] + GenPart_pdgId[prefsrLeps[1]]")
    df = df.Define("csSineCosThetaPhi", "wrem::csSineCosThetaPhi(genl, genlanti)")
    return df

def define_scale_tensor(df, clipWeight=10.0):
    # convert vector of scale weights to 3x3 tensor and clip weights to |weight|<10.
    df = df.Define("scaleWeights_tensor", f"wrem::makeScaleTensor(LHEScaleWeight, {clipWeight});")
    df = df.Define("scaleWeights_tensor_wnom", "auto res = scaleWeights_tensor; res = nominal_weight*res; return res;")

    return df

def make_scale_hist(df, axes, cols, hname=""):
    scaleHist = df.HistoBoost("qcdScale" if hname=="" else f"{hname}_qcdScale", axes, [*cols, "scaleWeights_tensor_wnom"], tensor_axes=scale_tensor_axes)
    return scaleHist

def pdf_info_map(dataset, pdfset):
    infoMap = pdfMap if dataset not in extended_pdf_datasets else pdfMapExtended

    if (pdfset != "nnpdf31" and dataset in only_central_pdf_datasets) or pdfset not in infoMap:
        raise ValueError(f"Skipping PDF {pdfset} for dataset {dataset}")
    return infoMap[pdfset]

def make_pdf_hists(df, dataset, axes, cols, pdfs, hname=""):
    res = []
    for pdf in pdfs:
        try:
            pdfInfo = pdf_info_map(dataset, pdf)
        except ValueError as e:
            logging.info(e)
            continue

        pdfName = pdfInfo["name"]
        tensorName = f"{pdfName}Weights_tensor"
        tensorASName = f"{pdfName}ASWeights_tensor"
        pdfHist = df.HistoBoost(pdfName if hname=="" else f"{hname}_{pdfName}", axes, [*cols, tensorName])

        alphaSHist = df.HistoBoost(f"alphaS002{pdfName}" if hname=="" else f"{hname}_alphaS002{pdfName}", axes, [*cols, tensorASName])
        res.extend([pdfHist, alphaSHist])
    return res

def define_pdf_columns(df, dataset, pdfs, noAltUnc):
    for i, pdf in enumerate(pdfs):
        try:
            pdfInfo = pdf_info_map(dataset, pdf)
        except ValueError as e:
            logging.info(e)
            return df

        pdfName = pdfInfo["name"]
        pdfBranch = pdfInfo["branch"]
        tensorName = f"{pdfName}Weights_tensor"
        tensorASName = f"{pdfName}ASWeights_tensor"
        entries = 1 if i != 0 and noAltUnc else pdfInfo["entries"]

        df = df.Define(tensorName, f"auto res = wrem::clip_tensor(wrem::vec_to_tensor_t<double, {entries}>({pdfBranch}), 10.); res = nominal_weight/nominal_pdf_cen*res; return res;")

        df = df.Define(tensorASName, "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; "
                f"res(0) = nominal_weight*nominal_pdf_cen*{pdfInfo['alphas'][0]}; "
                f"res(1) = nominal_weight*nominal_pdf_cen*{pdfInfo['alphas'][1]}; "
                "return wrem::clip_tensor(res, 10.)")
    return df

def define_weights_and_corrs(df, weight_expr, dataset_name, helpers, args):
    #TODO: organize this better
    if dataset_name in common.vprocs:
        df = df.Define("nominal_pdf_cen", pdf_central_weight(dataset_name, args.pdfs[0]))
        weight_expr = f"{weight_expr}*nominal_pdf_cen"
    df = define_prefsr_vars(df)

    if args.theory_corr and dataset_name in helpers:
        helper = helpers[dataset_name]
        df = define_theory_corr(df, weight_expr, helper, generators=args.theory_corr, 
                modify_central_weight=not args.theory_corr_alt_only)
    else:
        df = df.Define("nominal_weight", weight_expr)
    return df 

def pdf_central_weight(dataset, pdfset):
    pdfInfo = pdf_info_map(dataset, pdfset)
    pdfBranch = pdfInfo["branch"]
    return f"{pdfBranch}[0]"

def define_theory_corr(df, weight_expr, helpers, generators, modify_central_weight):
    for i, generator in enumerate(generators):
        if i == 0:
            if modify_central_weight and generator in helpers:
                df = df.Define("nominal_weight_uncorr", weight_expr)
            else:
                df = df.Define("nominal_weight", weight_expr)
                df = df.Alias("nominal_weight_uncorr", "nominal_weight")
        if generator not in helpers:
            continue

        helper = helpers[generator]

        if "Helicity" in generator:
            df = df.Define(f"{generator}Weight_tensor", helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "nominal_weight_uncorr"])
        else:
            df = df.Define(f"{generator}Weight_tensor", helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "nominal_weight_uncorr"])

        df = df.Define(f"{generator}CentralWeight", f"{generator}Weight_tensor(0)")

        if i == 0 and modify_central_weight:
            df = df.Alias("nominal_weight", f"{generator}CentralWeight")

    return df

def make_theory_corr_hists(df, name, axes, cols, helpers, generators, modify_central_weight, with_uncertainties=False):
    res = []
    
    for i, generator in enumerate(generators):
        if generator not in helpers:
            continue
        helper = helpers[generator]
        if i == 0 and modify_central_weight:
            nominal_uncorr = df.HistoBoost(f"{name}_uncorr", axes, [*cols, "nominal_weight_uncorr"])
            res.append(nominal_uncorr)
            res.append(df.HistoBoost("weight_uncorr", [hist.axis.Regular(100, -2, 2)], ["nominal_weight_uncorr"]))

        hist_name = f"{generator}Corr"
        if name != "nominal":
            hist_name = f"{name}_{hist_name}"

        if with_uncertainties:
            hist_name += "_unc"
            unc = df.HistoBoost(hist_name, axes, [*cols, "{generator}Weight_tensor"], tensor_axes=helper.tensor_axes)
            res.append(unc)
        else:
            nominal = df.HistoBoost(hist_name, axes, [*cols, f"{generator}CentralWeight"])
            res.append(nominal)

    return res

def scale_angular_moments(hist_moments_scales):
    # e.g. from arxiv:1708.00008 eq. 2.13, note A_0 is NOT the const term!
    scales = np.array([1., -10., 5., 10., 4., 4., 5., 5., 4.])

    hel_idx = hist_moments_scales.axes.name.index("helicity")
    scaled_vals = np.moveaxis(hist_moments_scales.view(flow=True), hel_idx, -1)*scales
    hnew = hist.Hist(*hist_moments_scales.axes, storage=hist_moments_scales._storage_type())
    hnew[...] = np.moveaxis(scaled_vals, -1, hel_idx) 
    return hnew

def replace_by_neighbors(vals, replace):
    if np.count_nonzero(replace) == vals.size:
        raise ValueError("Cannot replace all values with nearest non-zero neighbour")

    indices = ndimage.distance_transform_edt(replace, return_distances=False, return_indices=True)
    return vals[tuple(indices)]

def moments_to_angular_coeffs(hist_moments_scales, cutoff=1e-5):
    if hist_moments_scales.sum().value == 0:
       raise ValueError("Cannot make coefficients from empty hist")
    # broadcasting happens right to left, so move to rightmost then move back
    hel_ax = hist_moments_scales.axes["helicity"]
    hel_idx = hist_moments_scales.axes.name.index("helicity")
    vals = np.moveaxis(scale_angular_moments(hist_moments_scales).view(flow=True), hel_idx, -1) 
    
    # select constant term, leaving dummy axis for broadcasting
    unpol_idx = hel_ax.index(-1)
    norm_vals = vals[...,unpol_idx:unpol_idx+1].value
    norm_vals = np.where(np.abs(norm_vals) < cutoff, np.ones_like(norm_vals), norm_vals)

    # e.g. from arxiv:1708.00008 eq. 2.13, note A_0 is NOT the const term!
    offsets = np.array([0., 4., 0., 0., 0., 0., 0., 0., 0.])

    coeffs = vals / norm_vals + offsets

    # replace values in zero-xsec regions (otherwise A0 is spuriously set to 4.0 from offset)
    coeffs = np.where(np.abs(vals.value) < cutoff, np.full_like(vals, hist.accumulators.WeightedSum(0,0)), coeffs)
    coeffs = np.moveaxis(coeffs, -1, hel_idx)

    hist_coeffs_scales = hist.Hist(*hist_moments_scales.axes, storage = hist_moments_scales._storage_type(), name = "hist_coeffs_scales",
        data = coeffs
    )

    return hist_coeffs_scales

def qcdByHelicityLabels():
    coeffs = ["const"]+[f"a{i}" for i in range(8)]
    scaleVars = ["muRmuF", "muR", "muF"]
    return  [f"{var}_{coeff}{t}" for var in scaleVars for t in ["Up", "Down"] for coeff in coeffs]

def qcdScaleNames():
    # Exclude central and extreme variations
    shifts = ["muRmuFDown", "muRDown", "", "muFDown", "", "muFUp", "muRUp", "", "muRmuFUp"]
    return ["_".join(["QCDscale", s]) if s != "" else s for s in shifts]

def massWeightNames(matches=None, wlike=False):
    central=10
    nweights=21
    names = [f"massShift{int(abs(central-i)*10)}MeV{'Down' if i < central else 'Up'}" for i in range(nweights)]
    
    if wlike:
        # This is the PDG uncertainty
        names.extend(["massShift2p1MeVDown", "massShift2p1MeVUp"])

    # If name is "" it won't be stored
    return [x if not matches or any(y in x for y in matches) else "" for x in names]

def pdfNames(cardTool, pdf, skipFirst=True):
    size = 101
    names = cardTool.mirrorNames(f"pdf{{i}}{pdf}", size)
    if skipFirst:
        names[0] = ""
        names[size] = ""
    # TODO: This is probably not needed anymore, check with low PU
    if False and pdf == "NNPDF31":
        names[size-2] = "pdfAlphas002Up"
        names[size-1] = "pdfAlphas002Down"
        # Drop the mirrored alphaS variations
        names[size*2-2] = ""
        names[size*2-1] = ""
    return names

def pdfNamesAsymHessian(entries):
    pdfNames = [""] # Skip central weight
    pdfNames.extend(["pdf{i}{shift}".format(i=int(j/2), shift="Up" if j % 2 else "Down") for j in range(entries-1)])
    return pdfNames

def pdfSymmetricShifts(hdiff, axis_name):
    sq = hh.multiplyHists(hdiff, hdiff)
    ss = sq[{axis_name : hist.sum}]
    rss = hh.sqrtHist(ss)
    return rss, rss

def pdfAsymmetricShifts(hdiff, axis_name):
    # Assuming that the last axis is the syst axis
    # TODO: add some check to verify this
    def shiftHist(vals, varis, hdiff, axis_name):
        hnew = hdiff.copy()[{axis_name : hist.sum}]
        vals, varis = hh.multiplyWithVariance(vals, vals, varis, varis)
        ss = np.stack((np.sum(vals, axis=-1), np.sum(varis, axis=-1)), axis=-1)
        hnew[...] = ss
        return hh.sqrtHist(hnew)

    # The error sets are ordered up,down,up,down...
    upshift = shiftHist(hdiff.values(flow=True)[...,1::2], hdiff.variances(flow=True)[...,1::2], hdiff, axis_name)
    downshift = shiftHist(hdiff.values(flow=True)[...,2::2], hdiff.variances(flow=True)[...,2::2], hdiff, axis_name)
    return upshift, downshift

def hessianPdfUnc(h, axis_name="tensor_axis_0", uncType="symHessian", scale=1.):
    symmetric = uncType == "symHessian"
    diff = hh.addHists(h, -1*h[{axis_name : 0}])*scale
    shiftFunc = pdfSymmetricShifts if symmetric else pdfAsymmetricShifts
    rssUp, rssDown = shiftFunc(diff, axis_name)
    hUp = hh.addHists(h[{axis_name : 0}], rssUp)
    hDown = hh.addHists(h[{axis_name : 0}], -1*rssDown)
    return hUp, hDown
