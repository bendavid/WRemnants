import ROOT
import hist
import numpy as np
import copy
from utilities import boostHistHelpers as hh,common,logging
from wremnants import theory_corrections
from scipy import ndimage

logger = logging.child_logger(__name__)
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
        "alphasRange" : "001", # TODO: Check this
    },
    "ct18" : {
        # This has CT18 + CT18Z in it :-/
        "name" : "pdfCT18",
        "branch" : "LHEPdfWeightAltSet18",
        "combine" : "asymHessian",
        "entries" : 59,
        "alphas" : ["LHEPdfWeightAltSet18[59]", "LHEPdfWeightAltSet18[60]"],
        "alphasRange" : "002",
        "scale" : 1/1.645 # Convert from 90% CL to 68%
    },
    "mmht" : {
        "name" : "pdfMMHT",
        "branch" : "LHEPdfWeightAltSet19",
        "combine" : "asymHessian",
        "entries" : 51,
        "alphas" : ["LHEPdfWeightAltSet20[1]", "LHEPdfWeightAltSet20[2]"],
        "alphasRange" : "001",
    },
    "nnpdf30" : {
        "name" : "pdfNNPDF30",
        "branch" : "LHEPdfWeightAltSet13",
        "combine" : "symHessian",
        "entries" : 101,
        "alphas" : ["LHEPdfWeightAltSet15[0]", "LHEPdfWeightAltSet16[0]"],
        "alphasRange" : "001",
    },
}

pdfMapExtended = copy.deepcopy(pdfMap)
pdfMapExtended["ct18"]["branch"] = "LHEPdfWeightAltSet11"
pdfMapExtended["ct18"]["alphas"] = ["LHEPdfWeightAltSet11[59]", "LHEPdfWeightAltSet11[62]"]
pdfMapExtended["mmht"]["branch"] = "LHEPdfWeightAltSet13"
pdfMapExtended["mmht"]["alphas"] = ["LHEPdfWeightAltSet13[1]", "LHEPdfWeightAltSet13[2]"]
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
        "combine" : "asymHessian",
        "entries" : 65,
        "alphas" : ["LHEPdfWeightAltSet12[67]", "LHEPdfWeightAltSet12[70]"],
        # 66-71 - are LHAPDF ID 27500 = 27506, 27501 is 0.0116 and 27504 is 0.0120
        "alphasRange" : "002", 
    },
    "ct18z" : {
        # This has CT18 + CT18Z in it :-/
        "name" : "pdfCT18Z",
        "branch" : "LHEPdfWeightAltSet18",
        "combine" : "asymHessian",
        "entries" : 59,
        "first_entry" : 63,
        "alphas" : ["LHEPdfWeightAltSet18[122]", "LHEPdfWeightAltSet18[123]"],
        "alphasRange" : "002",
        "scale" : 1/1.645 # Convert from 90% CL to 68%
    },
    "atlasWZj20" : {
        "name" : "pdfATLASWZJ20",
        "branch" : "LHEPdfWeightAltSet19",
        "combine" : "asymHessian",
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

extended_pdf_datasets = [x for x in common.vprocs+common.vprocs_lowpu if not any(y in x for y in ["NNLOPS", "MiNLO"])]

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

def define_scale_tensor(df):
    # convert vector of scale weights to 3x3 tensor and clip weights to |weight|<10.
    df = df.Define("scaleWeights_tensor", f"wrem::makeScaleTensor(LHEScaleWeight, theory_weight_truncate);")
    df = df.Define("scaleWeights_tensor_wnom", "auto res = scaleWeights_tensor; res = nominal_weight*res; return res;")

    return df

def define_ew_vars(df):
    df = df.Define("ewLeptons", "wrem::ewLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi)")
    df = df.Define("ewPhotons", "wrem::ewPhotons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi)")
    df = df.Define('ewMll', '(ewLeptons[0]+ewLeptons[1]).mass()')
    df = df.Define('ewMlly', 'wrem::ewMLepPhos(ewLeptons, ewPhotons)')
    df = df.Define('ewLogDeltaM', 'log10(ewMlly-ewMll)')

    return df

def make_ew_binning(mass = 91.1535, width = 2.4932, initialStep = 0.1):
    maxVal = ROOT.Math.breitwigner_pdf(mass, width, mass)
    bins = [mass]
    currentMass = mass
    while currentMass - mass < 100:
        binSize = maxVal / ROOT.Math.breitwigner_pdf(currentMass, width, mass) * initialStep
        currentMass += binSize
        bins.append(currentMass)
        lowMass = 2*mass - currentMass
        if lowMass - binSize > 0:
            bins.insert(0, lowMass)
    bins.insert(0, 0.)
    return bins

def pdf_info_map(dataset, pdfset):
    infoMap = pdfMap if dataset not in extended_pdf_datasets else pdfMapExtended

    if "horace" in dataset or (pdfset != "nnpdf31" and dataset in only_central_pdf_datasets) or pdfset not in infoMap:
        raise ValueError(f"Skipping PDF {pdfset} for dataset {dataset}")
    return infoMap[pdfset]

def define_pdf_columns(df, dataset_name, pdfs, noAltUnc):
    if dataset_name not in common.vprocs_all or \
            "horace" in dataset_name or \
            "LHEPdfWeight" not in df.GetColumnNames():
        logger.warning(f"Did not find PDF weights for sample {dataset_name}! Using nominal PDF in sample")
        return df

    for i, pdf in enumerate(pdfs):
        try:
            pdfInfo = pdf_info_map(dataset_name, pdf)
        except ValueError as e:
            return df

        pdfName = pdfInfo["name"]
        pdfBranch = pdfInfo["branch"]
        tensorName = f"{pdfName}Weights_tensor"
        tensorASName = f"{pdfName}ASWeights_tensor"
        entries = 1 if i != 0 and noAltUnc else pdfInfo["entries"]
        start = 0 if "first_entry" not in pdfInfo else pdfInfo["first_entry"]

        df = df.Define(tensorName, f"auto res = wrem::clip_tensor(wrem::vec_to_tensor_t<double, {entries}>({pdfBranch}, {start}), theory_weight_truncate); res = nominal_weight/central_pdf_weight*res; return res;")

        if i == 0:
            tensorNameNominal = tensorName

        if pdfName == "pdfMSHT20":
            df = pdfBugfixMSHT20(df, tensorName)

        df = df.Define(tensorASName, "Eigen::TensorFixedSize<double, Eigen::Sizes<2>> res; "
                f"res(0) = nominal_weight/central_pdf_weight*{pdfInfo['alphas'][0]}; "
                f"res(1) = nominal_weight/central_pdf_weight*{pdfInfo['alphas'][1]}; "
                "return wrem::clip_tensor(res, theory_weight_truncate)")

    return df

def define_central_pdf_weight(df, dataset_name, pdf):
    try:
        pdfInfo = pdf_info_map(dataset_name, pdf)
    except ValueError as e:
        logger.warning(f"Did not find PDF {pdf} for sample {dataset_name}! Using nominal PDF in sample")
        return df.DefinePerSample("central_pdf_weight", "1.0")

    pdfName = pdfInfo["name"]
    pdfBranch = pdfInfo["branch"]
    #return df.Define("central_pdf_weight", f"std::clamp({pdfBranch}[0], -1*theory_weight_truncate, theory_weight_truncate)")
    return df.Define("central_pdf_weight", f"std::clamp<float>({pdfBranch}[0], -10., 10.)")

def define_theory_weights_and_corrs(df, dataset_name, helpers, args):
    #TODO: Reimplement this
    #if args.highptscales:
    #    weight_expr = f"{weight_expr}*MEParamWeightAltSet3[0]"
    df = define_prefsr_vars(df)
    df = define_ew_vars(df)

    df = df.DefinePerSample("theory_weight_truncate", "10.")
    df = define_central_pdf_weight(df, dataset_name, args.pdfs[0])
    df = define_theory_corr(df, dataset_name, helpers, generators=args.theoryCorr, 
            modify_central_weight=not args.theoryCorrAltOnly)
    df = df.Define("nominal_weight", "exp_weight*central_pdf_weight*theory_corr_weight")
    df = define_pdf_columns(df, dataset_name, args.pdfs, args.altPdfOnlyCentral)
        
    return df 

def define_theory_corr(df, dataset_name, helpers, generators, modify_central_weight):
    df = df.Define(f"nominal_weight_uncorr", "exp_weight*central_pdf_weight")

    sample_helpers = helpers[dataset_name] if dataset_name in helpers else None

    if not modify_central_weight or not sample_helpers or generators[0] not in sample_helpers:
        df = df.DefinePerSample("theory_corr_weight", "1.0")
        return df

    for i, generator in enumerate(generators):
        if generator not in sample_helpers:
            continue

        helper = sample_helpers[generator]

        if "Helicity" in generator:
            df = df.Define(f"{generator}Weight_tensor", helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "nominal_weight_uncorr"])
        elif 'ew' in generator:
            df = df.DefinePerSample(f"{generator}Dummy", "0.")
            multiplicative_weight = i != 0 and modify_central_weight
            df = df.Define(f"{generator}Weight_tensor", helper, ["ewMll", "ewLogDeltaM", f"{generator}Dummy", "chargeVgen", "nominal_weight_uncorr" if not multiplicative_weight else "nominal_weight"]) # multiplying with nominal QCD weight
        else:
            df = df.Define(f"{generator}Weight_tensor", helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "nominal_weight_uncorr"])

        if i == 0 and modify_central_weight:
            df = df.Define("theory_corr_weight", f"{generator}Weight_tensor(0)/nominal_weight_uncorr")

    return df

def make_theory_corr_hists(df, name, axes, cols, helpers, generators, modify_central_weight):
    res = []
    
    for i, generator in enumerate(generators):
        if generator not in helpers:
            continue
        
        if i == 0 and modify_central_weight:
            nominal_uncorr = df.HistoBoost(f"{name}_uncorr", axes, [*cols, "nominal_weight_uncorr"])
            res.append(nominal_uncorr)
            res.append(df.HistoBoost("weight_uncorr", [hist.axis.Regular(100, -2, 2)], ["nominal_weight_uncorr"]))

        hist_name = f"{name}_{generator}Corr"
        unc = df.HistoBoost(hist_name, axes, [*cols, f"{generator}Weight_tensor"], tensor_axes=helpers[generator].tensor_axes[-1:])
        res.append(unc)

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

def pdfNamesAsymHessian(entries, pdfset=""):
    pdfNames = ["pdf0"+pdfset.replace("pdf", "")] 
    pdfNames.extend([f"pdf{int((j+2)/2)}{pdfset.replace('pdf', '')}{'Up' if j % 2 else 'Down'}" for j in range(entries-1)])
    return pdfNames

def pdfNamesSymHessian(entries, pdfset=""):
    return [f"pdf{i+1}{pdfset.replace('pdf', '')}" for i in range(entries)]

def pdfSymmetricShifts(hdiff, axis_name):
    sq = hh.multiplyHists(hdiff, hdiff)
    ss = sq[{axis_name : hist.sum}]
    rss = hh.sqrtHist(ss)
    return rss, rss

def pdfAsymmetricShifts(hdiff, axis_name):
    # Assuming that the last axis is the syst axis
    # TODO: add some check to verify this
    def shiftHist(vals, varis, hdiff, axis_name):
        hnew = hdiff.copy()[{axis_name : 0}]
        vals, varis = hh.multiplyWithVariance(vals, vals, varis, varis)
        ss = np.stack((np.sum(vals, axis=-1), np.sum(varis, axis=-1)), axis=-1)
        hnew[...] = ss
        return hh.sqrtHist(hnew)

    ax = hdiff.axes[axis_name] 
    underflow = hdiff.axes[axis_name].traits.underflow
    overflow = hdiff.axes[axis_name].traits.overflow
    if type(ax) == hist.axis.StrCategory and all(["Up" in x or "Down" in x for x in ax][1:]):
        # Remove the overflow from the categorical axis
        end = int((ax.size-1)/2)
        upvals = hdiff[{axis_name : [x for x in ax if "Up" in x]}].values(flow=True)[...,:end]
        upvars = hdiff[{axis_name : [x for x in ax if "Up" in x]}].variances(flow=True)[...,:end]
        downvals = hdiff[{axis_name : [x for x in ax if "Down" in x]}].values(flow=True)[...,:end]
        downvars = hdiff[{axis_name : [x for x in ax if "Down" in x]}].variances(flow=True)[...,:end]
        if upvals.shape != downvals.shape:
            raise ValueError("Malformed PDF uncertainty hist! Expect equal number of up and down vars")
    else:
        end = ax.size+underflow
        upvals = hdiff.values(flow=True)[...,1+underflow:end:2]
        upvars = hdiff.variances(flow=True)[...,1+underflow:end:2]
        downvals = hdiff.values(flow=True)[...,2+underflow:end:2]
        downvars = hdiff.variances(flow=True)[...,2+underflow:end:2]

    # The error sets are ordered up,down,up,down...
    upshift = shiftHist(upvals, upvars, hdiff, axis_name)
    downshift = shiftHist(downvals, downvars, hdiff, axis_name)
    return upshift, downshift 

def hessianPdfUnc(h, axis_name="pdfVar", uncType="symHessian", scale=1.):
    underflow = h.axes[axis_name].traits.underflow
    symmetric = uncType == "symHessian"
    diff = hh.addHists(h, -1*h[{axis_name : 0}])*scale
    if diff.axes[axis_name].traits.overflow:
        diff[...,hist.overflow] = np.zeros_like(diff[{axis_name : 0}].view(flow=True))
    shiftFunc = pdfSymmetricShifts if symmetric else pdfAsymmetricShifts
    rssUp, rssDown = shiftFunc(diff, axis_name)
    hUp = hh.addHists(h[{axis_name : 0}], 1*rssUp)
    hDown = hh.addHists(h[{axis_name : 0}], -1*rssDown)
    return hUp, hDown

def pdfBugfixMSHT20(df , tensorPDFName):
    # There is a known bug in MSHT20 where member 15 and 16 are identical
    #   to fix this, one has to be mirrored:
    #   pdf(15) = pdf(0) - (pdf(15) - pdf(0))
    return df.Redefine(tensorPDFName, 
        f"auto res = {tensorPDFName};"
        f"res(15) = {tensorPDFName}(0) - ({tensorPDFName}(15) - {tensorPDFName}(0));"
        "return res")
        
