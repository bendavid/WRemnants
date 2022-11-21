import hist
import numpy as np
from utilities import boostHistHelpers as hh, common
from wremnants import theory_tools
import collections.abc
import logging

logger = common.child_logger(__name__)

def syst_transform_map(base_hist, hist_name):
    pdfInfo = theory_tools.pdfMapExtended 
    pdfNames = [pdfInfo[k]["name"] for k in pdfInfo.keys()]

    def pdfUnc(h, pdfName):
        key =  list(pdfInfo.keys())[list(pdfNames).index(pdfName)]
        unc = pdfInfo[key]["combine"]
        scale = pdfInfo[key]["scale"] if "scale" in pdfInfo[key] else 1.
        return theory_tools.hessianPdfUnc(h, "tensor_axis_0", unc, scale)

    def uncHist(unc):
        return unc if base_hist == "nominal" else f"{base_hist}_{unc}"

    transforms = {}
    transforms.update({pdf+"Up" : {"hist" : uncHist(pdf), "action" : lambda h,p=pdf: pdfUnc(h, p)[0]} for pdf in pdfNames})
    transforms.update({pdf+"Down" : {"hist" : uncHist(pdf), "action" : lambda h,p=pdf: pdfUnc(h, p)[1]} for pdf in pdfNames})
    transforms.update({
        "massShift100MeVDown" : {"hist" : "massWeight", "action" : lambda h: h[{"tensor_axis_0" : 0}]},
        "massShift100MeVUp" : {"hist" : "massWeight", "action" : lambda h: h[{"tensor_axis_0" : 20}]},
    })

    s = hist.tag.Slicer()
    transforms.update({
        "QCDscale_muRmuFUp" : {
            "action" : lambda h: h[{"muRfact" : 2.j, "muFfact" : 2.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRmuFDown" : {
            "action" : lambda h: h[{"muRfact" : 0.5j, "muFfact" : 0.5j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRUp" : {
            "action" : lambda h: h[{"muRfact" : 2.j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRDown" : {
            "action" : lambda h: h[{"muRfact" : 0.5j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muFUp" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 2.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muFDown" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 0.5j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_cen" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
    })

    def scetlibIdx(h, i):
        return h if not ("vars" in h.axes.name and h.axes["vars"].size > i) else h[{"vars" : i}]

    def projAx(hname):
        return [hname] if hname != "unrolled" else ["pt", "eta"]

    transforms.update({
        "resumFOScaleUp" : {
            "action" : lambda h: scetlibIdx(h, 2)},
        "resumFOScaleDown" : {
            "action" : lambda h: scetlibIdx(h, 1)},
        "resumLambdaDown" : {
            "action" : lambda h: scetlibIdx(h, 3)},
        "resumLambdaUp" : {
            "action" : lambda h: scetlibIdx(h, 4)},
        "resumTransitionMax" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(5,9), no_flow=["ptVgen"], do_min=False)},
        "resumTransitionMin" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(5,9), no_flow=["ptVgen"], do_min=True)},
        "resumScaleMax" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(9,44), no_flow=["ptVgen"], do_min=False)},
        "resumScaleMin" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(9,44), no_flow=["ptVgen"], do_min=True)},
    })
    for k,v in transforms.items():
        if any([x in k for x in ["QCDscale", "resum", "pdf"]]):
            v["procs"] = common.vprocs 
            if any([x in k for x in ["QCDscale", "resum", ]]):
                unc = "qcdScale" if "QCDscale" in k else "scetlibCorr_unc"
                v["hist"] = unc if base_hist == "nominal" else f"{base_hist}_{unc}"

    return transforms

def scale_helicity_hist_to_variations(scale_hist, sum_axes=[], rebinPtV=None):
    
    s = hist.tag.Slicer()
    # select nominal QCD scales, but keep the sliced axis at size 1 for broadcasting
    nom_scale_hist = scale_hist[{"muRfact" : s[1.j:1.j+1], "muFfact" : s[1.j:1.j+1]}]
    axisNames = scale_hist.axes.name
    # select nominal QCD scales and project down to nominal axes
    genAxes = ["ptVgen", "chargeVgen", "helicity"]
    nom_hist = nom_scale_hist[{"muRfact" : s[1.j], "muFfact" : s[1.j] }]
    for genAxis in genAxes:
        if genAxis in axisNames:
            nom_hist = nom_hist[{genAxis : s[::hist.sum]}]
    
    hasHelicityAxis = "helicity" in axisNames
    hasPtAxis = "ptVgen" in axisNames

    if rebinPtV and hasPtAxis:
        # Treat single bin array as a float
        array_rebin = isinstance(rebinPtV, collections.abc.Sequence)
        if array_rebin and len(rebinPtV) == 1:
            rebinPtV = rebinPtV[0]
            array_rebin = False

        if array_rebin:
            scale_hist = hh.rebinHist(scale_hist, "ptVgen", rebinPtV)
            nom_scale_hist = hh.rebinHist(nom_scale_hist, "ptVgen", rebinPtV)
        else:
            scale_hist = scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
            nom_scale_hist = nom_scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
    elif not hasPtAxis:
        raise ValueError("In scale_helicity_hist_to_variations: axis 'ptVgen' not found in histogram.")
            
    for axis in sum_axes:
        if axis in axisNames:
            scale_hist = scale_hist[{axis : s[::hist.sum]}]
            nom_scale_hist = nom_scale_hist[{axis : s[::hist.sum]}]
        else:
            logger.warning(f"In scale_helicity_hist_to_variations: axis '{axis}' not found in histogram.")
        
    # difference between a given scale and the nominal, plus the sum
    # this emulates the "weight if idx else nominal" logic and corresponds to the decorrelated
    # variations
    if scale_hist.name is None:
        out_name = "scale_helicity_variations" if hasHelicityAxis else "scale_vpt_variations" if hasPtAxis else "scale_vcharge_variations"
    else:
        out_name = scale_hist.name + "_variations"

    expd = scale_hist.ndim - nom_hist.ndim
    expandnom = np.expand_dims(nom_hist.view(flow=True), [-expd+i for i in range(expd)])
    systhist = scale_hist.view(flow=True) - nom_scale_hist.view(flow=True) + expandnom

    scale_variation_hist = hist.Hist(*scale_hist.axes, storage = scale_hist._storage_type(), 
                                     name = out_name, data = systhist)

    return scale_variation_hist


def uncertainty_hist_from_envelope(h, proj_ax, entries):
    hdown = hh.syst_min_or_max_env_hist(h, proj_ax, "vars", entries, no_flow=["ptVgen"], do_min=True)
    hup = hh.syst_min_or_max_env_hist(h, proj_ax, "vars", entries, no_flow=["ptVgen"], do_min=False)
    hnew = hist.Hist(*h.axes[:-1], common.down_up_axis, storage=h._storage_type())
    hnew[...,0] = hdown.view(flow=True)
    hnew[...,1] = hup.view(flow=True)
    return hnew

def define_mass_weights(df, isW, nominal_axes=None, nominal_cols=None):
    # nweights = 21 if isW else 23
    # from -100 to 100 MeV with 10 MeV increment
    nweights = 21
    df = df.Define("massWeight_tensor", f"wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight)")

    hist = None
    if nominal_axes and nominal_cols:
        df = df.Define("massWeight_tensor_wnom", "auto res = massWeight_tensor; res = nominal_weight*res; return res;")
        hist = df.HistoBoost("massWeight", nominal_axes, [*nominal_cols, "massWeight_tensor_wnom"])

    return df, hist
