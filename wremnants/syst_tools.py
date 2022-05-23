import hist
import numpy as np

def scale_helicity_hist_to_variations(scale_hist, sum_helicity=False, sum_ptV=False, rebinPtV=0):

    s = hist.tag.Slicer()
    # select nominal QCD scales, but keep the sliced axis at size 1 for broadcasting
    nom_scale_hist = scale_hist[{"muRfact" : s[1.j:1.j+1], "muFfact" : s[1.j:1.j+1]}]

    # select nominal QCD scales and project down to nominal axes
    nom_hist = nom_scale_hist[{"ptVgen" : s[::hist.sum], "chargeVgen" : s[::hist.sum], "helicity" : s[::hist.sum], "muRfact" : s[1.j], "muFfact" : s[1.j] }]

    if rebinPtV > 0:
        scale_hist = scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
        nom_scale_hist = nom_scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
    
    if sum_helicity:
        scale_hist = scale_hist[{"helicity" : s[::hist.sum]}]
        nom_scale_hist = nom_scale_hist[{"helicity" : s[::hist.sum]}]

    if sum_ptV:
        scale_hist = scale_hist[{"ptVgen" : s[::hist.sum]}]
        nom_scale_hist = nom_scale_hist[{"ptVgen" : s[::hist.sum]}]

    # difference between a given scale and the nominal, plus the sum
    # this emulates the "weight if idx else nominal" logic and corresponds to the decorrelated
    # variations
    if scale_hist.name is None:
        out_name = "scale_helicity_variations"
    else:
        out_name = scale_hist.name + "_variations"

    expd = scale_hist.ndim - nom_hist.ndim
    expandnom = np.expand_dims(nom_hist.view(flow=True), [-expd+i for i in range(expd)])
    systhist = scale_hist.view(flow=True) - nom_scale_hist.view(flow=True) + expandnom

    scale_variation_hist = hist.Hist(*scale_hist.axes, storage = scale_hist._storage_type(), 
                name = out_name, data = systhist)

    return scale_variation_hist

def define_mass_weights(df, isW, nominal_axes=None, nominal_cols=None):
    #nweights = 21 if isW else 23
    nweights = 21
    df = df.Define("massWeight_tensor", f"wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight)")

    hist = None
    if nominal_axes and nominal_cols:
        df = df.Define("massWeight_tensor_wnom", "auto res = massWeight_tensor; res = nominal_weight*res; return res;")
        hist = df.HistoBoost("massWeight", nominal_axes, [*nominal_cols, "massWeight_tensor_wnom"])

    return df, hist
