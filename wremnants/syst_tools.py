import hist
import numpy as np

def scale_helicity_hist_to_variations(scale_hist, sum_axis=[], rebinPtV=0):
    
    s = hist.tag.Slicer()
    # select nominal QCD scales, but keep the sliced axis at size 1 for broadcasting
    nom_scale_hist = scale_hist[{"muRfact" : s[1.j:1.j+1], "muFfact" : s[1.j:1.j+1]}]
    axisNames = [ax.name for ax in scale_hist.axes]
    hasHelicityAxis = "helicity" in axisNames
    # select nominal QCD scales and project down to nominal axes
    if hasHelicityAxis:
        nom_hist = nom_scale_hist[{"ptVgen" : s[::hist.sum], "chargeVgen" : s[::hist.sum], "helicity" : s[::hist.sum], "muRfact" : s[1.j], "muFfact" : s[1.j] }]
    else:
        nom_hist = nom_scale_hist[{"ptVgen" : s[::hist.sum], "chargeVgen" : s[::hist.sum], "muRfact" : s[1.j], "muFfact" : s[1.j] }]

    if rebinPtV > 0:
        if "ptVgen" in axisNames:
            scale_hist = scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
            nom_scale_hist = nom_scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
        else:
            raise ValueError("In scale_helicity_hist_to_variations: axis 'ptVgen' not found in histogram.")
            
    for axis in sum_axis:
        if axis in axisNames:
            scale_hist = scale_hist[{axis : s[::hist.sum]}]
            nom_scale_hist = nom_scale_hist[{axis : s[::hist.sum]}]
        else:
            raise ValueError(f"In scale_helicity_hist_to_variations: axis '{axis}' not found in histogram.")
        
    # difference between a given scale and the nominal, plus the sum
    # this emulates the "weight if idx else nominal" logic and corresponds to the decorrelated
    # variations
    if scale_hist.name is None:
        out_name = "scale_helicity_variations" if hasHelicityAxis else "scale_vpt_variations"
    else:
        out_name = scale_hist.name + "_variations"

    expd = scale_hist.ndim - nom_hist.ndim
    expandnom = np.expand_dims(nom_hist.view(flow=True), [-expd+i for i in range(expd)])
    systhist = scale_hist.view(flow=True) - nom_scale_hist.view(flow=True) + expandnom

    scale_variation_hist = hist.Hist(*scale_hist.axes, storage = scale_hist._storage_type(), 
                name = out_name, data = systhist)

    return scale_variation_hist


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
