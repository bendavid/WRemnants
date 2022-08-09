import os
import ROOT
import pathlib
import hist
import narf
import numpy as np
import lz4.frame
import pickle
from .correctionsTensor_helper import makeCorrectionsTensor
from wremnants import boostHistHelpers as hh,theory_tools
from wremnants.common import data_dir

def make_corr_helper_fromnp(filename=f"{data_dir}/N3LLCorrections/inclusive_{{process}}_pT.npz", isW=True):
    if isW:
        corrf_Wp = np.load(filename.format(process="Wp"), allow_pickle=True)
        corrf_Wm = np.load(filename.format(process="Wm"), allow_pickle=True)
        bins = corrf_Wp["bins"]
        axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
    else:
        corrf = np.load(filename.format(process="Z"), allow_pickle=True)
        bins = corrf["bins"]
        axis_charge = hist.axis.Regular(1, -1., 1., underflow=False, overflow=False, name = "charge")

    axis_syst = hist.axis.Regular(len(bins[0]) - 1, bins[0][0], bins[0][-1], 
                    name="systIdx", overflow=False, underflow=False)
    axis_mass = hist.axis.Variable(bins[1], name="mass")
    axis_y = hist.axis.Variable(bins[2], name="y")
    axis_pt = hist.axis.Regular(len(bins[-1]) - 1, bins[-1][0], bins[-1][-1], name="pT", underflow=False)

    corrh = hist.Hist(axis_mass, axis_y, axis_pt, axis_charge, axis_syst)
    if isW:
        corrh[...,1,:] = np.moveaxis(corrf_Wp["scetlibCorr3D_Wp"], 0, -1)
        corrh[...,0,:] = np.moveaxis(corrf_Wm["scetlibCorr3D_Wm"], 0, -1)
    else:
        corrh[...,0,:] = np.moveaxis(corrf["scetlibCorr3D_Z"], 0, -1)
    corrh[hist.underflow,...] = 1.
    corrh[hist.overflow,...] = 1.
    corrh[:,hist.underflow,...] = 1.
    corrh[:,hist.overflow,...] = 1.
    corrh[:,:,hist.overflow,...] = 1.

    return makeCorrectionsTensor(corrh, ROOT.wrem.TensorCorrectionsHelper, tensor_rank=1)

def make_corr_helper(filename, proc, histname):
    with lz4.frame.open(filename) as f:
        corr = pickle.load(f)
        corrh = corr[proc][histname]

    return makeCorrectionsTensor(corrh, ROOT.wrem.TensorCorrectionsHelper, tensor_rank=1)

def rebin_corr_hists(hists, ndim=-1):
    # Allow trailing dimensions to be different (e.g., variations)
    ndims = min([x.ndim for x in hists]) if ndim < 0 else ndim
    for i in range(ndims):
        # This is a workaround for now for the fact that MiNNLO has mass binning up to
        # Inf whereas SCETlib has 13 TeV
        if all([h.axes[i].size == 1 for h in hists]):
            continue
        hists = hh.rebinHistsToCommon(hists, i)
    return hists

# Assuming the 3 physics variable dimensions are first
def set_corr_ratio_flow(corrh):
    corrh[hist.underflow,...] = np.ones_like(corrh[0,...])
    corrh[hist.overflow,...] = np.ones_like(corrh[0,...])
    corrh[:,hist.underflow,...] = np.ones_like(corrh[:,0,...])
    corrh[:,hist.overflow,...] = np.ones_like(corrh[:,0,...])
    corrh[:,:,hist.overflow,...] = np.ones_like(corrh[:,:,0,...])
    return corrh

def make_corr_from_ratio(denom_hist, num_hist):
    denom_hist, num_hist = rebin_corr_hists([denom_hist, num_hist])

    corrh = hh.divideHists(num_hist, denom_hist)
    return set_corr_ratio_flow(corrh)

def make_corr_by_helicity(ref_helicity_hist, target_sigmaul, target_sigma4, ndim=3):
    ref_helicity_hist, target_sigmaul, target_sigma4 = rebin_corr_hists([ref_helicity_hist, target_sigmaul, target_sigma4], ndim)
    
    ref_coeffs = theory_tools.moments_to_angular_coeffs(ref_helicity_hist)


    target_a4_coeff = make_a4_coeff(target_sigma4, target_sigmaul)

    sigmaUL_ratio = hh.divideHists(target_sigmaul, ref_helicity_hist[{"helicity" : -1.j}])
    # This is because I haven't run W+ yet

    corr_ax = hist.axis.Boolean(name="corr")
    corr_coeffs = hist.Hist(*ref_coeffs.axes, corr_ax)
    # Corr = False is the uncorrected coeffs, corrected coeffs are scaled by the new sigma UL
    # and have the new A4
    corr_coeffs[...,False] = ref_coeffs.values(flow=True)
    corr_coeffs[...,True] = hh.multiplyHists(ref_coeffs, sigmaUL_ratio).values(flow=True)

    corr_coeffs[...,4.j,True] = target_a4_coeff[{"vars" : 0}].values()
    if corr_coeffs.axes["chargeVgen"].size == 2:
        corr_coeffs[...,1.j,:,:].view(flow=True)[...] = corr_coeffs[{"chargeVgen" : -1.j}].view(flow=True)
    return corr_coeffs

def read_scetlib_hist(path, pt_axis=None, nonsing="auto", flip_y_sign=False, charge=None):
    f = np.load(path, allow_pickle=True)
    var_axis = hist.axis.Integer(f["bins"][0][0], f["bins"][0][-1], name="vars", flow=False)
    mass_axis = hist.axis.Variable(f["bins"][1], name="mass")
    y_axis = hist.axis.Variable(f["bins"][2], name="y")
    
    if not pt_axis:
        pt_axis = hist.axis.Variable(f["bins"][3], name="pt")

    h = f["hist"]
    storage = hist.storage.Double()
    axes = [mass_axis,y_axis,pt_axis,var_axis]
    varax_idx = -1 
    vals = np.moveaxis(h, 0, varax_idx)

    if "hist_err" in f:
        err = f["hist_err"]
        storage = hist.storage.Weight()
        vals = np.stack((vals, np.moveaxis(err, 0, varax_idx)), axis=-1)

    if charge is not None:
        charge_args = (2, -2., 2.) if charge != 0 else (1, 0, 1) 
        charge_axis = hist.axis.Regular(*charge_args, flow=False, name = "charge")
        axes.insert(-1, charge_axis)
    
    scetlibh = hist.Hist(*axes, storage=storage)
    if charge is None:
        scetlibh[...] = vals
    else:
        scetlibh[...,charge_axis.index(charge),:] = vals

    if nonsing:
        if nonsing == "auto":
            nonsing = path.replace(".npz", "_nons.npz")
        nonsing = read_scetlib_hist(nonsing, pt_axis, False, 
                        flip_y_sign=flip_y_sign, charge=charge)
        scetlibh = scetlibh + nonsing
    
    if flip_y_sign:
        mid = y_axis.index(0)
        s = hist.tag.Slicer()
        scetlibh[{"y" : s[mid:]}] = scetlibh[{"y" : s[mid:]}].view()*-1

    return scetlibh 

def make_a4_coeff(sigma4_hist, ul_hist):
    return hh.divideHists(sigma4_hist, ul_hist, cutoff=0.0001)

def read_dyturbo_hist(filenames, path="", axis="pt"):
    isfile = list(filter(lambda x: os.path.isfile(x), ["/".join([path, f]) if path else f for f in filenames]))

    if not isfile:
        raise ValueError("Must pass in a valid file")

    hists = [read_dyturbo_file(f, axis) for f in isfile]
    return hh.sumHists(hists)

# Ignoring the scale unc for now
def read_matrixRadish_hist(filename, axname="pt"):
    data = read_text_data(filename)
    bins = list(set(data[:,0].flatten()))
    
    ax = hist.axis.Variable(bins, name=axname)
    h = hist.Hist(ax, storage=hist.storage.Weight())

    h[...] = data[:-1,1:3]
    return h*1/1000
    
def read_text_data(filename):
    data = []
    for line in open(filename).readlines():
        entry = line.split("#")[0]
        entry_data = [float(i.strip()) for i in entry.split()]
        if not entry_data:
            continue
        data.append(entry_data)
    return np.array(data, dtype=float)

def read_dyturbo_file(filename, axname="pt"):
    data = read_text_data(filename)
    # Last line is the total cross section
    bins = list(set(data[:-1,:2].flatten()))
    
    ax = hist.axis.Variable(bins, name=axname)
    h = hist.Hist(ax, storage=hist.storage.Weight())

    h[...] = data[:-1,2:4]
    return h*1/1000
