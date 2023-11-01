import os
import ROOT
import pathlib
import hist
import narf
import numpy as np
import lz4.frame
import pickle
import re
import glob
from .correctionsTensor_helper import makeCorrectionsTensor
from utilities import boostHistHelpers as hh, common, logging
from utilities.io_tools import input_tools
from wremnants import theory_tools

logger = logging.child_logger(__name__)

def valid_theory_corrections():
    corr_files = glob.glob(common.data_dir+"TheoryCorrections/*Corr*.pkl.lz4")
    matches = [re.match("(^.*)Corr[W|Z]\.pkl\.lz4", os.path.basename(c)) for c in corr_files]
    return [m[1] for m in matches if m]

def load_corr_helpers(procs, generators):
    corr_helpers = {}
    for proc in procs:
        corr_helpers[proc] = {}
        for generator in generators:
            fname = f"{common.data_dir}/TheoryCorrections/{generator}Corr{proc[0]}.pkl.lz4"
            if not os.path.isfile(fname):
                logger.warning(f"Did not find correction file for process {proc}, generator {generator}. No correction will be applied for this process!")
                continue
            helper_func = make_corr_helper if "Helicity" not in generator else make_corr_by_helicity_helper
            # Hack for now
            corr_hist_name = get_corr_name(generator)
            weighted_corr = generator in theory_tools.theory_corr_weight_map
            corr_helpers[proc][generator] = helper_func(fname, proc[0], corr_hist_name, weighted_corr)
    for generator in generators:
        if not any([generator in corr_helpers[proc] for proc in procs]):
            raise ValueError(f"Did not find correction for generator {generator} for any processes!")
    return corr_helpers


def make_corr_helper_fromnp(filename=f"{common.data_dir}/N3LLCorrections/inclusive_{{process}}_pT.npz", isW=True):
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

    return makeCorrectionsTensor(corrh)

def load_corr_hist(filename, proc, histname):
    with lz4.frame.open(filename) as f:
        corr = pickle.load(f)
        corrh = corr[proc][histname]
    return corrh

def make_corr_helper(filename, proc, histname, weighted_corr):
    corrh = load_corr_hist(filename, proc, histname)
    return makeCorrectionsTensor(corrh, weighted_corr=weighted_corr)

def make_corr_by_helicity_helper(filename, proc, histname, weighted_corr):
    corrh = load_corr_hist(filename, proc, histname)
    return makeCorrectionsTensor(corrh, ROOT.wrem.CentralCorrByHelicityHelper, tensor_rank=3)

def get_corr_name(generator):
    # Hack for now
    label = generator.replace("1D", "")
    return f"{label}_minnlo_ratio" if "Helicity" not in generator else f"{label.replace('Helicity', '')}_minnlo_coeffs"

def rebin_corr_hists(hists, ndim=-1, binning=None):
    # Allow trailing dimensions to be different (e.g., variations)
    ndims = min([x.ndim for x in hists]) if ndim < 0 else ndim
    if binning:
        try:
            hists = [h if not h else hh.rebinHistMultiAx(h, binning) for h in hists]
        except ValueError as e:
            logger.warning("Can't rebin axes to predefined binning")
        return hists

    for i in range(ndims):
        # This is a workaround for now for the fact that MiNNLO has mass binning up to
        # Inf whereas SCETlib has 13 TeV
        if all([h.axes[i].size == 1 for h in hists if h]):
            continue
        hists = hh.rebinHistsToCommon(hists, i)
    return hists

# Assuming the 3 physics variable dimensions are first
def set_corr_ratio_flow(corrh):
    # Probably there's a better way to do this...
    if corrh.axes[0].traits.underflow:
        corrh[hist.underflow,...] = np.ones_like(corrh[0,...].view(flow=True))
    if corrh.axes[0].traits.overflow:
        corrh[hist.overflow,...] = np.ones_like(corrh[0,...].view(flow=True))

    if corrh.axes[1].traits.underflow:
        corrh[:,hist.underflow,...] = np.ones_like(corrh[:,0,...].view(flow=True))
    if corrh.axes[1].traits.overflow:
        corrh[:,hist.overflow,...] = np.ones_like(corrh[:,0,...].view(flow=True))

    if corrh.axes[2].traits.underflow:
        corrh[:,:,hist.underflow,...] = np.ones_like(corrh[:,:,0,...].view(flow=True))
    if corrh.axes[2].traits.overflow:
        corrh[:,:,hist.overflow,...] = np.ones_like(corrh[:,:,0,...].view(flow=True))
    return corrh

def make_corr_from_ratio(denom_hist, num_hist, rebin=False):

    denom_hist, num_hist = rebin_corr_hists([denom_hist, num_hist], binning=rebin)


    corrh = hh.divideHists(num_hist, denom_hist, flow=False, by_ax_name=False)
    return set_corr_ratio_flow(corrh), denom_hist, num_hist

def make_corr_by_helicity(ref_helicity_hist, target_sigmaul, target_sigma4, coeff_hist=None, coeffs_from_hist=[], binning=None, ndim=3):
    ref_helicity_hist, target_sigmaul, target_sigma4 = rebin_corr_hists([ref_helicity_hist, target_sigmaul, target_sigma4], ndim, binning)

    apply_coeff_corr = coeff_hist is not None and coeffs_from_hist

    # broadcast back mass and helicity axes (or vars axis)
    sigmaUL_ratio = hh.divideHists(target_sigmaul, ref_helicity_hist[{"massVgen" : 0, "helicity" : -1.j}], 
                                   flow=False, by_ax_name=False).values()[np.newaxis,...,np.newaxis,:] \
                            if target_sigmaul else np.ones_like(ref_helicity_hist)[...,np.newaxis]

    ref_coeffs = theory_tools.moments_to_angular_coeffs(ref_helicity_hist)

    corr_ax = hist.axis.Boolean(name="corr")
    vars_ax = target_sigmaul.axes["vars"] if target_sigmaul else hist.axis.Regular(1 ,0, 1, name="vars")
    corr_coeffs = hist.Hist(*ref_coeffs.axes, corr_ax, vars_ax)
    # Corr = False is the uncorrected coeffs, corrected coeffs have the new A4
    # NOTE: the corrected coeffs are multiplied through by the sigmaUL correction, so that the 
    # new correction can be made as the ratio of the sum. To get the correct coeffs, this should
    # be divided back out
    corr_coeffs[...] = ref_coeffs.values()[...,np.newaxis,np.newaxis]

    if target_sigma4:
        target_a4_coeff = make_angular_coeff(target_sigma4, target_sigmaul)
        corr_coeffs[...,4.j,True,:] = target_a4_coeff.values()

    if apply_coeff_corr:
        for coeff in coeffs_from_hist:
            scale = -1 if coeff in [1, 4] else 1
            idx = complex(0, coeff)
            corr_coeffs[...,idx,True,:] = coeff_hist[{"helicity" : idx}].values()[...,np.newaxis]*scale
    # Scale by the sigmaUL correction, 
    corr_coeffs.values()[...,corr_coeffs.axes["corr"].index(True),:] *= sigmaUL_ratio

    corr_coeffs = set_corr_ratio_flow(corr_coeffs)
    return corr_coeffs

def make_angular_coeff(sigmai_hist, ul_hist):
    return hh.divideHists(sigmai_hist, ul_hist, cutoff=0.0001)

def read_combined_corrs(procNames, generator, corr_files, axes=[], absy=True, rebin={}):
    h = None
    if not corr_files:
        return h

    for procName in procNames:
        if procName[0] == "W":
            proc_files = list(filter(lambda x: procName[:2].lower() in x.lower(), corr_files))
            charge = 1 if procName[:2] == "Wp" else -1
        else:
            charge = 0
            proc_files = corr_files
        hproc = read_corr(generator, proc_files, charge, axes)
        h = hproc if not h else h+hproc

    if absy and "Y" in axes:
        h = hh.makeAbsHist(h, "Y")

    if rebin:
        h = hh.rebinHistMultiAx(h, rebin)

    return h

def read_corr(generator, corr_files, charge, axes=[]):
    if "scetlib" in generator:
        coeff=None
        if any("A4" in c for c in corr_files):
            coeff = "a4"
        if "dyturbo" in generator:
            scetlib_files = [x for x in corr_files if pathlib.Path(x).suffix == ".pkl"]
            if len(scetlib_files) != 2:
                raise ValueError(f"scetlib_dyturbo correction requires two SCETlib files (resummed and FO singular). Found {len(scetlib_files)}")
            if not any("nnlo_sing" in x for x in scetlib_files):
                raise ValueError("Must pass in a fixed order singular file")
            nnlo_sing_idx = 0 if "nnlo_sing" in scetlib_files[0] else 1
            resumf = scetlib_files[~nnlo_sing_idx]
            nnlo_singf = scetlib_files[nnlo_sing_idx]

            dyturbo_files = [x for x in corr_files if pathlib.Path(x).suffix == ".txt"]
            if len(dyturbo_files) != 1:
                raise ValueError("scetlib_dyturbo correction requires one DYTurbo file (fixed order contribution)")

            corrh = input_tools.read_matched_scetlib_dyturbo_hist(resumf, nnlo_singf, dyturbo_files[0], axes, charge=charge, coeff=coeff)
        else:
            corrh = input_tools.read_scetlib_hist(corr_files[0], charge=charge, nonsing=nons, flip_y_sign=coeff=="a4")
    else:
        if generator == "matrix_radish":
            h = input_tools.read_matrixRadish_hist(corr_file[0], axes[0])
        elif generator == "dyturbo":
            h = input_tools.read_dyturbo_hist(corr_files, axes=axes, charge=charge)

        vars_ax = h.axes["vars"] if "vars" in h.axes.name else hist.axis.StrCategory(["central"], name="vars") 
        hnD = hist.Hist(*h.axes, vars_ax)
        # Leave off the overflow, we won't use it anyway
        hnD[...] = np.reshape(h.values(), hnD.shape)
        corrh = hnD

    return corrh
