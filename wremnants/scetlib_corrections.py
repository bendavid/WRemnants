import ROOT
import pathlib
import hist
import narf
import numpy as np
from .correctionsTensor_helper import makeCorrectionsTensor

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def makeScetlibCorrHelper(filename=f"{data_dir}/N3LLCorrections/inclusive_{{process}}_pT.npz"):
    corrf_Wp = np.load(filename.format(process="Wp"), allow_pickle=True)
    corrf_Wm = np.load(filename.format(process="Wm"), allow_pickle=True)
    bins = corrf_Wp["bins"]
    axis_syst = hist.axis.Regular(len(bins[0]) - 1, bins[0][0], bins[0][-1], 
                    name="systIdx", overflow=False, underflow=False)
    axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
    axis_mass = hist.axis.Variable(bins[1], name="mass")
    axis_y = hist.axis.Variable(bins[2], name="y")
    axis_pt = hist.axis.Regular(len(bins[-1]) - 1, bins[-1][0], bins[-1][-1], name="pT", underflow=False)

    corrh = hist.Hist(axis_charge, axis_mass, axis_y, axis_pt, axis_syst)
    corrh[1,...] = np.moveaxis(corrf_Wp["scetlibCorr3D_Wp"], 0, -1)
    corrh[0,...] = np.moveaxis(corrf_Wm["scetlibCorr3D_Wm"], 0, -1)
    corrh[:,hist.underflow,...] = 1.
    corrh[:,hist.overflow,...] = 1.
    corrh[:,:,hist.underflow,...] = 1.
    corrh[:,:,hist.overflow,...] = 1.
    corrh[...,hist.overflow,:] = 1.

    return makeCorrectionsTensor(corrh, ROOT.wrem.TensorCorrectionsHelper, tensor_rank=1)

def readScetlibHist(path, pt_axis=None, add_nonsing=True, flip_y_sign=False):
    f = np.load(path, allow_pickle=True)
    var_axis = hist.axis.Integer(f["bins"][0][0], f["bins"][0][-1], name="vars")
    mass_axis = hist.axis.Variable(f["bins"][1], name="mass")
    y_axis = hist.axis.Variable(f["bins"][2], name="y")
    
    if not pt_axis:
        pt_axis = hist.axis.Variable(f["bins"][3], name="pt")    

    h = f["hist"]
    storage = hist.storage.Double()
    vals = h
    if "hist_err" in f:
        err = f["hist_err"]
        storage = hist.storage.Weight()
        vals = np.stack((h, err), axis=-1)

    scetlibh = hist.Hist(var_axis,mass_axis,y_axis,pt_axis, storage=storage)
    scetlibh[...,:h.shape[-1]] = vals

    if add_nonsing:
        nonsing = readScetlibHist(path.replace(".npz", "_nons.npz"), pt_axis, False)
        scetlibh = scetlibh + nonsing

    if flip_y_sign:
        mid = y_axis.index(0)
        s = hist.tag.Slicer()
        scetlibh[{"y" : s[mid:]}] = scetlibh[{"y" : s[mid:]}].view()*-1

    return scetlibh 

