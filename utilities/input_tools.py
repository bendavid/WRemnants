import lz4.frame
import pickle
import numpy
import hist
from utilities import boostHistHelpers as hh
import numpy as np
import logging
import os

def read_and_scale(fname, proc, histname):
    with lz4.frame.open(fname) as f:
        out = pickle.load(f)
        
    return load_and_scale(out, proc, histname)

def load_and_scale(res_dict, proc, histname):
    h = res_dict[proc]["output"][histname]
    scale = res_dict[proc]["dataset"]["xsec"]/res_dict[proc]["weight_sum"]
    return h*scale

def read_all_and_scale(fname, procs, histnames):
    with lz4.frame.open(fname) as f:
        out = pickle.load(f)

    hists = []
    for histname in histnames:
        h = load_and_scale(out, procs[0], histname)
        for proc in procs[1:]:
            h += load_and_scale(out, proc, histname)
        hists.append(h)

    return hists

def read_scetlib_hist(path, nonsing="auto", flip_y_sign=False, charge=None):
    if path[-4:] == ".npz":
        f = np.load(path, allow_pickle=True)
    elif path[-4:] == ".pkl":
        with open(path, "rb") as picklefile:
            f = pickle.load(picklefile)
    else:
        ValueError("File {path} is not a recognized file format")

    var_axis = hist.axis.Integer(f["bins"][0][0], f["bins"][0][-1], name="vars", flow=False)
    # Won't actually have overflow/underflow, but set to match MiNNLO
    mass_underflow = f["bins"][1][0] > 0.
    mass_overflow = f["bins"][1][-1] < 13000.
    mass_axis = hist.axis.Variable(f["bins"][1], name="mass", overflow=mass_overflow, underflow=mass_underflow)
    y_axis = hist.axis.Variable(f["bins"][2], name="y")
    
    # Use 0.1 here rather than 0, because the nonsingular behaves much better with a "cut" at > 0.1
    pt_underflow = f["bins"][3][0] > 0.1
    pt_axis = hist.axis.Variable(f["bins"][3], name="pt", underflow=pt_underflow)

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

    if nonsing and nonsing != "skip":
        if nonsing == "auto":
            nonsing = path.replace(*((".", "_nons.") if "sing" not in path else ("sing", "nons")))
        nonsingh = read_scetlib_hist(nonsing, nonsing="skip", flip_y_sign=flip_y_sign, charge=charge)
        scetlibh = hh.addHists(scetlibh, nonsingh)
    elif nonsing != "skip":
        logging.warning("Will not include nonsingular contribution!")
    
    if flip_y_sign:
        mid = y_axis.index(0)
        s = hist.tag.Slicer()
        scetlibh[{"y" : s[mid:]}] = scetlibh[{"y" : s[mid:]}].view()*-1

    return scetlibh 

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
    
    ax = hist.axis.Variable(bins, name=axname, underflow=not (bins[0] == 0 and "pt" in axname))
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
    
    ax = hist.axis.Variable(bins, name=axname, underflow=not (bins[0] == 0 and "pt" in axname))
    h = hist.Hist(ax, storage=hist.storage.Weight())

    h[...] = data[:-1,2:4]
    return h*1/1000
