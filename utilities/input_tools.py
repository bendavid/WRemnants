import lz4.frame
import pickle
import numpy
import hist
from utilities import boostHistHelpers as hh
import numpy as np
import logging
import os

def read_and_scale(fname, proc, histname, lumi):
    with lz4.frame.open(fname) as f:
        out = pickle.load(f)
        
    return load_and_scale(out, proc, histname, lumi)

def load_and_scale(res_dict, proc, histname, lumi=False):
    h = res_dict[proc]["output"][histname]
    scale = res_dict[proc]["dataset"]["xsec"]/res_dict[proc]["weight_sum"]
    if lumi:
        data_keys = [p for p in res_dict.keys() if "dataset" in res_dict[p] and res_dict[p]["dataset"]["is_data"]]
        lumi = sum([res_dict[p]["lumi"] for p in data_keys])*1000
        scale *= lumi
    return h*scale

def read_all_and_scale(fname, procs, histnames, lumi=False):
    with lz4.frame.open(fname) as f:
        out = pickle.load(f)

    hists = []
    for histname in histnames:
        h = load_and_scale(out, procs[0], histname, lumi)
        for proc in procs[1:]:
            h += load_and_scale(out, proc, histname, lumi)
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

def read_dyturbo_hist(filenames, path="", axes=("y", "pt")):
    print("Files", filenames)
    isfile = list(filter(lambda x: os.path.isfile(x), 
        [os.path.join(f, os.path.expanduser(f)) for f in filenames]))
    print("Valid Files", isfile)

    if not isfile:
        raise ValueError("Must pass in a valid file")


    hists = [read_dyturbo_file(f, axes) for f in isfile]
    if len(hists) > 1:
        print([h.sum() for h in hists])
        hists = hh.rebinHistsToCommon(hists, 0)

    return hh.sumHists(hists)

def distribution_to_hist(data):
    next_bin = data[1:,0]
    bin_width = next_bin-data[:-1,0]
    data[:-1,1:] = data[:-1,1:]*bin_width[:,np.newaxis]
    return data

# Ignoring the scale unc for now
def read_matrixRadish_hist(filename, axname="pt"):
    data = read_text_data(filename)
    # Multiply through by bin width
    data = distribution_to_hist(data)
    bins = np.unique(data[:,0])
    
    ax = hist.axis.Variable(bins, name=axname, underflow=not (bins[0] == 0 and "pt" in axname))
    var_ax = hist.axis.Integer(0, 3, name="vars", flow=False)
    h = hist.Hist(ax, var_ax, storage=hist.storage.Double())

    h[...] = data[:-1, np.array([1,3,5])]
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

def read_dyturbo_file(filename, axnames=("y", "pt")):
    data = read_text_data(filename)
    # 2 numbers per axis + result + error
    if data.shape[1] != len(axnames)*2+2:
        raise ValueError(f"Mismatch between number of axes advertised and found ({(data.shape[1]-2)/2} vs. {len(axnames)})")

    axes = []
    offset = True
    for i,name in enumerate(axnames):
        # Normally last line is the total cross section, also possible it isn't, so check the bin ranges
        offset = offset and data[-1,2*i] == data[0,2*i] and data[-1,2*i+1] == data[-2,2*i+1]
        bins = sorted(list(set(data[:len(data)-offset,2*i:2*i+2].flatten())))
        axes.append(hist.axis.Variable(bins, name=name, underflow=not (bins[0] == 0 and "pt" in name)))

    h = hist.Hist(*axes, storage=hist.storage.Weight())

    h[...] = np.reshape(data[:len(data)-offset,len(axes)*2:], (*h.axes.size, 2))
    return h*1/1000
