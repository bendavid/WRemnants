import lz4.frame
import pickle
import hist
from utilities import boostHistHelpers as hh,logging
import numpy as np
import os
import json
import hdf5plugin
import h5py
from narf import ioutils
import ROOT
import uproot

logger = logging.child_logger(__name__)

def read_and_scale_pkllz4(fname, proc, histname, calculate_lumi=False, scale=1):
    with lz4.frame.open(fname) as f:
        results = pickle.load(f)
        
    return load_and_scale(results, proc, histname, calculate_lumi, scale)

def read_hist_names(fname, proc):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
        if proc not in results:
            raise ValueError(f"Invalid process {proc}! No output found in file {fname}")
        return results[proc]["output"].keys()

def read_keys(fname):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
        return results.keys()

def read_xsec(fname, proc):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
        return results[proc]["dataset"]["xsec"]

def read_sumw(fname, proc):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
        return results[proc]["weight_sum"]

def read_and_scale(fname, proc, histname, calculate_lumi=False, scale=1, apply_xsec=True):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
            
        return load_and_scale(results, proc, histname, calculate_lumi, scale, apply_xsec)

def load_and_scale(res_dict, proc, histname, calculate_lumi=False, scale=1., apply_xsec=True):
    h = res_dict[proc]["output"][histname]
    if isinstance(h, ioutils.H5PickleProxy):
        h = h.get()
    if not res_dict[proc]["dataset"]["is_data"]:
        if apply_xsec:
            scale = res_dict[proc]["dataset"]["xsec"]/res_dict[proc]["weight_sum"]*scale
        if calculate_lumi:
            data_keys = [p for p in res_dict.keys() if "dataset" in res_dict[p] and res_dict[p]["dataset"]["is_data"]]
            lumi = sum([res_dict[p]["lumi"] for p in data_keys])*1000
            if not lumi:
                logger.warning("Did not find a data hist! Skipping calculate_lumi option")
                lumi = 1
            scale *= lumi
    return h*scale

def read_all_and_scale(fname, procs, histnames, lumi=False):
    h5file = h5py.File(fname, "r")
    results = ioutils.pickle_load_h5py(h5file["results"])

    hists = []
    for histname in histnames:
        h = load_and_scale(results, procs[0], histname, lumi)
        for proc in procs[1:]:
            h += load_and_scale(results, proc, histname, lumi)
        hists.append(h)

    return hists

def read_scetlib_hist(path, nonsing="none", flip_y_sign=False, charge=None):
    if path[-4:] == ".npz":
        f = np.load(path, allow_pickle=True)
    elif path[-4:] == ".pkl":
        with open(path, "rb") as picklefile:
            f = pickle.load(picklefile)
    else:
        ValueError("File {path} is not a recognized file format")

    if type(f["hist"]) == hist.Hist:
        scetlibh = f["hist"]
    else:
        var_axis = hist.axis.Integer(f["bins"][0][0], f["bins"][0][-1], name="vars", flow=False)
        mass_axis = hist.axis.Variable(f["bins"][1], name="Q", flow=False)
        y_axis = hist.axis.Variable(f["bins"][2], name="Y", flow=False)
        
        # Use 0.1 here rather than 0, because the nonsingular behaves much better with a "cut" at > 0.1
        pt_underflow = f["bins"][3][0] > 0.1
        pt_axis = hist.axis.Variable(f["bins"][3], name="qT", flow=False)

        h = f["hist"]
        storage = hist.storage.Double()
        axes = [mass_axis,y_axis,pt_axis,var_axis]

        varax_idx = -1 
        vals = np.moveaxis(h, 0, varax_idx)

        if "hist_err" in f:
            err = f["hist_err"]
            storage = hist.storage.Weight()
            vals = np.stack((vals, np.moveaxis(err, 0, varax_idx)), axis=-1)

        scetlibh = hist.Hist(*axes, storage=storage, data=vals)

    if charge is not None:
        scetlibh = add_charge_axis(scetlibh, charge)

    if nonsing and nonsing != "none":
        if nonsing == "auto":
            nonsing = path.replace(*((".", "_nons.") if "sing" not in path else ("sing", "nons")))
        nonsingh = read_scetlib_hist(nonsing, nonsing="none", flip_y_sign=flip_y_sign, charge=charge)
        # The overflow in the categorical axis breaks the broadcast
        # FIXME: Only set central for variations that aren't present
        if "vars" in nonsingh.axes.name and nonsingh.axes["vars"].size == 1:
            nonsingh = nonsingh[{"vars" : 0}]
        scetlibh = hh.addHists(scetlibh, nonsingh)
        logger.warning("Adding NLO nonsingular contribution!")
    elif nonsing != "none":
        logger.warning("Will not include nonsingular contribution!")
    
    if flip_y_sign:
        scetlibh = flip_hist_y_sign(scetlibh)

    return scetlibh 

def flip_hist_y_sign(h, yaxis="Y"):
    centers = h.axes[yaxis].centers
    scale = np.ones_like(centers)
    scale[centers < 0] *= -1
    h.values()[...] = h.values()*scale[(None if ax.name != yaxis else slice(0, ax.size) for ax in h.axes)]
    return h 

def read_dyturbo_pdf_hist(base_name, pdf_members, axes, charge=None):
    pdf_ax = hist.axis.StrCategory([f"pdf{i}" for i in range(pdf_members)], name="vars")
    pdf_hist = None
    
    for i in range(pdf_members):
        h = read_dyturbo_hist([base_name.format(i=i)], axes=axes, charge=charge)
        if not pdf_hist:
            pdf_hist = hist.Hist(*h.axes, pdf_ax, storage=h._storage_type())
        pdf_hist[...,i] = h.view(flow=True)
        
    return pdf_hist

def read_dyturbo_hist(filenames, path="", axes=("y", "pt"), charge=None, coeff=None):
    filenames = [os.path.expanduser(os.path.join(path, f)) for f in filenames]
    isfile = list(filter(lambda x: os.path.isfile(x), filenames))

    if not isfile:
        raise ValueError(f"Did not find any valid files in {filenames}")

    hists = [read_dyturbo_file(f, axes, charge, coeff) for f in isfile]
    if len(hists) > 1:
        hists = hh.rebinHistsToCommon(hists, 0)

    h = hh.sumHists(hists)

    if charge is not None and "charge" not in h.axes.name:
        charge_args = (2, -2., 2.) if charge != 0 else (1, 0, 1) 
        charge_axis = hist.axis.Regular(*charge_args, flow=False, name = "charge")
        hnew = hist.Hist(*h.axes, charge_axis, storage=h._storage_type())
        hnew[...,charge_axis.index(charge)] = h.view(flow=True)
        return hnew
    else:
        return h

def expand_dyturbo_filenames(path, basename, varname, pieces=["n3ll_born", "n2ll_ct", "n2lo_vj"], append=None):
    return [os.path.join(path, "_".join(filter(None, [basename, piece, varname, append]))+".txt") for piece in pieces]

def dyturbo_varnames():
    return ["mur{0}_muf{1}_mures{2}".format(i,j,k).replace("0", "H") for i in range(3) for j in range(3) for k in range(3) 
        if abs(i-j) < 2 and abs(i-k) < 2 and abs(j-k) < 2 and not (i == 1 and j == 1 and k == 1)]

def read_dyturbo_variations(path, basename, varnames, axes, pieces=["n3ll_born", "n2ll_ct", "n2lo_vj"], append=None, charge=None):
    central_files = expand_dyturbo_filenames(path, basename, "", pieces, append)
    centralh = read_dyturbo_hist(central_files, axes=axes, charge=charge)
    var_ax = hist.axis.Integer(0, len(varnames)+1, name="vars")
    varh = hist.Hist(*centralh.axes, var_ax, storage=centralh._storage_type())
    varh[...,0] = centralh.view(flow=True)
    for i,var in enumerate(varnames):
        filenames = expand_dyturbo_filenames(path, basename, var, pieces, append)
        varh[...,i+1] = read_dyturbo_hist(filenames, axes=axes, charge=charge).view(flow=True)
    return varh 

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

def read_dyturbo_file(filename, axnames=("Y", "qT"), charge=None, coeff=None):
    if filename.endswith(".root"):
        f = uproot.open(filename)
        hname = "_".join((["wgt", coeff] if coeff else ["s"])+[axnames[0].lower()]) 
        h = f[hname].to_hist()
        if coeff == "a4":
            h = -1*h
    else:
        data = read_text_data(filename)
        # 2 numbers per axis + result + error
        if data.shape[1] != len(axnames)*2+2:
            raise ValueError(f"Mismatch between number of axes advertised ({len(axnames)} ==> {axnames}) and found ({(data.shape[1]-2)/2})")

        axes = []
        offset = True
        for i,name in enumerate(axnames):
            # Normally last line is the total cross section, also possible it isn't, so check the bin ranges
            offset = offset and data[-1,2*i] == data[0,2*i] and data[-1,2*i+1] == data[-2,2*i+1]
            bins = sorted(list(set(data[:len(data)-offset,2*i:2*i+2].flatten())))
            axes.append(hist.axis.Variable(bins, name=name, underflow=not (bins[0] == 0 and "qT" in name)))

        h = hist.Hist(*axes, storage=hist.storage.Weight())
        h[...] = np.reshape(data[:len(data)-offset,len(axes)*2:], (*h.axes.size, 2))

    if charge is not None:
        h = add_charge_axis(h, charge)

    return h*1/1000

def add_charge_axis(h, charge):
    charge_args = (2, -2., 2.) if charge != 0 else (1, 0, 1) 
    charge_axis = hist.axis.Regular(*charge_args, flow=False, name = "charge")

    has_vars = h.axes.name[-1] == "vars"
    new_axes = (*h.axes, charge_axis) if not has_vars else (*h.axes[:-1], charge_axis, h.axes[-1])
    hnew = hist.Hist(*new_axes, storage=h._storage_type())
    if has_vars:
        hnew[...,charge_axis.index(charge),:] = h.view(flow=True)
    else:
        hnew[...,charge_axis.index(charge)] = h.view(flow=True)
    return hnew

def read_matched_scetlib_dyturbo_hist(scetlib_resum, scetlib_fo_sing, dyturbo_fo, axes=None, charge=None, fix_nons_bin0=True, coeff=None):
    hsing = read_scetlib_hist(scetlib_resum, charge=charge, flip_y_sign=coeff=="a4")
    hfo_sing = read_scetlib_hist(scetlib_fo_sing, charge=charge, flip_y_sign=coeff=="a4")

    if axes:
        newaxes = [*axes, "vars"]
        if charge is not None:
            newaxes.insert(-1, "charge")
        hfo_sing = hfo_sing.project(*newaxes)
        hsing = hsing.project(*newaxes)
    if all("pdf" in x for x in hsing.axes["vars"]) and hsing.axes["vars"].size > 1:
        logger.info("Reading PDF variations for DYTurbo")
        pdf_members = hsing.axes["vars"].size
        hfo = read_dyturbo_pdf_hist(dyturbo_fo, pdf_members=pdf_members, axes=axes if axes else hsing.axes.name[:-1], charge=charge)
    else:
        hfo = read_dyturbo_hist([dyturbo_fo], axes=axes if axes else hsing.axes.name[:-1], charge=charge, coeff=coeff)
    for ax in ["Y", "Q"]:
        if ax in set(hfo.axes.name).intersection(set(hfo_sing.axes.name)).intersection(set(hsing.axes.name)):
            hfo, hfo_sing, hsing = hh.rebinHistsToCommon([hfo, hfo_sing, hsing], ax)
    if "vars" in hfo.axes.name and hfo.axes["vars"].size != hfo_sing.axes["vars"].size:
        if hfo.axes["vars"].size == 1:
            hfo = hfo[{"vars" : 0}]
    hnonsing = hh.addHists(-1*hfo_sing, hfo, flow=False, by_ax_name=False)
    if fix_nons_bin0:
        # The 2 is for the WeightedSum
        res = np.zeros((*hnonsing[{"qT" : 0}].shape, 2))
        if "charge" in hnonsing.axes.name:
            hnonsing[...,0,:,:] = res
        else:
            hnonsing[...,0,:] = res
    # TODO: Validate
    if hnonsing.axes["vars"].size != hsing.axes["vars"].size:
        htmp_nonsing = hsing.copy()
        for var in htmp_nonsing.axes["vars"]:
            logger.warning(f"Did not find variation {var} for nonsingular! Assuming nominal")
            htmp_nonsing = hnonsing[{"vars" : var if var in hnonsing.axes["vars"] else 0}]
        hnonsing = htmp_nonsing
    return hh.addHists(hsing, hnonsing, by_ax_name=False)

def read_json(fIn):

    if not os.path.exists(fIn):
        logger.warning(f"File {fIn} not found")
        return False
    else:
        with open(fIn) as f: jsDict = json.load(f)
        return jsDict

def safeGetRootObject(fileObject, objectName, quitOnFail=True, silent=False, detach=True):
    obj = fileObject.Get(objectName)
    if obj == None:
        error_msg = f"Error getting {objectName} from file {fileObject.GetName()}"
        if not silent:
            logger.error(error_msg)
        if quitOnFail:
            raise IOError(error_msg)
        return None
    else:
        if detach:
            obj.SetDirectory(0)
        return obj
        
def safeOpenRootFile(fileName, quitOnFail=True, silent=False, mode="READ"):
    fileObject = ROOT.TFile.Open(fileName, mode)
    if not fileObject or fileObject.IsZombie():
        error_msg = f"Error when opening file {fileName}"
        if not silent:
            logger.error(error_msg)
        if quitOnFail:
            raise IOError(error_msg)
        else:
            return None
    elif not fileObject.IsOpen():
        error_msg = f"File {fileName} was not opened"
        if not silent:
            logger.error(error_msg)
        if quitOnFail:
            raise IOError(error_msg)
        else:
            return None
    else:
        return fileObject

def args_from_metadata(card_tool, arg):
    meta_data = card_tool.datagroups.getMetaInfo()
    if "args" not in meta_data.keys():
        raise IOError(f"The argument {arg} was not found in the metadata, maybe you run on an obsolete file.")
    elif arg not in meta_data["args"].keys():
        raise IOError(f"Did not find the argument {arg} in the meta_data dict. Maybe it is an outdated option")

    return meta_data["args"][arg]

def get_metadata(infile):
    results = None
    if infile.endswith(".pkl.lz4"):
        with lz4.frame.open(infile) as f:
            results = pickle.load(f)
    elif infile.endswith(".pkl"):
        with open(infile, "rb") as f:
            results = pickle.load(f)
    elif infile.endswith(".hdf5"):
        h5file = h5py.File(infile, "r")
        meta = h5file.get("results", h5file.get("meta", None))
        results = ioutils.pickle_load_h5py(meta) if meta else None

    if results is None:
        logger.warning("Failed to find results dict. Note that only pkl, hdf5, and pkl.lz4 file types are supported")
        return None

    return results["meta_info"] if "meta_info" in results else results["meta_data"]

def get_scetlib_config(infile):
    if infile.endswith(".pkl"):
        with open(infile, "rb") as f:
            results = pickle.load(f)
        return results["config"]
    else:
        raise ValueError("Expected scetlib output in pkl format")

def read_infile(input):
    # read histogramer input file(s)
    result = {}
    meta = []
    infiles = []
    if isinstance(input, list):
        for inpt in input:
            r, m, h = read_infile(inpt)
            result.update(r)
            meta += m
            infiles += h
        return result, meta, infiles
    
    logger.info(f"Load {input}")
    if input.endswith(".pkl.lz4"):
        with lz4.frame.open(input) as f:
            result = pickle.load(f)
    elif input.endswith(".hdf5"):
        h5file = h5py.File(input, "r")
        infiles = [h5file]
        result = ioutils.pickle_load_h5py(h5file["results"])
    else:
        raise ValueError("Unsupported file type")

    meta = result["meta_info"] if "meta_info" in result else result["meta_data"]

    return result, [meta], [infiles]

def read_dyturbo_angular_coeffs(dyturbof, boson=None, rebin=None, absy=True, add_axes=[]):
    if add_axes and not all(ax.size == 1 for ax in add_axes):
        raise ValueError("Can only add axes of size 1!")

    if type(dyturbof) == str:
        dyturbof = uproot.open(dyturbof)

    if not boson:
        boson = "Wp" if "wp" in dyturbof.file_path else ("Wm" if "wm" in dyturbof else "Z")

    sigma_ul = dyturbof["s_qt_vs_y"].to_hist()
    for ax,name in zip(sigma_ul.axes, ["qT", "Y"]):
        ax._ax.metadata["name"] = name

    if rebin:
        sigma_ul = hh.rebinHistMultiAx(sigma_ul, rebin)
    if absy:
        sigma_ul = hh.makeAbsHist(sigma_ul, "Y")

    charge_range = (2, -2, 2) if "w" in boson.lower() else (1, -1, 1)
    charge = 0 if boson.lower() == "z" else (-1 if "wm" in boson.lower() else 1)
    charge_ax = hist.axis.Regular(*charge_range, name='charge', flow=False)
    h = hist.Hist(*add_axes, *reversed(sigma_ul.axes), charge_ax, hist.axis.Regular(8, 0, 8, flow=False, name="helicity"), storage=sigma_ul.storage_type())
    for i in range(8):
        sigma_i = dyturbof[f"wgt_a{i}_y_qt"].to_hist()
        for ax,name in zip(sigma_i.axes, ["qT", "Y"]):
            ax._ax.metadata["name"] = name
        if rebin:
            sigma_i = hh.rebinHistMultiAx(sigma_i, rebin)
        if absy:
            sigma_i = hh.makeAbsHist(sigma_i, "Y")
        entry = (*[0]*len(add_axes), Ellipsis, charge_ax.index(charge), i)
        # qT, y order in DY turbo is reversed wrt MiNNLO
        h[entry] = hh.divideHists(sigma_i, sigma_ul).view().T

    return h

def read_mu_hist_combine_tau(minnlof, mu_sample, hist_name):
    hmu = read_and_scale(minnlof, mu_sample, hist_name, apply_xsec=False)
    sumw = read_sumw(minnlof, mu_sample)
    xsec = read_xsec(minnlof, mu_sample)
    
    tau_sample = mu_sample.replace("mu", "tau")
    htau = read_and_scale(minnlof, tau_sample, hist_name, apply_xsec=False)
    sumw += read_sumw(minnlof, tau_sample)
    return (hmu + htau)*xsec/sumw
