import lz4.frame
import pickle
import hist
from utilities import boostHistHelpers as hh,logging
import numpy as np
import os
import hdf5plugin
import h5py
from narf import ioutils

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

def read_and_scale(fname, proc, histname, calculate_lumi=False, scale=1):
    with h5py.File(fname, "r") as h5file:
        results = ioutils.pickle_load_h5py(h5file["results"])
            
        return load_and_scale(results, proc, histname, calculate_lumi, scale)

def load_and_scale(res_dict, proc, histname, calculate_lumi=False, scale=1.):
    h = res_dict[proc]["output"][histname]
    if isinstance(h, ioutils.H5PickleProxy):
        h = h.get()
    if not res_dict[proc]["dataset"]["is_data"]:
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
        mid = y_axis.index(0)
        s = hist.tag.Slicer()
        scetlibh[{"Y" : s[mid:]}] = scetlibh[{"Y" : s[mid:]}].view()*-1

    return scetlibh 

def read_dyturbo_pdf_hist(base_name, pdf_members, axes, charge=None):
    pdf_ax = hist.axis.StrCategory([f"pdf{i}" for i in range(pdf_members)], name="vars")
    pdf_hist = None
    
    for i in range(pdf_members):
        h = read_dyturbo_hist([base_name.format(i=i)], axes=axes, charge=charge)
        if not pdf_hist:
            pdf_hist = hist.Hist(*h.axes, pdf_ax, storage=h._storage_type())
        pdf_hist[...,i] = h.view(flow=True)
        
    return pdf_hist

def read_dyturbo_hist(filenames, path="", axes=("y", "pt"), charge=None):
    isfile = list(filter(lambda x: os.path.isfile(x), 
        [os.path.expanduser(os.path.join(path, f)) for f in filenames]))

    if not isfile:
        raise ValueError("Must pass in a valid file")

    hists = [read_dyturbo_file(f, axes, charge) for f in isfile]
    if len(hists) > 1:
        hists = hh.rebinHistsToCommon(hists, 0)

    h = hh.sumHists(hists)

    if charge is not None:
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

def read_dyturbo_file(filename, axnames=("Y", "qT"), charge=None):
    charge=None
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

def readImpacts(rtfile, group, sort=True, add_total=True, stat=0.0):
    histname = "nuisance_group_impact_nois" if group else "nuisance_impact_nois"
    impacts = rtfile[histname].to_hist()
    labels = np.array([impacts.axes[1].value(i) for i in range(impacts.axes[1].size)])
    total = rtfile["fitresults"][impacts.axes[0].value(0)+"_err"].array()[0]
    impacts = impacts.values()[0,:]
    if sort:
        order = np.argsort(impacts)
        impacts = impacts[order]
        labels = labels[order]
    if add_total:
        impacts = np.append(impacts, total)
        labels = np.append(labels, "Total")

    if stat > 0:
        idx = np.argwhere(labels == "stat")
        impacts[idx] = stat

    return impacts,labels

def read_matched_scetlib_dyturbo_hist(scetlib_resum, scetlib_fo_sing, dyturbo_fo, axes=None, charge=None, fix_nons_bin0=True):
    hsing = read_scetlib_hist(scetlib_resum, charge=charge)
    hfo_sing = read_scetlib_hist(scetlib_fo_sing, charge=charge)
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
        hfo = read_dyturbo_hist([dyturbo_fo], axes=axes if axes else hsing.axes.name[:-1], charge=charge)
    for ax in ["Y", "Q"]:
        if ax in set(hfo.axes.name).intersection(set(hfo_sing.axes.name)).intersection(set(hsing.axes.name)):
            hfo, hfo_sing, hsing = hh.rebinHistsToCommon([hfo, hfo_sing, hsing], ax)
    if "vars" in hfo.axes.name and hfo.axes["vars"].size != hfo_sing.axes["vars"].size:
        if hfo.axes["vars"].size == 1:
            hfo = hfo[{"vars" : 0}]
    hnonsing = hh.addHists(-1*hfo_sing, hfo)
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
    return hh.addHists(hsing, hnonsing)
