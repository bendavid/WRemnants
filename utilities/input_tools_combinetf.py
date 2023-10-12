import hist
from utilities import logging
import numpy as np
import os
import json
import hdf5plugin
import h5py
from narf import ioutils
import ROOT
import uproot
import pandas as pd

logger = logging.child_logger(__name__)

def get_fitresult(fitresult_filename):
    if fitresult_filename.endswith(".root"):
        return uproot.open(fitresult_filename)
    elif fitresult_filename.endswith(".hdf5"):
        return h5py.File(fitresult_filename, mode='r')
    else:
        raise IOError(f"Unknown format of fitresult {fitresult}")

def is_root_file(fileobject):
    return isinstance(fileobject, uproot.ReadOnlyDirectory)

def is_h5_file(fileobject):
    return isinstance(fileobject, h5py.File)

def get_poi_names(fitresult_file, poi_type="mu"):
    if is_h5_file(fitresult_file):
        names = get_poi_names_h5(fitresult_file, poi_type)
    elif is_root_file(fitresult_file):
        names = get_poi_names_root(fitresult_file, poi_type)
    else:
        raise IOError(f"Unknown format of fitresult {fitresult}")

    if len(names)==0:
        logger.warning('No free parameters found (neither signal strenght(s), nor W mass)')
        return [None]
    return names

def get_poi_names_h5(h5file, poi_type="mu"):
    outnames = h5file["outnames"][...].astype(str)
    names = np.array([])
    if poi_type is not None and poi_type in outnames:
        names = h5file[f"{poi_type}_names"][...].astype(str)

    if 'nois' in outnames:
        names = np.append(names, 'Wmass')

    return names

def get_poi_names_root(rtfile, poi_type="mu"):
    names = []
    if poi_type is not None and f'nuisance_impact_{poi_type}' in [k.replace(";1","") for k in rtfile.keys()]:
        impacts = rtfile[f'nuisance_impact_{poi_type}'].to_hist()
        names = [impacts.axes[0].value(i) for i in range(impacts.axes[0].size)]

    if 'nuisance_impact_nois' in [k.replace(";1","") for k in rtfile.keys()]:
        impacts = rtfile['nuisance_impact_nois'].to_hist()
        names.append('Wmass')

    return np.array(names)

def read_impacts_poi(fileobject, group, sort=True, add_total=True, stat=0.0, poi='Wmass', normalize=True):
    # read impacts of a single POI
    if is_h5_file(fileobject):
        impacts, labels, norm, total = read_impacts_poi_h5(fileobject, group, poi=poi)
    elif is_root_file(fileobject):
        impacts, labels, norm, total = read_impacts_poi_root(fileobject, group, poi=poi)
    else:
        raise IOError(f"Unknown format of fitresult {fitresult}")

    if sort:
        order = np.argsort(impacts)
        impacts = impacts[order]
        labels = labels[order]

    if add_total:
        impacts = np.append(impacts, total)
        labels = np.append(labels, "Total")

    if normalize:
        impacts /= norm

    if stat > 0:
        idx = np.argwhere(labels == "stat")
        impacts[idx] = stat

    return impacts, labels, norm

def read_impacts_poi_h5(h5file, group, poi='Wmass', skip_systNoConstraint=True):
    if poi is None:
        poi_type=None
    else:
        poi_type = poi.split("_")[-1] if poi else None
        poi_names = get_poi_names(h5file, poi_type)
        if poi=='Wmass':
            impact_hist_total = "nuisance_impact_nois"
            impact_hist = "nuisance_group_impact_nois" if group else "nuisance_impact_nois"
        elif poi in poi_names:
            impact_hist_total = f"nuisance_impact_{poi_type}"
            impact_hist = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"
        else:
            raise ValueError(f"Invalid POI: {poi}")

    if group:
        labels = h5file["hsystgroups"][...].astype(str)
        labels = np.append(labels, "stat")
    if not group:
        labels = h5file["hsysts"][...].astype(str)

    if poi_type is None:
        impacts = np.zeros_like(labels, dtype=int)
        total = 0.
        norm = 0.
    else:
        ipoi = 0 if poi=='Wmass' else poi_names.index(poi)
        impacts = h5file[impact_hist][...][ipoi]
        total = h5file[impact_hist_total][...][ipoi,ipoi]
        norm = h5file["x"][...][ipoi]

    if len(labels)+1 == len(impacts): 
        labels = np.append(labels, "binByBinStat")

    if skip_systNoConstraint:
        noConstraint = np.array([l in h5file["hsystsnoconstraint"][...].astype(str) for l in labels])
        labels = labels[~noConstraint]
        impacts = impacts[~noConstraint]

    return impacts, labels, norm, total

def read_impacts_poi_root(rtfile, group, poi='Wmass'):
    poi_type = poi.split("_")[-1] if poi else None
    poi_names = get_poi_names(rtfile, poi_type)
    if poi=='Wmass':
        impact_hist = "nuisance_group_impact_nois" if group else "nuisance_impact_nois"
    elif poi in poi_names:
        impact_hist = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"
    else:
        raise ValueError(f"Invalid POI: {poi}")

    if group:
        histname = impact_hist
    else:
        histname = "correlation_matrix_channelmu"
        
    h = rtfile[histname].to_hist()
    labels = np.array(list(h.axes["yaxis"]), dtype=object)

    if impact_hist not in rtfile:
        logger.warning("Did not find impact hist in file. Skipping!")
        return np.zeros_like(labels), labels, 1., 1.
    else:
        impacts = rtfile[impact_hist].to_hist()
        ipoi = 0 if poi=='Wmass' else poi_names.index(poi)
        total = rtfile["fitresults"][impacts.axes[0].value(ipoi)+"_err"].array()[0]
        norm = rtfile["fitresults"][impacts.axes[0].value(ipoi)].array()[0]
        impacts = impacts.values()[ipoi,:]

    return impacts, labels, norm, total

def decode_poi_bin(name, var):
    name_split = name.split(var)
    if len(name_split) == 1:
        return None
    else:
        return name_split[-1].split("_")[0]

def filter_poi_bins(names, gen_axes, selections={}, base_processes=[], flow=False):       
    if isinstance(gen_axes, str):
        gen_axes = [gen_axes]
    if isinstance(base_processes, str):
        base_processes = [base_processes]
    df = pd.DataFrame({"Name":names})
    for axis in gen_axes:
        df[axis] = df["Name"].apply(lambda x, a=axis: decode_poi_bin(x, a))
        if flow:
            # set underflow to -1, overflow to max bin number+1
            max_bin = pd.to_numeric(df[axis],errors='coerce').max()
            df[axis] = df[axis].apply(lambda x, iu=max_bin: -1 if x=="U" else iu+1 if x=="O" else int(x) if x is not None else None)
        else:
            # set underflow and overflow to None
            df[axis] = df[axis].apply(lambda x: None if x in ["U","O",None] else int(x))

    # filter out rows with NaNs
    mask = ~df.isna().any(axis=1)
    # select rows from base process
    if len(base_processes):
        mask = mask & df["Name"].apply(lambda x, p=base_processes: any([x.startswith(p) for p in base_processes]))
    # gen bin selections
    for k, v in selections.items():
        mask = mask & (df[k] == v)

    filtered_df = df[mask]

    return filtered_df.sort_values(list(gen_axes)).index.to_list()

def select_pois(df, gen_axes=[], selections={}, base_processes=[], flow=False):
    return df.iloc[filter_poi_bins(df["Name"].values, gen_axes, selections=selections, base_processes=base_processes, flow=flow)]

def read_impacts_pois(fileobject, poi_type, scale=1.0, group=True, uncertainties=None):
    # read impacts of all pois from type 'poi_type'
    if is_root_file(fileobject):
        res = read_impacts_pois_root(fileobject, poi_type, group, uncertainties)
    elif is_h5_file(fileobject):
        res = read_impacts_pois_h5(fileobject, poi_type, group, uncertainties)
    else:
        raise IOError(f"Unknown fitresult format for object {fileobject}")

    df = pd.DataFrame({"Name":res[0], "value":res[1], "err_total":res[2], **res[3]})

    if scale != 1:
        df["value"] /= scale
        df["err_total"] /= scale
        for u in res[3].keys():
            df[u] /= scale

    return df

def read_impacts_pois_h5(h5file, poi_type, group=True, uncertainties=None):
    names = h5file[f"{poi_type}_names"][...].astype(str)
    centrals = h5file[f"{poi_type}_outvals"][...]

    npoi = len(names)
    # make matrix between POIs only; assume POIs come first
    totals = np.sqrt(np.diagonal(h5file[f"{poi_type}_outcov"][:npoi,:npoi]))

    impact_hist = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"

    impacts = h5file[impact_hist][...]

    if group:
        labels = h5file["hsystgroups"][...].astype(str)
        labels = np.append(labels, "stat")
        if len(labels)+1 == impacts.shape[1]:
            labels = np.append(labels, "binByBinStat")
    else:
        labels = h5file["hsysts"][...].astype(str)

    logger.debug(f"Load ucertainties")
    # pick uncertainties
    if uncertainties is None:
        uncertainties = {f"err_{k}": impacts[:,i] for i, k in enumerate(labels)}
    else:
        uncertainties = {f"err_{k}": impacts[:,i] for i, k in enumerate(labels) if k in uncertainties}

    return names, centrals, totals, uncertainties

def read_impacts_pois_root(rtfile, poi_type, group=True, uncertainties=None):   
    histname = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"
    if f"{histname};1" not in rtfile.keys():
        raise RuntimeError(f"Histogram {histname};1 not found in fitresult file")
        return None
    impacts = rtfile[histname].to_hist()

    # process names
    names = [k for k in impacts.axes[0]]

    logger.debug(f"Load ucertainties")
    # pick uncertainties
    if uncertainties is None:
        uncertainties = {f"err_{k}": impacts.values()[:,i] for i, k in enumerate(impacts.axes[1])}
    else:
        uncertainties = {f"err_{k}": impacts.values()[:,i] for i, k in enumerate(impacts.axes[1]) if k in uncertainties}

    # measured central value
    fitresults = rtfile["fitresults;1"]
    centrals = [fitresults[n].array()[0] for n in names]

    # total uncertainties
    totals = [fitresults[n+"_err"].array()[0] for n in names]

    return names, centrals, totals, uncertainties

def select_covariance_pois(cov, names, gen_axes=[], selections={}, base_processes=[], flow=False):
    indices = filter_poi_bins(names, gen_axes, selections=selections, base_processes=base_processes, flow=flow)

    # make matrix between selected POIs only
    new_cov = hist.Hist(
        hist.axis.Regular(bins=len(indices), start=0.5, stop=len(indices)+0.5, underflow=False, overflow=False), 
        hist.axis.Regular(bins=len(indices), start=0.5, stop=len(indices)+0.5, underflow=False, overflow=False), 
        storage=hist.storage.Double())
    new_cov.view(flow=False)[...] = cov.view(flow=False)[indices, :][:, indices]

    return new_cov

def load_covariance_pois(fitresult, poi_type="mu"):   
    if is_root_file(fitresult):
        values, names = load_covariance_pois_root(fitresult, poi_type)
    elif is_h5_file(fitresult):
        values, names = load_covariance_pois_h5(fitresult, poi_type)
    else:
        raise IOError(f"Unknown fitresult format for object {fitresult}")

    cov = hist.Hist(
        hist.axis.Regular(bins=len(names), start=0.5, stop=len(names)+0.5, underflow=False, overflow=False), 
        hist.axis.Regular(bins=len(names), start=0.5, stop=len(names)+0.5, underflow=False, overflow=False), 
        storage=hist.storage.Double())
    cov.view(flow=False)[...] = values
    return cov, names

def load_covariance_pois_root(rtfile, poi_type="mu"):   
    matrix_key = f"covariance_matrix_channel{poi_type}"
    if matrix_key not in [c.replace(";1","") for c in rtfile.keys()]:
        IOError(f"Histogram {matrix_key} was not found in the fit results file!")
    hist2d = rtfile[matrix_key].to_hist()
    names = [n for n in hist2d.axes[0]]
    hcov = hist2d.values()
    return hcov, names

def load_covariance_pois_h5(h5file, poi_type="mu"):   
    matrix_key = f"{poi_type}_outcov"
    names_key = f"{poi_type}_names"
    if matrix_key not in h5file.keys():
        IOError(f"Matrix {matrix_key} was not found in the fit results file!")
    if names_key not in h5file.keys():
        IOError(f"Names {names_key} not found in the fit results file!")
    
    names = h5file[names_key][...].astype(str)
    npoi = len(names)
    # make matrix between POIs only; assume POIs come first
    hcov = h5file[f"{poi_type}_outcov"][:npoi,:npoi]
    return hcov, names

def get_theoryfit_data(fitresult, axes, base_processes = ["W"], poi_type="pmaskedexp", flow=False):
    logger.info(f"Prepare theory fit: load measured differential cross secction distribution and covariance matrix")

    cov, names = load_covariance_pois(fitresult, poi_type)
    df = read_impacts_pois(fitresult, poi_type, group=False, uncertainties=[])

    # select POIs 
    all_axes = [a for b in axes for a in b]
    cov = select_covariance_pois(cov, names, gen_axes=all_axes, base_processes=base_processes, flow=flow)

    # write out unfolded data as flat 1D hist for each channel
    data = []
    for a, p in zip(axes, base_processes):
        df_c = select_pois(df, a, base_processes=p, flow=flow)

        hist_xsec = hist.Hist(
            hist.axis.Regular(bins=len(df_c), start=0.5, stop=len(df_c)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
        hist_xsec.view(flow=False)[...] = np.stack([df_c["value"].values, (df_c["err_total"].values)**2], axis=-1)

        data.append(hist_xsec.values(flow=False).flatten())

    return data, cov


