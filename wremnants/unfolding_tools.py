from utilities import differential
from wremnants import syst_tools, theory_tools, logging
from copy import deepcopy
import hist
import numpy as np
import pandas as pd

logger = logging.child_logger(__name__)

def add_out_of_acceptance(datasets, group):
    # Copy datasets from specified group to make out of acceptance contribution
    datasets_ooa = []
    for dataset in datasets:
        if dataset.group == group:
            ds = deepcopy(dataset)

            ds.group = "Bkg"+ds.group
            ds.out_of_acceptance = True

            datasets_ooa.append(ds)

    return datasets + datasets_ooa

def define_gen_level(df, gen_level, dataset_name, mode="wmass"):
    # gen level definitions
    gen_levels = ["preFSR", "postFSR"]
    if gen_level not in gen_levels:
        raise ValueError(f"Unknown gen level '{gen_level}'! Supported gen level definitions are '{gen_levels}'.")

    modes = ["wmass", "wlike", "dilepton"]
    if mode not in modes:
        raise ValueError(f"Unknown mode '{mode}'! Supported modes are '{modes}'.")

    if gen_level == "preFSR":
        df = theory_tools.define_prefsr_vars(df)

        # needed for fiducial phase space definition
        df = df.Alias("lepGen", "genl")
        df = df.Alias("antilepGen", "genlanti")

        df = df.Alias("massVGen", "massVgen")
        df = df.Alias("ptVGen", "ptVgen")
        df = df.Alias("absYVGen", "absYVgen")

        if mode in ["wmass", "wlike"]:
            df = df.Define("mTWGen", "wrem::mt_2(genl.pt(), genl.phi(), genlanti.pt(), genlanti.phi())")   

        if mode == "wmass":
            df = df.Define("ptGen", "chargeVgen < 0 ? genl.pt() : genlanti.pt()")   
            df = df.Define("absEtaGen", "chargeVgen < 0 ? fabs(genl.eta()) : fabs(genlanti.eta())")
        elif mode == "wlike":
            df = df.Define("ptGen", "event % 2 == 0 ? genl.pt() : genlanti.pt()")
            df = df.Define("absEtaGen", "event % 2 == 0 ? fabs(genl.eta()) : fabs(genlanti.eta())")

    elif gen_level == "postFSR":

        df = df.Define("postFSRleps", "GenPart_status == 1 && (GenPart_statusFlags & 1) && (GenPart_pdgId >= 11 && GenPart_pdgId <= 16)")
        df = df.Define("postFSRantileps", "GenPart_status == 1 && (GenPart_statusFlags & 1) && (GenPart_pdgId <= -11 && GenPart_pdgId >= -16)")
        df = df.Define("postFSRlepIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRleps])")
        df = df.Define("postFSRantilepIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRantileps])")

        if mode in ["wmass", "wlike"]:
            df = df.Define("mTWGen", "wrem::mt_2(GenPart_pt[postFSRleps][postFSRlepIdx], GenPart_phi[postFSRleps][postFSRlepIdx], GenPart_pt[postFSRantileps][postFSRantilepIdx], GenPart_phi[postFSRantileps][postFSRantilepIdx])")   

        if mode == "wmass":
            if "Wplus" in dataset_name:
                idx = "postFSRantilepIdx" 
                muons = "postFSRantileps"
            else:
                idx = "postFSRlepIdx" 
                muons = "postFSRleps"

            df = df.Define("ptGen", f"GenPart_pt[{muons}][{idx}]")
            df = df.Define("absEtaGen", f"fabs(GenPart_eta[{muons}][{idx}])")                
        elif mode == "wlike":
            df = df.Define("ptGen", "event % 2 == 0 ? GenPart_pt[postFSRleps][postFSRlepIdx] : GenPart_pt[postFSRantileps][postFSRantilepIdx]")
            df = df.Define("absEtaGen", "event % 2 == 0 ? fabs(GenPart_eta[postFSRleps][postFSRlepIdx]) : fabs(GenPart_eta[postFSRantileps][postFSRantilepIdx])")    

        df = df.Define("lepGen", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRleps][postFSRlepIdx], GenPart_eta[postFSRleps][postFSRlepIdx], GenPart_phi[postFSRleps][postFSRlepIdx], GenPart_mass[postFSRleps][postFSRlepIdx])")
        df = df.Define("antilepGen", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRantileps][postFSRantilepIdx], GenPart_eta[postFSRantileps][postFSRantilepIdx], GenPart_phi[postFSRantileps][postFSRantilepIdx], GenPart_mass[postFSRantileps][postFSRantilepIdx])")
        df = df.Define("VGen", "ROOT::Math::PxPyPzEVector(lepGen)+ROOT::Math::PxPyPzEVector(antilepGen)")

        df = df.Define("massVGen", "VGen.mass()")
        df = df.Define("ptVGen", "VGen.pt()")
        df = df.Define("absYVGen", "fabs(VGen.Rapidity())")  
    
    if mode == "wlike":
        df = df.Define("qGen", "event % 2 == 0 ? -1 : 1")

    return df

def select_fiducial_space(df, accept=True, mode="wmass", pt_min=26, pt_max=55, mass_min=60, mass_max=120, mtw_min=0, selections=[]):
    # Define a fiducial phase space and either select events inside/outside
    # accept = True: select events in fiducial phase space 
    # accept = False: reject events in fiducial pahse space
    
    if mode == "wmass":
        selection = f"""
            (ptGen > {pt_min}) && (ptGen < {pt_max}) 
            && (absEtaGen < 2.4)"""
    elif mode in ["wlike", "dilepton"]:
        selection = f"""
            (fabs(lepGen.eta()) < 2.4) && (fabs(antilepGen.eta()) < 2.4) 
            && (lepGen.pt() > {pt_min}) && (antilepGen.pt() > {pt_min}) 
            && (lepGen.pt() < {pt_max}) && (antilepGen.pt() < {pt_max}) 
            && (massVGen > {mass_min}) && (massVGen < {mass_max})
            """
    else:
        raise NotImplementedError(f"No fiducial phase space definiton found for mode '{mode}'!") 

    if mtw_min > 0:
        selection += f" && (mTWGen > {mtw_min})"

    for sel in selections:
        logger.debug(f"Add selection {sel} for fiducial phase space")
        selection += f" && ({sel})"

    df = df.Define("fiducial", selection)

    if accept:
        df = df.Filter("fiducial")
    else:
        df = df.Filter("fiducial == 0")

    return df

def add_xnorm_histograms(results, df, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols):
    # add histograms before any selection
    df_xnorm = df
    df_xnorm = df_xnorm.DefinePerSample("exp_weight", "1.0")

    df_xnorm = theory_tools.define_theory_weights_and_corrs(df_xnorm, dataset_name, corr_helpers, args)

    df_xnorm = df_xnorm.Define("xnorm", "0.5")

    axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)

    xnorm_axes = [axis_xnorm, *unfolding_axes]
    xnorm_cols = ["xnorm", *unfolding_cols]
    
    results.append(df_xnorm.HistoBoost("xnorm", xnorm_axes, [*xnorm_cols, "nominal_weight"]))

    syst_tools.add_theory_hists(results, df_xnorm, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, xnorm_axes, xnorm_cols, base_name="xnorm")


### for plotting

def get_bin(name, var):
    name_split = name.split(var)
    if len(name_split) == 1:
        return -1
    else:
        return int(name_split[-1].split("_")[0])

def getProcessBins(name, axes=["qGen", "ptGen", "absEtaGen", "ptVGen", "absYVGen"], base_processes=["Zmumu","Wmunu", "Z", "W"]):
    res = {
        x: get_bin(name, x) if get_bin(name, x) else 0 for x in axes
    }
    res["name"] = name

    for p in base_processes:
        if p in name:
            res["proc"] = p
            break

    if not res.get("proc", False):
        res["proc"] = None
    
    return res

def matrix_poi(rfile, poi_type="mu", base_process=None, axes=None, keys=None):   
    if isinstance(poi_type, list):
        poi_type = poi_type[0]
    
    matrix = f"covariance_matrix_channel{poi_type}"

    if matrix not in [c.replace(";1","") for c in rfile.keys()]:
        logger.error(f"Histogram {matrix} was not found in the fit results file!")
        return

    hist2d = rfile[matrix].to_hist()

    # select signal parameters
    key = matrix.split("channel")[-1]
    xentries = [(i, hist2d.axes[0][i]) for i in range(len(hist2d.axes[0])) if hist2d.axes[0][i].endswith(key)]

    if base_process is not None:
        xentries = [x for x in xentries if base_process in x[1]]  

    if keys is not None:
        xentries = [v for v in filter(lambda x, keys=keys: all([f"_{k}_" in x[1] for k in keys]), xentries)]

    if axes is not None:
        if isinstance(axes, str):
            axes = [axes]
        
        # select specified axes
        xentries = [v for v in filter(lambda x, axes=axes: all([f"_{a}" in x[1] for a in axes]), xentries)]

        # sort them in the specified order
        xentries = sorted(xentries, key=lambda x, axes=axes: [get_bin(x[1], a) for a in axes], reverse=False)

    # make matrix between POIs only
    cov_mat = np.zeros((len(xentries), len(xentries)))
    for i, ia in enumerate(xentries):
        for j, ja in enumerate(xentries):
            cov_mat[i][j] = hist2d[ia[0], ja[0]]

    hist_cov = hist.Hist(
        hist.axis.Regular(bins=len(xentries), start=0.5, stop=len(xentries)+0.5, underflow=False, overflow=False), 
        hist.axis.Regular(bins=len(xentries), start=0.5, stop=len(xentries)+0.5, underflow=False, overflow=False), 
        storage=hist.storage.Double())
    hist_cov.view(flow=False)[...] = cov_mat

    return hist_cov

def get_results(rtfile, poi_type, scale=1.0, group=True, uncertainties=None, gen_axes=["qGen", "ptGen", "absEtaGen", "ptVGen", "absYVGen"]):
    if isinstance(poi_type, list):
        results = []
        for p in poi_type:
            result = get_results(rtfile, p, scale=scale, group=group, uncertainties=uncertainties)
            if result is not None:
                results.append(result)
        return pd.concat(results)

    results = []

    fitresult = rtfile["fitresults"]

    histname = f"nuisance_group_impact_{poi_type}" if group else f"nuisance_impact_{poi_type}"

    if f"{histname};1" not in rtfile.keys():
        logger.debug(f"Histogram {histname};1 not found in fitresult file")
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
    centrals = [fitresult[n].array()[0] for n in names]

    # total uncertainties
    totals = [fitresult[n+"_err"].array()[0] for n in names]

    df = pd.DataFrame({"Name":names, "value":centrals, "err_total":totals, **uncertainties})

    if scale != 1:
        df["value"] /= scale
        df["err_total"] /= scale
        for u in uncertainties.keys():
            df[u] /= scale

    # try to decode the name string into bin number
    for axis in gen_axes:
        df[axis] = df["Name"].apply(lambda x, a=axis: get_bin(x, a))

    df = df.sort_values(gen_axes, ignore_index=True)
    return df