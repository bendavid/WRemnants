from utilities import differential
from wremnants import syst_tools, theory_tools, logging
from copy import deepcopy
import hist
import numpy as np
import pandas as pd
import h5py
import uproot

logger = logging.child_logger(__name__)

def add_out_of_acceptance(datasets, group, newGroupName=None):
    # Copy datasets from specified group to make out of acceptance contribution
    datasets_ooa = []
    for dataset in datasets:
        if dataset.group == group:
            ds = deepcopy(dataset)

            if newGroupName is None:
                ds.group = ds.group+"OOA"
            else:
                ds.group = newGroupName
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
        df = df.Alias("massVGen", "massVgen")
        df = df.Alias("ptVGen", "ptVgen")
        df = df.Alias("absYVGen", "absYVgen")

        if mode in ["wmass", "wlike"]:
            df = df.Define("mTWGen", "wrem::mt_2(genl.pt(), genl.phi(), genlanti.pt(), genlanti.phi())")   

        if mode == "wmass":
            df = df.Define("ptGen", "chargeVgen < 0 ? genl.pt() : genlanti.pt()")   
            df = df.Define("absEtaGen", "chargeVgen < 0 ? fabs(genl.eta()) : fabs(genlanti.eta())")

        if mode in ["wlike", "dilepton"]:
            df = df.Define("ptGen", "event % 2 == 0 ? genl.pt() : genlanti.pt()")
            df = df.Define("absEtaGen", "event % 2 == 0 ? fabs(genl.eta()) : fabs(genlanti.eta())")
            df = df.Define("ptOtherGen", "event % 2 == 0 ? genlanti.pt() : genl.pt()")
            df = df.Define("absEtaOtherGen", "event % 2 == 0 ? genlanti.pt() : genl.pt()")

    elif gen_level == "postFSR":
        df = theory_tools.define_postfsr_vars(df, mode=mode)

        df = df.Alias("ptGen", f"postfsrLep_pt")
        df = df.Alias("absEtaGen", f"postfsrLep_absEta")           

        if mode in ["wmass", "wlike"]:
            df = df.Alias("mTWGen", "postfsrMT")   
   
        if mode in ["wlike", "dilepton"]:
            df = df.Alias("ptOtherGen", "postfsrOtherLep_pt")
            df = df.Alias("absEtaOtherGen", f"postfsrOtherLep_absEta")                

            df = df.Alias("massVGen", "postfsrMV")
            df = df.Define("absYVGen", "fabs(postfsrYV)")  

        df = df.Alias("ptVGen", "postfsrPTV")      

    if mode == "wlike":
        df = df.Define("qGen", "event % 2 == 0 ? -1 : 1")

    return df

def select_fiducial_space(df, select=True, accept=True, mode="wmass", pt_min=None, pt_max=None, mass_min=60, mass_max=120, mtw_min=0, selections=[]):
    # Define a fiducial phase space and if select=True, either select events inside/outside
    # accept = True: select events in fiducial phase space 
    # accept = False: reject events in fiducial pahse space
    
    if mode == "wmass":
        selection = "(absEtaGen < 2.4)"        
    elif mode == "wlike":
        selection = f"""
            (absEtaGen < 2.4) && (absEtaOtherGen < 2.4) 
            && (ptOtherGen > {pt_min}) && (ptOtherGen < {pt_max})
            && (massVGen > {mass_min}) && (massVGen < {mass_max})
            """
    elif mode == "dilepton":
        selection = f"""
            (absEtaGen < 2.4) && (absEtaOtherGen < 2.4) 
            && (ptGen > {pt_min}) && (ptOtherGen > {pt_min}) 
            && (ptGen < {pt_max}) && (ptOtherGen < {pt_max}) 
            && (massVGen > {mass_min}) && (massVGen < {mass_max})
            """
    else:
        raise NotImplementedError(f"No fiducial phase space definiton found for mode '{mode}'!") 

    if mtw_min > 0:
        selection += f" && (mTWGen > {mtw_min})"

    for sel in selections:
        logger.debug(f"Add selection {sel} for fiducial phase space")
        selection += f" && ({sel})"

    df = df.Define("acceptance", selection)

    if select and accept:
        df = df.Filter("acceptance")
    elif select :
        df = df.Filter("acceptance == 0")

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

