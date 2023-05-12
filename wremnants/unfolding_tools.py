from utilities import differential
from wremnants import syst_tools, theory_tools

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

        if mode == "wmass":
            df = df.Define("ptGen", "chargeVgen < 0 ? genl.pt() : genlanti.pt()")   
            df = df.Define("etaGen", "chargeVgen < 0 ? abs(genl.eta()) : abs(genlanti.eta())")
        else:
            # needed for fiducial phase space definition
            df = df.Alias("muGen", "genl")
            df = df.Alias("antimuGen", "genlanti")
            df = df.Alias("massVGen", "massVgen")

            if mode == "wlike":
                df = df.Define("ptGen", "event % 2 == 0 ? genl.pt() : genlanti.pt()")
                df = df.Define("etaGen", "event % 2 == 0 ? abs(genl.eta()) : abs(genlanti.eta())")

            elif mode == "dilepton":
                df = df.Alias("ptVGen", "ptVgen")
                df = df.Alias("yVGen", "yVgen")

    elif gen_level == "postFSR":

        df = df.Define("postFSRmus", "GenPart_status == 1 && (GenPart_statusFlags & 1) && GenPart_pdgId == 13")
        df = df.Define("postFSRantimus", "GenPart_status == 1 && (GenPart_statusFlags & 1) && GenPart_pdgId == -13")
        df = df.Define("postFSRmuIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRmus])")
        df = df.Define("postFSRantimuIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRantimus])")

        if mode == "wmass":
            if "Wplusmunu" in dataset_name:
                idx = "postFSRantimuIdx" 
                muons = "postFSRantimus"
            else:
                idx = "postFSRmuIdx" 
                muons = "postFSRmus"

            df = df.Define("ptGen", f"GenPart_pt[{muons}][{idx}]")
            df = df.Define("etaGen", f"abs(GenPart_eta[{muons}][{idx}])")                

        else:
            df = df.Define("muGen", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRmus][postFSRmuIdx], GenPart_eta[postFSRmus][postFSRmuIdx], GenPart_phi[postFSRmus][postFSRmuIdx], GenPart_mass[postFSRmus][postFSRmuIdx])")
            df = df.Define("antimuGen", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRantimus][postFSRantimuIdx], GenPart_eta[postFSRantimus][postFSRantimuIdx], GenPart_phi[postFSRantimus][postFSRantimuIdx], GenPart_mass[postFSRantimus][postFSRantimuIdx])")
            df = df.Define("VGen", "ROOT::Math::PxPyPzEVector(muGen)+ROOT::Math::PxPyPzEVector(antimuGen)")
            df = df.Define("massVGen", "VGen.mass()")

            if mode == "wlike":
                df = df.Define("ptGen", "event % 2 == 0 ? GenPart_pt[postFSRmus][postFSRmuIdx] : GenPart_pt[postFSRantimus][postFSRantimuIdx]")
                df = df.Define("etaGen", "event % 2 == 0 ? abs(GenPart_eta[postFSRmus][postFSRmuIdx]) : abs(GenPart_eta[postFSRantimus][postFSRantimuIdx])")    

            if mode == "dilepton":
                df = df.Define("ptVGen", "VGen.pt()")
                df = df.Define("yVGen", "VGen.Rapidity()")  
    
    if mode == "wlike":
        df = df.Define("qGen", "event % 2 == 0 ? -1 : 1")

    return df

def define_fiducial_space(df, mode="wmass", pt_min=26, pt_max=55, gen_var=None):
    # Define a fiducial phase space where all out of acceptance contribution is stored in a single bin 
    #   (instead of having several overflow/underflow bins in the gen axes). 

    if mode == "wmass":
        df = df.Define("fiducial", f"""
            (etaGen < 2.4) 
            && (ptGen > {pt_min}) 
            && (ptGen < {pt_max}) 
            """)
    elif mode == "wlike":
        df = df.Define("fiducial", f"""
            (abs(muGen.eta()) < 2.4) && (abs(antimuGen.eta()) < 2.4) 
            && (muGen.pt() > {pt_min}) && (antimuGen.pt() > {pt_min}) 
            && (muGen.pt() < {pt_max}) && (antimuGen.pt() < {pt_max}) 
            && (massVGen > 60) & (massVGen < 120)
            """)

    elif mode == "dilepton":
        df = df.Define("fiducial", f"""
            (abs(muGen.eta()) < 2.4) && (abs(antimuGen.eta()) < 2.4) 
            && (muGen.pt() > {pt_min}) && (antimuGen.pt() > {pt_min}) 
            && (muGen.pt() < {pt_max}) && (antimuGen.pt() < {pt_max}) 
            && (massVGen > 60) & (massVGen < 120)
            && (ptVGen < 100)
            """)

    return df        

def add_xnorm_histograms(results, df, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols):
    # add histograms before any selection
    df_xnorm = df
    df_xnorm = df_xnorm.DefinePerSample("exp_weight", "1.0")

    df_xnorm = theory_tools.define_theory_weights_and_corrs(df_xnorm, dataset_name, corr_helpers, args)

    df_xnorm = df_xnorm.DefinePerSample("count", "0.5")

    xnorm_axes = [*unfolding_axes, differential.axis_xnorm]
    xnorm_cols = [*unfolding_cols, "count"]
    
    results.append(df_xnorm.HistoBoost("xnorm", xnorm_axes, [*xnorm_cols, "nominal_weight"]))

    syst_tools.add_theory_hists(results, df_xnorm, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, xnorm_axes, xnorm_cols, base_name="xnorm")
