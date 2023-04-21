from utilities import differential
from wremnants import syst_tools, theory_tools

def define_gen_level(df, gen_level, dataset_name, mode="wmass"):
    # gen level definitions
    gen_levels = ["preFSR", "postFSR"]
    if gen_level not in gen_levels:
        raise ValueError(f"Unknown gen level '{gen_level}'! Supported gen level definitions are '{gen_levels}'.")

    modes = ["wmass", "wlike"]
    if mode not in modes:
        raise ValueError(f"Unknown mode '{mode}'! Supported modes are '{modes}'.")

    if gen_level == "preFSR":
        df = theory_tools.define_prefsr_vars(df)

        if mode == "wmass":
            df = df.Define("ptGen", "chargeVgen < 0 ? genl.pt() : genlanti.pt()")   
            df = df.Define("etaGen", "chargeVgen < 0 ? abs(genl.eta()) : abs(genlanti.eta())")
        elif mode == "wlike":
            df = df.Define("ptGen", "event % 2 == 0 ? genl.pt() : genlanti.pt()")
            df = df.Define("etaGen", "event % 2 == 0 ? abs(genl.eta()) : abs(genlanti.eta())")

    elif gen_level == "postFSR":

        if mode == "wmass":
            pdgId = -13 if "Wplusmunu" in dataset_name else 13
            df = df.Define("postFSR_mu", f"GenPart_status == 1 && (GenPart_statusFlags & 1) && GenPart_pdgId == {pdgId}")
            df = df.Define("postFSR_mu_idx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSR_mu])")
            df = df.Define("ptGen", "GenPart_pt[postFSR_mu][postFSR_mu_idx]")
            df = df.Define("etaGen", "abs(GenPart_eta[postFSR_mu][postFSR_mu_idx])")
        elif mode == "wlike":
            df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && GenPart_pdgId == 13")
            df = df.Define("postFSRantimuons", "GenPart_status == 1 && (GenPart_statusFlags & 1) && GenPart_pdgId == -13")
            df = df.Define("postFSRmuonIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRmuons])")
            df = df.Define("postFSRantimuonIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRantimuons])")
            df = df.Define("ptGen", "event % 2 == 0 ? GenPart_pt[postFSRmuons][postFSRmuonIdx] : GenPart_pt[postFSRantimuons][postFSRantimuonIdx]")
            df = df.Define("etaGen", "event % 2 == 0 ? abs(GenPart_eta[postFSRmuons][postFSRmuonIdx]) : abs(GenPart_eta[postFSRantimuons][postFSRantimuonIdx])")    
    
    if mode == "wlike":
        df = df.Define("qGen", "event % 2 == 0 ? -1 : 1")

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
