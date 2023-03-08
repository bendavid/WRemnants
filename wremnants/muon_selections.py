
from wremnants import muon_calibration
from utilities.common import background_MCprocs as bkgMCprocs

def apply_met_filters(df):
    df = df.Filter("Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter")

    return df

def select_veto_muons(df, nMuons=1):

    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")
    df = df.Filter(f"Sum(vetoMuons) == {nMuons}")

    return df

def select_good_muons(df, nMuons=1, use_trackerMuons=False, use_isolation=False):

    if use_trackerMuons:
        if dataset.group in bkgMCprocs:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_highPurity")
        else:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_innerTrackOriginalAlgo != 13 && Muon_innerTrackOriginalAlgo != 14 && Muon_highPurity")
    else:
        df = df.Define("Muon_category", "Muon_isGlobal && Muon_highPurity")

    goodMuonsSelection = "vetoMuons && Muon_mediumId && Muon_category"
    if use_isolation:
        # for w like we directly require isolated muons, for w we need non-isolated for qcd estimation
        goodMuonsSelection += " && Muon_pfRelIso04_all < 0.15"

    df = df.Define("goodMuons", goodMuonsSelection) 
    df = df.Filter(f"Sum(goodMuons) == {nMuons}")

    return df

def define_trigger_muons(df):

    # mu- for even event numbers, mu+ for odd event numbers
    df = df.Define("trigMuons_charge0", "event % 2 == 0 ? -1 : 1")
    df = df.Define("nonTrigMuons_charge0", "-trigMuons_charge0")

    df = df.Define("trigMuons", "goodMuons && Muon_correctedCharge == trigMuons_charge0")
    df = df.Define("nonTrigMuons", "goodMuons && Muon_correctedCharge == nonTrigMuons_charge0")

    df = muon_calibration.define_corrected_reco_muon_kinematics(df, "trigMuons", ["pt", "eta", "phi"])
    df = muon_calibration.define_corrected_reco_muon_kinematics(df, "nonTrigMuons", ["pt", "eta", "phi"])

    return df

def select_z_candidate(df, ptLow, ptHigh):

    df = df.Filter("Sum(trigMuons) == 1 && Sum(nonTrigMuons) == 1")
    df = df.Filter(f"trigMuons_pt0 > {ptLow} && trigMuons_pt0 < {ptHigh}")
    df = df.Filter(f"nonTrigMuons_pt0 > {ptLow} && nonTrigMuons_pt0 < {ptHigh}")

    df = df.Define("trigMuons_mom4", "ROOT::Math::PtEtaPhiMVector(trigMuons_pt0, trigMuons_eta0, trigMuons_phi0, wrem::muon_mass)")
    df = df.Define("nonTrigMuons_mom4", "ROOT::Math::PtEtaPhiMVector(nonTrigMuons_pt0, nonTrigMuons_eta0, nonTrigMuons_phi0, wrem::muon_mass)")
    df = df.Define("ll_mom4", "ROOT::Math::PxPyPzEVector(trigMuons_mom4)+ROOT::Math::PxPyPzEVector(nonTrigMuons_mom4)")
    df = df.Define("mll", "ll_mom4.mass()")

    df = df.Filter("mll >= 60. && mll < 120.")

    return df

def apply_triggermatching_muon(df, dataset, muon_eta, muon_phi):
    if dataset.group in bkgMCprocs:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    else:
        df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Filter(f"wrem::hasTriggerMatch({muon_eta},{muon_phi},TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
    
    return df

def veto_electrons(df):

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")
    df = df.Filter("Sum(vetoElectrons) == 0")
    
    return df

def select_standalone_muons(df, dataset, use_trackerMuons=False, muons="goodMuons", idx=0):

    from utilities.common import muonEfficiency_standaloneNumberOfValidHits as nHitsSA

    #standalone quantities, currently only in data and W/Z samples
    if dataset.group in bkgMCprocs:
        df = df.Alias(f"{muons}_SApt{idx}",  f"{muons}_pt{idx}")
        df = df.Alias(f"{muons}_SAeta{idx}", f"{muons}_eta{idx}")
        df = df.Alias(f"{muons}_SAphi{idx}", f"{muons}_phi{idx}")
    elif use_trackerMuons:
        # try to use standalone variables when possible
        df = df.Define(f"{muons}_SApt{idx}",  f"Muon_isStandalone[{muons}][{idx}] ? Muon_standalonePt[{muons}][{idx}] : {muons}_pt{idx}")
        df = df.Define(f"{muons}_SAeta{idx}", f"Muon_isStandalone[{muons}][{idx}] ? Muon_standaloneEta[{muons}][{idx}] : {muons}_eta{idx}")
        df = df.Define(f"{muons}_SAphi{idx}", f"Muon_isStandalone[{muons}][{idx}] ? Muon_standalonePhi[{muons}][{idx}] : {muons}_phi{idx}")
    else:
        df = df.Define(f"{muons}_SApt{idx}",  f"Muon_standalonePt[{muons}][{idx}]")
        df = df.Define(f"{muons}_SAeta{idx}", f"Muon_standaloneEta[{muons}][{idx}]")
        df = df.Define(f"{muons}_SAphi{idx}", f"Muon_standalonePhi[{muons}][{idx}]")
    
    # the next cuts are mainly needed for consistency with the reco efficiency measurement for the case with global muons
    # note, when SA does not exist this cut is still fine because of how we define these variables
    df = df.Filter(f"{muons}_SApt{idx} > 15.0 && wrem::deltaR2({muons}_SAeta{idx}, {muons}_SAphi{idx}, {muons}_eta{idx}, {muons}_phi{idx}) < 0.09")
    if nHitsSA > 0 and not use_trackerMuons and not dataset.group in bkgMCprocs:
        df = df.Filter(f"Muon_standaloneNumberOfValidHits[{muons}][{idx}] >= {nHitsSA}")

    return df
