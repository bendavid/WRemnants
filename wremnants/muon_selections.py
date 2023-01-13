
def select_veto_muons(df, is_w_like):

    nMuons = 2 if is_w_like else 1

    # n.b. charge = -99 is a placeholder for invalid track refit/corrections (mostly just from tracks below
    # the pt threshold of 8 GeV in the nano production)
    df = df.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_correctedCharge != -99")
    df = df.Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")

    df = df.Filter(f"Sum(vetoMuons) == {nMuons}")

    return df

def select_good_muons(df, is_w_like, use_trackerMuons):

    nMuons = 2 if is_w_like else 1

    if use_trackerMuons:
        if dataset.group in ["Top", "Diboson"]:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_highPurity")
        else:
            df = df.Define("Muon_category", "Muon_isTracker && Muon_innerTrackOriginalAlgo != 13 && Muon_innerTrackOriginalAlgo != 14 && Muon_highPurity")
    else:
        df = df.Define("Muon_category", "Muon_isGlobal")

    goodMuonsSelection = "vetoMuons && Muon_mediumId && Muon_category"
    if is_w_like:
        # for w like we directly require isolated muons, for w we need non-isolated for qcd estimation
        goodMuonsSelection += " && Muon_pfRelIso04_all < 0.15"

    df = df.Define("goodMuons", goodMuonsSelection) 
    df = df.Filter(f"Sum(goodMuons) == {nMuons}")

    return df

def veto_electrons(df):

    df = df.Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")
    df = df.Filter("Sum(vetoElectrons) == 0")
    
    return df

# def select_standalone_muons(df, is_w_like):
