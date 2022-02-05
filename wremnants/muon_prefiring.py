import ROOT
import pathlib

ROOT.gInterpreter.Declare('#include "muon_prefiring.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_prefiring_helpers(filename = data_dir + "/testMuonSF/L1MuonPrefiringParametriations_histograms.root", era = None):

    fin = ROOT.TFile.Open(filename);

    if era == "H":
        hparameters = fin.Get("L1prefiring_muonparam_2016H")
    elif era != "GToH":
        # BG should be like preVFP, but more data was used to derive corrections
        #hMuonPrefiringNew[era] = *(dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016preVFP")));
        hparameters = fin.Get("L1prefiring_muonparam_2016BG")
    else:
        hparameters = fin.Get("L1prefiring_muonparam_2016postVFP")

    netabins = hparameters.GetXaxis().GetNbins()

    hparameters_hotspot = fin.Get("L1prefiring_muonparam_2016_hotspot")

    helper = ROOT.wrem.muon_prefiring_helper(hparameters, hparameters_hotspot)

    # histograms are copied into the helper class so we can close the file without detaching them
    fin.Close()

    helper_stat = ROOT.wrem.muon_prefiring_helper_stat[netabins](helper)
    helper_syst = ROOT.wrem.muon_prefiring_helper_syst(helper)

    return helper, helper_stat, helper_syst
