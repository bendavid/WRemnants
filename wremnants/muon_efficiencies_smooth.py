import ROOT
import pathlib
import hist
import narf
import numpy as np
import boost_histogram as bh
import logging
import pickle
import lz4.frame

from utilities import common

ROOT.gInterpreter.Declare('#include "muon_efficiencies_smooth.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_efficiency_helpers_smooth(filename = data_dir + "/testMuonSF/allSmooth_GtoH.root",
                                        era = None, is_w_like = False, max_pt = np.inf):

    eradict = { "2016PreVFP" : "BtoF",
                "2016PostVFP" : "GtoH" }
    eratag = eradict[era]

    axis_eta_eff = None
    axis_pt_eff = None
    # categorical axes in python bindings always have an overflow bin, so use a regular
    # axis for the charge
    axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "SF charge")
    isoEff_types = ["iso", "isonotrig", "antiiso", "antiisonotrig"]
    allEff_types = ["reco", "tracking", "idip", "trigger"] + isoEff_types
    axis_allEff_type = hist.axis.StrCategory(allEff_types, name = "allEff_type")
    axis_nom_syst = hist.axis.Integer(0, 2, underflow = False, overflow =False, name = "nom-syst") # only one syst for now (and the nominal in the first bin)
    axis_down_up = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar") # for the stat variations, to avoid stacking on same axis

    charges = { -1. : "minus", 1. : "plus" }
    chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
    
    fin = ROOT.TFile.Open(filename)
    if fin is None or fin.IsZombie():
        print(f"Error: file {filename} was not opened correctly")
        quit()
        
    ## start with NOMI and SYST
    sf_syst = None
    ############
    # all steps, for the nominal and systematics
    for charge, charge_tag in charges.items():
        for eff_type in allEff_types:
            # for iso can use histogram for efficiency variation only in data (the only one we have for now)
            nameTag = "nomiAndAlt_onlyDataVar" if any(x in eff_type for x in ["iso", "antiiso"]) else "nomiAndAlt"
            chargeTag = charge_tag if eff_type in chargeDependentSteps else "both"
            hist_name = f"SF_{nameTag}_{eratag}_{eff_type}_{chargeTag}"
            hist_root = fin.Get(hist_name)
            #print(ROOT.AddressOf(hist_root))
            if hist_root is None:
                print(f"Error: {hist_name} onot found in file {filename}")
                quit()
            #print(f"syst: {eff_type} -> {hist_name}")

            hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt", "nomi-statUpDown-syst"])

            if sf_syst is None:
                axis_eta_eff = hist_hist.axes[0]
                axis_pt_eff = hist_hist.axes[1]
                # store all systs (currently only 1) with the nominal, for all efficiency steps
                sf_syst = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge, axis_allEff_type, axis_nom_syst, name = "sf_syst", storage = hist.storage.Weight())
            # this axis might change for different histograms, because of a different number of effStat variations
            axis_nomiAlt_eff = hist_hist.axes[2]
            # could use max_pt to remove some of the pt bins for the input histogram
            # extract nominal (first bin that is not underflow) and put in corresponding bin of destination (bin 0 is the first bin because no underflow)
            sf_syst.view(flow=True)[:, :, axis_charge.index(charge), axis_allEff_type.index(eff_type), 0] = hist_hist.view(flow=True)[:,:,1]
            # extract syst (last bin except overflow) and put in corresponding bin of destination (bin 1 is the second bin because no underflow)
            sf_syst.view(flow=True)[:, :, axis_charge.index(charge), axis_allEff_type.index(eff_type), 1] = hist_hist.view(flow=True)[:,:,axis_nomiAlt_eff.extent-2]

    # set overflow and underflow eta-pt bins equal to adjacent bins
    sf_syst.view(flow=True)[0, ...] = sf_syst.view(flow=True)[1, ...]
    sf_syst.view(flow=True)[axis_eta_eff.extent-1, ...] = sf_syst.view(flow=True)[axis_eta_eff.extent-2, ...]
    sf_syst.view(flow=True)[:, 0, ...] = sf_syst.view(flow=True)[:, 1, ...]
    sf_syst.view(flow=True)[:, axis_pt_eff.extent-1, ...] = sf_syst.view(flow=True)[:, axis_pt_eff.extent-2, ...]

    # might be convenient to store them for easier usage
    # outpkl = filename.replace(".root", ".pkl.lz4")
    # with lz4.frame.open(outpkl, "wb") as f:
    #     pickle.dump(sf_syst, f, protocol = pickle.HIGHEST_PROTOCOL)
    # print(f"Saved file {outpkl}")
    
    sf_syst_pyroot = narf.hist_to_pyroot_boost(sf_syst)
    
    # nomi and syst are stored in the same histogram, just use different helpers to override the () operator for now, until RDF is improved
    helper = ROOT.wrem.muon_efficiency_smooth_helper[str(is_w_like).lower(), type(sf_syst_pyroot)]( ROOT.std.move(sf_syst_pyroot) )
    helper_syst = ROOT.wrem.muon_efficiency_smooth_helper_syst[str(is_w_like).lower(), type(sf_syst_pyroot)]( helper )
    # define axis for syst variations with all steps
    axis_all = hist.axis.Integer(0, 5, underflow = False, overflow = False, name = "reco-tracking-idip-trigger-iso")
    helper_syst.tensor_axes = [axis_all]
    ##############
    ## now the EFFSTAT
    
    # for the stat histograms, an axis will contain the efficiency type, in some cases it might have a single bin (e.g. tracking, trigger)
    # but for isolation needs more, because of the different cases to handle inside the helper class (iso-antiiso with/without trigger)
    # could actually stack reco-idip-trig together, but better to have the helper deal with a single histogram at a time
    # then there will be a specialization of the class for isolation
    # in this way the name of the (single) effType axis can be used to activate some internal behaviour as needed (for example trigger won't have variations for not triggering lepton)
    
    effStat_manager = {"sf_reco": {"nPtEigenBins" : 0,
                                   "axisLabels" : ["reco"],
                                   "boostHist" : None,
                                   "helper" : None},
                       "sf_tracking": {"nPtEigenBins" : 0,
                                       "axisLabels" : ["tracking"],
                                       "boostHist" : None,
                                       "helper" : None},
                       "sf_idip": {"nPtEigenBins" : 0,
                                   "axisLabels" : ["idip"],
                                   "boostHist" : None,
                                   "helper" : None},
                       "sf_trigger": {"nPtEigenBins" : 0,
                                   "axisLabels" : ["trigger"],
                                   "boostHist" : None,
                                   "helper" : None},
                       "sf_iso_effData": {"nPtEigenBins" : 0,
                                          "axisLabels" : ["iso", "isonotrig", "antiiso", "antiisonotrig"],
                                          "boostHist" : None,
                                          "helper" : None},
                       "sf_iso_effMC": {"nPtEigenBins" : 0,
                                        "axisLabels" : ["iso", "isonotrig", "antiiso", "antiisonotrig"],
                                        "boostHist" : None,
                                        "helper" : None},
    }

    for effStatKey in effStat_manager.keys():
        down_nom_up_effStat_axis = None
        axis_eff_type = None
        for charge, charge_tag in charges.items():
            for eff_type in effStat_manager[effStatKey]["axisLabels"]:
                print(eff_type)
                nameTag = "nomiAndAlt"
                if "effData" in effStatKey:
                    nameTag += "_onlyDataVar"
                elif "effMC" in effStatKey:
                    nameTag += "_onlyMCVar"
                chargeTag = charge_tag if eff_type in chargeDependentSteps else "both"
                hist_name = f"SF_{nameTag}_{eratag}_{eff_type}_{chargeTag}"
                hist_root = fin.Get(hist_name)
                # print(f"stat: {effStatKey}|{eff_type} -> {hist_name}")
                # print(ROOT.AddressOf(hist_root))
                if hist_root is None:
                    print(f"Error: {hist_name} onot found in file {filename}")
                    quit()

                if down_nom_up_effStat_axis is None:
                    # Z axis has nominal, 2*N effStat (first N up then N down) and one syst in the last bin
                    nPtEigenBins = int(int(hist_root.GetNbinsZ() - 2)/ 2)
                    effStat_manager[effStatKey]["nPtEigenBins"] = nPtEigenBins
                    down_nom_up_effStat_axis = hist.axis.Regular(int(1+2*nPtEigenBins), -0.5-nPtEigenBins, 0.5+nPtEigenBins, underflow=False, overflow=False, name = "downNomUpVar")

                hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt", "nomi-statUpDown-syst"])

                if effStat_manager[effStatKey]["boostHist"] is None:
                    axis_eta_eff = hist_hist.axes[0]
                    axis_pt_eff = hist_hist.axes[1]
                    axis_eff_type = hist.axis.StrCategory(effStat_manager[effStatKey]["axisLabels"], name = f"{effStatKey}_eff_type")
                    effStat_manager[effStatKey]["boostHist"] = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge,
                                                                         axis_eff_type,
                                                                         down_nom_up_effStat_axis,
                                                                         name = effStatKey,
                                                                         storage = hist.storage.Weight())
                    
                # extract nominal (first bin that is not underflow)
                effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type.index(eff_type), down_nom_up_effStat_axis.index(0)] = hist_hist.view(flow=True)[:,:,1]
                # now extract the stat variations
                # note that for hist_hist the nominal bin of the third axis is number 1 when the underflow is active
                for iup in range(1, 1 + nPtEigenBins):
                    # up variations
                    effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type.index(eff_type), down_nom_up_effStat_axis.index(iup)] = hist_hist.view(flow=True)[:,:,1+iup]
                    # down variations
                    idown = -1 * iup
                    effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type.index(eff_type), down_nom_up_effStat_axis.index(idown)] = hist_hist.view(flow=True)[:,:,1+iup+nPtEigenBins]
                
        # set overflow and underflow equal to adjacent bins
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[0, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[1, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[axis_eta_eff.extent-1, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[axis_eta_eff.extent-2, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, 0, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, 1, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, axis_pt_eff.extent-1, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, axis_pt_eff.extent-2, ...]

        netabins = axis_eta_eff.size
        sf_stat_pyroot = narf.hist_to_pyroot_boost(effStat_manager[effStatKey]["boostHist"])
        if "sf_iso" in effStatKey:
            helper_stat = ROOT.wrem.muon_efficiency_smooth_helper_stat_iso[str(is_w_like).lower(), netabins, nPtEigenBins, type(sf_stat_pyroot)]( ROOT.std.move(sf_stat_pyroot) )
        else:
            helper_stat = ROOT.wrem.muon_efficiency_smooth_helper_stat[str(is_w_like).lower(), netabins, nPtEigenBins, type(sf_stat_pyroot)]( ROOT.std.move(sf_stat_pyroot) )
        # make new versions of these axes without overflow/underflow to index the tensor
        if isinstance(axis_eta_eff, bh.axis.Regular):
            axis_eta_eff_tensor = hist.axis.Regular(axis_eta_eff.size, axis_eta_eff.edges[0], axis_eta_eff.edges[-1], name = axis_eta_eff.name, overflow = False, underflow = False)
        elif isinstance(axis_eta_eff, bh.axis.Variable):
            axis_eta_eff_tensor = hist.axis.Variable(axis_eta_eff.edges, name = axis_eta_eff.name, overflow = False, underflow = False)
        axis_ptEigen_eff_tensor = hist.axis.Integer(0, effStat_manager[effStatKey]["nPtEigenBins"], underflow = False, overflow =False, name = "nPtEigenBins")    
        helper_stat.tensor_axes = [axis_eta_eff_tensor, axis_ptEigen_eff_tensor, axis_charge, axis_down_up]
        effStat_manager[effStatKey]["helper"] = helper_stat

    fin.Close()
    ####
    # return nomi, effsyst, and a dictionary with effStat to use them by name
    return helper, helper_syst, {k : effStat_manager[k]["helper"] for k in effStat_manager.keys()}
