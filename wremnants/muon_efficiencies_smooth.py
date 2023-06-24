import ROOT
import pathlib
import hist
import narf
import numpy as np
import boost_histogram as bh
import pickle
import lz4.frame
import pdb
import copy

from utilities import common, logging
logger = logging.child_logger(__name__)

narf.clingutils.Declare('#include "muon_efficiencies_smooth.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def cloneAxis(ax, overflow=False, underflow=False):
    if isinstance(ax, bh.axis.Regular):
        #logger.debug(f"Cloning regular axis {ax.name} with {ax.size} bins from {ax.edges[0]} to {ax.edges[-1]}, and underflow/overflow = {underflow}/{overflow}")
        newax = hist.axis.Regular(ax.size, ax.edges[0], ax.edges[-1], name=ax.name, overflow=overflow, underflow=underflow)
    elif isinstance(ax, bh.axis.Variable):
        #logger.debug(f"Cloning Variable axis {ax.name} with {len(ax.edges)-1} bins from {ax.edges[0]} to {ax.edges[-1]}, and underflow/overflow = {underflow}/{overflow}")
        newax = hist.axis.Variable(ax.edges, name=ax.name, overflow=overflow, underflow=underflow)
    else:
        logger.error(f"In cloneAxis(): only Regular and Variable axes are supported for now")
        quit()
    return newax

def make_muon_efficiency_helpers_smooth(filename = data_dir + "/testMuonSF/allSmooth_GtoH.root",
                                        era = None, is_w_like = False, max_pt = np.inf,
                                        isoEfficiencySmoothing = False,
                                        smooth3D=False):
    logger.debug(f"Make efficiency helper smooth")

    eradict = { "2016PreVFP" : "BtoF",
                "2016PostVFP" : "GtoH" }
    eratag = eradict[era]

    axis_eta_eff = None
    axis_pt_eff = None
    # categorical axes in python bindings always have an overflow bin, so use a regular
    # axis for the charge
    axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "SF charge")
    axis_charge_inclusive = hist.axis.Regular(1, -2., 2., underflow=False, overflow=False, name = "SF charge") # for isolation and effStat only
    isoEff_types = ["iso", "isonotrig", "antiiso", "antiisonotrig"]
    allEff_types = ["reco", "tracking", "idip", "trigger"] + isoEff_types
    #eff_types_3D = ["trigger"] + [x for x in isoEff_types] # if this is set to [] the code should revert back to all 2D SF
    eff_types_3D = [] if not smooth3D else ["trigger"] + [x for x in isoEff_types]
    eff_types_2D = [x for x in allEff_types if x not in eff_types_3D]
    axis_allEff_type = hist.axis.StrCategory(allEff_types, name = "allEff_type")
    axis_eff_type_2D = hist.axis.StrCategory(eff_types_2D, name = "eff_types_2D_etapt")
    axis_eff_type_3D = hist.axis.StrCategory(eff_types_3D, name = "eff_types_3D_etaptut")
    effSyst_decorrEtaEdges = [round(-2.4 + 0.1*i,1) for i in range(49)]
    Nsyst = 1 + (len(effSyst_decorrEtaEdges) - 1) # inclusive variation + all decorrelated bins
    axis_nom_syst = hist.axis.Integer(0, 1+Nsyst, underflow = False, overflow =False, name = "nom-systs") # nominal in first bin

    charges = { -1. : "minus", 1. : "plus" }
    chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
  
    fin = ROOT.TFile.Open(filename)
    if fin is None or fin.IsZombie():
        raise IOError(f"Error: file {filename} was not opened correctly")
        
    ## start with NOMI and SYST
    sf_syst_2D = None
    ############
    # all ut-independent steps, for the nominal and systematics
    for charge, charge_tag in charges.items():
        for eff_type in axis_eff_type_2D:
            # for iso can use histogram for efficiency variation only in data (the only one we have for now)
            nameTag = "nomiAndAlt_onlyDataVar" if (isoEfficiencySmoothing and any(x in eff_type for x in ["iso", "antiiso"])) else "nomiAndAlt"
            chargeTag = charge_tag if eff_type in chargeDependentSteps else "both"
            hist_name = f"SF_{nameTag}_{eratag}_{eff_type}_{chargeTag}"
            hist_root = fin.Get(hist_name)
            if hist_root is None:
                logger.info(f"Error: {hist_name} not found in file {filename}")
                quit()
            logger.debug(f"syst: {eff_type} -> {hist_name}")

            hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt", "nomi-statUpDown-syst"])

            if sf_syst_2D is None:
                axis_eta_eff = hist_hist.axes[0]
                axis_pt_eff = hist_hist.axes[1]
                # store all systs (currently only 1) with the nominal, for all efficiency steps
                sf_syst_2D = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge, axis_eff_type_2D, axis_nom_syst, name = "sf_syst_2D", storage = hist.storage.Weight())
            # this axis might change for different histograms, because of a different number of effStat variations
            axis_nomiAlt_eff = hist_hist.axes[2]
            # could use max_pt to remove some of the pt bins for the input histogram
            # extract nominal (first bin that is not underflow) and put in corresponding bin of destination (bin 0 is the first bin because no underflow)
            sf_syst_2D.view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type_2D.index(eff_type), 0] = hist_hist.view(flow=True)[:,:,1]
            # extract syst (last bin except overflow) and put in corresponding bin of destination (bin 1 is the second bin because no underflow)
            sf_syst_2D.view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type_2D.index(eff_type), 1] = hist_hist.view(flow=True)[:,:,axis_nomiAlt_eff.extent-2]
            for isyst in range(len(effSyst_decorrEtaEdges)-1):
                # first copy the nominal
                sf_syst_2D.view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type_2D.index(eff_type), 2+isyst] = hist_hist.view(flow=True)[:,:,1]
                # now update with actual syst all eta bins inside interval [effSyst_decorrEtaEdges[isyst], effSyst_decorrEtaEdges[isyst+1]]
                # add epsilon to ensure picking the bin on the right of the edge (for the right edge given by
                # effSyst_decorrEtaEdges[isyst+1]] the range selection in boost later on will stop at the left
                #edge of the chosen bin number, e.g. h[b:b+1] will pick the range containing the single bin b, unlike in ROOT
                # also sum 1 because sf_syst_2D.view(flow=True) will always put the underflow in the 0 bin index
                indexEtaLow = 1 + axis_eta_eff.index(effSyst_decorrEtaEdges[isyst] + 0.001) # add epsilon to ensure picking the bin on the right of the edge
                indexEtaHigh = 1 + axis_eta_eff.index(effSyst_decorrEtaEdges[isyst+1] + 0.001) 
                sf_syst_2D.view(flow=True)[indexEtaLow:indexEtaHigh, :, axis_charge.index(charge), axis_eff_type_2D.index(eff_type), 2+isyst] = hist_hist.view(flow=True)[indexEtaLow:indexEtaHigh, :,axis_nomiAlt_eff.extent-2]
                        
    # set overflow and underflow eta-pt bins equal to adjacent bins
    sf_syst_2D.view(flow=True)[0, ...] = sf_syst_2D.view(flow=True)[1, ...]
    sf_syst_2D.view(flow=True)[axis_eta_eff.extent-1, ...] = sf_syst_2D.view(flow=True)[axis_eta_eff.extent-2, ...]
    sf_syst_2D.view(flow=True)[:, 0, ...] = sf_syst_2D.view(flow=True)[:, 1, ...]
    sf_syst_2D.view(flow=True)[:, axis_pt_eff.extent-1, ...] = sf_syst_2D.view(flow=True)[:, axis_pt_eff.extent-2, ...]

    ## now proceed with loading 3D SF if existing, then create helpers as appropriate
    if len(eff_types_3D):

        logger.warning("Using SF in 3D for trigger/isolation")
        dict_SF3D = None
        # temporary file stored locally for tests
        fileSF3D = data_dir + "/testMuonSF/smoothSF3D.pkl.lz4"
        with lz4.frame.open(fileSF3D) as f3D:
            dict_SF3D = pickle.load(f3D)

        sf_syst_3D = None
        ############
        # all ut-dependent steps, for the nominal and systematics
        for charge, charge_tag in charges.items():
            for eff_type in axis_eff_type_3D:
                effTypeNameInHist = f"{eff_type}{charge_tag}" if eff_type == "trigger" else eff_type
                if "anti" in effTypeNameInHist:
                    effTypeNameInHist = effTypeNameInHist.replace("anti", "") # FIXME temporary until we have all
                hist_hist = dict_SF3D[f"smoothSF3D_{effTypeNameInHist}"]
                #logger.debug(f"{charge} {eff_type}: hist_hist = {hist_hist.axes}")
                if sf_syst_3D is None:
                    ## Redefine axes to force the existence of overflow bins
                    axis_eta_eff = cloneAxis(hist_hist.axes[0], overflow=True, underflow=True)
                    axis_pt_eff = cloneAxis(hist_hist.axes[1], overflow=True, underflow=True)
                    axis_ut_eff = cloneAxis(hist_hist.axes[2], overflow=True, underflow=True)
                    # ut axis is put last to keep same axes sorting as for 2D histograms created before
                    sf_syst_3D = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge,
                                           axis_eff_type_3D, axis_nom_syst, axis_ut_eff,
                                           name = "sf_syst_3D", storage = hist.storage.Weight())
                # could use max_pt to remove some of the pt bins for the input histogram
                # extract nominal (first bin that is not underflow) and put in corresponding bin of destination (bin 0 is the first bin because no underflow)
                ## Note: must call sf_syst_3D with flow=False because it has overflow bins but hist_hist does not
                sf_syst_3D.view(flow=False)[:, :, axis_charge.index(charge), axis_eff_type_3D.index(eff_type), 0, :] = hist_hist.view(flow=False)[:,:,:,0]
                # extract syst (last bin except overflow) and put in corresponding bin of destination (bin 1 is the second bin because no underflow)
                ##
                ## FIXME: for now using syst=nominal since we don't have syst (must copy from 2D version)
                ##
                sf_syst_3D.view(flow=False)[:, :, axis_charge.index(charge), axis_eff_type_3D.index(eff_type), 1, :] = hist_hist.view(flow=False)[:,:,:,0]  ## <---- FIX THIS !!!
                for isyst in range(len(effSyst_decorrEtaEdges)-1):
                    # first copy the nominal
                    sf_syst_3D.view(flow=False)[:, :, axis_charge.index(charge), axis_eff_type_3D.index(eff_type), 2+isyst, :] = hist_hist.view(flow=False)[:,:,:,0]
                    # now update with actual syst all eta bins inside interval [effSyst_decorrEtaEdges[isyst], effSyst_decorrEtaEdges[isyst+1]]
                    # add epsilon to ensure picking the bin on the right of the edge (for the right edge given by
                    # effSyst_decorrEtaEdges[isyst+1]] the range selection in boost later on will stop at the left
                    #edge of the chosen bin number, e.g. h[b:b+1] will pick the range containing the single bin b, unlike in ROOT
                    # also do not sum 1 because sf_syst_3D.view(flow=False) will consider 0 the first bin index (with flow=True instead 0 is the underflow but only if it exists, otherwise 0 is the first bin)
                    indexEtaLow = axis_eta_eff.index(effSyst_decorrEtaEdges[isyst] + 0.001) # add epsilon to ensure picking the bin on the right of the edge
                    indexEtaHigh = axis_eta_eff.index(effSyst_decorrEtaEdges[isyst+1] + 0.001) 
                    sf_syst_3D.view(flow=False)[indexEtaLow:indexEtaHigh, :, axis_charge.index(charge), axis_eff_type_3D.index(eff_type), 2+isyst, :] = hist_hist.view(flow=False)[indexEtaLow:indexEtaHigh, :, :, 0]
                    
        # set overflow and underflow eta-pt bins equal to adjacent bins
        sf_syst_3D.view(flow=True)[0, ...]                       = sf_syst_3D.view(flow=True)[1, ...]
        sf_syst_3D.view(flow=True)[axis_eta_eff.extent-1, ...]   = sf_syst_3D.view(flow=True)[axis_eta_eff.extent-2, ...]
        sf_syst_3D.view(flow=True)[:, 0, ...]                    = sf_syst_3D.view(flow=True)[:, 1, ...]
        sf_syst_3D.view(flow=True)[:, axis_pt_eff.extent-1, ...] = sf_syst_3D.view(flow=True)[:, axis_pt_eff.extent-2, ...]
        sf_syst_3D.view(flow=True)[..., 0]                       = sf_syst_3D.view(flow=True)[..., 1]
        sf_syst_3D.view(flow=True)[..., axis_ut_eff.extent-1]    = sf_syst_3D.view(flow=True)[..., axis_ut_eff.extent-2]

        #logger.error(f"axis_pt_eff.extent = {axis_pt_eff.extent}")
        #logger.error(f"axis_eta_eff.extent = {axis_eta_eff.extent}")
        #logger.error(f"axis_ut_eff.extent = {axis_ut_eff.extent}")
        ## try rebinning
        rebinUt = True
        if rebinUt:
            rebin = 60
            logger.error(f"Attention, rebinning ut by {rebin} as a test")
            sf_syst_3D = sf_syst_3D[{axis_ut_eff.name : hist.rebin(rebin)}] # try with and without rebinning
        logger.debug("")
        logger.debug(f"sf_syst_3D.shape = {sf_syst_3D.shape}") # this is currently (48, 205, 2, 5, 50, N_ut) for eta,pt,charge,effType,nomiAndSyst,ut
        #logger.debug("")
        #logger.debug(f"sf_syst_2D.axes = {sf_syst_2D.axes}")
        logger.debug("")
        logger.debug(f"sf_syst_3D.axes = {sf_syst_3D.axes}")
        
        sf_syst_2D_pyroot = narf.hist_to_pyroot_boost(sf_syst_2D)
        sf_syst_3D_pyroot = narf.hist_to_pyroot_boost(sf_syst_3D)
        # nomi and syst are stored in the same histogram, just use different helpers to override the () operator for now, until RDF is improved
        helper = ROOT.wrem.muon_efficiency_smooth3D_helper[str(is_w_like).lower(), Nsyst,
                                                           type(sf_syst_2D_pyroot),
                                                           type(sf_syst_3D_pyroot)](
                                                               ROOT.std.move(sf_syst_2D_pyroot),
                                                               ROOT.std.move(sf_syst_3D_pyroot)
                                                           )
        helper_syst = ROOT.wrem.muon_efficiency_smooth3D_helper_syst[str(is_w_like).lower(), Nsyst, type(sf_syst_2D_pyroot), type(sf_syst_3D_pyroot)]( helper )
        # define axis for syst variations with all steps
        axis_all = hist.axis.Integer(0, 5, underflow = False, overflow = False, name = "reco-tracking-idip-trigger-iso")
        axis_nsyst = hist.axis.Integer(0, Nsyst, underflow = False, overflow = False, name = "n_syst_variations")
        helper_syst.tensor_axes = [axis_all, axis_nsyst]
        #

    else:
        logger.warning("Using all SF in 2D")
        # case with only 2D histograms
        sf_syst_2D_pyroot = narf.hist_to_pyroot_boost(sf_syst_2D)
        # nomi and syst are stored in the same histogram, just use different helpers to override the () operator for now, until RDF is improved
        helper = ROOT.wrem.muon_efficiency_smooth_helper[str(is_w_like).lower(), Nsyst, type(sf_syst_2D_pyroot)]( ROOT.std.move(sf_syst_2D_pyroot) )
        helper_syst = ROOT.wrem.muon_efficiency_smooth_helper_syst[str(is_w_like).lower(), Nsyst, type(sf_syst_2D_pyroot)]( helper )
        # define axis for syst variations with all steps
        axis_all = hist.axis.Integer(0, 5, underflow = False, overflow = False, name = "reco-tracking-idip-trigger-iso")
        axis_nsyst = hist.axis.Integer(0, Nsyst, underflow = False, overflow = False, name = "n_syst_variations")
        helper_syst.tensor_axes = [axis_all, axis_nsyst]
        #
                       
    ##############
    ## now the EFFSTAT
    
    # for the stat histograms, an axis will contain the efficiency type, in some cases it might have a single bin (e.g. tracking, trigger)
    # but for isolation needs more, because of the different cases to handle inside the helper class (iso-antiiso with/without trigger)
    # could actually stack reco-idip-trig together, but better to have the helper deal with a single histogram at a time
    # then there will be a specialization of the class for isolation
    # in this way the name of the (single) effType axis can be used to activate some internal behaviour as needed (for example trigger won't have variations for not triggering lepton)
    
    effStat_manager = {"sf_reco": {"nPtEigenBins" : None,
                                   "nCharges" : None,
                                   "axisLabels" : ["reco"],
                                   "boostHist" : None,
                                   "helper" : None},
                       "sf_tracking": {"nPtEigenBins" : None,
                                       "nCharges" : None,
                                       "axisLabels" : ["tracking"],
                                       "boostHist" : None,
                                       "helper" : None},
                       "sf_idip": {"nPtEigenBins" : None,
                                   "nCharges" : None,
                                   "axisLabels" : ["idip"],
                                   "boostHist" : None,
                                   "helper" : None},
                       "sf_trigger": {"nPtEigenBins" : None,
                                      "nCharges" : None,
                                      "axisLabels" : ["trigger"],
                                      "boostHist" : None,
                                      "helper" : None},
                       }
    if not isoEfficiencySmoothing:
        effStat_manager["sf_iso"] = {"nPtEigenBins" : None,
                                     "nCharges" : None,
                                     "axisLabels" : ["iso", "isonotrig", "antiiso", "antiisonotrig"],
                                     "boostHist" : None,
                                     "helper" : None}
    else:
        effStat_manager["sf_iso_effData"] =  {"nPtEigenBins" : None,
                                              "nCharges" : None,
                                              "axisLabels" : ["iso", "isonotrig", "antiiso", "antiisonotrig"],
                                              "boostHist" : None,
                                              "helper" : None}
        effStat_manager["sf_iso_effMC"] = {"nPtEigenBins" : None,
                                           "nCharges" : None,
                                           "axisLabels" : ["iso", "isonotrig", "antiiso", "antiisonotrig"],
                                           "boostHist" : None,
                                           "helper" : None}
        
    for effStatKey in effStat_manager.keys():
        nom_up_effStat_axis = None
        axis_eff_type = None
        axis_charge_def = None
        for ic, (charge, charge_tag) in enumerate(charges.items()):
            for eff_type in effStat_manager[effStatKey]["axisLabels"]:
                if eff_type in chargeDependentSteps:
                    axis_charge_def = axis_charge
                    chargeTag = charge_tag
                else:
                    axis_charge_def = axis_charge_inclusive
                    chargeTag = "both"                    
                    if ic: continue # must continue after having set the charge axis
                nameTag = "nomiAndAlt"
                if "effData" in effStatKey:
                    nameTag += "_onlyDataVar"
                elif "effMC" in effStatKey:
                    nameTag += "_onlyMCVar"
                hist_name = f"SF_{nameTag}_{eratag}_{eff_type}_{chargeTag}"
                hist_root = fin.Get(hist_name)
                # logger.info(f"stat: {effStatKey}|{eff_type} -> {hist_name}")
                # logger.info(ROOT.AddressOf(hist_root))
                if hist_root is None:
                    logger.info(f"Error: {hist_name} not found in file {filename}")
                    quit()

                if nom_up_effStat_axis is None:
                    # Z axis has nominal, 2*N effStat (first N up then N down) and one syst in the last bin
                    nPtEigenBins = int(int(hist_root.GetNbinsZ() - 2)/ 2)
                    effStat_manager[effStatKey]["nPtEigenBins"] = nPtEigenBins
                    effStat_manager[effStatKey]["nCharges"] = 2 if eff_type in chargeDependentSteps else 1
                    nom_up_effStat_axis = hist.axis.Regular(int(1+nPtEigenBins), -0.5, 0.5+nPtEigenBins, underflow=False, overflow=False, name = "nomUpVar")

                hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt", "nomi-statUpDown-syst"])

                if effStat_manager[effStatKey]["boostHist"] is None:
                    axis_eta_eff = hist_hist.axes[0]
                    axis_pt_eff = hist_hist.axes[1]
                    axis_eff_type = hist.axis.StrCategory(effStat_manager[effStatKey]["axisLabels"], name = f"{effStatKey}_eff_type")
                    effStat_manager[effStatKey]["boostHist"] = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge_def,
                                                                         axis_eff_type,
                                                                         nom_up_effStat_axis,
                                                                         name = effStatKey,
                                                                         storage = hist.storage.Weight())
                    
                # extract nominal (first bin that is not underflow)
                effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, :, axis_charge_def.index(charge), axis_eff_type.index(eff_type), nom_up_effStat_axis.index(0)] = hist_hist.view(flow=True)[:,:,1]
                # now extract the stat variations
                # note that for hist_hist the nominal bin of the third axis is number 1 when the underflow is active
                for iup in range(1, 1 + nPtEigenBins):
                    # up variations (down are going to be created later by mirroring in CardTool.py
                    effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, :, axis_charge_def.index(charge), axis_eff_type.index(eff_type), nom_up_effStat_axis.index(iup)] = hist_hist.view(flow=True)[:,:,1+iup]
                
        # set overflow and underflow equal to adjacent bins
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[0, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[1, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[axis_eta_eff.extent-1, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[axis_eta_eff.extent-2, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, 0, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, 1, ...]
        effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, axis_pt_eff.extent-1, ...] = effStat_manager[effStatKey]["boostHist"].view(flow=True)[:, axis_pt_eff.extent-2, ...]

        netabins = axis_eta_eff.size
        ncharges = effStat_manager[effStatKey]["nCharges"]
        sf_stat_pyroot = narf.hist_to_pyroot_boost(effStat_manager[effStatKey]["boostHist"])
        if "sf_iso" in effStatKey:
            helper_stat = ROOT.wrem.muon_efficiency_smooth_helper_stat_iso[str(is_w_like).lower(), netabins, nPtEigenBins, ncharges, type(sf_stat_pyroot)]( ROOT.std.move(sf_stat_pyroot) )
        else:
            helper_stat = ROOT.wrem.muon_efficiency_smooth_helper_stat[str(is_w_like).lower(), netabins, nPtEigenBins, ncharges, type(sf_stat_pyroot)]( ROOT.std.move(sf_stat_pyroot) )
        # make new versions of these axes without overflow/underflow to index the tensor
        if isinstance(axis_eta_eff, bh.axis.Regular):
            axis_eta_eff_tensor = hist.axis.Regular(axis_eta_eff.size, axis_eta_eff.edges[0], axis_eta_eff.edges[-1], name = axis_eta_eff.name, overflow = False, underflow = False)
        elif isinstance(axis_eta_eff, bh.axis.Variable):
            axis_eta_eff_tensor = hist.axis.Variable(axis_eta_eff.edges, name = axis_eta_eff.name, overflow = False, underflow = False)
        axis_ptEigen_eff_tensor = hist.axis.Integer(0, effStat_manager[effStatKey]["nPtEigenBins"], underflow = False, overflow =False, name = "nPtEigenBins")    
        effStatTensorAxes = [axis_eta_eff_tensor, axis_ptEigen_eff_tensor, axis_charge_def]
        helper_stat.tensor_axes = effStatTensorAxes
        effStat_manager[effStatKey]["helper"] = helper_stat

    fin.Close()

    logger.debug(f"Return efficiency helper!")

    ####
    # return nomi, effsyst, and a dictionary with effStat to use them by name
    return helper, helper_syst, {k : effStat_manager[k]["helper"] for k in effStat_manager.keys()}
