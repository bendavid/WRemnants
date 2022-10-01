import ROOT
import pathlib
import hist
import narf
import numpy as np
import boost_histogram as bh

ROOT.gInterpreter.Declare('#include "muon_efficiencies.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def make_muon_efficiency_helpers(filename = data_dir + "/testMuonSF/scaleFactorProduct_29Sept2022_vertexWeight_OSchargeExceptTracking.root", era = None, is_w_like = False, max_pt = np.inf):

    eradict = { "2016PreVFP" : "BtoF",
                "2016PostVFP" : "GtoH" }
    eratag = eradict[era]

    axis_eta_eff = None
    axis_pt_eff = None
    # categorical axes in python bindings always have an overflow bin, so use a regular
    # axis for the charge
    axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "SF charge")
    axis_eff_type = hist.axis.StrCategory(["idip", "idip_trig", "iso_triggering", "antiiso_triggering", "iso_nontriggering", "antiiso_nontriggering"], name = "eff_type")
    axis_with_trigger = hist.axis.Boolean(name = "with_trigger")
    axis_nom_alt = hist.axis.Integer(0, 2, underflow = False, overflow =False, name = "nom-alt")

    sf_idip_trig_iso = None

    charges = { -1. : "Minus", 1. : "Plus" }
    eff_types = { "idip" : "noisoNotrig",
                 "idip_trig" : "noisoTrig",
                 "iso_triggering" : "isoOnly",
                 "antiiso_triggering" : "antiisoOnly",
                 "iso_nontriggering" : "isoNotrigOnly",
                 "antiiso_nontriggering" : "antiisoNotrigOnly",
                 }
    nom_systs = { 0 : "nominal", 1 : "dataAltSig" }

    fin = ROOT.TFile.Open(filename);

    # all the (idip_trig_iso) SF are belong to us
    for charge, charge_tag in charges.items():
        for eff_type, eff_type_tag in eff_types.items():
            for nom_syst, nom_syst_tag in nom_systs.items():
                hist_name = f"fullSF2D_{nom_syst_tag}_{eff_type_tag}{charge_tag}_{eratag}"
                hist_root = fin.Get(hist_name)
                #print(ROOT.AddressOf(hist_root))
                # might not have a charge specific version
                if not isinstance(hist_root, ROOT.TH1):
                    hist_name = f"fullSF2D_{nom_syst_tag}_{eff_type_tag}_{eratag}"
                    hist_root = fin.Get(hist_name)

                hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt"])

                if sf_idip_trig_iso is None:
                    axis_eta_eff = hist_hist.axes[0]
                    axis_pt_eff = hist_hist.axes[1]

                    sf_idip_trig_iso = hist.Hist(axis_eta_eff, axis_pt_eff, axis_charge, axis_eff_type, axis_nom_alt,  name = "sf_idip_trig_iso", storage = hist.storage.Weight())

                sf_idip_trig_iso.view(flow=True)[:, :, axis_charge.index(charge), axis_eff_type.index(eff_type), axis_nom_alt.index(nom_syst)] = hist_hist.view(flow=True)[...]

    # set overflow and underflow equal to adjacent bins
    sf_idip_trig_iso.view(flow=True)[0, ...] = sf_idip_trig_iso.view(flow=True)[1, ...]
    sf_idip_trig_iso.view(flow=True)[axis_eta_eff.extent-1, ...] = sf_idip_trig_iso.view(flow=True)[axis_eta_eff.extent-2, ...]
    sf_idip_trig_iso.view(flow=True)[:, 0, ...] = sf_idip_trig_iso.view(flow=True)[:, 1, ...]
    sf_idip_trig_iso.view(flow=True)[:, axis_pt_eff.extent-1, ...] = sf_idip_trig_iso.view(flow=True)[:, axis_pt_eff.extent-2, ...]

    sf_tracking = None
    sf_reco = None

    # tracking and reco SFs
    for nom_syst, nom_syst_tag in nom_systs.items():
        ## first tracking
        hist_name = f"fullSF2D_{nom_syst_tag}_tracking_{eratag}"
        hist_root = fin.Get(hist_name)
        hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt"])

        if sf_tracking is None:
            axis_eta_eff_tracking = hist_hist.axes[0]
            axis_pt_eff_tracking = hist_hist.axes[1]
            sf_tracking = hist.Hist(axis_eta_eff_tracking, axis_pt_eff_tracking, axis_charge, axis_nom_alt,  name = "sf_tracking", storage = hist.storage.Weight())

        # same SFs for + and - for now
        sf_tracking.view(flow=True)[:, :, axis_charge.index(-1.), axis_nom_alt.index(nom_syst)] = hist_hist.view(flow=True)[...]
        sf_tracking.view(flow=True)[:, :, axis_charge.index(1.), axis_nom_alt.index(nom_syst)] = hist_hist.view(flow=True)[...]

        ## FIXME: avoid code duplication?
        ## now again for reco
        hist_name = f"fullSF2D_{nom_syst_tag}_reco_{eratag}"
        hist_root = fin.Get(hist_name)
        hist_hist = narf.root_to_hist(hist_root, axis_names = ["SF eta", "SF pt"])

        if sf_reco is None:
            axis_eta_eff_reco = hist_hist.axes[0]
            axis_pt_eff_reco = hist_hist.axes[1]
            sf_reco = hist.Hist(axis_eta_eff_reco, axis_pt_eff_reco, axis_charge, axis_nom_alt,  name = "sf_reco", storage = hist.storage.Weight())

        # same SFs for + and - for now
        sf_reco.view(flow=True)[:, :, axis_charge.index(-1.), axis_nom_alt.index(nom_syst)] = hist_hist.view(flow=True)[...]
        sf_reco.view(flow=True)[:, :, axis_charge.index(1.), axis_nom_alt.index(nom_syst)] = hist_hist.view(flow=True)[...]

        
    # set overflow and underflow equal to adjacent bins
    sf_tracking.view(flow=True)[0, ...] = sf_tracking.view(flow=True)[1, ...]
    sf_tracking.view(flow=True)[axis_eta_eff_tracking.extent-1, ...] = sf_tracking.view(flow=True)[axis_eta_eff_tracking.extent-2, ...]
    sf_tracking.view(flow=True)[:, 0, ...] = sf_tracking.view(flow=True)[:, 1, ...]
    sf_tracking.view(flow=True)[:, axis_pt_eff_tracking.extent-1, ...] = sf_tracking.view(flow=True)[:, axis_pt_eff_tracking.extent-2, ...]
    # set overflow and underflow equal to adjacent bins
    sf_reco.view(flow=True)[0, ...] = sf_reco.view(flow=True)[1, ...]
    sf_reco.view(flow=True)[axis_eta_eff_reco.extent-1, ...] = sf_reco.view(flow=True)[axis_eta_eff_reco.extent-2, ...]
    sf_reco.view(flow=True)[:, 0, ...] = sf_reco.view(flow=True)[:, 1, ...]
    sf_reco.view(flow=True)[:, axis_pt_eff_reco.extent-1, ...] = sf_reco.view(flow=True)[:, axis_pt_eff_reco.extent-2, ...]

    
    fin.Close()

    # TODO implement a convenient way to use enums for the bool/category/nom-var axes such that
    # the C++ code is less error-prone?

    netabins = axis_eta_eff.size  # this is the same also for reco and tracking, should one check explicitly?

    # exclude pt bins outside of analysis range
    nptbins = 0
    nptbins_tracking = 0
    nptbins_reco = 0
    for edge in axis_pt_eff.edges:
        if edge < max_pt:
            nptbins += 1
    for edge in axis_pt_eff_tracking.edges:
        if edge < max_pt:
            nptbins_tracking += 1
    for edge in axis_pt_eff_reco.edges:
        if edge < max_pt:
            nptbins_reco += 1

    sf_idip_trig_iso_pyroot = narf.hist_to_pyroot_boost(sf_idip_trig_iso)
    sf_tracking_pyroot = narf.hist_to_pyroot_boost(sf_tracking)
    sf_reco_pyroot = narf.hist_to_pyroot_boost(sf_reco)

    helper = ROOT.wrem.muon_efficiency_helper[str(is_w_like).lower(), netabins, nptbins,
                                              type(sf_idip_trig_iso_pyroot),
                                              type(sf_tracking_pyroot),
                                              type(sf_reco_pyroot)](
                                                  ROOT.std.move(sf_idip_trig_iso_pyroot),
                                                  ROOT.std.move(sf_tracking_pyroot),
                                                  ROOT.std.move(sf_reco_pyroot)
                                              )

    helper_stat = ROOT.wrem.muon_efficiency_helper_stat[str(is_w_like).lower(), netabins, nptbins, type(sf_idip_trig_iso_pyroot), type(sf_tracking_pyroot), type(sf_reco_pyroot)](helper)
    helper_stat_tracking = ROOT.wrem.muon_efficiency_helper_singleStep_stat[str(is_w_like).lower(), netabins, nptbins_tracking, type(sf_idip_trig_iso_pyroot), type(sf_tracking_pyroot), type(sf_reco_pyroot)](helper)
    helper_stat_reco     = ROOT.wrem.muon_efficiency_helper_singleStep_stat[str(is_w_like).lower(), netabins, nptbins_reco,     type(sf_idip_trig_iso_pyroot), type(sf_tracking_pyroot), type(sf_reco_pyroot)](helper)

    #TODO make this a categorical once we can disable overflow
    axis_idiptrig_iso = hist.axis.Integer(0, 2, underflow = False, overflow = False, name = "idiptrig-iso") # for stat variations with idiptrig and iso
    axis_reco_tracking_idiptrig_iso = hist.axis.Integer(0, 4, underflow = False, overflow = False, name = "reco-tracking-idiptrig-iso") # for syst variations with all steps (idiptrig is one single piece)

    # make new versions of these axes without overflow/underflow to index the tensor
    if isinstance(axis_eta_eff, bh.axis.Regular):
        axis_eta_eff_tensor = hist.axis.Regular(axis_eta_eff.size, axis_eta_eff.edges[0], axis_eta_eff.edges[-1], name = axis_eta_eff.name, overflow = False, underflow = False)
    elif isinstance(axis_eta_eff, bh.axis.Variable):
        axis_eta_eff_tensor = hist.axis.Variable(axis_eta_eff.size, axis_eta_eff.edges, name = axis_eta_eff.name, overflow = False, underflow = False)

    # for pt need to additionally remove the out of range bins if any
    if isinstance(axis_pt_eff, bh.axis.Regular):
        axis_pt_eff_tensor = hist.axis.Regular(nptbins, axis_pt_eff.edges[0], axis_pt_eff.edges[nptbins], name = axis_pt_eff.name, overflow = False, underflow = False)
    elif isinstance(axis_pt_eff, bh.axis.Variable):
        axis_pt_eff_tensor = hist.axis.Variable(axis_pt_eff.edges[:nptbins+1], name = axis_pt_eff.name, overflow = False, underflow = False)
    # repeat for tracking
    if isinstance(axis_pt_eff_tracking, bh.axis.Regular):
        axis_pt_eff_tensor_tracking = hist.axis.Regular(nptbins_tracking, axis_pt_eff_tracking.edges[0], axis_pt_eff_tracking.edges[nptbins_tracking], name = axis_pt_eff_tracking.name, overflow = False, underflow = False)
    elif isinstance(axis_pt_eff, bh.axis.Variable):
        axis_pt_eff_tensor_tracking = hist.axis.Variable(axis_pt_eff_tracking.edges[:nptbins_tracking+1], name = axis_pt_eff_tracking.name, overflow = False, underflow = False)
    # now also reco
    if isinstance(axis_pt_eff_reco, bh.axis.Regular):
        axis_pt_eff_tensor_reco = hist.axis.Regular(nptbins_reco, axis_pt_eff_reco.edges[0], axis_pt_eff_reco.edges[nptbins_reco], name = axis_pt_eff_reco.name, overflow = False, underflow = False)
    elif isinstance(axis_pt_eff, bh.axis.Variable):
        axis_pt_eff_tensor_reco = hist.axis.Variable(axis_pt_eff_reco.edges[:nptbins_reco+1], name = axis_pt_eff_reco.name, overflow = False, underflow = False)

    # attach tensor axes to the stat helper
    helper_stat.tensor_axes = [axis_eta_eff_tensor, axis_pt_eff_tensor, axis_charge, axis_idiptrig_iso]
    helper_stat_tracking.tensor_axes = [axis_eta_eff_tensor, axis_pt_eff_tensor_tracking, axis_charge]
    helper_stat_reco.tensor_axes = [axis_eta_eff_tensor, axis_pt_eff_tensor_reco, axis_charge]

    helper_syst = ROOT.wrem.muon_efficiency_helper_syst[str(is_w_like).lower(), netabins, nptbins, type(sf_idip_trig_iso_pyroot), type(sf_tracking_pyroot), type(sf_tracking_pyroot)](helper)
    helper_syst.tensor_axes = [axis_reco_tracking_idiptrig_iso]

    return helper, helper_stat, helper_stat_tracking, helper_stat_reco, helper_syst
