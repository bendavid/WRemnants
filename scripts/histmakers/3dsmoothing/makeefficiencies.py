import logging
import lz4.frame
import pickle
import hdf5plugin
import h5py
import narf
import ROOT

#FIRST STEP IN EFFICIENCY CREATION

name = "rnd"

histlist = []

histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_idip{name}.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_triggerMC{name}.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_trigger{name}.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_trigger{name}error.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_isoMC{name}.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_iso{name}.hdf5")
histlist.append(f"mw_with_mu_eta_pt_scetlibCorr_iso{name}error.hdf5")

namelist = []
namelist.append("IDIP")
namelist.append("TriggerMC")
namelist.append("Trigger")
namelist.append("TriggerError")
namelist.append("IsoMC")
namelist.append("Iso")
namelist.append("IsoError")

histos = []

for idx,filename in enumerate(histlist) :
    h5file = h5py.File(filename, "r")
    results = narf.ioutils.pickle_load_h5py(h5file["results"]))
    plus = results["WplusmunuPostVFP"]["output"]
    minus = results["WminusmunuPostVFP"]["output"]
    plus = plus["nominal"]
    minus = minus["nominal"]
    if "error" in filename:
        plus = plus[:,:,1:2,:,:,:]
        plus = plus[:,:,sum,sum,sum,:]
        minus = minus[:,:,1:2,:,:,:]
        minus = minus[:,:,sum,sum,sum,:]
    else :
        plus = plus[:,:,1:2,:,:]
        plus = plus[:,:,sum,sum,sum]
        minus = minus[:,:,1:2,:,:]
        minus = minus[:,:,sum,sum,sum]
    plus=narf.hist_to_root(plus)
    plus.SetName(f"{namelist[idx]}Plus")
    minus=narf.hist_to_root(minus)
    minus.SetName(f"{namelist[idx]}Minus")
    histos.append(plus)
    histos.append(minus)

output_file = ROOT.TFile("makeefficiencies{name}.root","RECREATE")
output_file.cd()

for histo in histos:
    histo.Write()
output_file.Close()
