import logging
import lz4.frame
import pickle
import hdf5plugin
import h5py
import narf
import ROOT

#FIRST STEP IN EFFICIENCY CREATION

name = "rnf"

histlist = []

directory = "directory_name"

histlist.append(f"{directory}/mw_with_mu_eta_pt_scetlibCorr_idip{name}.hdf5")
histlist.append(f"{directory}/mw_with_mu_eta_pt_scetlibCorr_trigger{name}.hdf5")
histlist.append(f"{directory}/mw_with_mu_eta_pt_scetlibCorr_iso{name}.hdf5")

listoflists = [["nominal"],["nominal_smoothMC","nominal","nominal_smoothErr"],["nominal_smoothMC","nominal","nominal_smoothErr"]] 

namelist = []
namelist.append("IDIP")
namelist.append("TriggerMC")
namelist.append("Trigger")
namelist.append("TriggerError")
namelist.append("IsoMC")
namelist.append("Iso")
namelist.append("IsoError")

histos = []

counter = 0
for idx,filename in enumerate(histlist) :
    h5file = h5py.File(filename, "r")
    results = narf.ioutils.pickle_load_h5py(h5file["results"])
    plus = results["WplusmunuPostVFP"]["output"]
    minus = results["WminusmunuPostVFP"]["output"]
    for i in range(0,len(listoflists[idx])):
        plus = plus[listoflists[idx][i]].get()
        minus = minus[listoflists[idx][i]].get()
        if "Error" in namelist[counter]:
            plus = plus[:,:,1:2,:,:,:]
            plus = plus[:,:,sum,sum,sum,:]
            minus = minus[:,:,0:1,:,:,:]
            minus = minus[:,:,sum,sum,sum,:]
        else :
            plus = plus[:,:,1:2,:,:]
            plus = plus[:,:,sum,sum,sum]
            minus = minus[:,:,0:1,:,:]
            minus = minus[:,:,sum,sum,sum]
        plus=narf.hist_to_root(plus)
        plus.SetName(f"{namelist[counter]}Plus")
        minus=narf.hist_to_root(minus)
        minus.SetName(f"{namelist[counter]}Minus")
        histos.append(plus)
        histos.append(minus)
        counter = counter + 1

output_file = ROOT.TFile(f"makeefficiencies{name}.root","RECREATE")
output_file.cd()

for histo in histos:
    histo.Write()
output_file.Close()
