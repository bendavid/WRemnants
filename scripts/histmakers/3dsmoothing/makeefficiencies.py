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
        Plus = plus[listoflists[idx][i]].get()
        Minus = minus[listoflists[idx][i]].get()
        if "Error" in namelist[counter]:
            Plus = Plus[:,:,1:2,:,:,:]
            Plus = Plus[:,:,sum,sum,sum,:]
            Minus = Minus[:,:,0:1,:,:,:]
            Minus = Minus[:,:,sum,sum,sum,:]
        else :
            Plus = Plus[:,:,1:2,:,:]
            Plus = Plus[:,:,sum,sum,sum]
            Minus = Minus[:,:,0:1,:,:]
            Minus = Minus[:,:,sum,sum,sum]
        Plus=narf.hist_to_root(Plus)
        Plus.SetName(f"{namelist[counter]}Plus")
        Minus=narf.hist_to_root(Minus)
        Minus.SetName(f"{namelist[counter]}Minus")
        histos.append(Plus)
        histos.append(Minus)
        counter = counter + 1

output_file = ROOT.TFile(f"makeefficiencies{name}.root","RECREATE")
output_file.cd()

for histo in histos:
    histo.Write()
output_file.Close()
