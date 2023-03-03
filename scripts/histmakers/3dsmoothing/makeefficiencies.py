import lz4.frame
import pickle
import uproot
import narf
import ROOT

#FIRST STEP IN EFFICIENCY CREATION

name = "rnd"

infileidip = f"mw_with_mu_eta_pt_scetlibCorr_idip{name}.pkl.lz4"
infiletriggerMC = f"mw_with_mu_eta_pt_scetlibCorr_triggerMC{name}.pkl.lz4"
infiletrigger = f"mw_with_mu_eta_pt_scetlibCorr_trigger{name}.pkl.lz4"
infiletriggererror = f"mw_with_mu_eta_pt_scetlibCorr_trigger{name}error.pkl.lz4"
infileisoMC = f"mw_with_mu_eta_pt_scetlibCorr_isoMC{name}.pkl.lz4"
infileiso = f"mw_with_mu_eta_pt_scetlibCorr_iso{name}.pkl.lz4"
infileisoerror = f"mw_with_mu_eta_pt_scetlibCorr_iso{name}error.pkl.lz4"

with lz4.frame.open(infileidip) as f:
    resultsidip = pickle.load(f)
with lz4.frame.open(infiletriggerMC) as f:
    resultstriggerMC = pickle.load(f)
with lz4.frame.open(infiletrigger) as f:
    resultstrigger = pickle.load(f)
with lz4.frame.open(infiletriggererror) as f:
    resultstriggererror = pickle.load(f)
with lz4.frame.open(infileisoMC) as f:
    resultsisoMC = pickle.load(f)
with lz4.frame.open(infileiso) as f:
    resultsiso = pickle.load(f)
with lz4.frame.open(infileisoerror) as f:
    resultsisoerror = pickle.load(f)

outputidipplus = resultsidip["WplusmunuPostVFP"]["output"]
nominalidipplus = outputidipplus["nominal"]
miniidipplus = nominalidipplus[:,:,1:2,:,:]
miniidipplus = miniidipplus[:,:,sum,sum,sum]
outputidipminus = resultsidip["WminusmunuPostVFP"]["output"]
nominalidipminus = outputidipminus["nominal"]
miniidipminus = nominalidipminus[:,:,0:1,:,:]
miniidipminus = miniidipminus[:,:,sum,sum,sum]
outputtriggerMCplus = resultstriggerMC["WplusmunuPostVFP"]["output"]
nominaltriggerMCplus = outputtriggerMCplus["nominal"]
minitriggerMCplus = nominaltriggerMCplus[:,:,1:2,:,:]
minitriggerMCplus = minitriggerMCplus[:,:,sum,sum,sum]
outputtriggerMCminus = resultstriggerMC["WminusmunuPostVFP"]["output"]
nominaltriggerMCminus = outputtriggerMCminus["nominal"]
minitriggerMCminus = nominaltriggerMCminus[:,:,0:1,:,:]
minitriggerMCminus = minitriggerMCminus[:,:,sum,sum,sum]
outputtriggerplus = resultstrigger["WplusmunuPostVFP"]["output"]
nominaltriggerplus = outputtriggerplus["nominal"]
minitriggerplus = nominaltriggerplus[:,:,1:2,:,:]
minitriggerplus = minitriggerplus[:,:,sum,sum,sum]
outputtriggerminus = resultstrigger["WminusmunuPostVFP"]["output"]
nominaltriggerminus = outputtriggerminus["nominal"]
minitriggerminus = nominaltriggerminus[:,:,0:1,:,:]
minitriggerminus = minitriggerminus[:,:,sum,sum,sum]
outputtriggererrorplus = resultstriggererror["WplusmunuPostVFP"]["output"]
nominaltriggererrorplus = outputtriggererrorplus["nominal"]
minitriggererrorplus = nominaltriggererrorplus[:,:,1:2,:,:,:]
minitriggererrorplus = minitriggererrorplus[:,:,sum,sum,sum,:]
outputtriggererrorminus = resultstriggererror["WminusmunuPostVFP"]["output"]
nominaltriggererrorminus = outputtriggererrorminus["nominal"]
minitriggererrorminus = nominaltriggererrorminus[:,:,0:1,:,:,:]
minitriggererrorminus = minitriggererrorminus[:,:,sum,sum,sum,:]
outputisoMCplus = resultsisoMC["WplusmunuPostVFP"]["output"]
nominalisoMCplus = outputisoMCplus["nominal"]
miniisoMCplus = nominalisoMCplus[:,:,1:2,:,:]
miniisoMCplus = miniisoMCplus[:,:,sum,sum,sum]
outputisoMCminus = resultsisoMC["WminusmunuPostVFP"]["output"]
nominalisoMCminus = outputisoMCminus["nominal"]
miniisoMCminus = nominalisoMCminus[:,:,0:1,:,:]
miniisoMCminus = miniisoMCminus[:,:,sum,sum,sum]
outputisoplus = resultsiso["WplusmunuPostVFP"]["output"]
nominalisoplus = outputisoplus["nominal"]
miniisoplus = nominalisoplus[:,:,1:2,:,:]
miniisoplus = miniisoplus[:,:,sum,sum,sum]
outputisominus = resultsiso["WminusmunuPostVFP"]["output"]
nominalisominus = outputisominus["nominal"]
miniisominus = nominalisominus[:,:,0:1,:,:]
miniisominus = miniisominus[:,:,sum,sum,sum]
outputisoerrorplus = resultsisoerror["WplusmunuPostVFP"]["output"]
nominalisoerrorplus = outputisoerrorplus["nominal"]
miniisoerrorplus = nominalisoerrorplus[:,:,1:2,:,:,:]
miniisoerrorplus = miniisoerrorplus[:,:,sum,sum,sum,:]
outputisoerrorminus = resultsisoerror["WminusmunuPostVFP"]["output"]
nominalisoerrorminus = outputisoerrorminus["nominal"]
miniisoerrorminus = nominalisoerrorminus[:,:,0:1,:,:,:]
miniisoerrorminus = miniisoerrorminus[:,:,sum,sum,sum,:]

idipplus=narf.hist_to_root(miniidipplus)
idipplus.SetName("IDIPPlus")
idipminus=narf.hist_to_root(miniidipminus)
idipminus.SetName("IDIPMinus")
triggerMCplus=narf.hist_to_root(minitriggerMCplus)
triggerMCplus.SetName("TriggerMCPlus")
triggerMCminus=narf.hist_to_root(minitriggerMCminus)
triggerMCminus.SetName("TriggerMCMinus")
triggerplus=narf.hist_to_root(minitriggerplus)
triggerplus.SetName("TriggerPlus")
triggerminus=narf.hist_to_root(minitriggerminus)
triggerminus.SetName("TriggerMinus")
triggererrorplus=narf.hist_to_root(minitriggererrorplus)
triggererrorplus.SetName("TriggerErrorPlus")
triggererrorminus=narf.hist_to_root(minitriggererrorminus)
triggererrorminus.SetName("TriggerErrorMinus")
isoMCplus=narf.hist_to_root(miniisoMCplus)
isoMCplus.SetName("IsoMCPlus")
isoMCminus=narf.hist_to_root(miniisoMCminus)
isoMCminus.SetName("IsoMCMinus")
isoplus=narf.hist_to_root(miniisoplus)
isoplus.SetName("IsoPlus")
isominus=narf.hist_to_root(miniisominus)
isominus.SetName("IsoMinus")
isoerrorplus=narf.hist_to_root(miniisoerrorplus)
isoerrorplus.SetName("IsoErrorPlus")
isoerrorminus=narf.hist_to_root(miniisoerrorminus)
isoerrorminus.SetName("IsoErrorMinus")
output_file = ROOT.TFile("makeefficiencies{name}.root","RECREATE")
output_file.cd()
idipplus.Write()
idipminus.Write()
triggerMCplus.Write()
triggerMCminus.Write()
triggerplus.Write()
triggerminus.Write()
triggererrorplus.Write()
triggererrorminus.Write()
isoMCplus.Write()
isoMCminus.Write()
isoplus.Write()
isominus.Write()
isoerrorplus.Write()
isoerrorminus.Write()
