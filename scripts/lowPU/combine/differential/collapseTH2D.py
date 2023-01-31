
# script collapse a TH2D to TH1 with one bin
# necessary for masked channel combine

import sys
import ROOT
import copy

baseProc = "reco_mll_Zmumu_genBin"
norms = []
hists = []
fIn = ROOT.TFile(sys.argv[1])

# compute normalizations for each gen process
for i in range(0, 10):

    hName = "%s%d_mumu_xsec" % (baseProc, i)
    hIn = fIn.Get(hName)
    if hIn: norms.append(hIn.Integral())
    else: norms.append(0)
    


for key in fIn.GetListOfKeys():

    hName = key.GetName()
    if not "reco_mll" in hName: continue
    if "Zmumu_genBin" in hName:
        genBin = int(hName.split("Zmumu_genBin")[1].split("_")[0])
        norm = norms[genBin]
        h_ = fIn.Get(hName)
        h_.SetName(hName+"_tmp")
        h_.Scale(1./norm)
        h = ROOT.TH1D(hName, "", 1, 0, 1)
        h.SetBinContent(1, h_.Integral())
        hists.append(copy.deepcopy(h))
        #print(h, h.Integral())
        
    else:
        print(hName)
        h_ = fIn.Get(hName)
        h_.SetName(hName+"_tmp")
        h_.Scale(1./h_.Integral())
        h = ROOT.TH1D(hName, "", 1, 0, 1)
        h.SetBinContent(1, h_.Integral())
        hists.append(copy.deepcopy(h))
        #print(h, h.Integral())

fIn.Close()

fIn = ROOT.TFile(sys.argv[1], "RECREATE")
for h in hists: h.Write() 
fIn.Close()