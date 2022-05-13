
import sys,array,math,os,copy,fnmatch
from collections import OrderedDict

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import plotter
import functions

import lz4.frame
import pickle
import narf

from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z


def parseProc(histCfg, procName, syst="", rebin=1):

    axis = histCfg['axis']
    hNames = histCfg['name'].split(",")

    
    bhist = None
    for hName in hNames:
        
        bhist = groups.readProc(hName, procName, axis=axis)
        if bhist == None: continue
        label = "%s_%s" % (hName, procName)
        break
        
    #print(bhist)
    rhist = narf.hist_to_root(bhist)
    rhist.Rebin(rebin)
    rhist.SetName(label)

    #print("Get histogram %s, yield=%d" % (label, rhist.Integral()))
    return rhist
                        
 
	
	
if __name__ == "__main__":

    flavor = "mumu"
    #flavor = "ee"

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_differential_combineOutput.root" % fitcfg)
      
      
    sys.exit()  
    groups = datagroupsLowPU_Z("mz_lowPU_%s.pkl.lz4" % flavor)

    if flavor == "mumu":
    
        procs = ['EWK', 'TTbar', 'Zmumu']
        procs = ['DYee', 'DYmumu', 'DYtautau', 'WJets', 'WZTo3LNu', 'WWTo2L2Nu', 'ZZ', 'TTbar']
        data = 'SingleMuon'
        
    if flavor == "ee":
    
        procs = ['DYee', 'DYmumu', 'DYtautau', 'WJets', 'WZTo3LNu', 'WWTo2L2Nu', 'ZZ', 'TTbar']
        data = 'SingleElectron'    
        
    histCfg = {"name": "dilepton", "axis": "mll" }
    h_data = parseProc(histCfg, data)
    print("Data \t\t %.2f" % h_data.Integral())
    
    hTot = 0
    for i,proc in enumerate(procs):
    
        hist = parseProc(histCfg, proc,)
        hTot += hist.Integral()
        print("%s \t\t %.2f" % (proc, hist.Integral()))
        
    print("MC tot \t\t %.2f" % hTot)
    print("Data/MC ratio: \t %.3f" % (h_data.Integral()/hTot))