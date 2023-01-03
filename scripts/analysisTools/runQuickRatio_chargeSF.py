#!/usr/bin/env python3

import os, array, math
import argparse

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

#sys.path.append(os.getcwd() + "/plotUtils/")
#from utility import *
from scripts.analysisTools.plotUtils.utility import *

#mainPath = "/eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_Sept2022_binnedInPtEta_mass60to120"
#workingPoints = ["reco", "trigger", "iso"]
mainPath = "/eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_globalMuons_testByCharge"
#inputhPath = mainPath
inputhPath = "/home/m/mciprian/tnp/egm_tnp_analysis/localplots/results_globalMuons_testByCharge"
#mainPath = "/eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_globalMuons_ntuplesXYZ_1orMoreNvalidHitsStandalone"
workingPoints = ["tracking"]

tag = "plusOverMinus" # just something to name the ratio
elements = ["plus", "minus"]
ratioTitle = f"{elements[0]}/{elements[1]}"
asymTitle = f"({elements[0]}-{elements[1]})/({elements[0]}+{elements[1]})"
#hToPlot = ["SF2D_nominal", "EffData2D", "EffMC2D", "EffDataAltSig2D", "EffMCAltSig2D"]
#hToPlot = ["SF2D_nominal", "EffData2D", "EffMC2D", "EffDataAltSig2D"]
hToPlot = ["SF2D_nominal", "EffData2D", "EffMC2D"]
#hToPlot = ["SF2D_nominal"]
eras = ["GtoH"] # ["BtoF", "GtoH"]:

xyRange = ""
subFolder = ""
skipAsym = True

#xyRange = " --xRange 0.0 2.4"
#subFolder = "eta_0p0To2p4"

for wp in workingPoints:
    for era in eras:
        outdir = f"{mainPath}/compareEffAndSFbyCharge/{era}/{wp}/{subFolder}/"
        file1  = f"{inputhPath}/efficiencies_{era}/mu_{wp}_{elements[0]}/allEfficiencies_2D.root"
        file2  = f"{inputhPath}/efficiencies_{era}/mu_{wp}_{elements[1]}/allEfficiencies_2D.root"
        
        for n in hToPlot:
            h1 = f"{n}_{era}"
            
            basecmd = f"python w-mass-13TeV/makeRatioTH2.py {file1} {n} {file2} {n} -o {outdir} {xyRange} -p -n ratio_{tag}_{n}_{era} -x 'Muon #eta' -y 'Muon p_{{T}} (GeV)' -z '{era} {n} ratio::0.99,1.01' --skip1DPlot -t '{ratioTitle}' --palette -1 "
            print()
            print(basecmd)
            
            basecmd_asym = f"python w-mass-13TeV/makeRatioTH2.py {file1} {n} {file2} {n} -o {outdir} {xyRange} -a -n asymmetry_{tag}_{n}_{era} -x 'Muon #eta' -y 'Muon p_{{T}} (GeV)' -z '{era} {n} asymmetry::-0.01,0.01' --skip1DPlot -t '{asymTitle}' --palette -1 "
            if not skipAsym:
                print()
                print(basecmd_asym)
        
    print()
    print()
