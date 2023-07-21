#!/usr/bin/env python3

from wremnants.datasets.datagroups2016 import make_datagroups_2016
from wremnants import histselections as sel
#from wremnants import plot_tools,theory_tools,syst_tools
from utilities import boostHistHelpers as hh,common, logging

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist

import numpy as np
from utilities import input_tools

import lz4.frame

import argparse
import os
import shutil
import logging
import re

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

from scripts.analysisTools.plotUtils.utility import *

sys.path.append(os.getcwd())

# should really avoid having numbers by hand, but this started with less then 10 in all graphs together :) 

def runWithBkg():
    
    outdir = "/eos/user/m/mciprian/www/WMassAnalysis/fromMyWremnants/fitResults/theoryAgnostic/v2/WMass_pt_eta_statOnly/impactScalingWithRecoPtBinning/"
    createPlotDirAndCopyPhp(outdir)

    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D", "", 900, 800)

    xvals = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0]
    yvals_gdata = [ 6.6, 12.2, 16.4, 21.9, 26.2, 30.8]
    yvals_gmc   = [ 9.3, 12.8, 17.0, 22.7, 27.2, 31.9]
    yvals_gtot  = [11.4, 17.7, 24.0, 31.6, 37.7, 44.3]

    gdata = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gdata))
    gmc   = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gmc))
    gtot  = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gtot))
    glist = [gtot, gdata, gmc]
    leg = ["Total", "Data stat", "MC stat"]
    
    drawGraphCMS(glist, "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50", "impacts_mw_ptRecoBinWidth",
                 outdir, leg,
                 legendCoords="0.2,0.7,0.6,0.88;1", vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    xvals = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0]
    yvals_halfMC_gdata = [ 2.11, 9.14, 12.5, 17.4, 21.4, 25.9]
    yvals_halfMC_gmc   = [ 9.17, 13.1, 17.8, 24.7, 30.4, 36.6]
    yvals_halfMC_gtot  = [ 9.41, 15.9, 21.7, 30.3, 37.1, 44.8]
    
    gdata_halfMC = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gdata))
    gmc_halfMC   = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gmc))
    gtot_halfMC  = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gtot))
    glist_halfMC = [gtot_halfMC, gdata_halfMC, gmc_halfMC]
    tag = "(1/2 MC stat)"
    leg_halfMC = [f"{l} {tag}" for l in leg]
    
    drawGraphCMS(glist_halfMC, "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50", "impacts_mw_ptRecoBinWidth_halfMCstat",
                 outdir, leg_halfMC,
                 legendCoords="0.2,0.7,0.6,0.88;1", vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    # together
    drawGraphCMS([*glist, *glist_halfMC], "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50",
                 "impacts_mw_ptRecoBinWidth_overlaid",
                 outdir, [*leg, *leg_halfMC],
                 legendCoords="0.55,0.14,0.95,0.44;1",
                 vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 vecLineStyle=[1,1,1,2,2,2],
                 vecMarkerStyle=[20,20,20,25,25,25],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    xvals_traditional =       [  0.1, 0.25, 0.5,  1.0,  1.5,  2.0]
    yvals_traditional_gdata = [ 2.27, 2.37, 2.39, 2.40, 2.41, 2.42]
    yvals_traditional_gmc   = [ 2.50, 2.49, 2.51, 2.52, 2.52, 2.53]
    yvals_traditional_gtot  = [ 7.69, 9.01, 9.88, 10.9, 12.2, 9.18]
    
    gdata_traditional = ROOT.TGraph(len(xvals_traditional), array('d', xvals_traditional), array('d', yvals_traditional_gdata))
    gmc_traditional   = ROOT.TGraph(len(xvals_traditional), array('d', xvals_traditional), array('d', yvals_traditional_gmc))
    gtot_traditional  = ROOT.TGraph(len(xvals_traditional), array('d', xvals_traditional), array('d', yvals_traditional_gtot))
    glist_traditional = [gtot_traditional, gdata_traditional, gmc_traditional]
    tag = ""
    leg_traditional = [f"{l} {tag}" for l in leg]
    
    drawGraphCMS(glist_traditional, "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,15", "impacts_mw_ptRecoBinWidth_traditional",
                 outdir, leg_traditional,
                 legendCoords="0.2,0.7,0.6,0.88;1", vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

def runSignalOnly():
    outdir = "/eos/user/m/mciprian/www/WMassAnalysis/fromMyWremnants/fitResults/theoryAgnostic/v2/signalOnly/WMass_pt_eta_statOnly/impactScalingWithRecoPtBinning/"
    createPlotDirAndCopyPhp(outdir)

    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D", "", 900, 800)

    xvals =       [ 0.1,  0.2,   0.4,   0.6,   1.0,   1.2,   1.5,   2.0]
    yvals_gdata = [ 7.86, 10.70, 14.16, 16.59, 20.51, 22.21, 24.67, 29.40]
    yvals_gmc   = [ 7.56, 10.31, 13.66, 16.05, 19.88, 21.54, 23.95, 28.55]
    yvals_gtot  = [10.9,  14.9,  19.7,  23.1,  28.6 , 30.9,  34.4,  41.0]

    gdata = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gdata))
    gmc   = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gmc))
    gtot  = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_gtot))
    glist = [gtot, gdata, gmc]
    leg = ["Total", "Data stat", "MC stat"]
    
    drawGraphCMS(glist, "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50", "impacts_mw_ptRecoBinWidth_signalOnly",
                 outdir, leg,
                 legendCoords="0.2,0.7,0.6,0.88;1", vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    xvals =              [ 0.1,  0.2,   0.4,   0.6,   1.0,   1.2,   1.5,   2.0]
    yvals_halfMC_gdata = [ 5.72,  7.94, 10.85, 12.94, 16.50, 18.22, 20.30, 24.78]
    yvals_halfMC_gmc   = [ 7.76, 10.79, 14.75, 17.63, 22.52, 24.89, 27.73, 33.90]
    yvals_halfMC_gtot  = [ 9.64, 13.4,  18.3,  21.9,  27.9,  30.8,  34.3,  42.0]
    
    gdata_halfMC = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gdata))
    gmc_halfMC   = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gmc))
    gtot_halfMC  = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_halfMC_gtot))
    glist_halfMC = [gtot_halfMC, gdata_halfMC, gmc_halfMC]
    tag = "(1/2 MC stat)"
    leg_halfMC = [f"{l} {tag}" for l in leg]
    
    drawGraphCMS(glist_halfMC, "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50", "impacts_mw_ptRecoBinWidth_halfMCstat",
                 outdir, leg_halfMC,
                 legendCoords="0.2,0.7,0.6,0.88;1",
                 vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 vecLineStyle=[2,2,2],
                 vecMarkerStyle=[25,25,25],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    # together
    drawGraphCMS([*glist, *glist_halfMC], "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50",
                 "impacts_mw_ptRecoBinWidth_overlaid",
                 outdir, [*leg, *leg_halfMC],
                 legendCoords="0.55,0.14,0.95,0.44;1",
                 vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 vecLineStyle=[1,1,1,2,2,2],
                 vecMarkerStyle=[20,20,20,25,25,25],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)

    # compute true sigma from Lorenzo's calculations, using gdata
    yvals_true_gdata  = [0.0 for i in range(len(yvals_gdata))]
    for iy,y in enumerate(yvals_gdata):
        r = yvals_halfMC_gdata[iy]/y
        if r < 1./np.sqrt(2):
            print("Error for x = {xvals[iy]}: r = {r} is too small (must be >= 1./sqrt(2) )")
            fact = 0
        else:
            fact = r / np.sqrt(2 * r * r - 1.0)
        yvals_true_gdata[iy] = fact * y
        
    gdata_true = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals_true_gdata))

    drawGraphCMS([gdata_true, gdata, gdata_halfMC], "p_{T} bin width (GeV)", "Impact on m_{W} (MeV)::0,50",
                 "impacts_mw_ptRecoBinWidth_trueSensitivityData",
                 outdir, ["Data stat true", "Data stat (full MC stat)", f"Data stat {tag}"],
                 legendCoords="0.2,0.7,0.6,0.88;1",
                 vecMCcolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue],
                 vecLineStyle=[1,1,1],
                 vecMarkerStyle=[20,20,20],
                 passCanvas=canvas1D, graphDrawStyle="pl", legEntryStyle="PL",
                 skipLumi=True)
    

    rf = safeOpenFile(outdir+"impactMW_signalOnly_statOnly.root", mode="RECREATE")
    rf.cd()
    gdata.Write("gdata_fullMCstat")  
    gmc.Write("gmc_fullMCstat")    
    gtot.Write("gtot_fullMCstat")    
    gdata_halfMC.Write("gdata_halfMCstat")    
    gmc_halfMC.Write("gmc_halfMCstat")    
    gtot_halfMC.Write("gtot_halfMCstat")    
    rf.Close()
    
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    runSignalOnly()
