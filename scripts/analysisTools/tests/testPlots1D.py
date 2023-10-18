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
from utilities.io_tools import input_tools

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
from scripts.analysisTools.tests.cropNegativeTemplateBins import cropNegativeContent

def plotDistribution1D(hdata, hmc, datasets, outfolder_dataMC, canvas1Dshapes=None,
                       xAxisName="variable", plotName="variable_failIso_jetInclusive",
                       draw_both0_noLog1_onlyLog2=1, ratioPadYaxisTitle="Data/pred::0.9,1.1"):
    
    createPlotDirAndCopyPhp(outfolder_dataMC)
    if not canvas1Dshapes:
        canvas1Dshapes = ROOT.TCanvas("canvas1Dshapes","",700,800)
                
    legend = ROOT.TLegend(0.2,0.72,0.95,0.92)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(3)

    for d in datasets:
        if d == "Data":
            legend.AddEntry(hdata, "Data", "EP")
        else:
            cropNegativeContent(hmc[d])
            hmc[d].SetFillColor(colors_plots_[d])
            hmc[d].SetLineColor(ROOT.kBlack)
            hmc[d].SetMarkerSize(0)
            hmc[d].SetMarkerStyle(0)

    stack_1D = ROOT.THStack("stack_1D", "signal and backgrounds")
    hmcSortedKeys = sorted(hmc.keys(), key= lambda x: hmc[x].Integral())
    for i in hmcSortedKeys:
        stack_1D.Add(hmc[i])
    # reverse sorting for legend, first the ones with larger integral
    for i in list(reversed(hmcSortedKeys)):
        legend.AddEntry(hmc[i], legEntries_plots_[i], "LF")

    drawTH1dataMCstack(hdata, stack_1D, xAxisName, "Events", plotName,
                       outfolder_dataMC, legend, ratioPadYaxisNameTmp=ratioPadYaxisTitle, passCanvas=canvas1Dshapes,
                       lumi="16.8", drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True,
                       draw_both0_noLog1_onlyLog2=draw_both0_noLog1_onlyLog2)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1)
    parser.add_argument("outputfolder",   type=str, nargs=1)
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");
    parser.add_argument('-p','--processes', default=None, nargs='*', type=str,
                        help='Choose what processes to plot, otherwise all are done')
    parser.add_argument('--plot', nargs='+', type=str,
                        help='Choose what distribution to plot by name')
    parser.add_argument("-x", "--x-axis-name", dest="xAxisName", nargs='+', type=str, help="x axis name")
    args = parser.parse_args()
    
    logger = logging.setup_logger(os.path.basename(__file__), args.verbose)
    
    fname = args.inputfile[0]
    outdir = args.outputfolder[0]
    createPlotDirAndCopyPhp(outdir)
        
    ROOT.TH1.SetDefaultSumw2()

    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    groups = make_datagroups_2016(fname)
    datasets = groups.getNames()
    if args.processes is not None and len(args.processes):
        datasets = list(filter(lambda x: x in args.processes, datasets))
    logger.info(f"Will plot datasets {datasets}")

    for ip,p in enumerate(args.plot):

        groups.setNominalName(p)
        groups.loadHistsForDatagroups(p, syst="", procsToRead=datasets)

        histInfo = groups.getDatagroups()
        rootHists = {}
        
        for d in datasets:
            hnarf = histInfo[d].hists[p]
            rootHists[d] = narf.hist_to_root(hnarf)
            rootHists[d].SetName(f"{p}_{d}")

        hdata = rootHists["Data"]
        hmc = {d : rootHists[d] for d in datasets if d != "Data"}
        plotDistribution1D(hdata, hmc, datasets, outdir, canvas1Dshapes=canvas1D,
                           xAxisName=args.xAxisName[ip], plotName=p)

