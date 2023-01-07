#!/usr/bin/env python3

import os, re, array, math
import time
import argparse

## safe batch mode                                 
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

#sys.path.append(os.getcwd() + "/plotUtils/")
#from utility import *
from scripts.analysisTools.plotUtils.utility import *

sys.path.append(os.getcwd())

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("outputfolder", type=str, nargs=1)
    parser.add_argument("inputfileMC",   type=str, nargs=2, help="Input files for MC")
    parser.add_argument("labels",   type=str, nargs=2, help="Labels for efficiency plots")
    parser.add_argument(     "--rebin-y", dest="rebinY", default=1, type=int, help="To rebin y axis (pt)")
    parser.add_argument(     "--rebin-z", dest="rebinZ", default=1, type=int, help="To rebin z axis (eta)")
    parser.add_argument(     '--nContours', dest='nContours', default=51, type=int,
                             help='Number of contours in palette. Default is 51')
    parser.add_argument(     '--palette'  , dest='palette',   default=87, type=int,
                             help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette', action='store_true',
                             help='Inverte color ordering in palette')
    args = parser.parse_args()

    for il,l in enumerate(args.labels):
        args.labels[il] = l.replace(" ", "_") # safety thing since it is also used to name output files

    ROOT.TH1.SetDefaultSumw2()

    outdir = args.outputfolder[0]
    createPlotDirAndCopyPhp(outdir)
    
    f = safeOpenFile(args.inputfileMC[0])
    hmcfail3D = safeGetObject(f, "fail_mu_DY_postVFP")
    hmcpass3D = safeGetObject(f, "pass_mu_DY_postVFP")
    f.Close()

    f = safeOpenFile(args.inputfileMC[1])
    hmcfail3Dtest = safeGetObject(f, "fail_mu_DY_postVFP")
    hmcpass3Dtest = safeGetObject(f, "pass_mu_DY_postVFP")
    f.Close()

    hists = [hmcpass3D, hmcfail3D, hmcpass3Dtest, hmcfail3Dtest]
    for h in hists:
        if args.rebinY > 1: h.RebinY(args.rebinY)
        if args.rebinZ > 1: h.RebinZ(args.rebinZ)
    
    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 800, 800)

    print(f"{hmcpass3D.GetNbinsZ()} eta bins")
    print(f"{hmcpass3D.GetNbinsY()} pt  bins")
    
    hmcpass2D = hmcpass3D.Project3D("yze") # do y versus z which is pt versus eta 
    hmcpass2D.SetName("hmcpass2D")
    hmcfail2D = hmcfail3D.Project3D("yze") # do y versus z which is pt versus eta 
    hmcfail2D.SetName("hmcfail2D")
    hmcpass2Dtest = hmcpass3Dtest.Project3D("yze") # do y versus z which is pt versus eta 
    hmcpass2Dtest.SetName("hmcpass2Dtest")
    hmcfail2Dtest = hmcfail3Dtest.Project3D("yze") # do y versus z which is pt versus eta 
    hmcfail2Dtest.SetName("hmcfail2Dtest")

    effMC = copy.deepcopy(hmcpass2D.Clone(args.labels[0]))
    hmcfail2D.Add(hmcpass2D)
    effMC.Divide(hmcfail2D)
    effMC.SetTitle(f"{args.labels[0]}")

    effMCtest = copy.deepcopy(hmcpass2Dtest.Clone(args.labels[1]))
    hmcfail2Dtest.Add(hmcpass2Dtest)
    effMCtest.Divide(hmcfail2Dtest)
    effMCtest.SetTitle(f"{args.labels[1]}")

    drawCorrelationPlot(effMC, "muon #eta", "muon p_{T} (GeV)", f"MC efficiency",
                        f"effMC_{args.labels[0]}", plotLabel="ForceTitle", outdir=outdir,
                        smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                        draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                        nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)

    drawCorrelationPlot(effMCtest, "muon #eta", "muon p_{T} (GeV)", f"MC efficiency",
                        f"effMC_{args.labels[1]}", plotLabel="ForceTitle", outdir=outdir,
                        smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                        draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                        nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)

    effMC.Divide(effMCtest)
    effMC.SetTitle(f"{args.labels[0]} / {args.labels[1]}")
    drawCorrelationPlot(effMC, "muon #eta", "muon p_{T} (GeV)", f"MC efficiency ratio",
                        f"effRatioMC_{args.labels[0]}_over_{args.labels[1]}", plotLabel="ForceTitle", outdir=outdir,
                        smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                        draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                        nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
