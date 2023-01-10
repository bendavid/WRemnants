#!/usr/bin/env python3

# mainly fro tracking efficiencies and failing probes
# example:
# python w-mass-13TeV/compareTnpMass.py /home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_trackerMuons/tnp_tracking_data_vertexWeights1_oscharge0.root /home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_trackerMuons/tnp_tracking_mc_vertexWeights1_oscharge0.root plots/TnP/Steve_Marc_Raj/testTrackerMuons/tracking/ --zbin 1 3


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
    parser.add_argument("inputfileData", type=str, nargs=1, help="Input file for data")
    parser.add_argument("inputfileMC",   type=str, nargs=1, help="Input file for MC")
    parser.add_argument("outputfolder", type=str, nargs=1)
    #parser.add_argument("--hname", default=None, required=True, type=str, nargs=2, help="Name of histograms to get from the input files")
    parser.add_argument("-x", "--x-axis-name", dest="xAxisName", default="Invariant mass (GeV) ", help="x axis name")
    parser.add_argument(     "--rebin-x", dest="rebinX", default=1, type=int, help="To rebin x axis (mass)")
    parser.add_argument(     "--rebin-y", dest="rebinY", default=1, type=int, help="To rebin y axis (pt)")
    parser.add_argument(     "--rebin-z", dest="rebinZ", default=1, type=int, help="To rebin z axis (eta)")
    parser.add_argument(     "--ybin", type=int, nargs=2, default=[0, 0], help="Bins for y axis to plot, default is to do all")
    parser.add_argument(     "--zbin", type=int, nargs=2, default=[0, 0], help="Bins for z axis to plot, default is to do all")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    outdir = args.outputfolder[0]
    createPlotDirAndCopyPhp(outdir)
    
    f = safeOpenFile(args.inputfileData[0])
    hdata3D = safeGetObject(f, "fail_mu_RunGtoH")
    f.Close()
    
    f = safeOpenFile(args.inputfileMC[0])
    hmc3D = safeGetObject(f, "fail_mu_DY_postVFP")
    hmcalt3D = safeGetObject(f, "pass_mu_DY_postVFP_alt")
    f.Close()
    
    hmcTot3D = copy.deepcopy(hmc3D.Clone("all_mu_DY_postVFP"))
    hmcTot3D.Add(hmcalt3D)

    hists = [hmcTot3D, hmc3D, hdata3D]
    for h in hists:
        if args.rebinX > 1: h.RebinX(args.rebinX)
        if args.rebinY > 1: h.RebinY(args.rebinY)
        if args.rebinZ > 1: h.RebinZ(args.rebinZ)
    
    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 900, 800)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.04)
    canvas.cd()

    print(f"{hmcTot3D.GetNbinsZ()} eta bins")
    print(f"{hmcTot3D.GetNbinsY()} pt  bins")

    iymin = 1
    iymax = hmcTot3D.GetNbinsY()
    if args.ybin[0] > 0 and args.ybin[1] > 0:
        iymin,iymax = args.ybin

    izmin = 1
    izmax = hmcTot3D.GetNbinsZ()
    if args.zbin[0] > 0 and args.zbin[1] > 0:
        izmin,izmax = args.zbin
    
    for ieta in range(1, 1 + hmcTot3D.GetNbinsZ()):
        if not (izmin <= ieta <= izmax):
            continue
        for ipt in range(1, 1 + hmcTot3D.GetNbinsY()):
            if not (iymin <= ipt <= iymax):
                continue
            hmcTot = hmcTot3D.ProjectionX("hmcTot", ipt, ipt, ieta, ieta, "e")
            hmc = hmc3D.ProjectionX("hmc", ipt, ipt, ieta, ieta, "e")
            hdata = hdata3D.ProjectionX("hdata", ipt, ipt, ieta, ieta, "e")

            hmcTot.SetMarkerSize(0)
            #hmcTot.SetFillColor(ROOT.kGreen+2)
            hmcTot.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
            hmcTot.SetFillStyle(1001)

            hmc.SetLineColor(ROOT.kBlue)
            hmc.SetLineWidth(2)

            hdata.SetMarkerStyle(20)
            hdata.SetMarkerSize(1)

            hmcTotScale = hdata.Integral()/hmcTot.Integral()
            hmcTot.Scale(hmcTotScale)
            hmcScale = hdata.Integral()/hmc.Integral()
            hmc.Scale(hmcScale)

            miny, maxy =  getMinMaxMultiHisto([hdata, hmc, hmcTot], excludeEmpty=False, sumError=False)
            
            hmcTot.SetStats(0)
            hmcTot.SetMarkerSize(0)
            hmcTot.SetMarkerStyle(0)
            hmcTot.GetXaxis().SetTitle(args.xAxisName)
            hmcTot.GetXaxis().SetTitleOffset(1.3)
            hmcTot.GetXaxis().SetTitleSize(0.05)
            hmcTot.GetXaxis().SetLabelSize(0.04)
            hmcTot.GetYaxis().SetTitle("Data events / bin")
            hmcTot.GetYaxis().SetTitleOffset(1.25)
            hmcTot.GetYaxis().SetTitleSize(0.05)
            hmcTot.GetYaxis().SetLabelSize(0.04)
            hmcTot.GetYaxis().SetRangeUser(0, 1.25 * maxy)

            hmcTot.Draw("HE")
            hmc.Draw("HIST SAME")
            hdata.Draw("EP SAME")

            header = "{} < #eta < {}".format(round(hmcTot3D.GetZaxis().GetBinLowEdge(ieta),1), round(hmcTot3D.GetZaxis().GetBinUpEdge(ieta),1))
            header += "   ---   "
            header += "{} < p_{{T}} < {} GeV".format(round(hmcTot3D.GetYaxis().GetBinLowEdge(ipt),0), round(hmcTot3D.GetYaxis().GetBinUpEdge(ipt),0))
            
            leg = ROOT.TLegend(0.2, 0.78, 0.9, 0.9)
            leg.SetNColumns(3)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetFillColorAlpha(0,0.6)
            leg.SetBorderSize(0)
            leg.SetHeader(header)
            leg.AddEntry(hdata,  "Data", "EP")
            leg.AddEntry(hmc,    "MC (fail probes)"  , "L")
            leg.AddEntry(hmcTot, "MC (all probes)", "LF")
            leg.Draw("same")

            canvasName = f"failProbeMass_ieta_{ieta}_ipt_{ipt}"

            canvas.RedrawAxis("sameaxis")
            for ext in ["png","pdf"]:
                canvas.SaveAs(f"{outdir}/{canvasName}.{ext}")
