#!/usr/bin/env python3
import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

sys.path.append(os.getcwd())

logging.basicConfig(level=logging.INFO)

ROOT.gInterpreter.ProcessLine(".O3")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with TH2 histograms")
    parser.add_argument("outdir",   type=str, nargs=1, help="Output folder")
    parser.add_argument("-l", "--lumi",     type=str, default=None, help="Luminosity to print on canvas, by default it is not printed")
    parser.add_argument("--pt-range-projection", dest="ptRangeProjection", default=(0,-1), type=float, nargs=2, help="Pt range to select bins to use for 1D projection (for upper range remember that upper bin edge belongs to next bin in ROOT)")
    parser.add_argument("--wlike", dest="isWlike", action="store_true", help="Flag for W-like analysis")
    args = parser.parse_args()
           
    fname = args.rootfile[0]
    outdir = args.outdir[0]
    createPlotDirAndCopyPhp(outdir)
    
    ROOT.TH1.SetDefaultSumw2()

    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    cwide = ROOT.TCanvas("cwide","",2400,600)                      
    adjustSettings_CMS_lumi()
    
    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    processes = ["Data", "Wmunu", "Wtau", "Zmumu", "Ztautau", "Fake", "Top", "Diboson"]
    if args.isWlike:
        processes = ["Data", "Zmumu", "Ztautau", "Other"]
    charges = ["plus", "minus"]

    xAxisName = "Muon #eta"
    yAxisName = "Muon p_{T} (GeV)"

    colors = {"Wmunu"      : ROOT.kRed+2,
              "Zmumu"      : ROOT.kAzure+2,
              "Wtau"       : ROOT.kCyan+1,
              "Ztautau"    : ROOT.kSpring+9,
              "Top"        : ROOT.kGreen+2,
              "Diboson"    : ROOT.kViolet,
              "Fake"       : ROOT.kGray,
              "Other"      : ROOT.kGray}

    legEntries = {"Wmunu"      : "W#rightarrow#mu#nu",
                  "Zmumu"      : "Z#rightarrow#mu#mu",
                  "Wtau"       : "W#rightarrow#tau#nu",
                  "Ztautau"    : "Z#rightarrow#tau#tau",
                  "Top"        : "t quark",
                  "Diboson"    : "Diboson",
                  "Fake"       : "Multijet",
                  "Other"      : "Other"}

    for charge in charges:
    
        # read histograms
        nomihists = {}
        infile = safeOpenFile(fname)
        for proc in processes:
            nomihists[proc] = safeGetObject(infile, f"x_{proc}_{charge}", detach=True)
        infile.Close()

        hdata2D = nomihists["Data"]
        hmc2D = [nomihists[x] for x in nomihists.keys() if x != "Data"]

        outdir_dataMC = f"{outdir}dataMC_{charge}/"
        createPlotDirAndCopyPhp(outdir_dataMC)
        
        fShapesName = outdir_dataMC + "plots.root"
        fShapes = safeOpenFile(fShapesName, mode="RECREATE")
        
        hmc2D = sorted(hmc2D, key= lambda x: x.Integral()) # , reverse=True) 
        
        stack_eta = ROOT.THStack("stack_eta", "signal and backgrounds")
        stack_pt = ROOT.THStack("stack_pt", "signal and backgrounds")
        stack_unrolled = ROOT.THStack("stack_unrolled", "signal and backgrounds")

        ratio2D = hdata2D.Clone("dataOverMC2D")
        den2D = hdata2D.Clone("sigAndBkg2D")
        den2D.Reset("ICESM")

        hdata2D.SetMarkerColor(ROOT.kBlack)
        hdata2D.SetLineColor(ROOT.kBlack)
        #hdata2D.SetLineWidth(2)
        hdata2D.SetMarkerStyle(20)
        hdata2D.SetMarkerSize(1)
        hdata2D.SetTitle("")
        hdata_unrolled = unroll2Dto1D(hdata2D, newname=f"{hdata2D.GetName()}_unrolled")

        # for projections along eta
        ptRange = ""
        if args.ptRangeProjection[0] < args.ptRangeProjection[1]:
            lowPtbin = max(1, hdata2D.GetYaxis().FindFixBin(args.ptRangeProjection[0]))
            highPtbin = min(hdata2D.GetNbinsY(), hdata2D.GetYaxis().FindFixBin(args.ptRangeProjection[1])) # hdata2D.GetNbinsY()
            ptRange = "_%gTo%g" % (hdata2D.GetYaxis().GetBinLowEdge(lowPtbin), hdata2D.GetYaxis().GetBinLowEdge(1+highPtbin))
            ptRange = ptRange.replace(".","p")
        else:
            lowPtbin = 0
            highPtbin = -1

        hdata_eta = hdata2D.ProjectionX("data_eta",lowPtbin,highPtbin,"e")
        hdata_pt  = hdata2D.ProjectionY("data_pt",0,-1,"e")

        legend = ROOT.TLegend(0.2,0.72,0.95,0.92)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetNColumns(3)

        leg_unrolled = prepareLegend(legWidth=0.60, textSize=0.045, nColumns=5)
        
        hdata2D.Write("data2D")
        hdata_eta.Write()
        hdata_pt.Write()
        hdata_unrolled.Write()

        legend.AddEntry(hdata2D, "Data", "EP")
        leg_unrolled.AddEntry(hdata_unrolled, "Data", "EP")
        hmc_unrolled = []
        for h in hmc2D:
            h.SetTitle("")
            h.SetFillColor(colors[h.GetName().split("_")[1]])
            h.SetLineColor(ROOT.kBlack)
            stack_eta.Add(h.ProjectionX(f"{h.GetName()}_eta",lowPtbin,highPtbin,"e"))
            stack_pt.Add( h.ProjectionY(f"{h.GetName()}_pt",0,-1,"e"))
            den2D.Add(h)
            h.Write()
            hmc_unrolled.append(unroll2Dto1D(h, newname=f"{h.GetName()}_unrolled"))
            hmc_unrolled[-1].SetFillColor(colors[h.GetName().split("_")[1]])
            hmc_unrolled[-1].SetLineColor(ROOT.kBlack)
            hmc_unrolled[-1].SetMarkerSize(0)
            hmc_unrolled[-1].SetMarkerStyle(0)
            hmc_unrolled[-1].Write()
            stack_unrolled.Add(hmc_unrolled[-1])
        for i in sorted(range(len(hmc2D)), reverse=True):
            legend.AddEntry(hmc2D[i], legEntries[hmc2D[i].GetName().split("_")[1]], "F")
            leg_unrolled.AddEntry(hmc_unrolled[i], legEntries[hmc2D[i].GetName().split("_")[1]], "F")

        stack_eta.Write()
        stack_pt.Write()
        den2D.Write()

        ratio2D.Divide(den2D)
        ratio2D.Write()

        drawTH1dataMCstack(hdata_eta, stack_eta, "Muon #eta", "Events", "muon_eta" + ptRange,
                           outdir_dataMC, legend, ratioPadYaxisNameTmp="Data/pred::0.92,1.08", passCanvas=canvas1D, lumi=args.lumi,
                           drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
        )
        drawTH1dataMCstack(hdata_pt, stack_pt, "Muon p_{T} (GeV)", "Events", "muon_pt",
                           outdir_dataMC, legend, ratioPadYaxisNameTmp="Data/pred::0.92,1.08", passCanvas=canvas1D, lumi=args.lumi,
                           drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
        )

        
        ratio2D.SetTitle("data / (signal + background)")
        drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "Data/pred",
                            f"muon_eta_pt_dataMCratio", plotLabel="ForceTitle", outdir=outdir_dataMC,
                            palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)
        drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "Data/pred. statistical uncertainty",
                            f"muon_eta_pt_dataMCratio_absUncertainty", plotLabel="ForceTitle", outdir=outdir_dataMC,
                            palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True, plotError=True)


        allHists = hmc2D + [hdata2D]
        for h in allHists:
            if "Data" not in h.GetName():
                h.SetTitle(legEntries[h.GetName().split("_")[1]] + " " + charge)
            else:
                h.SetTitle(f"Data {charge}")
            drawCorrelationPlot(h, xAxisName, yAxisName, "Events",
                                f"muon_eta_pt_{h.GetName()}", plotLabel="ForceTitle", outdir=outdir_dataMC,
                                palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)

        #
        etabins = [round(hdata2D.GetXaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hdata2D.GetNbinsX())]
        ptbins = [round(hdata2D.GetYaxis().GetBinLowEdge(i), 1) for i in range(1,  2 + hdata2D.GetNbinsY())]
        recoBins = templateBinning(etabins, ptbins)
        nRecoBins = recoBins.NTotBins
        #following array is used to call function dressed2DfromFit()
        binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]        
        cnameUnroll = f"muon_etaPtUnrolled"
        XlabelUnroll = "unrolled template along #eta:  #eta #in [%.1f, %.1f]" % (recoBins.etaBins[0], recoBins.etaBins[-1])
        YlabelUnroll = "Events::%.2f,%.2f" % (0, 2.*hdata_unrolled.GetBinContent(hdata_unrolled.GetMaximumBin()))
        # to draw panels in the unrolled plots
        ptBinRanges = []
        for ipt in range(0,recoBins.Npt):
            #ptBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=recoBins.ptBins[ipt], ptmax=recoBins.ptBins[ipt+1]))
            ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(recoBins.ptBins[ipt]), ptmax=int(recoBins.ptBins[ipt+1])))
            
        drawTH1dataMCstack(hdata_unrolled, stack_unrolled, XlabelUnroll, YlabelUnroll, cnameUnroll, outdir_dataMC, leg_unrolled, ratioPadYaxisNameTmp="Data/pred::0.92,1.08",
                           passCanvas=cwide,
                           wideCanvas=True, leftMargin=0.05,rightMargin=0.02,lumi=args.lumi, 
                           drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta),
                           textForLines=ptBinRanges, etaptbinning=binning, noLegendRatio=True, textSize=0.04, textAngle=30)
