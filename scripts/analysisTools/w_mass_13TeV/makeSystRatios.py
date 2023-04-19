#!/usr/bin/env python3

# plot ratios of histogram variation over nominal
# example
# python w-mass-13TeV/makeSystRatios.py /scratch/mciprian/CombineStudies/WMass/scetlibCorr_nnpdf31_testSF_trashdebug/byHelicityPtCharge_correlateEffStatIsoByCharge/WMassCombineInput.root /eos/user/m/mciprian/www/WMassAnalysis/fromMyWremnants/Wmass_fit/scetlibCorr_nnpdf31_testSF_trashdebug/byHelicityPtCharge_correlateEffStatIsoByCharge/makeSystRatios/ -s ".*effStatTnP_trigger_eta20.*q1.*" -c plus --systPostfix effStatTnP_trigger_eta20plus

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

#sys.path.append(os.getcwd() + "/plotUtils/")
#from utility import *
from scripts.analysisTools.plotUtils.utility import *

def plotUnrolledHistogram(h, process, syst, outdir, canvas, hist2DforBins, yAxisTitle="syst/nomi",
                          errorBars=False, channelCharge=None,
                          canvasNamePrefix="", canvasFullName=None):

    canvas.cd()
    yTitleOffset = 0.65
    histTitle = f"process: {process},      {syst}" 
    if channelCharge:
        histTitle = f"Reco charge: {channelCharge} -- " + histTitle 
    h.SetTitle(histTitle)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitle("Unrolled muon #eta-p_{T} bin")
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetYaxis().SetTitleOffset(yTitleOffset)
    h.GetYaxis().SetTitle(yAxisTitle)
    h.GetYaxis().SetTickSize(0.01)

    ########
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    #h.SetMarkerStyle(20)
    #h.SetMarkerSize(0.9)
    h.SetMarkerStyle(1)
    #h.SetFillColor(ROOT.kGray+3)
    #h.SetFillColorAlpha(ROOT.kGray+3, 0.4)
    
    miny, maxy = getMinMaxHisto(h, sumError=True if errorBars else False)
    diff = maxy - miny
    #print(f"min,max = {miny}, {maxy}:   diff = {diff}")
    #print(f"bin min,max = {h.GetMinimumBin()}, {h.GetMaximumBin()}")
    miny -= diff * 0.1
    maxy += diff * 0.3
    #print(f"new min,max = {miny}, {maxy}")
    h.GetYaxis().SetRangeUser(miny, maxy)

    h.Draw("HE" if errorBars else "HIST")

    bintext = ROOT.TLatex()
    textSize = 0.04
    textAngle = 30
    bintext.SetTextSize(textSize)
    bintext.SetTextFont(42)
    bintext.SetTextAngle(textAngle)

    # to draw panels in the unrolled plots
    nptBins = int(hist2DforBins.GetNbinsY())
    etarange = float(hist2DforBins.GetNbinsX())

    vertline = ROOT.TLine(36, 0, 36, canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 larger hatches
    for i in range(1, nptBins): # do not need line at canvas borders
        vertline.DrawLine(etarange*i, miny, etarange*i, maxy)

    ptBinRanges = []
    for ipt in range(0, nptBins):
        ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                           ptmax=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+2))))
    offsetText = etarange / 4.0
    for i in range(0,len(ptBinRanges)): # we need nptBins texts
        bintext.DrawLatex(etarange*i + offsetText, maxy - 0.2*diff, ptBinRanges[i])        
        
    line = ROOT.TF1("horiz_line", "1", h.GetXaxis().GetBinLowEdge(1), h.GetXaxis().GetBinLowEdge(h.GetNbinsX()+1))
    line.SetLineWidth(1)
    line.SetLineColor(ROOT.kRed)
    line.Draw("Lsame")
        
    canvas.RedrawAxis("sameaxis")

    if canvasNamePrefix:
        canvasNamePrefix += "_"

    canvasName = canvasFullName if canvasFullName else f"{canvasNamePrefix}unrolled_{process}_{syst}"
    for ext in ["png","pdf"]:
        canvas.SaveAs(f"{outdir}{canvasName}.{ext}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input root file")
    parser.add_argument("outdir",   type=str, nargs=1, help="Folder for plots")
    parser.add_argument("-s", "--systematics",    type=str, default=".*pdf.*", help="Comma separated list of syst names or regular expressions to select systematics to make ratios with nominal")
    parser.add_argument("-p", "--processes",    type=str, default="Wmunu", help="Comma separated list of processes to plot (full name please)")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use 0 for a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--plot', dest='plot', default="unrolled", choices=["2D", "unrolled", "all"],  help='Plots to make')
    parser.add_argument(     '--plotNominal', dest='plotNominal' , default=False , action='store_true',   help='Plot nominal histogram')
    parser.add_argument(     '--plotSystOriginal', dest='plotSystOriginal' , default=False , action='store_true',   help='Plot syst histogram also before ratio with nominal')
    parser.add_argument(     '--plotStat', dest='plotStat' , default=False , action='store_true',   help='Plot histogram to show statistical uncertainty on each nominal process')
    parser.add_argument(     '--addErrorBars', dest='addErrorBars' , default=False , action='store_true',   help='Plot error bars on ratio plots for systematics (not recommended if not for test, but better to use --plotStat to check the stat uncertainty)')
    parser.add_argument(     '--statUncRatio', dest='statUncRatio' , default=False , action='store_true',   help='Add ratio of stat uncertainties between syst and nomi, so see if there are correlations with fluctuating bins in the ratio of yields')
    parser.add_argument(     '--source', dest='source' , type=str, default="wrem", choices=["wrem", "cmg"],  help='Select which tool was used to make histograms (cmg=CMGTools or wrem=WRemnants), the naming convention is different for processes and systematics')
    parser.add_argument('-c','--charge', dest='charge', default=None, choices=["plus", "minus"], type=str, help='For source=wrem, needs to specify the charge, since all histograms are in the same file')
    parser.add_argument(     '--systPostfix', type=str, default=None, help="Postfix for plot showing some syst variations")
    args = parser.parse_args()

    doUnrolled = (args.plot == "unrolled") or (args.plot == "all")
    do2D       = (args.plot == "2D") or (args.plot == "all")
    
    fname = args.rootfile[0]
    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

    isWrem = True if args.source == "wrem" else False
    if isWrem and not args.charge:
        print("For histograms from WRemnants the charge must be specified using -c [minus|plus].")
        quit()
    
    processes = args.processes.split(',')
    inputsysts = args.systematics.split(',')
    regexp_syst = re.compile(args.systematics.replace(',','|'))

    nominals = {p : None for p in processes}
    #print(nominals)
    ratios = {p : [] for p in processes}

    canvas = ROOT.TCanvas("canvas","",900,800) 
    canvas1D = ROOT.TCanvas("canvas","",800,900) 

    canvas_unroll = ROOT.TCanvas("canvas_unroll","",3000,800) 
    leftMargin = 0.06
    rightMargin = 0.01
    bottomMargin = 0.12
    canvas_unroll.SetTickx(1)
    canvas_unroll.SetTicky(1)
    canvas_unroll.cd()
    canvas_unroll.SetLeftMargin(leftMargin)
    canvas_unroll.SetRightMargin(rightMargin)
    canvas_unroll.cd()
    canvas_unroll.SetBottomMargin(bottomMargin)

    setTDRStyle() # this one removes the stat box
    
    f = safeOpenFile(fname)
    # get nominals
    for p in processes:
        if isWrem:
            nominals[p] = safeGetObject(f, f"x_{p}_{args.charge}")
        else:
            nominals[p] = safeGetObject(f, f"x_{p}")
        nominals[p].SetTitle(p)
        if args.plotNominal:
            if do2D:
                drawCorrelationPlot(nominals[p], "Muon #eta", "Muon p_{T} (GeV)", f"Events",
                                    nominals[p].GetName(), plotLabel="ForceTitle", outdir=outdir,
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                                    palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette,
                                    passCanvas=canvas, drawOption="COLZ0")
            if doUnrolled:
                nomi_unrolled = unroll2Dto1D(nominals[p], newname=f"unrolled_{nominals[p].GetName()}", cropNegativeBins=False)
                plotUnrolledHistogram(nomi_unrolled, p, nominals[p].GetName().replace("x_",""), outdir, canvas_unroll, nominals[p], yAxisTitle="Events", channelCharge=args.charge)
            
        if args.plotStat:
            h_relativeStatUnc = copy.deepcopy(nominals[p].Clone(f"relativeStatUnc_{p}_{args.charge}" if isWrem else f"relativeStatUnc_{p}"))
            for ib in range(1+h_relativeStatUnc.GetNcells()):
                value = 0.0 if nominals[p].GetBinContent(ib) == 0.0 else nominals[p].GetBinError(ib) / nominals[p].GetBinContent(ib)
                h_relativeStatUnc.SetBinContent(ib, value)
            if do2D:
                drawCorrelationPlot(h_relativeStatUnc, "Muon #eta", "Muon p_{T} (GeV)", f"Relative stat. uncertainty",
                                    h_relativeStatUnc.GetName(), plotLabel="ForceTitle", outdir=outdir,
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                                    palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette,
                                    passCanvas=canvas, drawOption="COLZ0")
            if doUnrolled:
                ratio_unrolled = unroll2Dto1D(h_relativeStatUnc, newname=f"unrolled_{h_relativeStatUnc.GetName()}", cropNegativeBins=False)
                plotUnrolledHistogram(ratio_unrolled, p, f"relStatUnc_{args.charge}", outdir, canvas_unroll, h_relativeStatUnc, yAxisTitle="Relative stat uncertainty", channelCharge=args.charge)

    systList = {p : [] for p in processes}
    systLeg  = {p : [] for p in processes}
    systList_eta = {p : [] for p in processes}
    systList_pt  = {p : [] for p in processes}
    for p in processes:
        systList[p].append(unroll2Dto1D(nominals[p], newname=f"unrolled_{nominals[p].GetName()}", cropNegativeBins=False))
        systLeg[p].append(p)
        systList_eta[p].append(nominals[p].ProjectionX(f"{nominals[p].GetName()}_eta", 1, nominals[p].GetNbinsY(), "e"))
        systList_pt[p].append( nominals[p].ProjectionY(f"{nominals[p].GetName()}_pt",  1, nominals[p].GetNbinsX(), "e"))
    # now make a full loop for systematics
    for k in f.GetListOfKeys():
        name = k.GetName()
        if isWrem and not name.endswith(args.charge): continue
        # check name also allowing for perfect matching
        if not any(x in name for x in inputsysts) and not regexp_syst.match(name): continue

        if isWrem:
            tokens = name.split("_")
            pname = f"{tokens[1]}"
            if pname not in processes: continue
            sname = f"{'_'.join(tokens[2:-1])}"
            systLeg[pname].append(sname)
            sname += f"_{args.charge}"
        else:
            tokens = name.split("__") # remove "x_" and name of nuisance
            #print(tokens)
            pname = tokens[0].lstrip("x_")
            if pname not in processes: continue
            sname = tokens[1]
            systLeg[pname].append(sname)

        alternate = f.Get(name)
        alternate.SetDirectory(0)    
        systList[pname].append(unroll2Dto1D(alternate, newname=f"unrolled_{alternate.GetName()}", cropNegativeBins=False))
        systList_eta[pname].append(alternate.ProjectionX(f"{alternate.GetName()}_eta", 1, alternate.GetNbinsY(), "e"))
        systList_pt[pname].append( alternate.ProjectionY(f"{alternate.GetName()}_pt",  1, alternate.GetNbinsX(), "e"))
        ratio = alternate.Clone(f"systRatio_{sname}")
        if args.plotSystOriginal:
            if do2D:
                drawCorrelationPlot(alternate, "Muon #eta", "Muon p_{T} (GeV)", "Events",
                                    f"{name}_Syst", plotLabel="ForceTitle", outdir=outdir,
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                                    palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette,
                                    passCanvas=canvas, drawOption="COLZ0")
            if doUnrolled:
                plotUnrolledHistogram(systList[pname][-1], pname, sname+"_Syst", outdir, canvas_unroll, ratio, errorBars=args.addErrorBars, yAxisTitle="Events", channelCharge=args.charge)
            
        ratio.SetDirectory(0)
        #ratios[pname].append(htmp.Divide(nominals[pname]))
        ratio.Divide(nominals[pname])
        ratio.SetTitle(f"syst: {sname}")
        if do2D:
            drawCorrelationPlot(ratio, "Muon #eta", "Muon p_{T} (GeV)", f"{pname}: syst / nominal",
                                name, plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                                palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette,
                                passCanvas=canvas, drawOption="COLZ0")
        if doUnrolled:
            ratio_unrolled = unroll2Dto1D(ratio, newname=f"unrolledRatio_{name}", cropNegativeBins=False)
            plotUnrolledHistogram(ratio_unrolled, pname, sname, outdir, canvas_unroll, ratio, errorBars=args.addErrorBars, channelCharge=args.charge, canvasFullName=ratio_unrolled.GetName())

            if args.statUncRatio:
                statUncNomi = nominals[pname].Clone(f"statUnc_{nominals[pname].GetName()}")
                fillTH2fromTH2part(statUncNomi, nominals[pname], fillWithError=True)
                statUncAlt = alternate.Clone(f"statUnc_{alternate.GetName()}")
                fillTH2fromTH2part(statUncAlt, alternate, fillWithError=True)
                ratioStatUnc = statUncAlt.Clone(f"ratioStatUncSystOverNomi_{sname}")
                ratioStatUnc.Divide(statUncNomi)
                ratioStatUnc_unrolled = unroll2Dto1D(ratioStatUnc, newname=f"unrolledStatUncRatio_{name}", cropNegativeBins=False)
                plotUnrolledHistogram(ratioStatUnc_unrolled, pname, f"{sname}_statUncRatioWithNomi", outdir, canvas_unroll, ratioStatUnc, errorBars=False, yAxisTitle="stat. uncertainty ratio: syst/nomi", channelCharge=args.charge)
    print()

    ptBinRanges = []
    for i in range(nominals[p].GetNbinsY()):
        ptBinRanges.append("") # keep dummy otherwise there's too much text here
    systPostfix = f"_{args.charge}"
    if args.systPostfix:
        systPostfix += f"_{args.systPostfix}"
    for p in processes:
        if len(systList[p]) > 11:
            print("Not running drawNTH1() function to draw curves, there are too many lines ({})".format(len(systList[p])))
        else:
            drawNTH1(systList[p], systLeg[p], "Unrolled eta-p_{T} bin", "Events", f"nominalAndSyst_{p}{systPostfix}", outdir,
                     leftMargin=0.06, rightMargin=0.01, labelRatioTmp="syst/nomi",
                     legendCoords="0.1,0.9,0.82,0.9;5", lowerPanelHeight=0.5, skipLumi=True, passCanvas=canvas_unroll,
                     drawVertLines="{a},{b}".format(a=nominals[p].GetNbinsY(),b=nominals[p].GetNbinsX()),
                     textForLines=ptBinRanges, transparentLegend=False,
                     onlyLineColor=True, noErrorRatioDen=True, useLineFirstHistogram=True, setOnlyLineRatio=True, lineWidth=1)
            #
            drawNTH1(systList_eta[p], systLeg[p], "Muon eta", "Events", f"nominalAndSyst_{p}{systPostfix}_projEta", outdir,
                     topMargin=0.14, leftMargin=0.16, rightMargin=0.05, labelRatioTmp="syst/nomi",
                     legendCoords="0.02,0.98,0.91,0.99;3", lowerPanelHeight=0.5, skipLumi=True, passCanvas=canvas1D,
                     transparentLegend=False,
                     onlyLineColor=True, noErrorRatioDen=True, useLineFirstHistogram=True, setOnlyLineRatio=True, lineWidth=2)
            #
            drawNTH1(systList_pt[p], systLeg[p], "Muon p_{T} (GeV)", "Events", f"nominalAndSyst_{p}{systPostfix}_projPt", outdir,
                     topMargin=0.14, leftMargin=0.16, rightMargin=0.05, labelRatioTmp="syst/nomi",
                     legendCoords="0.02,0.98,0.91,0.99;3", lowerPanelHeight=0.45, skipLumi=True, passCanvas=canvas1D,
                     transparentLegend=False,
                     onlyLineColor=True, noErrorRatioDen=True, useLineFirstHistogram=True, setOnlyLineRatio=True, lineWidth=2)
    print()
