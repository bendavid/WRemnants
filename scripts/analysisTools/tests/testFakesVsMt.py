#!/usr/bin/env python3

# recent examples from THn (change '-c plus' to '-c both' to do both charges in one go)
#
# python w-mass-13TeV/testFakesVsMt.py plots/testNanoAOD/WmassPlots/histograms_isoChargeMtPtEta_fakeRegion_deepMET_NanoV9/nominal/postVFP//plots_fakerate.root plots/testNanoAOD/WmassPlots/histograms_isoChargeMtPtEta_fakeRegion_deepMET_NanoV9/nominal/postVFP/testFakesVsMt_rebin_eta8_pt2/ --palette 87 --rebin-x 8 --rebin-y 2 --mt-bin-edges "0,10,20,30,40,50,60" --mt-nominal-range "0,40" --mt-value-correction 53.0 -c plus


# python w-mass-13TeV/testFakesVsMt.py plots/testNanoAOD/WmassPlots/histograms_isoChargeMtPtEta_fakeRegion/nominal/postVFP//plots_fakerate.root plots/testNanoAOD/WmassPlots/histograms_isoChargeMtPtEta_fakeRegion/nominal/postVFP/testFakesVsMt_rebin_eta4_pt2/ --palette 87 --rebin-x 4 --rebin-y 2 --mt-bin-edges "0,10,20,30,40,50,60" --mt-nominal-range "0,40" --mt-value-correction 55.0 -c plus

# older exmples

# example for 2 mT bins with border at 40 GeV, and validation in signal region (needs existing plots)
# python w-mass-13TeV/testFakesVsMt.py plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_lowIso//plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_highIso//plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/testFakesVsMt_0to40to60/ --palette 87 --rebin-x 3 --rebin-y 5 --mt-bin-edges 0,40,60 --test-file plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/fakeRateRegion_postVFP_plus_systTH3/postprocessing/distributions_signalRegion/plots.root

# example with correction from extrapolated linear fit (charge plus)
# python w-mass-13TeV/testFakesVsMt.py plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_lowIso/plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion//fakeRateRegion_postVFP_plus_highIso/plots_fakerate.root plots/testNanoAOD/WmassPlots_MtPtEta_fakeRegion/testFakesVsMt_chargePlus/ --palette 87 --rebin-x 3 --rebin-y 5 --mt-bin-edges "0,10,20,30,40,50,60" --mt-nominal-range "0,40" --mt-value-correction 55.0 --test-file plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/fakeRateRegion_postVFP_plus_systTH3/postprocessing/distributions_signalRegion/plots.root

import os, re, array, math
import time
import argparse

import narf
import wremnants
import hist
import lz4.frame, pickle
from wremnants.datasets.datagroups import datagroups2016
from wremnants import histselections as sel

from utilities import boostHistHelpers as hh,common

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
from scripts.analysisTools.tests.testPlots1D import plotDistribution1D

# will use as input the histogram mTStudyForFakes with eta-pt-charge-mt-passIso-hasJet, where mt is the full distribution

# function to plot 1D data/MC distributions for the iso/nJet bins we have in this study
def plotProjection1D(rootHists, datasets, outfolder_dataMC, canvas1Dshapes=None, chargeBin=1,
                     projectAxisToKeep=0, isoAxisRange=[1,1], jetAxisRange=[1,2],
                     xAxisName="variable", plotName="variable_failIso_jetInclusive", mTaboveThis=None,
                     rebinVariable=None):

    firstBinMt = 1
    lastBinMt = rootHists["Data"].GetAxis(3).GetNbins()
    if mTaboveThis:
        firstBinMt = rootHists["Data"].GetAxis(3).FindFixBin(mTaboveThis+0.001)

    hdata = None
    hmc = {}
    for d in datasets:
        rootHists[d].GetAxis(2).SetRange(chargeBin, chargeBin)
        rootHists[d].GetAxis(3).SetRange(firstBinMt, lastBinMt)
        rootHists[d].GetAxis(4).SetRange(isoAxisRange[0], isoAxisRange[1])
        rootHists[d].GetAxis(5).SetRange(jetAxisRange[0], jetAxisRange[1])
        if d == "Data":
            hdata = rootHists[d].Projection(projectAxisToKeep, "EO")
            hdata.SetName(f"{plotName}_{d}")
            if rebinVariable:
                hdata.Rebin(rebinVariable)
        else:
            hmc[d] = rootHists[d].Projection(projectAxisToKeep, "EO")
            hmc[d].SetName(f"{plotName}_{d}")
            cropNegativeContent(hmc[d])  # should I crop the 2D before projecting?
            hmc[d].SetFillColor(colors_plots_[d])
            hmc[d].SetLineColor(ROOT.kBlack)
            hmc[d].SetMarkerSize(0)
            hmc[d].SetMarkerStyle(0)
            if rebinVariable:
                hmc[d].Rebin(rebinVariable)

    plotDistribution1D(hdata, hmc, datasets,
                       outfolder_dataMC, canvas1Dshapes=canvas1Dshapes,
                       xAxisName=xAxisName, plotName=plotName)

def plotProjection1Dfrom3D(rootHists, datasets, outfolder_dataMC, canvas1Dshapes=None, chargeBin=1,
                           projectAxisToKeep=0, xAxisName="variable", plotName="variable_failIso_jetInclusive",
                           correctionFakes=None):

    hdata = None
    hmc = {}
    axisProj = "x" if projectAxisToKeep == 0 else "y" if projectAxisToKeep == 1 else "z"
    for d in datasets:
        rootHists[d].GetZaxis().SetRange(chargeBin, chargeBin)
        if d == "Data":
            hdata = rootHists[d].Project3D(f"{axisProj}eo")
            hdata.SetName(f"{plotName}_{d}")
        else:
            if correctionFakes and d == "Fake":
                # correct 2D before going to 1D
                hmc2D = rootHists[d].Project3D(f"yxe")
                cropNegativeContent(hmc2D)
                hmc2D.SetName(f"{plotName}_{d}_2Dtmp")                                                 
                print()
                print(f"PlotName: {plotName}")
                print(f"Fakes norm before corr: {hmc2D.Integral()}")
                scaleTH2byOtherTH2(hmc2D, correctionFakes, scaleUncertainty=True, zeroOutOfAcceptance=False)
                print(f"Fakes norm after corr : {hmc2D.Integral()}")
                print()
                if axisProj == "x":
                    hmc[d] = hmc2D.ProjectionX(f"{plotName}_{d}", 1, hmc2D.GetNbinsX(), "eo")
                elif axisProj == "y":
                    hmc[d] = hmc2D.ProjectionY(f"{plotName}_{d}", 1, hmc2D.GetNbinsY(), "eo")
                else:
                    print("Error in plotProjection1Dfrom3D: 'z' axis not permitted for TH2. Abort")
                    quit()
            else:
                hmc[d] = rootHists[d].Project3D(f"{axisProj}eo")
                hmc[d].SetName(f"{plotName}_{d}")
            #cropNegativeContent(hmc[d]) # should I crop the 2D before projecting?
            hmc[d].SetFillColor(colors_plots_[d])
            hmc[d].SetLineColor(ROOT.kBlack)
            hmc[d].SetMarkerSize(0)
            hmc[d].SetMarkerStyle(0)
        
    plotDistribution1D(hdata, hmc, datasets,
                       outfolder_dataMC, canvas1Dshapes=canvas1Dshapes,
                       xAxisName=xAxisName, plotName=plotName)


def drawAndFitFRF(h1, labelXtmp="xaxis", labelYtmp="yaxis",
                  canvasName="default", outdir="./",
                  rebinFactorX=0,
                  draw_both0_noLog1_onlyLog2=1,
                  leftMargin=0.15,
                  rightMargin=0.04,
                  labelRatioTmp="Rel.Unc.::0.5,1.5",
                  drawStatBox=False,
                  legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                  canvasSize="600,700",  # use X,Y to pass X and Y size     
                  lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                  drawLineLowerPanel="luminosity uncertainty::0.025", # if not empty, draw band at 1+ number after ::, and add legend with title
                  passCanvas=None,
                  lumi=None,
                  moreTextLatex="",
                  fitRange="0,40", # xmin and xmax
                  fitOptions="MSQ+",
                  evalAt=None,
                  mTshape=None,
                  fitPolDegree=3,
                  useBinnedCorr=False
):

    # moreTextLatex is used to write text with TLatex, and the syntax is textToWrite::x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)
    # multiple text lines can be defined in textToWrite separated by ";", e.g. line1;line2;line3
    
    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
        # case in which rebinFactorX is a list of bin edges
        else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    #h1.SetMarkerSize(0)

    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        diff = ymax - ymin
        ymin -= diff * 0.25
        ymax += diff * 0.3
        if ymin < 0: ymin = 0

    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)

    # draw with bands
    # h1.Draw("HIST")
    # h1err = h1.Clone("h1err")
    # h1err.SetFillColor(ROOT.kRed+2)
    # h1err.SetFillStyle(3001)  # 3001 is better than 3002 for pdf, while 3002 is perfect for png
    # h1err.Draw("E2same")
    h1.Draw("EP")
    
    xMinFit, xMaxFit = map(float, fitRange.split(','))
    #print(f"Fit range [{xMinFit}, {xMaxFit}] ({xMaxFit-xMinFit} wide)")

    # tensorflow fit cannot select a subrange to fit, so need to make the boost histogram directly from
    # a new ROOT histograms defined in the actual fit range
    fitBinLow  = h1.GetXaxis().FindFixBin(xMinFit + 0.001)
    fitBinHigh = h1.GetXaxis().FindFixBin(xMaxFit + 0.001)
    fitBinEdges = [round(h1.GetXaxis().GetBinLowEdge(i),2) for i in range(fitBinLow, fitBinHigh+1)]
    #print(fitBinEdges)
    h1fitRange = ROOT.TH1D(f"{h1.GetName()}_fitRange", "", len(fitBinEdges)-1, array("d", fitBinEdges))
    for ib in range(1, 1+h1fitRange.GetNbinsX()):
        binOriginalHisto = ib + fitBinLow - 1
        h1fitRange.SetBinContent(ib, h1.GetBinContent(binOriginalHisto))
        h1fitRange.SetBinError(ib, h1.GetBinError(binOriginalHisto))
    boost_hist = narf.root_to_hist(h1fitRange)

    # OLD ROOT BASED FIT
    #f1 = ROOT.TF1("f1",polN_scaled_global, xMinFit, xMaxFit, 3)
    ## TODO: set coefficient of x^1 as 0 to have the maximum at 0?
    #f1.SetParLimits(2, -10.0, 0.0)
    #realFitOptions = fitOptions
    #if "B" not in fitOptions:
    #    realFitOptions = "B" + fitOptions
    #fitres = h1.Fit("f1", realFitOptions, "", xMinFit, xMaxFit)

    # FIT WITH TENSORFLOW
    global polN_scaled_global
    fitModel = f"pol{fitPolDegree}" # just for legend, set here to make sure it is changed consistently with line below
    params = np.array([1.0] + [0.0 for i in range(fitPolDegree)])
    if polN_scaled_global == None: 
        polN_scaled_global = partial(polN_root_, xLowVal=xMinFit, xFitRange=xMaxFit-xMinFit, degree=fitPolDegree)

    fitres_tf1 = narf.fit_hist(boost_hist, polN_scaled_global, params)
    f1 = ROOT.TF1("f1",polN_scaled_global, xMinFit, xMaxFit, len(params))
    status = fitres_tf1["status"]
    covstatus = fitres_tf1["covstatus"]
    if status:
        print(f"\n\n-----> STATUS: fit had status {status}\n\n")
        if status not in badFitsID.keys():
            badFitsID[status] = 1
        else:
            badFitsID[status] += 1                                       
    if covstatus:
        print(f"\n\n-----> COV: fit had covstatus {covstatus}\n\n")
        if covstatus not in badCovMatrixID.keys():
            badCovMatrixID[covstatus] = 1
        else:
            badCovMatrixID[covstatus] += 1                                       

    f1.SetParameters( np.array( fitres_tf1["x"], dtype=np.dtype('d') ) )
    f1.SetLineWidth(3)
    #f1.SetLineStyle(9) # ROOT.kDashed == thin dashes, almost dotted
    f1.SetLineColor(ROOT.kBlue+1)
    # draw extrapolation
    f2 = ROOT.TF1("f2",polN_scaled_global, xMaxFit, h1.GetXaxis().GetBinLowEdge(1+h1.GetNbinsX()), len(params))
    for i in range(1+fitPolDegree):
        f2.SetParameter(i, f1.GetParameter(i))
    f2.SetLineColor(ROOT.kRed+2)
    f2.SetLineWidth(3)
    f2.SetLineStyle(9)

    xvalMin = h1.GetXaxis().GetBinLowEdge(1)
    xvalMax = h1.GetXaxis().GetBinUpEdge(h1.GetNbinsX())
    funcFullRange = ROOT.TF1("funcFullRange", polN_scaled_global, xvalMin, 200, len(params)) # define this one until "infinity"
    funcFullRange.SetParameters( np.array( fitres_tf1["x"], dtype=np.dtype('d') ) )
    hband = ROOT.TH1D("hband", "", 400, xvalMin, xvalMax)
    hband.SetStats(0)
    hband.SetFillColor(ROOT.kGray)
    #hband.SetFillStyle(3001)
    for ib in range(1, hband.GetNbinsX()+1):
        xval = hband.GetBinCenter(ib)
        val = f1.Eval(xval)
        #print(f"val({xval}) = {val}")
        hband.SetBinContent(ib, val)

    npar = len(params)
    # diagonalize and get eigenvalues and eigenvectors
    e, v = np.linalg.eigh(fitres_tf1["cov"])
    # store all variations for faster access below
    altParameters = np.array([np.zeros(npar, dtype=np.dtype('d'))] * (npar * 2), dtype=np.dtype('d'))
    #print(altParameters)
    for ivar in range(npar):
        shift = np.sqrt(e[ivar]) * v[:, ivar]
        altParameters[ivar]      = fitres_tf1["x"] + shift
        altParameters[ivar+npar] = fitres_tf1["x"] - shift

    tf1_func_alt = ROOT.TF1()
    tf1_func_alt.SetName("tf1_func_alt")
    funcFullRange.Copy(tf1_func_alt)
    tf1_func_alt.SetLineWidth(2)

    for ib in range(1, hband.GetNbinsX()+1):
        xval = hband.GetBinCenter(ib)
        err = 0.0
        for ivar in range(npar):
            # set parameters for a given hessian
            # make band using only up variations, the result for down ones should be very symmetric 
            tf1_func_alt.SetParameters(altParameters[ivar]) # this is for Up variations
            funcVal = tf1_func_alt.Eval(xval)
            diff = funcVal - hband.GetBinContent(ib)
            err += diff * diff
        err = math.sqrt(err)
        hband.SetBinError(ib, err)

    hband.Draw("E2 SAME")
    f1.Draw("L SAME")
    if xMaxFit < xvalMax:
        f2.Draw("L SAME")
    h1.Draw("EP SAME")
        
    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    # resize if one line is not needed 
    if xMaxFit >= xvalMax:
        ly1 = ly2 - 0.75 * (ly2 - ly1)
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
        
    leg.SetNColumns(nColumnsLeg)
    leg.AddEntry(h1,"Measurement","PE")
    leg.AddEntry(f1,f"Fit {fitModel} in [{int(xMinFit)}, {int(xMaxFit)}]","L")
    if xMaxFit < xvalMax:
        leg.AddEntry(f2,f"Extrapolation","L")
    leg.AddEntry(hband,"Fit uncertainty","F")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        h1.SetStats(0)

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        ratio = h1.Clone("ratio")
        den_noerr = h1.Clone("den_noerr")
        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        ratio.Divide(den_noerr)
        ratio.SetFillColor(ROOT.kGray+1)
        #den_noerr.SetFillColor(ROOT.kGray)
        frame.Draw()
        ratio.SetMarkerSize(0)
        ratio.SetMarkerStyle(0) # important to remove dots at y = 1
        ratio.Draw("E2same")

        line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineColor(ROOT.kRed)
        line.SetLineWidth(1)
        line.Draw("Lsame")

        if drawLineLowerPanel:
            legEntry,yline = drawLineLowerPanel.split('::')
            line2 = ROOT.TF1("horiz_line_2",str(1+float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line3 = ROOT.TF1("horiz_line_3",str(1-float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line2.SetLineColor(ROOT.kBlue)
            line2.SetLineWidth(1)
            line2.Draw("Lsame")
            line3.SetLineColor(ROOT.kBlue)
            line3.SetLineWidth(1)
            line3.Draw("Lsame")
            x1leg2 = 0.2 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.5 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.25 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.35 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.AddEntry(line2,legEntry,"L")
            leg2.Draw("same")

        
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)

    # when deriving the average mT at which to evaluate the FRF, use the average of the FRF over the mT distribution
    # essentially integral of FRF(mT)*G(mT)*dmT, where G is the normalized shape of the mT distribution 

    averageFRF = 0.0
    if mTshape:
        ## get integral in restricted range to normalize the mT PDF
        #mTshape.GetXaxis().SetRange(1, mTshape.GetXaxis().FindFixBin(xvalMax-0.001))
        integralMt = mTshape.Integral()
        for ib in range(1, 1+mTshape.GetNbinsX()):
            mTbinVal = mTshape.GetBinCenter(ib)
            # evaluate FRF as 0 when mT is too large, the function would be negative most likely,
            # but if positive it should be set to 0 to avoid issues when the fitted parabola increases versus mT 
            # perhaps the parabola should be forced to have negative concavity
            #valFRF = 0.0 if mTbinVal > xvalMax else funcFullRange.Eval(mTbinVal)
            valFRF = h1.GetBinContent(max(1, min(h1.GetXaxis().FindFixBin(mTbinVal), h1.GetNbinsX()))) if useBinnedCorr else funcFullRange.Eval(mTbinVal) 
            averageFRF += max(0.0, valFRF) * mTshape.GetBinContent(ib)
        averageFRF /= integralMt
        
    if evalAt:
        return averageFRF, funcFullRange.Eval(evalAt)
    else:
        return averageFRF, 0.0
        
################################################################

# need to define functions as global below to avoid them being deleted out of fitTurnOnTF
polN_scaled_global = None
badFitsID = {}
badCovMatrixID = {}


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1)
    parser.add_argument("outputfolder",   type=str, nargs=1)
    parser.add_argument("-x", "--x-axis-name", dest="xAxisName", default="Muon #eta", help="x axis name")
    parser.add_argument("-y", "--y-axis-name", dest="yAxisName", default="Muon p_{T} (GeV)", help="y axis name")
    parser.add_argument("-z", "--z-axis-name", dest="zAxisName", default="m_{T} (GeV)", help="z axis name")
    parser.add_argument("-c", "--charge", default="plus", choices=["plus", "minus", "both"], help="charge")
    parser.add_argument("--mt-bin-edges", dest="mtEdges", default="0,10,20,30,40,50,60", type=str, help="Comma-separated list of bin edges for mT")
    parser.add_argument("--mt-nominal-range", dest="mtNominalRange", default="0,40", type=str, help="Comma-separated list of 2 bin edges for mT, representing the nominal range, used to derive the correction using also option --mt-value-correction")
    parser.add_argument("--mt-fit-range", dest="mtFitRange", default="0,50", type=str, help="Comma-separated list of 2 bin edges for mT, representing the fit range, might be the same as --mt-nominal-range but not necessarily")
    parser.add_argument(     '--fit-pol-degree'  , dest='fitPolDegree',      default=3, type=int, help='Degree for polynomial used in the fits')
    parser.add_argument("--integral-mt-method", dest="integralMtMethod", default="sideband", choices=["sideband", "fullRange"], type=str, help="How to integrate mT distribution to derive the FRF correction (default is 'sideband', orthogonal to the nominal fake region, 'fullRange' uses the full range)")
    parser.add_argument(     '--use-binned-correction', dest='useBinnedCorrection', action='store_true',   help='Use binned FRF to derive the correction (mainly for tests)')
    parser.add_argument(     '--jet-cut', dest='jetCut', action='store_true',   help='Use jet cut to derive the FRF (sample will be more QCD enriched but might boas the FRF)')
    parser.add_argument(     "--rebin-x", dest="rebinEta", default=1, type=int, help="To rebin x axis (eta)")
    parser.add_argument(     "--rebin-y", dest="rebinPt", default=1, type=int, help="To rebin y axis (pt)")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--absolute-eta', dest='doAbsEta' , action='store_true',   help='Do study using histograms folded into absolute value of pseudorapidity')
    parser.add_argument(     '--skip-plot-2D', dest='skipPlot2D' , action='store_true',   help='skip some 2D plots with FRF vs eta-Mt, pt-mT, and so on')
    parser.add_argument("--postfix", default="", type=str, help="Postfix for folder name")
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");
    args = parser.parse_args()
    
    logger = common.setup_color_logger(os.path.basename(__file__), args.verbose)
    if args.doAbsEta:
        logger.error("Option --absolute-eta not implemented correctly yet. Abort")
        quit()
    
    ROOT.TH1.SetDefaultSumw2()

    if args.charge == "both":
        logger.warning("Running both charges together with -c both is currently deprecated")
        logger.warning("There is some unexpected rebinning versus mt to be fixed")
        quit()
        charges = ["plus", "minus"]
    else:
        charges = [args.charge]

    subFolder = f"/{args.integralMtMethod}Integral_fit{args.mtFitRange.replace(',','to')}_pol{args.fitPolDegree}"
    if args.rebinEta > 1:
        subFolder += f"_rebinEta{args.rebinEta}"
    if args.rebinPt > 1:
        subFolder += f"_rebinPt{args.rebinPt}"
    if args.jetCut:
        subFolder += "_1orMoreJetFRF"
    else:
        subFolder += "_jetInclusiveFRF"
    if args.postfix:
        subFolder += f"_{args.postfix}"
    subFolder += "/"
    mainOutputFolder = args.outputfolder[0] + subFolder
        
    hFakerateFactorCorrection = {x: None for x in charges}
    hFRF_proj = {x: None for x in charges}

    etaLabel = "#eta" if not args.doAbsEta else "|#eta|"
    
    xAxisName = args.xAxisName
    yAxisName = args.yAxisName
    zAxisName = args.zAxisName

    groups = datagroups2016(args.inputfile[0], applySelection=False)
    datasets = groups.getNames()
    datasetsNoQCD = list(filter(lambda x: x != "QCD", datasets)) # exclude QCD MC if present
    datasetsNoFakes = list(filter(lambda x: x != "Fake", datasets)) 
    datasetsNoQCDFakes = list(filter(lambda x: x not in ["QCD", "Fake"], datasets))
    logger.info(f"Will plot datasets {datasets}")
    inputHistName = "mTStudyForFakes"
    groups.setNominalName(inputHistName)
    groups.loadHistsForDatagroups(inputHistName, syst="", procsToRead=datasets, applySelection=False)
    histInfo = groups.getDatagroups()
    rootHists = {d: None for d in datasets}
    for d in datasets:
        #print(d)
        hnarf = histInfo[d][inputHistName]
        rootHists[d] = narf.hist_to_root(hnarf) # this is a THnD with eta-pt-charge-mt-passIso
        
    ########
    ########
    # data-MC already done in the management of the groups above
    # should check that the subtraction by hand yields the same results and uncertainties
    # the other processes are still needed to make other plots with data and MC, like mT in the different regions
    histo_fakes = copy.deepcopy(rootHists["Fake"])    

    # get the standard signal region with the data-driven fakes,
    # it will be needed to make the plot with/without the corrections
    inputHistName = "nominal"
    groups.setNominalName(inputHistName)
    groups.setSelectOp(sel.histWmass_passMT_passIso)
    groups.setSelectOp(sel.fakeHistABCD, processes=["Fake"])
    groups.loadHistsForDatagroups(inputHistName, syst="", procsToRead=datasetsNoQCD)
    histInfo_nominal = groups.getDatagroups()
    rootHists_nominal = {d: None for d in datasetsNoQCD}
    for d in rootHists_nominal.keys():
        hnarf = histInfo_nominal[d][inputHistName]
        rootHists_nominal[d] = narf.hist_to_root(hnarf) # this is a TH3D with eta-pt-charge
    #######
    
    adjustSettings_CMS_lumi()    
    canvas = ROOT.TCanvas("canvas","",800,800)
    canvas1D = ROOT.TCanvas("canvas1D","",800,700)
    canvas1Dshapes = ROOT.TCanvas("canvas1Dshapes","",700,800)

    axisVar = {0 : ["muon_eta", "#eta"],
               1 : ["muon_pt",  "p_{T} (GeV)"],
               3 : ["mT", "m_{T} (GeV)"]
    }
               
    for charge in charges:

        outfolder = f"{mainOutputFolder}/{charge}/"
        createPlotDirAndCopyPhp(outfolder)

        outfolder_dataMC = f"{outfolder}/shapesDataMC/"
        createPlotDirAndCopyPhp(outfolder_dataMC)

        # bin number from root histogram
        chargeBin = 1 if charge == "minus" else 2

        # plot mT, eta, pt in some regions iso-nJet regions, for checks
        # don't plot fakes here
        for xbin in axisVar.keys():
            rebinVariable = 2 if xbin == 3 else None 
            #if xbin == 3:
            #    rebinVariable=2 # mT has 1 GeV original binning, but for plot 2 GeV might be better
            plotProjection1D(rootHists, datasetsNoFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                             projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1], plotName=f"{axisVar[xbin][0]}_passIso_1orMoreJet",
                             isoAxisRange=[2,2], jetAxisRange=[2,2], rebinVariable=rebinVariable)
            plotProjection1D(rootHists, datasetsNoFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                             projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1], plotName=f"{axisVar[xbin][0]}_passIso_jetInclusive",
                             isoAxisRange=[2,2], jetAxisRange=[1,2], rebinVariable=rebinVariable)
            plotProjection1D(rootHists, datasetsNoFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                             projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1], plotName=f"{axisVar[xbin][0]}_failIso_1orMoreJet",
                             isoAxisRange=[1,1], jetAxisRange=[2,2], rebinVariable=rebinVariable)
            plotProjection1D(rootHists, datasetsNoFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                             projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1], plotName=f"{axisVar[xbin][0]}_failIso_jetInclusive",
                             isoAxisRange=[1,1], jetAxisRange=[1,2], rebinVariable=rebinVariable)
            # signal region adding mT cut too
            plotProjection1D(rootHists, datasetsNoFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                             projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1], plotName=f"{axisVar[xbin][0]}_passIso_jetInclusive_passMt",
                             isoAxisRange=[2,2], jetAxisRange=[1,2], mTaboveThis=40, rebinVariable=rebinVariable)
            
        ###################################
        ###################################
        ###   Now the actual study on fakes
        ###################################
        ###################################

        # select events with at least one jet (or jet inclusive depending on what you need to do)
        # jetInclusive: GetAxis(5).SetRange(1, 2)
        # 1orMoreJet:   GetAxis(5).SetRange(2, 2)
        
        # set charge from charge axis
        histo_fakes.GetAxis(2).SetRange(chargeBin, chargeBin)
        histo_fakes.GetAxis(4).SetRange(2, 2) # passIso, equivalent to lowIso
        histo_fakes.GetAxis(5).SetRange(2 if args.jetCut else 1, 2)
        # now get a TH3
        histoPassIso = histo_fakes.Projection(0, 1, 3, "E")
        histoPassIso.SetName("fakes_passIso")
        histoPassIso.SetTitle("fakes_passIso")
        histo_fakes.GetAxis(2).SetRange(chargeBin, chargeBin)
        histo_fakes.GetAxis(4).SetRange(1, 1) # FailIso
        histo_fakes.GetAxis(5).SetRange(2 if args.jetCut else 1, 2)
        histoFailIso = histo_fakes.Projection(0, 1, 3, "E")
        histoFailIso.SetName("fakes_failIso")
        histoFailIso.SetTitle("fakes_failIso")

        # to get the mean only in the desired mt range, get projection using only mt bins in the high mt region
        # this is always done in the jet inclusive case
        mtThreshold = float(args.mtNominalRange.split(",")[-1])
        histo_fakes.GetAxis(2).SetRange(chargeBin, chargeBin)
        histo_fakes.GetAxis(3).SetRange(histo_fakes.GetAxis(3).FindFixBin(mtThreshold+0.001),
                                        histo_fakes.GetAxis(3).GetNbins()) # high mT
        histo_fakes.GetAxis(4).SetRange(1, 1) # FailIso
        histo_fakes.GetAxis(5).SetRange(1, 2) # jet inclusive
        histoPassMtFailIso = histo_fakes.Projection(0, 1, 3, "E")
        histoPassMtFailIso.SetName("fakes_passMt_failIso_jetInclusive")
        histoPassMtFailIso.SetTitle("fakes_passMt_failIso_jetInclusive")

        # don't rebin mT here (default is 1 GeV),
        # otherwise the new binning might no longer be consistent with the one passed to --mt-bin-edges
        histoPassIso.Rebin3D(args.rebinEta, args.rebinPt, 1)
        histoFailIso.Rebin3D(args.rebinEta, args.rebinPt, 1)
        histoPassMtFailIso.Rebin3D(args.rebinEta, args.rebinPt, 1)

        cropNegativeContent(histoPassIso)
        cropNegativeContent(histoFailIso)
        cropNegativeContent(histoPassMtFailIso)

        mtEdges = [round(int(x),1) for x in args.mtEdges.split(',')] 
        nMtBins = len(mtEdges) -1
        ratio = []
        for imt in range(nMtBins):
            lowEdge = mtEdges[imt]
            highEdge = mtEdges[imt+1]
            binStart = histoPassIso.GetZaxis().FindFixBin(lowEdge)
            binEnd = histoPassIso.GetZaxis().FindFixBin(highEdge+0.001) - 1 # bin up edges belong to "next" bin
            h2PassIso = getTH2fromTH3(histoPassIso, f"pt_eta_mt{lowEdge}to{highEdge}_passIso", binStart, binEnd)
            h2FailIso = getTH2fromTH3(histoFailIso, f"pt_eta_mt{lowEdge}to{highEdge}_failIso", binStart, binEnd)

            h2PassIso.SetTitle("Low Iso: m_{T} #in [%d, %d]" % (lowEdge, highEdge))
            if not args.skipPlot2D:
                drawCorrelationPlot(h2PassIso,
                                    xAxisName,
                                    yAxisName,
                                    "Events (data - MC)",
                                    h2PassIso.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)

            h2FailIso.SetTitle("High Iso: m_{T} #in [%d, %d]" % (lowEdge, highEdge))
            if not args.skipPlot2D:
                drawCorrelationPlot(h2FailIso,
                                    xAxisName,
                                    yAxisName,
                                    "Events (data - MC)",
                                    h2FailIso.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
                
            ratio.append(h2PassIso.Clone(f"fakerateFactor_mt{lowEdge}to{highEdge}"))
            ratio[imt].SetTitle("m_{T} #in [%d, %d]" % (lowEdge, highEdge))
            ratio[imt].Divide(h2FailIso)
            if not args.skipPlot2D:
                drawCorrelationPlot(ratio[imt],
                                    xAxisName,
                                    yAxisName,
                                    "fakerate factor: N(iso) / N(not-iso)::0,3",
                                    ratio[imt].GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)


        if args.mtNominalRange:
            lowEdge, highEdge = map(int, args.mtNominalRange.split(','))
            binStart = histoPassIso.GetZaxis().FindFixBin(lowEdge)
            binEnd = histoPassIso.GetZaxis().FindFixBin(highEdge+0.001) - 1 # bin up edges belong to "next" bin
            h2PassIso = getTH2fromTH3(histoPassIso, f"pt_eta_nominalmt{lowEdge}to{highEdge}_passIso", binStart, binEnd)
            h2FailIso = getTH2fromTH3(histoFailIso, f"pt_eta_nominalmt{lowEdge}to{highEdge}_failIso", binStart, binEnd)
            cropNegativeContent(h2PassIso)
            cropNegativeContent(h2FailIso)
            nominalFakerateFactor = h2PassIso.Clone(f"nominalFakerateFactor_mt{lowEdge}to{highEdge}")
            nominalFakerateFactor.SetTitle("m_{T} #in [%d, %d]" % (lowEdge, highEdge))        
            nominalFakerateFactor.Divide(h2FailIso)

        # rounding with python is sometimes tricky, add an epsilon to get the desired number
        etaLow = round(0.01 + h2PassIso.GetXaxis().GetBinLowEdge(1), 1)
        etaHigh = round(0.01 + h2PassIso.GetXaxis().GetBinLowEdge(1+h2PassIso.GetNbinsX()), 1)
        ptLow = round(0.01 + h2PassIso.GetYaxis().GetBinLowEdge(1), 1)
        ptHigh = round(0.01 + h2PassIso.GetYaxis().GetBinLowEdge(1+h2PassIso.GetNbinsY()), 1)

        hFakerateFactorCorrection[charge] = ROOT.TH2D(f"fakerateFactorCorrection_{charge}", "m_{T} > %d GeV" % int(args.mtNominalRange.split(',')[1]),
                                                      h2PassIso.GetNbinsX(), round(etaLow,1), round(etaHigh,1),
                                                      h2PassIso.GetNbinsY(), round(ptLow,1), round(ptHigh,1))
        
        # now preparing a summary for each pt bin
        ptCentralBin = h2PassIso.GetYaxis().FindFixBin(39.5)
        for ipt in range(1, 1+h2PassIso.GetNbinsY()):
            ptBinLow = int(h2PassIso.GetYaxis().GetBinLowEdge(ipt))
            ptBinHigh = int(h2PassIso.GetYaxis().GetBinLowEdge(ipt+1))
            fakerateFactor_vs_etaMt = ROOT.TH2D("fakerateFactor_vs_etaMt_pt%dto%d" % (ptBinLow, ptBinHigh),
                                                "Muon p_{T} #in [%d, %d] GeV" % (ptBinLow, ptBinHigh),
                                                h2PassIso.GetNbinsX(), round(etaLow,1), round(etaHigh,1),
                                                nMtBins, array("d", mtEdges)
                                               )

            outfolder1D = outfolder + "fakerateFactor_fits_pt%dto%d/" % (ptBinLow, ptBinHigh)
            createPlotDirAndCopyPhp(outfolder1D)


            for ieta in range(1, 1+fakerateFactor_vs_etaMt.GetNbinsX()):

                etaBinLowNoRound = fakerateFactor_vs_etaMt.GetXaxis().GetBinLowEdge(ieta)
                etaBinHighNoRound = fakerateFactor_vs_etaMt.GetXaxis().GetBinLowEdge(ieta+1)
                etaBinLow =  round(0.01 + etaBinLowNoRound, 1)
                etaBinHigh = round(0.01 + etaBinHighNoRound, 1)
                # print(f"ieta = {ieta}    {ptBinLow} < pt < {ptBinHigh}     {etaBinLow} < eta < {etaBinHigh}    {etaBinLow} < etaNoRound < {etaBinHigh}")
                hFRfactorVsMt = ROOT.TH1D(f"hFRfactorVsMt_ieta{ieta}_pt{ptBinLow}to{ptBinHigh}",
                                          "%.1f < %s < %.1f, p_{T} #in [%d, %d] GeV" % (etaBinLow, etaLabel, etaBinHigh, ptBinLow, ptBinHigh),
                                          nMtBins, array("d", mtEdges))

                # to make easier computation of correction factor below
                hTmp = []

                for imt in range(1, 1+fakerateFactor_vs_etaMt.GetNbinsY()):
                    binContent = ratio[imt-1].GetBinContent(ieta, ipt)
                    binError = ratio[imt-1].GetBinError(ieta, ipt)
                    fakerateFactor_vs_etaMt.SetBinContent(ieta, imt, binContent)
                    fakerateFactor_vs_etaMt.SetBinError(ieta, imt, binError)
                    hFRfactorVsMt.SetBinContent(imt, binContent)
                    hFRfactorVsMt.SetBinError(  imt, binError)
                    if nMtBins == 2:
                        hTmp.append(ROOT.TH1D(f"hTmp{imt}","",1,0,1))
                        hTmp[imt-1].SetBinContent(1, max(0.0, binContent))
                        hTmp[imt-1].SetBinError(  1, binError)

                textLatex = "%.1f < %s < %.1f;p_{T} #in [%d, %d] GeV::0.2,0.3,0.1,0.045" % (etaBinLow, etaLabel, etaBinHigh, ptBinLow, ptBinHigh)
                if nMtBins > 2:
                    ## get high mt value to evaluate the correction, can use the mean of the mt distribution for each etapt bin 
                    if args.integralMtMethod == "sideband":
                        projMt = histoPassMtFailIso.ProjectionZ(f"projZ_{histoPassMtFailIso.GetName()}", ieta, ieta, ipt, ipt, "eo")
                    elif args.integralMtMethod == "fullRange":
                        ## or maybe in the full range to get the actual average
                        projMt = histoFailIso.ProjectionZ(f"projZ_{histoFailIso.GetName()}", ieta, ieta, ipt, ipt, "eo")
                    cropNegativeContent(projMt)
                    meanMtInThisEtaPtBin = projMt.GetMean()                        
                    valFRFhighMt, valFRFhighMtAtMean = drawAndFitFRF(hFRfactorVsMt, zAxisName, "Fakerate factor: N(iso) / N(not-iso)",
                                                                     hFRfactorVsMt.GetName(),
                                                                     outfolder1D, lowerPanelHeight=0.0, passCanvas=canvas1D,
                                                                     moreTextLatex=textLatex,
                                                                     legendCoords="0.64,0.96,0.69,0.93", fitRange=args.mtFitRange,
                                                                     fitOptions="MSQ+",
                                                                     evalAt=meanMtInThisEtaPtBin, mTshape=projMt,
                                                                     fitPolDegree=args.fitPolDegree,
                                                                     useBinnedCorr=args.useBinnedCorrection)
                    print(f"{valFRFhighMt}, {valFRFhighMtAtMean}")
                    if valFRFhighMt < 0:
                        printLine(marker=" ")
                        printLine()
                        logger.warning(f"Warning: ieta = {ieta},   ipt = {ipt},   FRF = {valFRFhighMt}")
                        #print(f"Warning: ieta = {ieta},   ipt = {ipt},   FRF(mt={meanMtInThisEtaPtBin}) = {valFRFhighMt}")
                        logger.warning("Setting FRF to 0.01!")
                        printLine()
                        printLine(marker=" ")
                        valFRFhighMt = 0.01
                    hFakerateFactorCorrection[charge].SetBinContent(ieta, ipt, valFRFhighMt / nominalFakerateFactor.GetBinContent(ieta, ipt))
                elif nMtBins == 2:
                    hTmp[1].Divide(hTmp[0])
                    hFakerateFactorCorrection[charge].SetBinContent(ieta, ipt, hTmp[1].GetBinContent(1))
                    hFakerateFactorCorrection[charge].SetBinError(  ieta, ipt, hTmp[1].GetBinError(1))
                    drawSingleTH1(hFRfactorVsMt, zAxisName, "Fakerate factor: N(iso) / N(not-iso)", hFRfactorVsMt.GetName(),
                                  outfolder1D, lowerPanelHeight=0.0, passCanvas=canvas1D, moreTextLatex=textLatex,
                                  legendCoords="0.64,0.96,0.77,0.93")

            if not args.skipPlot2D:
                drawCorrelationPlot(fakerateFactor_vs_etaMt,
                                    xAxisName,
                                    zAxisName,
                                    "fakerate factor: N(iso) / N(not-iso)",
                                    fakerateFactor_vs_etaMt.GetName(), plotLabel="ForceTitle", outdir=outfolder,
                                    draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                    invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)

        if args.integralMtMethod == "sideband":
            if args.jetCut:
                corrRange = "0,1.0"
            else:
                corrRange = "0.2,1.2"
        else:
            minCorr = hFakerateFactorCorrection[charge].GetBinContent(hFakerateFactorCorrection[charge].GetMinimumBin())
            maxCorr = hFakerateFactorCorrection[charge].GetBinContent(hFakerateFactorCorrection[charge].GetMaximumBin())
            diff = max(abs(maxCorr-1.0), abs(minCorr-1.0))
            corrRange = f"{1.0-diff},{1.0+diff}"
        drawCorrelationPlot(hFakerateFactorCorrection[charge],
                            xAxisName,
                            yAxisName,
                            f"QCD template correction::{corrRange}",
                            hFakerateFactorCorrection[charge].GetName(), plotLabel="ForceTitle", outdir=outfolder,
                            draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                            invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
        # also plot average vs pt and eta
        # currently doesn't plot as I expect
        # hFRF_profileX = hFakerateFactorCorrection[charge].ProfileX(hFakerateFactorCorrection[charge].GetName()+"_profileX",
        #                                                            1, hFakerateFactorCorrection[charge].GetNbinsY())
        # hFRF_profileY = hFakerateFactorCorrection[charge].ProfileY(hFakerateFactorCorrection[charge].GetName()+"_profileY",
        #                                                            1, hFakerateFactorCorrection[charge].GetNbinsX())
        # drawTH1(hFRF_profileX, xAxisName, f"QCD template correction", hFRF_profileX.GetName(), outfolder, drawStatBox=False,
        #         plotTitleLatex=f"Charge {charge}")
        # drawTH1(hFRF_profileY, yAxisName, f"QCD template correction", hFRF_profileY.GetName(), outfolder, drawStatBox=False,
        #         plotTitleLatex=f"Charge {charge}")
        if hFakerateFactorCorrection[charge].GetNbinsX() > 1 and hFakerateFactorCorrection[charge].GetNbinsY() > 1:
            drawCorrelationPlot(hFakerateFactorCorrection[charge],
                                xAxisName,
                                yAxisName,
                                f"QCD template correction::{corrRange}",
                                f"{hFakerateFactorCorrection[charge].GetName()}_smooth", plotLabel="ForceTitle", outdir=outfolder, smoothPlot=True,
                                draw_both0_noLog1_onlyLog2=1, nContours=args.nContours, palette=args.palette,
                                invertePalette=args.invertePalette, passCanvas=canvas, skipLumi=True)
        else:
            if hFakerateFactorCorrection[charge].GetNbinsY() == 1:
                hFRF_proj[charge] = hFakerateFactorCorrection[charge].ProjectionX(hFakerateFactorCorrection[charge].GetName()+"_projX", 1, 1, "e")
                axisName1D = xAxisName
            else:
                hFRF_proj[charge] = hFakerateFactorCorrection[charge].ProjectionY(hFakerateFactorCorrection[charge].GetName()+"_projY", 1, 1, "e")
                axisName1D = yAxisName
            drawTH1(hFRF_proj[charge], axisName1D, f"QCD template correction::{corrRange}", hFRF_proj[charge].GetName(), outfolder, drawStatBox=False,
                    plotTitleLatex=f"Charge {charge}")

        # plot 1D eta and pt in the signal region with the correction on the fakes, all from the nominal histogram
        for xbin in axisVar.keys():
            # only eta and pt, a bit hardcoded ...
            if xbin > 1:
                continue
            plotProjection1Dfrom3D(rootHists_nominal, datasetsNoQCDFakes, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                                   projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1],
                                   plotName=f"{axisVar[xbin][0]}_passIso_jetInclusive_passMt_noFakes",
                                   correctionFakes=None)
            plotProjection1Dfrom3D(rootHists_nominal, datasetsNoQCD, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                                   projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1],
                                   plotName=f"{axisVar[xbin][0]}_passIso_jetInclusive_passMt_FakesNoCorr",
                                   correctionFakes=None)
            plotProjection1Dfrom3D(rootHists_nominal, datasetsNoQCD, outfolder_dataMC, canvas1Dshapes=canvas1Dshapes, chargeBin=chargeBin,
                                   projectAxisToKeep=xbin, xAxisName=axisVar[xbin][1],
                                   plotName=f"{axisVar[xbin][0]}_passIso_jetInclusive_passMt_correctFakes",
                                   correctionFakes=hFakerateFactorCorrection[charge])
            
        outFile = outfolder + "fakerateFactorMtBasedCorrection_vsEtaPt.root"
        fout = safeOpenFile(outFile, mode="RECREATE")
        hFakerateFactorCorrection[charge].Write()
        if hFRF_proj[charge]:
            hFRF_proj[charge].Write()
        print()
        print(f"Saving FRF correction vs eta-pt in file\n{outFile}\nfor charge {charge}")
        print()
        fout.Close()

        print("-"*30)
        print("### Bad fit status")
        print("### | status | number of bad bins|")
        for key in sorted(badFitsID.keys()):
            print(f"{key}  {badFitsID[key]}")
        print("-"*30)
        print("### Bad covariance matrix status")
        print("### | status | number of bad bins|")
        for key in sorted(badCovMatrixID.keys()):
            print(f"{key}  {badCovMatrixID[key]}")
        print("-"*30)
        

    if len(charges) == 2:
        outFile = mainOutputFolder + "/fakerateFactorMtBasedCorrection_vsEtaPt.root"
        fout = safeOpenFile(outFile, mode="RECREATE")
        for charge in charges:
            hFakerateFactorCorrection[charge].Write()
            hFRF_proj[charge].Write()
        print()
        print(f"Saving FRF correction vs eta-pt in file {outFile}")
        print()
        fout.Close()
