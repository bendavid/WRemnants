#!/usr/bin/env python3

# Latest commands (check input file name)
# python w-mass-13TeV/smoothLeptonScaleFactors.py /eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_08Oct2022_binnedInPtEta_mass60to120/allSFs.root /eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_08Oct2022_binnedInPtEta_mass60to120/smoothLeptonScaleFactors/ -s iso

import os, re, array, math
import argparse
from copy import *

import numpy as np
import tensorflow as tf
import hist
import narf

from functools import partial

import utilitiesCMG
utilities = utilitiesCMG.util()

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

from make2DEffStatVariations import effStatVariations

import wremnants

# for a quick summary at the end
badFitsID_data = {}
badFitsID_mc = {}
badFitsID_sf = {}
badCovMatrixID_data = {}
badCovMatrixID_mc = {}
badCovMatrixID_sf = {}

# need to define functions as global below to avoid them being deleted out of fitTurnOnTF
pol3_tf_scaled = None
pol4_tf_scaled = None
erfPol2_tf_scaled = None
#############################################
# some functions for tensorflow

# for root to be consistent with functions passed to tensorflow
def pol3_root(xvals, parms, xLowVal = 0.0, xFitRange = 1.0):
    xscaled = (xvals[0] - xLowVal) / xFitRange
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2 + parms[3]*xscaled**3

def pol4_root(xvals, parms, xLowVal = 0.0, xFitRange = 1.0):
    xscaled = (xvals[0] - xLowVal) / xFitRange 
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2 + parms[3]*xscaled**3 + parms[4]*xscaled**4

#by hand version
def pol3_tf(xvals, parms):
    xscaled = (xvals[0] - xvals[0][0])/(xvals[0][-1] - xvals[0][0])
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2 + parms[3]*xscaled**3

def pol4_tf(xvals, parms):
    xscaled = (xvals[0] - xvals[0][0])/(xvals[0][-1] - xvals[0][0])
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2 + parms[3]*xscaled**3 + parms[4]*xscaled**4

def pol5_tf(xvals, parms):
    xscaled = (xvals[0] - xvals[0][0])/(xvals[0][-1] - xvals[0][0])
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2 + parms[3]*xscaled**3 + parms[4]*xscaled**4 + parms[5]*xscaled**5

#using tensforflow functions
def pol3_tf_v2(xvals, parms):
    coeffs = tf.unstack(tf.reverse(parms, axis=[0]))
    return tf.math.polyval(coeffs, xvals[0])

def erf_tf(xvals, parms):
    return parms[0] * (1.0 + tf.math.erf( (xvals[0] - parms[1]) / parms[2] ))

def erfPol2_tf(xvals, parms, xLowVal = 0.0, xFitRange = 1.0):
    xscaled = (xvals[0] - xLowVal) / xFitRange
    return  (parms[0] + parms[1] * xscaled + parms[2] * xscaled**2) * (1.0 + tf.math.erf( (xvals[0] - parms[3]) / parms[4] ))

def erfRatio_tf(xvals, parms):
    return parms[0] * (1.0 + tf.math.erf( (xvals[0] - parms[1]) / parms[2] )) / (1.0 + tf.math.erf( (xvals[0] - parms[3]) / parms[4] ))
                       
#############################################

def getReducedChi2andLabel(func):
    if func.GetNDF():
        reducedChi2 = func.GetChisquare() / func.GetNDF()
        lineChi2 = "#chi^{{2}} = {chi2:.2g} / {ndf}".format(chi2=func.GetChisquare(),ndf=int(func.GetNDF()))
    else:
        reducedChi2 = 0
        lineChi2 = "BAD! #chi^{{2}} = {chi2:.2g} / {ndf}".format(chi2=func.GetChisquare(),ndf=int(func.GetNDF()))

    return float(reducedChi2),lineChi2


def copyHisto(h1, h2, copyError=False):

    if h1.GetDimension() != h2.GetDimension():
        print(f"Error in copyHisto(): histograms have different dimensions. Dim(h1)={h1.GetDimension()}  Dim(h2)={h2.GetDimension()}. Exit")
        quit()

    if h1.GetDimension() == 1:
        for ix in range(h2.GetNbinsX()+2):
            if copyError:
                h1.SetBinContent(ix, h2.GetBinContent(ix, iy))
                h1.SetBinError(ix, h2.GetBinError(ix, iy))
            else:
                h1.SetBinContent(ix, h2.GetBinError(ix, iy))
                h1.SetBinError(ix, 0.0)
    elif h1.GetDimension() == 2:
        for ix in range(h2.GetNbinsX()+2):
            for iy in range(h2.GetNbinsY()+2):
                if copyError:
                    h1.SetBinContent(ix,iy, h2.GetBinError(ix, iy))
                    h1.SetBinError(ix,iy, 0.0)
                else:
                    h1.SetBinContent(ix,iy, h2.GetBinContent(ix, iy))
                    h1.SetBinError(ix,iy, h2.GetBinError(ix, iy))

    else:
        print("Error in copyHisto(): function not implemented for dimension > 2. Exit")
        quit()        

def make1Dhist(namePrefix, h2D, ptbins, step):
    hpt = {}
    for x in range(1, h2D.GetNbinsX()+1):
        binID = x-1
        hpt[binID] = ROOT.TH1D(f"{namePrefix}_{binID}",
                               "%s: %.4g < #eta < %.4g" % (step, h2D.GetXaxis().GetBinLowEdge(x), h2D.GetXaxis().GetBinLowEdge(x+1)),
                               len(ptbins)-1, array('d',ptbins)
        )
        for y in range(1, h2D.GetNbinsY()+1):
            hpt[binID].SetBinContent(y, h2D.GetBinContent(x,y))         
            hpt[binID].SetBinError(y, h2D.GetBinError(x,y))
    return hpt

def smoothSpline(hist, key, outname, channel="mu", drawFit=True):

    outdir = "{out}spline/".format(out=outname)
    createPlotDirAndCopyPhp(outdir)

    canvas = ROOT.TCanvas("canvas_%s" % key,"",700,700)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.06)
    canvas.cd()                           

    setTDRStyle()

    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)

    hist.GetXaxis().SetTitle("%s p_{T} (GeV)" % ("electron" if channel == "el" else "muon"))
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetTitle("Correction")
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.04)
    miny, maxy = getMinMaxHisto(hist, sumError=True)
    offset = 0.1 * (maxy - miny)
    miny -= offset
    maxy += offset
    hist.SetStats(0)
    hist.Draw("EP")

    # TSpline
    xval = []
    yval = []
    for i in range(1,hist.GetNbinsX()+1):
        xval.append(hist.GetBinCenter(i))
        yval.append(hist.GetBinContent(i))
    #spl = ROOT.TSpline3("spline3",array('d',xval),array('d',yval),len(xval),"b1e1")
    spl = ROOT.TSpline3(hist,"b1e1")
    spl.SetLineWidth(2)
    spl.SetLineColor(ROOT.kRed+1)
    spl.Draw("pclsame")

    leg = ROOT.TLegend(0.5, 0.2, 0.9, 0.3 if isIso else 0.4)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(spl, "spline", 'LF')

    for ext in ["pdf","png"]:
        canvas.SaveAs("{out}CorrectionFactorVsPt_{ch}_eta{b}.{ext}".format(out=outdir,ch=channel,b=key,ext=ext))            

    return spl

def smoothSomeFile(fname,hname,outname,outfilename,channel,widthPt):

    tfile = ROOT.TFile.Open(fname)        
    hist = tfile.Get(hname)
    if (hist == 0):
        print(f"Error: could not retrieve hist from input file {fname}. Exit")
        quit()
    hist.SetDirectory(0)
    tfile.Close()

    etabins = hist.GetXaxis().GetXbins()
    ptbins  = hist.GetYaxis().GetXbins()
    hist.SetTitle("")
    histSmooth = ROOT.TH2D("histSmooth","",
                           len(etabins)-1,array('d',etabins),
                           int(math.ceil(hist.GetYaxis().GetBinLowEdge(1+hist.GetNbinsY()) - hist.GetYaxis().GetBinLowEdge(1))/widthPt), # bins of 0.2 GeV
                           hist.GetYaxis().GetBinLowEdge(1),hist.GetYaxis().GetBinLowEdge(1+hist.GetNbinsY())        
                           )

    histpt = {}
    for x in range(1,hist.GetNbinsX()+1):
        bin = x-1
        histpt[bin] = ROOT.TH1D("histpt_{b}".format(b=str(bin)),
                                "%.4g <= #eta < %.4g" % (hist.GetXaxis().GetBinLowEdge(x), hist.GetXaxis().GetBinLowEdge(x+1)),
                                len(ptbins)-1,array('d',ptbins)
                                )
        for y in range(1,hist.GetNbinsY()+1):
            histpt[bin].SetBinContent(y,hist.GetBinContent(x,y))         
            histpt[bin].SetBinError(y,hist.GetBinError(x,y))
            
    bestFit_hist = {}
    for key in histpt:

        spline = smoothSpline(histpt[key],key,outname,channel=channel)
        for ipt in range(1,histSmooth.GetNbinsY()+1):
            ptval = histSmooth.GetYaxis().GetBinCenter(ipt)
            histSmooth.SetBinContent(key+1,ipt, spline.Eval(ptval))

    
    setTDRStyle()
    drawCorrelationPlot(hist,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Original value::-0.2,1.2",
                        "original_{hn}".format(hn=hname),"ForceTitle",outname,palette=55)
    drawCorrelationPlot(histSmooth,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Smoothed value::-0.2,1.2",
                        "smooth_{hn}".format(hn=hname),"ForceTitle",outname,palette=55)

    tfile = ROOT.TFile.Open(outname+outfilename,'recreate')
    histSmooth.Write()
    hist.Write("histOriginal")
    tfile.Close()
    print()
    print(f"Created file {outname+outfilename}")
    print()
    return 0
    
############

def fitTurnOnTF(hist, key, outname, mc, channel="el", hist_chosenFunc=0, drawFit=True, 
                step=None,
                fitRange=None,
                hist_reducedChi2=0,
                hist_FuncParam_vs_eta=0,
                hist_FuncCovMatrix_vs_eta=0,  # TH3, eta on x and cov matrix in yz
                charge = "both",
                etabins = []
):

    doingSF = True if mc == "SF" else False
        
    chargeText = ""
    if charge == "plus": chargeText = "positive"
    if charge == "minus": chargeText = "negative"
    
    originalMaxPt = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())

    outdir = "{out}{mc}/".format(out=outname,mc=mc)
    createPlotDirAndCopyPhp(outdir)

    canvas = ROOT.TCanvas(f"canvas_{mc}_{key}","",700,700)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.06)
    canvas.cd()                           

    setTDRStyle()

    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)

    hist.GetXaxis().SetTitle(f"{chargeText} muon p_{{T}} (GeV)")
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)

    if fitRange != None:
        if fitRange[0] >= 0 and fitRange[1] >= 0:
            hist.GetXaxis().SetRangeUser(fitRange[0], fitRange[1])
        elif fitRange[0] >= 0:
            hist.GetXaxis().SetRangeUser(fitRange[0], hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX()))
        elif fitRange[1] >= 0:
            hist.GetXaxis().SetRangeUser(hist.GetXaxis().GetBinLowEdge(1),fitRange[1])

    if mc == "SF":
        hist.GetYaxis().SetTitle("Data/MC scale factor")
    else:
        hist.GetYaxis().SetTitle("{mc} efficiency".format(mc=mc))
    hist.GetYaxis().SetTitleOffset(1.4)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.04)
    miny,maxy = getMinMaxHisto(hist, sumError=True)
    offset = 0.1 * (maxy - miny)
    upOffset = offset * (2.5 if doingSF else 3.4)
    miny -= offset
    maxy += upOffset
    hist.GetYaxis().SetRangeUser(miny, maxy)
    hist.SetStats(0)
    hist.Draw("EP")

    maxFitRange = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())
    minFitRange = hist.GetXaxis().GetBinLowEdge(1)
    if fitRange != None:
        if fitRange[1] > 0:
            maxFitRange = fitRange[1]
        if fitRange[0] > 0:
            minFitRange = fitRange[0]
    ##
    ## TODO
    ## Unlike root, tensorflow fits use the bin centers to run the actual fit, and there is no concept of fit range
    ## so, if one wants to fit a subset of the histogram range one needs to pass the slice of the boost histogram
    ## currently this is not implemented
        
            
    ###################
    # fits
    ####################
    # Erf defined here: https://root.cern.ch/doc/v608/namespaceTMath.html#a44e82992cba4684280c72679f0f39adc
    # Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x 
    boost_hist = narf.root_to_hist(hist)

    ###############################################################
    fitFunction = None
    fitres_TF = None
    defaultFunc = ""
    badFitsID = None
    badCovMatrixID = None
    
    xFitRange = maxFitRange - minFitRange

    if doingSF:
        #binCenter1 = hist.GetXaxis().GetBinCenter(1)
        #xHistRange = hist.GetXaxis().GetBinCenter(hist.GetNbinsX()) - binCenter1 # use bin centers for consistency with xvals in tensorflow fits
        global pol3_tf_scaled
        if pol3_tf_scaled == None:
            pol3_tf_scaled = partial(pol3_root, xLowVal=minFitRange, xFitRange=xFitRange)
        params = np.array([1.0, 0.0, 0.0, 0.0])
        res_tf1_pol3 = narf.fit_hist(boost_hist, pol3_tf_scaled, params)
        tf1_pol3 = ROOT.TF1("tf1_pol3", pol3_tf_scaled, minFitRange, maxFitRange, len(params))
        tf1_pol3.SetParameters( np.array( res_tf1_pol3["x"], dtype=np.dtype('d') ) )
        tf1_pol3.SetLineWidth(3)
        tf1_pol3.SetLineColor(ROOT.kRed+2)

        # tf1_pol3_test = ROOT.TF1("tf1_pol3_test", pol3_tf_scaled, minFitRange, maxFitRange, len(params))
        # tf1_pol3_test.SetParameters( np.array( res_tf1_pol3["x"], dtype=np.dtype('d') ) )
        # tf1_pol3_test.SetLineStyle(ROOT.kDashed)
        # tf1_pol3_test.SetLineWidth(5)
        # tf1_pol3_test.SetLineColor(ROOT.kBlue)
        # fitopt = "FMBRQS+" # add FM if using Minuit   
        # hist.Fit(tf1_pol3_test, fitopt)
        
        fitres_TF = {"pol3_tf" : res_tf1_pol3}
        fitFunction = {
            "pol3_tf" : {
                "func" : tf1_pol3,
                "leg"  : "Pol3",
            }
        }
        defaultFunc = "pol3_tf"
        badFitsID = badFitsID_sf
        badCovMatrixID = badCovMatrixID_sf
    else:
        # global erfPol2_tf_scaled
        # erfPol2_tf_scaled = partial(erfPol2_tf, xLowVal=minFitRange, xFitRange=xFitRange)
        # params = np.array([1.0, 0.0, 0.0, 35.0, 3.0])
        # res_tf1_erfPol2 = narf.fit_hist(boost_hist, erfPol2_tf_scaled, params)
        # tf1_erfPol2 = ROOT.TF1("tf1_erfPol2", erfPol2_tf_scaled, minFitRange, maxFitRange, len(params))
        # tf1_erfPol2.SetParameters( np.array( res_tf1_erfPol2["x"], dtype=np.dtype('d') ) )
        # tf1_erfPol2.SetLineWidth(2)
        # tf1_erfPol2.SetLineStyle(ROOT.kDashed)
        # tf1_erfPol2.SetLineColor(ROOT.kGreen+2)

        tf1_erf = ROOT.TF1("tf1_erf","[0] * (1.0 + TMath::Erf((x-[1])/[2]))", minFitRange, maxFitRange)
        res_tf1_erf = narf.fit_hist(boost_hist, erf_tf, np.array([1.0, 35.0, 3.0]))
        tf1_erf.SetParameters( np.array( res_tf1_erf["x"], dtype=np.dtype('d') ) )
        tf1_erf.SetLineWidth(2)
        tf1_erf.SetLineStyle(ROOT.kDashed)
        tf1_erf.SetLineColor(ROOT.kBlue)

        global pol4_tf_scaled
        if pol4_tf_scaled == None:
            pol4_tf_scaled = partial(pol4_root, xLowVal=minFitRange, xFitRange=xFitRange)        
        params = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
        res_tf1_pol4 = narf.fit_hist(boost_hist, pol4_tf_scaled, params)
        tf1_pol4 = ROOT.TF1("tf1_pol4", pol4_tf_scaled, minFitRange, maxFitRange, len(params))
        tf1_pol4.SetParameters( np.array( res_tf1_pol4["x"], dtype=np.dtype('d') ) )
        tf1_pol4.SetLineWidth(3)
        tf1_pol4.SetLineColor(ROOT.kRed+2)
        #
        fitres_TF = {"erf" : res_tf1_erf,
                     "pol4_tf" : res_tf1_pol4,
                     #"erfPol2_tf" : res_tf1_erfPol2
        }
        fitFunction = {
            "pol4_tf" : {
                "func" : tf1_pol4,
                "leg"  : "Pol4",
            },
            "erf" : {
                "func" : tf1_erf,
                "leg"  : "Erf",
            },
            # "erfPol2_tf" : {
            #     "func" : tf1_erfPol2,
            #     "leg"  : "ErfPol2",
            # },
        }
        defaultFunc = "pol4_tf"
        if mc == "MC":
            badFitsID = badFitsID_mc
            badCovMatrixID = badCovMatrixID_mc
        else:
            badFitsID = badFitsID_data
            badCovMatrixID = badCovMatrixID_data

    ###############################################################
            
    for fr in fitres_TF.keys():
        status = fitres_TF[fr]["status"]
        covstatus = fitres_TF[fr]["covstatus"]
        if status:
            print(f"Function {fr} had status {status}")
            if key not in badFitsID.keys():
                badFitsID[key] = {}
            badFitsID[key][fr] = status
        if covstatus:
            print(f"Function {fr} had covstatus {covstatus}")
            if key not in badCovMatrixID.keys():
                badCovMatrixID[key] = {}
            badCovMatrixID[key][fr] = covstatus

    npar = fitFunction[defaultFunc]["func"].GetNpar()
    
    if hist_FuncCovMatrix_vs_eta:
        for i in range(npar):
            for j in range(npar):
                hist_FuncCovMatrix_vs_eta.SetBinContent(key+1, i+1, j+1, fitres_TF[defaultFunc]["cov"][i][j])

    if hist_FuncParam_vs_eta:
        # key is the eta bin number, but starts from 0, so add 1
        for ip in range(npar):
            hist_FuncParam_vs_eta.SetBinContent(key+1, ip+1, fitres_TF[defaultFunc]["x"][ip])
            hist_FuncParam_vs_eta.SetBinError(  key+1, ip+1, math.sqrt(fitres_TF[defaultFunc]["cov"][ip][ip]))

    hband = ROOT.TH1D("hband", "", int((maxFitRange-minFitRange)/0.2), minFitRange, maxFitRange)
    hband.SetStats(0)
    hband.SetFillColor(ROOT.kGray)
    #hband.SetFillStyle(3001)
    for ib in range(1, hband.GetNbinsX()+1):
        pt = hband.GetBinCenter(ib)
        hband.SetBinContent(ib, fitFunction[defaultFunc]["func"].Eval(pt))
    # eigen decomposition to plot alternate curves
    #
    # can pass full histogram, eta bin to use is set below
    systCalc = ROOT.wrem.EtaPtCorrelatedEfficiency(hist_FuncCovMatrix_vs_eta, hist_FuncParam_vs_eta, minFitRange, maxFitRange)
    systCalc.setSmoothingFunction(defaultFunc)    
    vecUp = ROOT.std.vector["double"]()
    vecUp = systCalc.DoEffSyst(key+1)
    #systCalc.setEigenShift(-1.0); # shift down
    #vecDown = systCalc.DoEffSyst(key+1)

    ## some functions cannot be cloned, although in python it might work
    ##
    #tf1_func_alt = copy.deepcopy(fitFunction[defaultFunc]["func"].Clone("tf1_func_alt"))
    #tf1_func_alt = copy.deepcopy(fitFunction[defaultFunc]["func"])
    tf1_func_alt = ROOT.TF1()
    tf1_func_alt.SetName("tf1_func_alt")
    fitFunction[defaultFunc]["func"].Copy(tf1_func_alt)
    tf1_func_alt.SetLineWidth(1)
    for ib in range(1, hband.GetNbinsX()+1):
        pt = hband.GetBinCenter(ib)
        err2 = 0.0
        for ivar in range(npar):
            startIndex = npar*ivar
            # set parameters for a given hessian
            tf1_func_alt.SetParameters(np.array([vecUp[i] for i in range(startIndex, startIndex+npar)], dtype=np.dtype('d')))
            diff = tf1_func_alt.Eval(pt) - hband.GetBinContent(ib)
            err2 += diff * diff 
        hband.SetBinError(ib, math.sqrt(err2))

    hband.Draw("E4SAME")
    # redraw to have them on top
    hist.Draw("EPSAME")
    for f in fitFunction.keys():
        fitFunction[f]["func"].Draw("LSAME")
        
    colors_alt = [ROOT.kCyan+1, ROOT.kGreen+2, ROOT.kBlue, ROOT.kSpring+9, ROOT.kOrange+2, ROOT.kPink, ROOT.kViolet, ROOT.kMagenta]
    for ivar in range(npar):
        startIndex = npar*ivar
        # set parameters for a given hessian
        tf1_func_alt.SetParameters(np.array([vecUp[i] for i in range(startIndex, startIndex+npar)], dtype=np.dtype('d')))
        #tf1_func_alt.SetLineStyle(ROOT.kDotted)
        tf1_func_alt.SetLineColor(colors_alt[ivar])
        tf1_func_alt.DrawCopy("LSAME") # can also not draw these lines, too busy plot otherwise
        #tf1_func_alt.SetParameters(np.array([vecDown[i] for i in range(startIndex, startIndex+npar)], dtype=np.dtype('d')))
        #tf1_func_alt.DrawCopy("LSAME") # can also not draw these lines, too busy plot otherwise
        
    #######
    nFits = len(fitFunction.keys())

    upLeg = 0.9
    downLeg = max(0.6, upLeg - 0.06 * (nFits + 1)) # +1 to include the eror band
    leftLeg = 0.2
    rightLeg = 0.9
        
    leg = ROOT.TLegend(leftLeg, downLeg, rightLeg, upLeg)
    leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(0,0.6)
    
    # for chosen function
    reducedChi2 = 0.0
    for f in fitFunction.keys():
        chi2 = fitres_TF[f]["loss_val"]
        ndof = int(hist.GetNbinsX() - fitFunction[f]["func"].GetNpar())
        legEntry = fitFunction[f]["leg"] + f"   #chi^{{2}} = {round(chi2,1)} / {ndof}"
        chi2_nsigma = abs(chi2 - ndof) / math.sqrt(2.0 * ndof) 
        if chi2_nsigma > 3.0:
            legEntry += " BAD!"
        leg.AddEntry(fitFunction[f]["func"], legEntry, 'L')  
        if f == defaultFunc:
            reducedChi2 = chi2/ndof
    leg.AddEntry(hband, "Model uncertainty", 'F')
    leg.Draw('same')
    canvas.RedrawAxis("sameaxis")

    setTDRStyle()
    ROOT.gStyle.SetOptTitle(1)  # use histogram title with binning as canvas title

    if hist_reducedChi2:
        hist_reducedChi2.Fill(reducedChi2)
    if hist_chosenFunc:
        hist_chosenFunc.Fill(defaultFunc, 1)

    # lat = ROOT.TLatex()
    # line = ""
    # lineChi2 = ""
    # lat.SetNDC();
    # lat.SetTextSize(0.045);
    # lat.SetTextFont(42);
    # lat.SetTextColor(1);        
    # xmin = 0.20 
    # yhi = 0.85
    # lat.DrawLatex(xmin, yhi, line);
    # lat.DrawLatex(xmin, yhi-0.05, lineChi2);
    
    tmpch = ""
    if charge != "both":
        tmpch = "_" + charge
    for ext in ["pdf","png"]:
        if mc == "SF":
            canvas.SaveAs("{out}sf_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,ch=channel,b=key,charge=tmpch,ext=ext))            
        else:
            canvas.SaveAs("{out}eff{mc}_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,charge=tmpch,ext=ext))                            
    
    return fitFunction[defaultFunc]["func"]

######

minmaxSF = {"trigger"      : "0.65,1.15",
            "idip"         : "0.95,1.01",
            "iso"          : "0.975,1.025",
            "antiiso"      : "0.6,1.25",
            "isonotrig"    : "0.97,1.03",
            "antiisonotrig": "0.6,1.25",
            "tracking"     : "0.98,1.01",
            "reco"         : "0.94,1.02",
}

    
        
if __name__ == "__main__":
            
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile',  type=str, nargs=1, help='input root file with TH2')
    parser.add_argument('outdir', type=str, nargs=1, help='output directory to save things')
    parser.add_argument('-c','--charge', default='both', choices=['both', 'plus', 'minus'], type=str, help='Plus or minus if the efficiencies were derived separately for each charge. If both, assumes no charge splitting in the inputs')
    parser.add_argument('-e','--era',  dest='era',     default='GtoH', choices=['GtoH'], type=str, help='Efficiency era')
    parser.add_argument('-s','--step', dest='step', default='', choices=list(minmaxSF.keys()), required=True, help='Working point to smooth')
    parser.add_argument('-r','--pt-fit-range', dest='ptFitRange', type=float, nargs=2, default=[-1, -1], help='Pt range fo the fit: pass two values for min and max. If one of them (or both) is negative, the corresponding histogram range is used')
    parser.add_argument('-w','--width-pt',     dest='widthPt',default='0.2', type=float, help='Pt bin width for the smoothed histogram')
    parser.add_argument(     '--set-max-pt-histo',     dest='setMaxPtHisto', default='-1.0', type=float, help='Set upper pt for output histograms. If negative use default max from input histograms')
    #parser.add_argument(     '--save-TF1', dest='saveTF1',action="store_true", default=False, help='Save TF1 as well, not just TH2 with many bins')
    parser.add_argument(    '--input-hist-names', dest='inputHistNames', default='', type=str, help='Pass comma separated list of 3  names, for eff(data),eff(MC),SF, to be used instead of the default names')
    parser.add_argument(     '--palette'  , dest='palette',      default=87, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--skip-eff', dest='skipEff', action="store_true", default=False, help='Skip efficiencies and do only SF (to save time and if one only wants to smooth SF directly)')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    channel = "mu"
    charge = "" if args.charge == "both" else "positive" if args.charge == "plus" else "negative"
    lepton = f"{charge} muon"

    outname = args.outdir[0]
    addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
    outname += f"{args.era}/{channel}_{args.step}_{args.charge}/"
    createPlotDirAndCopyPhp(outname)
    
    outfilename = f"smoothedSFandEffi_{args.step}_{args.era}_{args.charge}.root"

    #########################################    
    #########################################    

    datahistname = f"effData_{args.step}_{args.era}_{args.charge}"
    mchistname   = f"effMC_{args.step}_{args.era}_{args.charge}"
    sfhistname   = f"SF2D_nominal_{args.step}_{args.era}_{args.charge}"

    if len(args.inputHistNames):
        datahistname,mchistname,sfhistname = args.inputHistNames.split(",")
    
    tfile = safeOpenFile(args.inputfile[0])
    hsf =   safeGetObject(tfile, sfhistname)
    if args.skipEff:
        hdata = copy.deepcopy(hsf.Clone("data_dummy"))
        hdata.Reset("ICESM")
        hmc = copy.deepcopy(hdata.Clone("mc_dummy"))
    else:
        hdata = safeGetObject(tfile, datahistname)
        hmc =   safeGetObject(tfile, mchistname)
    tfile.Close()
        
    etabins = [round(hdata.GetXaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hdata.GetNbinsX())]
    ptbins =  [round(hdata.GetYaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hdata.GetNbinsY())]
        
    dummybins = [-0.5, 0.5, 1.5, 2.5]
    title = "parameters: [0]*TMath::Erf((x-[1])/[2])"
    hist_FuncParam_vs_eta_data = ROOT.TH2D("hist_FuncParam_vs_eta_data",
                                          title,
                                          len(etabins)-1, array('d',etabins),
                                          len(dummybins)-1, array('d',dummybins))    
    hist_FuncParam_vs_eta_mc   = ROOT.TH2D("hist_FuncParam_vs_eta_mc",
                                          title,
                                          len(etabins)-1, array('d',etabins),
                                          len(dummybins)-1, array('d',dummybins))

    hist_FuncCovMatrix_vs_eta_data = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_data",
                                              "Covariance matrix: eta on X",
                                              len(etabins)-1, array('d',etabins),
                                              len(dummybins)-1, array('d',dummybins),
                                              len(dummybins)-1, array('d',dummybins))
    hist_FuncCovMatrix_vs_eta_mc = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_mc",
                                            "Covariance matrix: eta on X",
                                            len(etabins)-1, array('d',etabins),
                                            len(dummybins)-1, array('d',dummybins),
                                            len(dummybins)-1, array('d',dummybins))

    dummybins_sf = [-0.5, 0.5, 1.5, 2.5, 3.5]
    title_sf = "parameters: cheb3"
    hist_FuncParam_vs_eta_sf = ROOT.TH2D("hist_FuncParam_vs_eta_sf",
                                         title_sf,
                                         len(etabins)-1, array('d',etabins),
                                         len(dummybins_sf)-1, array('d',dummybins_sf))    
    hist_FuncCovMatrix_vs_eta_sf = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_sf",
                                             "Covariance matrix: eta on X",
                                             len(etabins)-1, array('d',etabins),
                                             len(dummybins_sf)-1, array('d',dummybins_sf),
                                             len(dummybins_sf)-1, array('d',dummybins_sf))
    
    # utility histogram to show what function will be used based on some quality criterium (usually we choose Erf, but can be pol3 or something else)
    hist_chosenFunc = ROOT.TH1D("chosenFitFunc", "Best fit function for each #eta bin for data or MC", 3, 0, 3)
    hist_chosenFunc.GetXaxis().SetBinLabel(1, "tf1_erf")
    hist_chosenFunc.GetXaxis().SetBinLabel(2, "tf1_pol2")
    hist_chosenFunc.GetXaxis().SetBinLabel(3, "tf1_pol3")
    # 
    hist_chosenFunc_SF = ROOT.TH1D("chosenFitFunc_SF", "Best fit function for each #eta bin for SF", 3, 0, 3)
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(1, "tf1_cheb3")
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(2, "tf1_pol2")
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(3, "tf1_pol3")

    hist_reducedChi2_data = ROOT.TH1D("reducedChi2_data", "Reduced #chi^{2}", 50, 0, 5) # will not have many entries (~100 depending on how many eta bins)
    hist_reducedChi2_data.StatOverflows() # use underflow and overflow to compute mean and RMS
    hist_reducedChi2_MC = ROOT.TH1D("reducedChi2_MC", "Reduced #chi^{2}", 50, 0, 5) # will not have many entries (~100 depending on how many eta bins)
    hist_reducedChi2_MC.StatOverflows() # use underflow and overflow to compute mean and RMS
    hist_reducedChi2_sf = ROOT.TH1D("reducedChi2_sf", "Reduced #chi^{2}", 50, 0, 5) # will not have many entries (~100 depending on how many eta bins)
    hist_reducedChi2_sf.StatOverflows() # use underflow and overflow to compute mean and RMS

    ######################
    # to make ratio
    ######################
    ratioData =  ROOT.TH2D("dataEfficiencyRatio","Original/smooth Data efficiency ratio",
                           len(etabins)-1, array('d',etabins),
                           len(ptbins)-1, array('d',ptbins)
                           )
    ratioMC =  ROOT.TH2D("mcEfficiencyRatio","Original/smooth MC efficiency ratio",
                           len(etabins)-1, array('d',etabins),
                           len(ptbins)-1, array('d',ptbins)
                           )
    ratioSF =  ROOT.TH2D("scaleFactorRatio","Original/smooth scale factor ratio",
                         len(etabins)-1, array('d',etabins),
                         len(ptbins)-1, array('d',ptbins)
                         )

    copyHisto(ratioData,hdata)
    copyHisto(ratioMC,hmc)
    copyHisto(ratioSF,hsf)

    ######################
    # to make ratio
    ######################
    pullData =  ROOT.TH2D("dataEfficiencyPull","Original/smooth Data efficiency pull",
                           len(etabins)-1, array('d',etabins),
                           len(ptbins)-1, array('d',ptbins)
                           )
    pullMC =  ROOT.TH2D("mcEfficiencyPull","Original/smooth MC efficiency pull",
                           len(etabins)-1, array('d',etabins),
                           len(ptbins)-1, array('d',ptbins)
                           )
    pullSF =  ROOT.TH2D("scaleFactorPull","Original/smooth scale factor pull",
                         len(etabins)-1, array('d',etabins),
                         len(ptbins)-1, array('d',ptbins)
                         )

    copyHisto(pullData, hdata)
    copyHisto(pullMC, hmc)
    copyHisto(pullSF, hsf)
    errData = copy.deepcopy(hdata.Clone("errData"))
    errMC   = copy.deepcopy(hmc.Clone("errMC"))
    errSF   = copy.deepcopy(hsf.Clone("errSF"))
    copyHisto(errData, hdata, copyError=True)
    copyHisto(errMC, hmc, copyError=True)
    copyHisto(errSF, hsf, copyError=True)
    
    #############
    # these will be used to check the smoothed efficiency
    ###############
    maxPtHistoData = hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY())
    maxPtHistoMC = hmc.GetYaxis().GetBinLowEdge(1+hmc.GetNbinsY())
    if args.setMaxPtHisto > 0.0:
        maxPtHistoData = args.setMaxPtHisto
        maxPtHistoMC = args.setMaxPtHisto

    nFinePtBins = int(math.ceil(maxPtHistoData - hdata.GetYaxis().GetBinLowEdge(1))/args.widthPt)
    minPtHistoData = hdata.GetYaxis().GetBinLowEdge(1)
    
    hdataSmoothCheck = ROOT.TH2D("hdataSmoothCheck","Data smoothed efficiency",
                                 len(etabins)-1, array('d',etabins),
                                 nFinePtBins, minPtHistoData, maxPtHistoData)
                                 
    hmcSmoothCheck = ROOT.TH2D("hmcSmoothCheck","MC smoothed efficiency",
                               len(etabins)-1, array('d',etabins),
                               nFinePtBins, minPtHistoData, maxPtHistoData)
                               
    hsfSmoothCheck = ROOT.TH2D("hsfSmoothCheck","Data/MC smoothed scale factor",
                               len(etabins)-1, array('d',etabins),
                               nFinePtBins, minPtHistoData, maxPtHistoData)
                               
    hdataSmoothCheck_origBinPt = ROOT.TH2D("hdataSmoothCheck_origBinPt","Data smoothed efficiency",
                                           len(etabins)-1, array('d',etabins),
                                           len(ptbins)-1, array('d',ptbins)
                                           )
    hmcSmoothCheck_origBinPt = ROOT.TH2D("hmcSmoothCheck_origBinPt","MC smoothed efficiency",
                                         len(etabins)-1, array('d',etabins),
                                         len(ptbins)-1, array('d',ptbins)
                                         )
    hsfSmoothCheck_origBinPt = ROOT.TH2D("hsfSmoothCheck_origBinPt","Data/MC smoothed scale factor",
                                         len(etabins)-1, array('d',etabins),
                                         len(ptbins)-1, array('d',ptbins)
                                         )

    # hmc and hdata have eta on X and pt on Y
    # we select slices at constant eta and fit along pt with some function

    label = args.step + (args.charge if args.charge != "both" else "")
    hmcpt = make1Dhist("hmcpt", hmc, ptbins, label)
    hdatapt = make1Dhist("hdatapt", hdata, ptbins, label)
    hsfpt = make1Dhist("hsfpt", hsf, ptbins, label)

    bestFit_MC = {}
    bestFit_Data = {}
    if not args.skipEff:
        ###########################
        # first MC
        ###########################
        for key in hmcpt:
            bestFitFunc = fitTurnOnTF(hmcpt[key], key, outname, "MC",channel=channel,
                                      hist_chosenFunc=hist_chosenFunc, 
                                      step=args.step,
                                      fitRange=args.ptFitRange,
                                      hist_reducedChi2=hist_reducedChi2_MC,
                                      hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_mc,
                                      hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_mc,
                                      charge=args.charge,
                                      etabins=etabins)
            bestFit_MC["smoothFunc_MC_ieta%d" % key] = bestFitFunc
            for ipt in range(1, hmcSmoothCheck.GetNbinsY()+1):
                ptval = hmcSmoothCheck.GetYaxis().GetBinCenter(ipt)
                hmcSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
            for ipt in range(1, hmcSmoothCheck_origBinPt.GetNbinsY()+1):
                ptval = hmcSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
                hmcSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

        ###########################
        # now data
        ###########################
        for key in hdatapt:

            bestFitFunc = fitTurnOnTF(hdatapt[key],key,outname, "Data",channel=channel,hist_chosenFunc=hist_chosenFunc, 
                                      step=args.step,
                                      fitRange=args.ptFitRange,
                                      hist_reducedChi2=hist_reducedChi2_data,
                                      hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_data,
                                      hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_data,
                                      charge=args.charge,
                                      etabins=etabins)
            bestFit_Data["smoothFunc_Data_ieta%d" % key] = bestFitFunc
            for ipt in range(1,hdataSmoothCheck.GetNbinsY()+1):
                ptval = hdataSmoothCheck.GetYaxis().GetBinCenter(ipt)
                #hdataSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
                hdataSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
            for ipt in range(1,hdataSmoothCheck_origBinPt.GetNbinsY()+1):
                ptval = hdataSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
                #hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
                hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

    # ###########################
    # # now SF
    # ###########################
    bestFit_SF = {}
    for key in hsfpt:
        
        bestFitFunc = fitTurnOnTF(hsfpt[key],key,outname, "SF",channel=channel,hist_chosenFunc=hist_chosenFunc_SF, 
                                  step=args.step,
                                  fitRange=args.ptFitRange, hist_reducedChi2=hist_reducedChi2_sf,
                                  hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_sf,
                                  hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_sf,
                                  charge=args.charge, etabins=etabins)
        bestFit_SF["smoothFunc_SF_ieta%d" % key] = bestFitFunc
        for ipt in range(1,hsfSmoothCheck.GetNbinsY()+1):
            ptval = hsfSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #hsfSmoothCheck.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hsfSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        for ipt in range(1,hsfSmoothCheck_origBinPt.GetNbinsY()+1):
            ptval = hsfSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
            #hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, a0 + a1 * ptval + a2 * ptval * ptval)
            hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))

    #################################
    # start to make plots
    #################################
    zaxisRange = ""
    zaxisRangeSF = "::" + minmaxSF[args.step]

    canvas = ROOT.TCanvas("canvas","",700,625)

    # plot original histograms
    drawCorrelationPlot(hmc,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency%s" % zaxisRange,
                        "inputEfficiency_MC","",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hdata,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency%s" % zaxisRange,
                        "inputEfficiency_Data","",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hsf,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor%s" % zaxisRangeSF,
                        "inputScaleFactor","",outname,palette=args.palette,passCanvas=canvas)
    # now the new ones

    # make a sanity check plot: fill eta-pt with smoothed efficiency
    drawCorrelationPlot(hmcSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC smoothed efficiency%s" % zaxisRange,
                        "smoothEfficiency_MC","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hdataSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data smoothed efficiency%s" % zaxisRange,
                        "smoothEfficiency_Data","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hsfSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC smoothed scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactorDirectly","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    # scale factor: data/MC
    scaleFactor = ROOT.TH2D("scaleFactor","Scale factor from smoothed efficiencies",
                            len(etabins)-1, array('d',etabins),
                            int(math.ceil((maxPtHistoData - hdata.GetYaxis().GetBinLowEdge(1))/args.widthPt)),
                            hdata.GetYaxis().GetBinLowEdge(1),maxPtHistoData)                            

    copyHisto(scaleFactor, hdataSmoothCheck)
    scaleFactor.Divide(hmcSmoothCheck)
    scaleFactor.SetMinimum(scaleFactor.GetBinContent(scaleFactor.GetMinimumBin()))
    scaleFactor.SetMaximum(scaleFactor.GetBinContent(scaleFactor.GetMaximumBin()))
    drawCorrelationPlot(scaleFactor,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactor","ForceTitle",outname,palette=args.palette,passCanvas=canvas)

    #################################
    # plot also with oiginal binning
    ################################
    
    # divide before drawing the denominator, whose axis settings are modified by drawCorrelationPlot and seems to affect the ratio as well if divided afterwards
    ratioData.Divide(hdataSmoothCheck_origBinPt)
    ratioMC.Divide(hmcSmoothCheck_origBinPt)
    
    drawCorrelationPlot(hmcSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC smoothed efficiency%s" % zaxisRange,
                        "smoothEfficiency_MC_origBinPt","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hdataSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data smoothed efficiency%s" % zaxisRange,
                        "smoothEfficiency_Data_origBinPt","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hsfSmoothCheck_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC smoothed scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactorDirectly_origBinPt","ForceTitle",outname,palette=args.palette,passCanvas=canvas)

    # scale factor: data/MC
    scaleFactor_origBinPt = ROOT.TH2D("scaleFactor_origBinPt","Scale factor",
                                      len(etabins)-1, array('d',etabins),
                                      len(ptbins)-1, array('d',ptbins)
                                      )
    copyHisto(scaleFactor_origBinPt, hdataSmoothCheck_origBinPt)
    scaleFactor_origBinPt.Divide(hmcSmoothCheck_origBinPt)
    scaleFactor_origBinPt.SetMinimum(scaleFactor_origBinPt.GetBinContent(scaleFactor_origBinPt.GetMinimumBin()))
    scaleFactor_origBinPt.SetMaximum(scaleFactor_origBinPt.GetBinContent(scaleFactor_origBinPt.GetMaximumBin()))

    # to make ratio, divide before passing to function, to avoid changes in the histogram
    ratioSF.Divide(scaleFactor_origBinPt)

    #scaleFactor_origBinPt.GetZaxis().SetTitle("Data/MC scale factor")
    drawCorrelationPlot(scaleFactor_origBinPt,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactor_origBinPt","ForceTitle",outname,palette=args.palette,passCanvas=canvas)


    ######################
    # finally SF(smooth)/SF(original)
    ######################
    drawCorrelationPlot(ratioData,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency ratio (original/smooth)::0.98,1.02",
                        "dataEfficiencyRatio","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(ratioMC,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency ratio (original/smooth)::0.98,1.02",
                        "mcEfficiencyRatio","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(ratioSF,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"scale factor ratio (original/smooth)::0.98,1.02",
                        "scaleFactorRatio","ForceTitle",outname,palette=args.palette,passCanvas=canvas)

    pullData.Add(hdataSmoothCheck_origBinPt, -1.0)
    pullData.Divide(errData)
    pullMC.Add(hmcSmoothCheck_origBinPt, -1.0)
    pullMC.Divide(errMC)
    pullSF.Add(hsfSmoothCheck_origBinPt, -1.0)
    pullSF.Divide(errSF)
    drawCorrelationPlot(pullData,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data eff. pull (original-smooth)/err::-5.0,5.0",
                        "dataEfficiencyPull","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(pullMC,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC eff. pull (original-smooth)/err::-5.0,5.0",
                        "mcEfficiencyPull","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(pullSF,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"scale factor pull (original-smooth)/err::-5.0,5.0",
                        "scaleFactorPull","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    # ######################
    # # See the difference between smoothing Data and MC efficiency and taking the ratio or smoothing directly the efficiency ratio
    # ######################    
    ratioSF_smoothNumDen_smoothRatio = ROOT.TH2D("ratioSF_smoothNumDen_smoothRatio","SF ratio: smooth eff or ratio directly",
                                                 len(etabins)-1, array('d',etabins),
                                                 len(ptbins)-1, array('d',ptbins)
    )

    copyHisto(ratioSF_smoothNumDen_smoothRatio,scaleFactor_origBinPt)
    ratioSF_smoothNumDen_smoothRatio.Divide(hsfSmoothCheck_origBinPt)
    drawCorrelationPlot(ratioSF_smoothNumDen_smoothRatio,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),
                        "SF ratio: smooth eff or ratio directly::0.98,1.02",
                        "ratioSF_smoothNumDen_smoothRatio","ForceTitle",outname,palette=args.palette,passCanvas=canvas)

    ##############
    # Erf[x] parameter and error for data and MC efficiency
    drawCorrelationPlot(hist_FuncParam_vs_eta_data,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
                        "parameter value",
                        "hist_FuncParam_vs_eta_data","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)
    
    drawCorrelationPlot(hist_FuncParam_vs_eta_mc,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
                        "parameter value",
                        "hist_FuncParam_vs_eta_mc","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)

    drawCorrelationPlot(hist_FuncParam_vs_eta_data,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
                        "parameter uncertainty",
                        "hist_FuncParamError_vs_eta_data","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)

    drawCorrelationPlot(hist_FuncParam_vs_eta_mc,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
                        "parameter uncertainty",
                        "hist_FuncParamError_vs_eta_mc","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)

    drawCorrelationPlot(hist_FuncParam_vs_eta_sf,"{lep} #eta".format(lep=lepton),"Pol3 parameter number",
                        "parameter value",
                        "hist_FuncParam_vs_eta_sf","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)

    drawCorrelationPlot(hist_FuncParam_vs_eta_sf,"{lep} #eta".format(lep=lepton),"Pol3 parameter number",
                        "parameter uncertainty",
                        "hist_FuncParamError_vs_eta_sf","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)

    
    c = ROOT.TCanvas("c","",700,700)
    c.SetTickx(1)
    c.SetTicky(1)
    c.cd()
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.06)
    c.cd()
    def quickPlotTH1(c, h, outname,channel, postfix=""):
        c.cd()
        h.Scale(1./h.Integral())
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetXaxis().SetLabelSize(0.06)
        #h.GetXaxis().LabelsOption("v")
        h.GetYaxis().SetTitle("Fraction of events")
        h.GetYaxis().SetTitleOffset(1.15)
        h.GetYaxis().SetTitleSize(0.05)
        h.GetYaxis().SetLabelSize(0.04)
        h.Draw("HIST")
        for ext in ["png","pdf"]:
            c.SaveAs(f"{outname}bestFitFunction{postfix}_{channel}.{ext}")
    if not args.skipEff:
        quickPlotTH1(c, hist_chosenFunc, outname, channel)
    quickPlotTH1(c, hist_chosenFunc_SF, outname, channel, postfix="SF")
        
    # now the chi2 histogram
    # put overflow in last bin
    # but get number of entries before (same as integral since we used Fill() without weights)
    entries_MC = hist_reducedChi2_MC.GetEntries()
    entries_data = hist_reducedChi2_data.GetEntries()
    entries_sf = hist_reducedChi2_sf.GetEntries()

    lastBin = hist_reducedChi2_data.GetNbinsX()
    hist_reducedChi2_data.SetBinContent(lastBin, hist_reducedChi2_data.GetBinContent(lastBin) + hist_reducedChi2_data.GetBinContent(1+lastBin) )
    hist_reducedChi2_data.SetBinError(lastBin, math.sqrt( hist_reducedChi2_data.GetBinError(lastBin)*hist_reducedChi2_data.GetBinError(lastBin) 
                                                     + hist_reducedChi2_data.GetBinError(1+lastBin)*hist_reducedChi2_data.GetBinError(1+lastBin) ) 
                                 )
    lastBin = hist_reducedChi2_MC.GetNbinsX()
    hist_reducedChi2_MC.SetBinContent(lastBin, hist_reducedChi2_MC.GetBinContent(lastBin) + hist_reducedChi2_MC.GetBinContent(1+lastBin) )
    hist_reducedChi2_MC.SetBinError(lastBin, math.sqrt( hist_reducedChi2_MC.GetBinError(lastBin)*hist_reducedChi2_MC.GetBinError(lastBin) 
                                                     + hist_reducedChi2_MC.GetBinError(1+lastBin)*hist_reducedChi2_MC.GetBinError(1+lastBin) ) 
                                 )
    lastBin = hist_reducedChi2_sf.GetNbinsX()
    hist_reducedChi2_sf.SetBinContent(lastBin, hist_reducedChi2_sf.GetBinContent(lastBin) + hist_reducedChi2_sf.GetBinContent(1+lastBin) )
    hist_reducedChi2_sf.SetBinError(lastBin, math.sqrt( hist_reducedChi2_sf.GetBinError(lastBin)*hist_reducedChi2_sf.GetBinError(lastBin) 
                                                     + hist_reducedChi2_sf.GetBinError(1+lastBin)*hist_reducedChi2_sf.GetBinError(1+lastBin) ) 
                                 )

    tmpmin,tmpmax = getMinMaxHisto(hist_reducedChi2_MC, sumError=True)
    tmpmin1,maxY = getMinMaxHisto(hist_reducedChi2_data, sumError=True)
    tmpmin2,maxY2 = getMinMaxHisto(hist_reducedChi2_sf, sumError=True)
    maxY = max([tmpmax,maxY,maxY2])
    maxY = 1.5 * maxY

    hist_reducedChi2_data.GetXaxis().SetTitleOffset(1.2)
    hist_reducedChi2_data.SetLineWidth(2)
    hist_reducedChi2_data.SetLineColor(ROOT.kRed+2)
    hist_reducedChi2_data.SetFillColor(ROOT.kRed+2)
    hist_reducedChi2_data.SetFillStyle(3003)
    hist_reducedChi2_data.GetXaxis().SetTitleSize(0.05)
    hist_reducedChi2_data.GetXaxis().SetLabelSize(0.06)
    hist_reducedChi2_data.GetYaxis().SetTitleOffset(1.15)
    hist_reducedChi2_data.GetYaxis().SetTitleSize(0.05)
    hist_reducedChi2_data.GetYaxis().SetLabelSize(0.04)
    hist_reducedChi2_data.GetYaxis().SetTitle("Events")
    hist_reducedChi2_data.GetYaxis().SetRangeUser(0,maxY)
    hist_reducedChi2_data.GetXaxis().SetTitle("#chi^{2} / NDF")
    hist_reducedChi2_data.Draw("HE")
    lat = ROOT.TLatex()
    xmin = 0.20 
    yhi = 0.85
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(62);
    lat.DrawLatex(xmin,yhi,"")
    lat.DrawLatex(xmin,yhi-0.05,"entries")
    lat.DrawLatex(xmin,yhi-0.1,"mean")
    lat.DrawLatex(xmin,yhi-0.15,"rms")
    line1 = "= {0}".format(int(entries_data))
    line2 = "= {:.2f}".format(hist_reducedChi2_data.GetMean())
    line3 = "= {:.2f}".format(hist_reducedChi2_data.GetStdDev())
    lat.SetTextFont(42);
    lat.SetTextColor(ROOT.kRed+2);
    xmin = xmin + 0.2 
    lat.DrawLatex(xmin,yhi,"data")
    lat.DrawLatex(xmin,yhi-0.05,line1)
    lat.DrawLatex(xmin,yhi-0.1,line2)
    lat.DrawLatex(xmin,yhi-0.15,line3)
    # now MC
    hist_reducedChi2_MC.SetLineWidth(2)
    hist_reducedChi2_MC.SetLineColor(ROOT.kBlack)
    hist_reducedChi2_MC.SetFillColor(ROOT.kGray)
    #hist_reducedChi2_MC.SetFillStyle(3004)
    hist_reducedChi2_MC.Draw("HE SAME")
    #######################
    # redraw some stuff that might be covered by FillColor
    histCopy = hist_reducedChi2_MC.DrawCopy("HE SAME")
    histCopy.SetFillColor(0)
    histCopy.SetFillStyle(0)
    hist_reducedChi2_data.Draw("HE SAME")
    c.RedrawAxis("sameaxis")
    #################
    line1 = "= {0}".format(int(entries_MC))
    line2 = "= {:.2f}".format(hist_reducedChi2_MC.GetMean())
    line3 = "= {:.2f}".format(hist_reducedChi2_MC.GetStdDev())
    lat.SetTextColor(ROOT.kBlack);
    xmin = xmin + 0.2
    lat.DrawLatex(xmin,yhi,"MC")
    lat.DrawLatex(xmin,yhi-0.05,line1)
    lat.DrawLatex(xmin,yhi-0.1,line2)
    lat.DrawLatex(xmin,yhi-0.15,line3)
    # now SF
    hist_reducedChi2_sf.SetLineWidth(2)
    hist_reducedChi2_sf.SetLineColor(ROOT.kBlue)
    hist_reducedChi2_sf.Draw("HE SAME")
    #######################
    # redraw some stuff that might be covered by FillColor
    c.RedrawAxis("sameaxis")
    #################
    line1 = "= {0}".format(int(entries_sf))
    line2 = "= {:.2f}".format(hist_reducedChi2_sf.GetMean())
    line3 = "= {:.2f}".format(hist_reducedChi2_sf.GetStdDev())
    lat.SetTextColor(ROOT.kBlue);
    xmin = xmin + 0.2
    lat.DrawLatex(xmin,yhi,"SF")
    lat.DrawLatex(xmin,yhi-0.05,line1)
    lat.DrawLatex(xmin,yhi-0.1,line2)
    lat.DrawLatex(xmin,yhi-0.15,line3)
    
    for ext in ["png","pdf"]:
        c.SaveAs("{out}reducedChi2_{ch}.{ext}".format(out=outname,ch=channel,ext=ext))

    # before saving things, assign an uncertainty from the original input (will need to devise a better way to estimate them)
    # following histograms have same eta-pt binning
    for ix in range(1,hdataSmoothCheck.GetNbinsX()+1):
        for iy in range(1,hdataSmoothCheck.GetNbinsY()+1):
            ieta = hdata.GetXaxis().FindFixBin(hdataSmoothCheck.GetXaxis().GetBinCenter(ix))
            ipt = hdata.GetYaxis().FindFixBin(hdataSmoothCheck.GetYaxis().GetBinCenter(iy))
            hdataSmoothCheck.SetBinError(ix,iy,hdata.GetBinError(ieta,ipt))            
            hmcSmoothCheck.SetBinError(ix,iy,hmc.GetBinError(ieta,ipt))
            scaleFactor.SetBinError(ix,iy,hsf.GetBinError(ieta,ipt))

    # now I also smooth the scale factors versus eta
    xarray = array('d', scaleFactor.GetXaxis().GetXbins())
    # don't know why I can't get y binning as I did for X axis. Darn you ROOT!
    # yarray = array('d', hist2d.GetYaxis().GetXbins())
    # print yarray
    tmparray = []
    for i in range(1,scaleFactor.GetNbinsY()+2):
        tmparray.append(round(scaleFactor.GetYaxis().GetBinLowEdge(i),4))
    yarray = array('d', tmparray)

    # xarray might not be uniform, so if we want more granularity, we must build the new binning bin by bin
    binSplitFactor = 3  # 3 is good enough, with more splitting, the smooth affects only a narrower strip between two bins
    newxarray = []
    for i in range(len(xarray)-1):
        width = xarray[i+1]-xarray[i]
        for j in range(binSplitFactor):
            newxarray.append(round(xarray[i] + float(j)*width/float(binSplitFactor), 4))
    newxarray.append(xarray[-1])
    #print newxarray   # I suggest you print once in life to see what happens

    # do also a manual smoothing interpolating with a line (so modyfing the two sub-bins at the borders of the original bin)
    scaleFactor_etaInterpolated = ROOT.TH2D("scaleFactor_etaInterpolated","",len(newxarray)-1, array('d',newxarray), len(yarray)-1, yarray)
    # now fill it from the input (which is coarser)
    for ix in range(1,1+scaleFactor_etaInterpolated.GetNbinsX()):
        for iy in range(1,1+scaleFactor_etaInterpolated.GetNbinsY()):
            xval = scaleFactor_etaInterpolated.GetXaxis().GetBinCenter(ix)
            yval = scaleFactor_etaInterpolated.GetYaxis().GetBinCenter(iy)
            hist2xbin = scaleFactor.GetXaxis().FindFixBin(xval)
            hist2ybin = scaleFactor.GetYaxis().FindFixBin(yval)
            scaleFactor_etaInterpolated.SetBinContent(ix,iy,scaleFactor.GetBinContent(hist2xbin,hist2ybin))
            scaleFactor_etaInterpolated.SetBinError(ix,iy,scaleFactor.GetBinError(hist2xbin,hist2ybin))
            # now interpolate eta
            if ix == 1 or ix == scaleFactor_etaInterpolated.GetNbinsX(): continue # do not modify the outer bins
            etabinID = ix%binSplitFactor  # which sub-bin inside the bin (if binSplitFactor=3, can be 1,2,0 for first, second and third sub-bin)
            thisVal = scaleFactor.GetBinContent(hist2xbin,hist2ybin)
            otherVal = 0
            if  etabinID == 1:   
                # if sub-bin on the left, take this -1/3 of the difference between this and the previous (computed in the central sub-bin)      
                otherVal = scaleFactor.GetBinContent(hist2xbin-1,hist2ybin)
                val = thisVal - 1. * (thisVal - otherVal) / 3.
                scaleFactor_etaInterpolated.SetBinContent(ix,iy,val)
            elif etabinID == 0:
                # if sub-bin on the right, take this -1/3 of the difference between this and the following (computed in the central sub-bin)
                otherVal = scaleFactor.GetBinContent(hist2xbin+1,hist2ybin)
                val = thisVal - 1. * (thisVal - otherVal) / 3.
                scaleFactor_etaInterpolated.SetBinContent(ix,iy,val)

    drawCorrelationPlot(scaleFactor_etaInterpolated,
                        "{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactor_etaInterpolated","ForceTitle",
                        outname,palette=args.palette,passCanvas=canvas)

    ###
    ### Do eigen deconvolution on scale factors or efficiencies as appropriate
    ### Should check that covariance matrices are well defined everywhere for this to make sense
    sf3D = None
    # print()
    # print("Running eigen decomposition on scale factors ...")
    # print()
    # sf3D = effStatVariations(outname+"/eigenDecomposition/", hist_FuncCovMatrix_vs_eta_sf, hist_FuncParam_vs_eta_sf,
    #                          nFinePtBins, minPtHistoData, maxPtHistoData, smoothFunction="pol3_tf", suffix="SF", palette=args.palette)
    # sf3D.SetTitle("Nominal in first Z bin, eigen vars elsewhere")
    print()
    if not args.skipEff:
        ## TODO
        pass    
    
    ###########################
    # Now save things
    ###########################
    tfile = ROOT.TFile.Open(outname+outfilename,'recreate')
    hdataSmoothCheck.Write()
    hmcSmoothCheck.Write()
    hsfSmoothCheck.Write()
    scaleFactor.Write()
    scaleFactor_etaInterpolated.Write()
    hsf.Write("scaleFactorOriginal")
    hdata.Write("efficiencyDataOriginal")
    hmc.Write("efficiencyMCOriginal")
    hist_FuncParam_vs_eta_data.Write("hist_FuncParam_vs_eta_data")
    hist_FuncParam_vs_eta_mc.Write("hist_FuncParam_vs_eta_mc")
    hist_FuncParam_vs_eta_sf.Write("hist_FuncParam_vs_eta_sf")
    hist_FuncCovMatrix_vs_eta_data.Write("hist_FuncCovMatrix_vs_eta_data")
    hist_FuncCovMatrix_vs_eta_mc.Write("hist_FuncCovMatrix_vs_eta_mc")
    hist_FuncCovMatrix_vs_eta_sf.Write("hist_FuncCovMatrix_vs_eta_sf")
    # if args.saveTF1:
    #     for key in bestFit_MC:
    #         bestFit_MC[key].Write(key)
    #     for key in bestFit_Data:
    #         bestFit_Data[key].Write(key)
    if sf3D:
        sf3D.Write("SF2D_withEigenVars")
    tfile.Close()
    print()
    print(f"Created file {outname+outfilename}")
    print()

    print("="*30)
    print("Summary of bad fits (Erf for data/MC and pol3 for SF)")
    print("="*30)
    print("### Bad fit status (Data/MC/SF,  key,  fitstatus)")
    for key in sorted(badFitsID_data.keys()):
        print(f"DATA  {key}  {badFitsID_data[key]}")
    for key in sorted(badFitsID_mc.keys()):
        print(f"MC    {key}  {badFitsID_mc[key]}")
    for key in sorted(badFitsID_sf.keys()):
        print(f"SF    {key}  {badFitsID_sf[key]}")
    print("-"*30)
    print("### Bad covariance matrix status (Data/MC/SF,  key,  covquality).")
    for key in sorted(badCovMatrixID_data.keys()):
        print(f"DATA  {key}  {badCovMatrixID_data[key]}")
    for key in sorted(badCovMatrixID_mc.keys()):
        print(f"MC  {key}  {badCovMatrixID_mc[key]}")
    for key in sorted(badCovMatrixID_sf.keys()):
        print(f"SF  {key}  {badCovMatrixID_sf[key]}")
