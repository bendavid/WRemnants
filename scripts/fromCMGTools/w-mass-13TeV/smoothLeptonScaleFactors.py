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
pol2_tf_scaled = None
pol3_tf_scaled = None
pol4_tf_scaled = None
erfPol2_tf_scaled = None
#############################################
# some functions for tensorflow

# for root to be consistent with functions passed to tensorflow
def pol2_root(xvals, parms, xLowVal = 0.0, xFitRange = 1.0):
    xscaled = (xvals[0] - xLowVal) / xFitRange
    return parms[0] + parms[1]*xscaled + parms[2]*xscaled**2

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
    
############

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

#########

def fitTurnOnTF(hist, key, outname, mc, channel="el", hist_chosenFunc=0, drawFit=True, 
                step=None,
                fitRange=None,
                hist_reducedChi2=None,
                #hist_FuncParam_vs_eta=None,
                #hist_FuncCovMatrix_vs_eta=None,  # TH3, eta on x and cov matrix in yz
                charge = "both",
                etabins = [],
                widthPtSmooth=0.2,
                hist_nomiAndAlt_etapt=None,
                histAlt = None
):

    doingSF = True if mc == "SF" else False
        
    chargeText = ""
    if charge == "plus": chargeText = "positive"
    if charge == "minus": chargeText = "negative"
    
    originalMaxPt = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())

    outdir = "{out}{mc}/".format(out=outname,mc=mc)
    createPlotDirAndCopyPhp(outdir)
    adjustSettings_CMS_lumi()
    
    leftMargin = 0.15
    rightMargin = 0.04
    canvas = ROOT.TCanvas(f"canvas_{mc}_{key}","",700,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()                           

    setTDRStyle()
    lowerPanelHeight = 0.3 # for bottom panel with pulls or residuals
    canvas.SetBottomMargin(lowerPanelHeight)
    pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
    pad2.SetTopMargin(1-lowerPanelHeight)
    pad2.SetRightMargin(rightMargin)
    pad2.SetLeftMargin(leftMargin)
    pad2.SetFillColor(0)
    pad2.SetGridy(1)
    pad2.SetFillStyle(0)
    
    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)
    
    hist.GetXaxis().SetLabelSize(0)
    hist.GetXaxis().SetTitle("")
    # hist.GetXaxis().SetTitle(f"{chargeText} muon p_{{T}} (GeV)")
    # hist.GetXaxis().SetTitleOffset(1.1)
    # hist.GetXaxis().SetTitleSize(0.05)
    # hist.GetXaxis().SetLabelSize(0.04)

    maxFitRange = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())
    minFitRange = hist.GetXaxis().GetBinLowEdge(1)
    if fitRange != None:
        if fitRange[1] > 0:
            maxFitRange = fitRange[1]
        if fitRange[0] > 0:
            minFitRange = fitRange[0]
    hist.GetXaxis().SetRangeUser(minFitRange, maxFitRange)

    if mc == "SF":
        hist.GetYaxis().SetTitle("Data/MC scale factor")
    else:
        hist.GetYaxis().SetTitle("{mc} efficiency".format(mc=mc))
    hist.GetYaxis().SetTitleOffset(1.45)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.04)
    if histAlt:
        histAlt.SetStats(0)
        histAlt.SetLineColor(ROOT.kGray+2)
        histAlt.SetMarkerColor(ROOT.kGray+2)
        histAlt.SetMarkerStyle(ROOT.kOpenCircle)
        histAlt.SetMarkerSize(1)
        miny,maxy = getMinMaxMultiHisto([hist, histAlt], sumError=True)
    else:
        miny,maxy = getMinMaxHisto(hist, sumError=True)
    offset = 0.1 * (maxy - miny)
    upOffset = offset * (2.5 if doingSF else 3.5)
    miny -= offset
    maxy += upOffset
    hist.GetYaxis().SetRangeUser(miny, maxy)
    hist.SetStats(0)
    hist.Draw("EP")
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

        if step == "tracking":

            global pol2_tf_scaled
            if pol2_tf_scaled == None:
                pol2_tf_scaled = partial(pol2_root, xLowVal=minFitRange, xFitRange=xFitRange)
            params = np.array([1.0, 0.0, 0.0])
            res_tf1_pol2 = narf.fit_hist(boost_hist, pol2_tf_scaled, params)
            tf1_pol2 = ROOT.TF1("tf1_pol2", pol2_tf_scaled, minFitRange, maxFitRange, len(params))
            tf1_pol2.SetParameters( np.array( res_tf1_pol2["x"], dtype=np.dtype('d') ) )
            tf1_pol2.SetLineWidth(3)
            tf1_pol2.SetLineColor(ROOT.kRed+2)

            # tf1_pol2_test = ROOT.TF1("tf1_pol2_test", pol2_tf_scaled, minFitRange, maxFitRange, len(params))
            # tf1_pol2_test.SetParameters( np.array( res_tf1_pol2["x"], dtype=np.dtype('d') ) )
            # tf1_pol2_test.SetLineStyle(ROOT.kDashed)
            # tf1_pol2_test.SetLineWidth(5)
            # tf1_pol2_test.SetLineColor(ROOT.kBlue)
            # fitopt = "FMBRQS+" # add FM if using Minuit   
            # hist.Fit(tf1_pol2_test, fitopt)

            fitres_TF = {"pol2_tf" : res_tf1_pol2}
            fitFunction = {
                "pol2_tf" : {
                    "func" : tf1_pol2,
                    "leg"  : "Pol2",
                }
            }
            defaultFunc = "pol2_tf"
        else:
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
        ###
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
        ## security check on chi2, to make sure TF1 for plotting is consistent with actual fit    
        manualChi2 = 0.0
        for ib in range(1, 1+hist.GetNbinsX()):
            pt = hist.GetXaxis().GetBinCenter(ib)
            item = (hist.GetBinContent(ib) - fitFunction[fr]["func"].Eval(pt))/hist.GetBinError(ib)
            manualChi2 += item * item
        fitChi2 = fitres_TF[fr]["loss_val"]
        if abs(manualChi2 - fitChi2) > 0.01:
            print(f"{fr}: manual/fit chi2  = {manualChi2}/{fitChi2}")
            
    npar = fitFunction[defaultFunc]["func"].GetNpar()
    
    # if hist_FuncCovMatrix_vs_eta:
    #     for i in range(npar):
    #         for j in range(npar):
    #             hist_FuncCovMatrix_vs_eta.SetBinContent(key+1, i+1, j+1, fitres_TF[defaultFunc]["cov"][i][j])

    # if hist_FuncParam_vs_eta:
    #     # key is the eta bin number, but starts from 0, so add 1
    #     for ip in range(npar):
    #         hist_FuncParam_vs_eta.SetBinContent(key+1, ip+1, fitres_TF[defaultFunc]["x"][ip])
    #         hist_FuncParam_vs_eta.SetBinError(  key+1, ip+1, math.sqrt(fitres_TF[defaultFunc]["cov"][ip][ip]))

    hband = ROOT.TH1D("hband", "", int(math.ceil((maxFitRange-minFitRange))/widthPtSmooth), minFitRange, maxFitRange)
    hband.SetStats(0)
    hband.SetFillColor(ROOT.kGray)
    #hband.SetLineColor(fitFunction[defaultFunc]["func"].GetLineColor())
    #hband.SetFillStyle(3001)
    for ib in range(1, hband.GetNbinsX()+1):
        pt = hband.GetBinCenter(ib)
        val = fitFunction[defaultFunc]["func"].Eval(pt)
        hband.SetBinContent(ib, val)
        hist_nomiAndAlt_etapt.SetBinContent(key+1, ib, 1, val) # assuming the pt binning is the same, which it should, although the code should be made more robust
    #
    # eigen decomposition to plot alternate curves, using the C++ helper
    #
    ## can pass full histogram, eta bin to use is set below
    #systCalc = ROOT.wrem.EtaPtCorrelatedEfficiency(hist_FuncCovMatrix_vs_eta, hist_FuncParam_vs_eta, minFitRange, maxFitRange)
    #systCalc.setSmoothingFunction(defaultFunc)    
    #vecUp = ROOT.std.vector["double"]()
    #vecUp = systCalc.DoEffSyst(key+1)
    #systCalc.setEigenShift(-1.0); # shift down
    #vecDown = systCalc.DoEffSyst(key+1)

    #can also try the eigen decomposition in the python way with numpy
    # e, v = np.linalg.eigh(fitres_TF[defaultFunc]["cov"])
    # altpar_i = fitrs["x"] +/- np.sqrt(e[i])*v[:, i]
    # where altpar_i is a full set of parameters
    
    ## some functions cannot be cloned, although in python it might work
    ##
    #tf1_func_alt = copy.deepcopy(fitFunction[defaultFunc]["func"].Clone("tf1_func_alt"))
    #tf1_func_alt = copy.deepcopy(fitFunction[defaultFunc]["func"])

    # diagonalize and get eigenvalues and eigenvectors
    e, v = np.linalg.eigh(fitres_TF[defaultFunc]["cov"])
    # store all variations for faster access below
    altParameters = np.array([np.zeros(npar, dtype=np.dtype('d'))] * (npar * 2), dtype=np.dtype('d'))
    #print(altParameters)
    for ivar in range(npar):
        shift = np.sqrt(e[ivar]) * v[:, ivar]
        altParameters[ivar]      = fitres_TF[defaultFunc]["x"] + shift
        altParameters[ivar+npar] = fitres_TF[defaultFunc]["x"] - shift

        tf1_func_alt = ROOT.TF1()
    tf1_func_alt.SetName("tf1_func_alt")
    fitFunction[defaultFunc]["func"].Copy(tf1_func_alt)
    tf1_func_alt.SetLineWidth(2)
    for ib in range(1, hband.GetNbinsX()+1):
        pt = hband.GetBinCenter(ib)
        err = 0.0
        for ivar in range(npar):
            # set parameters for a given hessian
            tf1_func_alt.SetParameters(altParameters[ivar]) # this is for Up variations, Down ones could not be the mirror image
            funcVal = tf1_func_alt.Eval(pt)
            diff = funcVal - hband.GetBinContent(ib)
            err += diff * diff
            # now fill TH3, also with down variations
            hist_nomiAndAlt_etapt.SetBinContent(key+1, ib, 2+ivar, funcVal)
            # repeat for Down variations
            tf1_func_alt.SetParameters(altParameters[ivar+npar])
            funcVal = tf1_func_alt.Eval(pt)
            hist_nomiAndAlt_etapt.SetBinContent(key+1, ib, 2+ivar+npar, funcVal)
        err = math.sqrt(err)
        hband.SetBinError(ib, err)
        hist_nomiAndAlt_etapt.SetBinError(key+1, ib, 2+ivar+npar, err)

    hband.Draw("E4SAME")
    # redraw to have them on top
    #for f in fitFunction.keys():
    #    fitFunction[f]["func"].Draw("LSAME")
    hist.Draw("EPSAME")
    if histAlt:
        histAlt.Draw("EPSAME")

    ## draw curves if needed on top of the band
    # colors_alt = [ROOT.kBlue, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+2, ROOT.kGreen+2, ROOT.kViolet, ROOT.kSpring+9, ROOT.kPink]
    # for ivar in range(npar):
    #     #tf1_func_alt.SetLineStyle(ROOT.kDotted)
    #     tf1_func_alt.SetLineColor(colors_alt[ivar])
    #     # set parameters for a given hessian
    #     tf1_func_alt.SetParameters(altParameters[ivar])
    #     tf1_func_alt.DrawCopy("LSAME") # can also not draw these lines, too busy plot otherwise
    #     #tf1_func_alt.SetParameters(altParameters[ivar+npar])
    #     #tf1_func_alt.DrawCopy("LSAME") # can also not draw these lines, too busy plot otherwise

    #######
    nFits = len(fitFunction.keys())

    upLeg = 0.92
    downLeg = max(0.5, upLeg - 0.06 * (nFits + 1)) # +1 to include the eror band as a different entry
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
        # chi2_nsigma = abs(chi2 - ndof) / math.sqrt(2.0 * ndof) 
        # if chi2_nsigma > 3.0:
        #    legEntry += " BAD!"
        chi2prob = ROOT.TMath.Prob(chi2, ndof)
        if chi2prob < 0.05:
            perc_chi2prob = 100.0 * chi2prob
            sign = "="
            if perc_chi2prob < 0.1:
                perc_chi2prob = 0.1
                sign = "<"
            legEntry += " (prob {} {}%)".format(sign, round(perc_chi2prob,1))
        leg.AddEntry(fitFunction[f]["func"], legEntry, 'L')  
        #leg.AddEntry(fitFunction[f]["func"], legEntry, 'L')  
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

    # Now the bottom panel
    pad2.Draw()
    pad2.cd()

    frame = hist.Clone("frame")
    frame.SetTitle("")
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)
    frame.Reset("ICES")
    frame.GetXaxis().SetRangeUser(minFitRange, maxFitRange)
    frame.GetYaxis().SetNdivisions(5)
    denLabel = fitFunction[defaultFunc]["leg"]
    frame.GetYaxis().SetTitle(f"X / {denLabel}") # or Pulls
    frame.GetYaxis().SetTitleOffset(1.45)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().SetTitle(f"{chargeText} muon p_{{T}} (GeV)")

    # to draw Pulls
    # pulls = []
    # defaultIndex = 0
    # for ifunc,f in enumerate(fitFunction.keys()):
    #     pull = copy.deepcopy(hist.Clone(f"pull_{f}"))
    #     pull.Reset("ICESM")
    #     pull.SetLineColor(fitFunction[f]["func"].GetLineColor())
    #     pull.SetLineWidth(2)
    #     pull.SetLineStyle(fitFunction[f]["func"].GetLineStyle())
    #     if f == defaultFunc:
    #         pull.SetFillColor(fitFunction[f]["func"].GetLineColor())
    #         defaultIndex = ifunc
    #     for ib in range(1, 1 + pull.GetNbinsX()):
    #         pull.SetBinContent(ib, (hist.GetBinContent(ib) - fitFunction[f]["func"].Eval(hist.GetBinCenter(ib))) / hist.GetBinError(ib))
    #     pulls.append(pull)

    # miny, maxy = getMinMaxMultiHisto(pulls, excludeEmpty=True, sumError=False,
    #                                        excludeUnderflow=True, excludeOverflow=True)
    # diffy = maxy - miny
    # offset = 0.05
    # miny -= offset * diffy
    # maxy += offset * diffy
    # frame.GetYaxis().SetRangeUser(miny, maxy)
    # frame.Draw()    
    # for ip,p in enumerate(ratios):
    #     drawOpt = "F" if ip == defaultIndex else "HIST"         
    #     p.Draw(f"{drawOpt}SAME")    

    ratios = []
    defaultIndex = 0
    den_noerr = copy.deepcopy(hband.Clone("den_noerr"))
    den = copy.deepcopy(hband.Clone("den"))
    den.SetTitle("")    
    for iBin in range (1, den_noerr.GetNbinsX()+1):
        den_noerr.SetBinError(iBin, 0.)
    den.Divide(den_noerr)
    for ifunc,f in enumerate(fitFunction.keys()):
        if f == defaultFunc:
            continue
        ratio = copy.deepcopy(hband.Clone(f"ratio_{f}"))
        ratio.Reset("ICESM")
        ratio.SetMarkerSize(0)
        ratio.SetMarkerStyle(0)
        ratio.SetLineColor(fitFunction[f]["func"].GetLineColor())
        ratio.SetLineWidth(2)
        ratio.SetFillColor(0)
        ratio.SetLineStyle(fitFunction[f]["func"].GetLineStyle())
        for ib in range(1, 1 + ratio.GetNbinsX()):
            xval = ratio.GetBinCenter(ib)
            ratio.SetBinContent(ib, fitFunction[f]["func"].Eval(xval)/ den_noerr.GetBinContent(ib))
        ratios.append(ratio)
    # now data which is less granular
    dataRatio = copy.deepcopy(hist.Clone("dataRatio"))
    for ib in range(1, 1 + dataRatio.GetNbinsX()):
        ibinDen = hband.GetXaxis().FindFixBin(dataRatio.GetBinCenter(ib))
        denVal = hband.GetBinContent(ibinDen)
        dataRatio.SetBinContent(ib, dataRatio.GetBinContent(ib) / denVal)
        dataRatio.SetBinError(  ib, dataRatio.GetBinError(ib)   / denVal)
        
    miny, maxy = getMinMaxMultiHisto(ratios+[den, dataRatio], excludeEmpty=True, sumError=True,
                                     excludeUnderflow=True, excludeOverflow=True)

    diffy = maxy - miny
    offset = 0.1
    miny -= offset * diffy
    maxy += offset * diffy
    frame.GetYaxis().SetRangeUser(miny, maxy)
    frame.Draw()
    den.Draw("E4SAME")
    denOnlyLine = copy.deepcopy(den.Clone("denOnlyLine"))
    denOnlyLine.SetFillColor(0)
    denOnlyLine.SetLineWidth(1)
    denOnlyLine.SetLineColor(fitFunction[defaultFunc]["func"].GetLineColor())
    denOnlyLine.Draw("HIST SAME")
    for ip,p in enumerate(ratios):
        drawOpt = "C" # like HIST but smooth curve through points
        p.Draw(f"{drawOpt}SAME")    
    dataRatio.Draw("EP SAME")
    pad2.RedrawAxis("sameaxis")

    #---------------------------------
    
    tmpch = ""
    if charge != "both":
        tmpch = "_" + charge
    for ext in ["pdf","png"]:
        if mc == "SF":
            canvas.SaveAs("{out}sf_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,ch=channel,b=key,charge=tmpch,ext=ext))            
        else:
            canvas.SaveAs("{out}eff{mc}_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,charge=tmpch,ext=ext))                            
    
    return fitFunction[defaultFunc]["func"]

############################################################################

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
    parser.add_argument(    '--input-hist-names', dest='inputHistNames', default='', type=str, help='Pass comma separated list of 3  names, for eff(data),eff(MC),SF, to be used instead of the default names')
    parser.add_argument(    '--input-hist-names-alt', dest='inputHistNamesAlt', default='', type=str, help='Pass comma separated list of 2  names for alternate variations, for eff(data),SF, to be used instead of the default names')
    parser.add_argument(     '--palette'  , dest='palette',      default=87, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--skip-eff', dest='skipEff', action="store_true", default=False, help='Skip efficiencies and do only SF (to save time and if one only wants to smooth SF directly)')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    ### TODO:
    ### options to set fit range and some histogram definition in the fit function are not always consistent.
    ### also, one should slice the boost histogram used for the fit, but it is not implemented yet
    ### If one uses the full range everything is fine, otherwise it should be checked
    if any(x > 0 for x in args.ptFitRange):
        print("WARNING: using a restricted fit range is currently not implemented.")
        quit()
    
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

    datahistnameAlt = f"effData_altSig_{args.step}_{args.era}_{args.charge}"
    sfhistnameAlt   = f"SF2D_dataAltSig_{args.step}_{args.era}_{args.charge}"
    if len(args.inputHistNamesAlt):
        datahistnameAlt,sfhistnameAlt = args.inputHistNamesAlt.split(",")
        
    tfile = safeOpenFile(args.inputfile[0])
    hsf =   safeGetObject(tfile, sfhistname)
    hsfAlt = safeGetObject(tfile, sfhistnameAlt)
    if args.skipEff:
        hdata = copy.deepcopy(hsf.Clone("data_dummy"))
        hdata.Reset("ICESM")
        hmc = copy.deepcopy(hdata.Clone("mc_dummy"))
        hdataAlt = copy.deepcopy(hdata.Clone("dataAlt_dummy"))
    else:
        hdata = safeGetObject(tfile, datahistname)
        hmc =   safeGetObject(tfile, mchistname)
        hdataAlt = safeGetObject(tfile, datahistnameAlt)
    tfile.Close()
        
    etabins = [round(hdata.GetXaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hdata.GetNbinsX())]
    ptbins =  [round(hdata.GetYaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hdata.GetNbinsY())]
        
    # dummybins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
    # title = "parameters: Sum_{i=0,..,4} a_{i}*[(x-min)/xrange]**i"
    # hist_FuncParam_vs_eta_data = ROOT.TH2D("hist_FuncParam_vs_eta_data",
    #                                       title,
    #                                       len(etabins)-1, array('d',etabins),
    #                                       len(dummybins)-1, array('d',dummybins))    
    # hist_FuncParam_vs_eta_mc   = ROOT.TH2D("hist_FuncParam_vs_eta_mc",
    #                                       title,
    #                                       len(etabins)-1, array('d',etabins),
    #                                       len(dummybins)-1, array('d',dummybins))

    # hist_FuncCovMatrix_vs_eta_data = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_data",
    #                                           "Covariance matrix: eta on X",
    #                                           len(etabins)-1, array('d',etabins),
    #                                           len(dummybins)-1, array('d',dummybins),
    #                                           len(dummybins)-1, array('d',dummybins))
    # hist_FuncCovMatrix_vs_eta_mc = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_mc",
    #                                         "Covariance matrix: eta on X",
    #                                         len(etabins)-1, array('d',etabins),
    #                                         len(dummybins)-1, array('d',dummybins),
    #                                         len(dummybins)-1, array('d',dummybins))

    # dummybins_sf = [-0.5, 0.5, 1.5, 2.5, 3.5]
    # title_sf = "parameters: Sum_{i=0,..,3} a_{i}*[(x-min)/xrange]**i"
    # hist_FuncParam_vs_eta_sf = ROOT.TH2D("hist_FuncParam_vs_eta_sf",
    #                                      title_sf,
    #                                      len(etabins)-1, array('d',etabins),
    #                                      len(dummybins_sf)-1, array('d',dummybins_sf))    
    # hist_FuncCovMatrix_vs_eta_sf = ROOT.TH3D("hist_FuncCovMatrix_vs_eta_sf",
    #                                          "Covariance matrix: eta on X",
    #                                          len(etabins)-1, array('d',etabins),
    #                                          len(dummybins_sf)-1, array('d',dummybins_sf),
    #                                          len(dummybins_sf)-1, array('d',dummybins_sf))

    # utility histogram to show what function will be used (usually we choose pol4 or pol3 for efficiencies of scale factors)
    hist_chosenFunc = ROOT.TH1D("chosenFitFunc", "Best fit function for each #eta bin for data or MC", 3, 0, 3)
    hist_chosenFunc.GetXaxis().SetBinLabel(1, "erf")
    hist_chosenFunc.GetXaxis().SetBinLabel(2, "tf1_pol3")
    hist_chosenFunc.GetXaxis().SetBinLabel(3, "pol4_tf")
    # 
    hist_chosenFunc_SF = ROOT.TH1D("chosenFitFunc_SF", "Best fit function for each #eta bin for SF", 3, 0, 3)
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(1, "tf1_cheb3")
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(2, "tf1_pol2")
    hist_chosenFunc_SF.GetXaxis().SetBinLabel(3, "pol3_tf")

    hist_reducedChi2_data = ROOT.TH1D("reducedChi2_data", "Reduced #chi^{2}", 25, 0, 5) # will not have many entries (~100 depending on how many eta bins)
    hist_reducedChi2_data.StatOverflows() # use underflow and overflow to compute mean and RMS
    hist_reducedChi2_MC = ROOT.TH1D("reducedChi2_MC", "Reduced #chi^{2}", 25, 0, 5) # will not have many entries (~100 depending on how many eta bins)
    hist_reducedChi2_MC.StatOverflows() # use underflow and overflow to compute mean and RMS
    hist_reducedChi2_sf = ROOT.TH1D("reducedChi2_sf", "Reduced #chi^{2}", 25, 0, 5) # will not have many entries (~100 depending on how many eta bins)
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
    maxPtHisto = hdata.GetYaxis().GetBinLowEdge(1+hdata.GetNbinsY())
    if args.setMaxPtHisto > 0.0:
        maxPtHisto = args.setMaxPtHisto

    nFinePtBins = int(math.ceil(maxPtHisto - hdata.GetYaxis().GetBinLowEdge(1))/args.widthPt)
    minPtHisto = hdata.GetYaxis().GetBinLowEdge(1)

    # get these directly by projecting the TH3 filled below
    #
    # hdataSmoothCheck = ROOT.TH2D("hdataSmoothCheck","Data smoothed efficiency",
    #                              len(etabins)-1, array('d',etabins),
    #                              nFinePtBins, minPtHisto, maxPtHisto)
                                 
    # hmcSmoothCheck = ROOT.TH2D("hmcSmoothCheck","MC smoothed efficiency",
    #                            len(etabins)-1, array('d',etabins),
    #                            nFinePtBins, minPtHisto, maxPtHisto)
                               
    # hsfSmoothCheck = ROOT.TH2D("hsfSmoothCheck","Data/MC smoothed scale factor",
    #                            len(etabins)-1, array('d',etabins),
    #                            nFinePtBins, minPtHisto, maxPtHisto)
                               
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

    ## prepare histograms with nominal eta-pt efficiency or scale factor, and other variations on the Z axis
    ## will save both up and down variations, so have 1+2*nPar bins
    ## assume we have 48 eta bins to use constructor of TH3
    hist_effData_nomiAndAlt_etapt = ROOT.TH3D("hist_effData_nomiAndAlt_etapt",
                                              "Smooth nominal and alternate data efficiency",
                                              48,-2.40,2.40,
                                              nFinePtBins, minPtHisto, maxPtHisto,
                                              11,0.5,11.5)
    hist_effMC_nomiAndAlt_etapt = ROOT.TH3D("hist_effMC_nomiAndAlt_etapt",
                                              "Smooth nominal and alternate MC efficiency",
                                              48,-2.40,2.40,
                                              nFinePtBins, minPtHisto, maxPtHisto,
                                              11,0.5,11.5)
    nBinsSF = 7 if args.step == "tracking" else 9 # fits are done with pol2 for tracking (3 parameters) and pol3 otherwise (4 parameters), then this is 1 + 2 * nBinsSF
    hist_SF_nomiAndAlt_etapt = ROOT.TH3D("hist_SF_nomiAndAlt_etapt",
                                         "Smooth nominal and alternate scale factor",
                                         48,-2.40,2.40,
                                         nFinePtBins, minPtHisto, maxPtHisto,
                                         nBinsSF,0.5,0.5+nBinsSF)

    # hmc and hdata have eta on X and pt on Y
    # we select slices at constant eta and fit along pt with some function

    label = args.step + (args.charge if args.charge != "both" else "")
    hsfpt = make1Dhist("hsfpt", hsf, ptbins, label)
    hsfptAlt = make1Dhist("hsfptAlt", hsfAlt, ptbins, label)

    outfolder_eigenVars = f"{outname}/eigenDecomposition/" 

    if not args.skipEff:
        hmcpt = make1Dhist("hmcpt", hmc, ptbins, label)
        ###########################
        # first MC
        ###########################
        for key in hmcpt:
            bestFitFunc = fitTurnOnTF(hmcpt[key], key, outname, "MC",channel=channel,
                                      hist_chosenFunc=hist_chosenFunc, 
                                      step=args.step,
                                      fitRange=args.ptFitRange,
                                      hist_reducedChi2=hist_reducedChi2_MC,
                                      #hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_mc,
                                      #hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_mc,
                                      charge=args.charge,
                                      etabins=etabins,
                                      widthPtSmooth=args.widthPt,
                                      hist_nomiAndAlt_etapt=hist_effMC_nomiAndAlt_etapt
            )
            #for ipt in range(1, hmcSmoothCheck.GetNbinsY()+1):
            #    ptval = hmcSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #    hmcSmoothCheck.SetBinContent(key+1, ipt, bestFitFunc.Eval(ptval))
            for ipt in range(1, hmcSmoothCheck_origBinPt.GetNbinsY()+1):
                ptval = hmcSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
                hmcSmoothCheck_origBinPt.SetBinContent(key+1, ipt, bestFitFunc.Eval(ptval))
        hmcSmoothCheck = getTH2fromTH3(hist_effMC_nomiAndAlt_etapt, "hmcSmoothCheck", 1, 1)            
        hmcSmoothCheck.SetTitle("Smooth MC efficiency")
        ###########################
        # now data
        ###########################
        hdatapt = make1Dhist("hdatapt", hdata, ptbins, label)
        hdataptAlt = make1Dhist("hdataptAlt", hdataAlt, ptbins, label)
        for key in hdatapt:

            bestFitFunc = fitTurnOnTF(hdatapt[key],key,outname, "Data",channel=channel,hist_chosenFunc=hist_chosenFunc, 
                                      step=args.step,
                                      fitRange=args.ptFitRange,
                                      hist_reducedChi2=hist_reducedChi2_data,
                                      #hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_data,
                                      #hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_data,
                                      charge=args.charge,
                                      etabins=etabins,
                                      widthPtSmooth=args.widthPt,
                                      hist_nomiAndAlt_etapt=hist_effData_nomiAndAlt_etapt,
                                      histAlt=hdataptAlt[key]
            )
            #for ipt in range(1,hdataSmoothCheck.GetNbinsY()+1):
            #    ptval = hdataSmoothCheck.GetYaxis().GetBinCenter(ipt)
            #    hdataSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
            for ipt in range(1,hdataSmoothCheck_origBinPt.GetNbinsY()+1):
                ptval = hdataSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
                hdataSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        hdataSmoothCheck = getTH2fromTH3(hist_effData_nomiAndAlt_etapt, "hdataSmoothCheck", 1, 1)
        hdataSmoothCheck.SetTitle("Smooth data efficiency")

    # ###########################
    # # now SF
    # ###########################
    for key in hsfpt:
        
        bestFitFunc = fitTurnOnTF(hsfpt[key],key,outname, "SF",channel=channel,hist_chosenFunc=hist_chosenFunc_SF, 
                                  step=args.step,
                                  fitRange=args.ptFitRange,
                                  hist_reducedChi2=hist_reducedChi2_sf,
                                  #hist_FuncParam_vs_eta=hist_FuncParam_vs_eta_sf,
                                  #hist_FuncCovMatrix_vs_eta=hist_FuncCovMatrix_vs_eta_sf,
                                  charge=args.charge,
                                  etabins=etabins,
                                  widthPtSmooth=args.widthPt,
                                  hist_nomiAndAlt_etapt=hist_SF_nomiAndAlt_etapt,
                                  histAlt=hsfptAlt[key]
        )
        #for ipt in range(1,hsfSmoothCheck.GetNbinsY()+1):
        #    ptval = hsfSmoothCheck.GetYaxis().GetBinCenter(ipt)
        #    hsfSmoothCheck.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
        for ipt in range(1,hsfSmoothCheck_origBinPt.GetNbinsY()+1):
            ptval = hsfSmoothCheck_origBinPt.GetYaxis().GetBinCenter(ipt)
            hsfSmoothCheck_origBinPt.SetBinContent(key+1,ipt, bestFitFunc.Eval(ptval))
    hsfSmoothCheck = getTH2fromTH3(hist_SF_nomiAndAlt_etapt, "hsfSmoothCheck", 1, 1)
    hsfSmoothCheck.SetTitle("Smooth scale factor")

    #################################
    # start to make plots
    #################################
    zaxisRange = ""
    zaxisRangeSF = "::" + minmaxSF[args.step]

    canvas = ROOT.TCanvas("canvas","",700,625)

    # plot eigen variations
    if not args.skipEff:
        nvars = int((hist_effMC_nomiAndAlt_etapt.GetNbinsZ() - 1) / 2)
        for iv in range(nvars):
            binVar = 2 + iv
            # MC
            hvar = getTH2fromTH3(hist_effMC_nomiAndAlt_etapt, f"effStatVar_p{iv}_effMC",binVar, binVar)
            hvar.SetTitle(f"Eff. stat. nuisance p{iv}")
            hvar.Add(hmcSmoothCheck, -1.0)
            drawCorrelationPlot(hvar, "{lep} #eta".format(lep=lepton), "{lep} p_{{T}} [GeV]".format(lep=lepton), "Alternate - nominal (MC efficiency)",
                                hvar.GetName(), "ForceTitle", outfolder_eigenVars,
                                palette=args.palette, passCanvas=canvas)
            # data
            hvar = getTH2fromTH3(hist_effData_nomiAndAlt_etapt, f"effStatVar_p{iv}_effData",binVar, binVar)
            hvar.SetTitle(f"Eff. stat. nuisance p{iv}")
            hvar.Add(hdataSmoothCheck, -1.0)
            drawCorrelationPlot(hvar, "{lep} #eta".format(lep=lepton), "{lep} p_{{T}} [GeV]".format(lep=lepton), "Alternate - nominal (data efficiency)",
                                hvar.GetName(), "ForceTitle", outfolder_eigenVars,
                                palette=args.palette, passCanvas=canvas)

    # SF
    nvars = int((hist_SF_nomiAndAlt_etapt.GetNbinsZ() - 1) / 2)
    for iv in range(nvars):
        binVar = 2 + iv
        hvar = getTH2fromTH3(hist_SF_nomiAndAlt_etapt, f"effStatVar_p{iv}_SF",binVar, binVar)
        hvar.SetTitle(f"Eff. stat. nuisance p{iv}")
        hvar.Add(hsfSmoothCheck, -1.0)
        drawCorrelationPlot(hvar, "{lep} #eta".format(lep=lepton), "{lep} p_{{T}} [GeV]".format(lep=lepton), "Alternate - nominal (scale factor)",
                            hvar.GetName(), "ForceTitle", outfolder_eigenVars,
                            palette=args.palette, passCanvas=canvas)


    # plot original histograms
    drawCorrelationPlot(hmc,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC efficiency%s" % zaxisRange,
                        "inputEfficiency_MC","",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hdata,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data efficiency%s" % zaxisRange,
                        "inputEfficiency_Data","",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hsf,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC scale factor%s" % zaxisRangeSF,
                        "inputScaleFactor","",outname,palette=args.palette,passCanvas=canvas)
    # now the new ones
    if not args.skipEff:
        drawCorrelationPlot(hmcSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"MC smoothed efficiency%s" % zaxisRange,
                            "smoothEfficiency_MC","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
        drawCorrelationPlot(hdataSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data smoothed efficiency%s" % zaxisRange,
                            "smoothEfficiency_Data","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    drawCorrelationPlot(hsfSmoothCheck,"{lep} #eta".format(lep=lepton),"{lep} p_{{T}} [GeV]".format(lep=lepton),"Data/MC smoothed scale factor%s" % zaxisRangeSF,
                        "smoothScaleFactorDirectly","ForceTitle",outname,palette=args.palette,passCanvas=canvas)
    # scale factor: data/MC
    scaleFactor = ROOT.TH2D("scaleFactor","Scale factor from smoothed efficiencies",
                            len(etabins)-1, array('d',etabins),
                            int(math.ceil((maxPtHisto - hdata.GetYaxis().GetBinLowEdge(1))/args.widthPt)),
                            hdata.GetYaxis().GetBinLowEdge(1),maxPtHisto)                            

    if not args.skipEff:
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
    # drawCorrelationPlot(hist_FuncParam_vs_eta_data,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
    #                     "parameter value",
    #                     "hist_FuncParam_vs_eta_data","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)
    
    # drawCorrelationPlot(hist_FuncParam_vs_eta_mc,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
    #                     "parameter value",
    #                     "hist_FuncParam_vs_eta_mc","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)

    # drawCorrelationPlot(hist_FuncParam_vs_eta_data,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
    #                     "parameter uncertainty",
    #                     "hist_FuncParamError_vs_eta_data","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)

    # drawCorrelationPlot(hist_FuncParam_vs_eta_mc,"{lep} #eta".format(lep=lepton),"Erf[x] parameter number",
    #                     "parameter uncertainty",
    #                     "hist_FuncParamError_vs_eta_mc","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)
    # # now SF
    # drawCorrelationPlot(hist_FuncParam_vs_eta_sf,"{lep} #eta".format(lep=lepton),"Pol3 parameter number",
    #                     "parameter value",
    #                     "hist_FuncParam_vs_eta_sf","ForceTitle",outname,draw_both0_noLog1_onlyLog2=1,palette=args.palette,passCanvas=canvas)

    # drawCorrelationPlot(hist_FuncParam_vs_eta_sf,"{lep} #eta".format(lep=lepton),"Pol3 parameter number",
    #                     "parameter uncertainty",
    #                     "hist_FuncParamError_vs_eta_sf","ForceTitle",outname,draw_both0_noLog1_onlyLog2=2,palette=args.palette,passCanvas=canvas,plotError=True)

    
    c = ROOT.TCanvas("c","",700,700)
    c.SetTickx(1)
    c.SetTicky(1)
    c.cd()
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.06)
    c.cd()
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
    hist_reducedChi2_data.Draw("HIST")
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
    hist_reducedChi2_MC.Draw("HIST SAME")
    #######################
    # redraw some stuff that might be covered by FillColor
    histCopy = hist_reducedChi2_MC.DrawCopy("HIST SAME")
    histCopy.SetFillColor(0)
    histCopy.SetFillStyle(0)
    hist_reducedChi2_data.Draw("HIST SAME")
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
    hist_reducedChi2_sf.Draw("HIST SAME")
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
    # for ix in range(1,hdataSmoothCheck.GetNbinsX()+1):
    #     for iy in range(1,hdataSmoothCheck.GetNbinsY()+1):
    #         ieta = hdata.GetXaxis().FindFixBin(hdataSmoothCheck.GetXaxis().GetBinCenter(ix))
    #         ipt = hdata.GetYaxis().FindFixBin(hdataSmoothCheck.GetYaxis().GetBinCenter(iy))
    #         #if not args.skipEff:
    #         #    hdataSmoothCheck.SetBinError(ix,iy,hdata.GetBinError(ieta,ipt))            
    #         #    hmcSmoothCheck.SetBinError(ix,iy,hmc.GetBinError(ieta,ipt))
    #         scaleFactor.SetBinError(ix,iy,hsf.GetBinError(ieta,ipt))
    
    ###########################
    # Now save things
    ###########################
    hist_postfix = f"_{args.step}_{args.charge}"
    tfile = ROOT.TFile.Open(outname+outfilename,'recreate')
    #hdataSmoothCheck.Write()
    #hmcSmoothCheck.Write()
    #hsfSmoothCheck.Write()
    if not args.skipEff:
        scaleFactor.Write("SF_fromSmoothEfficiencyRatio" + hist_postfix)
    hsf.Write("SF_original" + hist_postfix)
    hdata.Write("effData_original" + hist_postfix)
    hmc.Write("effMC_original" + hist_postfix)
    hist_effData_nomiAndAlt_etapt.Write("effData_nomiAndAlt" + hist_postfix)
    hist_effMC_nomiAndAlt_etapt.Write("effMC_nomiAndAlt" + hist_postfix)
    hist_SF_nomiAndAlt_etapt.Write("SF_nomiAndAlt" + hist_postfix)
    #hist_FuncParam_vs_eta_data.Write("hist_FuncParam_vs_eta_data")
    #hist_FuncParam_vs_eta_mc.Write("hist_FuncParam_vs_eta_mc")
    #hist_FuncParam_vs_eta_sf.Write("hist_FuncParam_vs_eta_sf")
    #hist_FuncCovMatrix_vs_eta_data.Write("hist_FuncCovMatrix_vs_eta_data")
    #hist_FuncCovMatrix_vs_eta_mc.Write("hist_FuncCovMatrix_vs_eta_mc")
    #hist_FuncCovMatrix_vs_eta_sf.Write("hist_FuncCovMatrix_vs_eta_sf")
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
