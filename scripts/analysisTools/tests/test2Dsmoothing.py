#!/usr/bin/env python3

import os, re, array, math
import argparse
from copy import *

import numpy as np
import tensorflow as tf
import hist
import boost_histogram as bh
import narf
import narf.fitutils
import subprocess

from functools import partial

#import utilitiesCMG
#utilities = utilitiesCMG.util()

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from scripts.analysisTools.plotUtils.utility import *

import wremnants

def polN_2d(xvals, parms,
            xLowVal = 0.0, xFitRange = 1.0, degreeX=2,
            yLowVal = 0.0, yFitRange = 1.0, degreeY=3):

    parms2d = tf.reshape(parms, [degreeY+1, degreeX+1])
    
    xscaled = (xvals[0] - xLowVal) / xFitRange
    yscaled = (xvals[1] - yLowVal) / yFitRange
    ret = 0.0
    for dx in range(1+degreeX):
        powx = xscaled**dx if dx else 1.0
        for dy in range(1+degreeY):
            powy = yscaled**dy if dy else 1.0
            ret += (parms2d[dy][dx] * powx) * powy
    return ret


if __name__ == "__main__":
            
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile',  type=str, nargs=1, help='input root file with histogram')
    parser.add_argument('outdir', type=str, nargs=1, help='output directory to save things')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    outdir = args.outdir[0]
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    outdir += os.path.basename(args.inputfile[0]).replace(".root","")
    outdir += "/"
    createPlotDirAndCopyPhp(outdir)

    adjustSettings_CMS_lumi()
    
    leftMargin = 0.15
    rightMargin = 0.04
    bottomMargin = 0.12
    canvas = ROOT.TCanvas(f"canvas","",800,700)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetBottomMargin(bottomMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()                           
    
    sfhistname = "SF3D_nominal_isolation" # uT - eta - pt
    tfile = safeOpenFile(args.inputfile[0])
    hsf =   safeGetObject(tfile, sfhistname)
    hsf.GetXaxis().SetRange(2, hsf.GetNbinsX()-1) # remove extreme bins for now, they extend up to infinity

    # 
    polnx = 2
    polny = 3

    isTest = 0

    if isTest:
        # just a test    
        hsf.GetYaxis().SetRange(1, 1)
        htest = copy.deepcopy(hsf.Project3D("zxe").Clone("htest"))
        htest.SetTitle(f"Test")
        print("Preparing test histogram")
        relativeBinUncertainty = 0.1
        if isTest == 1:
            # here when I start with only 3 non zero parameters, removing the pt (i.e. "y") dependence, I get a bias on the other
            # parameters, which should be 0 but are about 1e-3 or 1e-4. This is with relative uncertainty in each bin of 10%,
            # if I reduce it to 0.1% the bias gets smaller (parameters which should be 0 are returned as ~1e-14
            ptLow = hsf.GetZaxis().GetBinLowEdge(1)
            utLow = hsf.GetXaxis().GetBinLowEdge(2)
            ptRange = hsf.GetZaxis().GetBinLowEdge(1+hsf.GetNbinsZ()) - ptLow
            utRange = hsf.GetXaxis().GetBinLowEdge(hsf.GetNbinsX()) - utLow
            polN_2d_scaled = partial(polN_2d,
                                     xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
                                     yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)
            arr = [2.0, 0.05, 0.2] + [0.0 for i in range(9)]
            params = np.array(arr)
            for ix in range(1,1+htest.GetNbinsX()):
                bincx = htest.GetXaxis().GetBinCenter(ix)
                for iy in range(1,1+htest.GetNbinsY()):
                    bincy = htest.GetYaxis().GetBinCenter(iy)
                    htest.SetBinContent(ix, iy, polN_2d_scaled([bincx, bincy], params))
                    htest.SetBinError(ix, iy, relativeBinUncertainty * htest.GetBinContent(ix, iy))
            ## reset parameters to see if the fit can find the original ones
            arr = [1.0 for x in range((polnx+1)*(polny+1))]
            params = np.array(arr)
        elif isTest == 2:
            ptLow = hsf.GetZaxis().GetBinCenter(1)
            utLow = hsf.GetXaxis().GetBinCenter(2)
            ptRange = hsf.GetZaxis().GetBinCenter(hsf.GetNbinsZ()) - ptLow
            utRange = hsf.GetXaxis().GetBinCenter(hsf.GetNbinsX()-1) - utLow
            polN_2d_scaled = partial(polN_2d,
                                     xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
                                     yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)
            for ix in range(1,1+htest.GetNbinsX()):
                bincx = (htest.GetXaxis().GetBinCenter(ix) - utLow) / utRange
                for iy in range(1,1+htest.GetNbinsY()):
                    #bincy = (htest.GetYaxis().GetBinCenter(iy) - ptLow) / ptRange
                    htest.SetBinContent(ix, iy, 2.0 + 0.05 * bincx + 0.2 * bincx * bincx)
                    htest.SetBinError(ix, iy, relativeBinUncertainty * htest.GetBinContent(ix, iy))
            arr = [1.0 for x in range((polnx+1)*(polny+1))]
            params = np.array(arr)
                    
        # draw and fit
        drawCorrelationPlot(htest, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Scale factor",
                            htest.GetName(), "ForceTitle", outdir,
                            palette=87, passCanvas=canvas)
        
        boost_hist = narf.root_to_hist(htest)
        print("Test fit")
        params = np.array(arr)
        res_polN_2d = narf.fitutils.fit_hist(boost_hist, polN_2d_scaled, params)
        status = res_polN_2d["status"]
        covstatus = res_polN_2d["covstatus"]
        print(res_polN_2d["x"])
        print(f"status/covstatus = {status}/{covstatus}")
        quit()
    ## END OF ISTEST

    ptLow = hsf.GetZaxis().GetBinCenter(1)
    ptHigh = hsf.GetZaxis().GetBinCenter(hsf.GetNbinsZ())
    ptRange = ptHigh - ptLow
    utLow = hsf.GetXaxis().GetBinCenter(2) # remove first bin, whose range is too extreme
    utHigh = hsf.GetXaxis().GetBinCenter(hsf.GetNbinsX()-1) # remove last bin, whose range is too extreme    
    utRange = utHigh - utLow  
    polN_2d_scaled = partial(polN_2d,
                             xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
                             yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)

    # set initial parameters, starting with constant at unit
    arr = [1.0] + [0.0 for x in range((polnx+1)*(polny+1)-1)]
    ##    
    for ieta in range(1, 1 + hsf.GetNbinsY()):
    # for ieta in range(1, 2):
        print(f"Going to fit eta bin {ieta}")
        hsf.GetYaxis().SetRange(ieta, ieta)
        h = hsf.Project3D("zxe")
        h.SetName(f"{hsf.GetName()}_eta{ieta}")
        etaLow = round(hsf.GetYaxis().GetBinLowEdge(ieta), 1)
        etaHigh = round(hsf.GetYaxis().GetBinLowEdge(ieta+1), 1)
        etaRange = f"{etaLow} < #eta < {etaHigh}"
        h.SetTitle(etaRange)
        hpull = copy.deepcopy(h.Clone(f"{h.GetName()}_pull2Dfit"))
        hpull.Reset("ICESM")
        drawCorrelationPlot(h, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Scale factor",
                            h.GetName(), "ForceTitle", outdir,
                            palette=87, passCanvas=canvas)
        boost_hist = narf.root_to_hist(h)
        params = np.array(arr)
        res_polN_2d = narf.fitutils.fit_hist(boost_hist, polN_2d_scaled, params)
        status = res_polN_2d["status"]
        covstatus = res_polN_2d["covstatus"]
        postfit_params = res_polN_2d['x']
        print(f"postfit_params = {postfit_params}")
        print(f"status/covstatus = {status}/{covstatus}")
        if status != 0 or covstatus != 0:
            print("BAD FIT!!!")
            quit()
        for ix in range(1,1+h.GetNbinsX()):
            bincx = h.GetXaxis().GetBinCenter(ix)
            for iy in range(1,1+h.GetNbinsY()):
                bincy = h.GetYaxis().GetBinCenter(iy)
                fitval = polN_2d_scaled([bincx, bincy], postfit_params)
                hpull.SetBinContent(ix, iy, (fitval - h.GetBinContent(ix, iy)) / h.GetBinError(ix, iy))
                #hpull.SetBinError(ix, iy, relativeBinUncertainty * htest.GetBinContent(ix, iy))
        drawCorrelationPlot(hpull, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Pull: (fit - meas)/meas_unc::-5,5",
                            hpull.GetName(), "ForceTitle", outdir,
                            palette=87, nContours=20, passCanvas=canvas)

        # 1 GeV for recoil and pt
        utNbins = int(utRange+0.001)
        ptNbins = int(ptRange+0.001)
        hfit = ROOT.TH2D(f"fit2D_ieta{ieta}", etaRange, utNbins, utLow, utHigh, ptNbins, ptLow, ptHigh)
        print(f"Creating smooth histogram ...")
        for ix in range(1,1+hfit.GetNbinsX()):
            bincx = hfit.GetXaxis().GetBinCenter(ix)
            for iy in range(1,1+hfit.GetNbinsY()):
                bincy = hfit.GetYaxis().GetBinCenter(iy)
                fitval = polN_2d_scaled([bincx, bincy], postfit_params)
                hfit.SetBinContent(ix, iy, fitval)
        drawCorrelationPlot(hfit, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Smooth scale factor",
                            hfit.GetName(), "ForceTitle", outdir,
                            palette=87, passCanvas=canvas)

        print("-"*30)
