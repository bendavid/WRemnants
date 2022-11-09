#!/usr/bin/env python3

# Latest commands (check input file name)
# python w-mass-13TeV/smoothLeptonScaleFactors.py /eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_08Oct2022_binnedInPtEta_mass60to120/allSFs.root /eos/user/m/mciprian/www/WMassAnalysis/TnP/egm_tnp_analysis/results_08Oct2022_binnedInPtEta_mass60to120/smoothLeptonScaleFactors/ -s iso

import ROOT, os, sys, re, array, math
import argparse

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

# for a quick summary at the end
badFitsID_data = {}
badFitsID_mc = {}
badFitsID_sf = {}
badCovMatrixID_data = {}
badCovMatrixID_mc = {}
badCovMatrixID_sf = {}

def getReducedChi2andLabel(func):
    if func.GetNDF():
        reducedChi2 = func.GetChisquare() / func.GetNDF()
        lineChi2 = "#chi^{{2}} = {chi2:.2g} / {ndf}".format(chi2=func.GetChisquare(),ndf=int(func.GetNDF()))
    else:
        reducedChi2 = 0
        lineChi2 = "BAD! #chi^{{2}} = {chi2:.2g} / {ndf}".format(chi2=func.GetChisquare(),ndf=int(func.GetNDF()))

    return float(reducedChi2),lineChi2


def copyHisto(h1, h2):

    if h1.GetDimension() != h2.GetDimension():
        print(f"Error in copyHisto(): histograms have different dimensions. Dim(h1)={h1.GetDimension()}  Dim(h2)={h2.GetDimension()}. Exit")
        quit()

    if h1.GetDimension() == 1:
        for ix in range(h2.GetNbinsX()+2):
                h1.SetBinContent(ix, h2.GetBinContent(ix, iy))
                h1.SetBinError(ix, h2.GetBinError(ix, iy))
    elif h1.GetDimension() == 2:
        for ix in range(h2.GetNbinsX()+2):
            for iy in range(h2.GetNbinsY()+2):
                h1.SetBinContent(ix,iy, h2.GetBinContent(ix, iy))
                h1.SetBinError(ix,iy, h2.GetBinError(ix, iy))
    else:
        print("Error in copyHisto(): function not implemented for dimension > 2. Exit")
        quit()        


def fitTurnOn(hist, key, outname, mc, channel="el", hist_chosenFunc=0, drawFit=True, 
              step=None,
              fitRange=None,
              hist_reducedChi2=0,
              hist_FuncParam_vs_eta=0,
              hist_FuncCovMatrix_vs_eta=0,  # TH3, eta on x and cov matrix in yz
              charge = "both",
              etabins = []
              ):
    
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2","Migrad")
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(2)

    chargeText = ""
    if charge == "plus": chargeText = "positive"
    if charge == "minus": chargeText = "negative"

    # some old options, might be removed
    doingSF = True if mc == "SF" else False
    forcePol3 = False
    forceCheb3Always = True if doingSF else False
    forceErfAlways = False if doingSF else True
    
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
    hist.GetXaxis().SetTitleOffset(1.2)
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
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.04)
    miny,maxy = getMinMaxHisto(hist, sumError=True)
    offset = 0.1 * (maxy - miny)
    miny -= offset
    maxy += offset
    hist.GetYaxis().SetRangeUser(miny, maxy)
    hist.SetStats(0)
    hist.Draw("EP")

    fitopt = "FMBQS+" # add FM if using Minuit   
    maxFitRange = hist.GetXaxis().GetBinLowEdge(1+hist.GetNbinsX())
    minFitRange = hist.GetXaxis().GetBinLowEdge(1)
    if fitRange != None:
        if "R" not in fitopt:
            fitopt = "R" + fitopt
        if fitRange[1] > 0:
            maxFitRange = fitRange[1]
        if fitRange[0] > 0:
            minFitRange = fitRange[0]

    ###################
    # fits
    ####################
    # Erf defined here: https://root.cern.ch/doc/v608/namespaceTMath.html#a44e82992cba4684280c72679f0f39adc
    # Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x 
    tf1_erf = ROOT.TF1("tf1_erf", "[0]*TMath::Erf((x-[1])/[2])", minFitRange, maxFitRange) 
    tf1_erf.SetParameter(0,1.0)
    tf1_erf.SetParameter(1,38)
    tf1_erf.SetParameter(2,3.0)
    tf1_erf.SetLineWidth(2)
    tf1_erf.SetLineColor(ROOT.kOrange+1)

    tf1_pol1 = ROOT.TF1("tf1_pol1","pol1",minFitRange,maxFitRange)
    tf1_pol1.SetLineWidth(2)
    tf1_pol1.SetLineColor(ROOT.kGreen+2)

    tf1_pol2 = ROOT.TF1("tf1_pol2","pol2",minFitRange,maxFitRange)
    tf1_pol2.SetLineWidth(2)
    tf1_pol2.SetLineColor(ROOT.kRed+2)

    tf1_pol3 = ROOT.TF1("tf1_pol3","pol3",minFitRange,maxFitRange)
    tf1_pol3.SetLineWidth(2)
    tf1_pol3.SetLineColor(ROOT.kCyan+1)
    tf1_pol3.SetParameters(1.0, 0.0, 0.0, 0.0);

    tf1_cheb3 = ROOT.TF1("tf1_cheb3","cheb3",minFitRange,maxFitRange)
    tf1_cheb3.SetLineWidth(2)
    #tf1_cheb3.SetLineStyle(2)
    tf1_cheb3.SetLineColor(ROOT.kRed+2)
    tf1_cheb3.SetParameters(1.0, 0.0, 0.0, 0.0);
    tf1_cheb3.SetParLimits(0,  0.0, 2.0);
    tf1_cheb3.SetParLimits(1, -0.5, 0.5);
    tf1_cheb3.SetParLimits(2, -0.5, 0.5);
    tf1_cheb3.SetParLimits(3, -0.5, 0.5);

    tf1_cheb2 = ROOT.TF1("tf1_cheb2","cheb2",minFitRange,maxFitRange)
    tf1_cheb2.SetLineWidth(3)
    tf1_cheb2.SetLineStyle(2)
    tf1_cheb2.SetLineColor(ROOT.kBlue)
    tf1_cheb2.SetParameters(1.0, 0.0, 0.0);
    tf1_cheb2.SetParLimits(0,  0.0, 2.0);
    tf1_cheb2.SetParLimits(1, -0.5, 0.5);
    tf1_cheb2.SetParLimits(2, -0.5, 0.5);

    # hardcoded manual tuning for Erf parameters to fix failing fits
    # was based on efficiencies for SMP-18-012, need to be redone for wmass, but the idea is the one below
    # It doesn't need to be a very optimized setup, fits fail randomly sometimes so the quickest solution is the best
    # if step == "XXX":
    #     if charge == "both":
    #         if mc == "Data":
    #             if any(key == x for x in [13,14,19,21,26,28,33,34]):
    #                 tf1_erf.SetParameter(1,32)
    #             #elif any(key == x for x in [9,23,28]):
    #             #    tf1_erf.SetParameter(1,27)
    #             #    tf1_erf.SetParameter(2,2.0)
    #         elif mc == "MC":
    #             if any(key == x for x in [19,20,26,27,28,34]):
    #                 tf1_erf.SetParameter(1,32)
    #                 #tf1_erf.SetParameter(2,5.0)
    #             elif any(key == x for x in [14]):
    #                 tf1_erf.SetParameter(1,30)
    #     elif charge == "plus":
    #         # following commented was when using uncertainty from TnP txt file
    #         # if mc == "Data":
    #         #     if any(key == x for x in [12,14,22,25,35,40,6,7]):
    #         #         tf1_erf.SetParameter(1,32)
    #         # elif mc == "MC":
    #         #     # following commented was used with point uncertainty equal to stat+syst from TnP
    #         #     #if any(key == x for x in [12,14,22,25,33,35,40]):
    #         #     #    tf1_erf.SetParameter(1,32)
    #         #     if any(key == x for x in [12,14,22,25,33,35]):
    #         #         tf1_erf.SetParameter(1,32)
    #         if mc == "Data":
    #             if any(key == x for x in [12,14,22,25,35,40,6,7, 20, 32, 41, 46]):
    #                 tf1_erf.SetParameter(1,32)
    #         elif mc == "MC":
    #             # following commented was used with point uncertainty equal to stat+syst from TnP
    #             #if any(key == x for x in [12,14,22,25,33,35,40]):
    #             #    tf1_erf.SetParameter(1,32)
    #             if any(key == x for x in [12,14,22,25,33,35, 8, 18, 30, 32]):
    #                 tf1_erf.SetParameter(1,32)

    #     elif charge == "minus":
    #         if mc == "Data":
    #             # following commented was when using uncertainty from TnP txt file
    #             #if any(key == x for x in [2, 18,20,29,39,6, 9]):
    #             #    tf1_erf.SetParameter(1,32)
    #             if any(key == x for x in [2, 18, 29,39,6, 9, 4, 23, 32, 45]):
    #                 # key 9 actually has an evident drop, should be RPC transition 
    #                 # http://mciprian.web.cern.ch/mciprian/wmass/13TeV/scaleFactors_Final/muon/muFullData_trigger_fineBin_noMu50_MINUS/Data/effVsPt_Data_mu_eta9_minus.png
    #                 # maybe here I should force pol3, or assign a large uncertainty
    #                 tf1_erf.SetParameter(1,32)
    #             # elif any(key == x for x in [2]):
    #             #     tf1_erf.SetParameter(1,32)
    #         elif mc == "MC":                    
    #             # following commented was when using uncertainty from TnP txt file
    #             # if any(key == x for x in [0,14,18,20,29,44,47]):
    #             #     tf1_erf.SetParameter(1,32)
    #             # # elif any(key == x for x in [0]):
    #             # #     tf1_erf.SetParameter(1,32)
    #             if any(key == x for x in [0, 18,20,29,44,47, 4, 6, 12, 14, 23, 32, 41, 45]):
    #                 tf1_erf.SetParameter(1,32) 
    #             elif any(key == x for x in [17]):
    #                 tf1_erf.SetParameter(1,34)
    #                 tf1_erf.SetParameter(2,5.0)
                
    erf_fitresPtr = None
    cheb3_fitresPtr = None
    # fit and draw (if required)
    # for SF draw cb3 after pol3, since they might overlap
    if doingSF:
        cheb3_fitresPtr = hist.Fit(tf1_cheb3,fitopt)
        cheb2_fitresPtr = hist.Fit(tf1_cheb2,fitopt)
        pol3_fitresPtr = hist.Fit(tf1_pol3,fitopt+"0") # do not draw pol3 here, but keep until all the mess below is cleaned
    else:
        erf_fitresPtr = hist.Fit(tf1_erf,fitopt)
        pol3_fitresPtr = hist.Fit(tf1_pol3,fitopt)
    #spl = ROOT.TSpline3(hist,"b1e1")
    #spl.SetLineWidth(2)
    #spl.SetLineColor(ROOT.kRed+3)
    #spl.Draw("pclsame")

    fitresPtr = cheb2_fitresPtr if doingSF else erf_fitresPtr
    if fitresPtr != None:
        fitstatus = int(fitresPtr)
        #print "fit status: ", str(fitstatus)
        #print "fit status: ", str(fitresPtr.Status())
        # status is 0 if all is ok (if option M was used, might be 4000 in case the improve command of Minuit failed, but it is ok)
        # without M the fit sometimes fails and should be tuned by hand (option M does it in some case)
        if fitstatus != 0 and fitstatus != 4000: 
            print(f"##### WARNING: FIT HAS STATUS --> {str(fitstatus)}") 
            if mc == "Data":
                badFitsID_data[key] = fitstatus
            elif mc == "MC":
                badFitsID_mc[key] = fitstatus
            else:
                badFitsID_sf[key] = fitstatus
            if doingSF and int(cheb2_fitresPtr) == 0:
                print("##### Cheb2 had good status 0")
        cov = fitresPtr.GetCovarianceMatrix()
        #cov.Print()
        #print "%s" % str(fitresPtr.CovMatrix(1,1))
        #print "Covariance matrix status = ", str(fitresPtr.CovMatrixStatus())
        # covariance matrix status code using Minuit convention : =0 not calculated, =1 approximated, =2 made pos def , =3 accurate 
        if fitresPtr.CovMatrixStatus() != 3: 
            print(f"##### WARNING: COVARIANCE MATRIX HAS QUALITY --> {str(fitresPtr.CovMatrixStatus())}")
            if mc == "Data":
                badCovMatrixID_data[key] = fitresPtr.CovMatrixStatus()
            elif mc == "MC":
                badCovMatrixID_mc[key] = fitresPtr.CovMatrixStatus()
            else:
                badCovMatrixID_sf[key] = fitresPtr.CovMatrixStatus()
            if doingSF and cheb2_fitresPtr.CovMatrixStatus() == 3:
                print("##### Cheb2 had good covariance matrix with quality 3")
                
        if hist_FuncCovMatrix_vs_eta:
            # erf has 3 parameters, pol3 and cheb3 have 4, should be made more general
            npar = 4 if doingSF else 3
            for i in range(npar):
                for j in range(npar):
                    hist_FuncCovMatrix_vs_eta.SetBinContent(key+1,i+1,j+1,fitresPtr.CovMatrix(i,j))

    quit()
    upLeg = 0.4
    downLeg = 0.2
    leftLeg = 0.5
    rightLeg = 0.9
    if mc == "SF":
        upLeg = 0.85
        downLeg = 0.65
        leftLeg = 0.6
        rightLeg = 0.9
        
    leg = ROOT.TLegend(leftLeg, downLeg, rightLeg, upLeg)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    legEntry = {}
    legEntry[tf1_erf.GetName()]  = "Erf[x]"
    legEntry[tf1_pol3.GetName()] = "pol3"
    legEntry[tf1_cheb3.GetName()] = "cheb3"
    legEntry[tf1_cheb2.GetName()] = "cheb2"
    legEntry[tf1_pol1.GetName()] = "pol1"
    legEntry[tf1_pol2.GetName()] = "pol2"

    if doingSF:
        leg.AddEntry(tf1_cheb3,  legEntry[tf1_cheb3.GetName()], 'LF')
        leg.AddEntry(tf1_cheb2,  legEntry[tf1_cheb2.GetName()], 'LF')  
    else:
        leg.AddEntry(tf1_erf,  legEntry[tf1_erf.GetName()], 'LF')
        leg.AddEntry(tf1_pol3, legEntry[tf1_pol3.GetName()], "LF")
    #leg.AddEntry(spl,"spline", "LF")
    leg.Draw('same')

    canvas.RedrawAxis("sameaxis")

    setTDRStyle()
    ROOT.gStyle.SetOptTitle(1)  # use histogram title with binning as canvas title

    fit_erf =  hist.GetFunction(tf1_erf.GetName()) if not doingSF else None
    fit_pol3 = hist.GetFunction(tf1_pol3.GetName())
    fit_cheb3 = hist.GetFunction(tf1_cheb3.GetName()) if doingSF else None
    fit_cheb2 = hist.GetFunction(tf1_cheb2.GetName()) if doingSF else None
    
    functions = {}
    functions[tf1_erf.GetName()] = fit_erf
    functions[tf1_pol3.GetName()] = fit_pol3
    functions[tf1_cheb3.GetName()] = fit_cheb3
    functions[tf1_cheb2.GetName()] = fit_cheb2
    
    lat = ROOT.TLatex()
    line = ""
    lineChi2 = ""
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.SetTextColor(1);

    chi2 = 1000000.0
    funcMinChi2 = None
    for name in functions.keys():
        if name == "tf1_pol3": continue
        f = functions[name]
        #print "Name: %s func %s" % (name, f)                                                                                                                    
        if f == None: continue
        if f.GetNDF() == 0: continue
        if f.GetChisquare() < chi2:
            chi2 = f.GetChisquare()
            funcMinChi2 = f
    if funcMinChi2 == None:
        print("ERROR: no function found with at least 1 degree of freedom")
        quit()
    #print("Function %s had the best Chi2/Ndof: %.3f/%d among non-pol3" % (funcMinChi2.GetName(),funcMinChi2.GetChisquare(),funcMinChi2.GetNDF()))
    #print("pol3 Chi2/Ndof: %.3f/%d" % (fit_pol3.GetChisquare(),fit_pol3.GetNDF()))
    if forcePol3:
        retFunc = fit_pol3
    else:
        #print(funcMinChi2.GetName())
        nChi2Sigma = abs(funcMinChi2.GetChisquare()-funcMinChi2.GetNDF())/math.sqrt(2.0*funcMinChi2.GetNDF())  # Chi2 variance is 2*Ndof
        nChi2Sigma_pol3 = abs(fit_pol3.GetChisquare()-fit_pol3.GetNDF())/math.sqrt(2.0*fit_pol3.GetNDF()) if fit_pol3.GetNDF() else 999
        # pol3 will generally fit very well also in case of weird points
        # for good looking points, pol3 might be better because it can change curvature, while other functions cannot (which would be more physical)
        # allow non-pol3 fit to have Chi2 within 2 standard deviation from the expected one
        # in this case choose that value, otherwise use the one closer to expected Chisquare
        if nChi2Sigma < 5.0:
            retFunc = funcMinChi2
            #print("HERE")
        elif nChi2Sigma_pol3 < nChi2Sigma:
            retFunc = fit_pol3
        else:
            retFunc = funcMinChi2
            
    line = "Best fit: " + legEntry[retFunc.GetName()] + (" (forced)" if forcePol3 else "")
    if forceErfAlways and retFunc.GetName() != fit_erf.GetName():
            retFunc = fit_erf
            line = "Best fit: Erf[x] (forced)"
    if forceCheb3Always and retFunc.GetName() != fit_cheb3.GetName():
            retFunc = fit_cheb3
            line = "Best fit: cheb3 (forced)"
            
    reducedChi2, lineChi2 = getReducedChi2andLabel(retFunc)
    if hist_reducedChi2:
        hist_reducedChi2.Fill(reducedChi2)
    if hist_chosenFunc:
        hist_chosenFunc.Fill(retFunc.GetName(), 1)
    
    xmin = 0.20 
    yhi = 0.85
    lat.DrawLatex(xmin, yhi, line);
    lat.DrawLatex(xmin, yhi-0.05, lineChi2);
    tmpch = ""
    if charge != "both":
        tmpch = "_" + charge
    for ext in ["pdf","png"]:
        if mc == "SF":
            canvas.SaveAs("{out}sf_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,ch=channel,b=key,charge=tmpch,ext=ext))            
        else:
            canvas.SaveAs("{out}eff{mc}_pt_{ch}_eta{b}{charge}.{ext}".format(out=outdir,mc=mc,ch=channel,b=key,charge=tmpch,ext=ext))                            

    if hist_FuncParam_vs_eta:
        # key is the eta bin number, but starts from 0, so add 1
        if doingSF:
            for ip in range(retFunc.GetNpar()):
                hist_FuncParam_vs_eta.SetBinContent(key+1, ip+1, retFunc.GetParameter(ip))
                hist_FuncParam_vs_eta.SetBinError(  key+1, ip+1, retFunc.GetParError(ip))
        else:
            if retFunc.GetName() == tf1_erf.GetName():
                for ip in range(retFunc.GetNpar()):
                    hist_FuncParam_vs_eta.SetBinContent(key+1, ip+1, retFunc.GetParameter(ip))
                    hist_FuncParam_vs_eta.SetBinError(  key+1, ip+1, retFunc.GetParError(ip))
    
    # some last moment hand-fix for some bins if needed
    # if False:
    #     message = f">>> Warning: returning pol3 for key {key} in {mc} as interpolating function"
    #     if charge == "plus":
    #         if (mc == "MC" and any(key == x for x in [9,38])) or (mc == "Data" and any(key == x for x in [38])):
    #             print(message)
    #             return fit_pol3
    #     elif charge == "minus":        
    #         if (mc == "MC" and any(key == x for x in [9,38])) or (mc == "Data" and any(key == x for x in [9])):
    #             print(message)
    #             return fit_pol3

    return retFunc


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

    tf = ROOT.TFile.Open(fname)        
    hist = tf.Get(hname)
    if (hist == 0):
        print(f"Error: could not retrieve hist from input file {fname}. Exit")
        quit()
    hist.SetDirectory(0)
    tf.Close()

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

    tf = ROOT.TFile.Open(outname+outfilename,'recreate')
    histSmooth.Write()
    hist.Write("histOriginal")
    tf.Close()
    print()
    print(f"Created file {outname+outfilename}")
    print()
    return 0
    
############
    
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
    parser.add_argument(     '--save-TF1', dest='saveTF1',action="store_true", default=False, help='Save TF1 as well, not just TH2 with many bins')
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
    
    tf = safeOpenFile(args.inputfile[0])
    hdata = safeGetObject(tf, datahistname)
    hmc =   safeGetObject(tf, mchistname)
    hsf =   safeGetObject(tf, sfhistname)
    tf.Close()
        
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
    # let's select an error function 

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

    label = args.step + (args.charge if args.charge != "both" else "")
    hmcpt = make1Dhist("hmcpt", hmc, ptbins, label)
    hdatapt = make1Dhist("hdatapt", hdata, ptbins, label)
    hsfpt = make1Dhist("hsfpt", hsf, ptbins, args.step)

    bestFit_MC = {}
    bestFit_Data = {}
    if not args.skipEff:
        ###########################
        # first MC
        ###########################
        for key in hmcpt:
            bestFitFunc = fitTurnOn(hmcpt[key], key, outname, "MC",channel=channel,
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

            bestFitFunc = fitTurnOn(hdatapt[key],key,outname, "Data",channel=channel,hist_chosenFunc=hist_chosenFunc, 
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
        
        bestFitFunc = fitTurnOn(hsfpt[key],key,outname, "SF",channel=channel,hist_chosenFunc=hist_chosenFunc_SF, 
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
    if args.skipEff:
        print()
        print("Running eigen decomposition ...")
        print()
        sf3D = effStatVariations(outname+"/eigenDecomposition/", hist_FuncCovMatrix_vs_eta_sf, hist_FuncParam_vs_eta_sf,
                                 nFinePtBins, minPtHistoData, maxPtHistoData, smoothFunction="cheb3", suffix="SF", palette=args.palette)
        sf3D.SetTitle("Nominal in first Z bin, eigen vars for cheb3 elsewhere")
        print()
        
    ###########################
    # Now save things
    ###########################
    tf = ROOT.TFile.Open(outname+outfilename,'recreate')
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
    if args.saveTF1:
        for key in bestFit_MC:
            bestFit_MC[key].Write(key)
        for key in bestFit_Data:
            bestFit_Data[key].Write(key)
    if sf3D:
        sf3D.Write("SF2D_withEigenVars")
    tf.Close()
    print()
    print(f"Created file {outname+outfilename}")
    print()

    print("="*30)
    print("Summary of bad fits (Erf for data/MC and pol3 for SF)")
    print("="*30)
    print("### Bad fit status (Data/MC,  key,  fitstatus)")
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
