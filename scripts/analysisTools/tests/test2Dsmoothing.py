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
import pickle
import lz4.frame
import time
from functools import partial

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
logger = logging.setup_logger(__file__, 3, False)

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
            ret += parms2d[dy][dx] * powx * powy
    return ret


# def runTests(args, test=1):
    
#     outdir = args.outdir[0]
#     addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
#     outdir += os.path.basename(args.inputfile[0]).replace(".root","")
#     outdir += "/"
#     createPlotDirAndCopyPhp(outdir)
#     inputfile = args.inputfile[0]
    
#     adjustSettings_CMS_lumi()
    
#     leftMargin = 0.15
#     rightMargin = 0.04
#     bottomMargin = 0.12
#     canvas = ROOT.TCanvas(f"canvas","",800,700)
#     canvas.SetTickx(1)
#     canvas.SetTicky(1)
#     canvas.cd()
#     canvas.SetLeftMargin(leftMargin)
#     canvas.SetBottomMargin(bottomMargin)
#     canvas.SetRightMargin(rightMargin)
#     canvas.cd()                           
    
#     sfhistname = args.histname[0] #"SF3D_nominal_isolation" # uT - eta - pt
#     tfile = safeOpenFile(inputfile)
#     hsf =   safeGetObject(tfile, sfhistname)
#     hsf.GetXaxis().SetRange(2, hsf.GetNbinsX()-1) # remove extreme bins for now, they extend up to infinity

#     # 
#     polnx = args.polDegree[0]
#     polny = args.polDegree[1]

#     isTest = 0

#     if isTest:
#         # just a test    
#         hsf.GetYaxis().SetRange(1, 1)
#         htest = copy.deepcopy(hsf.Project3D("zxe").Clone("htest"))
#         htest.SetTitle(f"Test")
#         print("Preparing test histogram")
#         relativeBinUncertainty = 0.1
#         if isTest == 1:
#             # here when I start with only 3 non zero parameters, removing the pt (i.e. "y") dependence, I get a bias on the other
#             # parameters, which should be 0 but are about 1e-3 or 1e-4. This is with relative uncertainty in each bin of 10%,
#             # if I reduce it to 0.1% the bias gets smaller (parameters which should be 0 are returned as ~1e-14
#             ptLow = hsf.GetZaxis().GetBinLowEdge(1)
#             utLow = hsf.GetXaxis().GetBinLowEdge(2)
#             ptRange = hsf.GetZaxis().GetBinLowEdge(1+hsf.GetNbinsZ()) - ptLow
#             utRange = hsf.GetXaxis().GetBinLowEdge(hsf.GetNbinsX()) - utLow
#             polN_2d_scaled = partial(polN_2d,
#                                      xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
#                                      yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)
#             arr = [2.0, 0.05, 0.2] + [0.0 for i in range(9)]
#             params = np.array(arr)
#             for ix in range(1,1+htest.GetNbinsX()):
#                 bincx = htest.GetXaxis().GetBinCenter(ix)
#                 for iy in range(1,1+htest.GetNbinsY()):
#                     bincy = htest.GetYaxis().GetBinCenter(iy)
#                     htest.SetBinContent(ix, iy, polN_2d_scaled([bincx, bincy], params))
#                     htest.SetBinError(ix, iy, relativeBinUncertainty * htest.GetBinContent(ix, iy))
#             ## reset parameters to see if the fit can find the original ones
#             arr = [1.0 for x in range((polnx+1)*(polny+1))]
#             params = np.array(arr)
#         elif isTest == 2:
#             ptLow = hsf.GetZaxis().GetBinCenter(1)
#             utLow = hsf.GetXaxis().GetBinCenter(2)
#             ptRange = hsf.GetZaxis().GetBinCenter(hsf.GetNbinsZ()) - ptLow
#             utRange = hsf.GetXaxis().GetBinCenter(hsf.GetNbinsX()-1) - utLow
#             polN_2d_scaled = partial(polN_2d,
#                                      xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
#                                      yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)
#             for ix in range(1,1+htest.GetNbinsX()):
#                 bincx = (htest.GetXaxis().GetBinCenter(ix) - utLow) / utRange
#                 for iy in range(1,1+htest.GetNbinsY()):
#                     #bincy = (htest.GetYaxis().GetBinCenter(iy) - ptLow) / ptRange
#                     htest.SetBinContent(ix, iy, 2.0 + 0.05 * bincx + 0.2 * bincx * bincx)
#                     htest.SetBinError(ix, iy, relativeBinUncertainty * htest.GetBinContent(ix, iy))
#             arr = [1.0 for x in range((polnx+1)*(polny+1))]
#             params = np.array(arr)
                    
#         # draw and fit
#         drawCorrelationPlot(htest, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Scale factor",
#                             htest.GetName(), "ForceTitle", outdir,
#                             palette=87, passCanvas=canvas)
        
#         boost_hist = narf.root_to_hist(htest)
#         print("Test fit")
#         params = np.array(arr)
#         res_polN_2d = narf.fitutils.fit_hist(boost_hist, polN_2d_scaled, params)
#         status = res_polN_2d["status"]
#         covstatus = res_polN_2d["covstatus"]
#         print(res_polN_2d["x"])
#         print(f"status/covstatus = {status}/{covstatus}")
#         quit()
#     ## END OF ISTEST


def runSmoothing(inputfile, histname, outdir, step, args):

    #inputfile = args.inputfile[0]
    #histname = args.histname[0]
    #outdir = args.outdir[0]
    #step = args.step[0]
    outdirNew = outdir
    addStringToEnd(outdirNew, "/", notAddIfEndswithMatch=True)
    outdirNew += os.path.basename(inputfile).replace(".root","")
    outdirNew += "/"
    createPlotDirAndCopyPhp(outdirNew)

    logger.info(inputfile)
    logger.info(histname)
    logger.info(step)

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
    
    sfhistname = histname #"SF3D_nominal_isolation" # uT - eta - pt
    tfile = safeOpenFile(inputfile)
    hsf =   safeGetObject(tfile, sfhistname)
    uT_binOffset = 1 # 1 to exclude first bin, use 0 to use all
    utEdges  = [round(hsf.GetXaxis().GetBinLowEdge(i), 1) for i in range(1+uT_binOffset, 2+hsf.GetNbinsX()-uT_binOffset)] # remove extreme bins for now
    etaEdges = [round(hsf.GetYaxis().GetBinLowEdge(i), 1) for i in range(1, 2+hsf.GetNbinsY())]
    ptEdges  = [round(hsf.GetZaxis().GetBinLowEdge(i), 1) for i in range(1, 2+hsf.GetNbinsZ())]
    if uT_binOffset == 0:
        # redefine last uT bins with a sensible range (twice the width of the second to last bin)
        stretch = 2.0
        utEdges[0] = utEdges[1] - stretch * (utEdges[2] - utEdges[1])
        utEdges[-1] = utEdges[-2] - stretch * (utEdges[-3] - utEdges[-2])
        name = hsf.GetName()
        hsf.SetName(f"{hsf.GetName()}_AAA")
        hsfTmp = ROOT.TH3D(name, hsf.GetTitle(),
                           len(utEdges)-1, array('d', utEdges),
                           len(etaEdges)-1, array('d', etaEdges),
                           len(ptEdges)-1, array('d', ptEdges))
        ROOT.wrem.fillTH3fromTH3part(hsfTmp, hsf) # fill new histogram with original values, they have same number of bins
        hsf = hsfTmp # replace the input histogram, to use it later
        logger.warning(f"Setting first and last uT edges to {utEdges[0]} and {utEdges[-1]}")
        
    hsf.GetXaxis().SetRange(1+uT_binOffset, hsf.GetNbinsX()-uT_binOffset) # remove extreme bins for now, they extend up to infinity

    # 
    polnx = args.polDegree[0]
    polny = args.polDegree[1]

    # for consistency with how tensorflow reads the histogram to fit, it is better to use bin centers to define the range.
    # It is then fine even if the fit function extends outside this range
    ptLow = hsf.GetZaxis().GetBinCenter(1)
    ptHigh = hsf.GetZaxis().GetBinCenter(hsf.GetNbinsZ())
    ptRange = ptHigh - ptLow
    utLow = hsf.GetXaxis().GetBinCenter(1+uT_binOffset) # remove first bin, whose range is too extreme
    utHigh = hsf.GetXaxis().GetBinCenter(hsf.GetNbinsX()-uT_binOffset) # remove last bin, whose range is too extreme    
    utRange = utHigh - utLow  
    polN_2d_scaled = partial(polN_2d,
                             xLowVal=utLow, xFitRange=utRange, degreeX=polnx,
                             yLowVal=ptLow, yFitRange=ptRange, degreeY=polny)
    ## save actual histogram edges
    ptEdgeLow = hsf.GetZaxis().GetBinLowEdge(1)
    ptEdgeHigh = hsf.GetZaxis().GetBinLowEdge(1 + hsf.GetNbinsZ())
    utEdgeLow = hsf.GetXaxis().GetBinLowEdge(1+uT_binOffset) # remove first bin, whose range is too extreme
    utEdgeHigh = hsf.GetXaxis().GetBinLowEdge(hsf.GetNbinsX()+1-uT_binOffset) # remove last bin, whose range is too extreme    
    # 1 GeV for recoil and pt
    utNbins = int(utEdgeHigh - utEdgeLow + 0.001) # 1 GeV width
    ptNbins = 5 * int(ptEdgeHigh - ptEdgeLow + 0.001) # 0.2 GeV width

    nEtaBins = hsf.GetNbinsY()
    
    # to store the pull versus eta-pt for bins of uT, for some plots
    hpull_utEtaPt = ROOT.TH3D("hpull_utEtaPt", "",
                              len(utEdges)-1, array('d', utEdges),
                              len(etaEdges)-1, array('d', etaEdges),
                              len(ptEdges)-1, array('d', ptEdges))

    # create final boost histogram with smooth SF, eta-pt-ut-ivar
    axis_eta = hist.axis.Regular(nEtaBins, etaEdges[0], etaEdges[-1], name = "eta", overflow = False, underflow = False)
    axis_pt  = hist.axis.Regular(ptNbins,  ptEdgeLow,   ptEdgeHigh,   name = "pt",  overflow = False, underflow = False)
    axis_ut  = hist.axis.Regular(utNbins,  utEdgeLow,   utEdgeHigh,   name = "ut",  overflow = False, underflow = False)
    axis_var = hist.axis.Integer(0, 1+(1+polnx)*(1+polny), underflow = False, overflow =False, name = "nomi-eigenVars")

    histSF3D = hist.Hist(axis_eta, axis_pt, axis_ut, axis_var,
                         name = f"smoothSF3D_{step}",
                         storage = hist.storage.Weight())
    
    # set initial parameters, starting with constant at unit
    arr = [1.0] + [0.0 for x in range((polnx+1)*(polny+1)-1)]
    ##
    postfix = f"{args.postfix}_" if len(args.postfix) else ""
    etaBinsToRun = args.eta if len(args.eta) else range(1, 1 + hsf.GetNbinsY())

    hpull1D_uTpT = ROOT.TH1D("hpull1D", "", 20, -5, 5)
    hpullSummary_eta_mean  = ROOT.TH1D("hpullSummary_eta_mean",  "Pull distribution mean",  nEtaBins, etaEdges[0], etaEdges[-1])
    hpullSummary_eta_sigma = ROOT.TH1D("hpullSummary_eta_sigma", "Pull distribution width", nEtaBins, etaEdges[0], etaEdges[-1])
    
    for ieta in etaBinsToRun:
    # for ieta in range(1, 2):
        hsf.GetYaxis().SetRange(ieta, ieta)
        h = hsf.Project3D("zxe")
        h.SetName(f"{hsf.GetName()}_eta{ieta}")
        etaLow = round(hsf.GetYaxis().GetBinLowEdge(ieta), 1)
        etaCenter = hsf.GetYaxis().GetBinCenter(ieta)
        etaHigh = round(hsf.GetYaxis().GetBinLowEdge(ieta+1), 1)
        etaRange = f"{etaLow} < #eta < {etaHigh}"
        eta_index = axis_eta.index(etaCenter)
        logger.warning(f"Going to fit eta bin {ieta}, {etaRange}")
        h.SetTitle(etaRange)
        hpull = copy.deepcopy(h.Clone(f"{h.GetName()}_pull2Dfit"))
        hpull.Reset("ICESM")
        hpull1D_uTpT.Reset("ICESM")
        drawCorrelationPlot(h, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Scale factor",
                            h.GetName(), "ForceTitle", outdirNew,
                            palette=87, passCanvas=canvas)
                
        boost_hist = narf.root_to_hist(h)
        params = np.array(arr)
        res_polN_2d = narf.fitutils.fit_hist(boost_hist, polN_2d_scaled, params)
        status = res_polN_2d["status"]
        covstatus = res_polN_2d["covstatus"]
        postfit_params = res_polN_2d['x']
        npar = len(postfit_params)
        fitChi2 = round(res_polN_2d["loss_val"], 1)
        ndof = h.GetNbinsX()*h.GetNbinsY() - npar 
        logger.info(f"postfit_params = {postfit_params}")
        logger.info(f"status/covstatus = {status}/{covstatus}")
        if status != 0 or covstatus != 0:
            logger.error("BAD FIT!!!")
            quit()

        # get chi2    
        chi2text = f"#chi^{{2}} = {fitChi2}/{ndof}"
        chi2prob = ROOT.TMath.Prob(fitChi2, ndof)
        if chi2prob < 0.05:
            perc_chi2prob = 100.0 * chi2prob
            sign = "="
            if perc_chi2prob < 0.1:
                perc_chi2prob = 0.1
                sign = "<"
            chi2text += " (prob {} {}%)".format(sign, round(perc_chi2prob,1))
            
        # do it with boost histograms
        for ix in range(1,1+h.GetNbinsX()):
            bincx = h.GetXaxis().GetBinCenter(ix)
            for iy in range(1,1+h.GetNbinsY()):
                bincy = h.GetYaxis().GetBinCenter(iy)
                fitval = polN_2d_scaled([bincx, bincy], postfit_params)
                pull = (fitval - h.GetBinContent(ix, iy)) / h.GetBinError(ix, iy)
                hpull.SetBinContent(ix, iy, pull)
                hpull_utEtaPt.SetBinContent(ix, ieta, iy, pull)
                hpull1D_uTpT.Fill(pull)

        hpullSummary_eta_mean.SetBinContent(ieta, hpull1D_uTpT.GetMean())
        hpullSummary_eta_mean.SetBinError(ieta, hpull1D_uTpT.GetMeanError())
        hpullSummary_eta_sigma.SetBinContent(ieta, hpull1D_uTpT.GetStdDev())
        hpullSummary_eta_sigma.SetBinError(ieta, hpull1D_uTpT.GetStdDevError())
                
        drawCorrelationPlot(hpull, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Pull: (fit - meas)/meas_unc::-5,5",
                            hpull.GetName(), "ForceTitle", outdirNew,
                            palette=87, nContours=20, passCanvas=canvas)

        hfit_alt = []
        # diagonalize and get eigenvalues and eigenvectors
        logger.info("Diagonalizing covariance matrix ...")
        e, v = np.linalg.eigh(res_polN_2d["cov"])
        postfit_params_alt = np.array([np.zeros(npar, dtype=np.dtype('d'))] * (npar * 2), dtype=np.dtype('d'))
        for ivar in range(npar):
            shift = np.sqrt(e[ivar]) * v[:, ivar]
            postfit_params_alt[ivar]      = postfit_params + shift
            postfit_params_alt[ivar+npar] = postfit_params - shift

        logger.info(f"Creating smooth histogram ...")
        boost_hist_smooth = hist.Hist(axis_ut, axis_pt,
                                      name = f"htmp_fit2D_ieta{ieta}",
                                      storage = hist.storage.Weight())

        dtype = tf.float64
        xvals = [tf.constant(center, dtype=dtype) for center in boost_hist_smooth.axes.centers]
        boost_hist_smooth.values()[...] = polN_2d_scaled(xvals, postfit_params)
        # copy into final histograms with 3D smoothed SF
        # first swap pt and ut axes from uT-pT to pT-uT by projecting onto itself with reshuffled axes
        #axes = [axis_pt.name, axis_ut.name]
        #hswap = boost_hist_smooth.project(*axes)
        histSF3D.values()[eta_index, :, :, 0] = boost_hist_smooth.values().T # hswap.values()[:,:]
        # convert to root for plotting
        hfit = narf.hist_to_root(boost_hist_smooth)
        hfit.SetName(f"fit2D_ieta{ieta}_{postfix}")
        hfit.SetTitle(f"{etaRange}   {chi2text}")
        drawCorrelationPlot(hfit, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", f"Smooth scale factor",
                            hfit.GetName(), "ForceTitle", outdirNew,
                            palette=87, passCanvas=canvas)

        for ivar in range(npar):
            boost_hist_smooth.values()[...] = polN_2d_scaled(xvals, postfit_params_alt[ivar])
            histSF3D.values()[eta_index, :, :, ivar+1] = boost_hist_smooth.values().T
            hfit_alt.append(narf.hist_to_root(boost_hist_smooth))
            hfit_alt[ivar].SetName(f"fit2D_ieta{ieta}_eigen{ivar}_{postfix}")
            hfit_alt[ivar].SetTitle(f"{etaRange}: eigen {ivar}")            

        if args.plotEigenVar:
            outdir_eigen = outdirNew + f"eigenDecomposition/eta_{ieta}/"
            createPlotDirAndCopyPhp(outdir_eigen)
            for ivar in range(npar):
                hratio = copy.deepcopy(hfit_alt[ivar].Clone(f"ratioSF_hfit_alt[ivar].GetName()"))
                hratio.Divide(hfit)
                drawCorrelationPlot(hratio, "Projected recoil u_{T} (GeV)", "Muon p_{T} (GeV)", "Alternate / nominal SF ratio",
                                    hratio.GetName(), "ForceTitle", outdir_eigen,
                                    palette=87, passCanvas=canvas)

        logger.info("-"*30)

    if not len(args.eta):
        outdir_pull = outdirNew + f"/pulls_etapt_utBins/"
        createPlotDirAndCopyPhp(outdir_pull)
        for iut in range(1, 1 + hpull_utEtaPt.GetNbinsX()):
            hpull_utEtaPt.GetXaxis().SetRange(iut, iut)
            hpulletapt = hpull_utEtaPt.Project3D("zye")
            hpulletapt.SetName(f"hpull_etapt_ut{iut}")
            utBinLow = hpull_utEtaPt.GetXaxis().GetBinLowEdge(iut)
            utBinHigh = hpull_utEtaPt.GetXaxis().GetBinLowEdge(iut+1)
            hpulletapt.SetTitle(f"{utBinLow} < u_{{T}} < {utBinHigh}")
            drawCorrelationPlot(hpulletapt, "Muon #eta", "Muon p_{T} (GeV)", "Pull: (fit - meas)/meas_unc::-5,5",
                                hpulletapt.GetName(), "ForceTitle", outdir_pull,
                                palette=87, nContours=20, passCanvas=canvas)
        drawSingleTH1(hpullSummary_eta_mean, "Muon #eta", "Mean of pulls in u_{T}-p_{T} plane::-1, 1",
                      hpullSummary_eta_mean.GetName(), outdir_pull, legendCoords=None,
                      lowerPanelHeight=0.0, drawLineTopPanel=0.0,
                      passCanvas=canvas, skipLumi=True)
        drawSingleTH1(hpullSummary_eta_sigma, "Muon #eta", "Width of pulls in u_{T}-p_{T} plane::0, 2",
                      hpullSummary_eta_sigma.GetName(), outdir_pull, legendCoords=None,
                      lowerPanelHeight=0.0, drawLineTopPanel=1.0,
                      passCanvas=canvas, skipLumi=True)


    return histSF3D


if __name__ == "__main__":
            
    parser = argparse.ArgumentParser()
    #parser.add_argument('inputfile',  type=str, nargs=1, help='input root file with histogram')
    #parser.add_argument('histname',  type=str, nargs=1, help='Histogram name to read from the file')
    parser.add_argument('outdir', type=str, nargs=1, help='output directory to save things')
    #parser.add_argument('step', type=str, nargs=1, help='Step, also to name output histogram (should be parsed from input file though)')
    parser.add_argument('--eta', type=int, nargs="*", default=[], help='Select some eta bins (ID goes from 1 to Neta')
    parser.add_argument('--polDegree', type=int, nargs=2, default=[2, 3], help='Select degree of polynomial for 2D smoothing')
    parser.add_argument('--plotEigenVar', action="store_true", help='Plot eigen variations (it actually produces histogram ratios alt/nomi)')
    parser.add_argument('-p', '--postfix', type=str, default="", help='Postfix for plot names (can be the step name)')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    
    work = []
    work.append(["/home/m/mciprian/isolation3DSFUT.root",     "SF3D_nominal_isolation",     "iso"])
    work.append(["/home/m/mciprian/isonotrigger3DSFVQT.root", "SF3D_nominal_isonotrigger",  "isonotrig"])
    work.append(["/home/m/mciprian/triggerplus3DSFUT.root",   "SF3D_nominal_trigger_plus",  "triggerplus"])
    work.append(["/home/m/mciprian/triggerminus3DSFUT.root",  "SF3D_nominal_trigger_minus", "triggerminus"])

    outdir = args.outdir[0]
    
    resultDict = {}
    for w in work:
        inputfile, histname, step = w
        ret = runSmoothing(inputfile, histname, outdir, step, args)
        resultDict[ret.name] = ret
        
    outfile = outdir + "smoothSF3D.pkl.lz4"
    logger.info(f"Going to store histograms {resultDict.keys()} in file {outfile}")
    time0 = time.time()
    with lz4.frame.open(outfile, 'wb') as f:
        pickle.dump(resultDict, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info(f"Output saved: {time.time()-time0}")
