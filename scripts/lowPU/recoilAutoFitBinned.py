
import sys,array,math,os,copy

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import functions
import plotter

mirror = False

def mirrorx(h):

    hnew = h.Clone("%s_mirrorx" % h.GetName())
    for i in range(1, hnew.GetNbinsX()+1):
        idx = hnew.GetNbinsX()-i+1
        hnew.SetBinContent(i, h.GetBinContent(idx))
        hnew.SetBinError(i, h.GetBinError(idx))

    #sys.exit()
    return hnew



def doFitMC_3Gauss(tag, name, hOut, qTbin, qTmin, qTmax, mean1_, mean2_, mean3_, norm1_, norm2_, sigma1_, sigma2_, sigma3_):

    hist = fIn.Get(name)
    hist.Rebin(1)
    
    if qTmin < 30: xMin,xMax = -100, 100
    elif qTmin < 60: xMin,xMax = -150, 50
    elif qTmin < 150: xMin,xMax = -200, 50
    else: xMin,xMax = -500, 0
    if "perp" in name: xMin,xMax = -100, 100
    if "_qT" in name: xMin,xMax = -100, 100
    
    if mirror:
        hist = mirrorx(hist)
        xMax = -xMax
        xMin = -xMin
        tmp = xMax
        xMax = xMin
        xMin = tmp
    
    #hist = hist.Rebin(1)
    hist.Scale(1./hist.Integral())
    
    qT = 0.5*(qTmin + qTmax)
    cfg['xmin'], cfg['xmax'] = xMin, xMax
    cfg['ymax'] = 1.75*hist.GetMaximum()
    cfg['ymax'] = 1e3
    cfg['ymin'] = 1e-6

    recoil = ROOT.RooRealVar("recoil", "Recoil parallel (GeV)", 0, -500, 500)
    mean1 = ROOT.RooRealVar("mean1", "", mean1_, -500, 500)
    mean2 = ROOT.RooRealVar("mean2", "", mean2_, -500, 500)
    mean3 = ROOT.RooRealVar("mean3", "", mean3_, -500, 500)
    sigma1 = ROOT.RooRealVar("sigma1", "", sigma1_, 0.1, 100)
    sigma2 = ROOT.RooRealVar("sigma2", "", sigma2_, 0.1, 100)
    sigma3 = ROOT.RooRealVar("sigma3", "", sigma3_, 0.1, 100)
    norm1 = ROOT.RooRealVar("norm1", "", norm1_, 0, 1)
    norm2 = ROOT.RooRealVar("norm2", "", norm2_, 0, 1)
    gauss1 = ROOT.RooGaussian("gauss1", "gauss1", recoil, mean1, sigma1)
    gauss2 = ROOT.RooGaussian("gauss2", "gauss2", recoil, mean2, sigma2)
    gauss3 = ROOT.RooGaussian("gauss3", "gauss3", recoil, mean3, sigma3)

    pdf = ROOT.RooAddPdf("pdf", '', ROOT.RooArgList(gauss1, gauss2, gauss3), ROOT.RooArgList(norm1, norm2))
    rdh = ROOT.RooDataHist("rdh", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(hist))
    

    if "perp" in name:
    
        #mean1.setVal(0)
        #mean1.setConstant(ROOT.kTRUE)
        #mean2.setVal(0)
        #mean2.setConstant(ROOT.kTRUE)
        #mean3.setVal(0)
        #mean3.setConstant(ROOT.kTRUE)
        
        norm1.setVal(0.5)
        norm1.setConstant(ROOT.kTRUE)
        
        norm2.setVal(0.25)
        norm2.setConstant(ROOT.kTRUE)

    # ML fit
    fitRes = pdf.fitTo(rdh, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save(ROOT.kTRUE)) # ROOT.RooFit.Extended(ROOT.kTRUE), 
   
    #mean2.setVal(mean1.getVal())
    #mean3.setVal(mean1.getVal())    


    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=0)
    
    
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    dummyT.Draw("HIST")
    
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    plt = recoil.frame()
    rdh.plotOn(plt)
    gauss1.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(norm1.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    gauss2.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Normalization(norm2.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    gauss3.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Normalization(1.-norm1.getVal()-norm2.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    pdf.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
    histpull = plt.pullHist()
    chi2 = plt.chiSquare() 

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.20, 0.85, label)
    latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin, qTmax))
    latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist.GetMean(), hist.GetMeanError()))
    latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist.GetRMS(), hist.GetRMSError()))
        
    latex.DrawLatex(0.60, 0.85, "#mu_{1} = %.3f #pm %.3f" % (mean1.getVal(), mean1.getError()))
    latex.DrawLatex(0.60, 0.80, "#mu_{2} = %.3f #pm %.3f" % (mean2.getVal(), mean2.getError()))
    latex.DrawLatex(0.60, 0.75, "#mu_{3} = %.3f #pm %.3f" % (mean3.getVal(), mean3.getError()))
    latex.DrawLatex(0.60, 0.70, "#sigma_{1} = %.3f #pm %.3f" % (sigma1.getVal(), sigma1.getError()))
    latex.DrawLatex(0.60, 0.65, "#sigma_{2} = %.3f #pm %.3f" % (sigma2.getVal(), sigma2.getError()))
    latex.DrawLatex(0.60, 0.60, "#sigma_{3} = %.3f #pm %.3f" % (sigma3.getVal(), sigma3.getError()))
    latex.DrawLatex(0.60, 0.55, "N_{1} = %.3f #pm %.3f" % (norm1.getVal(), norm1.getError()))
    latex.DrawLatex(0.60, 0.50, "N_{2} = %.3f #pm %.3f" % (norm2.getVal(), norm2.getError()))
    latex.DrawLatex(0.60, 0.45, "#chi^{2}/ndof = %.3f" % chi2)

    plt.Draw("SAME")
    plotter.auxRatio()
       
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
    plt = recoil.frame()
    plt.addPlotable(histpull, "P")
    plt.Draw("SAME")
    dummyL.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s/%03d_recoil.png" % (outDir, tag, qTbin))
    #canvas.SaveAs("%s/%03d_recoil_%s.pdf" % (outDir, qTbin, name))

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)   

    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()
    
    
    hOut.SetBinContent(qTbin, 1, 1, mean1.getVal())
    hOut.SetBinContent(qTbin, 1, 2, mean1.getError())
    hOut.SetBinContent(qTbin, 2, 1, sigma1.getVal())
    hOut.SetBinContent(qTbin, 2, 2, sigma1.getError())
    hOut.SetBinContent(qTbin, 3, 1, norm1.getVal())
    hOut.SetBinContent(qTbin, 3, 2, norm1.getError())
    
    hOut.SetBinContent(qTbin, 4, 1, mean2.getVal())
    hOut.SetBinContent(qTbin, 4, 2, mean2.getError())
    hOut.SetBinContent(qTbin, 5, 1, sigma2.getVal())
    hOut.SetBinContent(qTbin, 5, 2, sigma2.getError())
    hOut.SetBinContent(qTbin, 6, 1, norm2.getVal())
    hOut.SetBinContent(qTbin, 6, 2, norm2.getError())
    
    hOut.SetBinContent(qTbin, 7, 1, mean3.getVal())
    hOut.SetBinContent(qTbin, 7, 2, mean3.getError())
    hOut.SetBinContent(qTbin, 8, 1, sigma3.getVal())
    hOut.SetBinContent(qTbin, 8, 2, sigma3.getError())
    hOut.SetBinContent(qTbin, 9, 1, (1. - norm1.getVal() - norm2.getVal()))
    hOut.SetBinContent(qTbin, 9, 2, 0)
    
    return mean1.getVal(), mean2.getVal(), mean3.getVal(), norm1.getVal(), norm2.getVal(), sigma1.getVal(), sigma2.getVal(), sigma3.getVal()
 



def doFitBkg_3Gauss(tag, name, qTbin, qTmin, qTmax, mean1_, mean2_, mean3_, norm1_, norm2_, sigma1_, sigma2_, sigma3_):

    hist = fIn.Get(name)
    
    #hist.Rebin(2)
    hist = functions.Rebin(hist, bins_perp, binWidth=False)
    
    if qTmin < 30: xMin,xMax = -100, 100
    elif qTmin < 60: xMin,xMax = -150, 50
    elif qTmin < 150: xMin,xMax = -200, 50
    else: xMin,xMax = -500, 0
    if "perp" in name: xMin,xMax = -300, 300
    if "_qT" in name: xMin,xMax = -200, 200
    
    #hist = hist.Rebin(1)
    hist_norm = hist.Integral()
    hist.Scale(1./hist_norm)    
    
    if mirror:
        hist = mirrorx(hist)
        xMax = -xMax
        xMin = -xMin
        tmp = xMax
        xMax = xMin
        xMin = tmp

    
    qT = 0.5*(qTmin + qTmax)
    


    cfg['xmin'], cfg['xmax'] = xMin, xMax
    cfg['ymax'] = 1.75*hist.GetMaximum()
    cfg['ymax'] = 1e5
    cfg['ymin'] = 1e-5
    


    recoil = ROOT.RooRealVar("recoil", "Recoil parallel (GeV)", 0, -500, 500)
    mean1 = ROOT.RooRealVar("mean1", "", mean1_, -500, 500)
    mean2 = ROOT.RooRealVar("mean2", "", mean2_, -500, 500)
    mean3 = ROOT.RooRealVar("mean3", "", mean3_, -500, 500)
    sigma1 = ROOT.RooRealVar("sigma1", "", sigma1_, 0.1, 100)
    sigma2 = ROOT.RooRealVar("sigma2", "", sigma2_, 0.1, 100)
    sigma3 = ROOT.RooRealVar("sigma3", "", sigma3_, 0.1, 100)
    norm1 = ROOT.RooRealVar("norm1", "", norm1_, 0, 1)
    norm2 = ROOT.RooRealVar("norm2", "", norm2_, 0, 1)
    gauss1 = ROOT.RooGaussian("gauss1", "gauss1", recoil, mean1, sigma1)
    gauss2 = ROOT.RooGaussian("gauss2", "gauss2", recoil, mean2, sigma2)
    gauss3 = ROOT.RooGaussian("gauss3", "gauss3", recoil, mean3, sigma3)
    

    
   
    pdf = ROOT.RooAddPdf("pdf", '', ROOT.RooArgList(gauss1, gauss2, gauss3), ROOT.RooArgList(norm1, norm2))
    rdh = ROOT.RooDataHist("rdh", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(hist))

    if "perp" in name:
        mean1.setVal(0)
        mean1.setConstant(ROOT.kTRUE)
        mean2.setVal(0)
        mean2.setConstant(ROOT.kTRUE)
        mean3.setVal(0)
        mean3.setConstant(ROOT.kTRUE)
        
        
    fitRes = pdf.fitTo(rdh, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save(ROOT.kTRUE))
    
    #mean2.setVal(mean1.getVal())
    #mean3.setVal(mean1.getVal())

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=0)
    
    
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    dummyT.Draw("HIST")
    
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    plt = recoil.frame()
    rdh.plotOn(plt)
    pdf.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
    histpull = plt.pullHist()
    chi2 = plt.chiSquare() 

        
        
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.20, 0.85, label)
    latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin, qTmax))
    latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist.GetMean(), hist.GetMeanError()))
    latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist.GetRMS(), hist.GetRMSError()))
        
    latex.DrawLatex(0.60, 0.85, "#mu_{1} = %.3f #pm %.3f" % (mean1.getVal(), mean1.getError()))
    latex.DrawLatex(0.60, 0.80, "#mu_{2} = %.3f #pm %.3f" % (mean2.getVal(), mean2.getError()))
    latex.DrawLatex(0.60, 0.75, "#mu_{3} = %.3f #pm %.3f" % (mean3.getVal(), mean3.getError()))
    latex.DrawLatex(0.60, 0.70, "#sigma_{1} = %.3f #pm %.3f" % (sigma1.getVal(), sigma1.getError()))
    latex.DrawLatex(0.60, 0.65, "#sigma_{2} = %.3f #pm %.3f" % (sigma2.getVal(), sigma2.getError()))
    latex.DrawLatex(0.60, 0.60, "#sigma_{3} = %.3f #pm %.3f" % (sigma3.getVal(), sigma3.getError()))
    latex.DrawLatex(0.60, 0.55, "N_{1} = %.3f #pm %.3f" % (norm1.getVal(), norm1.getError()))
    latex.DrawLatex(0.60, 0.50, "N_{2} = %.3f #pm %.3f" % (norm2.getVal(), norm2.getError()))
    latex.DrawLatex(0.60, 0.45, "#chi^{2}/ndof = %.3f" % chi2)

    plt.Draw("SAME")
    plotter.auxRatio()
       
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
    
    plt = recoil.frame()
    plt.addPlotable(histpull, "P")

    plt.Draw("SAME")
        
    dummyL.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s/%03d_recoil.png" % (outDir, tag, qTbin))
    #canvas.SaveAs("%s/%03d_recoil_%s.pdf" % (outDir, qTbin, name))

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)   

    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()
    
    return mean1.getVal(), mean2.getVal(), mean3.getVal(), norm1.getVal(), norm2.getVal(), sigma1.getVal(), sigma2.getVal(), sigma3.getVal(), hist_norm



def doFitData_3Gauss(tag, name, hOut, qTbin, qTmin, qTmax, mean1_, mean2_, mean3_, norm1_, norm2_, sigma1_, sigma2_, sigma3_, mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg, hist_norm_bkg):

    hist = fIn.Get(name)
    hist = functions.Rebin(hist, bins_perp, binWidth=False)

    hist_bkg = fIn.Get(name.replace("data_", "bkg_"))
    hist_bkg = functions.Rebin(hist_bkg, bins_perp, binWidth=False)
    
    hist.Add(hist_bkg, -1)

    
    if qTmin < 30: xMin,xMax = -100, 100
    elif qTmin < 60: xMin,xMax = -150, 50
    elif qTmin < 150: xMin,xMax = -200, 50
    else: xMin,xMax = -500, 0
    if "perp" in name: xMin,xMax = -100, 100
    if "_qT" in name: xMin,xMax = -100, 100
    
    hist_norm_data = hist.Integral()
    hist.Scale(1./hist_norm_data)    

 
    recoil = ROOT.RooRealVar("recoil", "Recoil parallel (GeV)", 0, -500, 500)
    cfg['xmin'], cfg['xmax'] = xMin, xMax
    cfg['ymax'] = 1.75*hist.GetMaximum()
    cfg['ymax'] = 1e3
    cfg['ymin'] = 1e-6
    
    


    ## signal
    mean1 = ROOT.RooRealVar("mean1", "", mean1_, -500, 500)
    mean2 = ROOT.RooRealVar("mean2", "", mean2_, -500, 500)
    mean3 = ROOT.RooRealVar("mean3", "", mean3_, -500, 500)
    sigma1 = ROOT.RooRealVar("sigma1", "", sigma1_, 0.1, 100)
    sigma2 = ROOT.RooRealVar("sigma2", "", sigma2_, 0.1, 100)
    sigma3 = ROOT.RooRealVar("sigma3", "", sigma3_, 0.1, 100)
    norm1 = ROOT.RooRealVar("norm1", "", norm1_, 0, 1)
    norm2 = ROOT.RooRealVar("norm2", "", norm2_, 0, 1)
    gauss1 = ROOT.RooGaussian("gauss1", "gauss1", recoil, mean1, sigma1)
    gauss2 = ROOT.RooGaussian("gauss2", "gauss2", recoil, mean2, sigma2)
    gauss3 = ROOT.RooGaussian("gauss3", "gauss3", recoil, mean3, sigma3)
    pdf = ROOT.RooAddPdf("pdf_sig", '', ROOT.RooArgList(gauss1, gauss2, gauss3), ROOT.RooArgList(norm1, norm2))

    rdh = ROOT.RooDataHist("rdh", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(hist))
    

    if "perp" in name:
        #mean1.setVal(0)
        #mean1.setConstant(ROOT.kTRUE)
        #mean2.setVal(0)
        #mean2.setConstant(ROOT.kTRUE)
        #mean3.setVal(0)
        #mean3.setConstant(ROOT.kTRUE)
        
        
        norm1.setVal(0.5)
        norm1.setConstant(ROOT.kTRUE)
        
        norm2.setVal(0.25)
        norm2.setConstant(ROOT.kTRUE)

    fitRes = pdf.fitTo(rdh, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save(ROOT.kTRUE))   
       
    #mean2.setVal(mean1.getVal())
    #mean3.setVal(mean1.getVal())    

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=0)
    canvas.cd()
    padT.Draw()
    padT.cd()
    dummyT.Draw("HIST")
    
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    plt = recoil.frame()
    rdh.plotOn(plt)
    gauss1.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(norm1.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    gauss2.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Normalization(norm2.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    gauss3.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Normalization(1.-norm1.getVal()-norm2.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
    pdf.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
    histpull = plt.pullHist()
    chi2 = plt.chiSquare() 
  
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.20, 0.85, label)
    latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin, qTmax))
    latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist.GetMean(), hist.GetMeanError()))
    latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist.GetRMS(), hist.GetRMSError()))

    latex.DrawLatex(0.60, 0.85, "#mu_{1} = %.3f #pm %.3f" % (mean1.getVal(), mean1.getError()))
    latex.DrawLatex(0.60, 0.80, "#mu_{2} = %.3f #pm %.3f" % (mean2.getVal(), mean2.getError()))
    latex.DrawLatex(0.60, 0.75, "#mu_{3} = %.3f #pm %.3f" % (mean3.getVal(), mean3.getError()))
    latex.DrawLatex(0.60, 0.70, "#sigma_{1} = %.3f #pm %.3f" % (sigma1.getVal(), sigma1.getError()))
    latex.DrawLatex(0.60, 0.65, "#sigma_{2} = %.3f #pm %.3f" % (sigma2.getVal(), sigma2.getError()))
    latex.DrawLatex(0.60, 0.60, "#sigma_{3} = %.3f #pm %.3f" % (sigma3.getVal(), sigma3.getError()))
    latex.DrawLatex(0.60, 0.55, "N_{1} = %.3f #pm %.3f" % (norm1.getVal(), norm1.getError()))
    latex.DrawLatex(0.60, 0.50, "N_{2} = %.3f #pm %.3f" % (norm2.getVal(), norm2.getError()))
    latex.DrawLatex(0.60, 0.45, "#chi^{2}/ndof = %.3f" % chi2)

    plt.Draw("SAME")
    plotter.auxRatio()
       
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
    plt = recoil.frame()
    plt.addPlotable(histpull, "P")
    plt.Draw("SAME")
    dummyL.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s/%03d_recoil.png" % (outDir, tag, qTbin))
    #canvas.SaveAs("%s/%03d_recoil_%s.pdf" % (outDir, qTbin, name))

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)   

    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()
    
    hOut.SetBinContent(qTbin, 1, 1, mean1.getVal())
    hOut.SetBinContent(qTbin, 1, 2, mean1.getError())
    hOut.SetBinContent(qTbin, 2, 1, sigma1.getVal())
    hOut.SetBinContent(qTbin, 2, 2, sigma1.getError())
    hOut.SetBinContent(qTbin, 3, 1, norm1.getVal())
    hOut.SetBinContent(qTbin, 3, 2, norm1.getError())
    
    hOut.SetBinContent(qTbin, 4, 1, mean2.getVal())
    hOut.SetBinContent(qTbin, 4, 2, mean2.getError())
    hOut.SetBinContent(qTbin, 5, 1, sigma2.getVal())
    hOut.SetBinContent(qTbin, 5, 2, sigma2.getError())
    hOut.SetBinContent(qTbin, 6, 1, norm2.getVal())
    hOut.SetBinContent(qTbin, 6, 2, norm2.getError())
    
    hOut.SetBinContent(qTbin, 7, 1, mean3.getVal())
    hOut.SetBinContent(qTbin, 7, 2, mean3.getError())
    hOut.SetBinContent(qTbin, 8, 1, sigma3.getVal())
    hOut.SetBinContent(qTbin, 8, 2, sigma3.getError())
    hOut.SetBinContent(qTbin, 9, 1, (1. - norm1.getVal() - norm2.getVal()))
    hOut.SetBinContent(qTbin, 9, 2, 0)
    
    return mean1.getVal(), mean2.getVal(), mean3.getVal(), norm1.getVal(), norm2.getVal(), sigma1.getVal(), sigma2.getVal(), sigma3.getVal()
 
 

if __name__ == "__main__":

    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrection/fits_Z_auto_RoccoR_lowPU_v0_binned/"
    fIn = ROOT.TFile("/eos/user/j/jaeyserm/analysis/lowPU/recoil_lowPU_RoccoR_lowPU_v0.root")


    cfg = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : -300,
        'xmax'              : 100,
        'ymin'              : 1e-1,
        'ymax'              : 1,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)",
        'ytitle'            : "Events / 1 GeV",
        
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Pull",
        'yminR'             : -3.8,
        'ymaxR'             : 3.8,
    }   
    

   
    recoil_qTbins = list(range(0, 30, 1)) + list(range(30, 50, 2)) + list(range(50, 100, 5)) + list(range(100, 150, 10)) + [150, 175, 200, 10000]
    
    bins_perp = [-200, -100, -75, -50, -40, -35, -30, -25, -20, -15, -14, -13, - -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100, 200]
    bins_perp = range(-200, 202, 2)

    # PERP MC
    h_perp_mc = ROOT.TH3D("h_perp_mc", "", 100, 0, 100, 10, 0, 10, 10, 0, 10) # qt, nominal, error
    ##mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.4, 0.4, 8, 12, 6
    ##mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.5, 0.1, 8.5, 10, 6
    ##mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.5, 0.2, 8, 12, 6
    mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.5, 0.4, 8, 5, 18
    for i in range(1, len(recoil_qTbins)):
        continue
        mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = doFitMC_3Gauss("mc_perp", "mc_perp_bin%d" % i, h_perp_mc, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__)
    

    # PARA MC
    h_para_mc = ROOT.TH3D("h_para_mc", "", 100, 0, 100, 10, 0, 10, 10, 0, 10) # qt, nominal, error
    ##mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = -4.58011e-01, -4.45388e-01, -3.04321e-01, 5.97111e-01, 2.05346e-01, 7.58247e+00, 4.00142e+00, 1.32090e+01
    mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0.2, 0.2, 0.2, 0.5, 0.4, 8, 5, 18
    for i in range(1, len(recoil_qTbins)):
        continue
        mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = doFitMC_3Gauss("mc_para", "mc_para_qT_bin%d" % i, h_para_mc, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__)
    
    # PARA DATA    
    h_para_data = ROOT.TH3D("h_para_data", "", 100, 0, 100, 10, 0, 10, 10, 0, 10) # qt, nominal, error
    ##mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg = 0, 0, 0, 0.3, 0.7, 3, 15, 8
    ##mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 4.21743e-01, 4.01514e-02, 5.61924e+00, 1.78429e+01, 9.60333e+00
    mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg = 1.189, 1.416, 1.081, 0.3, 0.7, 5.810, 32.229, 7.987
    mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.5, 0.4, 10, 5, 18
    for i in range(1, len(recoil_qTbins)):
        continue
        mean1__bkg, mean2__bkg, mean3__bkg,  norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg, hist_norm_bkg = doFitBkg_3Gauss("bkg_para", "bkg_para_bin%d" % i, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg)
        mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = doFitData_3Gauss("data_para", "data_para_qT_bin%d" % i, h_para_data, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__, mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg, hist_norm_bkg)

    #sys.exit()
    # PERP DATA    
    h_perp_data = ROOT.TH3D("h_perp_data", "", 100, 0, 100, 10, 0, 10, 10, 0, 10) # qt, nominal, error
    mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg = 0, 0, 0, 0.3, 0.7, 3, 15, 8
    #mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 4.21743e-01, 4.01514e-02, 5.61924e+00, 1.78429e+01, 9.60333e+00
    mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg = 0, 0, 0, 0.3, 0.0, 3, 15, 8
    mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = 0, 0, 0, 0.5, 0.2, 10, 10, 20
    for i in range(1, len(recoil_qTbins)):
        #continue
        mean1__bkg, mean2__bkg, mean3__bkg,  norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg, hist_norm_bkg = doFitBkg_3Gauss("bkg_perp", "bkg_perp_bin%d" % i, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg)
        mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__ = doFitData_3Gauss("data_perp", "data_perp_bin%d" % i, h_perp_data, i, recoil_qTbins[i-1], recoil_qTbins[i], mean1__, mean2__, mean3__, norm1__, norm2__, sigma1__, sigma2__, sigma3__, mean1__bkg, mean2__bkg, mean3__bkg, norm1__bkg, norm2__bkg, sigma1__bkg, sigma2__bkg, sigma3__bkg, hist_norm_bkg)
        
      
    fOut = ROOT.TFile("/eos/user/j/jaeyserm/analysis/lowPU/recoil_lowPU_fitParam_RoccoR_lowPU_v0_mirror.root", "RECREATE")
    h_para_mc.Write()
    h_perp_mc.Write()
    h_para_data.Write()
    h_perp_data.Write()
    fOut.Close()      
    fIn.Close()