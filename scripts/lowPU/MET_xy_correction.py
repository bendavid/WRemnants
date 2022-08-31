
import sys,array,math,os,copy,json
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import functions
import plotter


from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z
from wremnants.datasets.datagroups import datagroups2016

import lz4.frame
import pickle
import narf
import numpy as np


def readProc(datagroups, hName, procName):

    label = "%s_%s" % (hName, procName)
    datagroups.setHists(hName, "", label=label, procsToRead=[procName], selectSignal=False)
    bhist = datagroups.groups[procName][label]
    return bhist 


def makePlot(hist_data, hist_mc, fOut, xLabel, npv, outDir_):

    hist_data.Scale(1./hist_data.Integral())
    hist_mc.Scale(1./hist_mc.Integral())
    
    hist_data.SetLineColor(ROOT.kBlack)
    hist_data.SetLineWidth(2)
    
    hist_mc.SetLineColor(ROOT.kRed)
    hist_mc.SetLineWidth(2)
    
    ## sigmas
    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : -100,
        'xmax'              : 100,
        'ymin'              : 1e-5,
        'ymax'              : 1e0,
            
        'xtitle'            : xLabel,
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    hist_data.Draw("SAME")
    hist_mc.Draw("SAME")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.030)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    
    latex.DrawLatex(0.20, 0.90, "NPV = %d" % npv)
    
    latex.DrawLatex(0.20, 0.85, "Data")
    latex.DrawLatex(0.20, 0.82, "Mean %.2f #pm %.2f GeV" % (hist_data.GetMean(), hist_data.GetMeanError()))
    
    latex.DrawLatex(0.60, 0.85, "#color[2]{MC}")
    latex.DrawLatex(0.60, 0.82, "#color[2]{Mean %.2f #pm %.2f GeV}" % (hist_mc.GetMean(), hist_mc.GetMeanError()))
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s.png" % (outDir_, fOut))
    canvas.SaveAs("%s/%s.pdf" % (outDir_, fOut))
    canvas.Delete()

def METxyCorrection(direction = "x", corrType="uncorr", polyOrderData=-1, polyOrderMC=-1):

    
    # corrType uncorr or corr_xy
    outDir_ = "%s/%s" % (outDir, corrType)
    functions.prepareDir(outDir_, False)

    
    procs = ["DYmumu"] + bkgs
    
    b_data = readProc(datagroups, "MET%s_%s_npv" % (direction, corrType), "SingleMuon" if "mu" in flavor else "SingleElectron")
    b_mc = None
    for proc in procs:
    
        b = readProc(datagroups, "MET%s_%s_npv" % (direction, corrType), proc)
        if b_mc == None: b_mc = b
        else: b_mc += b

    h_data = narf.hist_to_root(b_data)
    h_mc = narf.hist_to_root(b_mc)
    
    g_data = ROOT.TGraphErrors()
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1)
    g_data.SetMarkerColor(ROOT.kBlack)
    
    g_mc = ROOT.TGraphErrors()
    g_mc.SetLineColor(ROOT.kRed)
    g_mc.SetMarkerStyle(21)
    g_mc.SetMarkerSize(1)
    g_mc.SetMarkerColor(ROOT.kRed)
    
    outDict = {"data": { 'nominal': {} }, "mc": { 'nominal': {} }}
    for iBin in range(1, npv_max):

        hist_data = h_data.ProjectionY("data_bin%d" % iBin, iBin, iBin)
        hist_mc = h_mc.ProjectionY("mc_bin%d" % iBin, iBin, iBin)
        
        mean_data, mean_data_err = hist_data.GetMean(), hist_data.GetMeanError()
        mean_mc, mean_mc_err = hist_mc.GetMean(), hist_mc.GetMeanError()
        print(iBin, mean_mc, mean_mc_err)
        
        makePlot(hist_data, hist_mc, "npv_%d_%d_MET%s" % (iBin, iBin, direction), "MET %s (GeV)" % direction, iBin, outDir_)
        
        g_data.SetPoint(iBin-1, iBin, mean_data)
        g_data.SetPointError(iBin-1, 0, mean_data_err)
        g_mc.SetPoint(iBin-1, iBin, mean_mc)
        g_mc.SetPointError(iBin-1, 0, mean_mc_err)
  
    if polyOrderData > 0:
        fit_data = ROOT.TF1("fit_data_%s" % direction, "pol%d" % polyOrderData, 0, 10)
        result = g_data.Fit(fit_data.GetName(), "NSE", "", 0, 10) 
        fit_data.SetLineColor(ROOT.kBlack)
        fit_data.GetXaxis().SetRangeUser(0, 10)
        fit_data.SetLineWidth(2)
        for iP in range(0, polyOrderData+1): outDict['data']['nominal']['p%d' % iP] = fit_data.GetParameter(iP)
        outDict['data']['nominal']['polyOrder'] = polyOrderData
    
    if polyOrderMC > 0:
        fit_mc = ROOT.TF1("fit_mc_%s" % direction, "pol%d" % polyOrderMC, 0, npv_fit_max)
        result = g_mc.Fit(fit_mc.GetName(), "NSE", "", 0, npv_fit_max) 
        fit_mc.SetLineColor(ROOT.kRed)
        fit_mc.SetLineStyle(ROOT.kDashed)
        fit_mc.GetXaxis().SetRangeUser(0, npv_max)
        fit_mc.SetLineWidth(2)
        for iP in range(0, polyOrderMC+1): outDict['mc']['nominal']['p%d' % iP] = fit_mc.GetParameter(iP)
        outDict['mc']['nominal']['polyOrder'] = polyOrderMC
  
    ## sigmas
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : npv_max,
        'ymin'              : -5,
        'ymax'              : 5,
            
        'xtitle'            : "NPV",
        'ytitle'            : "#LT MET_{%s} #GT (Gev)" % direction,
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    g_data.Draw("PE SAME")
    if polyOrderData > 0: fit_data.Draw("L SAME")
    g_mc.Draw("PE SAME")
    if polyOrderMC > 0: fit_mc.Draw("L SAME")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    
    if flavor == "mumu": label = "Z #rightarrow #mu^{#plus}#mu^{#minus}"
    elif flavor == "ee": label = "Z #rightarrow e^{#plus}e^{#minus}"
    elif flavor == "mu": label = "W #rightarrow #mu"
    elif flavor == "e": label = "W #rightarrow e"
    else: label = ""
    latex.DrawLatex(0.20, 0.90, label)
    latex.DrawLatex(0.20, 0.85, "#LT MET_{%s} #GT, %s" % (direction, met))

    leg = ROOT.TLegend(.2, 0.7, .5, .80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(g_data, "Data", "LP")
    leg.AddEntry(g_mc, "MC", "LP")
    leg.Draw("SAME")
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/MET_%s_npv.png" % (outDir_, direction))
    canvas.SaveAs("%s/MET_%s_npv.pdf" % (outDir_, direction))
    canvas.Delete()
    
    return outDict
    
    
if __name__ == "__main__":

    met = "RawPFMET" # RawPFMET DeepMETReso
    flavor = "ee" # mu, e, mumu, ee
    npv_max, npv_fit_max = 10, 10

    # DATA For electron channels!
    
    ####################################################################
    datagroups = datagroupsLowPU_Z("lowPU_%s_%s.pkl.lz4" % (flavor, met), flavor=flavor)
    bkgs = ["DYmumu", "DYee", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]

    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/METxy_correction/METxy_%s_%s/" % (flavor, met)
    fOut = "wremnants/data/lowPU/MET_xy_corr_coeff_%s_%s.json" % (flavor, met)
    functions.prepareDir(outDir, True)
    

  
    dictout = {}
    dictX = METxyCorrection(direction="x", corrType="corr_lep", polyOrderData=1, polyOrderMC=1)
    dictY = METxyCorrection(direction="y", corrType="corr_lep", polyOrderData=1, polyOrderMC=1)
    dictout['x'] = dictX
    dictout['y'] = dictY
    jsOut = json.dumps(dictout, indent = 4)
    with open(fOut, "w") as outfile: outfile.write(jsOut)
    os.system("cp %s %s" % (fOut, outDir)) # make copy to web dir

    
    METxyCorrection(direction="x", corrType="corr_xy", polyOrderData=1, polyOrderMC=1)
    METxyCorrection(direction="y", corrType="corr_xy", polyOrderData=1, polyOrderMC=1)