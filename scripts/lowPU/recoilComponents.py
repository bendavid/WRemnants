
import sys,array,math,os,copy,json
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import functions
import plotter


from wremnants.datasets.datagroupsLowPU import make_datagroups_lowPU

import lz4.frame
import pickle
import narf
import numpy as np


def readProc(datagroups, hName, procName):

    label = "%s_%s" % (hName, procName)
    datagroups.setHists(hName, "", label=label, procsToRead=[procName], selectSignal=False)
    bhist = datagroups.groups[procName][label]
    return bhist 



def makePlot(tag="para_qT"):

    corr = "corr_xy" # corr_xy uncorr corr_lep

    b_data = readProc(datagroups, "recoil_%s_%s_qTbinned" % (corr, tag), data)
    b_sig = readProc(datagroups, "recoil_%s_%s_qTbinned" % (corr, tag), sig)
    b_bkgs = None
    for bkg in bkgs:
    
        b = readProc(datagroups, "recoil_%s_%s_qTbinned" % (corr, tag), bkg)
        if b_bkgs == None: b_bkgs = b
        else: b_bkgs += b

    h_data = narf.hist_to_root(b_data)
    h_sig = narf.hist_to_root(b_sig)
    h_bkgs = narf.hist_to_root(b_bkgs)
    
    rebin=1

    g_data_mean = ROOT.TGraphErrors()
    g_data_mean.SetLineColor(ROOT.kBlack)
    g_data_mean.SetMarkerStyle(20)
    g_data_mean.SetMarkerSize(1)
    g_data_mean.SetMarkerColor(ROOT.kBlack)
    
    g_sig_mean = ROOT.TGraphErrors()
    g_sig_mean.SetLineColor(ROOT.kRed)
    g_sig_mean.SetMarkerStyle(21)
    g_sig_mean.SetMarkerSize(1)
    g_sig_mean.SetMarkerColor(ROOT.kRed)
    
    g_data_bkgs_mean = ROOT.TGraphErrors()
    g_data_bkgs_mean.SetLineColor(ROOT.kBlue)
    g_data_bkgs_mean.SetMarkerStyle(21)
    g_data_bkgs_mean.SetMarkerSize(1)
    g_data_bkgs_mean.SetMarkerColor(ROOT.kBlue)
    
    g_data_sigma = ROOT.TGraphErrors()
    g_data_sigma.SetLineColor(ROOT.kBlack)
    g_data_sigma.SetMarkerStyle(20)
    g_data_sigma.SetMarkerSize(1)
    g_data_sigma.SetMarkerColor(ROOT.kBlack)
    
    g_sig_sigma = ROOT.TGraphErrors()
    g_sig_sigma.SetLineColor(ROOT.kRed)
    g_sig_sigma.SetMarkerStyle(21)
    g_sig_sigma.SetMarkerSize(1)
    g_sig_sigma.SetMarkerColor(ROOT.kRed)
    
    g_data_bkgs_sigma = ROOT.TGraphErrors()
    g_data_bkgs_sigma.SetLineColor(ROOT.kBlue)
    g_data_bkgs_sigma.SetMarkerStyle(21)
    g_data_bkgs_sigma.SetMarkerSize(1)
    g_data_bkgs_sigma.SetMarkerColor(ROOT.kBlue)
    
    g_data_response = ROOT.TGraphErrors()
    g_data_response.SetLineColor(ROOT.kBlack)
    g_data_response.SetMarkerStyle(20)
    g_data_response.SetMarkerSize(1)
    g_data_response.SetMarkerColor(ROOT.kBlack)
    
    g_sig_response = ROOT.TGraphErrors()
    g_sig_response.SetLineColor(ROOT.kRed)
    g_sig_response.SetMarkerStyle(21)
    g_sig_response.SetMarkerSize(1)
    g_sig_response.SetMarkerColor(ROOT.kRed)
    
    g_data_bkgs_response = ROOT.TGraphErrors()
    g_data_bkgs_response.SetLineColor(ROOT.kBlue)
    g_data_bkgs_response.SetMarkerStyle(21)
    g_data_bkgs_response.SetMarkerSize(1)
    g_data_bkgs_response.SetMarkerColor(ROOT.kBlue)  
  
    for iBin in range(1, len(recoil_qTbins)):
    
        qT = 0.5*(recoil_qTbins[iBin-1]+recoil_qTbins[iBin])
        hist_data = h_data.ProjectionY("data_bin%d" % iBin, iBin, iBin)
        hist_sig = h_sig.ProjectionY("sig_bin%d" % iBin, iBin, iBin)
        hist_bkgs = h_bkgs.ProjectionY("data_bkgs_bin%d" % iBin, iBin, iBin)
        
        hist_data.Rebin(rebin)
        hist_sig.Rebin(rebin)
        hist_bkgs.Rebin(rebin)
        
        mean_data, mean_data_err = hist_data.GetMean(), hist_data.GetMeanError()
        sigma_data, sigma_data_err = hist_data.GetRMS(), hist_data.GetRMSError()
        g_data_mean.SetPoint(iBin-1, qT, mean_data)
        g_data_mean.SetPointError(iBin-1, 0, mean_data_err)
        g_data_sigma.SetPoint(iBin-1, qT, sigma_data)
        g_data_sigma.SetPointError(iBin-1, 0, sigma_data_err)

        mean_sig, mean_sig_err = hist_sig.GetMean(), hist_sig.GetMeanError()
        sigma_sig, sigma_sig_err = hist_sig.GetRMS(), hist_sig.GetRMSError()
        g_sig_mean.SetPoint(iBin-1, qT, mean_sig)
        g_sig_mean.SetPointError(iBin-1, 0, mean_sig_err)
        g_sig_sigma.SetPoint(iBin-1, qT, sigma_sig)
        g_sig_sigma.SetPointError(iBin-1, 0, sigma_sig_err)
        
        hist_data.Add(hist_bkgs, -1)
        mean_data_bkgs, mean_data_bkgs_err = hist_data.GetMean(), hist_data.GetMeanError()
        sigma_data_bkgs, sigma_data_bkgs_err = hist_data.GetRMS(), hist_data.GetRMSError()
        g_data_bkgs_mean.SetPoint(iBin-1, qT, mean_data_bkgs)
        g_data_bkgs_mean.SetPointError(iBin-1, 0, mean_data_bkgs_err)
        g_data_bkgs_sigma.SetPoint(iBin-1, qT, sigma_data_bkgs)
        g_data_bkgs_sigma.SetPointError(iBin-1, 0, sigma_data_bkgs_err)
        
        

        g_data_response.SetPoint(iBin-1, qT, mean_data/qT)
        g_data_response.SetPointError(iBin-1, 0, mean_data_err/qT)
        g_sig_response.SetPoint(iBin-1, qT, mean_sig/qT)
        g_sig_response.SetPointError(iBin-1, 0, mean_sig_err/qT)
  
  
    if flavor == "mumu": label = "Z #rightarrow #mu^{#plus}#mu^{#minus}"
    elif flavor == "ee": label = "Z #rightarrow e^{#plus}e^{#minus}"
    elif flavor == "mu": label = "W #rightarrow #mu"
    elif flavor == "e": label = "W #rightarrow e"
    else: label = ""
    
    if tag == "para_qT": tag_ = "U_{#parallel} + q_{T}"
    elif tag == "para": tag_ = "U_{#parallel}"
    elif tag == "perp": tag_ = "U_{#perp}  "

    yMin, yMax = -5, 15
    if tag == "para": yMin, yMax = -100, 50

    ## mean
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : 100,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "#LT #mu(%s) #GT" % tag_,
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    
    line = ROOT.TLine(0, 0, 100, 0)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineWidth(3)
    line.SetLineStyle(1)
    line.Draw("SAME")
    
    g_data_mean.Draw("PE SAME")
    #g_data_bkgs_mean.Draw("PE SAME")
    g_sig_mean.Draw("PE SAME")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.90, label)
    latex.DrawLatex(0.20, 0.85, "#LT #mu(%s) #GT, %s" % (tag_, met))

    leg = ROOT.TLegend(.2, 0.7, .5, .80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(g_data_mean, "Data", "LP")
    #leg.AddEntry(g_data_bkgs, "Data #minus bkgs.", "LP")
    leg.AddEntry(g_sig_mean, "MC", "LP")
    leg.Draw("SAME")
    
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s_mean.png" % (outDir, tag))
    canvas.SaveAs("%s/%s_mean.pdf" % (outDir, tag))
    canvas.Delete()

    
    
    ## sigma
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : 100,
        'ymin'              : 5,
        'ymax'              : 20,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "#LT #sigma(%s) #GT" % tag_,
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    g_data_sigma.Draw("PE SAME")
    #g_data_bkgs_sigma.Draw("PE SAME")
    g_sig_sigma.Draw("PE SAME")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.90, label)
    latex.DrawLatex(0.20, 0.85, "#LT #sigma(%s) #GT, %s" % (tag_, met))

    leg = ROOT.TLegend(.2, 0.7, .5, .80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(g_data_sigma, "Data", "LP")
    #leg.AddEntry(g_data_bkgs_sigma, "Data #minus bkgs.", "LP")
    leg.AddEntry(g_sig_sigma, "MC", "LP")
    leg.Draw("SAME")
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s_sigma.png" % (outDir, tag))
    canvas.SaveAs("%s/%s_sigma.pdf" % (outDir, tag))
    canvas.Delete()
    
    
    ## response
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : 100,
        'ymin'              : 0,
        'ymax'              : 2,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "#LT #mu(%s) #GT" % tag_,
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    
    line = ROOT.TLine(0, 0, 100, 0)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineWidth(3)
    line.SetLineStyle(1)
    line.Draw("SAME")
    
    g_data_response.Draw("PE SAME")
    #g_data_bkgs_mean.Draw("PE SAME")
    g_sig_response.Draw("PE SAME")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.90, label)
    latex.DrawLatex(0.20, 0.85, "#LT #mu(%s) #GT, %s" % (tag_, met))

    leg = ROOT.TLegend(.2, 0.7, .5, .80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(g_data_response, "Data", "LP")
    #leg.AddEntry(g_data_bkgs, "Data #minus bkgs.", "LP")
    leg.AddEntry(g_sig_response, "MC", "LP")
    leg.Draw("SAME")
    
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s_response.png" % (outDir, tag))
    canvas.SaveAs("%s/%s_response.pdf" % (outDir, tag))
    canvas.Delete()

    
    
if __name__ == "__main__":

    met = "RawPFMET" # RawPFMET DeepMETReso
    flavor = "mumu" # mu, e, mumu, ee

    recoil_qTbins = list(range(0, 50, 1)) + list(range(50, 70, 2)) + list(range(70, 100, 5)) + list(range(100, 150, 10)) + [150, 200, 300, 10000]
    
    ####################################################################
    datagroups = make_datagroups_lowPU("lowPU_%s_%s.pkl.lz4" % (flavor, met), flavor=flavor)
    data, sig, bkgs = "SingleMuon", "DYmumu", ['EWK', 'TTbar'] 
    #bkgs = ["DYee", "DYtautau", "TTTo2L2Nu", "TTToSemiLeptonic", "ZZ", "WZTo3LNu", "WWTo2L2Nu", "WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToENu", "WminusJetsToENu", "WplusJetsToTauNu", "WminusJetsToTauNu"]
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrection/RecoilComponents_%s_%s/" % (flavor, met)
    functions.prepareDir(outDir, True)
    

  
    makePlot(tag="para_qT")
    makePlot(tag="para")
    makePlot(tag="perp")