
import sys,array,math,os,copy,fnmatch
from collections import OrderedDict
import ctypes

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import lz4.frame
import pickle
import narf
import hist
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z

import plotter
import functions

import scripts.lowPU.config as lowPUcfg




def parseHist(fIn, hName, bins):

    h = ROOT.TH1D(hName, "", len(bins)-1, array.array('d', bins))
    
    hIn = fIn.Get(hName)
    for iReco in range(0, len(bins)-1):
    
        
        err = ctypes.c_double()
        integral = hIn.IntegralAndError(iReco*mZbins+1, (iReco+1)*mZbins, err)
        
        h.SetBinContent(iReco+1, integral)
        h.SetBinError(iReco+1, err)
        #print(iReco, integral, err)

    return h
    
    
def parseHist1D(fIn, hName, bins):

    h = ROOT.TH1D(hName, "", len(bins)-1, array.array('d', bins))
    
    hIn = fIn.Get(hName)
    for iReco in range(0, len(bins)-1):
    
        
        err = ctypes.c_double()
        integral = hIn.IntegralAndError(iReco*mZbins+1, (iReco+1)*mZbins, err)
        
        h.SetBinContent(iReco+1, integral)
        h.SetBinError(iReco+1, err)
        #print(iReco, integral, err)

    return h
    
def doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.06):

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_output.root" % fitcfg)
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Combine/Z/fit_%s" % fitcfg
    functions.prepareDir(outDir, remove=False)
    
    resTree = fIn.Get("fitresults")
    
    if flavor == "mumu":
        dataLabel = "Data (#mu^{#plus}#mu^{#minus})"
        signalLabel = "DY#rightarrow#mu^{#plus}#mu^{#minus}"
        
    if flavor == "ee":
        dataLabel = "Data (e^{#plus}e^{#minus})"
        signalLabel = "DY#rightarrowe^{#plus}e^{#minus}"
        
    if fitcfg == "mumu": label = "Muon channel"
    elif fitcfg == "ee": label = "Electron channel"
    else: label = "Muon + electron channels"

    
    
    sigUncs = [0.077181, 0.171347, 0.100207, 0.119719, 0.115287, 0.108075]
    sigUncs = [0.035796, 0.072002, 0.050658, 0.059626, 0.057034, 0.060837]
    sigUncs = [0.041657, 0.052776, 0.039676, 0.050306, 0.052252, 0.065191]
    


    
    leg = ROOT.TLegend(.17, 0.57, .9, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.031)
    leg.SetNColumns(2)
    leg.SetMargin(0.7*leg.GetMargin())
    
    
    # get data
    h_data = parseHist(fIn, "obs;1", bins_recoil_reco)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)       
    h_data.Scale(1, "width")

    
    st = ROOT.THStack()
    st.SetName("stack")
    
    # backgrounds
    h_top = parseHist(fIn, "expproc_TTbar_%s;1" % fitmode, bins_recoil_reco)
    h_top.SetFillColor(ROOT.TColor.GetColor(222, 90, 106))
    h_top.SetLineColor(ROOT.kBlack)
    h_top.SetLineWidth(1)
    h_top.SetLineStyle(1)
    h_top.Scale(1, "width")
    st.Add(h_top)
    
        
    h_ewk = parseHist(fIn, "expproc_EWK_%s;1" % fitmode, bins_recoil_reco)
    h_ewk.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
    h_ewk.SetLineColor(ROOT.kBlack)
    h_ewk.SetLineWidth(1)
    h_ewk.SetLineStyle(1)
    h_ewk.Scale(1, "width")
    st.Add(h_ewk)
    
    
    # sum of backgrounds
    h_tot = parseHist(fIn, "expfull_%s;1" % fitmode, bins_recoil_reco)
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetFillColor(0)
    h_tot.SetLineWidth(2)
    h_tot.Scale(1, "width")



    # ratio
    h_ratio = h_data.Clone("ratio")
    h_tot_noerr = h_tot.Clone("h_tot_noerr")
    for i in range(0, h_tot_noerr.GetNbinsX()+1): h_tot_noerr.SetBinError(i, 0)
    h_ratio.Divide(h_tot_noerr)
    
    # ratio err
    h_tot_err = None
    for iGenBin in range(1, len(lowPUcfg.bins_recoil_gen)):

        h_err = parseHist(fIn, "expproc_DY_genBin%d_%s;1" % (iGenBin, fitmode), bins_recoil_reco)
        
        # append the signal strength uncertainties
        if fitmode == "postfit":
            for iBin in range(0, h_err.GetNbinsX()+1):
            
                nom = h_err.GetBinContent(iBin)
                err = h_err.GetBinError(iBin)
                err_sig = nom*sigUncs[iGenBin-1]
                err_new = math.sqrt(err**2 + err_sig**2)
                h_err.SetBinError(iBin, err_new)
                #if iGenBin != 2: h_err.SetBinError(iBin, 0)
                print(iBin, nom, err, err_sig, err_new)

        if h_tot_err == None: h_tot_err = h_err
        else: h_tot_err.Add(h_err)
        
     
    
    #h_err = parseHist("expproc_Top_%s;1" % mode)
    #h_tot_err.Add(h_err)
    #h_err = parseHist("expproc_EWK_%s;1" % mode)
    #h_tot_err.Add(h_err)
    
    h_tot_err.Scale(1, "width")
    h_tot_err.SetFillColor(ROOT.kBlack)
    h_tot_err.SetMarkerSize(0)
    h_tot_err.SetLineWidth(0)
    h_tot_err.SetFillStyle(3004)
    
    for i in range(1, 10):
        print(h_tot_err.GetBinContent(i), h_tot_err.GetBinError(i))
    
    '''
    h_tot_err = parseHist("expfull_%s;1" % mode)
    h_tot_err.Scale(1, "width")
    h_tot_err.SetFillColor(ROOT.kBlack)
    h_tot_err.SetMarkerSize(0)
    h_tot_err.SetLineWidth(0)
    h_tot_err.SetFillStyle(3004)
    '''

    h_tot_err_ratio = h_tot_err.Clone("ratio_err")
    h_tot_err_ratio_denom = h_tot_err.Clone("ratio_err_denom")
    #for i in range(0, h_tot_err_ratio_denom.GetNbinsX()+1): h_tot_err_ratio_denom.SetBinError(i, 0)
    h_tot_err_ratio.Divide(h_tot_err_ratio_denom)
    
    
    leg.AddEntry(h_data, dataLabel, "PE")
    leg.AddEntry(h_tot_err, "Stat. + Syst. Unc. (%s)" % fitmode, "F")
    leg.AddEntry(h_top, "Top", "F")
    leg.AddEntry(h_ewk, "EWK (#tau#tau, diboson)", "F")
    
    
    

    

    
    # gen signals
    sigColors = [ROOT.kOrange, ROOT.kOrange+1, ROOT.kOrange+2, ROOT.kOrange+4, ROOT.kOrange-1, ROOT.kOrange+5]
    for iGenBin in range(1, len(lowPUcfg.bins_recoil_gen)):

        h_gen_sig = parseHist(fIn, "expproc_DY_genBin%d_%s;1" % (iGenBin, fitmode), bins_recoil_reco)
        h_gen_sig.SetFillColor(sigColors[iGenBin-1])
        h_gen_sig.SetLineColor(ROOT.kBlack)
        h_gen_sig.SetLineWidth(1)
        h_gen_sig.SetLineStyle(1)
        h_gen_sig.Scale(1, "width")

        
        leg.AddEntry(h_gen_sig, "%s, p_{T}^{Z, gen} #in [%d, %d] GeV" % (signalLabel, lowPUcfg.bins_recoil_gen[iGenBin-1], lowPUcfg.bins_recoil_gen[iGenBin]), "F")
        st.Add(h_gen_sig)
          
          


    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : 150,
        'ymin'              : 1e0,
        'ymax'              : 1e7, # 3e6
            
        'xtitle'            : "p_{T}^{Z} (GeV)",
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Data/MC",
        
        'yminR'             : 2.-ratio,
        'ymaxR'             : ratio
    }
    
    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## top panel
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    dummyT.Draw("HIST")
        
    st.Draw("HIST SAME")
    h_tot_err.Draw("SAME E2")
    h_tot.Draw("HIST SAME")
    h_data.Draw("PE SAME")
    leg.Draw("SAME")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.86, label)
    
    plotter.auxRatio()  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()  



    ## bottom panel
    canvas.cd()
    padB.Draw()
    padB.SetFillStyle(0)
    padB.cd()
    padB.SetGrid()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")
   
    h_ratio.Draw("P SAME")
    h_tot_err_ratio.Draw("SAME E2")

    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/%s_qT_reco.png" % (outDir, fitmode))
    canvas.SaveAs("%s/%s_qT_reco.pdf" % (outDir, fitmode))
    canvas.Close()    
 

def doTransverseMassPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.06):

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_wlike_combineOutput.root" % fitcfg)
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Combine/Z/Wlike_%s" % fitcfg
    functions.prepareDir(outDir, remove=False)
    
    resTree = fIn.Get("fitresults")
    
    if flavor == "mumu":
        dataLabel = "Data (#mu^{#plus}#mu^{#minus})"
        signalLabel = "DY#rightarrow#mu^{#plus}#mu^{#minus}"
        
    if flavor == "ee":
        dataLabel = "Data (e^{#plus}e^{#minus})"
        signalLabel = "DY#rightarrowe^{#plus}e^{#minus}"
        
    if fitcfg == "mumu": label = "W#minuslike, muon channel"
    elif fitcfg == "ee": label = "Electron channel"
    else: label = "Muon + electron channels"

    mT_bins = [0, 10, 15, 20, 25, 30, 35] + list(range(40, 100, 2)) + [100, 102, 104, 106, 108, 110, 115, 120, 125, 130, 140, 160, 200]
    ratio=1.08

    leg = ROOT.TLegend(.17, 0.65, .9, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetNColumns(2)
    leg.SetMargin(0.7*leg.GetMargin())
    
    
    # get data
    h_data = fIn.Get("obs;1")
    h_data = functions.Rebin(h_data, mT_bins)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)       

    st = ROOT.THStack()
    st.SetName("stack")
    
    # backgrounds
    h_ewk = fIn.Get("expproc_EWK_%s;1" % fitmode)
    h_ewk = functions.Rebin(h_ewk, mT_bins)
    h_ewk.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
    h_ewk.SetLineColor(ROOT.kBlack)
    h_ewk.SetLineWidth(1)
    h_ewk.SetLineStyle(1)
    st.Add(h_ewk)
    
    h_top = fIn.Get("expproc_TTbar_%s;1" % fitmode)
    h_top = functions.Rebin(h_top, mT_bins)
    h_top.SetFillColor(ROOT.TColor.GetColor(222, 90, 106))
    h_top.SetLineColor(ROOT.kBlack)
    h_top.SetLineWidth(1)
    h_top.SetLineStyle(1)
    st.Add(h_top)
    
    h_dy = fIn.Get("expproc_DYmumu_%s;1" % fitmode)
    h_dy = functions.Rebin(h_dy, mT_bins)
    h_dy.SetFillColor(ROOT.TColor.GetColor(248, 206, 104))
    h_dy.SetLineColor(ROOT.kBlack)
    h_dy.SetLineWidth(1)
    h_dy.SetLineStyle(1)
    h_dy.Scale(1.026)
    st.Add(h_dy)
    
    # sum of backgrounds
    h_tot = fIn.Get("expfull_%s;1" % fitmode)
    h_tot = functions.Rebin(h_tot, mT_bins)
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetFillColor(0)
    h_tot.SetLineWidth(2)
    h_tot.Scale(1.026)


    # ratio
    h_ratio = h_data.Clone("ratio")
    h_tot_noerr = h_tot.Clone("h_tot_noerr")
    for i in range(0, h_tot_noerr.GetNbinsX()+1): h_tot_noerr.SetBinError(i, 0)
    h_ratio.Divide(h_tot_noerr)
    
    # ratio err
    h_tot_err = fIn.Get("expbkg_%s;1" % fitmode)
    h_tot_err = functions.Rebin(h_tot_err, mT_bins)
    h_tot_err.SetFillColor(ROOT.kBlack)
    h_tot_err.SetMarkerSize(0)
    h_tot_err.SetLineWidth(0)
    h_tot_err.SetFillStyle(3004)
    

    h_tot_err_ratio = h_tot_err.Clone("ratio_err")
    h_tot_err_ratio_denom = h_tot_err.Clone("ratio_err_denom")
    #for i in range(0, h_tot_err_ratio_denom.GetNbinsX()+1): h_tot_err_ratio_denom.SetBinError(i, 0)
    h_tot_err_ratio.Divide(h_tot_err_ratio_denom)
    
    
    leg.AddEntry(h_data, dataLabel, "PE")
    leg.AddEntry(h_tot_err, "Stat. + Syst. Unc. (%s)" % fitmode, "F")
    leg.AddEntry(h_top, "Top", "F")
    leg.AddEntry(h_ewk, "EWK (#tau#tau, diboson)", "F")
    leg.AddEntry(h_dy, "DY #rightarrow #mu^{#plus}#mu^{#minus}", "F")

    ## mass variations
    if True:
    
        groups = datagroupsLowPU_Z("mz_lowPU_%s.pkl.lz4" % flavor)
        bhist = groups.readProc("mt_massWeight", "DYmumu")
        hist_nom = narf.hist_to_root(bhist[:, 10])
        hist_plus_1 = narf.hist_to_root(bhist[:, 20])
        hist_minus_1 = narf.hist_to_root(bhist[:, 0])
        hist_plus_2 = narf.hist_to_root(bhist[:, 15])
        hist_minus_2 = narf.hist_to_root(bhist[:, 5])
        
        hist_nom = functions.Rebin(hist_nom, mT_bins)
        hist_plus_1 = functions.Rebin(hist_plus_1, mT_bins)
        hist_minus_1 = functions.Rebin(hist_minus_1, mT_bins)
        hist_plus_2 = functions.Rebin(hist_plus_2, mT_bins)
        hist_minus_2 = functions.Rebin(hist_minus_2, mT_bins)

          
        
        hist_plus_1.SetLineColor(ROOT.kRed)
        hist_plus_1.SetFillColor(0)
        hist_plus_1.SetLineWidth(2)
        hist_minus_1.SetLineColor(ROOT.kRed)
        hist_minus_1.SetFillColor(0)
        hist_minus_1.SetLineWidth(2)
        
        hist_plus_2.SetLineColor(ROOT.kBlue)
        hist_plus_2.SetFillColor(0)
        hist_plus_2.SetLineWidth(2)
        hist_minus_2.SetLineColor(ROOT.kBlue)
        hist_minus_2.SetFillColor(0)
        hist_minus_2.SetLineWidth(2)
        
        hist_plus_err_1 = hist_plus_1.Clone("hist_plus_err_1")
        hist_plus_err_1.Divide(hist_nom)
        hist_minus_err_1 = hist_minus_1.Clone("hist_minus_err_1")
        hist_minus_err_1.Divide(hist_nom)
        
        hist_plus_err_2 = hist_plus_2.Clone("hist_plus_err_2")
        hist_plus_err_2.Divide(hist_nom)
        hist_minus_err_2 = hist_minus_2.Clone("hist_minus_err_2")
        hist_minus_err_2.Divide(hist_nom)
        
        
        leg.AddEntry(hist_plus_1, "m_{Z} #pm 100 MeV", "L")
    


    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : 200,
        'ymin'              : 1e-1,
        'ymax'              : 1e6, # 3e6
            
        'xtitle'            : "m_{T} (GeV)",
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Data/MC",
        
        'yminR'             : 2.-ratio,
        'ymaxR'             : ratio
    }
    
    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## top panel
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    dummyT.Draw("HIST")
        
    st.Draw("HIST SAME")
    h_tot_err.Draw("SAME E2")
    h_tot.Draw("HIST SAME")
    h_data.Draw("PE SAME")
    leg.Draw("SAME")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.86, label)
    
    plotter.auxRatio()  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()  



    ## bottom panel
    canvas.cd()
    padB.Draw()
    padB.SetFillStyle(0)
    padB.cd()
    padB.SetGrid()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")
   
    #h_ratio.Draw("P SAME")
    h_tot_err_ratio.Draw("SAME E2")
    
    hist_plus_err_1.Draw("HIST SAME")
    hist_minus_err_1.Draw("HIST SAME")
    #hist_plus_err_2.Draw("HIST SAME")
    #hist_minus_err_2.Draw("HIST SAME")

    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/%s_mT.png" % (outDir, fitmode))
    canvas.SaveAs("%s/%s_mT.pdf" % (outDir, fitmode))
    canvas.Close()    
     
    
if __name__ == "__main__":
    
    mZbins = 60 # make auto
    
    bins_recoil_reco = lowPUcfg.bins_recoil_reco
    bins_recoil_reco[-1] = 150
    
    
    doTransverseMassPlot()
    
    doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.03)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="prefit", ratio=1.03)
    
    doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="postfit", ratio=1.11)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="postfit", ratio=1.11)