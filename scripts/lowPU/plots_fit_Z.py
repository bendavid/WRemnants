
import sys,array,math,os,copy,fnmatch
from collections import OrderedDict
import ctypes

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import plotter
import functions

import scripts.lowPU.config as lowPUcfg


def parseHist(fIn, hName):

    h = ROOT.TH1D(hName, "", len(bins_recoil_reco)-1, array.array('d', bins_recoil_reco))
    
    hIn = fIn.Get(hName)
    for iReco in range(0, len(bins_recoil_reco)-1):
    
        
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
    sigUncs = [0.036836, 0.065440, 0.045241, 0.053807, 0.050354, 0.055627]

    
    leg = ROOT.TLegend(.17, 0.57, .9, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.031)
    leg.SetNColumns(2)
    leg.SetMargin(0.7*leg.GetMargin())
    
    
    # get data
    h_data = parseHist(fIn, "obs;1")
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)       
    h_data.Scale(1, "width")

    
    st = ROOT.THStack()
    st.SetName("stack")
    
    # backgrounds
    h_top = parseHist(fIn, "expproc_TTbar_%s;1" % fitmode)
    h_top.SetFillColor(ROOT.TColor.GetColor(222, 90, 106))
    h_top.SetLineColor(ROOT.kBlack)
    h_top.SetLineWidth(1)
    h_top.SetLineStyle(1)
    h_top.Scale(1, "width")
    st.Add(h_top)
    
        
    h_ewk = parseHist(fIn, "expproc_EWK_%s;1" % fitmode)
    h_ewk.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
    h_ewk.SetLineColor(ROOT.kBlack)
    h_ewk.SetLineWidth(1)
    h_ewk.SetLineStyle(1)
    h_ewk.Scale(1, "width")
    st.Add(h_ewk)
    
    
    # sum of backgrounds
    h_tot = parseHist(fIn, "expfull_%s;1" % fitmode)
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

        h_err = parseHist(fIn, "expproc_DY_genBin%d_%s;1" % (iGenBin, fitmode))
        
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

        h_gen_sig = parseHist(fIn, "expproc_DY_genBin%d_%s;1" % (iGenBin, fitmode))
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
    
    
if __name__ == "__main__":
    
    mZbins = 60 # make auto
    
    bins_recoil_reco = lowPUcfg.bins_recoil_reco
    bins_recoil_reco[-1] = 150
    
    
    doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.03)
    doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="prefit", ratio=1.03)
    
    doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="postfit", ratio=1.11)
    doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="postfit", ratio=1.11)