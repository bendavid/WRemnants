
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

sys.path.insert(0, "scripts/lowPU/")
import plotter
import functions

import scripts.lowPU.config as lowPUcfg




def doPlot(flavor="mu", charge="plus", syst="", ratio=1.06):

    
    
    met = "RawPFMET"
    lumi = "1p0"
    mode = "binned"
    fIn = ROOT.TFile("/scratch/jaeyserm/CombineStudies_Wmass/LowPU_Wmass_mu_RawPFMET_lumi1p0.root")
    
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Projection/Templates/%s_%s_%s" % (met, mode, charge)
    functions.prepareDir(outDir, remove=False)
    
    mTbins = list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200]
    mTbinRange = [1, 82] # assume 81 mT bins
    if charge == "minus" and fitcfg == "combined": mTbinRange = [82, 82+81]
    

    dataLabel = "Data (#mu^{#%s})" % charge
    signalLabel = "W^{#%s} #rightarrow #mu^{#%s}#nu" % (charge, charge)
 

    leg = ROOT.TLegend(.45, 0.6, .95, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    #leg.SetNColumns(2)
    leg.SetMargin(0.7*leg.GetMargin())
    

    
    h_base = fIn.Get("mT_corr_rec_W%sJetsToMuNu_%s;1" % (charge, charge))
    if h_base.ClassName() == "TH2D": h_base = h_base.ProjectionX("h_base", 1, 1)

    h_base.Scale(1., "width")
    h_base.SetLineColor(ROOT.kBlack)
    h_base.SetLineWidth(2)
    h_base.SetLineStyle(1)
    
    hUp = fIn.Get("mT_corr_rec_W%sJetsToMuNu_%sUp_%s;1" % (charge, syst, charge))
    if hUp.ClassName() == "TH2D": hUp = hUp.ProjectionX("hUp", 1, 1)
    hDw = fIn.Get("mT_corr_rec_W%sJetsToMuNu_%sDown_%s;1" % (charge, syst, charge))
    if hDw.ClassName() == "TH2D": hDw = hDw.ProjectionX("hDw", 1, 1)
    hUp.Scale(1., "width")
    hDw.Scale(1., "width")
    
    hUp_ratio = h_base.Clone("hUp_ratio")
    hUp_ratio.Divide(hUp)
    hUp_ratio.SetLineColor(ROOT.kRed)
    hUp_ratio.SetLineWidth(2)
    hUp_ratio.SetLineStyle(1)
    
    hDw_ratio = h_base.Clone("hDw_ratio")
    hDw_ratio.Divide(hDw)
    hDw_ratio.SetLineColor(ROOT.kBlue)
    hDw_ratio.SetLineWidth(2)
    hDw_ratio.SetLineStyle(1)
    
    # mass shift
    hUp_mass = fIn.Get("mT_corr_rec_W%sJetsToMuNu_%sUp_%s;1" % (charge, "massShift100MeV", charge))
    if hUp_mass.ClassName() == "TH2D": hUp_mass = hUp_mass.ProjectionX("hUp_mass", 1, 1)
    hDw_mass = fIn.Get("mT_corr_rec_W%sJetsToMuNu_%sDown_%s;1" % (charge, "massShift100MeV", charge))
    if hDw_mass.ClassName() == "TH2D": hDw_mass = hDw_mass.ProjectionX("hDw_mass", 1, 1)
    hUp_mass.Scale(1., "width")
    hDw_mass.Scale(1., "width")
    hUp_mass_ratio = h_base.Clone("hUp_mass_ratio")
    hDw_mass_ratio = h_base.Clone("hDw_mass_ratio")
    hUp_mass_ratio.Divide(hUp_mass)
    hDw_mass_ratio.Divide(hDw_mass)
    hUp_mass_ratio.SetLineColor(ROOT.kGray+1)
    hUp_mass_ratio.SetLineWidth(2)
    hUp_mass_ratio.SetLineStyle(1)
    hDw_mass_ratio.SetLineColor(ROOT.kGray+1)
    hDw_mass_ratio.SetLineWidth(2)
    hDw_mass_ratio.SetLineStyle(1)
    
    
    leg.SetHeader(syst)
    leg.AddEntry(h_base, "Base template", "L")
    leg.AddEntry(hUp_ratio, "Up variation", "L")
    leg.AddEntry(hDw_ratio, "Down variation", "L")
    leg.AddEntry(hDw_mass_ratio, "100 MeV mass variation", "L")
    
    for i in range(1, h_base.GetNbinsX()+1):
        print(h_base.GetBinCenter(i), h_base.GetBinContent(i), hUp_ratio.GetBinContent(i), hDw_ratio.GetBinContent(i))


    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 40,
        'xmax'              : 200,
        'ymin'              : 10,
        'ymax'              : 1e6,
            
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
        

    h_base.Draw("HIST SAME")
    hUp.Draw("HIST SAME")
    hDw.Draw("HIST SAME")
    leg.Draw("SAME")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.4, 0.86, label)
    
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
    hUp_mass_ratio.Draw("HIST SAME")
    hDw_mass_ratio.Draw("HIST SAME")
    hUp_ratio.Draw("HIST SAME")
    hDw_ratio.Draw("HIST SAME")


    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/%s.png" % (outDir, syst))
    canvas.SaveAs("%s/%s.pdf" % (outDir, syst))
    canvas.Close()    
     
    

if __name__ == "__main__":
    
    mZbins = 60 # make auto
    
    bins_recoil_reco = lowPUcfg.bins_recoil_reco
    bins_recoil_reco[-1] = 150
    
    #doPlot(flavor="mu", charge="plus", syst="massShift100MeV", ratio=1.03)
    
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_0", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_1", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_2", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_3", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_4", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="recoilSyst_target_para_5", ratio=1.03)
    
    
    doPlot(flavor="mu", charge="plus", syst="CMS_scale_m_ieta0", ratio=1.01)
    
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q0muF", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q0muR", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q0muRmuF", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q1muF", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q1muR", ratio=1.03)
    #doPlot(flavor="mu", charge="plus", syst="QCDscale_q1muRmuF", ratio=1.03)
    
    #doPlot(flavor="mu", charge="plus", syst="pdf100NNPDF31", ratio=1.03)

    
    #recoilSyst_target_para_1
    #doTransverseMassPlot(flavor="mu", charge="minus", fitcfg="combined", fitmode="prefit", ratio=1.08)
    
    #doTransverseMassPlot(flavor="mu", charge="plus", fitcfg="combined", fitmode="postfit", ratio=1.08)
    #doTransverseMassPlot(flavor="mu", charge="minus", fitcfg="combined", fitmode="postfit", ratio=1.08)
    
    #doTransverseMassPlot()
    
    #doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.06)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="prefit", ratio=1.03)
    
    #doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="postfit", ratio=1.11)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="postfit", ratio=1.11)
    
    
    #extractUncertainties()
    
    
    #plotImpactsVsPt(impact_type = "mu")
    #plotImpactsVsPt(impact_type = "abs")
    #plotImpactsVsPt(impact_type = "norm")
