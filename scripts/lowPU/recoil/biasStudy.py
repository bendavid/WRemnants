
import sys,array,math,os,copy,fnmatch
from collections import OrderedDict

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, "scripts/lowPU/")
sys.path.insert(0, "scripts/lowPU/recoil/")
import functions
import plotter

import lz4.frame
import pickle
import narf
import numpy as np

from wremnants.datasets.datagroupsLowPU import datagroupsLowPU
from wremnants.datasets.datagroups import datagroups2016

def doPlot(tag, xMin=-100, xMax=100, yMin=1, yMax=-1, label="", xLabel="", yRatio = 1.3, rebin=1, logPos=[.20, 0.65, .5, .88]):

    hist = narf.hist_to_root(functions.readBoostHistProc(groups, tag, ["Zmumu"]))
    hist_qTrw = narf.hist_to_root(functions.readBoostHistProc(groups, tag+"_qTrw", ["Zmumu"]))
    
    hist = functions.Rebin(hist, rebin)
    hist_qTrw = functions.Rebin(hist_qTrw, rebin)
    
    cfg = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax if yMax > 0 else max([hist.GetMaximum(), hist_qTrw.GetMaximum()])    ,
        
        'xtitle'            : xLabel,
        'ytitle'            : "Events / 1 GeV",
        
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : 1-(yRatio-1),
        'ymaxR'             : yRatio,
    }   
    
    leg = ROOT.TLegend(logPos[0], logPos[1], logPos[2], logPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    leg.SetHeader(label)

    hist.SetLineColor(ROOT.kRed)
    hist.SetFillColor(0)
    hist.SetLineWidth(2)
    
    hist_qTrw.SetLineColor(ROOT.kBlue)
    hist_qTrw.SetFillColor(0)
    hist_qTrw.SetLineWidth(2)   
    
    leg.AddEntry(hist, "Unreweighted", "L")
    leg.AddEntry(hist_qTrw, "Reweighted", "L")

    h_ratio = hist.Clone("h_uncorr_ratio")
    h_ratio.Divide(hist_qTrw)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetLineWidth(2)


    
    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
    canvas.cd()
    padT.Draw()
    padT.cd()
    dummyT.Draw("HIST")
    
    hist.Draw("HIST SAME")
    hist_qTrw.Draw("HIST SAME")
    leg.Draw()
    plotter.auxRatio()
       

    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")
    
    h_ratio.Draw("HIST SAME E")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s.png" % (outDir, tag))
    canvas.SaveAs("%s/%s.pdf" % (outDir, tag))

    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()
    


    
if __name__ == "__main__":

    # compare recoil, MET, mT, lep pT, ... for MC unreweighted and qT reweighted to data


    met = "RawPFMET" # RawPFMET DeepMETReso
    flavor = "mumu" # mumu, ee
    


    outDir = "/eos/user/j/jaeyserm/www/recoil/biasStudy/highPU_zmumu"
    functions.prepareDir(outDir, False)

    #groups = datagroupsLowPU("lowPU_%s_%s_nnpdf31.pkl.lz4" % (flavor, met), flavor=flavor)
    groups = datagroups2016("mz_wlike_with_mu_eta_pt_%s_nnpdf31.pkl.lz4" % (met))
    
    recoil_qTbins = list(range(0, 30, 1)) + list(range(30, 50, 2)) + list(range(50, 70, 5)) + list(range(70, 120, 10)) + [120, 150, 300]
    
    doPlot("qT", xMin=0, xMax=100, yMin=1e3, yMax= 1e6, label="", xLabel="q_{T} (GeV)", yRatio=1.15, rebin=recoil_qTbins, logPos=[.60, 0.65, .9, .88])
    
    doPlot("transverseMass", xMin=0, xMax=120, yMin=1e3, yMax= 1e6, label="", xLabel="m_{T} (GeV)", yRatio=1.03, rebin=2)
    #doPlot("lep_pT", xMin=26, xMax=55, yMin=1e3, yMax= 1e4, label="", xLabel="p_{T}(l) (GeV)", yRatio=1.03)
    
    doPlot("recoil_corr_xy_magn", xMin=0, xMax=100, yMin=1e3, yMax= 1e6, label="", xLabel="|U| (GeV)", yRatio=1.15, logPos=[.60, 0.65, .9, .88])
    doPlot("recoil_corr_xy_para", xMin=-150, xMax=50, yMin=1e3, yMax= 1e6, label="", xLabel="U_{#parallel} (GeV)", yRatio=1.07, rebin=4)
    doPlot("recoil_corr_xy_para_qT", xMin=-50, xMax=50, yMin=1e4, yMax= 1e7, label="", xLabel="U_{#parallel} + q_{T} (GeV)", yRatio=1.03)
    doPlot("recoil_corr_xy_perp", xMin=-50, xMax=50, yMin=1e4, yMax= 1e7, label="", xLabel="U_{#perp}   (GeV)", yRatio=1.03)
    doPlot("MET_corr_xy_pt", xMin=0, xMax=50, yMin=10, yMax= 1e6, label="", xLabel="MET p_{T} (GeV)", yRatio=1.03, rebin=2)
 
 
 
