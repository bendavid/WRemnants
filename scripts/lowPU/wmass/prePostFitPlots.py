
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

sys.path.insert(0, "scripts/lowPU/")
import plotter
import functions

import scripts.lowPU.config as lowPUcfg




def parseHist(hIn, bins, mTbinRange):

    h = ROOT.TH1D(hIn.GetName() + "parsed", "", len(bins)-1, array.array('d', bins))
    iBinTarget = 1
    for iBin in range(mTbinRange[0], mTbinRange[1]):
   
        err = ctypes.c_double()
        integral = hIn.IntegralAndError(iBin, iBin, err)
        h.SetBinContent(iBinTarget, integral)
        h.SetBinError(iBinTarget, err)
        #print(iBin, integral, err)
        iBinTarget += 1

    h.Scale(1, "width")
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
    
def doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.08):

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_differential_combineOutput.root" % fitcfg)
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Combine/Z/fit_%s" % fitcfg
    functions.prepareDir(outDir, remove=False)
    
    resTree = fIn.Get("fitresults")
    resTree.GetEntry(0)
    
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
    
    # extract the uncertainties
    name = "nuisance_group_impact_pmaskedexp"
    #name = "nuisance_group_impact_pmaskedexpnorm"
    #name = "nuisance_group_impact_mu"
        
    impact_hist = fIn.Get(name)
    xIdx = range(1, impact_hist.GetNbinsX()+1)
    xLabels = [impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)]    
    xIdx = [x for _,x in sorted(zip(xLabels,xIdx))]
    sigUncs = []
    for i,iBin in enumerate(xIdx):
        xLabel = xLabels[iBin-1]
        #print(xLabel)
        #print(getattr(resTree, "%s" % xLabel))
        err = getattr(resTree, "%s_err" % xLabel) / getattr(resTree, "%s" % xLabel)
        #err = getattr(resTree, "%s_err" % xLabel.replace("_mu", "_pmaskedexpnorm")) / getattr(resTree, "%s" % xLabel.replace("_mu", "_pmaskedexpnorm"))
        #err = getattr(resTree, "%s_err" % xLabel)
        sigUncs.append(err)
    
    #print(sigUncs)
    #sys.exit()

    
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
            for iBin in range(1, h_err.GetNbinsX()+1):
            
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
        print("iGenBin", iGenBin, "%d" % h_gen_sig.Integral())
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
 

def doTransverseMassPlot(flavor="mu", charge="plus", fitcfg="combined", fitmode="prefit", ratio=1.06):

    outDir = "/eos/user/j/jaeyserm/www/wmass/studies_mT/lowPU/"
    
    met = "RawPFMET"
    lumi = "1p0"
    mode = "binned"
    
    
    ## get data
    fIn_ = ROOT.TFile("/scratch/jaeyserm/CombineStudies_Wmass_mT/lowPU_mu_RawPFMET_lumi1p0_statOnly.root")
    fIn_.ls()
    h_data = copy.deepcopy(fIn_.Get("mT_corr_rec_SingleMuon_%s" % charge))
    #h_data = copy.deepcopy(h_data.ProjectionX("h_data"))
    h_data.Scale(1, "width")
    #sys.exit()
   
    
    
    #fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/lowPU_Wmass/%s/eta_inclusive_recoilOnly/LowPU_Wmass_mu_%s_lumi%s.root" % (mode, met, lumi))
    #fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/lowPU_Wmass/LowPU_Wmass_mu_%s_lumi%s.root" % (met, lumi))
    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/Wmass_mT/lowPU_mu_%s_lumi%s.root" % (met, lumi))
    
    mTbins = list(range(40, 110, 1)) + [110, 112, 114, 116, 118, 120, 125, 130, 140, 160, 180, 200]
    mTbinRange = [1, 82] # assume 81 mT bins
    if charge == "minus" and fitcfg == "combined": mTbinRange = [82, 82+81]
    

    dataLabel = "Data (#mu^{#%s})" % charge
    signalLabel = "W^{#%s} #rightarrow #mu^{#%s}#nu" % (charge, charge)
 
    #if fitcfg == "combined": label = "Combined fit, W^{#plus} + W^{#minus}"
    #else: label = "Single fit, W^{#%s}" % fitcfg
    #label += " (recoil %s)" % mode 
    label = "W^{#%s}, %s (%s)" % (charge, met, fitmode)


    leg = ROOT.TLegend(.40, 0.65, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetNColumns(2)
    
    leg_ratio = ROOT.TLegend(.18, 0.05, .40, .90)
    leg_ratio.SetBorderSize(0)
    leg_ratio.SetFillStyle(0)
    leg_ratio.SetTextSize(0.1)
    leg_ratio.SetMargin(0.7*leg.GetMargin())
    
    
    # get data
    #h_data = fIn.Get("obs;1")
    #h_data = parseHist(h_data, mTbins, mTbinRange)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)       

    st = ROOT.THStack()
    st.SetName("stack")
    
    # EWK backgrouns
    bkgs_ewk = ["WJetsToTauNu", "DY", "VV"] 
    h_ewk = None
    for bkg_ewk in bkgs_ewk:
        tmp = fIn.Get("expproc_%s_%s;1" % (bkg_ewk, fitmode))
        if h_ewk == None: h_ewk = tmp
        else: h_ewk.Add(tmp)
    h_ewk.SetName("EWK")
    h_ewk = parseHist(h_ewk, mTbins, mTbinRange)
    h_ewk.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
    h_ewk.SetLineColor(ROOT.kBlack)
    h_ewk.SetLineWidth(1)
    h_ewk.SetLineStyle(1)
    st.Add(h_ewk)

    h_top = fIn.Get("expproc_TTbar_%s;1" % fitmode)
    h_top = parseHist(h_top, mTbins, mTbinRange)
    h_top.SetFillColor(ROOT.TColor.GetColor(222, 90, 106))
    h_top.SetLineColor(ROOT.kBlack)
    h_top.SetLineWidth(1)
    h_top.SetLineStyle(1)
    st.Add(h_top)
    
    
 
    
    h_w = fIn.Get("expproc_WJetsToMuNu_%s;1" % (fitmode))
    h_w = parseHist(h_w, mTbins, mTbinRange)
    h_w.SetFillColor(ROOT.TColor.GetColor(248, 206, 104))
    h_w.SetLineColor(ROOT.kBlack)
    h_w.SetLineWidth(1)
    h_w.SetLineStyle(1)
    h_w.Scale(1.026)
    st.Add(h_w)
    
    # fakes
    h_fake = fIn.Get("expproc_Fake_%s;1" % fitmode)
    h_fake = parseHist(h_fake, mTbins, mTbinRange)
    h_fake.SetFillColor(ROOT.kGray)
    h_fake.SetLineColor(ROOT.kBlack)
    h_fake.SetLineWidth(1)
    h_fake.SetLineStyle(1)
    st.Add(h_fake)
    
    # sum of backgrounds
    h_tot = fIn.Get("expfull_%s;1" % fitmode)
    h_tot = parseHist(h_tot, mTbins, mTbinRange)
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetFillColor(0)
    h_tot.SetLineWidth(2)
    #h_tot.Scale(1.026)


    # ratio
    h_ratio = h_data.Clone("ratio")
    h_tot_noerr = h_tot.Clone("h_tot_noerr")
    for i in range(0, h_tot_noerr.GetNbinsX()+1): h_tot_noerr.SetBinError(i, 0)
    h_ratio.Divide(h_tot_noerr)
    h_ratio.SetMarkerSize(0.7)
    
    # ratio err
    h_tot_err = fIn.Get("expbkg_%s;1" % fitmode)
    h_tot_err = parseHist(h_tot_err, mTbins, mTbinRange)
    h_tot_err.SetFillColor(ROOT.kBlack)
    h_tot_err.SetMarkerSize(0)
    h_tot_err.SetLineWidth(0)
    h_tot_err.SetFillStyle(3004)
    

    h_tot_err_ratio = h_tot_err.Clone("ratio_err")
    h_tot_err_ratio_denom = h_tot_err.Clone("ratio_err_denom")
    #for i in range(0, h_tot_err_ratio_denom.GetNbinsX()+1): h_tot_err_ratio_denom.SetBinError(i, 0)
    h_tot_err_ratio.Divide(h_tot_err_ratio_denom)
    #h_tot_err_ratio.Scale(1, "width")
    

    leg.SetHeader(label)
    leg.AddEntry(h_data, dataLabel, "PE")
    leg.AddEntry(h_tot_err, "Stat. + Syst. Unc.", "F")
    leg.AddEntry(h_top, "Top", "F")
    leg.AddEntry(h_ewk, "EWK (#tau, VV, DY)", "F")
    leg.AddEntry(h_fake, "Nonprompt", "F")
    leg.AddEntry(h_w, signalLabel, "F")
    

    

    ## mass variations
    if True:
    
        hist_nom = fIn_.Get("mT_corr_rec_WJetsToMuNu_%s" % (charge))
        hist_dw = fIn_.Get("mT_corr_rec_WJetsToMuNu_massShift100MeVDown_%s" % (charge))
        hist_up = fIn_.Get("mT_corr_rec_WJetsToMuNu_massShift100MeVUp_%s" % (charge))
        
        hist_dw.Divide(hist_nom)
        hist_up.Divide(hist_nom)
        
        hist_dw = copy.deepcopy(hist_dw)
        hist_up = copy.deepcopy(hist_up)
          
        
        hist_dw.SetLineColor(ROOT.kBlue)
        hist_dw.SetFillColor(0)
        hist_dw.SetLineWidth(2)
        hist_up.SetLineColor(ROOT.kBlue)
        hist_up.SetFillColor(0)
        hist_up.SetLineWidth(2)
        
        leg_ratio.AddEntry(hist_up, "m_{W} #pm 100 MeV", "L")
    


    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 40,
        'xmax'              : 200,
        'ymin'              : 1,
        'ymax'              : 1e7,
            
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

    hist_dw.Draw("HIST SAME")
    hist_up.Draw("HIST SAME")
   
    h_ratio.Draw("P SAME")
    h_tot_err_ratio.Draw("SAME E2")
    leg_ratio.Draw()


    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/mT_%s_%sCharge_%sFit.png" % (outDir, fitmode, charge, fitcfg))
    canvas.SaveAs("%s/mT_%s_%sCharge_%sFit.pdf" % (outDir, fitmode, charge, fitcfg))
    canvas.Close()    
     
    
def extractUncertainties(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.06):

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_differential_combineOutput.root" % fitcfg)
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

    
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)
    
    name = "nuisance_group_impact_mu"
    impact_hist = fIn.Get(name) # nuisance_group_impact_nois
    xLabels = sorted([impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)])

    print("%s\t%s\t%s\t%s" % ("        ", "mu", "s_abs", "s_norm"))
    for i, xLabel in enumerate(xLabels):
    
        err_mu = getattr(tree, "%s_err" % xLabel)
        err_abs = getattr(tree, "%s_err" % xLabel.replace("_mu", "_pmaskedexp")) / getattr(tree, "%s" % xLabel.replace("_mu", "_pmaskedexp"))
        err_norm = getattr(tree, "%s_err" % xLabel.replace("_mu", "_pmaskedexpnorm")) / getattr(tree, "%s" % xLabel.replace("_mu", "_pmaskedexpnorm"))
        print("%s\t%.2f\t%.2f\t%.2f" % (xLabel.replace("_mu", ""), 100.*err_mu, 100.*err_abs, 100.*err_norm))
    
   
def plotImpactsVsPt(flavor="mumu", fitcfg="mumu", impact_type = "mu"):

    fIn = ROOT.TFile("/home/j/jaeyserm/combine/CMSSW_10_6_20/src/LowPU/LowPU_Z%s_differential_combineOutput.root" % fitcfg)
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Combine/Z/fit_%s" % fitcfg
    functions.prepareDir(outDir, remove=False)
    
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)
    
    recoil_gen = [0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 150]
    
    
    if impact_type == "abs": 
        name = "nuisance_group_impact_pmaskedexp"
        suffix = "abs"
        lab = "#sigma_{abs}"
    elif impact_type == "norm": 
        name = "nuisance_group_impact_pmaskedexpnorm"
        suffix = "norm"
        lab = "#sigma_{norm}"
    else: 
        name = "nuisance_group_impact_mu"
        suffix = "mu"
        lab = "#mu"
        
    impact_hist = fIn.Get(name)
    xIdx = range(1, impact_hist.GetNbinsX()+1)
    xLabels = [impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)]
    yLabels = [impact_hist.GetYaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsY()+1)]
    
    xIdx = [x for _,x in sorted(zip(xLabels,xIdx))]

    scale = 100
    rounded = 2
    dfs = []
    graphs = {}
    
    for k,yLabel in enumerate(yLabels):
        
        g = ROOT.TGraph()
        g.SetName(yLabel)
        g.SetTitle(yLabel)
        for i,iBin in enumerate(xIdx): g.SetPoint(i, 0.5*(recoil_gen[i]+recoil_gen[i+1]), 100.*impact_hist.GetBinContent(iBin, k+1))
       
        graphs[yLabel] = g
        
    # total
    g = ROOT.TGraph()
    g.SetName("total")
    g.SetTitle("Total")
    for i,iBin in enumerate(xIdx):
        xLabel = xLabels[iBin-1]
        if impact_type == "abs": err = getattr(tree, "%s_err" % xLabel.replace("_mu", "_pmaskedexp")) / getattr(tree, "%s" % xLabel.replace("_mu", "_pmaskedexp"))
        elif impact_type == "norm": err = getattr(tree, "%s_err" % xLabel.replace("_mu", "_pmaskedexpnorm")) / getattr(tree, "%s" % xLabel.replace("_mu", "_pmaskedexpnorm"))
        else: err = getattr(tree, "%s_err" % xLabel)
        g.SetPoint(i, 0.5*(recoil_gen[i]+recoil_gen[i+1]), 100.*err)
    graphs["total"] = g
    
    ######################################
    c = ROOT.TCanvas("c", "c", 1000, 750)
    c.SetTopMargin(0.06)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.14)
    c.SetLogy()
    
    # dummy
    dummy = ROOT.TH1D("h", "h", 1, 0, 200)
    dummy.GetXaxis().SetTitle("p_{T}^{GEN}(Z) (GeV)")
    dummy.GetXaxis().SetRangeUser(0, 200)
    dummy.GetXaxis().SetTitleFont(43)
    dummy.GetXaxis().SetTitleSize(40)
    dummy.GetXaxis().SetLabelFont(43)
    dummy.GetXaxis().SetLabelSize(35)
    dummy.GetXaxis().SetTitleOffset(1.1*dummy.GetXaxis().GetTitleOffset())
    dummy.GetXaxis().SetLabelOffset(1.2*dummy.GetXaxis().GetLabelOffset())
    dummy.GetYaxis().SetTitle("Impact on %s (%%)" % lab)
    dummy.GetYaxis().SetRangeUser(0.1, 10)
    dummy.GetYaxis().SetTitleFont(43)
    dummy.GetYaxis().SetTitleSize(40)
    dummy.GetYaxis().SetLabelFont(43)
    dummy.GetYaxis().SetLabelSize(35)
    dummy.GetYaxis().SetTitleOffset(1.4*dummy.GetYaxis().GetTitleOffset())
    dummy.GetYaxis().SetLabelOffset(1.4*dummy.GetYaxis().GetLabelOffset())
    dummy.Draw("HIST")
    
    leg = ROOT.TLegend(.7, 0.30, .95, .8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    #leg.SetNColumns(3)
    leg.SetMargin(0.7*leg.GetMargin())

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kMagenta, ROOT.kCyan, ROOT.kGray, ROOT.kYellow, ROOT.kSpring]
    graphs_plot = ["total", "QCDscale", "CMS_lumi_lowPU", "stat", "CMS_recoil_stat", "pdfNNPDF31", "binByBinStat", "CMS_lepton_eff", "CMS_bkg_norm", "CMS_prefire17"]
    for i, grName in enumerate(graphs_plot):
    
        g = graphs[grName]
        g.SetLineWidth(2)
        g.SetLineColor(colors[i])
        g.SetMarkerStyle(8)
        g.SetMarkerSize(1)
        g.SetMarkerColor(colors[i])
        g.Draw("SAME LP")
        leg.AddEntry(g, g.GetTitle(), "LP")

    leg.Draw("SAME")
    c.Modify()
    c.Update()
    
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(30)
    tr = 1.-c.GetRightMargin()
    latex.DrawLatex(tr, 0.94, "199 pb^{#minus1} (13 TeV)")
    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.15, 0.985, "#bf{CMS} #scale[0.7]{#it{Preliminary}}")
    
    
    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.70, 0.88, "DY #rightarrow #mu^{#plus}#mu^{#minus}")

    c.SaveAs("%s/impacts_%s.png" % (outDir, suffix))
    c.SaveAs("%s/impacts_%s.pdf" % (outDir, suffix))
    c.Close()  
 
if __name__ == "__main__":
    
    mZbins = 60 # make auto
    
    bins_recoil_reco = lowPUcfg.bins_recoil_reco
    bins_recoil_reco[-1] = 150
    
    doTransverseMassPlot(flavor="mu", charge="plus", fitcfg="combined", fitmode="prefit", ratio=1.11)
    doTransverseMassPlot(flavor="mu", charge="minus", fitcfg="combined", fitmode="prefit", ratio=1.11)
    
    doTransverseMassPlot(flavor="mu", charge="plus", fitcfg="combined", fitmode="postfit", ratio=1.11)
    doTransverseMassPlot(flavor="mu", charge="minus", fitcfg="combined", fitmode="postfit", ratio=1.11)
    
    #doTransverseMassPlot()
    
    #doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="prefit", ratio=1.06)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="prefit", ratio=1.03)
    
    #doRecoilPlot(flavor="mumu", fitcfg="mumu", fitmode="postfit", ratio=1.11)
    #doRecoilPlot(flavor="ee", fitcfg="ee", fitmode="postfit", ratio=1.11)
    
    
    #extractUncertainties()
    
    
    #plotImpactsVsPt(impact_type = "mu")
    #plotImpactsVsPt(impact_type = "abs")
    #plotImpactsVsPt(impact_type = "norm")
