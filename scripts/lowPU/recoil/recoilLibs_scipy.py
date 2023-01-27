
import sys,array,math,os,copy,decimal,random
import numpy as np
import ctypes
import json
import time
import threading

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import functions
import plotter
import recoilFunctions as rf

import lz4.frame
import pickle
import narf
import pickle
import tensorflow as tf

import fitter
import hist
import boost_histogram as bh

s = hist.tag.Slicer()

mirror = False

__topRight__ = "199 pb^{#minus1} (13 TeV)"

__min_recoil__, __max_recoil__ = -500, 500
__min_mean__, __max_mean__ = -50, 500
__min_sigma__, __max_sigma__ = 0.1, 100

def addParam(jsOut, p, func, vals, bnds=None):
    if p in jsOut:
        raise Exception("Parameter %d already in dictionary" % p)
    jsOut[p] = {}
    jsOut[p]['func'] = func
    jsOut[p]['funcName'] = None
    jsOut[p]['nParams'] = len(vals)
    for i in range(0, len(vals)):
        jsOut[p]["p%d"%i] = vals[i]
        jsOut[p]["p%d_err"%i] = 0
        jsOut[p]["p%d_isCte"%i] = False
        jsOut[p]["p%d_bnds" % i] = (None, None) if bnds==None else bnds[i]


def diagonalize(cov_matrix, parms_nom, sigma=1):

    parms_nom = np.array(parms_nom)
    npars = len(parms_nom)
    eig_vals, eig_vec = np.linalg.eig(cov_matrix)
    eig_vec_inv = np.linalg.inv(eig_vec) # = transposed
    
    parms_nom_eig = np.dot(parms_nom, eig_vec) # parameters in diagonalized space
    parms_pert = []
    for i, eig in enumerate(eig_vals): 


        parms_dnom_eig = copy.deepcopy(parms_nom_eig)
        if eig < 0:
            print("Eigenvalue negative for parIdx=%d, value=%f" % (i, eig))
        else:
            parms_dnom_eig[i] += sigma*math.sqrt(eig)
        parms_dnom = np.dot(parms_dnom_eig, eig_vec_inv)

        parms_pert.append(parms_dnom)

    return parms_pert
    



def doFitMultiGauss_mt(bhist, comp, fitCfg, label, outDir_, binning_qT, propagate=[], doFit=True, rebin=1, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, bkgCfg={}, yRatio = 1.15, recoilLow=-100, recoilHigh=100):

    average_eval_PDF = False
 
    cfgPlot = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)" if comp == "para" else "Recoil U_{#perp}   (GeV)",
        'ytitle'            : "Events" ,
        
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }
    
    axis_name = "recoil_perp" if comp == "perp" else "recoil_para_qT"
    #axis_name = "recoil_para_qT"
    qTlow, qThigh = min(binning_qT), max(binning_qT)
    bhist =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)], axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]  #"recoil_perp": s[complex(0,recoilLow):complex(0,recoilHigh)
    #h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}] 
    
    
    bins_qT = bhist.axes[0].edges
    bins_recoil = bhist.axes[1].edges
    centers_qT = bhist.axes[0].centers
    centers_recoil = bhist.axes[1].centers
    centers_recoil_tf = tf.constant(centers_recoil, dtype=tf.float64)
    bins_recoil_tf = tf.constant(bins_recoil, dtype=tf.float64)
    ##bins_widths = [bins_recoil[i]-bins_recoil[i-1] for i in range(1, len(bins_recoil))]
    
    functions.prepareDir(outDir_, True)
    func_name = fitCfg['func_name']
    func_parms = []
    if "func_parms_init" in fitCfg: func_parms = fitCfg['func_parms_init'] # initial values array
    else: # get the parameters from the normal dict structure
        for k in range(0, fitCfg['nParams']):
            for l in range(0, fitCfg["p%d"%k]['nParams']):
                func_parms.append(fitCfg["p%d"%k]["p%d"%l])
    func_model = getattr(rf, func_name)
    
    from multiprocessing import Process, Queue

    def singleFit(iBin, bhist, func_parms, cfg, ret):
    
        average_eval_PDF = False
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        qT_tf = tf.constant(qT, dtype=tf.float64)
        qTlow, qThigh = binning_qT[iBin-1], binning_qT[iBin]
        print("############# Recoil bin %d [%.2f, %.2f]" % (iBin, qTlow, qThigh))
        
        h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h =  h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum
        yield_, yield_err = h.sum().value, math.sqrt(h.sum().variance)
   
        start = time.time()
        res = fitter.fit_hist(h, func_model, np.array(func_parms), max_iter = 5, edmtol = 1e-5, mode = "nll") # nll chisq_normalized
        end = time.time()
        
        dt = end-start
        ret[0] = dt
            
        parms_vals = res['x']
        fit_status = res['status']
        cov_status = res['covstatus']
        hess_eigvals = res['hess_eigvals']
        cov_matrix = res['cov']
        print("Done", dt)
        return
        #func_parms = [parms_vals[i] if propagate[i] else func_parms[i] for i in range(0, len(func_parms))]
        #if propagate: func_parms = parms_vals
            
        #print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
        #print(parms_vals)
        
        if average_eval_PDF:
            pdf_values_bins = func_model([qT_tf], parms_vals).numpy()
            pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]
                
        else:
            pdf_values = func_model([centers_recoil_tf], parms_vals).numpy()
        
        import plotter
        plotter.cfg = cfg
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        
        '''
        hist_root = narf.hist_to_root(h) 
        hist_root.SetName("hist_%d" % random.randint(1e4, 1e5))
        #hist_root.Scale(1.0, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")

        g_pdf = ROOT.TGraphErrors()
        g_pdf.SetName("pdf_%d" % random.randint(1e4, 1e5))
        for i in range(0, len(centers_recoil)): 
            g_pdf.SetPoint(i, centers_recoil[i], pdf_values[i]*yield_)
        g_pdf.SetLineColor(ROOT.kBlue)
        g_pdf.SetLineWidth(3)
        g_pdf.Draw("L SAME")

        histRatio = hist_root.Clone("pdf_ratio%d" % random.randint(1e4, 1e5))
        histRatio.Reset("ACE")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        chi2, histRatio = ratioHist_tgraph(hist_root, histRatio, g_pdf)
        

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, label)
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTlow, qThigh))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        if doFit:
            latex.DrawLatex(0.20, 0.55, "#color[%d]{Fit status = %d}" % (2 if fit_status != 0 else 1, fit_status))
            latex.DrawLatex(0.20, 0.50, "#color[%d]{Covariance status = %d}" % (2 if cov_status != 0 else 1, cov_status))
 
        latex.SetTextSize(0.035)
         '''
        
        # save output
        outDict[iBin] = {}
        outDict[iBin]['yield'] = yield_
        outDict[iBin]['yield_err'] = yield_err

        for i in range(0, len(func_parms)):
            val, err = parms_vals[i], math.sqrt(cov_matrix[i][i]) if doFit and cov_status == 0 and cov_matrix[i][i] >= 0 else -2
            outDict[iBin]["p%d"%i], outDict[iBin]["p%d_err"%i] = val, err
            #latex.DrawLatex(0.7, 0.87-i*0.04, "p_{%d} = %.3f #pm %.3f" % (i, val, err))
        
        #plotter.auxRatio()
        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        dummyB.Draw("HIST")
        #histRatio.Draw("SAME E0")
        #dummyL.Draw("SAME")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, iBin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, iBin))
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
     
        return
        outDict[iBin]['chi2'] = chi2   
        ret[0] = arg
        
        
        
        
        
    
        
    
    threads = []
    outDict = {}
    for iBin in range(1, len(binning_qT)):
        #if iBin > 5: break
        outDict[iBin] = [None]*1
        #thread = myThread(1, "Thread-1", 1)
        thread = Process(target=singleFit, args=(iBin, bhist, func_parms, cfgPlot, outDict[iBin],))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()
        
    print("Done")    
        
    quit()   
        
    outDict = {}
    outDict['func_name'] = func_name
               
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()
    hist_root_tot = None
    pdf_tot = [0]*len(centers_recoil)
    yields, yields_err = [], []
    for iBin in range(1, len(binning_qT)):
    
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        qT_tf = tf.constant(qT, dtype=tf.float64)
        qTlow, qThigh = binning_qT[iBin-1], binning_qT[iBin]
        print("############# Recoil bin %d [%.2f, %.2f]" % (iBin, qTlow, qThigh))
        continue
        #h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)], axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]
        h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h =  h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum
        #h = resetErrors(h)
        yield_, yield_err = h.sum().value, math.sqrt(h.sum().variance)
        yields.append(yield_)
        yields_err.append(yield_err)
        if doFit:
            print(" -> Start fitting")
            start = time.time()
            res = fitter.fit_hist(h, func_model, np.array(func_parms), max_iter = 5, edmtol = 1e-5, mode = "nll") # nll chisq_normalized
            end = time.time()
            
            parms_vals = res['x']
            fit_status = res['status']
            cov_status = res['covstatus']
            hess_eigvals = res['hess_eigvals']
            cov_matrix = res['cov']
            func_parms = [parms_vals[i] if propagate[i] else func_parms[i] for i in range(0, len(func_parms))]
            #if propagate: func_parms = parms_vals
            
            print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
            print(parms_vals)

            if average_eval_PDF:
                pdf_values_bins = func_model([qT_tf], parms_vals).numpy()
                pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]
                
            else:
                pdf_values = func_model([centers_recoil_tf], parms_vals).numpy()
        else:
            parms_vals = func_parms
            if average_eval_PDF:
                pdf_values_bins = func_model([qT_tf, bins_recoil_tf], parms_vals).numpy()
                pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]
            else:
                pdf_values = func_model([qT_tf, centers_recoil_tf], parms_vals).numpy()
            #pdf_values = func_model([qT_tf, centers_recoil_tf], parms_vals).numpy()
            cov_status, fit_status = -1, -1   
            # errors?

        plotter.cfg = cfgPlot
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        
        hist_root = narf.hist_to_root(h) 
        #hist_root.Scale(1.0, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")
        if hist_root_tot == None: hist_root_tot = hist_root.Clone("hist_root_tot")
        else: hist_root_tot.Add(hist_root)

        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers_recoil)): 
            g_pdf.SetPoint(i, centers_recoil[i], pdf_values[i]*yield_)
            pdf_tot[i] += pdf_values[i]*yield_
        g_pdf.SetLineColor(ROOT.kBlue)
        g_pdf.SetLineWidth(3)
        g_pdf.Draw("L SAME")

        histRatio = hist_root.Clone("ratio")
        histRatio.Reset("ACE")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        chi2, histRatio = ratioHist_tgraph(hist_root, histRatio, g_pdf)
        

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, label)
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTlow, qThigh))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        if doFit:
            latex.DrawLatex(0.20, 0.55, "#color[%d]{Fit status = %d}" % (2 if fit_status != 0 else 1, fit_status))
            latex.DrawLatex(0.20, 0.50, "#color[%d]{Covariance status = %d}" % (2 if cov_status != 0 else 1, cov_status))
 
        latex.SetTextSize(0.035)
        
        
        # save output
        outDict[iBin] = {}
        outDict[iBin]['yield'] = yield_
        outDict[iBin]['yield_err'] = yield_err

        for i in range(0, len(func_parms)):
            val, err = parms_vals[i], math.sqrt(cov_matrix[i][i]) if doFit and cov_status == 0 and cov_matrix[i][i] >= 0 else -2
            outDict[iBin]["p%d"%i], outDict[iBin]["p%d_err"%i] = val, err
            latex.DrawLatex(0.7, 0.87-i*0.04, "p_{%d} = %.3f #pm %.3f" % (i, val, err))
        
        plotter.auxRatio()
        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        dummyB.Draw("HIST")
        histRatio.Draw("SAME E0")
        dummyL.Draw("SAME")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, iBin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, iBin))
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
     

        outDict[iBin]['chi2'] = chi2
        g_chi2.SetPoint(iBin-1, qT, chi2)
        g_yields.SetPoint(iBin-1, qT, yield_)
        g_yields.SetPointError(iBin-1, 0, yield_err)
        
    quit()
    outDict['qTbins'] = binning_qT
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=min(binning_qT), xMax=max(binning_qT))
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=min(binning_qT), xMax=max(binning_qT))
    
    
    # uncertainties
    doUnc = 'nStatVars' in fitCfg
    if doUnc:
        nStatVars = fitCfg['nStatVars']
        hist_root_tot_unc = hist_root_tot.Clone("hist_root_tot_unc")
        hist_root_tot_unc.SetFillColor(ROOT.kBlack)
        hist_root_tot_unc.SetMarkerSize(0)
        hist_root_tot_unc.SetLineWidth(0)
        hist_root_tot_unc.SetFillStyle(3004)    
        hist_root_tot_unc.Reset("ICE")
        
        hist_root_tot_unc_ratio = hist_root_tot.Clone("hist_root_tot_unc_ratio")
        hist_root_tot_unc_ratio.SetFillColor(ROOT.kBlack)
        hist_root_tot_unc_ratio.SetMarkerSize(0)
        hist_root_tot_unc_ratio.SetLineWidth(0)
        hist_root_tot_unc_ratio.SetFillStyle(3004)    
        hist_root_tot_unc_ratio.Reset("ICE") # reset errors

        for nStat in range(0, nStatVars):
            nStatName = "stat%d" % nStat
            #if nStat > 0: continue
            func_parms = []
            for k in range(0, fitCfg['nParams']):
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                    
            pdf_tot_var = [0]*len(centers_recoil)
            for iBin in range(1, len(binning_qT)):
            
                qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
                #print(qT)
                qT_tf = tf.constant(qT, dtype=tf.float64)
                
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var[i] += pdf_values_var[i]*yields[iBin-1] # iBin starts from 1
                        
            for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1):
                pdf_tot_eval = pdf_tot[i-1]
                
                if pdf_tot_eval > 0: 
                    err = abs(pdf_tot[i-1] - pdf_tot_var[i-1])
                else:
                    err = 0
                #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinContent(i), pdf_eval_unc, pdf_eval, err)
                hist_root_tot_unc_ratio.SetBinError(i, hist_root_tot_unc_ratio.GetBinError(i) + err*err) # add in quadrature
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i) + err*err)
                #print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])
           # quit() 
                    

        for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1): 
            pdf_tot_eval = pdf_tot[i-1]
            if pdf_tot_eval > 0: 
                hist_root_tot_unc_ratio.SetBinError(i, (hist_root_tot_unc_ratio.GetBinError(i)**0.5)/pdf_tot_eval)
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i)**0.5)
            hist_root_tot_unc_ratio.SetBinContent(i, 1)
            hist_root_tot_unc.SetBinContent(i, pdf_tot_eval)
            #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinError(i), histUnc.GetBinError(i))
            #if i == 100: print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])
            
            
            
    
    # global plot
    plotter.cfg = cfgPlot
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
        
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.SetMarkerStyle(20)
    hist_root_tot.SetMarkerColor(ROOT.kBlack)
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.Draw("PE SAME")

    g_pdf = ROOT.TGraphErrors()
    for i in range(0, len(centers_recoil)): 
        g_pdf.SetPoint(i, centers_recoil[i], pdf_tot[i])
    g_pdf.SetLineColor(ROOT.kBlue)
    g_pdf.SetLineWidth(3)
    if doUnc: hist_root_tot_unc.Draw("E2 SAME") 
    g_pdf.Draw("L SAME")
        
    histRatio = hist_root_tot.Clone("ratio")
    histRatio.Reset("ACE")
    histRatio.SetMarkerStyle(8)
    histRatio.SetMarkerSize(0.7)
    histRatio.SetMarkerColor(ROOT.kBlack)
    histRatio.SetLineColor(ROOT.kBlack)
    histRatio.SetLineWidth(1)
    chi2, histRatio = ratioHist_tgraph(hist_root_tot, histRatio, g_pdf)
        

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.85, label)
    latex.DrawLatex(0.20, 0.80, "#chi^{2} = %.3f" % chi2)

    plotter.auxRatio()
    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    dummyB.Draw("HIST")
    histRatio.Draw("SAME E0")
    if doUnc: hist_root_tot_unc_ratio.Draw("E2 SAME")
    dummyL.Draw("SAME")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/global.png" % outDir_)
    canvas.SaveAs("%s/global.pdf" % outDir_)
    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()



def doFitMultiGauss_plot(bhist, comp, fitCfg, procLabel, metLabel, outDir_, binning_qT, rebin=1, yMin=1e-6, yMax=1e4, bkgCfg={}, yRatio = 1.15, recoilLow=-100, recoilHigh=100, plotSignal=False):
    
    average_eval_PDF = False
 
    cfgPlot = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : recoilLow,
        'xmax'              : recoilHigh,
        'ymin'              : yMin,
        'ymax'              : yMax,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)" if comp == "para" else "Recoil U_{#perp}   (GeV)",
        'ytitle'            : "Events" ,
        
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }
    
    axis_name = "recoil_perp" if comp == "perp" else "recoil_para_qT"
    qTlow, qThigh = min(binning_qT), max(binning_qT)
    bhist =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)], axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]
    bhist = bhist[{axis_name: s[::bh.rebin(rebin)]}] # induces overflows?
    
    doUnc = 'nStatVars' in fitCfg
    bins_qT = bhist.axes[0].edges
    bins_recoil = bhist.axes[1].edges
    centers_qT = bhist.axes[0].centers
    centers_recoil = bhist.axes[1].centers
    centers_recoil_tf = tf.constant(centers_recoil, dtype=tf.float64)
    bins_recoil_tf = tf.constant(bins_recoil, dtype=tf.float64)
    
    functions.prepareDir(outDir_, True)
    func_name = fitCfg['func_name']
    func_parms = []
    for k in range(0, fitCfg['nParams']):
        for l in range(0, fitCfg["p%d"%k]['nParams']):
            func_parms.append(fitCfg["p%d"%k]["p%d"%l])
            #print("n", fitCfg["p%d"%k]["p%d"%l])
    func_model = getattr(rf, func_name)

               
    # backgrounds
    args = ()
    args_nobkg = ()
    if len(bkgCfg) > 0:
        args_, args_nobkg_ = {}, {}
        for i,bkg in enumerate(bkgCfg["procs"]):
            parms_cond = [] # construct array of conditional pdf parameters
            for k in range(0, bkgCfg['parms'][i]['nParams']):
                for l in range(0, bkgCfg['parms'][i]["p%d"%k]['nParams']):
                    parms_cond.append(tf.constant(bkgCfg['parms'][i]["p%d"%k]["p%d"%l], dtype=tf.float64))
            args_['parms_%s'%bkg] = parms_cond
            args_nobkg_['parms_%s'%bkg] = parms_cond
            yield_fracs = [] # construct array of yield fractions
            for iBin in range(1, len(binning_qT)):
                yield_fracs.append(bkgCfg['norms'][i]*bkgCfg['yields'][i]["%d"%iBin]["yield"]/bkgCfg['data_yields']["%d"%iBin]['yield'])
            args_['fracs_%s'%bkg] = tf.constant(yield_fracs, dtype=tf.float64)
            args_nobkg_['fracs_%s'%bkg] = tf.constant([0.0]*len(binning_qT), dtype=tf.float64)
        args = (args_)
        args_nobkg = (args_nobkg_)
               
    g_chi2 = ROOT.TGraphErrors()
    hist_root_tot = None
    pdf_tot, sig_pdf_tot = [0]*len(centers_recoil), [0]*len(centers_recoil)
    yields, yields_err, pdf_fracs = [], [], []
    for iBin in range(1, len(binning_qT)):
    
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        qT_tf = tf.constant(qT, dtype=tf.float64)
        qTlow, qThigh = binning_qT[iBin-1], binning_qT[iBin]
        print("############# Recoil bin %d [%.2f, %.2f]" % (iBin, qTlow, qThigh))

        h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h =  h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum
        #h = h[{axis_name: s[::bh.rebin(rebin)]}]
        yield_, yield_err = h.sum().value, math.sqrt(h.sum().variance)
        yields.append(yield_)
        yields_err.append(yield_err)

        if average_eval_PDF:
            pdf_values_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args).numpy()
            pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]
        else:
            pdf_values = func_model([qT_tf, centers_recoil_tf], func_parms, args).numpy()
        pdf_frac_range = sum(pdf_values)
        pdf_fracs.append(pdf_frac_range)
        
        hist_root = narf.hist_to_root(h)
        hist_root.Scale(1.0, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        if hist_root_tot == None: hist_root_tot = hist_root.Clone("hist_root_tot")
        else: hist_root_tot.Add(hist_root)

        if doUnc or False:
            nStatVars = fitCfg['nStatVars']
            hist_root_tot_unc = hist_root_tot.Clone("hist_root_tot_unc_%d"%iBin)
            hist_root_tot_unc.SetFillColor(18)
            hist_root_tot_unc.SetMarkerSize(0)
            hist_root_tot_unc.SetLineWidth(0)
            hist_root_tot_unc.Reset("ICE")
            
            hist_root_tot_unc_ratio = hist_root_tot.Clone("hist_root_tot_unc_ratio_%d"%iBin)
            hist_root_tot_unc_ratio.SetFillColor(18)
            hist_root_tot_unc_ratio.SetMarkerSize(0)
            hist_root_tot_unc_ratio.SetLineWidth(0)
            hist_root_tot_unc_ratio.Reset("ICE") # reset errors
            
            for nStat in range(0, nStatVars):
            
                ### PLUS VARIATION
                pdf_tot_var_up = [0]*len(centers_recoil)
                nStatName = "stat%d_p" % nStat
                func_parms = []
                for k in range(0, fitCfg['nParams']):
                    for l in range(0, fitCfg["p%d"%k]['nParams']):
                        func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                        
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms, args).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var_up[i] += pdf_values_var[i]*yields[iBin-1]/pdf_fracs[iBin-1] # iBin starts from 1
                         
                
                ### MINUS VARIATION
                pdf_tot_var_dw = [0]*len(centers_recoil)
                nStatName = "stat%d_m" % nStat
                func_parms = []
                for k in range(0, fitCfg['nParams']):
                    for l in range(0, fitCfg["p%d"%k]['nParams']):
                        func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms, args).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var_dw[i] += pdf_values_var[i]*yields[iBin-1]/pdf_fracs[iBin-1] # iBin starts from 1     
                   
                # make average of plus and minus variation  
                for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1):
                    pdf_tot_eval = pdf_values[i-1]*yield_/pdf_frac_range
                    if pdf_tot_eval > 0:
                        err_up = abs(pdf_tot_eval - pdf_tot_var_up[i-1])
                        err_dw = abs(pdf_tot_eval - pdf_tot_var_dw[i-1])
                        err = 0.5*(err_up+err_dw) # max(err_up, err_dw) 
                        #err = 0.5*abs(pdf_tot_var_up[i-1] - pdf_tot_var_dw[i-1])
                        #err = err_up
                        #err = abs(pdf_tot_eval - pdf_tot_var[i-1])
                    else:
                        err = 0
                    hist_root_tot_unc_ratio.SetBinError(i, hist_root_tot_unc_ratio.GetBinError(i) + err*err) # add in quadrature
                    hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i) + err*err)


            for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1): 
                pdf_tot_eval = pdf_values[i-1]*yield_/pdf_frac_range
                if pdf_tot_eval > 0: 
                    hist_root_tot_unc_ratio.SetBinError(i, (hist_root_tot_unc_ratio.GetBinError(i)**0.5)/pdf_tot_eval)
                    hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i)**0.5)
                hist_root_tot_unc_ratio.SetBinContent(i, 1)
                hist_root_tot_unc.SetBinContent(i, pdf_tot_eval)



        if plotSignal:
            if average_eval_PDF:
                pdf_values_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args_nobkg).numpy()
                sig_pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]
            else:
                sig_pdf_values = func_model([qT_tf, centers_recoil_tf], func_parms, args_nobkg).numpy()
            sig_pdf_frac_range = sum(sig_pdf_values)
           
        plotter.cfg = cfgPlot
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        #if doUnc: hist_root_tot_unc.Draw("E4 SAME")
        hist_root.Draw("PE SAME")
        

        yield_ /= rebin # reduce yield to accomodate for binning
        g_pdf = ROOT.TGraphErrors()
        g_sig_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers_recoil)): 
            g_pdf.SetPoint(i, centers_recoil[i], pdf_values[i]*yield_/pdf_frac_range)
            pdf_tot[i] += pdf_values[i]*yield_/pdf_frac_range
            if plotSignal:
                g_sig_pdf.SetPoint(i, centers_recoil[i], sig_pdf_values[i]*yield_/sig_pdf_frac_range)
                sig_pdf_tot[i] += sig_pdf_values[i]*yield_/sig_pdf_frac_range
        
        if plotSignal:
            g_sig_pdf.SetLineColor(ROOT.kBlue)
            g_sig_pdf.SetLineWidth(3)
            g_sig_pdf.Draw("L SAME")
        
        g_pdf.SetLineColor(ROOT.kRed)
        g_pdf.SetLineWidth(3)
        g_pdf.Draw("L SAME")

        histRatio = hist_root.Clone("ratio")
        histRatio.Reset("ACE")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        chi2, histRatio = ratioHist_tgraph(hist_root, histRatio, g_pdf)
        

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, procLabel)
        latex.DrawLatex(0.20, 0.80, metLabel)
        latex.DrawLatex(0.20, 0.75, "q_{T} = [%.1f, %.1f] GeV" % (qTlow, qThigh))
        latex.DrawLatex(0.20, 0.70, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.65, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.60, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.55, "#chi^{2}/ndof = %.3f" % chi2)

        
        padT.RedrawAxis()
        padT.RedrawAxis("G")
        plotter.auxRatio()
        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        padB.SetTickx()
        padB.SetTicky()
        dummyB.Draw("HIST")
        #if doUnc: hist_root_tot_unc_ratio.Draw("E4 SAME")
        dummyL.Draw("SAME")
        histRatio.Draw("SAME E0")
        padB.RedrawAxis()
        padB.RedrawAxis("G")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, iBin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, iBin))
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
        g_chi2.SetPoint(iBin-1, qT, chi2)
        
    plotChi2(g_chi2, "%s/chi2" % outDir_, procLabel, metLabel, xMin=min(binning_qT), xMax=max(binning_qT))
    
    
    # uncertainties
    if doUnc:
        nStatVars = fitCfg['nStatVars']
        hist_root_tot_unc = hist_root_tot.Clone("hist_root_tot_unc")
        hist_root_tot_unc.SetFillColor(18)
        hist_root_tot_unc.SetMarkerSize(0)
        hist_root_tot_unc.SetLineWidth(0)
        hist_root_tot_unc.Reset("ICE")
        
        hist_root_tot_unc_ratio = hist_root_tot.Clone("hist_root_tot_unc_ratio")
        hist_root_tot_unc_ratio.SetFillColor(18)
        hist_root_tot_unc_ratio.SetMarkerSize(0)
        hist_root_tot_unc_ratio.SetLineWidth(0)
        hist_root_tot_unc_ratio.Reset("ICE") # reset errors
        
        vx, vy = array.array('d'), array.array('d')
        exl, eyl, exh, eyh = array.array('d'), array.array('d'), array.array('d'), array.array('d')
        for k in range(0, len(centers_recoil)):
            vx.append(centers_recoil[k])
            vy.append(1)
            exl.append(0)
            eyl.append(0)
            exh.append(0)
            eyh.append(0)

        for nStat in range(0, nStatVars):
        
            #if nStat != 1: continue
            ### PLUS VARIATION
            pdf_tot_var_up = [0]*len(centers_recoil)
            nStatName = "stat%d_p" % nStat
            func_parms = []
            for k in range(0, fitCfg['nParams']):
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                    if nStat == 1: print("p", fitCfg[nStatName]["p%d"%k]["p%d"%l])
                    
            
            for iBin in range(1, len(binning_qT)):
                qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
                qT_tf = tf.constant(qT, dtype=tf.float64)
                
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms, args).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var_up[i] += pdf_values_var[i]*yields[iBin-1]/pdf_fracs[iBin-1] # iBin starts from 1
                     
            
            ### MINUS VARIATION
            pdf_tot_var_dw = [0]*len(centers_recoil)
            nStatName = "stat%d_m" % nStat
            func_parms = []
            for k in range(0, fitCfg['nParams']):
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                    if nStat == 1: print("m", fitCfg[nStatName]["p%d"%k]["p%d"%l])
            
            for iBin in range(1, len(binning_qT)):
                qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
                qT_tf = tf.constant(qT, dtype=tf.float64)
                
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms, args).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms, args).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var_dw[i] += pdf_values_var[i]*yields[iBin-1]/pdf_fracs[iBin-1] # iBin starts from 1     
               
            # make average of plus and minus variation  
            #pdf_tot_var = pdf_tot_var_up
            pdf_tot_var = []
            for i in range(0, len(centers_recoil)):
                pdf_tot_var.append(0.5*(pdf_tot_var_dw[i] + pdf_tot_var_up[i]))
                if nStat == 1:
                    print(i, centers_recoil[i], pdf_tot[i], pdf_tot_var_dw[i], pdf_tot_var_up[i])
            
            for i in range(0, len(centers_recoil)): 
                eyl[i] += (pdf_tot[i] - pdf_tot_var_dw[i])**2
                eyh[i] += (pdf_tot[i] - pdf_tot_var_up[i])**2
   
               
            for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1):
                pdf_tot_eval = pdf_tot[i-1]
                if pdf_tot_eval > 0:
                    err_up = abs(pdf_tot_eval - pdf_tot_var_up[i-1])
                    err_dw = abs(pdf_tot_eval - pdf_tot_var_dw[i-1])
                    err = 0.5*(err_up+err_dw) # max(err_up, err_dw) 
                    #err = 0.5*abs(pdf_tot_var_up[i-1] - pdf_tot_var_dw[i-1])
                    #err = err_up
                    #err = abs(pdf_tot_eval - pdf_tot_var[i-1])
                else:
                    err = 0
                #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinContent(i), pdf_eval_unc, pdf_eval, err)
                hist_root_tot_unc_ratio.SetBinError(i, hist_root_tot_unc_ratio.GetBinError(i) + err*err) # add in quadrature
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i) + err*err)
                #print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])


        for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1): 
            pdf_tot_eval = pdf_tot[i-1]
            if pdf_tot_eval > 0: 
                hist_root_tot_unc_ratio.SetBinError(i, (hist_root_tot_unc_ratio.GetBinError(i)**0.5)/pdf_tot_eval)
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i)**0.5)
            hist_root_tot_unc_ratio.SetBinContent(i, 1)
            hist_root_tot_unc.SetBinContent(i, pdf_tot_eval)
            #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinError(i), histUnc.GetBinError(i))
            #print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])
       
        for i in range(0, len(centers_recoil)): 
            eyl[i] = eyl[i]**0.5
            eyh[i] = eyh[i]**0.5
        hist_root_tot_unc_ratio_asym = ROOT.TGraphAsymmErrors(len(centers_recoil), vx, vy, exl, exh, eyl, eyh);
        hist_root_tot_unc_ratio_asym = hist_root_tot.Clone("hist_root_tot_unc_ratio")
        hist_root_tot_unc_ratio_asym.SetFillColor(18)
        hist_root_tot_unc_ratio_asym.SetMarkerSize(0)
        hist_root_tot_unc_ratio_asym.SetLineWidth(0)
            
    
    # global plot
    plotter.cfg = cfgPlot
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
    if doUnc: hist_root_tot_unc.Draw("E4 SAME")    
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.SetMarkerStyle(20)
    hist_root_tot.SetMarkerColor(ROOT.kBlack)
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.Draw("PE SAME")
    

    g_pdf = ROOT.TGraphErrors()
    g_sig_pdf = ROOT.TGraphErrors()
    for i in range(0, len(centers_recoil)): 
        g_pdf.SetPoint(i, centers_recoil[i], pdf_tot[i])
        if plotSignal:
            g_sig_pdf.SetPoint(i, centers_recoil[i], sig_pdf_tot[i])
    g_pdf.SetLineColor(ROOT.kRed)
    g_pdf.SetLineWidth(3)
    
    if plotSignal:
        g_sig_pdf.SetLineColor(ROOT.kBlue)
        g_sig_pdf.SetLineWidth(3)
        g_sig_pdf.Draw("L SAME")
    
    g_pdf.Draw("L SAME")
        
    histRatio = hist_root_tot.Clone("ratio")
    histRatio.Reset("ACE")
    histRatio.SetMarkerStyle(8)
    histRatio.SetMarkerSize(0.7)
    histRatio.SetMarkerColor(ROOT.kBlack)
    histRatio.SetLineColor(ROOT.kBlack)
    histRatio.SetLineWidth(1)
    chi2, histRatio = ratioHist_tgraph(hist_root_tot, histRatio, g_pdf)
        

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.85, procLabel)
    latex.DrawLatex(0.20, 0.80, metLabel)
    latex.DrawLatex(0.20, 0.75, "#chi^{2}/ndof = %.3f" % chi2)

    padT.RedrawAxis()
    padT.RedrawAxis("G")
    plotter.auxRatio()
    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    padB.SetTickx()
    padB.SetTicky()
    dummyB.Draw("HIST")
    if doUnc: hist_root_tot_unc_ratio.Draw("E4 SAME")
    #hist_root_tot_unc_ratio_asym.Draw("E4 SAME")
    dummyL.Draw("SAME")
    histRatio.Draw("SAME E0")
    padB.RedrawAxis()
    padB.RedrawAxis("G")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/global.png" % outDir_)
    canvas.SaveAs("%s/global.pdf" % outDir_)
    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()

 
def makeFunction(funcJs, i):

    cfg = funcJs["p%d"%i]
    func = ROOT.TF1(cfg['funcName'], getattr(rf, cfg['funcName']+"_")(), 0, 500)
    for j in range(0, cfg['nParams']):
        func.SetParameter(j, cfg["p%d"%j])
    return func
    
 
def doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_, binning_qT, yMin=1e-6, yMax=1e4, bkgCfg={}, funcJs=None, yRatio = 1.15, recoilLow=-100, recoilHigh=100, rebin=1):
    
    average_eval_PDF = False
 
    cfgPlot = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : recoilLow,
        'xmax'              : recoilHigh,
        'ymin'              : yMin,
        'ymax'              : yMax,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)" if comp == "para" else "Recoil U_{#perp}   (GeV)",
        'ytitle'            : "Events" ,
        
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }
    
    axis_name = "recoil_perp" if comp == "perp" else "recoil_para_qT"
    qTlow, qThigh = min(binning_qT), max(binning_qT)
    bhist =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)], axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]
    
    
    bins_qT = bhist.axes[0].edges
    bins_recoil = bhist.axes[1].edges
    centers_qT = bhist.axes[0].centers
    centers_recoil = bhist.axes[1].centers
    centers_recoil_tf = tf.constant(centers_recoil, dtype=tf.float64)
    bins_recoil_tf = tf.constant(bins_recoil, dtype=tf.float64)
    ##bins_widths = [bins_recoil[i]-bins_recoil[i-1] for i in range(1, len(bins_recoil))]
    
    print(len(centers_recoil))
    #quit()
    
    functions.prepareDir(outDir_, True)
    func_name = fitCfg['func_name']
    func_parms_vals = fitCfg['func_parms_vals']
    func_parms_cfg = fitCfg['func_parms_cfg']
    func_parms_func = [] # TF1s
    for i,cfg in enumerate(func_parms_cfg):
        if cfg != 2: func_parms_func.append(None)
        else: func_parms_func.append(makeFunction(funcJs, i))
        
    func_model = getattr(rf, func_name)
    outDict = {}
    outDict['func_name'] = func_name
    
    func_parms = func_parms_vals # initial conditions
    

    # backgrounds
    args = ()
    args_ = {}
    if len(bkgCfg) > 0:
        
        for i,bkg in enumerate(bkgCfg["procs"]):
            parms_cond = [] # construct array of conditional pdf parameters
            for k in range(0, bkgCfg['parms'][i]['nParams']):
                for l in range(0, bkgCfg['parms'][i]["p%d"%k]['nParams']):
                    parms_cond.append(tf.constant(bkgCfg['parms'][i]["p%d"%k]["p%d"%l], dtype=tf.float64))
            args_['parms_%s'%bkg] = parms_cond
            yield_fracs = [] # construct array of yield fractions
            for iBin in range(1, len(binning_qT)):
                yield_fracs.append(bkgCfg['norms'][i]*bkgCfg['yields'][i]["%d"%iBin]["yield"]/bkgCfg['data_yields']["%d"%iBin]['yield'])
            args_['fracs_%s'%bkg] = tf.constant(yield_fracs, dtype=tf.float64)
                
    args = (args_)

    g_chi2 = ROOT.TGraphErrors()
    hist_root_tot = None
    pdf_tot = [0]*len(centers_recoil)
    yields, yields_err = [], []
    for iBin in range(1, len(binning_qT)):
    
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        qT_tf = tf.constant(qT, dtype=tf.float64)
        qTlow, qThigh = binning_qT[iBin-1], binning_qT[iBin]
        print("############# Recoil bin %d [%.2f, %.2f]" % (iBin, qTlow, qThigh))

        # set the parameters
        for i,func in enumerate(func_parms_func):
            if func == None: continue
            func_parms[i] = func.Eval(qT)
            
        h = bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h = h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum  
        h = h[{axis_name: s[::bh.rebin(rebin)]}]
        yield_, yield_err = h.sum().value, math.sqrt(h.sum().variance)
        yields.append(yield_)
        yields_err.append(yield_err)
        
        args['qT'] = qT_tf # always propagate qT, in case it is used
            
        print(" -> Start fitting")
        start = time.time()
        res = fitter.fit_hist(h, func_model, np.array(func_parms), max_iter = 5, edmtol = 1e-5, mode = "nll", xargs=args, sumw2=True) # nll chisq_normalized
        end = time.time()
            
        parms_vals = res['x']
        fit_status = res['status']
        cov_status = res['covstatus']
        hess_eigvals = res['hess_eigvals']
        cov_matrix = res['cov']
        func_parms = [parms_vals[i] if func_parms_cfg[i]==1 else func_parms_vals[i] for i in range(0, len(func_parms))]
        #if propagate: func_parms = parms_vals
            
        print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
        print(parms_vals)
        

        if average_eval_PDF:
            pdf_values_bins = func_model([qT_tf], parms_vals, args).numpy()
            pdf_values = [0.5*(pdf_values_bins[i-1]+pdf_values_bins[i]) for i in range(1, len(pdf_values_bins))]      
        else:
            pdf_values = func_model([centers_recoil_tf], parms_vals, args).numpy()
            
        pdf_frac_range = sum(pdf_values)

        plotter.cfg = cfgPlot
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        
        hist_root = narf.hist_to_root(h) 
        hist_root.Scale(1.0, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")
        if hist_root_tot == None: hist_root_tot = hist_root.Clone("hist_root_tot")
        else: hist_root_tot.Add(hist_root)

        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers_recoil)): 
            g_pdf.SetPoint(i, centers_recoil[i], pdf_values[i]*yield_/pdf_frac_range)
            pdf_tot[i] += pdf_values[i]*yield_/pdf_frac_range
        g_pdf.SetLineColor(ROOT.kRed)
        g_pdf.SetLineWidth(3)
        g_pdf.Draw("L SAME")

        histRatio = hist_root.Clone("ratio")
        histRatio.Reset("ACE")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        chi2, histRatio = ratioHist_tgraph(hist_root, histRatio, g_pdf)
        

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, procLabel)
        latex.DrawLatex(0.20, 0.80, metLabel)
        latex.DrawLatex(0.20, 0.75, "q_{T} = [%.1f, %.1f] GeV" % (qTlow, qThigh))
        latex.DrawLatex(0.20, 0.70, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.65, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.60, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.55, "#chi^{2}/ndof = %.3f" % chi2)
        latex.DrawLatex(0.20, 0.50, "#color[%d]{Fit/Cov. status = %d/%d}" % (2 if fit_status!=0 or cov_status!=0 else 1, fit_status, cov_status))
 
        latex.SetTextSize(0.035)
        
        # save output
        outDict[iBin] = {}
        outDict[iBin]['yield'] = yield_
        outDict[iBin]['yield_err'] = yield_err

        for i in range(0, len(func_parms)):
            val, err = parms_vals[i], math.sqrt(cov_matrix[i][i]) if cov_status == 0 and cov_matrix[i][i] >= 0 else -2
            outDict[iBin]["p%d"%i], outDict[iBin]["p%d_err"%i] = val, err
            latex.DrawLatex(0.6, 0.86-i*0.04, "p_{%d} = %.2e #pm %.2e" % (i, val, err))
        
        plotter.auxRatio()
        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        padB.SetTickx()
        padB.SetTicky()
        dummyB.Draw("HIST")
        dummyL.Draw("SAME")
        histRatio.Draw("SAME E0")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, iBin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, iBin))
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
     

        outDict[iBin]['chi2'] = chi2
        g_chi2.SetPoint(iBin-1, qT, chi2)


    outDict['qTbins'] = binning_qT
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, procLabel, metLabel, xMin=min(binning_qT), xMax=max(binning_qT))
    
    
    
    # uncertainties
    doUnc = 'nStatVars' in fitCfg
    if doUnc:
        nStatVars = fitCfg['nStatVars']
        hist_root_tot_unc = hist_root_tot.Clone("hist_root_tot_unc")
        hist_root_tot_unc.SetFillColor(ROOT.kBlack)
        hist_root_tot_unc.SetMarkerSize(0)
        hist_root_tot_unc.SetLineWidth(0)
        hist_root_tot_unc.SetFillStyle(3004)    
        hist_root_tot_unc.Reset("ICE")
        
        hist_root_tot_unc_ratio = hist_root_tot.Clone("hist_root_tot_unc_ratio")
        hist_root_tot_unc_ratio.SetFillColor(ROOT.kBlack)
        hist_root_tot_unc_ratio.SetMarkerSize(0)
        hist_root_tot_unc_ratio.SetLineWidth(0)
        hist_root_tot_unc_ratio.SetFillStyle(3004)    
        hist_root_tot_unc_ratio.Reset("ICE") # reset errors

        for nStat in range(0, nStatVars):
            nStatName = "stat%d" % nStat
            #if nStat > 0: continue
            func_parms = []
            for k in range(0, fitCfg['nParams']):
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    func_parms.append(fitCfg[nStatName]["p%d"%k]["p%d"%l])
                    
            pdf_tot_var = [0]*len(centers_recoil)
            for iBin in range(1, len(binning_qT)):
            
                qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
                qT_tf = tf.constant(qT, dtype=tf.float64)
                
                if average_eval_PDF:
                    pdf_values_var_bins = func_model([qT_tf, bins_recoil_tf], func_parms).numpy()
                    pdf_values_var = [0.5*(pdf_values_var_bins[i-1]+pdf_values_var_bins[i]) for i in range(1, len(pdf_values_var_bins))]
                else:
                    pdf_values_var = func_model([qT_tf, centers_recoil_tf], func_parms).numpy()

                for i in range(0, len(centers_recoil)): 
                    pdf_tot_var[i] += pdf_values_var[i]*yields[iBin-1] # iBin starts from 1
                        
            for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1):
                pdf_tot_eval = pdf_tot[i-1]
                
                if pdf_tot_eval > 0: 
                    err = abs(pdf_tot[i-1] - pdf_tot_var[i-1])
                else:
                    err = 0
                #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinContent(i), pdf_eval_unc, pdf_eval, err)
                hist_root_tot_unc_ratio.SetBinError(i, hist_root_tot_unc_ratio.GetBinError(i) + err*err) # add in quadrature
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i) + err*err)
                #print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])
           # quit() 
                    

        for i in range(1, hist_root_tot_unc_ratio.GetNbinsX()+1): 
            pdf_tot_eval = pdf_tot[i-1]
            if pdf_tot_eval > 0: 
                hist_root_tot_unc_ratio.SetBinError(i, (hist_root_tot_unc_ratio.GetBinError(i)**0.5)/pdf_tot_eval)
                hist_root_tot_unc.SetBinError(i, hist_root_tot_unc.GetBinError(i)**0.5)
            hist_root_tot_unc_ratio.SetBinContent(i, 1)
            hist_root_tot_unc.SetBinContent(i, pdf_tot_eval)
            #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinError(i), histUnc.GetBinError(i))
            #if i == 100: print(i,hist_root_tot_unc_ratio.GetBinCenter(i), centers_recoil[i-1], pdf_tot[i-1], pdf_tot_var[i-1], err, pdf_tot[i-1]/pdf_tot_var[i-1])
            
            
            
    
    # global plot
    plotter.cfg = cfgPlot
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=1)
        
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
        
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.SetMarkerStyle(20)
    hist_root_tot.SetMarkerColor(ROOT.kBlack)
    hist_root_tot.SetLineColor(ROOT.kBlack)
    hist_root_tot.Draw("PE SAME")

    g_pdf = ROOT.TGraphErrors()
    for i in range(0, len(centers_recoil)):
        g_pdf.SetPoint(i, centers_recoil[i], pdf_tot[i])
    g_pdf.SetLineColor(ROOT.kRed)
    g_pdf.SetLineWidth(3)
    if doUnc: hist_root_tot_unc.Draw("E2 SAME") 
    g_pdf.Draw("L SAME")
        
    histRatio = hist_root_tot.Clone("ratio")
    histRatio.Reset("ACE")
    histRatio.SetMarkerStyle(8)
    histRatio.SetMarkerSize(0.7)
    histRatio.SetMarkerColor(ROOT.kBlack)
    histRatio.SetLineColor(ROOT.kBlack)
    histRatio.SetLineWidth(1)
    chi2, histRatio = ratioHist_tgraph(hist_root_tot, histRatio, g_pdf)
        

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.85, procLabel)
    latex.DrawLatex(0.20, 0.80, metLabel)
    latex.DrawLatex(0.20, 0.75, "#chi^{2}/ndof = %.3f" % chi2)

    plotter.auxRatio()
    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    padB.SetTickx()
    padB.SetTicky()
    dummyB.Draw("HIST")
    if doUnc: hist_root_tot_unc_ratio.Draw("E2 SAME")
    dummyL.Draw("SAME")
    histRatio.Draw("SAME E0")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/global.png" % outDir_)
    canvas.SaveAs("%s/global.pdf" % outDir_)
    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()



def doFitMultiGauss_scipy(fInName, proc, comp, fitCfg, label, outDir_, recoil_qTbins, ratio=True, funcJs=None, doFit=True, excludeBins=[], qTstartBin = 1, qTmin=0, qTmax=300, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, propagate=True, bkgCfg={}, orderSigmas=False, yRatio = 1.15):
    
    cfgPlot = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)" if comp == "para" else "Recoil U_{#perp}   (GeV)",
        'ytitle'            : "Events" ,
        
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio" if ratio else "Pull",
        'yminR'             : (1-(yRatio-1)) if ratio else -2.5,
        'ymaxR'             : yRatio if ratio else 2.5,
    }
    
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")
    hIn_y = hIn.ProjectionY("hIn_y")
    
    recoil_bins = [hIn_y.GetBinLowEdge(i) for i in range(1, hIn_y.GetNbinsX()+2)]
    centers = [hIn_y.GetBinCenter(i) for i in range(1, hIn_y.GetNbinsX()+1)]

    centers_tf = tf.constant(centers, dtype=tf.float64)
    
 
    functions.prepareDir(outDir_, True) # remove and recreate dir
    nGauss = len(fitCfg['mean'])
    dof_param = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(3)
    dof = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(0) + dof_param # if cfg=0 or 3, add DOF


    func_model = rf.func_4gauss
    parms_vals, parms_names = [], []
       
    # construct gauss
    for iGauss in range(1, nGauss+1):
        parms_vals.append(fitCfg['sigma'][iGauss-1])
        parms_names.append("sigma%d" % (iGauss))

    #for iGauss in range(1, nGauss):
    #    parms_vals.append(fitCfg['norm'][iGauss-1])
    #    parms_names.append("norm%d" % (iGauss))

    mean_ptrs, sigma_ptrs, norm_ptrs, gauss_ptrs = [], [], [], []
       
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()

   
         
    outDict = {}
    outDict['nGauss'] = nGauss
    recoil_qTbins_ = [0]
 
    if qTstartBin >= len(recoil_qTbins) or qTstartBin < 1: sys.exit("qTstartBin not in range (should be at least 1)")
    if recoil_qTbins[qTstartBin-1] < qTmin: sys.exit("qTstartBin is lower than the minimum range defined")
    if recoil_qTbins[qTstartBin] > qTmax: sys.exit("qTstartBin is higher than the maximum range defined")
    qTbin = qTstartBin
    direction = -1
    ###if qTstartBin == 1: direction = 1
    while True:

        if qTbin >= len(recoil_qTbins): break
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        qTmin_, qTmax_ = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        if qT > qTmax: break
        
        recoil_qTbins_.append(recoil_qTbins[qTbin])

        qTbins_ = list(functions.drange(qTmin_, qTmax_, 0.5))
        #qTbins_ = list(range(qTmin_, qTmax_, 1))
        hist = None
        for b in qTbins_:
            iBin = hIn_x.FindBin(b)
            h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
            if hist == None: 
                hist = copy.deepcopy(h)
                hist.SetName("hist")
            else: hist.Add(h)

        
        ##yield_err = ctypes.c_double(1.)
        ##yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        ##yield_err = yield_err.value

        h = narf.root_to_hist(hist, axis_names=["recoil_bins"]) 
        
        
        yield_ = h.sum().value
        yield_err = math.sqrt(h.sum().variance)
        
        if doFit:
        
            start = time.time()
            res = fitter.fit_hist(h, func_model, np.array(parms_vals), max_iter = 5, edmtol = 1e-5, mode = "chisq_normalized")
            end = time.time()
            
            parms_vals = res['x'] # update the parameters
            fit_status = res['status']
            cov_status = res['covstatus']
            hess_eigvals = res['hess_eigvals']
            cov_matrix = res['cov']
            
            print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
            pdf_values = func_model([centers_tf], parms_vals).numpy() # eager evaluation
         
        else:
        
            func = getattr(rf, funcJs['func_name'])
            nParams = funcJs['nParams']
            parms_vals = [funcJs['p%d'%iPar] for iPar in range(0, nParams)]
            print(parms_vals)
            
            pdf_values = func_model([centers_tf, qT], parms_vals).numpy() # eager evaluation, conditional pdf
            
            quit()
            
            cov_status = -1
            fit_status = -1
         
        if direction > 0 and qTbin == qTstartBin:
            qTbin += direction
            continue

        ## plot
        plotter.cfg = cfgPlot
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(line=1 if ratio else 0)
        
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        
        hist.SetLineColor(ROOT.kBlack)
        hist.SetMarkerStyle(20)
        hist.SetMarkerColor(ROOT.kBlack)
        hist.SetLineColor(ROOT.kBlack)
        hist.Draw("PE SAME")

        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers)): g_pdf.SetPoint(i, centers[i], pdf_values[i]*yield_)
        g_pdf.SetLineColor(ROOT.kBlue)
        g_pdf.SetLineWidth(3)
        g_pdf.Draw("L SAME")
        
        histRatio = hist.Clone("ratio")
        histRatio.Reset("ACE")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        chi2, histRatio = ratioHist_tgraph(hist, histRatio, g_pdf)
        

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, label)
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin_, qTmax_))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist.GetMean(), hist.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist.GetRMS(), hist.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        latex.DrawLatex(0.20, 0.55, "#color[%d]{Fit status = %d}" % (2 if fit_status != 0 else 1, fit_status))
        latex.DrawLatex(0.20, 0.50, "#color[%d]{Covariance status = %d}" % (2 if cov_status != 0 else 1, cov_status))
 
        latex.SetTextSize(0.035)
        
        # save output
        outDict[qTbin] = {}
        outDict[qTbin]['yield'] = yield_
        outDict[qTbin]['yield_err'] = yield_err
        

        # take remaining parameters as initial values
        nPar = 0
        for iGauss in range(1, nGauss+1):
            tag = 'mean%d' % iGauss
            if tag in parms_names: 
                idx = parms_names.index(tag)
                val, err = parms_vals[idx], math.sqrt(cov_matrix[idx][idx]) if doFit and cov_status == 0 and cov_matrix[idx][idx] >= 0 else -2
            else:
                val, err = fitCfg['mean'][iGauss-1], -1
            outDict[qTbin][tag], outDict[qTbin][tag+"_err"] = val, err
            latex.DrawLatex(0.7, 0.87-nPar*0.04, "#mu_{%d} = %.3f #pm %.3f" % (iGauss, val, err))
            nPar += 1
        for iGauss in range(1, nGauss+1):
            tag = 'sigma%d' % iGauss
            if tag in parms_names: 
                idx = parms_names.index(tag)
                val, err = parms_vals[idx], math.sqrt(cov_matrix[idx][idx]) if doFit and cov_status == 0 and cov_matrix[idx][idx] >= 0 else -2
            else:
                val, err = fitCfg['sigma'][iGauss-1], -1
            outDict[qTbin][tag], outDict[qTbin][tag+"_err"] = val, err
            latex.DrawLatex(0.7, 0.87-nPar*0.04, "#sigma_{%d} = %.3f #pm %.3f" % (iGauss, val, err))
            nPar += 1
        for iGauss in range(1, nGauss):
            tag = 'norm%d' % iGauss
            if tag in parms_names: 
                idx = parms_names.index(tag)
                val, err = parms_vals[idx], math.sqrt(cov_matrix[idx][idx]) if doFit and cov_status == 0 and cov_matrix[idx][idx] >= 0 else -2
            else:
                val, err = fitCfg['norm'][iGauss-1], -1
            outDict[qTbin][tag], outDict[qTbin][tag+"_err"] = val, err
            latex.DrawLatex(0.7, 0.87-nPar*0.04, "n_{%d} = %.3f #pm %.3f" % (iGauss, val, err))
            nPar += 1


        plotter.auxRatio()
        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        dummyB.Draw("HIST")
        histRatio.Draw("SAME E0")
        dummyL.Draw("SAME")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, qTbin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, qTbin))
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
     
        
        
        
        
        outDict[qTbin]['chi2'] = chi2
        g_chi2.SetPoint(qTbin-1, qT, chi2)
        g_yields.SetPoint(qTbin-1, qT, yield_)
        g_yields.SetPointError(qTbin-1, 0, yield_err)
        
        
        # update bin condition
        if qTmin_ <= qTmin or qTmin_ == 0:
            qTbin = qTstartBin - 1 # minus 1, redo the fit at qTstartBin for initial conditions
            direction = +1
        
        qTbin += direction # update bin 
        
        
        
    outDict['qTbins'] = recoil_qTbins_
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)

    with open("%s/results.json" % outDir_) as f: out = json.load(f)
    doGlobalPlot(fInName, proc, comp, out, recoil_qTbins, outDir_, label, funcJs=funcJs, ratio=ratio, bkgCfg=bkgCfg, qTmin=qTmin, qTmax=qTmax, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, yRatio=yRatio)
 



def doGlobalPlot(fInName, proc, comp, jsIn, binning_qT, outDir_, label, funcJs=None, bkgCfg={}, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, yRatio = 1.15, recoilLow=-100, recoilHigh=100):
    

    ptrs = {}
    garbage = [] # need to store the variables for memory issues
    means_arglist, sigmas_arglist, norms_arglist = [], [], []
    
    totYield = 0
    totHist = None
    for qTbin in range(1, len(recoil_qTbins)):
    
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        if qT > qTmax or qT < qTmin: continue
        totYield += jsIn[str(qTbin)]['yield']
    

    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")

    pdfs_arglist_tot, norms_arglist_tot = ROOT.RooArgList(), ROOT.RooArgList()
    pdfs_arglist_tot.setName("pdfs_arglist_tot")
    norms_arglist_tot.setName("norms_arglist_tot")
    for qTbin in range(1, len(recoil_qTbins)):
    
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        if qT > qTmax or qT < qTmin: continue
        qTmin_, qTmax_ = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        
        qTbins_ = list(functions.drange(qTmin_, qTmax_, 0.5))
        #qTbins_ = list(range(qTmin_, qTmax_, 1))
        hist = None
        for b in qTbins_:
            iBin = hIn_x.FindBin(b)
            h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
            if hist == None: 
                hist = copy.deepcopy(h)
                hist.SetName("hist")
            else: hist.Add(h)


        yield_ = hist.Integral()
        if totHist == None: totHist = hist.Clone("totHist")
        else: totHist.Add(hist)

        pdfs, norms = ROOT.RooArgList(), ROOT.RooArgList()
        pdfs.setName("pdfs_bin%d" % qTbin)
        norms.setName("norms_bin%d" % qTbin)
        for iGauss in range(1, nGauss+1):
            mean_ = ROOT.RooRealVar("mean%d_bin%d" % (iGauss, qTbin), "", jsIn[str(qTbin)]['mean%d' % iGauss], __min_mean__, __max_mean__)
            sigma_ = ROOT.RooRealVar("sigma%d_bin%d" % (iGauss, qTbin), "", jsIn[str(qTbin)]['sigma%d' % iGauss], __min_sigma__, __max_sigma__)
            gauss = ROOT.RooGaussian("gauss%d_bin%d" % (iGauss, qTbin), "", recoil, mean_, sigma_)
            pdfs.addOwned(gauss)
            garbage.append(mean_)
            garbage.append(sigma_)
            garbage.append(gauss)

        for iGauss in range(1, nGauss):
            norm_ = ROOT.RooRealVar("norm%d_bin%d" % (iGauss, qTbin), "", jsIn[str(qTbin)]['norm%d' % iGauss], 0, 1)
            norms.addOwned(norm_)
            garbage.append(norm_)
        
        garbage.append(pdfs)
        garbage.append(norms)
        pdf_sig = ROOT.RooAddPdf("pdf_sig_bin%d" % (qTbin), '', pdfs, norms)
        ptrs[pdf_sig.GetName()] = pdf_sig
        
        # include backgrounds (if configured)
        pdfs_arglist_, norms_arglist_ = ROOT.RooArgList(), ROOT.RooArgList()
        pdfs_arglist_.setName("pdfs_arglist_bin%d" % (qTbin))
        norms_arglist_.setName("norms_arglist_bin%d" % (qTbin))
        for bkg in bkgCfg:
            with open(bkgCfg[bkg]['filename']) as f: jsIn_bkg = json.load(f)
            pdfs_bkgs_arglist, norms_bkgs_arglist = ROOT.RooArgList(), ROOT.RooArgList()
            pdfs_bkgs_arglist.setName("pdfs_bkgs_arglist_bin%d_%s" % (qTbin, bkg))
            norms_bkgs_arglist.setName("norms_bkgs_arglist_bin%d_%s" % (qTbin, bkg))
            for iGauss in range(1, jsIn_bkg['nGauss']+1):
                mean = ROOT.RooRealVar("mean%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn_bkg[str(qTbin)]["mean%d" % iGauss], __min_mean__, __max_mean__)
                sigma = ROOT.RooRealVar("sigma%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn_bkg[str(qTbin)]["sigma%d" % iGauss], __min_sigma__, __max_sigma__)
                mean.setConstant(ROOT.kTRUE)
                sigma.setConstant(ROOT.kTRUE)
                gauss = ROOT.RooGaussian("gauss%d_bin%d_%s" % (iGauss, qTbin, bkg), "", recoil, mean, sigma)
                garbage.append(mean)
                garbage.append(sigma)
                garbage.append(gauss)
                pdfs_bkgs_arglist.addOwned(gauss)
            for iGauss in range(1, jsIn_bkg['nGauss']):
                norm = ROOT.RooRealVar("norm%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn_bkg[str(qTbin)]["norm%d" % iGauss], 0, 1)
                norm.setConstant(ROOT.kTRUE)
                norms_bkgs_arglist.addOwned(norm)
                garbage.append(norm)
            pdf_bkg = ROOT.RooAddPdf("bin%d_%s" % (qTbin, bkg), '', pdfs_bkgs_arglist, norms_bkgs_arglist)
            pdfs_arglist_.addOwned(pdf_bkg)
            norm_bkg = ROOT.RooRealVar("relyield_bin%d_%s" % (qTbin, bkg), "", (jsIn_bkg[str(qTbin)]['yield']*bkgCfg[bkg]['norm'])/yield_, 0, 1)
            if not bkgCfg[bkg]['float']: norm_bkg.setConstant(ROOT.kTRUE)
            norms_arglist_.addOwned(norm_bkg)
            garbage.append(pdfs_bkgs_arglist)
            garbage.append(norms_bkgs_arglist)
        
        garbage.append(pdfs_arglist_)
        garbage.append(norms_arglist_)
        pdfs_arglist_.addClone(pdf_sig) # append signal as the last one
        pdf_tot_ = ROOT.RooAddPdf("pdf_tot_bin%d" % (qTbin), '', pdfs_arglist_, norms_arglist_)
        garbage.append(pdf_tot_)
        ptrs[pdf_tot_.GetName()] = pdf_tot_
        
        norm_bkg = ROOT.RooRealVar("totyield_bin%d" % (qTbin), "", jsIn[str(qTbin)]['yield']/totYield)
        norms_arglist_tot.addOwned(norm_bkg)
        pdfs_arglist_tot.addOwned(pdf_tot_)
        garbage.append(norm_bkg)
        
        

    pdf_tot = ROOT.RooAddPdf("pdf_tot", '', pdfs_arglist_tot, norms_arglist_tot)
    rdh_tot = ROOT.RooDataHist("rdh_tot", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(totHist))
    
    doUnc = funcJs != None and 'nStatVars' in funcJs # uncertainties from the refit
    pdf_tot_vars = []
    if doUnc:
    
        histUnc = totHist.Clone("histUnc")
        histUnc.SetFillColor(ROOT.kBlack)
        histUnc.SetMarkerSize(0)
        histUnc.SetLineWidth(0)
        histUnc.SetFillStyle(3004)    
        histUnc.Reset("ICE")
        
        histRatioUnc = totHist.Clone("histRatioUnc")
        histRatioUnc.SetFillColor(ROOT.kBlack)
        histRatioUnc.SetMarkerSize(0)
        histRatioUnc.SetLineWidth(0)
        histRatioUnc.SetFillStyle(3004)    
        histRatioUnc.Reset("ICE")

        nStatVars = funcJs['nStatVars']
        for nStat in range(0, nStatVars):
            nStatName = "stat%d" % nStat
            #if nStat != 0: continue
            # get the functions
            funcs_mean, funcs_sigma, funcs_norm = [], [], []
            for iGauss in range(1, nGauss+1):
                func_mean = ROOT.TF1("mean%d_%s" % (iGauss, nStatName), funcJs['mean%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['mean%d' % iGauss]['nParams']): func_mean.SetParameter(iParam, funcJs[nStatName]['mean%d' % iGauss]['p%d' % iParam])
                funcs_mean.append(func_mean)
                func_sigma = ROOT.TF1("sigma%d_%s" % (iGauss, nStatName), funcJs['sigma%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['sigma%d' % iGauss]['nParams']): func_sigma.SetParameter(iParam, funcJs[nStatName]['sigma%d' % iGauss]['p%d' % iParam])
                funcs_sigma.append(func_sigma)
            for iGauss in range(1, nGauss):
                func_norm = ROOT.TF1("norm%d_%s" % (iGauss, nStatName), funcJs['norm%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['norm%d' % iGauss]['nParams']): func_norm.SetParameter(iParam, funcJs[nStatName]['norm%d' % iGauss]['p%d' % iParam])
                funcs_norm.append(func_norm)
                    
            pdfs_arglist_tot_var, norms_arglist_tot_var = ROOT.RooArgList(), ROOT.RooArgList()
            pdfs_arglist_tot_var.setName("pdfs_arglist_tot_%s" % nStatName)
            norms_arglist_tot_var.setName("norms_arglist_tot_%s" % nStatName)
            for qTbin in range(1, len(recoil_qTbins)):
            
                qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
                if qT > qTmax or qT < qTmin: continue
                qTmin_, qTmax_ = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
                
                pdfs, norms = ROOT.RooArgList(), ROOT.RooArgList()
                pdfs.setName("pdfs_bin%d_%s" % (qTbin, nStatName))
                norms.setName("norms_bin%d_%s" % (qTbin, nStatName))
                for iGauss in range(1, nGauss+1):
                    mean_ = ROOT.RooRealVar("mean%d_bin%d_%s" % (iGauss, qTbin, nStatName), "", funcs_mean[iGauss-1].Eval(qT), __min_mean__, __max_mean__)
                    sigma_ = ROOT.RooRealVar("sigma%d_bin%d_%s" % (iGauss, qTbin, nStatName), "", funcs_sigma[iGauss-1].Eval(qT), __min_sigma__, __max_sigma__)
                    gauss = ROOT.RooGaussian("gauss%d_bin%d_%s" % (iGauss, qTbin, nStatName), "", recoil, mean_, sigma_)
                    pdfs.addOwned(gauss)
                    garbage.append(mean_)
                    garbage.append(sigma_)
                    garbage.append(gauss)

                for iGauss in range(1, nGauss):
                    norm_ = ROOT.RooRealVar("norm%d_bin%d_%s" % (iGauss, qTbin, nStatName), "", funcs_norm[iGauss-1].Eval(qT), 0, 1)
                    norms.addOwned(norm_)
                    garbage.append(norm_)

                garbage.append(pdfs)
                garbage.append(norms)
                pdf_sig_var = ROOT.RooAddPdf("pdf_sig_bin%d_%s" % (qTbin, nStatName), '', pdfs, norms)
                ptrs[pdf_sig_var.GetName()] = pdf_sig_var
                
                # include backgrounds (if configured)
                pdfs_arglist_var_, norms_arglist_var_ = ROOT.RooArgList(), ROOT.RooArgList()
                pdfs_arglist_var_.setName("pdfs_arglist_bin%d_%s" % (qTbin, nStatName))
                norms_arglist_var_.setName("norms_arglist_bin%d_%s" % (qTbin, nStatName))
                for bkg in bkgCfg:
                    with open(bkgCfg[bkg]['filename']) as f: jsIn_bkg = json.load(f)
                    pdfs_bkgs_arglist, norms_bkgs_arglist = ROOT.RooArgList(), ROOT.RooArgList()
                    pdfs_bkgs_arglist.setName("pdfs_bkgs_arglist_bin%d_%s_%s" % (qTbin, bkg, nStatName))
                    norms_bkgs_arglist.setName("norms_bkgs_arglist_bin%d_%s_%s" % (qTbin, bkg, nStatName))
                    for iGauss in range(1, jsIn_bkg['nGauss']+1):
                        mean = ROOT.RooRealVar("mean%d_bin%d_%s_%s" % (iGauss, qTbin, bkg, nStatName), "", jsIn_bkg[str(qTbin)]["mean%d" % iGauss], __min_mean__, __max_mean__)
                        sigma = ROOT.RooRealVar("sigma%d_bin%d_%s_%s" % (iGauss, qTbin, bkg, nStatName), "", jsIn_bkg[str(qTbin)]["sigma%d" % iGauss], __min_sigma__, __max_sigma__)
                        mean.setConstant(ROOT.kTRUE)
                        sigma.setConstant(ROOT.kTRUE)
                        gauss = ROOT.RooGaussian("gauss%d_bin%d_%s_%s" % (iGauss, qTbin, bkg, nStatName), "", recoil, mean, sigma)
                        garbage.append(mean)
                        garbage.append(sigma)
                        garbage.append(gauss)
                        pdfs_bkgs_arglist.addOwned(gauss)
                    for iGauss in range(1, jsIn_bkg['nGauss']):
                        norm = ROOT.RooRealVar("norm%d_bin%d_%s_%s" % (iGauss, qTbin, bkg, nStatName), "", jsIn_bkg[str(qTbin)]["norm%d" % iGauss], 0, 1)
                        norm.setConstant(ROOT.kTRUE)
                        norms_bkgs_arglist.addOwned(norm)
                        garbage.append(norm)
                    pdf_bkg = ROOT.RooAddPdf("bin%d_%s_%s" % (qTbin, bkg, nStatName), '', pdfs_bkgs_arglist, norms_bkgs_arglist)
                    pdfs_arglist_.addOwned(pdf_bkg)
                    norm_bkg = ROOT.RooRealVar("relyield_bin%d_%s_%s" % (qTbin, bkg, nStatName), "", (jsIn_bkg[str(qTbin)]['yield']*bkgCfg[bkg]['norm'])/yield_, 0, 1)
                    if not bkgCfg[bkg]['float']: norm_bkg.setConstant(ROOT.kTRUE)
                    norms_arglist_.addOwned(norm_bkg)
                    garbage.append(pdfs_bkgs_arglist)
                    garbage.append(norms_bkgs_arglist)
                
                garbage.append(pdfs_arglist_var_)
                garbage.append(norms_arglist_var_)
                pdfs_arglist_var_.addClone(pdf_sig_var) # append signal as the last one
                pdf_tot_var_ = ROOT.RooAddPdf("pdf_tot_bin%d_%s" % (qTbin, nStatName), '', pdfs_arglist_var_, norms_arglist_var_)
                garbage.append(pdf_tot_var_)
                ptrs[pdf_tot_var_.GetName()] = pdf_tot_var_
                
                norm_bkg_var = ROOT.RooRealVar("totyield_bin%d_%s" % (qTbin, nStatName), "", jsIn[str(qTbin)]['yield']/totYield)
                norms_arglist_tot_var.addOwned(norm_bkg_var)
                pdfs_arglist_tot_var.addOwned(pdf_tot_var_)
                garbage.append(norm_bkg_var)
                
            pdf_tot_var = ROOT.RooAddPdf("pdf_tot_%s" % nStatName, '', pdfs_arglist_tot_var, norms_arglist_tot_var)
            pdf_tot_vars.append(pdf_tot_var)
            garbage.append(pdf_tot_var) 
            
            tmp = ROOT.RooArgSet(recoil)
            for i in range(1, histRatioUnc.GetNbinsX()+1):

                recoil.setVal(histRatioUnc.GetBinCenter(i))
                pdf_eval_unc = pdf_tot_var.getVal(tmp)
                pdf_eval = pdf_tot.getVal(tmp)
                if(pdf_eval > 0): 
                    err = abs(pdf_eval_unc - pdf_eval) 
                    #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinContent(i), pdf_eval_unc, pdf_eval, err)
                    histRatioUnc.SetBinError(i, histRatioUnc.GetBinError(i) + err*err)
                    histUnc.SetBinError(i, histUnc.GetBinError(i) + err*err)
                else: err = 0
    
        tmp = ROOT.RooArgSet(recoil)
        for i in range(1, histRatioUnc.GetNbinsX()+1): 
            recoil.setVal(histRatioUnc.GetBinCenter(i))
            pdf_eval = pdf_tot.getVal(tmp)
            if pdf_eval > 0: 
                histRatioUnc.SetBinError(i, histRatioUnc.GetBinError(i)**0.5/pdf_eval)
                histUnc.SetBinError(i, histUnc.GetBinError(i)**0.5)
            histRatioUnc.SetBinContent(i, 1)
            histUnc.SetBinContent(i, pdf_eval)
            #print(i, histRatioUnc.GetBinCenter(i), histRatioUnc.GetBinError(i), histUnc.GetBinError(i))
     

    cfgPlot = {

        'logy'              : True,
        'logx'              : False,
    
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
        
        'xtitle'            : "Recoil U_{#parallel} (GeV)" if comp == "para" else "Recoil U_{#perp}   (GeV)",
        'ytitle'            : "Events" ,
        
        'topRight'          : __topRight__,
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
        
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio" if ratio else "Pull",
        'yminR'             : (1-(yRatio-1)) if ratio else -2.5,
        'ymaxR'             : yRatio if ratio else 2.5,
    }
    
    
    # plot
    plotter.cfg = cfgPlot
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(line=1 if ratio else 0)
        
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetTickx()
    padT.SetTicky()
    padT.SetGrid()
    dummyT.Draw("HIST")
        
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    plt = recoil.frame()
    rdh_tot.plotOn(plt)
    if doUnc: histUnc.Draw("E2 SAME")
    pdf_tot.plotOn(plt) # ROOT.RooFit.Normalization(1, ROOT.RooAbsReal.NumEvent)
    
    #pdf_tot_vars[0].plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
    histpull = plt.pullHist()
    chi2 = plt.chiSquare()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.85, label)
    latex.DrawLatex(0.20, 0.80, "#chi^{2} = %.3f" % chi2)

    plt.Draw("SAME")
    plotter.auxRatio()

    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    dummyB.Draw("HIST")

    if ratio:
        histRatio = totHist.Clone("ratio")
        histRatio.SetMarkerStyle(8)
        histRatio.SetMarkerSize(0.7)
        histRatio.SetMarkerColor(ROOT.kBlack)
        histRatio.SetLineColor(ROOT.kBlack)
        histRatio.SetLineWidth(1)
        norm = ROOT.RooArgSet(recoil)
        histRatio = ratioHist(totHist, histRatio, pdf_tot, recoil, norm)
        histRatio.Draw("SAME E0")
        if doUnc: histRatioUnc.Draw("E20 SAME")
              
    else:
        plt = recoil.frame()
        plt.addPlotable(histpull, "P")
        plt.Draw("SAME")

    dummyL.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/global.png" % outDir_)
    canvas.SaveAs("%s/global.pdf" % outDir_)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)   
    dummyB.Delete()
    dummyT.Delete()
    padT.Delete()
    padB.Delete()    


def plotChi2(g_chi2, fOut, procLabel, metLabel, xMin=0, xMax=100, yMin=0, yMax=3):

    cfg_chi2 = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "#chi^{2}/ndof",
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    
    g_chi2.SetLineColor(ROOT.kRed)
    g_chi2.SetMarkerStyle(8)
    g_chi2.SetMarkerSize(1)
    g_chi2.SetMarkerColor(ROOT.kRed)
    
    f_chi2 = ROOT.TF1("f_chi2", "[0]*x + [1]", 0, 30)
    f_chi2.SetParameters(0, 1)
    f_chi2.SetLineColor(ROOT.kBlue)
    f_chi2.SetLineWidth(2)
    g_chi2.Fit("f_chi2", "NSE", "", 0, 30)
    
    plotter.cfg = cfg_chi2
    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTickx()
    canvas.SetTicky()
    dummy = plotter.dummy()
    canvas.cd()
    dummy.Draw("HIST")
    g_chi2.Draw("SAME LP")
    f_chi2.Draw("SAME L")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.88, procLabel)
    latex.DrawLatex(0.20, 0.84, metLabel)
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s.png" % fOut)
    canvas.SaveAs("%s.pdf" % fOut)
    canvas.Delete()


def ratioHist(hNom, hRatio, pdf, recoil, norm):

    hNom = hNom.Clone("test")
    hNom.Scale(1./hNom.Integral())
    for i in range(1, hRatio.GetNbinsX()+1):
            
        recoil.setVal(hRatio.GetBinLowEdge(i))
        pdfLow = pdf.getVal(norm)
        recoil.setVal(hRatio.GetBinLowEdge(i+1))
        pdfHigh = pdf.getVal(norm)
        pdf_eval = 0.5*(pdfLow + pdfHigh)
        y, y_err = 0, 0
        binWidth = hRatio.GetBinLowEdge(i+1) - hRatio.GetBinLowEdge(i)
        if(pdf_eval > 0): 
            y = hNom.GetBinContent(i)/pdf_eval/binWidth
            y_err = hNom.GetBinError(i)/hNom.GetBinContent(i)/binWidth if hNom.GetBinContent(i) > 0 else 0
        hRatio.SetBinContent(i, y)
        hRatio.SetBinError(i, y_err)
        #print(i, hRatio.GetBinLowEdge(i), hRatio.GetBinCenter(i), hRatio.GetBinLowEdge(i+1), hRatio.GetBinContent(i), hNom.GetBinContent(i), pdfLow, pdfHigh, pdf_eval)
    return hRatio

def ratioHist_tgraph(hNom, hRatio, pdf):

    chi2 = 0
    nBins = 0
    for i in range(1, hRatio.GetNbinsX()+1):
            
        xLow , xHigh = hRatio.GetBinLowEdge(i), hRatio.GetBinLowEdge(i+1)
        pdfLow, pdfHigh = pdf.Eval(xLow), pdf.Eval(xHigh)
        pdfVal = 0.5*(pdfLow + pdfHigh)
        pdfVal = pdf.Eval(hRatio.GetBinCenter(i))
        y, y_err = 0, 0
        if(pdfVal > 0): 
            y = hNom.GetBinContent(i)/pdfVal
            y_err = hNom.GetBinError(i)/hNom.GetBinContent(i) if hNom.GetBinContent(i) > 0 else 0
        hRatio.SetBinContent(i, y)
        hRatio.SetBinError(i, y_err)
    
        if hNom.GetBinContent(i) > 0:
            chi2 += (((pdfVal - hNom.GetBinContent(i)) / hNom.GetBinError(i)) ** 2)
            nBins += 1
            
    chi2 /= nBins      
    return chi2, hRatio
    

def plotYields(g_yields, fOut, procLabel, metLabel, xMin=0, xMax=100, yMin=0, yMax=3):

    cfg_chi2 = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 0,
        'ymax'              : 1.3*max(g_yields.GetY()),
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "Event yield",
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    
    g_yields.SetLineColor(ROOT.kRed)
    g_yields.SetMarkerStyle(8)
    g_yields.SetMarkerSize(1)
    g_yields.SetMarkerColor(ROOT.kRed)
    
    
    plotter.cfg = cfg_chi2
    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTickx()
    canvas.SetTicky()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    g_yields.Draw("SAME LP")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.88, procLabel)
    latex.DrawLatex(0.20, 0.84, metLabel)
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s.png" % fOut)
    canvas.SaveAs("%s.pdf" % fOut)
    canvas.Delete()






def plotParameter(param, jsIn, outDir, binning_qT, procLabel, metLabel, yMin=0, yMax=30, yTitle="", yRatio=1.15):

    xMin, xMax = min(binning_qT), max(binning_qT)

    f_base = ROOT.TF1("f_base", jsIn[param]['func'], 0, 500)
    for i in range(0, jsIn[param]['nParams']): f_base.SetParameter(i, jsIn[param]['p%d' % i])
    f_unc = ROOT.TF1("f_unc", jsIn[param]['func'], 0, 500)
    for i in range(0, jsIn[param]['nParams']): f_unc.SetParameter(i, jsIn[param]['p%d' % i])
    
    f_base.SetLineColor(ROOT.kRed)
    f_base.SetLineWidth(3)

    g = ROOT.TGraphErrors()
    #g.SetLineColor(ROOT.kBlack)
    #g.SetMarkerStyle(20)
    #g.SetMarkerSize(0.9)
    #g.SetMarkerColor(ROOT.kBlack)
    #g.SetLineColor(ROOT.kBlack)
    
           
    iPoint = 0
    for qTbin in range(1, len(binning_qT)):
        x = 0.5*(binning_qT[qTbin-1] + binning_qT[qTbin])
        x_err = 0.
        y = f_base.Eval(x)
        
        # loop over all stat uncertainties
        y_err = 0
        y_err_up, y_err_dw = 0, 0
        if 'nStatVars' in jsIn:
            for iStat in range(0, jsIn['nStatVars']):
                for i in range(0, jsIn[param]['nParams']): 
                    if param in jsIn['stat%d_p' % iStat]: f_unc.SetParameter(i, jsIn['stat%d_p' % iStat][param]['p%d' % i])
                    else : f_unc.SetParameter(i, jsIn[param]['p%d' % i])
                up = f_unc.Eval(x)
                s = up - y
                y_err_up += s**2
                
            for iStat in range(0, jsIn['nStatVars']):
                for i in range(0, jsIn[param]['nParams']): 
                    if param in jsIn['stat%d_m' % iStat]: f_unc.SetParameter(i, jsIn['stat%d_m' % iStat][param]['p%d' % i])
                    else : f_unc.SetParameter(i, jsIn[param]['p%d' % i])
                up = f_unc.Eval(x)
                s = up - y
                y_err_dw += s**2    
          
        y_err_dw = y_err_dw**(0.5)  
        y_err_up =y_err_up**(0.5)  
        y_err = 0.5*(y_err_dw + y_err_up) # y_err**(0.5)    
        #print(x, y, y_err)
        
        g.SetPoint(iPoint, x, y)   
        g.SetPointError(iPoint, 0, y_err)
        #print(x, y)
        iPoint += 1


    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : yTitle,
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
            
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }         

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()

    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
  
    
    g.SetFillColor(18)
    g.Draw("3SAME")
    f_base.Draw("L SAME")
    padT.RedrawAxis()
    padT.RedrawAxis("G")
    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.2, 0.85, procLabel)
    latex.DrawLatex(0.2, 0.80, metLabel)

    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    padB.SetTickx()
    padB.SetTicky()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")    

    padB.RedrawAxis()
    padB.RedrawAxis("G")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s.png" % (outDir, param))
    canvas.SaveAs("%s/%s.pdf" % (outDir, param))
    canvas.Delete()
    
    
    
def parameterizeGauss(jsIn, jsOut, comp, param, fitFuncName, iParams, outDir, binning_qT, procLabel, metLabel, cParams=[], excludeBins=[], fitMin=0, fitMax=200, xMin=0, xMax=150, yMin=0, yMax=30, yTitle="", bnds=None, doFit=True, cutOffMin=-99999, cutOffMax=99999, yRatio = 1.3, fitOpts="NS"):

    xMin, xMax = min(binning_qT), max(binning_qT)

    fitF = getattr(rf, fitFuncName + "_")()
    fit = ROOT.TF1("fit", fitF, xMin, xMax) if fitF != "" else None
    if cParams == []: cParams = [False]*len(iParams)
    for i, iPar in enumerate(iParams):
        if cParams[i]: fit.FixParameter(i, iPar)
        else: fit.SetParameter(i, iPar)


    g = ROOT.TGraphErrors()
    g.SetLineColor(ROOT.kBlack)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.55)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
           
    iPoint = 0
    for qTbin in range(1, len(binning_qT)):
        if not str(qTbin) in jsIn:
            continue
        if qTbin in excludeBins: continue
        x = 0.5*(binning_qT[qTbin-1] + binning_qT[qTbin])
        x_err = 0.5*(binning_qT[qTbin] - binning_qT[qTbin-1])
        if not param in jsIn[str(qTbin)]:
            print("WARNING: param %s not found for qTbin %d" % (param, qTbin))
            continue
        y = jsIn[str(qTbin)][param]
        y_err = jsIn[str(qTbin)][param + "_err"]
        if y_err < 0: y_err = 0
        
        #if x < fitMin or x > fitMax: continue
        if y < cutOffMin: continue
        if y > cutOffMax: continue
        #if y_err <= 0: continue
        g.SetPoint(iPoint, x, y)   
        g.SetPointError(iPoint, x_err, y_err)
        iPoint += 1

    if doFit and fit:
        result = g.Fit(fit.GetName(), fitOpts, "", fitMin, fitMax) 
        chisq = fit.GetChisquare()/fit.GetNDF()
        fit.SetLineColor(ROOT.kRed)
        fit.SetLineWidth(2)

            
        cov = result.GetCorrelationMatrix()       
        values = result.GetConfidenceIntervals(0.68, False)
        uncBand = ROOT.TGraphErrors()
        uncBand_ratio = ROOT.TGraphErrors()
        for i in range(len(values)):
            uncBand.SetPoint(i, g.GetX()[i], fit.Eval(g.GetX()[i]))
            uncBand.SetPointError(i, 0, values[i])
            uncBand_ratio.SetPoint(i, g.GetX()[i], 1)
            uncBand_ratio.SetPointError(i, 0, values[i])

        g_ratio = ROOT.TGraphErrors()
        g_ratio.SetLineColor(ROOT.kBlack)
        g_ratio.SetMarkerStyle(20)
        g_ratio.SetMarkerSize(0.55)
        g_ratio.SetMarkerColor(ROOT.kBlack)
        g_ratio.SetLineColor(ROOT.kBlack)
        
        iPoint = 0
        for qTbin in range(1, len(binning_qT)):
            if not str(qTbin) in jsIn:
                continue
            if qTbin in excludeBins: continue
            x = 0.5*(binning_qT[qTbin-1] + binning_qT[qTbin])
            x_err = 0.5*(binning_qT[qTbin] - binning_qT[qTbin-1])
            if not param in jsIn[str(qTbin)]: continue
            y = jsIn[str(qTbin)][param]
            y_err = jsIn[str(qTbin)][param + "_err"]
            y_fit = fit.Eval(x)
            
            #if x < fitMin or x > fitMax: continue
            if y < cutOffMin: continue
            if y > cutOffMax: continue
            
            if y_err > 0 and abs(y) != 0 and y_fit > 0: 
                ratio = y / y_fit
                ratio_err = y_err/abs(y)
            else:
                ratio = 1.0
                ratio_err = 0.0
            g_ratio.SetPoint(iPoint, x, ratio)
            g_ratio.SetPointError(iPoint, 0, ratio_err)
            iPoint += 1
          

    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : yTitle,
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
            
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }         

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")

    if doFit and fit:
        uncBand.SetFillColor(18)
        uncBand.Draw("3SAME")
    g.Draw("PE SAME")
    fit.Draw("L SAME")  
    padT.RedrawAxis()
    padT.RedrawAxis("G")
    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.2, 0.85, procLabel)
    latex.DrawLatex(0.2, 0.80, metLabel)
    latex.DrawLatex(0.2, 0.75, "#chi^{2}/ndof = %.2f" % chisq)
    
    latex.SetTextSize(0.035)
    for i in range(0, len(iParams)):
        latex.DrawLatex(0.6, 0.86-i*0.04, "a_{%d} = %.2e #pm %.2e" % (i, fit.GetParameter(i), fit.GetParError(i)))
        
       
        
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    padB.SetTickx()
    padB.SetTicky()
    dummyB.Draw("HIST")
        
    if doFit and fit:
        uncBand_ratio.SetFillColor(18)
        uncBand_ratio.Draw("3SAME")
    dummyL.Draw("SAME")
    g_ratio.Draw("PE0 SAME")
    
    padB.RedrawAxis()
    padB.RedrawAxis("G")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s.png" % (outDir, param))
    canvas.SaveAs("%s/%s.pdf" % (outDir, param))
    canvas.Delete()
    
    jsOut[param] = {}
    jsOut[param]["func"] = fitF
    jsOut[param]["funcName"] = fitFuncName
    jsOut[param]["nParams"] = len(iParams)
    for i, iPar in enumerate(iParams):
        jsOut[param]["p%d" % i] = fit.GetParameter(i)
        jsOut[param]["p%d_err" % i] = fit.GetParError(i)
        jsOut[param]["p%d_isCte" % i] = cParams[i]
        jsOut[param]["p%d_bnds" % i] = (None, None) if bnds==None else bnds[i]
    
    
  
 






def plotCovarianceMatrix(cov, outDir):

    import matplotlib.pyplot as plt

    s = len(cov)    
    fig, ax1 = plt.subplots(figsize=(s, s))
    pos = ax1.imshow(cov, cmap='Blues', interpolation='nearest')
    fig.colorbar(pos, ax=ax1)
    
    plt.savefig("%s/cov_matrix.png" % outDir)
    plt.savefig("%s/cov_matrix.pdf" % outDir)



    def correlation_from_covariance(covariance):
        v = np.sqrt(np.diag(covariance))
        outer_v = np.outer(v, v)
        correlation = covariance / outer_v
        correlation[covariance == 0] = 0
        return correlation

    corr = correlation_from_covariance(cov)
    fig, ax1 = plt.subplots(figsize=(s, s))
    pos = ax1.imshow(corr, cmap='Blues', interpolation='nearest', vmin=-1, vmax=1)
    fig.colorbar(pos, ax=ax1)
    
    plt.savefig("%s/corr_matrix.png" % outDir)
    plt.savefig("%s/corr_matrix.pdf" % outDir)


def combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg = {}, recoilLow=-100, recoilHigh=100, chisq_refit=False, rebin=1, outDir=""):
    
    # TODO: implement rebinning of binning_qT
    # TODO: what if fit failed, not invertible...
    # TODO: overflow?!
   
    axis_name = "recoil_perp" if comp == "perp" else "recoil_para_qT" # 
    qTlow, qThigh = min(binning_qT), max(binning_qT)
    h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)], axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]  #"recoil_perp": s[complex(0,recoilLow):complex(0,recoilHigh)
    h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
    h = h[{axis_name: s[::bh.rebin(rebin)]}]

    bins_qT = h.axes[0].edges
    bins_recoil = h.axes[1].edges
    centers_qT = h.axes[0].centers
    centers_recoil = h.axes[1].centers
    centers_recoil_tf = tf.constant(centers_recoil, dtype=tf.float64)

    func_name = fitCfg['func_name']
    print("Set initial parameters:")
    func_parms = []
    bnds = []
    for k in range(0, fitCfg['nParams']):
        for l in range(0, fitCfg["p%d"%k]['nParams']):
            func_parms.append(fitCfg["p%d"%k]["p%d"%l])
            bnds.append(fitCfg["p%d"%k]["p%d_bnds"%l] if "p%d_bnds"%l in fitCfg["p%d"%k] else (None, None))
            print("  -> p%d, p%d: %.2f" % (k, l, fitCfg["p%d"%k]["p%d"%l]))
    func_model = getattr(rf, func_name)

   
 
    # backgrounds
    args = ()
    if len(bkgCfg) > 0:
        args_ = {}
        for i,bkg in enumerate(bkgCfg["procs"]):
            parms_cond = [] # construct array of conditional pdf parameters
            for k in range(0, bkgCfg['parms'][i]['nParams']):
                for l in range(0, bkgCfg['parms'][i]["p%d"%k]['nParams']):
                    parms_cond.append(tf.constant(bkgCfg['parms'][i]["p%d"%k]["p%d"%l], dtype=tf.float64))
            args_['parms_%s'%bkg] = parms_cond
            yield_fracs = [] # construct array of yield fractions
            for iBin in range(1, len(binning_qT)):
                yield_fracs.append(bkgCfg['norms'][i]*bkgCfg['yields'][i]["%d"%iBin]["yield"]/bkgCfg['data_yields']["%d"%iBin]['yield'])
            args_['fracs_%s'%bkg] = tf.constant(yield_fracs, dtype=tf.float64)
        args = (args_)

    outDict = {}
    outDict['nParams'] = fitCfg['nParams']
    outDict['func_name'] = fitCfg['func_name']

               
    print(" -> Start NLL fitting")
    start = time.time()
    res = fitter.fit_hist(h, func_model, np.array(func_parms), max_iter = 5, edmtol = 1e-5, mode = "nll", xargs=args, bnds=bnds, norm_axes=[1], sumw2=True) # nll chisq_normalized
    end = time.time()
            
    parms_vals = res['x']
    fit_status = res['status']
    cov_status = res['covstatus']
    hess_eigvals = res['hess_eigvals']
    cov_matrix = res['cov']
    
    # plot cov matrix
    plotCovarianceMatrix(cov_matrix, outDir)

    print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
    print(parms_vals)
    print("cov")
    print(cov_matrix)
    
    if chisq_refit:
        print(" -> Start CHISQ fitting")
        start = time.time()
        res = fitter.fit_hist(h, func_model, np.array(parms_vals), max_iter = 5, edmtol = 1e-5, mode = "chisq_normalized", xargs=args) # nll chisq_normalized
        end = time.time()
                
        parms_vals = res['x']
        fit_status = res['status']
        cov_status = res['covstatus']
        hess_eigvals = res['hess_eigvals']
        cov_matrix = res['cov']
            
        print(" -> Fit ended, time=%.3f s, status=%d, cov_status=%d" % (end-start, fit_status, cov_status))
        print(parms_vals)
    
    # save output parms to dict
    iiter = 0
    for k in range(0, fitCfg['nParams']):
        outDict["p%d"%k] = {}
        outDict["p%d"%k]['nParams'] = fitCfg["p%d"%k]['nParams']
        outDict["p%d"%k]['func'] = fitCfg["p%d"%k]['func']
        for l in range(0, fitCfg["p%d"%k]['nParams']):
            val, err = parms_vals[iiter], math.sqrt(cov_matrix[iiter][iiter]) if cov_status == 0 and cov_matrix[iiter][iiter] >= 0 else -2
            outDict["p%d"%k]["p%d"%l], outDict["p%d"%k]["p%d_err"%l], outDict["p%d"%k]["p%d_err"%l] = val, err, fitCfg["p%d"%k]["p%d_bnds"%l] if "p%d_bnds"%l in fitCfg["p%d"%k] else (None, None)
            iiter+=1
                        
    # diagonalize covariance matrix
    if cov_status == 0:
        outDict['nStatVars'] = len(parms_vals)
        variations = diagonalize(cov_matrix, parms_vals, sigma=1)
        for iPert, parms_pert in enumerate(variations):
            outDict['stat%d_p' % iPert] = {}
            iiter = 0
            for k in range(0, fitCfg['nParams']):
                outDict['stat%d_p' % iPert]["p%d"%k] = {}
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    outDict['stat%d_p' % iPert]["p%d"%k]["p%d"%l] = parms_pert[iiter]
                    iiter+=1
    
        variations = diagonalize(cov_matrix, parms_vals, sigma=-1)
        for iPert, parms_pert in enumerate(variations):
            outDict['stat%d_m' % iPert] = {}
            iiter = 0
            for k in range(0, fitCfg['nParams']):
                outDict['stat%d_m' % iPert]["p%d"%k] = {}
                for l in range(0, fitCfg["p%d"%k]['nParams']):
                    outDict['stat%d_m' % iPert]["p%d"%k]["p%d"%l] = parms_pert[iiter]
                    iiter+=1

    return outDict



def parameterizeMean(bhist, comp, outDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=15, yRatio=1.15, recoilLow=-100, recoilHigh=100):
    
    xMin, xMax = min(binning_qT), max(binning_qT)
    
    fitFuncName = fitCfg['func_name']
    iParams = fitCfg['parms_init']
    cParams = fitCfg['parms_cte']
    
    axis_name = "recoil_perp" if comp == "perp" else "recoil_para_qT" # 
    
    fitF = getattr(rf, fitFuncName + "_")()
    fit = ROOT.TF1("fit", fitF, xMin, xMax)
    for i, iPar in enumerate(iParams):
        if cParams[i]: fit.FixParameter(i, iPar)
        else: fit.SetParameter(i, iPar)


    g = ROOT.TGraphErrors()
    g.SetLineColor(ROOT.kBlack)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.55)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    
    g_yields = ROOT.TGraphErrors()
    
    means, means_err = [], []
    yieldsDict = {}
    for iBin in range(1, len(binning_qT)):
    
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        qTlow, qThigh = binning_qT[iBin-1], binning_qT[iBin]
        #print("############# Recoil bin %d [%.2f, %.2f]" % (iBin, qTlow, qThigh))

        h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h =  h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum
        h =  h[{axis_name: s[complex(0,recoilLow):complex(0,recoilHigh)]}]
        h_root = narf.hist_to_root(h) 
        
        g.SetPoint(iBin-1, qT, h_root.GetMean())
        g.SetPointError(iBin-1, 0, h_root.GetMeanError())
        
        means.append(h_root.GetMean())
        means_err.append(h_root.GetMeanError())
      
        yield_, yield_err = h.sum().value, math.sqrt(h.sum().variance)
        yieldsDict[iBin] = {}
        yieldsDict[iBin]['yield'] = yield_
        yieldsDict[iBin]['yield_err'] = yield_err
        
        g_yields.SetPoint(iBin-1, qT, yield_)
        g_yields.SetPointError(iBin-1, 0, yield_err)
           


    result = g.Fit(fit.GetName(), "NSE", "", xMin, xMax) 
    chisq = fit.GetChisquare()/fit.GetNDF()
    fit.SetLineColor(ROOT.kRed)
    #fit.GetXaxis().SetRangeUser(0, 200)
    fit.SetLineWidth(2)

        
    cov = result.GetCorrelationMatrix()       
    values = result.GetConfidenceIntervals(0.68, False)
    uncBand = ROOT.TGraphErrors()
    uncBand_ratio = ROOT.TGraphErrors()
    for i in range(len(values)):
        uncBand.SetPoint(i, g.GetX()[i], fit.Eval(g.GetX()[i]))
        uncBand.SetPointError(i, 0, values[i])
        uncBand_ratio.SetPoint(i, g.GetX()[i], 1)
        uncBand_ratio.SetPointError(i, 0, values[i])

    # ratio 
    g_ratio = ROOT.TGraphErrors()
    g_ratio.SetLineColor(ROOT.kBlack)
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(0.55)
    g_ratio.SetMarkerColor(ROOT.kBlack)
    g_ratio.SetLineColor(ROOT.kBlack)
            
    iPoint = 0
    for iBin in range(1, len(binning_qT)):
    
        qT = 0.5*(binning_qT[iBin-1] + binning_qT[iBin])
        y = means[iBin-1]
        y_err = means_err[iBin-1]
        y_fit = fit.Eval(qT)

        if y_err > 0: 
            ratio = y / y_fit
            ratio_err = y_err/abs(y)
        else:
            ratio = 1.0
            ratio_err = 0.0
        g_ratio.SetPoint(iPoint, qT, ratio)
        g_ratio.SetPointError(iPoint, 0, ratio_err)
        iPoint += 1  

        print(qT, means[iBin-1], y_fit, means[iBin-1]/y_fit)
    
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "<U_{#parallel}> (GeV)" if comp=="para" else "<U_{#perp}  > (GeV)",
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
            
        'ratiofraction'     : 0.25,
        'ytitleR'           : "Ratio",
        'yminR'             : (1-(yRatio-1)),
        'ymaxR'             : yRatio,
    }         

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")

    uncBand.SetFillColor(ROOT.kGray)
    uncBand.Draw("3SAME")
    g.Draw("PE SAME")
    fit.Draw("L SAME")  
    padT.RedrawAxis()
    padT.RedrawAxis("G")
    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.2, 0.85, procLabel)
    latex.DrawLatex(0.2, 0.80, metLabel)
    latex.DrawLatex(0.2, 0.75, "#chi^{2}/ndof = %.2f" % chisq)
    
    latex.SetTextSize(0.035)
    for i in range(0, len(iParams)):
        latex.DrawLatex(0.6, 0.86-i*0.04, "a_{%d} = %.2e #pm %.2e" % (i, fit.GetParameter(i), fit.GetParError(i)))
        
       
        
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    padB.SetGrid()
    padB.SetTickx()
    padB.SetTicky()
    dummyB.Draw("HIST")
        
    #if doFit and fit:
    uncBand_ratio.SetFillColor(18)
    uncBand_ratio.Draw("3SAME")
    dummyL.Draw("SAME")
    g_ratio.Draw("PE0 SAME")
    padB.RedrawAxis()
    padB.RedrawAxis("G")
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/parametric_mean.png" % outDir)
    canvas.SaveAs("%s/parametric_mean.pdf" % outDir)
    canvas.Delete()
    
    plotYields(g_yields, "%s/yields" % outDir, procLabel, metLabel, xMin=min(binning_qT), xMax=max(binning_qT))
    
    outDict = {}
    outDict["procNames"] = procNames
    outDict["func"] = fitF
    outDict["funcName"] = fitFuncName
    outDict["nParams"] = len(iParams)
    for i, iPar in enumerate(iParams):
        outDict["p%d" % i] = fit.GetParameter(i)
        outDict["p%d_err" % i] = fit.GetParError(i)
        outDict["p%d_isCte" % i] = cParams[i]
     
     
    with open("%s/parametric_mean.json" % outDir, "w") as outfile: json.dump(outDict, outfile, indent=4)
    with open("%s/yields.json" % outDir, "w") as outfile: json.dump(yieldsDict, outfile, indent=4)
    

def export(exportCfg, fOut, jsInF, jsInF_mean=None):

    jsIn = functions.loadJSON(jsInF)
    if jsInF_mean != None: jsIn_mean = functions.loadJSON(jsInF_mean)
    
    nGauss = len(exportCfg['mean'])
    cfgOut = {}
    cfgOut['nGauss'] = nGauss
    for i in range(0, nGauss):
    
        # mean
        if isinstance(exportCfg['mean'][i], str): cfg = jsIn[exportCfg['mean'][i]]
        else:
            cfg = {}
            cfg['func'] = "[0]"
            cfg['nParams'] = 1
            cfg['p0'] = exportCfg['mean'][i]
            cfg['p0_err'] = -1
            cfg['p0_isCte'] = True
        if jsInF_mean != None: # add parameterized mean to mean
            func_mean = jsIn_mean['func']
            func_nParams = jsIn_mean['nParams']
            p_base = cfg['nParams']
            for iPar in reversed(range(0, func_nParams)):
                p = p_base + iPar
                func_mean = func_mean.replace("[%d]"%iPar, "[%d]"%p) # shift the parameters
                cfg['p%d'%p] = jsIn_mean['p%d'%iPar]
                cfg['p%d_err'%p] = jsIn_mean['p%d_err'%iPar]
                cfg['p%d_isCte'%p] = jsIn_mean['p%d_isCte'%iPar]
                cfg['nParams'] += 1
            cfg['func'] = "(%s) + (%s)" % (cfg['func'], func_mean)
        cfgOut['mean%d'%(i+1)] = cfg
        
        # sigma
        if isinstance(exportCfg['sigma'][i], str): cfg = jsIn[exportCfg['sigma'][i]]
        else:
            cfg = {}
            cfg['func'] = "[0]"
            cfg['nParams'] = 1
            cfg['p0'] = exportCfg['sigma'][i]
            cfg['p0_err'] = -1
            cfg['p0_isCte'] = True
        cfgOut['sigma%d'%(i+1)] = cfg
        
        # norm
        if i == nGauss-1: continue
        if isinstance(exportCfg['norm'][i], str): cfg = jsIn[exportCfg['norm'][i]]
        else:
            cfg = {}
            cfg['func'] = "[0]"
            cfg['nParams'] = 1
            cfg['p0'] = exportCfg['norm'][i]
            cfg['p0_err'] = -1
            cfg['p0_isCte'] = True
        cfgOut['norm%d'%(i+1)] = cfg
   

    # statistical variations
    if 'nStatVars' in jsIn:
        nStatVars = jsIn['nStatVars']
        cfgOut['nStatVars'] = nStatVars
        for nStat in range(0, nStatVars):
            for statVar in ["p", "m"]:
                statLabel = "stat%d_%s" % (nStat, statVar)
                cfgOut[statLabel] = {}
                for i in range(0, nGauss):
                    # mean
                    if isinstance(exportCfg['mean'][i], str): cfg = jsIn[statLabel][exportCfg['mean'][i]]
                    else:
                        cfg = {}
                        cfg['p0'] = exportCfg['mean'][i]
                    if jsInF_mean != None: # add parameterized mean to mean
                        func_mean = jsIn_mean['func']
                        func_nParams = jsIn_mean['nParams']
                        p_base = cfg['nParams']
                        for iPar in reversed(range(0, func_nParams)):
                            p = p_base + iPar
                            func_mean = func_mean.replace("[%d]"%iPar, "[%d]"%p) # shift the parameters
                            cfg['p%d'%p] = jsIn_mean['p%d'%iPar]
                            cfg['p%d_err'%p] = jsIn_mean['p%d_err'%iPar]
                            cfg['p%d_isCte'%p] = jsIn_mean['p%d_isCte'%iPar]
                            cfg['nParams'] += 1
                        cfg['func'] = "(%s) + (%s)" % (cfg['func'], func_mean)
                    cfgOut[statLabel]['mean%d'%(i+1)] = cfg
                    
                    # sigma
                    if isinstance(exportCfg['sigma'][i], str): cfg = jsIn[statLabel][exportCfg['sigma'][i]]
                    else:
                        cfg = {}
                        cfg['p0'] = exportCfg['sigma'][i]
      
                    cfgOut[statLabel]['sigma%d'%(i+1)] = cfg
                    
                    # norm
                    if i == nGauss-1: continue
                    if isinstance(exportCfg['norm'][i], str): cfg = jsIn[statLabel][exportCfg['norm'][i]]
                    else:
                        cfg = {}
                        cfg['p0'] = exportCfg['norm'][i]
                    cfgOut[statLabel]['norm%d'%(i+1)] = cfg   
            
    functions.writeJSON(fOut, cfgOut)
    