
import sys,array,math,os,copy,decimal
import numpy as np
import ctypes
import json
import time

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

#sys.path.insert(0, "/home/j/jaeyserm/analysis/WRemnants/env/lib/python3.10/site-packages")
import zfit
import zfit.z.numpy  # numpy-like backend
from zfit import z  # backend

import fitter
import hist
import boost_histogram as bh



mirror = False

__topRight__ = "199 pb^{#minus1} (13 TeV)"

__min_recoil__, __max_recoil__ = -500, 500
__min_mean__, __max_mean__ = -50, 500
__min_sigma__, __max_sigma__ = 0.1, 100

def getParamIdx(name, params):
    
    for i in range(0, len(params)): 
       
        if name == params[i].GetName(): return i
            
    return -1
    


def diagonalize(fitRes, verbose=False):

    cov = fitRes.covarianceMatrix() # covariance matrix: diagonal elements are the variances of each of the parameters
    if verbose: cov.Print()
    params = fitRes.floatParsFinal()
    npars = params.getSize()
    sigma = 1
    
    nom = np.zeros(npars)
    for i in range(0, npars): nom[i] = params[i].getVal()

        
    cov_ = np.zeros((npars, npars))
    for i in range(0, npars):
        for j in range(0, npars):
            cov_[i][j] = cov[i][j]
         
        
    # eigenvalues and eigenvectors
    eig_vals, eig_vec = np.linalg.eig(cov_)
    eig_vec_inv = np.linalg.inv(eig_vec)
        
    #print("COVARIANCE")
    #print(cov_)
    #print("EIGENVECTORS")
    #print(eig_vec)
    #print("EIGENVECTORS INV")
    #print(eig_vec_inv)
    #print("EIGENVALUES")
    #print(eig_vals)
    #print("NOMINAL")
    #print(nom)
                
    ret = []
    for iVar in range(0, npars):
            
        #print("******** vary", iVar+1, "********")
            
        dnom = copy.deepcopy(eig_vals)
        for i, eig in enumerate(dnom): 
           
            if i == iVar: dnom[i] = sigma*math.sqrt(eig)
            else: dnom[i] = 0
               
        #print("PERTURBED NON-ROTATED")
        #print(dnom)
            
        # rotate perturbation back to nominal base
        dnom = np.dot(eig_vec, dnom)
        nom_pert = nom + dnom
        
        ret.append(nom_pert)
        
        #print("%f %f %f" % (nom_pert[0], nom_pert[1], nom_pert[2]))
        #print("PERTURBED ROTATED")
        #print(nom_pert)
     
    return ret

def mirrorx(h):

    hnew = h.Clone("%s_mirrorx" % h.GetName())
    for i in range(1, hnew.GetNbinsX()+1):
        idx = hnew.GetNbinsX()-i+1
        hnew.SetBinContent(i, h.GetBinContent(idx))
        hnew.SetBinError(i, h.GetBinError(idx))

    #sys.exit()
    return hnew


 
def doFitMultiGauss_scipy_bhist(bhist, comp, fitCfg, label, outDir_, recoil_qTbins, ratio=True, funcJs=None, excludeBins=[], qTstartBin = 1, qTmin=0, qTmax=300, doFit=True, rebin=1, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, propagate=True, bkgCfg={}, orderSigmas=False, yRatio = 1.15):
    
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
    
    recoil_bins = bhist.axes[1].edges
    centers = bhist.axes[1].centers
    qT_bins = bhist.axes[0].edges
    

    
    functions.prepareDir(outDir_, True) # remove and recreate dir
    nGauss = len(fitCfg['mean'])
    dof_param = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(3)
    dof = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(0) + dof_param # if cfg=0 or 3, add DOF

    func_model = rf.func_4gauss
    parms = []
       
    # construct gauss
    mean_ptrs, sigma_ptrs, norm_ptrs, gauss_ptrs = [], [], [], []
    for iGauss in range(1, nGauss+1):
        parms.append(fitCfg['sigma'][iGauss-1])

    import hist
    s = hist.tag.Slicer()
                
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()        
    outDict = {}
    outDict['nGauss'] = nGauss
    recoil_qTbins_ = [0]
 
    if qTstartBin >= len(recoil_qTbins) or qTstartBin < 1: sys.exit("qTstartBin not in range (should be at least 1)")
    if recoil_qTbins[qTstartBin-1] < qTmin: sys.exit("qTstartBin is lower than the minimum range defined")
    if recoil_qTbins[qTstartBin] > qTmax: sys.exit("qTstartBin is higher than the maximum range defined")
    qTbin = qTstartBin
    direction = 1 if qTstartBin == 1 else -1

    while True:
        
        if qTbin >= len(recoil_qTbins): break
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        qTlow, qThigh = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        if qTlow >= qTmax: break
        
        print("##### DO BIN", qTbin, qT, qTlow, "->",  qThigh)

        
        h =  bhist[{"qTbinned": s[complex(0,qTlow):complex(0,qThigh)]}]
        h =  h[{"qTbinned": s[0:hist.overflow:hist.sum]}] # remove the overflow in the sum
        
        
        print(h)
        
       
        yield_ = h.sum().value
        yield_err = math.sqrt(h.sum().variance)
        res = fitter.fit_hist(h, func_model, np.array(parms), max_iter = 5, edmtol = 1e-5, mode = "nll")

        parms = res['x']
        pdf_values = [yield_*func_model(np.array([rec]), np.array(parms)).numpy() for rec in centers]
        
        print(pdf_values)
        print(res)
        
        # get the PDF values
        
        
            
        fitValid = True
        
        
        if direction > 0 and qTbin == qTstartBin and not qTstartBin == 1:  # avoid plotting/processing twice the start bin
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
        
        hist_root = narf.hist_to_root(h)
        hist_root.Scale(1., "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")

        
        #pdf_values = zfit.run(model.pdf(centers, norm_range=recoil_binned))
        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers)): g_pdf.SetPoint(i, centers[i], pdf_values[i])
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
        #if dof > 0 and yield_ > 0:
        #    latex.DrawLatex(0.20, 0.55, "#color[%d]{Cov. matrix quality = %d}" % (2 if fitRes.covQual() != 3 else 1, fitRes.covQual()))
        #    latex.DrawLatex(0.20, 0.50, "#color[%d]{MINUIT status = %d}" % (2 if fitRes.status() != 0 else 1, fitRes.status()))
 
        latex.SetTextSize(0.035)
        
        # save output
        outDict[qTbin] = {}
        outDict[qTbin]['yield'] = yield_
        outDict[qTbin]['yield_err'] = yield_err
        
        for i, s in enumerate(mean_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04, "#mu_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            #if fitCfg['mean_cfg'][i] >= 0 else mean_ptrs[abs(fitCfg['mean_cfg'][i])-1].value())) 
            outDict[qTbin]['mean%d' % (i+1)], outDict[qTbin]['mean%d_err' % (i+1)] = val, err 
        for i, s in enumerate(sigma_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs), "#sigma_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['sigma%d' % (i+1)], outDict[qTbin]['sigma%d_err' % (i+1)] = val, err
        for i, s in enumerate(norm_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs)*2, "n_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['norm%d' % (i+1)], outDict[qTbin]['norm%d_err' % (i+1)] = val, err
        

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
        

        recoil_qTbins_.append(recoil_qTbins[qTbin])
        if qTlow == qTmin:
            qTbin = qTstartBin # redo the fit at qTstartBin for initial conditions
            if qTstartBin != 1: qTbin -= 1 # subtract 1, as 1 is added below with direction
            direction = +1
        qTbin += direction # update bin 
        sys.exit()
    
    outDict['qTbins'] = recoil_qTbins_
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)

    with open("%s/results.json" % outDir_) as f: outDict = json.load(f)
    doGlobalPlot(bhist, comp, outDict, recoil_qTbins_, outDir_, label, funcJs=funcJs, ratio=ratio, rebin=rebin, bkgCfg=bkgCfg, qTmin=qTmin, qTmax=qTmax, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, yRatio=yRatio)

    
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
    
    sys.exit()
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")
    hIn_y = hIn.ProjectionY("hIn_y")
    
    recoil_bins = [hIn_y.GetBinLowEdge(i) for i in range(1, hIn_y.GetNbinsX()+2)]
    centers = [hIn_y.GetBinCenter(i) for i in range(1, hIn_y.GetNbinsX()+1)]
    bins = zfit.binned.VariableBinning(recoil_bins, name="recoil_bins")
    recoil = zfit.Space('recoil_bins', limits=(min(recoil_bins), max(recoil_bins)))
    recoil_binned = zfit.Space('recoil_bins', binning=bins)
    
 
    functions.prepareDir(outDir_, True) # remove and recreate dir
    nGauss = len(fitCfg['mean'])
    dof_param = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(3)
    dof = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(0) + dof_param # if cfg=0 or 3, add DOF

    pdf = getattr(rf, "model_4gauss_")
    parms = []
    
    


    # construct gauss
    mean_ptrs, sigma_ptrs, norm_ptrs, gauss_ptrs = [], [], [], []
    for iGauss in range(1, nGauss+1):
        if fitCfg['mean_cfg'][iGauss-1] >= 0:
            mean = zfit.Parameter("mean%d"%iGauss, fitCfg['mean'][iGauss-1], __min_mean__, __max_mean__, step_size=0.001)
            if fitCfg['mean_cfg'][iGauss-1] in [1, 2]: mean.floating = False
        else:
            dependent_tensor = mean_ptrs[abs(fitCfg['mean_cfg'][iGauss-1]-1)]
            mean = zfit.ComposedParameter("mean%d"%iGauss, dependent_tensor)
        sigma = zfit.Parameter("sigma%d"%iGauss, fitCfg['sigma'][iGauss-1], __min_sigma__, __max_sigma__, step_size=0.001)
        if fitCfg['sigma_cfg'][iGauss-1] in [1, 2]: sigma.floating = False
        gauss = zfit.pdf.Gauss(obs=recoil, mu=mean, sigma=sigma, name="gauss%d"%iGauss)
        gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss, recoil_binned)
        mean_ptrs.append(mean)
        sigma_ptrs.append(sigma)
        gauss_ptrs.append(gauss_binned)
    
    # identical norms for each gauss
    for iGauss in range(1, nGauss):
        norm = zfit.Parameter("norm%d"%iGauss, fitCfg['norm'][iGauss-1], 0, 1, step_size=0.001)
        if fitCfg['norm_cfg'][iGauss-1] in [1, 2]: norm.floating = False
        norm_ptrs.append(norm)

    
    # make the TF1s 
    if dof_param > 0 and funcJs == None: sys.exit("Parametric fit expected but input json not defined!")
    funcs = {}
    if funcJs != None:
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcJs:
                funcs["mean%d" % iGauss] = ROOT.TF1("mean%d" % iGauss, funcJs['mean%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['mean%d' % iGauss]['nParams']): funcs["mean%d" % iGauss].SetParameter(iParam, funcJs['mean%d' % iGauss]['p%d' % iParam])
            if "sigma%d" % iGauss in funcJs:
                funcs["sigma%d" % iGauss] = ROOT.TF1("sigma%d" % iGauss, funcJs['sigma%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['sigma%d' % iGauss]['nParams']): funcs["sigma%d" % iGauss].SetParameter(iParam, funcJs['sigma%d' % iGauss]['p%d' % iParam])

        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcJs:
                funcs["norm%d" % iGauss] = ROOT.TF1("norm%d" % iGauss, funcJs['norm%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['norm%d' % iGauss]['nParams']): funcs["norm%d" % iGauss].SetParameter(iParam, funcJs['norm%d' % iGauss]['p%d' % iParam])
       
       
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()
    model = zfit.pdf.BinnedSumPDF(gauss_ptrs, fracs=norm_ptrs, obs=recoil_binned, name="gauss")

    
    zfit.settings.set_seed(1)
    minimizer = zfit.minimize.Minuit() # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1() # tol=1e-4, verbosity=5, gradient=True
    
    
    h = narf.root_to_hist(hIn_y, axis_names=["recoil_bins"])
    data_ = zfit.data.BinnedData.from_hist(h)
    data_vals = data_.values()
    data_vars = data_.variances()

    
    # sampler trick: do initialization, replace data later
    somerandomdata = np.array(model.sample(1000).values())
    sampler = model.create_sampler(1000)
    sampler.holder.variances = data_vars
    sampler.sample_holder.assign(somerandomdata)
    

    start = time.time()
    print("##### START BinnedNLL")
    nll = zfit.loss.BinnedNLL(model=model, data=sampler)
    #nll = zfit.loss.BinnedNLL(model=model, data=data_)
    end = time.time()
    print("##### END BinnedNLL, DURATION", (end-start))
    
    start = time.time()
    print("##### START minimize")
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### END minimize, DURATION", (end-start))
         
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
        print(qTbin)
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

        
        yield_err = ctypes.c_double(1.)
        yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        yield_err = yield_err.value

        h = narf.root_to_hist(hist, axis_names=["recoil_bins"]) 
        data__ = zfit.data.BinnedData.from_hist(h)
        data_vals = data__.values()
        data_vars = data__.variances()


        # set the variable according to the functional forms
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcs:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].set_value(funcs["mean%d" % iGauss].Eval(qT))
            elif not propagate:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].set_value(fitCfg['mean'][iGauss-1])
            if "sigma%d" % iGauss in funcs: 
                sigma_ptrs[iGauss-1].set_value(funcs["sigma%d" % iGauss].Eval(qT))
            elif not propagate:    
                sigma_ptrs[iGauss-1].set_value(fitCfg['sigma'][iGauss-1])
            
        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcs:
                norm_ptrs[iGauss-1].set_value(funcs["norm%d" % iGauss].Eval(qT))
            elif not propagate:
                norm_ptrs[iGauss-1].set_value(fitCfg['norm'][iGauss-1])
        
        
        
        errCalc = "hesse" # minos, hesse

        fitValid = True
        if dof > 0 and yield_ > 0:
            
            start = time.time()
            print("##### START BinnedNLL")
            sampler.sample_holder.assign(data_vals)
            sampler.holder.variances = data_vars
            
            #data_.holder.values = data_vals
            #data_.holder.variances = data_vars
            
            # create new loss per point
            #nll = zfit.loss.BinnedNLL(model=model, data=data__) #BinnedChi2 BinnedNLL # , options={"subtr_const" : False}
            #nll_ = nll.create_new()
            
            # use initial nll, but update the data
            #if not nll:
            #    nll = zfit.loss.BinnedNLL(model=model, data=data__)
            #else: 
            #    nll = nll.create_new(data=data__)
            end = time.time()
            print("##### END BinnedNLL, DURATION", (end-start))

            
            
            start = time.time()
            print("##### START MIMIMIZATION")
            result = minimizer.minimize(loss=nll)
            end = time.time()
            print("##### END MIMIMIZATION, DURATION", (end-start))
            params = result.params
            print(result)
            print(params)
            #fitValid = (fitRes.covQual() == 3 and fitRes.status() == 0)
            
            if errCalc == "hesse":
                start = time.time()
                print("##### START HESSE")
                param_errors = result.hesse(name='hesse')
                end = time.time()
                print("##### END HESSE, DURATION", (end-start))
            
            if errCalc == "minos":
                start = time.time()
                print("##### START ERRORS")
                param_errors, _ = result.errors(name='hesse')
                end = time.time()
                print("##### END HESSE, DURATION", (end-start))
                for p in param_errors:
                    print(param_errors[p]['upper'], param_errors[p]['lower'])
                    param_errors[p]['error'] = 0.5*(abs(param_errors[p]['lower']) + abs(param_errors[p]['upper']))
            
        else: fitValid = True
        
      
        
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
        
        hist_root = narf.hist_to_root(h)
        #hist_root.Scale(1./yield_, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")

        pdf_values = zfit.run(model.pdf(centers, norm_range=recoil_binned))
        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers)): g_pdf.SetPoint(i, centers[i], pdf_values[i]*yield_)
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
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin_, qTmax_))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        #if dof > 0 and yield_ > 0:
        #    latex.DrawLatex(0.20, 0.55, "#color[%d]{Cov. matrix quality = %d}" % (2 if fitRes.covQual() != 3 else 1, fitRes.covQual()))
        #    latex.DrawLatex(0.20, 0.50, "#color[%d]{MINUIT status = %d}" % (2 if fitRes.status() != 0 else 1, fitRes.status()))
 
        latex.SetTextSize(0.035)
        
        # save output
        outDict[qTbin] = {}
        outDict[qTbin]['yield'] = yield_
        outDict[qTbin]['yield_err'] = yield_err
        
        for i, s in enumerate(mean_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04, "#mu_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            #if fitCfg['mean_cfg'][i] >= 0 else mean_ptrs[abs(fitCfg['mean_cfg'][i])-1].value())) 
            outDict[qTbin]['mean%d' % (i+1)], outDict[qTbin]['mean%d_err' % (i+1)] = val, err 
        for i, s in enumerate(sigma_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs), "#sigma_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['sigma%d' % (i+1)], outDict[qTbin]['sigma%d_err' % (i+1)] = val, err
        for i, s in enumerate(norm_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs)*2, "n_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['norm%d' % (i+1)], outDict[qTbin]['norm%d_err' % (i+1)] = val, err
        

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
        continue
        
    
    outDict['qTbins'] = recoil_qTbins_
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)

    with open("%s/results.json" % outDir_) as f: out = json.load(f)
    doGlobalPlot(fInName, proc, comp, out, recoil_qTbins, outDir_, label, funcJs=funcJs, ratio=ratio, bkgCfg=bkgCfg, qTmin=qTmin, qTmax=qTmax, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, yRatio=yRatio)
 

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
 


def doFitMultiGauss_zfit(fInName, proc, comp, fitCfg, label, outDir_, recoil_qTbins, ratio=True, funcJs=None, excludeBins=[], qTstartBin = 1, qTmin=0, qTmax=300, doFit=True, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, propagate=True, bkgCfg={}, orderSigmas=False, yRatio = 1.15):
    
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
    bins = zfit.binned.VariableBinning(recoil_bins, name="recoil_bins")
    recoil = zfit.Space('recoil_bins', limits=(min(recoil_bins), max(recoil_bins)))
    recoil_binned = zfit.Space('recoil_bins', binning=bins)
    
 
    functions.prepareDir(outDir_, True) # remove and recreate dir
    nGauss = len(fitCfg['mean'])
    dof_param = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(3)
    dof = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(0) + dof_param # if cfg=0 or 3, add DOF





    # construct gauss
    mean_ptrs, sigma_ptrs, norm_ptrs, gauss_ptrs = [], [], [], []
    for iGauss in range(1, nGauss+1):
        if fitCfg['mean_cfg'][iGauss-1] >= 0:
            mean = zfit.Parameter("mean%d"%iGauss, fitCfg['mean'][iGauss-1], __min_mean__, __max_mean__, step_size=0.001)
            if fitCfg['mean_cfg'][iGauss-1] in [1, 2]: mean.floating = False
        else:
            dependent_tensor = mean_ptrs[abs(fitCfg['mean_cfg'][iGauss-1]-1)]
            mean = zfit.ComposedParameter("mean%d"%iGauss, dependent_tensor)
        sigma = zfit.Parameter("sigma%d"%iGauss, fitCfg['sigma'][iGauss-1], __min_sigma__, __max_sigma__, step_size=0.001)
        if fitCfg['sigma_cfg'][iGauss-1] in [1, 2]: sigma.floating = False
        gauss = zfit.pdf.Gauss(obs=recoil, mu=mean, sigma=sigma, name="gauss%d"%iGauss)
        gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss, recoil_binned)
        mean_ptrs.append(mean)
        sigma_ptrs.append(sigma)
        gauss_ptrs.append(gauss_binned)
    
    # identical norms for each gauss
    for iGauss in range(1, nGauss):
        norm = zfit.Parameter("norm%d"%iGauss, fitCfg['norm'][iGauss-1], 0, 1, step_size=0.001)
        if fitCfg['norm_cfg'][iGauss-1] in [1, 2]: norm.floating = False
        norm_ptrs.append(norm)

    
    # make the TF1s 
    if dof_param > 0 and funcJs == None: sys.exit("Parametric fit expected but input json not defined!")
    funcs = {}
    if funcJs != None:
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcJs:
                funcs["mean%d" % iGauss] = ROOT.TF1("mean%d" % iGauss, funcJs['mean%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['mean%d' % iGauss]['nParams']): funcs["mean%d" % iGauss].SetParameter(iParam, funcJs['mean%d' % iGauss]['p%d' % iParam])
            if "sigma%d" % iGauss in funcJs:
                funcs["sigma%d" % iGauss] = ROOT.TF1("sigma%d" % iGauss, funcJs['sigma%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['sigma%d' % iGauss]['nParams']): funcs["sigma%d" % iGauss].SetParameter(iParam, funcJs['sigma%d' % iGauss]['p%d' % iParam])

        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcJs:
                funcs["norm%d" % iGauss] = ROOT.TF1("norm%d" % iGauss, funcJs['norm%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['norm%d' % iGauss]['nParams']): funcs["norm%d" % iGauss].SetParameter(iParam, funcJs['norm%d' % iGauss]['p%d' % iParam])
       
       
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()
    model = zfit.pdf.BinnedSumPDF(gauss_ptrs, fracs=norm_ptrs, obs=recoil_binned, name="gauss")

    
    zfit.settings.set_seed(1)
    minimizer = zfit.minimize.Minuit() # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1() # tol=1e-4, verbosity=5, gradient=True
    
    
    h = narf.root_to_hist(hIn_y, axis_names=["recoil_bins"])
    data_ = zfit.data.BinnedData.from_hist(h)
    data_vals = data_.values()
    data_vars = data_.variances()

    
    # sampler trick: do initialization, replace data later
    somerandomdata = np.array(model.sample(1000).values())
    sampler = model.create_sampler(1000)
    sampler.holder.variances = data_vars
    sampler.sample_holder.assign(somerandomdata)
    

    start = time.time()
    print("##### START BinnedNLL")
    nll = zfit.loss.BinnedNLL(model=model, data=sampler)
    #nll = zfit.loss.BinnedNLL(model=model, data=data_)
    end = time.time()
    print("##### END BinnedNLL, DURATION", (end-start))
    
    start = time.time()
    print("##### START minimize")
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### END minimize, DURATION", (end-start))
         
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
        print(qTbin)
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

        
        yield_err = ctypes.c_double(1.)
        yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        yield_err = yield_err.value

        h = narf.root_to_hist(hist, axis_names=["recoil_bins"]) 
        data__ = zfit.data.BinnedData.from_hist(h)
        data_vals = data__.values()
        data_vars = data__.variances()


        # set the variable according to the functional forms
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcs:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].set_value(funcs["mean%d" % iGauss].Eval(qT))
            elif not propagate:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].set_value(fitCfg['mean'][iGauss-1])
            if "sigma%d" % iGauss in funcs: 
                sigma_ptrs[iGauss-1].set_value(funcs["sigma%d" % iGauss].Eval(qT))
            elif not propagate:    
                sigma_ptrs[iGauss-1].set_value(fitCfg['sigma'][iGauss-1])
            
        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcs:
                norm_ptrs[iGauss-1].set_value(funcs["norm%d" % iGauss].Eval(qT))
            elif not propagate:
                norm_ptrs[iGauss-1].set_value(fitCfg['norm'][iGauss-1])
        
        
        
        errCalc = "hesse" # minos, hesse

        fitValid = True
        if dof > 0 and yield_ > 0:
            
            start = time.time()
            print("##### START BinnedNLL")
            sampler.sample_holder.assign(data_vals)
            sampler.holder.variances = data_vars
            
            #data_.holder.values = data_vals
            #data_.holder.variances = data_vars
            
            # create new loss per point
            #nll = zfit.loss.BinnedNLL(model=model, data=data__) #BinnedChi2 BinnedNLL # , options={"subtr_const" : False}
            #nll_ = nll.create_new()
            
            # use initial nll, but update the data
            #if not nll:
            #    nll = zfit.loss.BinnedNLL(model=model, data=data__)
            #else: 
            #    nll = nll.create_new(data=data__)
            end = time.time()
            print("##### END BinnedNLL, DURATION", (end-start))

            
            
            start = time.time()
            print("##### START MIMIMIZATION")
            result = minimizer.minimize(loss=nll)
            end = time.time()
            print("##### END MIMIMIZATION, DURATION", (end-start))
            params = result.params
            print(result)
            print(params)
            #fitValid = (fitRes.covQual() == 3 and fitRes.status() == 0)
            
            if errCalc == "hesse":
                start = time.time()
                print("##### START HESSE")
                param_errors = result.hesse(name='hesse')
                end = time.time()
                print("##### END HESSE, DURATION", (end-start))
            
            if errCalc == "minos":
                start = time.time()
                print("##### START ERRORS")
                param_errors, _ = result.errors(name='hesse')
                end = time.time()
                print("##### END HESSE, DURATION", (end-start))
                for p in param_errors:
                    print(param_errors[p]['upper'], param_errors[p]['lower'])
                    param_errors[p]['error'] = 0.5*(abs(param_errors[p]['lower']) + abs(param_errors[p]['upper']))
            
        else: fitValid = True
        
      
        
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
        
        hist_root = narf.hist_to_root(h)
        #hist_root.Scale(1./yield_, "width")
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        hist_root.Draw("PE SAME")

        pdf_values = zfit.run(model.pdf(centers, norm_range=recoil_binned))
        g_pdf = ROOT.TGraphErrors()
        for i in range(0, len(centers)): g_pdf.SetPoint(i, centers[i], pdf_values[i]*yield_)
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
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin_, qTmax_))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist_root.GetMean(), hist_root.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist_root.GetRMS(), hist_root.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        #if dof > 0 and yield_ > 0:
        #    latex.DrawLatex(0.20, 0.55, "#color[%d]{Cov. matrix quality = %d}" % (2 if fitRes.covQual() != 3 else 1, fitRes.covQual()))
        #    latex.DrawLatex(0.20, 0.50, "#color[%d]{MINUIT status = %d}" % (2 if fitRes.status() != 0 else 1, fitRes.status()))
 
        latex.SetTextSize(0.035)
        
        # save output
        outDict[qTbin] = {}
        outDict[qTbin]['yield'] = yield_
        outDict[qTbin]['yield_err'] = yield_err
        
        for i, s in enumerate(mean_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04, "#mu_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            #if fitCfg['mean_cfg'][i] >= 0 else mean_ptrs[abs(fitCfg['mean_cfg'][i])-1].value())) 
            outDict[qTbin]['mean%d' % (i+1)], outDict[qTbin]['mean%d_err' % (i+1)] = val, err 
        for i, s in enumerate(sigma_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs), "#sigma_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['sigma%d' % (i+1)], outDict[qTbin]['sigma%d_err' % (i+1)] = val, err
        for i, s in enumerate(norm_ptrs):
            val = zfit.run(s)
            err = param_errors[s]['error'] if s in param_errors else 0
            latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs)*2, "n_{%d} = %.3f #pm %.3f" % (i+1, val, err))
            outDict[qTbin]['norm%d' % (i+1)], outDict[qTbin]['norm%d_err' % (i+1)] = val, err
        

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
        continue
        
    
    outDict['qTbins'] = recoil_qTbins_
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(outDict, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)

    with open("%s/results.json" % outDir_) as f: out = json.load(f)
    doGlobalPlot(fInName, proc, comp, out, recoil_qTbins, outDir_, label, funcJs=funcJs, ratio=ratio, bkgCfg=bkgCfg, qTmin=qTmin, qTmax=qTmax, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, yRatio=yRatio)
 
def doFitMultiGauss(fInName, proc, comp, fitCfg, label, outDir_, recoil_qTbins, ratio=True, funcJs=None, excludeBins=[], qTstartBin = 1, qTmin=0, qTmax=300, doFit=True, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, propagate=True, bkgCfg={}, orderSigmas=False, yRatio = 1.15):

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
    
    #fIn = ROOT.TFile(fIn_)
    
    functions.prepareDir(outDir_, True) # remove and recreate dir
    recoil = ROOT.RooRealVar("recoil", "", 0, __min_recoil__, __max_recoil__) 
    nGauss = len(fitCfg['mean'])
    dof_param = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(3)
    dof = (fitCfg['mean_cfg'] + fitCfg['sigma_cfg'] + fitCfg['norm_cfg']).count(0) + dof_param # if cfg=0 or 3, add DOF

    # construct gauss
    pdfs, norms = ROOT.RooArgList(), ROOT.RooArgList()
    mean_ptrs, sigma_ptrs, norm_ptrs, gauss_ptrs = [], [], [], []
    for iGauss in range(1, nGauss+1):
        if fitCfg['mean_cfg'][iGauss-1] >= 0:
            mean = ROOT.RooRealVar("mean%d" % iGauss, "", fitCfg['mean'][iGauss-1], __min_mean__, __max_mean__)
            if fitCfg['mean_cfg'][iGauss-1] in [1, 2]: mean.setConstant(ROOT.kTRUE)
        else:
            mean = ROOT.RooFormulaVar("mean%d" % iGauss, "@0", ROOT.RooArgList(mean_ptrs[abs(fitCfg['mean_cfg'][iGauss-1])-1]))  
        sigma = ROOT.RooRealVar("sigma%d" % iGauss, "", fitCfg['sigma'][iGauss-1], __min_sigma__, __max_sigma__)
        if fitCfg['sigma_cfg'][iGauss-1] in [1, 2]: sigma.setConstant(ROOT.kTRUE)     
        gauss = ROOT.RooGaussian("gauss%d" % iGauss, "", recoil, mean, sigma)
        pdfs.addOwned(gauss)
        mean_ptrs.append(mean)
        sigma_ptrs.append(sigma)
        gauss_ptrs.append(gauss)
    
    # identical norms for each gauss
    for iGauss in range(1, nGauss):
        norm = ROOT.RooRealVar("norm%d" % iGauss, "", fitCfg['norm'][iGauss-1], 0, 1)
        if fitCfg['norm_cfg'][iGauss-1] in [1, 2]: norm.setConstant(ROOT.kTRUE)
        norms.addOwned(norm)
        norm_ptrs.append(norm)
        
        
    # make the TF1s 
    if dof_param > 0 and funcJs == None: sys.exit("Parametric fit expected but input json not defined!")
    funcs = {}
    if funcJs != None:
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcJs:
                funcs["mean%d" % iGauss] = ROOT.TF1("mean%d" % iGauss, funcJs['mean%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['mean%d' % iGauss]['nParams']): funcs["mean%d" % iGauss].SetParameter(iParam, funcJs['mean%d' % iGauss]['p%d' % iParam])
            if "sigma%d" % iGauss in funcJs:
                funcs["sigma%d" % iGauss] = ROOT.TF1("sigma%d" % iGauss, funcJs['sigma%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['sigma%d' % iGauss]['nParams']): funcs["sigma%d" % iGauss].SetParameter(iParam, funcJs['sigma%d' % iGauss]['p%d' % iParam])

        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcJs:
                funcs["norm%d" % iGauss] = ROOT.TF1("norm%d" % iGauss, funcJs['norm%d' % iGauss]['func'], 0, 500)
                for iParam in range(0, funcJs['norm%d' % iGauss]['nParams']): funcs["norm%d" % iGauss].SetParameter(iParam, funcJs['norm%d' % iGauss]['p%d' % iParam])
        
    pdf_sig = ROOT.RooAddPdf("pdf_sig", '', pdfs, norms)
    g_chi2 = ROOT.TGraphErrors()
    g_yields = ROOT.TGraphErrors()
    
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")

    out = {}
    out['nGauss'] = nGauss
    recoil_qTbins_ = [0]
 
    if qTstartBin >= len(recoil_qTbins) or qTstartBin < 1: sys.exit("qTstartBin not in range (should be at least 1)")
    if recoil_qTbins[qTstartBin-1] < qTmin: sys.exit("qTstartBin is lower than the minimum range defined")
    if recoil_qTbins[qTstartBin] > qTmax: sys.exit("qTstartBin is higher than the maximum range defined")
    qTbin = qTstartBin
    direction = -1
    ###if qTstartBin == 1: direction = 1
    while True:
        print(qTbin)
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

        
        yield_err = ctypes.c_double(1.)
        yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        
        # set the variable according to the functional forms
        for iGauss in range(1, nGauss+1):
            if "mean%d" % iGauss in funcs:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].setVal(funcs["mean%d" % iGauss].Eval(qT))
            elif not propagate:
                if not (fitCfg['mean_cfg'][iGauss-1] < 0 and iGauss != 1):
                    mean_ptrs[iGauss-1].setVal(fitCfg['mean'][iGauss-1])
            if "sigma%d" % iGauss in funcs: 
                sigma_ptrs[iGauss-1].setVal(funcs["sigma%d" % iGauss].Eval(qT))
            elif not propagate:    
                sigma_ptrs[iGauss-1].setVal(fitCfg['sigma'][iGauss-1])
            
        for iGauss in range(1, nGauss):
            if "norm%d" % iGauss in funcs:
                norm_ptrs[iGauss-1].setVal(funcs["norm%d" % iGauss].Eval(qT))
            elif not propagate:
                norm_ptrs[iGauss-1].setVal(fitCfg['norm'][iGauss-1])
            
        # include backgrounds (if configured)
        pdfs_arglist, norms_arglist = ROOT.RooArgList(), ROOT.RooArgList()
        pdfs_arglist.setName("pdfs_arglist_bin%d" % (qTbin))
        norms_arglist.setName("norms_arglist_bin%d" % (qTbin))
        bkg_ptrs = []
        for bkg in bkgCfg:
            if bkgCfg[bkg]['subtract']:
                hIn_bkg = fIn.Get("%s_%s" % (bkg, comp))
                for b in qTbins_:
                    iBin = hIn_x.FindBin(b)
                    h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
                    hist.Add(h, -1.*bkgCfg[bkg]['norm'])                    
            else:
                with open(bkgCfg[bkg]['filename']) as f: jsIn = json.load(f)
                pdfs_bkgs_arglist, norms_bkgs_arglist = ROOT.RooArgList(), ROOT.RooArgList()
                pdfs_bkgs_arglist.setName("pdfs_bkgs_arglist_bin%d_%s" % (qTbin, bkg))
                norms_bkgs_arglist.setName("norms_bkgs_arglist_bin%d_%s" % (qTbin, bkg))
                for iGauss in range(1, jsIn['nGauss']+1):
                    mean = ROOT.RooRealVar("mean%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn[str(qTbin)]["mean%d" % iGauss], __min_mean__, __max_mean__)
                    sigma = ROOT.RooRealVar("sigma%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn[str(qTbin)]["sigma%d" % iGauss], __min_sigma__, __max_sigma__)
                    mean.setConstant(ROOT.kTRUE)
                    sigma.setConstant(ROOT.kTRUE)
                    gauss = ROOT.RooGaussian("gauss%d_bin%d_%s" % (iGauss, qTbin, bkg), "", recoil, mean, sigma)
                    bkg_ptrs.append(mean)
                    bkg_ptrs.append(sigma)
                    bkg_ptrs.append(gauss)
                    pdfs_bkgs_arglist.addOwned(gauss)
                for iGauss in range(1, jsIn['nGauss']):
                    norm = ROOT.RooRealVar("norm%d_bin%d_%s" % (iGauss, qTbin, bkg), "", jsIn[str(qTbin)]["norm%d" % iGauss], 0, 1)
                    norm.setConstant(ROOT.kTRUE)
                    norms_bkgs_arglist.addOwned(norm)
                    bkg_ptrs.append(norm)
                pdf_bkg = ROOT.RooAddPdf("bin%d_%s" % (qTbin, bkg), '', pdfs_bkgs_arglist, norms_bkgs_arglist)
                pdfs_arglist.addOwned(pdf_bkg)
                norm_bkg = ROOT.RooRealVar("relyield_bin%d_%s" % (qTbin, bkg), "", (jsIn[str(qTbin)]['yield']*bkgCfg[bkg]['norm'])/yield_, 0, 1)
                if not bkgCfg[bkg]['float']: norm_bkg.setConstant(ROOT.kTRUE)
                norms_arglist.addOwned(norm_bkg)
                bkg_ptrs.append(pdfs_bkgs_arglist)
                bkg_ptrs.append(norms_bkgs_arglist)
        pdfs_arglist.addClone(pdf_sig) # append signal as the last one
        pdf_tot = ROOT.RooAddPdf("pdf_tot", '', pdfs_arglist, norms_arglist)
        
        rdh = ROOT.RooDataHist("rdh", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(hist, ROOT.kFALSE))
        if orderSigmas:
            for iGauss in range(2, nGauss+1): sigma_ptrs[iGauss-1].setRange(sigma_ptrs[iGauss-2].getVal(), __max_sigma__)
                
        
        if dof > 0 and yield_ > 0:
            fitRes = pdf_tot.fitTo(rdh, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save(ROOT.kTRUE))
            fitValid = (fitRes.covQual() == 3 and fitRes.status() == 0)
        else: fitValid = True
        
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
        padT.SetTickx()
        padT.SetTicky()
        dummyT.Draw("HIST")
        
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
        plt = recoil.frame()
        if dof > 0 and yield_ > 0:
            if fitValid: pdf_tot.plotOn(plt, ROOT.RooFit.VisualizeError(fitRes, 1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kGreen+1), ROOT.RooFit.Normalization(yield_, ROOT.RooAbsReal.NumEvent))
            else: pdf_tot.plotOn(plt, ROOT.RooFit.VisualizeError(fitRes, 1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Normalization(yield_, ROOT.RooAbsReal.NumEvent))
        rdh.plotOn(plt)
        colors = [ROOT.kBlue, ROOT.kGreen+1, ROOT.kOrange, ROOT.kCyan, ROOT.kViolet, ROOT.kGray, ROOT.kGray, ROOT.kGray, ROOT.kGray]
        for i, g in enumerate(gauss_ptrs):
            n = norm_ptrs[i].getVal() if i < len(gauss_ptrs)-1 else (1.-sum([x.getVal() for x in norm_ptrs]))
            g.plotOn(plt, ROOT.RooFit.LineColor(colors[i]), ROOT.RooFit.Normalization(n*yield_, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineStyle(ROOT.kDashed))
        pdf_tot.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed)) # ROOT.RooFit.Normalization(1, ROOT.RooAbsReal.NumEvent)
        histpull = plt.pullHist()
        chi2 = plt.chiSquare()
        pdf_sig.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed)) # ROOT.RooFit.Normalization(1, ROOT.RooAbsReal.NumEvent)
 
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.DrawLatex(0.20, 0.85, label)
        latex.DrawLatex(0.20, 0.80, "q_{T} = [%.1f, %.1f] GeV" % (qTmin_, qTmax_))
        latex.DrawLatex(0.20, 0.75, "Mean = %.3f #pm %.3f" % (hist.GetMean(), hist.GetMeanError()))
        latex.DrawLatex(0.20, 0.70, "RMS = %.3f #pm %.3f" % (hist.GetRMS(), hist.GetRMSError()))
        latex.DrawLatex(0.20, 0.65, "Yield = %.3f #pm %.3f" % (yield_, yield_err.value))
        latex.DrawLatex(0.20, 0.60, "#chi^{2} = %.3f" % chi2)
        if dof > 0 and yield_ > 0:
            latex.DrawLatex(0.20, 0.55, "#color[%d]{Cov. matrix quality = %d}" % (2 if fitRes.covQual() != 3 else 1, fitRes.covQual()))
            latex.DrawLatex(0.20, 0.50, "#color[%d]{MINUIT status = %d}" % (2 if fitRes.status() != 0 else 1, fitRes.status()))
 
        latex.SetTextSize(0.035)
        for i, s in enumerate(mean_ptrs): latex.DrawLatex(0.7, 0.87-i*0.04, "#mu_{%d} = %.3f #pm %.3f" % (i+1, s.getVal(), s.getError() if fitCfg['mean_cfg'][i] >= 0 else mean_ptrs[abs(fitCfg['mean_cfg'][i])-1].getError())) 
        for i, s in enumerate(sigma_ptrs): latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs), "#sigma_{%d} = %.3f #pm %.3f" % (i+1, s.getVal(), s.getError()))
        for i, s in enumerate(norm_ptrs): latex.DrawLatex(0.7, 0.87-i*0.04-0.04*len(mean_ptrs)*2, "n_{%d} = %.3f #pm %.3f" % (i+1, s.getVal(), s.getError()))
        

        plt.Draw("SAME")
        plotter.auxRatio()

        canvas.cd()
        padB.Draw()
        padB.cd()
        padB.SetGrid()
        dummyB.Draw("HIST")

        
        if ratio:
            histRatio = hist.Clone("ratio")
            histRatio.SetMarkerStyle(8)
            histRatio.SetMarkerSize(0.7)
            histRatio.SetMarkerColor(ROOT.kBlack)
            histRatio.SetLineColor(ROOT.kBlack)
            histRatio.SetLineWidth(1)
            norm = ROOT.RooArgSet(recoil)
            histRatio = ratioHist(hist, histRatio, pdf_tot, recoil, norm)
            histRatio.Draw("SAME E0")
        else:
            plt = recoil.frame()
            plt.addPlotable(histpull, "P")
            plt.Draw("SAME")
  
        dummyL.Draw("SAME")
        canvas.Modify()
        canvas.Update()
        canvas.Draw()
        canvas.SaveAs("%s/%03d_recoil.png" % (outDir_, qTbin))
        #canvas.SaveAs("%s/%03d_recoil.pdf" % (outDir_, qTbin))
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)   
        dummyB.Delete()
        dummyT.Delete()
        padT.Delete()
        padB.Delete()
        
        # save output
        out[qTbin] = {}
        out[qTbin]['yield'] = yield_
        out[qTbin]['yield_err'] = yield_err.value
        
        
        for i, s in enumerate(mean_ptrs):
            out[qTbin]['mean%d' % (i+1)], out[qTbin]['mean%d_err' % (i+1)] = s.getVal(), s.getError() if fitCfg['mean_cfg'][i] >= 0 else mean_ptrs[abs(fitCfg['mean_cfg'][i])-1].getError()
        
        for i, s in enumerate(sigma_ptrs):
            out[qTbin]['sigma%d' % (i+1)], out[qTbin]['sigma%d_err' % (i+1)] = s.getVal(), s.getError()
        
        for i, n in enumerate(norm_ptrs):
            out[qTbin]['norm%d' % (i+1)], out[qTbin]['norm%d_err' % (i+1)] = n.getVal(), n.getError()

        out[qTbin]['chi2'] = chi2
        g_chi2.SetPoint(qTbin-1, qT, chi2)
        g_yields.SetPoint(qTbin-1, qT, yield_)
        g_yields.SetPointError(qTbin-1, 0, yield_err)
        
        # update bin condition
        if qTmin_ <= qTmin or qTmin_ == 0:
            qTbin = qTstartBin - 1 # minus 1, redo the fit at qTstartBin for initial conditions
            direction = +1
        
        qTbin += direction # update bin 
        continue
    
    out['qTbins'] = recoil_qTbins_
    with open("%s/results.json" % outDir_, "w") as outfile: json.dump(out, outfile, indent=4)
    plotChi2(g_chi2, "%s/chi2" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)
    plotYields(g_yields, "%s/yields" % outDir_, label, "U_{#perp}" if comp == "perp" else "U_{#parallel}", xMin=qTmin, xMax=qTmax)

    with open("%s/results.json" % outDir_) as f: out = json.load(f)
    doGlobalPlot(fInName, proc, comp, out, recoil_qTbins, outDir_, label, funcJs=funcJs, ratio=ratio, bkgCfg=bkgCfg, qTmin=qTmin, qTmax=qTmax, xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, yRatio=yRatio)


def doGlobalPlot(fInName, proc, comp, jsIn, recoil_qTbins, outDir_, label, funcJs=None, ratio=True, qTmin=0, qTmax=50, bkgCfg={}, xMin=-100, xMax=100, yMin=1e-6, yMax=1e4, yRatio = 1.3):
    
    nGauss = jsIn['nGauss']
    recoil = ROOT.RooRealVar("recoil", "", 0, __min_recoil__, __max_recoil__)

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


def plotChi2(g_chi2, fOut, label, recoilComp, xMin=0, xMax=100, yMin=0, yMax=3):

    cfg_chi2 = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : "q_{T} (GeV)",
        'ytitle'            : "#chi^{2}",
            
        'topRight'          : __topRight__, 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    } 
    
    g_chi2.SetLineColor(ROOT.kRed)
    g_chi2.SetMarkerStyle(8)
    g_chi2.SetMarkerSize(1)
    g_chi2.SetMarkerColor(ROOT.kRed)
    
    f_chi2 = ROOT.TF1("f_chi2", "[0]*x + [1]", 0, 50)
    f_chi2.SetParameters(0, 1)
    f_chi2.SetLineColor(ROOT.kBlue)
    f_chi2.SetLineWidth(2)
    g_chi2.Fit("f_chi2", "NSE", "", 0, 50)
    
    plotter.cfg = cfg_chi2
    canvas = plotter.canvas()
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    g_chi2.Draw("SAME LP")
    f_chi2.Draw("SAME L")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.88, label)
    latex.DrawLatex(0.20, 0.82, recoilComp)
    
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
    

def plotYields(g_yields, fOut, label, recoilComp, xMin=0, xMax=100, yMin=0, yMax=3):

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
    dummy = plotter.dummy()   
    canvas.cd()
    dummy.Draw("HIST")
    g_yields.Draw("SAME LP")
    plotter.aux()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.20, 0.88, label)
    latex.DrawLatex(0.20, 0.82, recoilComp)
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s.png" % fOut)
    canvas.SaveAs("%s.pdf" % fOut)
    canvas.Delete()






def doPlot(comp, param, jsIn, outDir, recoil_qTbins, label, xMin=0, xMax=150, yMin=0, yMax=30, yTitle=""):


    f_base = ROOT.TF1("f_base", jsIn[param]['func'], 0, 500)
    for i in range(0, jsIn[param]['nParams']): f_base.SetParameter(i, jsIn[param]['p%d' % i])
    f_unc = ROOT.TF1("f_unc", jsIn[param]['func'], 0, 500)
    for i in range(0, jsIn[param]['nParams']): f_unc.SetParameter(i, jsIn[param]['p%d' % i])
    
    f_base.SetLineColor(ROOT.kBlack)

    g = ROOT.TGraphErrors()
    #g.SetLineColor(ROOT.kBlack)
    #g.SetMarkerStyle(20)
    #g.SetMarkerSize(0.9)
    #g.SetMarkerColor(ROOT.kBlack)
    #g.SetLineColor(ROOT.kBlack)
    
           
    iPoint = 0
    for qTbin in range(1, len(recoil_qTbins)):
        x = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        x_err = 0.
        y = f_base.Eval(x)
        
        # loop over all stat uncertainties
        y_err = 0
        for iStat in range(0, jsIn['nStatVars']):
            for i in range(0, jsIn[param]['nParams']): 
                if param in jsIn['stat%d' % iStat]: f_unc.SetParameter(i, jsIn['stat%d' % iStat][param]['p%d' % i])
                else : f_unc.SetParameter(i, jsIn[param]['p%d' % i])
            up = f_unc.Eval(x)
            s = up - y
            y_err += s**2
        y_err = y_err**(0.5)    
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
        'ytitleR'           : "Pull",
        'yminR'             : -2.5, # 0.88
        'ymaxR'             : 2.5, # 1.12
    }         

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
  
    
    g.SetFillColor(ROOT.kGreen+1)
    g.Draw("3SAME")
    f_base.Draw("L SAME")
    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.85, label)
    latex.DrawLatex(0.18, 0.80, "U_{#parallel}" if comp == "para" else "U_{#perp}")

    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
        


    line = ROOT.TLine(cfg['xmin'], 0, cfg['xmax'], 0)
    line.SetLineColor(ROOT.kBlue+2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s_refit.png" % (outDir, param))
    canvas.SaveAs("%s/%s_refit.pdf" % (outDir, param))
    canvas.Delete()
    
    
    
def doFit(jsIn, jsOut, comp, param, fitFuncName, iParams, outDir, recoil_qTbins, label, cParams=[], excludeBins=[], fitMin=0, fitMax=200, xMin=0, xMax=150, yMin=0, yMax=30, yTitle="", doFit=True, cutOffMin=-99999, cutOffMax=99999, yRatio = 1.15):

    fitF = getattr(rf, fitFuncName + "_")()
    fit = ROOT.TF1("fit", fitF, 0, 500) if fitF != "" else None
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
    for qTbin in range(1, len(recoil_qTbins)):
        if not str(qTbin) in jsIn:
            #print("WARNING: qTbin %s not found, skip it" % qTbin)
            continue
        if qTbin in excludeBins: continue
        x = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        x_err = 0.5*(recoil_qTbins[qTbin] - recoil_qTbins[qTbin-1])
        if not param in jsIn[str(qTbin)]:
            print("WARNING: param %s not found for qTbin %d" % (param, qTbin))
            continue
        y = jsIn[str(qTbin)][param]
        y_err = jsIn[str(qTbin)][param + "_err"]
        
        if y < cutOffMin: continue
        if y > cutOffMax: continue
        g.SetPoint(iPoint, x, y)   
        g.SetPointError(iPoint, x_err, y_err)
        iPoint += 1

  
    if doFit and fit:
        result = g.Fit(fit.GetName(), "NSE", "", fitMin, fitMax) 
        fit.SetLineColor(ROOT.kRed)
        #fit.GetXaxis().SetRangeUser(0, 200)
        fit.SetLineWidth(2)

        
        cov = result.GetCorrelationMatrix()
        cov.Print()
        
        values = result.GetConfidenceIntervals(0.68, False)
        uncBand = ROOT.TGraphErrors()
        uncBand_ratio = ROOT.TGraphErrors()
        for i in range(len(values)):
            uncBand.SetPoint(i, g.GetX()[i], fit.Eval(g.GetX()[i]))
            uncBand.SetPointError(i, 0, values[i])
            uncBand_ratio.SetPoint(i, g.GetX()[i], 0)
            uncBand_ratio.SetPointError(i, 0, values[i]/fit.Eval(g.GetX()[i]))

        # ratio 
        g_ratio = ROOT.TGraphErrors()
        g_ratio.SetLineColor(ROOT.kBlack)
        g_ratio.SetMarkerStyle(20)
        g_ratio.SetMarkerSize(0.55)
        g_ratio.SetMarkerColor(ROOT.kBlack)
        g_ratio.SetLineColor(ROOT.kBlack)
            
        iPoint = 0
        for qTbin in range(1, len(recoil_qTbins)):
            if not str(qTbin) in jsIn:
                #print("WARNING: qTbin %s not found, skip it" % qTbin)
                continue
            x = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
            x_err = 0.5*(recoil_qTbins[qTbin] - recoil_qTbins[qTbin-1])
            if not param in jsIn[str(qTbin)]: continue
            y = jsIn[str(qTbin)][param]
            y_err = jsIn[str(qTbin)][param + "_err"]
            y_fit = fit.Eval(x)
            #print(x, y, y_fit, y/y_fit)
            if y_err > 0: 
                pull = (y-y_fit)/y_err
                ratio = y / y_fit
                ratio_err = y_err/abs(y)
            else:
                pull = 0
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
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
  
    if doFit and fit:
        uncBand.SetFillColor(ROOT.kGreen+1)
        uncBand.Draw("3SAME")
    g.Draw("PE SAME")
    if fit: fit.Draw("L SAME")  

    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.85, label)
    latex.DrawLatex(0.18, 0.80, "U_{#parallel}" if comp == "para" else "U_{#perp}")
    
    for i in range(0, len(iParams)):
        latex.DrawLatex(0.65, 0.85-i*0.05, "a_{%d} = %.3f #pm %.3f" % (i, fit.GetParameter(i), fit.GetParError(i)))
        
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
        
    if doFit and fit:
        #uncBand_ratio.SetFillColor(ROOT.kGreen+1)
        #uncBand_ratio.Draw("3SAME")
        g_ratio.Draw("PE SAME")

    line = ROOT.TLine(cfg['xmin'], 1, cfg['xmax'], 1)
    line.SetLineColor(ROOT.kBlue+2)
    line.SetLineWidth(2)
    line.Draw("SAME")

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
     
    
    
def doFit_old(jsIn, jsOut, comp, param, fitF, iParams, outDir, recoil_qTbins, label, cParams=[], excludeBins=[], fitMin=0, fitMax=200, xMin=0, xMax=150, yMin=0, yMax=30, yTitle="", doFit=True, cutOffMin=-99999, cutOffMax=99999):

    fit = ROOT.TF1("fit", fitF, 0, 500) if fitF != "" else None
    if cParams == []: cParams = [False]*len(iParams)
    for i, iPar in enumerate(iParams):
        if cParams[i]: fit.FixParameter(i, iPar)
        else: fit.SetParameter(i, iPar)


    g = ROOT.TGraphErrors()
    g.SetLineColor(ROOT.kBlack)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
           
    iPoint = 0
    for qTbin in range(1, len(recoil_qTbins)):
        if not str(qTbin) in jsIn:
            #print("WARNING: qTbin %s not found, skip it" % qTbin)
            continue
        if qTbin in excludeBins: continue
        x = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        x_err = 0.5*(recoil_qTbins[qTbin] - recoil_qTbins[qTbin-1])
        if not param in jsIn[str(qTbin)]:
            print("WARNING: param %s not found for qTbin %d" % (param, qTbin))
            continue
        y = jsIn[str(qTbin)][param]
        y_err = jsIn[str(qTbin)][param + "_err"]
        
        if y < cutOffMin: continue
        if y > cutOffMax: continue
        g.SetPoint(iPoint, x, y)   
        g.SetPointError(iPoint, x_err, y_err)
        iPoint += 1

  
    if doFit and fit:
        result = g.Fit(fit.GetName(), "NSE", "", fitMin, fitMax) 
        fit.SetLineColor(ROOT.kRed)
        #fit.GetXaxis().SetRangeUser(0, 200)
        fit.SetLineWidth(2)

        
        cov = result.GetCorrelationMatrix()
        cov.Print()
        
        values = result.GetConfidenceIntervals(0.68, False)
        uncBand = ROOT.TGraphErrors()
        uncBand_ratio = ROOT.TGraphErrors()
        for i in range(len(values)):
            uncBand.SetPoint(i, g.GetX()[i], fit.Eval(g.GetX()[i]))
            uncBand.SetPointError(i, 0, values[i])
            uncBand_ratio.SetPoint(i, g.GetX()[i], 0)
            uncBand_ratio.SetPointError(i, 0, values[i]/fit.Eval(g.GetX()[i]))

        # ratio 
        g_ratio = ROOT.TGraphErrors()
        g_ratio.SetLineColor(ROOT.kBlack)
        g_ratio.SetMarkerStyle(20)
        g_ratio.SetMarkerSize(0.9)
        g_ratio.SetMarkerColor(ROOT.kBlack)
        g_ratio.SetLineColor(ROOT.kBlack)
            
        iPoint = 0
        for qTbin in range(1, len(recoil_qTbins)):
            if not str(qTbin) in jsIn:
                #print("WARNING: qTbin %s not found, skip it" % qTbin)
                continue
            x = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
            x_err = 0.5*(recoil_qTbins[qTbin] - recoil_qTbins[qTbin-1])
            if not param in jsIn[str(qTbin)]: continue
            y = jsIn[str(qTbin)][param]
            y_err = jsIn[str(qTbin)][param + "_err"]
            y_fit = fit.Eval(x)
            #print(x, y, y_fit, y/y_fit)
            if y_err > 0: pull = (y-y_fit)/y_err
            else: pull = 0
            g_ratio.SetPoint(iPoint, x, pull)
            #g_ratio.SetPointError(iPoint, x_err, y_err/abs(y_fit))
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
        'ytitleR'           : "Pull",
        'yminR'             : -2.5, # 0.88
        'ymaxR'             : 2.5, # 1.12
    }         

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetTickx()
    padT.SetTicky()
    dummyT.Draw("HIST")
  
    if doFit and fit:
        uncBand.SetFillColor(ROOT.kGreen+1)
        uncBand.Draw("3SAME")
    g.Draw("PE SAME")
    if fit: fit.Draw("L SAME")  

    plotter.auxRatio()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.85, label)
    latex.DrawLatex(0.18, 0.80, "U_{#parallel}" if comp == "para" else "U_{#perp}")
    
    for i in range(0, len(iParams)):
        latex.DrawLatex(0.65, 0.85-i*0.05, "a_{%d} = %.3f #pm %.3f" % (i, fit.GetParameter(i), fit.GetParError(i)))
        
    ## BOTTOM PAD ##
    canvas.cd()
    padB.Draw()
    padB.cd()
    dummyB.Draw("HIST")
        
    if doFit and fit:
        uncBand_ratio.SetFillColor(ROOT.kGreen+1)
        uncBand_ratio.Draw("3SAME")
        g_ratio.Draw("PE SAME")

    line = ROOT.TLine(cfg['xmin'], 0, cfg['xmax'], 0)
    line.SetLineColor(ROOT.kBlue+2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/%s.png" % (outDir, param))
    canvas.SaveAs("%s/%s.pdf" % (outDir, param))
    canvas.Delete()
    
    jsOut[param] = {}
    jsOut[param]["func"] = fitF
    jsOut[param]["nParams"] = len(iParams)
    for i, iPar in enumerate(iParams):
        jsOut[param]["p%d" % i] = fit.GetParameter(i)
        jsOut[param]["p%d_err" % i] = fit.GetParError(i)
        jsOut[param]["p%d_isCte" % i] = cParams[i]


  
 
def combinedFit(fInName, proc, comp, fitCfg, jsIn, recoil_qTbins, qTmin=0, qTmax=50, bkgCfg={}):
    
    nGauss = jsIn['nGauss']
    recoil = ROOT.RooRealVar("recoil", "", 0, __min_recoil__, __max_recoil__)
    

    cats = ROOT.RooCategory("category", "") # for each qT bin, define category
    hists = ROOT.std.map("string, RooDataHist*")() # container holding all RooDataHists
    ptrs = {}
    fitFunctions = {}
    garbage = [] # need to store the variables for memory issues
    means_arglist, sigmas_arglist, norms_arglist = [], [], []
    for iGauss in range(1, nGauss+1):
        # prepare variables
        mean_arglist, sigma_arglist = ROOT.RooArgList(), ROOT.RooArgList()
        mean_arglist.setName("gauss%d_mean_arglist" % iGauss)
        sigma_arglist.setName("gauss%d_sigma_arglist" % iGauss)
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            tmp = 'sigma%d_p%d' % (iGauss, iParam)
            ptrs[tmp] = ROOT.RooRealVar(tmp, "", jsIn['sigma%d' % iGauss]['p%d' % iParam], -100, 100) 
            ptrs[tmp].SetName(tmp)
            if jsIn['sigma%d' % iGauss]['p%d_isCte' % iParam]: ptrs[tmp].setConstant(ROOT.kTRUE)
            sigma_arglist.addOwned(ptrs[tmp])
            print(jsIn['sigma%d' % iGauss]['p%d' % iParam])

        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            
            if fitCfg['mean_cfg'][iGauss-1] >= 0:
                ptrs[tmp] = ROOT.RooRealVar(tmp, "", jsIn['mean%d' % iGauss]['p%d' % iParam], -100, 100)
                if jsIn['mean%d' % iGauss]['p%d_isCte' % iParam]: ptrs[tmp].setConstant(ROOT.kTRUE)
            else:
                #mean = ROOT.RooFormulaVar("mean%d" % iGauss, "@0", ROOT.RooArgList(mean_ptrs[abs(fitCfg['mean_cfg'][iGauss-1])-1]))  
                ptrs[tmp] = ROOT.RooFormulaVar(tmp, "@0", ROOT.RooArgList(ptrs['mean%d_p%d' % (abs(fitCfg['mean_cfg'][iGauss-1]), iParam)]))
                

            ptrs[tmp].SetName(tmp)
            mean_arglist.addOwned(ptrs[tmp])

        means_arglist.append(mean_arglist)
        sigmas_arglist.append(sigma_arglist)

        # construct the RooFormulaVar
        fitF = jsIn['sigma%d' % iGauss]['func'].replace("Exp", "ABC").replace("x", "{0}").replace("[", "x[").replace("ABC", "Exp")
        fitFunctions['sigma%d' % iGauss] = fitF
        fitF = jsIn['mean%d' % iGauss]['func'].replace("Exp", "ABC").replace("x", "{0}").replace("[", "x[").replace("ABC", "Exp")
        fitFunctions['mean%d' % iGauss] = fitF
   

    for iGauss in range(1, nGauss):
        # prepare variables
        norm_arglist = ROOT.RooArgList("gauss%d_mean_norm" % iGauss)
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            tmp = 'norm%d_p%d' % (iGauss, iParam)
            ptrs[tmp] = ROOT.RooRealVar(tmp, "", jsIn['norm%d' % iGauss]['p%d' % iParam], -100, 100)
            if jsIn['norm%d' % iGauss]['p%d_isCte' % iParam]: ptrs[tmp].setConstant(ROOT.kTRUE)
            norm_arglist.addOwned(ptrs[tmp])
        norms_arglist.append(norm_arglist)
        
        # construct the RooFormulaVar
        fitF = jsIn['norm%d' % iGauss]['func'].replace("Exp", "ABC").replace("x", "{0}").replace("[", "x[").replace("ABC", "Exp")
        fitFunctions['norm%d' % iGauss] = fitF

   

    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")

    out = {}
    out['nGauss'] = nGauss
    nBins = 0
    pdf_tot = ROOT.RooSimultaneous("pdf_tot", "", cats) # total pdf, containing all the categories
    for qTbin in range(1, len(recoil_qTbins)):
    
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        if qT > qTmax or qT < qTmin: continue
        qTmin_, qTmax_ = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        nBins += 1
        
        qTbins_ = list(functions.drange(qTmin_, qTmax_, 0.5))
        hist = None
        for b in qTbins_:
            iBin = hIn_x.FindBin(b)
            h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
            if hist == None: 
                hist = copy.deepcopy(h)
                hist.SetName("hist")
            else: hist.Add(h)

        yield_ = hist.Integral()
        rdh = ROOT.RooDataHist("", "", ROOT.RooArgList(recoil), ROOT.RooFit.Import(hist))
        rdh.SetName("rdh_bin%d" % (qTbin))
        catIDx = "cat_bin%s" % (qTbin)
        hists.insert(ROOT.std.pair("string, RooDataHist*")(catIDx, rdh))
        cats.defineType(catIDx, qTbin)    
        
        print("----------------------------------->", rdh.numEntries(), rdh.sumEntries())
        ptrs[rdh.GetName()] = rdh

        pdfs, norms = ROOT.RooArgList(), ROOT.RooArgList()
        pdfs.setName("pdfs_bin%d" % qTbin)
        norms.setName("norms_bin%d" % qTbin)        
        for iGauss in range(1, nGauss+1):
        
            mean_ = ROOT.RooFormulaVar("mean%d_bin%d" % (iGauss, qTbin), fitFunctions['mean%d' % iGauss].format(qT), means_arglist[iGauss-1])
            sigma_ = ROOT.RooFormulaVar("sigma%d_bin%d" % (iGauss, qTbin), fitFunctions['sigma%d' % iGauss].format(qT), sigmas_arglist[iGauss-1])
            gauss = ROOT.RooGaussian("gauss%d_bin%d" % (iGauss, qTbin), "", recoil, mean_, sigma_)
            pdfs.addOwned(gauss)
            garbage.append(mean_)
            garbage.append(sigma_)
            garbage.append(gauss)
 
        for iGauss in range(1, nGauss):
            norm_ = ROOT.RooFormulaVar("norm%d_bin%d" % (iGauss, qTbin), fitFunctions['norm%d' % iGauss].format(qT), norms_arglist[iGauss-1])
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
        #garbage.append(pdf_tot_)
        ptrs[pdf_tot_.GetName()] = pdf_tot_
        pdf_tot_.Print()
        pdf_tot.addPdf(pdf_tot_, "cat_bin%s" % (qTbin))
        
        print("----------------------------------->", rdh.numEntries(), rdh.sumEntries())


    #recoil.setBins(600) # seems to be necessary ??
    rdh_tot = ROOT.RooDataHist("rdh_tot", "", ROOT.RooArgList(recoil), cats, hists) # total RDH   , ROOT.RooFit.Minos(ROOT.kTRUE) , , ROOT.RooFit.Hesse(ROOT.kTRUE) 
    # ROOT.RooFit.SumW2Error(ROOT.kTRUE) --> errors calculated according to inverse Hessian method (see https://root-forum.cern.ch/t/using-roofit-with-weighted-data/20207/2)
    # ROOT.RooFit.Minos(ROOT.kTRUE) does not work with weighted values
    #fitRes = pdf_tot.fitTo(rdh_tot,  ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kTRUE)) # ROOT.RooFit.SumW2Error(ROOT.kTRUE), # , ROOT.RooFit.Minos(ROOT.kTRUE)
    
    
  
    fitRes = pdf_tot.fitTo(rdh_tot,  ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    #fitRes = pdf_tot.fitTo(rdh_tot,  ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.SumW2Error(ROOT.kFALSE) , ROOT.RooFit.Minos(ROOT.kTRUE)) # ROOT.RooFit.AsymptoticError(ROOT.kTRUE) ROOT.RooFit.SumW2Error(ROOT.kTRUE) , ROOT.RooFit.Minos(ROOT.kTRUE)
    
    ## https://root.cern.ch/doc/master/classRooAbsPdf.html#a5f79f16f4a26a19c9e66fb5c080f59c5  ROOT.RooFit.AsymptoticError(ROOT.kTRUE) , ROOT.RooFit.Strategy(2)
    fitValid = (fitRes.covQual() == 3 and fitRes.status() == 0)
    print("****************************")
    print("FIT STATUS")
    print("Covariance Quality = %d" % fitRes.covQual())
    print("Fit status = %d" % fitRes.status())
    print("****************************")
    
    cov = fitRes.covarianceMatrix()
    cov.Print()

    #
    
    floatingParams = fitRes.floatParsFinal()
    jsOut = copy.deepcopy(jsIn)
    for iGauss in range(1, nGauss+1):
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            tmp = 'sigma%d_p%d' % (iGauss, iParam)
            p = getParamIdx(tmp, floatingParams)
            print(p)
            if p == -1: val, err = ptrs[tmp].getVal(), 0
            else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            jsOut['sigma%d' % iGauss]['p%d' % iParam] = val
            jsOut['sigma%d' % iGauss]['p%d_err' % iParam] = err
            print("%s \t %.3e" % (tmp, val))
          
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            p = getParamIdx(tmp, floatingParams)
            if p == -1: val, err = ptrs[tmp].getVal(), 0
            else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            if fitCfg['mean_cfg'][iGauss-1] < 0:
                idx = abs(fitCfg['mean_cfg'][iGauss-1])
                tmp_1 = 'mean%d_p%d' % (idx, iParam)
                p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                if p == -1: val, err = ptrs[tmp_1].getVal(), 0
                else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            jsOut['mean%d' % iGauss]['p%d' % iParam] = val
            jsOut['mean%d' % iGauss]['p%d_err' % iParam] = err
            #print("%s \t %.3e" % (tmp, val))
    
    for iGauss in range(1, nGauss):
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            tmp = 'norm%d_p%d' % (iGauss, iParam)
            p = getParamIdx(tmp, floatingParams)
            if p == -1: val, err = ptrs[tmp].getVal(), 0
            else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            jsOut['norm%d' % iGauss]['p%d' % iParam] = val
            jsOut['norm%d' % iGauss]['p%d_err' % iParam] = err
            #print("%s \t %.3e" % (tmp, val))
    
    fIn.Close()
    
    # diagonalize covariance matrix and store perturbations
    variations = diagonalize(fitRes, verbose=False)
    jsOut['nStatVars'] = len(variations)
    for iVar, var in enumerate(variations):
        
        #print(iVar, var)
        jsOut['stat%d' % iVar] = {}
        t = jsOut['stat%d' % iVar]
        
        for iGauss in range(1, nGauss+1):
            t['sigma%d' % iGauss] = {}
            for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
                tmp = 'sigma%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['sigma%d' % iGauss]['p%d' % iParam] = val
            
            t['mean%d' % iGauss] = {}  
            for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
                tmp = 'mean%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                if fitCfg['mean_cfg'][iGauss-1] < 0:
                    idx = abs(fitCfg['mean_cfg'][iGauss-1])
                    tmp_1 = 'mean%d_p%d' % (idx, iParam)
                    p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                    if p == -1: val = ptrs[tmp_1].getVal()
                    else: val = var[p]
                    #print(p, val)
                t['mean%d' % iGauss]['p%d' % iParam] = val

        for iGauss in range(1, nGauss):
            t['norm%d' % iGauss] = {}
            for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
                tmp = 'norm%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['norm%d' % iGauss]['p%d' % iParam] = val
 
    return jsOut
 



def combinedFit_zlib(fInName, proc, comp, fitCfg, jsIn, recoil_qTbins, qTmin=0, qTmax=50, bkgCfg={}):
        
    # prepare file and binning
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")
    hIn_y = hIn.ProjectionY("hIn_y")
    
    recoil_bins = [hIn_y.GetBinLowEdge(i) for i in range(1, hIn_y.GetNbinsX()+2)]
    centers = [hIn_y.GetBinCenter(i) for i in range(1, hIn_y.GetNbinsX()+1)]
    bins = zfit.binned.VariableBinning(recoil_bins, name="recoil_bins")
    recoil = zfit.Space('recoil_bins', limits=(min(recoil_bins), max(recoil_bins)))
    #recoil_binned = zfit.Space('recoil_bins', binning=bins)
    recoil_binned = zfit.data.BinnedData.from_hist(narf.root_to_hist(hIn_y, axis_names=["recoil_bins"]) ).space
    nGauss = jsIn['nGauss']


    # construct Gauss
    fitFunctions, params_qT, params = {}, {}, {}
    for iGauss in range(1, nGauss+1):

        params_ = []
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            tmp = 'sigma%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            param = zfit.Parameter(tmp, float(jsIn['sigma%d' % iGauss]['p%d' % iParam]), step_size=0.001)
            if jsIn['sigma%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            params_.append(param)
        params["sigma_gauss%d" % iGauss] = params_
    
        params_ = []
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            if fitCfg['mean_cfg'][iGauss-1] >= 0:
                param = zfit.Parameter(tmp, float(jsIn['mean%d' % iGauss]['p%d' % iParam]), step_size=0.001)
                if jsIn['mean%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            else:
                dependent_tensor = params_[abs(fitCfg['mean_cfg'][iGauss-1]-1)]
                param = zfit.ComposedParameter(tmp, dependent_tensor)
            params_.append(param)
        params["mean_gauss%d" % iGauss] = params_

    for iGauss in range(1, nGauss):
        params_ = []
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            tmp = 'norm%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            param = zfit.Parameter(tmp, float(jsIn['norm%d' % iGauss]['p%d' % iParam]), step_size=0.001)
            if jsIn['norm%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            params_.append(param)
        params["norm_gauss%d" % iGauss] = params_
        




    models, datas = [], []
    nll_simultaneous = None
    nBins = 0
    for qTbin in range(1, len(recoil_qTbins)):
    
        if qTbin >= len(recoil_qTbins): break
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        qTlow, qThigh = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        if qTlow >= qTmax: break
        print("##### DO BIN", qTbin, qT, qTlow, "->",  qThigh)
        nBins += 1
        
        qTp = zfit.Parameter("qT%s" % nBins, qT)
        qTp.floating = False
        params_qT["qT%s" % nBins] = qTp
        
       
        qTbins_ = list(functions.drange(qTlow, qThigh, 0.5))
        hist = None
        for b in qTbins_:
            iBin = hIn_x.FindBin(b)
            h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
            if hist == None: 
                hist = copy.deepcopy(h)
                hist.SetName("hist")
            else: hist.Add(h)

        
        yield_err = ctypes.c_double(1.)
        yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        yield_err = yield_err.value

        h = narf.root_to_hist(hist, axis_names=["recoil_bins"]) 
        data__ = zfit.data.BinnedData.from_hist(h)
        data_vals = data__.values()
        data_vars = data__.variances()

        gauss_pdfs, gauss_norms = [], []      
        for iGauss in range(1, nGauss+1):
        
            # sigma
            func = getattr(rf, jsIn['sigma%d' % iGauss]['funcName'])
            params_ = [qTp] + params["sigma_gauss%d" % iGauss]
            p_sigma = zfit.ComposedParameter('sigma_gauss%d_qT%d' % (iGauss, qTbin), func, params=params_)

            # mean
            func = getattr(rf, jsIn['mean%d' % iGauss]['funcName'])
            params_ = [qTp] + params["mean_gauss%d" % iGauss]
            p_mean = zfit.ComposedParameter('mean_gauss%d_qT%d' % (iGauss, qTbin), func, params=params_)
            
            # construct gauss
            gauss = zfit.pdf.Gauss(obs=recoil, mu=p_mean, sigma=p_sigma, name="gauss%d"%iGauss)
            gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss, recoil_binned)
            gauss_pdfs.append(gauss_binned)
            
            # norm
            if iGauss == nGauss: continue
            func = getattr(rf, jsIn['norm%d' % iGauss]['funcName'])
            params_ = [qTp] + params["norm_gauss%d" % iGauss]
            p_norm = zfit.ComposedParameter('norm_gauss%d_qT%d' % (iGauss, qTbin), func, params=params_)
            gauss_norms.append(p_norm)
       
        model = zfit.pdf.BinnedSumPDF(gauss_pdfs, fracs=gauss_norms, obs=recoil_binned, name="nll_gauss%d" % iGauss)
        #print(model.values())
        #sys.exit()
        models.append(model)
        datas.append(data__)
        
        #print(model)
        
        
        #print(data__.space.binning)
        #print(model.space.binning)
        #self.space = space
        #self.values = values
        #self.variances = variances
 
        
        continue
        print("Construct NLL")
        nll = zfit.loss.BinnedNLL(model=model, data=data__, options={"subtr_const" : False})
        #nlls.append(nll)
        print("Done")
        
        #if nll_simultaneous == None: nll_simultaneous = nll
        #else: nll_simultaneous += nll

        sys.exit()
        break
      
            
            
            

        
        '''
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
        #garbage.append(pdf_tot_)
        ptrs[pdf_tot_.GetName()] = pdf_tot_
        pdf_tot_.Print()
        pdf_tot.addPdf(pdf_tot_, "cat_bin%s" % (qTbin))
        
        print("----------------------------------->", rdh.numEntries(), rdh.sumEntries())  
        '''
    
    # minimize
    
    start = time.time()
    print("##### START BinnedNLL")
    nll = zfit.loss.BinnedNLL(model=models, data=datas) # , options={"subtr_const" : False}
    end = time.time()
    print("##### END BinnedNLL, DURATION", (end-start))

    start = time.time()
    print("##### START MIMIMIZATION")
    zfit.settings.set_seed(1)
    minimizer = zfit.minimize.Minuit(verbosity=7) # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1(verbosity=7) # tol=1e-4, verbosity=5, gradient=True
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### END MIMIMIZATION, DURATION", (end-start))
    print(result.params)

            

    errCalc = "hesse"
    if errCalc == "hesse":
        start = time.time()
        print("##### START HESSE")
        param_errors = result.hesse(name='hesse')
        end = time.time()
        print("##### END HESSE, DURATION", (end-start))
            
    if errCalc == "minos":
        start = time.time()
        print("##### START ERRORS")
        param_errors, _ = result.errors(name='minos', method='minuit_minos')
        end = time.time()
        print("##### END HESSE, DURATION", (end-start))
        for p in param_errors:
            print(param_errors[p]['upper'], param_errors[p]['lower'])
            param_errors[p]['error'] = 0.5*(abs(param_errors[p]['lower']) + abs(param_errors[p]['upper']))

    print(param_errors)
    
    
    
    
    jsOut = copy.deepcopy(jsIn)
    for iGauss in range(1, nGauss+1):
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            param = params["sigma_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['sigma%d' % iGauss]['p%d' % iParam] = val
            jsOut['sigma%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_sigma%d \t %.3e" % (iGauss, iParam, val))
            
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            param = params["mean_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['mean%d' % iGauss]['p%d' % iParam] = val
            jsOut['mean%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_mean%d \t %.3e" % (iGauss, iParam, val))
          
        '''
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            p = getParamIdx(tmp, floatingParams)
            if p == -1: val, err = ptrs[tmp].getVal(), 0
            else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            if fitCfg['mean_cfg'][iGauss-1] < 0:
                idx = abs(fitCfg['mean_cfg'][iGauss-1])
                tmp_1 = 'mean%d_p%d' % (idx, iParam)
                p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                if p == -1: val, err = ptrs[tmp_1].getVal(), 0
                else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            jsOut['mean%d' % iGauss]['p%d' % iParam] = val
            jsOut['mean%d' % iGauss]['p%d_err' % iParam] = err
            #print("%s \t %.3e" % (tmp, val))
        '''
    for iGauss in range(1, nGauss):
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            param = params["norm_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['norm%d' % iGauss]['p%d' % iParam] = val
            jsOut['norm%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_norm%d \t %.3e" % (iGauss, iParam, val))
    
    fIn.Close()
    
    '''
    # diagonalize covariance matrix and store perturbations
    variations = diagonalize(fitRes, verbose=False)
    jsOut['nStatVars'] = len(variations)
    for iVar, var in enumerate(variations):
        
        #print(iVar, var)
        jsOut['stat%d' % iVar] = {}
        t = jsOut['stat%d' % iVar]
        
        for iGauss in range(1, nGauss+1):
            t['sigma%d' % iGauss] = {}
            for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
                tmp = 'sigma%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['sigma%d' % iGauss]['p%d' % iParam] = val
            
            t['mean%d' % iGauss] = {}  
            for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
                tmp = 'mean%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                if fitCfg['mean_cfg'][iGauss-1] < 0:
                    idx = abs(fitCfg['mean_cfg'][iGauss-1])
                    tmp_1 = 'mean%d_p%d' % (idx, iParam)
                    p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                    if p == -1: val = ptrs[tmp_1].getVal()
                    else: val = var[p]
                    #print(p, val)
                t['mean%d' % iGauss]['p%d' % iParam] = val

        for iGauss in range(1, nGauss):
            t['norm%d' % iGauss] = {}
            for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
                tmp = 'norm%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['norm%d' % iGauss]['p%d' % iParam] = val
    
    '''
    return jsOut
 





def combinedFit_zlib_cond(fInName, proc, comp, fitCfg, jsIn, recoil_qTbins, qTmin=0, qTmax=50, bkgCfg={}):
    
    # prepare file and binning
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")
    hIn_y = hIn.ProjectionY("hIn_y")
    
    
    # construct 2D PDF
    # condition it according to qT
    recoil_qTbins += [100] ## FIX
    
    qt_bins = zfit.binned.VariableBinning(recoil_qTbins, name="qt_space_binned")
    qt_space_binned = zfit.Space('qt_space_binned', binning=qt_bins)
    qt_space = zfit.Space('qt_space', limits=(min(recoil_qTbins), max(recoil_qTbins)))
    
    

    recoil_bins_list = [hIn_y.GetBinLowEdge(i) for i in range(1, hIn_y.GetNbinsX()+2)]
    #centers = [hIn_y.GetBinCenter(i) for i in range(1, hIn_y.GetNbinsX()+1)]
    recoil_bins = zfit.binned.VariableBinning(recoil_bins_list, name="recoil_space_binned")
    recoil_space_binned = zfit.Space('recoil_space_binned', binning=recoil_bins)
    recoil_space = zfit.Space('recoil_space', limits=(-500, 500))
    #recoil_space_binned = zfit.data.BinnedData.from_hist(narf.root_to_hist(hIn_y, axis_names=["recoil_bins"])).space
    nGauss = jsIn['nGauss']
    

    # construct Gauss
    fitFunctions, params_qT, params = {}, {}, {}
    for iGauss in range(1, nGauss+1):

        params_ = []
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            tmp = 'sigma%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            param = zfit.Parameter(tmp, float(jsIn['sigma%d' % iGauss]['p%d' % iParam]), step_size=0.001)
            if jsIn['sigma%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            params_.append(param)
        params["sigma_gauss%d" % iGauss] = params_
    
        params_ = []
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            if fitCfg['mean_cfg'][iGauss-1] >= 0:
                param = zfit.Parameter(tmp, float(jsIn['mean%d' % iGauss]['p%d' % iParam]), step_size=0.001)
                if jsIn['mean%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            else:
                dependent_tensor = params_[abs(fitCfg['mean_cfg'][iGauss-1]-1)]
                param = zfit.ComposedParameter(tmp, dependent_tensor)
            params_.append(param)
        params["mean_gauss%d" % iGauss] = params_

    for iGauss in range(1, nGauss):
        params_ = []
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            tmp = 'norm%d_p%d' % (iGauss, iParam)
            tmp_ = 'p%d' % iParam
            param = zfit.Parameter(tmp, float(jsIn['norm%d' % iGauss]['p%d' % iParam]), step_size=0.001)
            if jsIn['norm%d' % iGauss]['p%d_isCte' % iParam]: param.floating = False
            params_.append(param)
        params["norm_gauss%d" % iGauss] = params_
        



    qT = zfit.Parameter("qT", 5)
    gauss_pdfs, gauss_norms = [], []      
    for iGauss in range(1, nGauss+1):
        
        # sigma
        func = getattr(rf, jsIn['sigma%d' % iGauss]['funcName'])
        params_ = [qT] + params["sigma_gauss%d" % iGauss]
        p_sigma = zfit.ComposedParameter('sigma_gauss%d' % (iGauss), func, params=params_)

        # mean
        func = getattr(rf, jsIn['mean%d' % iGauss]['funcName'])
        params_ = [qT] + params["mean_gauss%d" % iGauss]
        p_mean = zfit.ComposedParameter('mean_gauss%d' % (iGauss), func, params=params_)
            
        # construct gauss
        gauss = zfit.pdf.Gauss(obs=recoil_space, mu=p_mean, sigma=p_sigma, name="gauss%d"%iGauss)
        gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss, recoil_space_binned)
        gauss_pdfs.append(gauss)
            
        # norm
        if iGauss == nGauss: continue
        func = getattr(rf, jsIn['norm%d' % iGauss]['funcName'])
        params_ = [qT] + params["norm_gauss%d" % iGauss]
        p_norm = zfit.ComposedParameter('norm_gauss%d' % (iGauss), func, params=params_)
        gauss_norms.append(p_norm)

    #model = zfit.pdf.BinnedSumPDF(gauss_pdfs, fracs=gauss_norms, name="recoil_pdf") # obs=recoil_space, 
    model = zfit.pdf.SumPDF(gauss_pdfs, fracs=gauss_norms, name="recoil_pdf") # obs=recoil_space, 
    cond = {}
    cond[qT] = qt_space
    model_cond = zfit.pdf.ConditionalPDFV1(model, cond, name="recoil_pdf_cond")
    model_cond_binned = model_cond.to_binned(recoil_space_binned*qt_space_binned)
    
    # the uncerlying PDFs must have binned obs, to be compatible with the data
    
    print(model_cond_binned)
    
    #def to_binned(self, space, *, extended=None, norm=None):
    #model_cond.toBinned()
    #model_cond_binned = zfit.pdf.BinnedFromUnbinnedPDF(model_cond, recoil_space)
    
    #sys.exit()





    ################################
    models, datas = [], []
    hist_sliced = ROOT.TH2D("hist_sliced", "", len(recoil_bins_list)-1, min(recoil_bins_list), max(recoil_bins_list), len(recoil_qTbins)-1, min(recoil_qTbins), max(recoil_qTbins))
    for qTbin in range(1, len(recoil_qTbins)):
        for recoilBin in range(1, len(recoil_bins_list)):
        
            hist_sliced.SetBinContent(recoilBin, qTbin, hIn.GetBinContent(qTbin, recoilBin))
            
        #continue
    
        if qTbin >= len(recoil_qTbins): break
        qT = 0.5*(recoil_qTbins[qTbin-1] + recoil_qTbins[qTbin])
        qTlow, qThigh = recoil_qTbins[qTbin-1], recoil_qTbins[qTbin]
        if qTlow >= qTmax: break
        print("##### DO BIN", qTbin, qT, qTlow, "->",  qThigh)
 
       
        qTbins_ = list(functions.drange(qTlow, qThigh, 0.5))
        hist = None
        for b in qTbins_:
            iBin = hIn_x.FindBin(b)
            h = hIn.ProjectionY("tmp_%d" % iBin, iBin, iBin)
            if hist == None: 
                hist = copy.deepcopy(h)
                hist.SetName("hist")
            else: hist.Add(h)

        
        yield_err = ctypes.c_double(1.)
        yield_ = hist.IntegralAndError(0, hist.GetNbinsX() + 1, yield_err)
        yield_err = yield_err.value

        h = narf.root_to_hist(hist, axis_names=["recoil_bins"]) 
        data__ = zfit.data.BinnedData.from_hist(h)
        datas.append(data__)
   
    
    h = narf.root_to_hist(hist_sliced, axis_names=["recoil_space", "qt_space"]) 
    print(h)
    
    data__ = zfit.data.BinnedData.from_hist(h)
    
    data__unbinned = data__.to_unbinned()

    #sys.exit()

    #models = [model_cond_binned]*(len(recoil_qTbins)-1)
    #nll = zfit.loss.BinnedNLL(model=model_cond_binned, data=data__) # , options={"subtr_const" : False}
    start = time.time()
    print("##### START UnbinnedNLL")
    nll = zfit.loss.UnbinnedNLL(model=model_cond, data=data__unbinned) # , options={"subtr_const" : False}
    end = time.time()
    print("##### END UnbinnedNLL, DURATION", (end-start))
    print(nll)
    
    start = time.time()
    print("##### START MIMIMIZATION")
    minimizer = zfit.minimize.Minuit(verbosity=10) # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1(verbosity=10) # tol=1e-4, verbosity=5, gradient=True
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### END MIMIMIZATION, DURATION", (end-start))

    
    
    
    print(result)
    print(result.params)
    
    
    '''
    
    
    sys.exit()
    
    minimizer = zfit.minimize.Minuit(verbosity=7) # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1(verbosity=7) # tol=1e-4, verbosity=5, gradient=True
    result = minimizer.minimize(nll)
    
    
    
    
    
    sys.exit()
    
    # minimize
    
    start = time.time()
    print("##### START BinnedNLL")
    nll = zfit.loss.BinnedNLL(model=models, data=datas) # , options={"subtr_const" : False}
    end = time.time()
    print("##### END BinnedNLL, DURATION", (end-start))

    start = time.time()
    print("##### START MIMIMIZATION")
    zfit.settings.set_seed(1)
    minimizer = zfit.minimize.Minuit(verbosity=7) # mode=0, gradient=True tol=1e-4, verbosity=5, gradient=True
    #minimizer = zfit.minimize.ScipyTrustConstrV1(verbosity=7) # tol=1e-4, verbosity=5, gradient=True
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### END MIMIMIZATION, DURATION", (end-start))
    print(result.params)

            

    errCalc = "hesse"
    if errCalc == "hesse":
        start = time.time()
        print("##### START HESSE")
        param_errors = result.hesse(name='hesse')
        end = time.time()
        print("##### END HESSE, DURATION", (end-start))
            
    if errCalc == "minos":
        start = time.time()
        print("##### START ERRORS")
        param_errors, _ = result.errors(name='minos', method='minuit_minos')
        end = time.time()
        print("##### END HESSE, DURATION", (end-start))
        for p in param_errors:
            print(param_errors[p]['upper'], param_errors[p]['lower'])
            param_errors[p]['error'] = 0.5*(abs(param_errors[p]['lower']) + abs(param_errors[p]['upper']))

    print(param_errors)
    
    '''
    
    
    err = 0
    
    jsOut = copy.deepcopy(jsIn)
    for iGauss in range(1, nGauss+1):
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            param = params["sigma_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            #err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['sigma%d' % iGauss]['p%d' % iParam] = val
            jsOut['sigma%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_sigma%d \t %.3e" % (iGauss, iParam, val))
            
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            param = params["mean_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            #err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['mean%d' % iGauss]['p%d' % iParam] = val
            jsOut['mean%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_mean%d \t %.3e" % (iGauss, iParam, val))
          
        '''
        for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
            tmp = 'mean%d_p%d' % (iGauss, iParam)
            p = getParamIdx(tmp, floatingParams)
            if p == -1: val, err = ptrs[tmp].getVal(), 0
            else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            if fitCfg['mean_cfg'][iGauss-1] < 0:
                idx = abs(fitCfg['mean_cfg'][iGauss-1])
                tmp_1 = 'mean%d_p%d' % (idx, iParam)
                p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                if p == -1: val, err = ptrs[tmp_1].getVal(), 0
                else: val, err = floatingParams[p].getVal(), floatingParams[p].getError()
            jsOut['mean%d' % iGauss]['p%d' % iParam] = val
            jsOut['mean%d' % iGauss]['p%d_err' % iParam] = err
            #print("%s \t %.3e" % (tmp, val))
        '''
    for iGauss in range(1, nGauss):
        for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
            param = params["norm_gauss%d" % iGauss][iParam]
            val = zfit.run(param)
            #err = param_errors[param]['error'] if param in param_errors else 0
            jsOut['norm%d' % iGauss]['p%d' % iParam] = val
            jsOut['norm%d' % iGauss]['p%d_err' % iParam] = err
            print("gauss%d_norm%d \t %.3e" % (iGauss, iParam, val))
    
    fIn.Close()
    
    '''
    # diagonalize covariance matrix and store perturbations
    variations = diagonalize(fitRes, verbose=False)
    jsOut['nStatVars'] = len(variations)
    for iVar, var in enumerate(variations):
        
        #print(iVar, var)
        jsOut['stat%d' % iVar] = {}
        t = jsOut['stat%d' % iVar]
        
        for iGauss in range(1, nGauss+1):
            t['sigma%d' % iGauss] = {}
            for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
                tmp = 'sigma%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['sigma%d' % iGauss]['p%d' % iParam] = val
            
            t['mean%d' % iGauss] = {}  
            for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
                tmp = 'mean%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                if fitCfg['mean_cfg'][iGauss-1] < 0:
                    idx = abs(fitCfg['mean_cfg'][iGauss-1])
                    tmp_1 = 'mean%d_p%d' % (idx, iParam)
                    p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                    if p == -1: val = ptrs[tmp_1].getVal()
                    else: val = var[p]
                    #print(p, val)
                t['mean%d' % iGauss]['p%d' % iParam] = val

        for iGauss in range(1, nGauss):
            t['norm%d' % iGauss] = {}
            for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
                tmp = 'norm%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['norm%d' % iGauss]['p%d' % iParam] = val
    
    '''
    return jsOut
 






def combinedFit_scipy(fInName, proc, comp, fitCfg, jsIn, recoil_qTbins, qTmin=0, qTmax=50, bkgCfg={}):
    
    # prepare file and binning
    fIn = ROOT.TFile(fInName)
    hIn = fIn.Get("%s_%s" % (proc, comp))
    hIn_x = hIn.ProjectionX("hIn_x")
    hIn_y = hIn.ProjectionY("hIn_y")
    

    recoil_qTbins_ = []
    for qT in recoil_qTbins:
        if qT <= qTmax: recoil_qTbins_.append(qT)
        else: break
    recoil_qTbins = recoil_qTbins_
    recoil_bins_list = [hIn_y.GetBinLowEdge(i) for i in range(1, hIn_y.GetNbinsX()+2)]
    nGauss = jsIn['nGauss']
    

    ################################
    hist_sliced = ROOT.TH2D("hist_sliced", "", len(recoil_bins_list)-1, min(recoil_bins_list), max(recoil_bins_list), len(recoil_qTbins)-1, min(recoil_qTbins), max(recoil_qTbins))
    for qTbin in range(1, len(recoil_qTbins)):
        for recoilBin in range(1, len(recoil_bins_list)):
            hist_sliced.SetBinContent(recoilBin, qTbin, hIn.GetBinContent(qTbin, recoilBin))
            hist_sliced.SetBinError(recoilBin, qTbin, hIn.GetBinError(qTbin, recoilBin))
            #print(recoilBin, qTbin, hIn.GetBinContent(qTbin, recoilBin), hIn.GetBinError(qTbin, recoilBin))
  
    h = narf.root_to_hist(hist_sliced, axis_names=["recoil_space", "qt_space"]) 
    print(h)    
    
    # initial parameters
    sigma_init = []
    for iGauss in range(1, nGauss+1):
        for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
            #if iParam == 0: continue
            sigma_init.append(float(jsIn['sigma%d' % iGauss]['p%d' % iParam]))
            
    
        
   
    init_parms = np.array(sigma_init)
    #init_parms = np.array([12, 0, 4, 0, 6, 0, 8, 0]) # lin
    #init_parms = np.array([12, 0, 0, 4, 0, 0, 6, 0, 0, 8, 0, 0]) # quadr
    #init_parms = np.array([12, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 0, 8, 0, 0, 0]) # cube
    init_parms = np.array([0, 1, 12, 0, 1, 4, 0, 1, 6, 0, 1, 8., 0.08, 0.20, 0.25]) # power   p0*tf.math.pow(qT, p1) + p2
    

    
    func_name = "mc_perp"
    func_model = getattr(rf, func_name)


    # do refit
    start = time.time()
    res = fitter.fit_hist(h, func_model, init_parms, norm_axes=[0], mode = "chisq_normalized", max_iter = 5, edmtol = 1e-5)
    dt_min = time.time() - start
    
    print(res)
    print("Duration", dt_min)

    fit_parms = res['x']
    
    
    # save all parameters to dict
    jsOut = {}
    jsOut["func_name"] = func_name
    jsOut["nParams"] = len(fit_parms)
    for iPar, par in enumerate(fit_parms):
        jsOut["p%d"%iPar] = par
    

    
    fIn.Close()
    

    '''
    # diagonalize covariance matrix and store perturbations
    variations = diagonalize(fitRes, verbose=False)
    
    for iVar, var in enumerate(variations):
        
        #print(iVar, var)
        jsOut['stat%d' % iVar] = {}
        t = jsOut['stat%d' % iVar]
        
        for iGauss in range(1, nGauss+1):
            t['sigma%d' % iGauss] = {}
            for iParam in range(0, jsIn['sigma%d' % iGauss]['nParams']):
                tmp = 'sigma%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['sigma%d' % iGauss]['p%d' % iParam] = val
            
            t['mean%d' % iGauss] = {}  
            for iParam in range(0, jsIn['mean%d' % iGauss]['nParams']):
                tmp = 'mean%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                if fitCfg['mean_cfg'][iGauss-1] < 0:
                    idx = abs(fitCfg['mean_cfg'][iGauss-1])
                    tmp_1 = 'mean%d_p%d' % (idx, iParam)
                    p = getParamIdx(tmp_1, floatingParams) # try first Gauss
                    if p == -1: val = ptrs[tmp_1].getVal()
                    else: val = var[p]
                    #print(p, val)
                t['mean%d' % iGauss]['p%d' % iParam] = val

        for iGauss in range(1, nGauss):
            t['norm%d' % iGauss] = {}
            for iParam in range(0, jsIn['norm%d' % iGauss]['nParams']):
                tmp = 'norm%d_p%d' % (iGauss, iParam)
                p = getParamIdx(tmp, floatingParams)
                if p == -1: val = ptrs[tmp].getVal()
                else: val = var[p]
                t['norm%d' % iGauss]['p%d' % iParam] = val
    
    '''

    return jsOut
 

