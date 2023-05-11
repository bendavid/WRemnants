
import sys,array,math,os
import numpy as np
import ctypes
import json
import copy

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


import utils
import plotter
import recoilLibs_scipy as rls
from wremnants.datasets.datagroups import datagroups2016

import lz4.frame
import narf



     

def ewk_perp():

    procs = ["EWK"]
    tag = "ewk_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    recoilMin, recoilMax = -150, 150
    procNames = groups.getProcNames(to_expand=procs)
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    rebin=2
    
    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    if True:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-2, yMax=2)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "ewk_perp"
    fitCfg['func_parms_vals'] = [40, 15, 5]
    fitCfg['func_parms_cfg'] = [0, 0, 0]
    fitCfg['func_parms_nparms'] = [3, 0, 0] # sigma, mean, norm 
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin, yRatio=1.25, sumw2=False) 
        return 
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "cpol2", [60, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "cpol4", [60, 0, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, cutOffMax=80, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "cpol4", [40, 0, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, cutOffMax=80, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")
        
        rls.addParam(jsOut, "p3", "[0]", [0.03])
        rls.addParam(jsOut, "p4", "[0]", [0.30])
     
        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return
        

    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "ewk_perp_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin, sumw2=False)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, statUnc=False, rebin=rebin, yRatio=1.25)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")  

def ewk_para():

    procs = ["EWK"]
    tag = "ewk_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    bhist = readProc(groups, "recoil_corr_xy_para_qTbinned", procs)
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]
    
    recoilMin, recoilMax = -150, 150
    rebin=2
    if True:
        
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.19057e+00, -2.37865e-01, 4.02425e-02, -1.34963e-03, 2.19717e-05, -1.73520e-07, 5.31624e-10]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=0, yMax=50)
        return
    
    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    fitCfg = {} 
    fitCfg['func_name'] = "ewk_para"
    fitCfg['func_parms_vals'] = [40, 15, 5, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1]
    fitCfg['func_parms_nparms'] = [3, 0, 0] # sigma, mean, norm 
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25, sumw2=False)
        #return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", transform=True)

        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{2} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{3} (GeV)", transform=True)


        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{1} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{2} (GeV)", transform=True)

        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{3} (GeV)", transform=True)


        rls.addParam(jsOut, "p6", "[0]", [0.03])
        rls.addParam(jsOut, "p7", "[0]", [0.35])

     
        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return
        

    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        fitCfg = jsIn
        fitCfg['func_name'] = "ewk_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, sumw2=False)
        with open("%s/results_refit.json" % baseDir, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        
        
        with open("%s/results_refit.json" % baseDir) as f: jsIn = json.load(f)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin, statUnc=False, yRatio=1.25)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")  
          
   
def ttbar_perp():

    procs = ["TTbar"]
    tag = "ttbar_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    recoilMin, recoilMax = -150, 150
    rebin = 2
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    
    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    if True:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-2, yMax=2)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "ttbar_perp"
    fitCfg['func_parms_vals'] = [120, 60, 40]
    fitCfg['func_parms_cfg'] = [1, 1, 1]
    fitCfg['func_parms_nparms'] = [3, 0, 0] # sigma, mean, norm 
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if True:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25)
        #return
        
    outDir_param_v0 = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "cpol2", [60, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "cpol2", [60, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, cutOffMax=80, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "cpol2", [40, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, cutOffMax=80, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")
        
        rls.addParam(jsOut, "p3", "[0]", [0.3])
        rls.addParam(jsOut, "p4", "[0]", [0.15])
     
        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        #return
        
    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "ttbar_perp_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, statUnc=False, yRatio=1.25)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
        #rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        #rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
   
def ttbar_para():

    procs = ["TTbar"]
    tag = "ttbar_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    rebin = 2 
    recoilMin, recoilMax = -150, 150
    bhist = readProc(groups, "recoil_corr_xy_para_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]
    # 0.5 GeV needed for the mean/yields
    if True:
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.99723e+00, -1.47337e-01, 1.61576e-02, -1.42718e-04, -9.91942e-07, 1.47306e-08, -3.50397e-11]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=50, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25)
        return

    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    fitCfg = {} 
    fitCfg['func_name'] = "ttbar_para"
    fitCfg['func_parms_vals'] = [120, 60, 40, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1]
    fitCfg['func_parms_nparms'] = [3, 3, 0] # sigma, mean, norm 
    
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin, yRatio=1.25)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}

        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)", transform=True, fitOpts="NS W")

        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yTitle = "#sigma_{2} (GeV)", transform=True, fitOpts="NS W")
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=200, xMax=200, yTitle = "#sigma_{3} (GeV)", transform=True, fitOpts="NS W")

        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-100, yTitle = "#mu_{1} (GeV)", transform=True, fitOpts="NS W")
        
        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-100, yTitle = "#mu_{2} (GeV)", transform=True, fitOpts="NS W")     

        fitF, params, cParams = "cpol4", [9.90024e-02, 9.34622e-01, 0, 0, 0], [False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-100, yTitle = "#mu_{3} (GeV)", transform=True, fitOpts="NS W")     

        rls.addParam(jsOut, "p6", "[0]", [0.3])
        rls.addParam(jsOut, "p7", "[0]", [0.2])
     
        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return
        

    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "ttbar_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
      
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, statUnc=False, rebin=rebin, yRatio=1.25)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=200, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=200, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")

 

def zmumu_perp():

    procs = ["Zmumu"]
    tag = "zmumu_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "DY #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    recoilMin, recoilMax = -150, 150
    
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)
    
    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-2, yMax=2)
        return


    fitCfg = {} 
    fitCfg['func_name'] = "dy_perp"
    fitCfg['func_parms_vals'] = [30, 20, 3, 5, 10]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    
    fitCfg['func_parms_vals'] = [30, 20, 15, 10, 5]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        #return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        jsIn = utils.loadJSON("%s/results.json" % outDir_fits_v0)
        
        #fitF, params, cParams = "cpol4", [7.64598e+00, 1.05385e-01, 9.46922e-04, -9.91530e-06, 2.21708e-08], [False, False, False, False, False]
        fitF, params, cParams = "linear", [1, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", cutOffMin=30)

        fitF, params, cParams = "cpol5", [18, 1, 0, 0, 0, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [10, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [6, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [4, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{5} (GeV)", transform=True)
       
        # add norm parameters
        rls.addParam(jsOut, "p5", "[0]", [0.0005])
        rls.addParam(jsOut, "p6", "[0]", [0.05])
        rls.addParam(jsOut, "p7", "[0]", [0.40])
        rls.addParam(jsOut, "p8", "[0]", [0.45])

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return


    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        #jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        #fitCfg = jsIn
        #fitCfg['func_name'] = "dy_perp_cond"
        #jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        #utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
      
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        #rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=150, yTitle = "#sigma_{1} (GeV)", transform=False)
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)", transform=True)
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)", transform=True)
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)", transform=True)
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{5} (GeV)", transform=True)
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", transform=False)
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", transform=False)
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", transform=False)
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", transform=False)
            
    # do a refit absorbing correlations
    if False:

        outDir_refit = "%s/params_refit_v1" % baseDir
        outDir_refit_fits = "%s/fits_refit_v1" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        #rls.printVars(fitCfg)
        #quit()
        fitCfg['func_name'] = "dy_perp_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-1
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-1
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-1
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-1
        
        fitCfg['p0']['p1'] = fitCfg['p0']['p2']
        fitCfg['p1']['p1'] = fitCfg['p1']['p2']
        fitCfg['p2']['p1'] = fitCfg['p2']['p2']
        fitCfg['p3']['p1'] = fitCfg['p3']['p2']
        
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        jsOut['func_name'] = "dy_perp_cond"
        jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
        jsOut['p1']['nParams'] = jsOut['p1']['nParams']+1
        jsOut['p2']['nParams'] = jsOut['p2']['nParams']+1
        jsOut['p3']['nParams'] = jsOut['p3']['nParams']+1

        jsOut['p0']['p2'] = jsOut['p0']['p1']
        jsOut['p1']['p2'] = jsOut['p1']['p1']
        jsOut['p2']['p2'] = jsOut['p2']['p1']
        jsOut['p3']['p2'] = jsOut['p3']['p1']
        
        jsOut['p0']['p1'] = jsOut['p0']['p0']*(2.5287E+02)
        jsOut['p1']['p1'] = jsOut['p1']['p0']*(-5.8640E+01)
        jsOut['p2']['p1'] = jsOut['p2']['p0']*(3.2202E+00)
        jsOut['p3']['p1'] = jsOut['p3']['p0']*(5.2319E+01)

        for i in range(jsOut['nStatVars']):
        
            jsOut['stat%d_p'%i]['p0']['p2'] = jsOut['stat%d_p'%i]['p0']['p1']
            jsOut['stat%d_p'%i]['p1']['p2'] = jsOut['stat%d_p'%i]['p1']['p1']
            jsOut['stat%d_p'%i]['p2']['p2'] = jsOut['stat%d_p'%i]['p2']['p1']
            jsOut['stat%d_p'%i]['p3']['p2'] = jsOut['stat%d_p'%i]['p3']['p1']
            
            jsOut['stat%d_p'%i]['p0']['p1'] = jsOut['stat%d_p'%i]['p0']['p0']*(2.5287E+02)
            jsOut['stat%d_p'%i]['p1']['p1'] = jsOut['stat%d_p'%i]['p1']['p0']*(-5.8640E+01)
            jsOut['stat%d_p'%i]['p2']['p1'] = jsOut['stat%d_p'%i]['p2']['p0']*(3.2202E+00)
            jsOut['stat%d_p'%i]['p3']['p1'] = jsOut['stat%d_p'%i]['p3']['p0']*(5.2319E+01)
            
            jsOut['stat%d_m'%i]['p0']['p2'] = jsOut['stat%d_m'%i]['p0']['p1']
            jsOut['stat%d_m'%i]['p1']['p2'] = jsOut['stat%d_m'%i]['p1']['p1']
            jsOut['stat%d_m'%i]['p2']['p2'] = jsOut['stat%d_m'%i]['p2']['p1']
            jsOut['stat%d_m'%i]['p3']['p2'] = jsOut['stat%d_m'%i]['p3']['p1']
            
            jsOut['stat%d_m'%i]['p0']['p1'] = jsOut['stat%d_m'%i]['p0']['p0']*(2.5287E+02)
            jsOut['stat%d_m'%i]['p1']['p1'] = jsOut['stat%d_m'%i]['p1']['p0']*(-5.8640E+01)
            jsOut['stat%d_m'%i]['p2']['p1'] = jsOut['stat%d_m'%i]['p2']['p0']*(3.2202E+00)
            jsOut['stat%d_m'%i]['p3']['p1'] = jsOut['stat%d_m'%i]['p3']['p0']*(5.2319E+01)
            

        utils.writeJSON("%s/results_refit_v1.json" % baseDir, jsOut)
  
        
        jsIn = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")

             
   
    # export and uncertainties
    if True:
    
        exportCfg = {}
        exportCfg["mean"] = [0, 0, 0, 0, 0]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3", "p4"]
        exportCfg["norm"] = ["p5", "p6", "p7", "p8"]
        exportCfg["transform"] = ["p1", "p2", "p3", "p4"]
        
        rls.export(exportCfg, "%s/recoil_zmumu_perp.json" % outCfgDir, "%s/results_refit.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        outDir_unc = "%s/uncertainties/" % baseDir
        utils.mkdir(outDir_unc, False)
        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)


def zmumu_para():

    procs = ["Zmumu"]
    tag = "zmumu_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "DY #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    
    bhist = readProc(groups, "recoil_corr_xy_para_qTbinned", procs)
    recoilMin, recoilMax = -150, 150
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pw_poly8_poly1"
        fitCfg['parms_init'] = [2.98078e+01, -2.50056e-02, 5.03451e-01, -2.86102e-02, 7.96316e-04, -8.18177e-07, -2.47077e-07, -2.73957e-09,  1.20691e-10]
        fitCfg['parms_init'] = [3.75252e+01,  -5.38110e-02, 6.33144e-01, -3.52892e-02, 7.99257e-04, 1.14259e-06, -2.67102e-07, -2.08546e-09, 1.41637e-10, -9.49462e-13]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=25)
        return
        
          

    fitCfg = {} 
    fitCfg['func_name'] = "dy_para"
    fitCfg['func_parms_vals'] = [30, 20, 15, 10, 5, 0, 0, 0, 0, 0.45] 
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    
    fitCfg['func_parms_vals'] = [30, 20, 15, 10, 5, 0, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    
    funcJs = utils.loadJSON("%s/results_refit.json" % baseDir.replace("para", "perp"))
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, funcJs=funcJs)
        #return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        excludeRegions=[[20, 35]]
        excludeRegions = []
        fitF, params, cParams = "linear", [50, 1], [False, False] 
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", cutOffMin=40)

        fitF, params, cParams = "cpol5", [20, 0, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[[0, 1]], fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)", transform=True, cutOffMin=20)
        
        fitF, params, cParams = "cpol5", [12, 0, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [7, 2.59440e-02, 7.84935e-03, -1.98100e-04, 2.00988e-06, -7.35763e-09], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [4, 3.04269e-02, 5.73096e-03, -9.34397e-05, 5.75346e-07, -1.24140e-09], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{5} (GeV)", transform=True)
        
       
        fitF, params, cParams = "cpol5", [1, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)", transform=True, cutOffMin=-5, cutOffMax=10)
        
        fitF, params, cParams = "cpol5", [1, 5.31943e-01, -1.50493e-02, 1.77988e-04, -8.64181e-07, 1.45325e-09], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)", transform=True, cutOffMin=0,  cutOffMax=10)
        
        fitF, params, cParams = "cpol5", [1, 5.31943e-01, -1.50493e-02, 1.77988e-04, -8.64181e-07, 1.45325e-09], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[[20,40]], fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)", transform=True, cutOffMin=0)

        fitF, params, cParams = "cpol5", [-4.04853e-01, 8.51386e-01, -3.64734e-02, 7.97315e-04, -7.86549e-06, 2.78457e-08], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p8", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[[20,40]], fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{4} (GeV)", transform=True, cutOffMin=-2,  cutOffMax=8)

        # add norm parameters
        rls.addParam(jsOut, "p9", "[0]", [0.00017167274118001641]) 
        rls.addParam(jsOut, "p10", "[0]", [0.0186479272549545])  
        rls.addParam(jsOut, "p11", "[0]", [0.26194731500718077])  
        rls.addParam(jsOut, "p12", "[0]", [0.5487736947147043])  

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return



    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        #jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        fitCfg = jsIn
        fitCfg['func_name'] = "dy_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, statUnc=False)

        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)", transform=False)
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)", transform=True)
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)", transform=True)
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)", transform=True)
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{5} (GeV)", transform=True)
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{1} (GeV)", transform=True)
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{2} (GeV)", transform=True)
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{3} (GeV)", transform=True)
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{4} (GeV)", transform=True)
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", transform=False)
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", transform=False)
        rls.plotParameter("p11", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", transform=False)
        rls.plotParameter("p12", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", transform=False)
   
 

    # do a refit absorbing correlations
    if False:

        outDir_refit = "%s/params_refit_v1" % baseDir
        outDir_refit_fits = "%s/fits_refit_v1" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.printVars(fitCfg)
        #quit()
        fitCfg['func_name'] = "dy_para_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-2
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        jsOut['func_name'] = "dy_para_cond"
        jsOut['p0']['nParams'] = jsOut['p0']['nParams']+2
        jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
        jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
        jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

        jsOut['p0']['p1'] = jsOut['p0']['p0']*3.6786E-02
        jsOut['p0']['p2'] = jsOut['p0']['p0']*3.3376E+00
        jsOut['p1']['p1'] = jsOut['p1']['p0']*2.8404E+01
        jsOut['p1']['p2'] = jsOut['p9']['p0']*5.6922E+01
        jsOut['p2']['p1'] = jsOut['p2']['p0']*1.2515E+01
        jsOut['p2']['p2'] = jsOut['p9']['p0']*8.6612E+01
        jsOut['p3']['p1'] = jsOut['p3']['p0']*5.2845E+02
        jsOut['p3']['p2'] = jsOut['p9']['p0']*1.1960E+02

        for i in range(jsOut['nStatVars']):
        
            jsOut['stat%d_p'%i]['p0']['p1'] = jsOut['stat%d_p'%i]['p0']['p0']*3.6786E-02
            jsOut['stat%d_p'%i]['p0']['p2'] = jsOut['stat%d_p'%i]['p0']['p0']*3.3376E+00
            jsOut['stat%d_p'%i]['p1']['p1'] = jsOut['stat%d_p'%i]['p1']['p0']*2.8404E+01
            jsOut['stat%d_p'%i]['p1']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*5.6922E+01
            jsOut['stat%d_p'%i]['p2']['p1'] = jsOut['stat%d_p'%i]['p2']['p0']*1.2515E+01
            jsOut['stat%d_p'%i]['p2']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*8.6612E+01
            jsOut['stat%d_p'%i]['p3']['p1'] = jsOut['stat%d_p'%i]['p3']['p0']*5.2845E+02
            jsOut['stat%d_p'%i]['p3']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*1.1960E+02
            
            jsOut['stat%d_m'%i]['p0']['p1'] = jsOut['stat%d_m'%i]['p0']['p0']*3.6786E-02
            jsOut['stat%d_m'%i]['p0']['p2'] = jsOut['stat%d_m'%i]['p0']['p0']*3.3376E+00
            jsOut['stat%d_m'%i]['p1']['p1'] = jsOut['stat%d_m'%i]['p1']['p0']*2.8404E+01
            jsOut['stat%d_m'%i]['p1']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*5.6922E+01
            jsOut['stat%d_m'%i]['p2']['p1'] = jsOut['stat%d_m'%i]['p2']['p0']*1.2515E+01
            jsOut['stat%d_m'%i]['p2']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*8.6612E+01
            jsOut['stat%d_m'%i]['p3']['p1'] = jsOut['stat%d_m'%i]['p3']['p0']*5.2845E+02
            jsOut['stat%d_m'%i]['p3']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*1.1960E+02
        
        utils.writeJSON("%s/results_refit_v1.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=150, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")

             
  
 
    # export and uncertainties
    if True:
    
        exportCfg = {}
        exportCfg["mean"] = ["p5", "p5", "p6", "p7", "p8"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3", "p4"]
        exportCfg["norm"] = ["p9", "p10", "p11", "p12"]
        exportCfg["transform"] = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]
        
        rls.export(exportCfg, "%s/recoil_zmumu_para.json" % outCfgDir, "%s/results_refit.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        #outDir_unc = "%s/uncertainties/" % baseDir
        #utils.mkdir(outDir_unc, False)
        #fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        #yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        #rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)


def zmumu_sumet():

    procs = ["Zmumu"]
    tag = "zmumu_sumet"
    comp = "sumet"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "DY #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    
    bhist = readProc(groups, "RawMET_sumEt_sqrt_qTbinned", procs)
    recoilMin, recoilMax = 0, 100
    
    binning_qT = list(utils.drange(0, 20, 1)) + [20]
    print(bhist)

    fitCfg = {} 
    fitCfg['func_name'] = "sumet"
    fitCfg['func_parms_vals'] = [11, 5, 7, 10]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 1]  # 0=float, 1=propagate, 2=TF1
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if True:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if True:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [4.54931e+00, 1, -4.25497e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [4.66409e-01, 1, 3.44635e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [5.18556e-01, 1, 15], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)")
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)")
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{5} (GeV)")
        
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        fitF, params, cParams = "quadratic", [1, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=20, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        fitF, params, cParams = "quadratic", [1, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=20, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        fitF, params, cParams = "quadratic", [1, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=20, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
        #fitF, params, cParams = "linear", [0, 1], [False, False]
        fitF, params, cParams = "quadratic", [1, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p8", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=20, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{4} (GeV)")
        
        #fitF, params, cParams = "linear", [0, 1], [False, False]
        fitF, params, cParams = "quadratic", [1, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p9", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=20, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{4} (GeV)")
        
        
    
        # add norm parameters
        rls.addParam(jsOut, "p10", "[0]", [0.1])
        rls.addParam(jsOut, "p11", "[0]", [0.15])
        rls.addParam(jsOut, "p12", "[0]", [0.5])
        rls.addParam(jsOut, "p13", "[0]", [0.24])

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return



    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "dy_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        quit()
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=150, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
 

    # do a refit absorbing correlations
    if True:

        outDir_refit = "%s/params_refit_v1" % baseDir
        outDir_refit_fits = "%s/fits_refit_v1" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.printVars(fitCfg)
        #quit()
        fitCfg['func_name'] = "dy_para_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-2
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        jsOut['func_name'] = "dy_para_cond"
        jsOut['p0']['nParams'] = jsOut['p0']['nParams']+2
        jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
        jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
        jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

        jsOut['p0']['p1'] = jsOut['p0']['p0']*3.6786E-02
        jsOut['p0']['p2'] = jsOut['p0']['p0']*3.3376E+00
        jsOut['p1']['p1'] = jsOut['p1']['p0']*2.8404E+01
        jsOut['p1']['p2'] = jsOut['p9']['p0']*5.6922E+01
        jsOut['p2']['p1'] = jsOut['p2']['p0']*1.2515E+01
        jsOut['p2']['p2'] = jsOut['p9']['p0']*8.6612E+01
        jsOut['p3']['p1'] = jsOut['p3']['p0']*5.2845E+02
        jsOut['p3']['p2'] = jsOut['p9']['p0']*1.1960E+02

        for i in range(jsOut['nStatVars']):
        
            jsOut['stat%d_p'%i]['p0']['p1'] = jsOut['stat%d_p'%i]['p0']['p0']*3.6786E-02
            jsOut['stat%d_p'%i]['p0']['p2'] = jsOut['stat%d_p'%i]['p0']['p0']*3.3376E+00
            jsOut['stat%d_p'%i]['p1']['p1'] = jsOut['stat%d_p'%i]['p1']['p0']*2.8404E+01
            jsOut['stat%d_p'%i]['p1']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*5.6922E+01
            jsOut['stat%d_p'%i]['p2']['p1'] = jsOut['stat%d_p'%i]['p2']['p0']*1.2515E+01
            jsOut['stat%d_p'%i]['p2']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*8.6612E+01
            jsOut['stat%d_p'%i]['p3']['p1'] = jsOut['stat%d_p'%i]['p3']['p0']*5.2845E+02
            jsOut['stat%d_p'%i]['p3']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*1.1960E+02
            
            jsOut['stat%d_m'%i]['p0']['p1'] = jsOut['stat%d_m'%i]['p0']['p0']*3.6786E-02
            jsOut['stat%d_m'%i]['p0']['p2'] = jsOut['stat%d_m'%i]['p0']['p0']*3.3376E+00
            jsOut['stat%d_m'%i]['p1']['p1'] = jsOut['stat%d_m'%i]['p1']['p0']*2.8404E+01
            jsOut['stat%d_m'%i]['p1']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*5.6922E+01
            jsOut['stat%d_m'%i]['p2']['p1'] = jsOut['stat%d_m'%i]['p2']['p0']*1.2515E+01
            jsOut['stat%d_m'%i]['p2']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*8.6612E+01
            jsOut['stat%d_m'%i]['p3']['p1'] = jsOut['stat%d_m'%i]['p3']['p0']*5.2845E+02
            jsOut['stat%d_m'%i]['p3']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*1.1960E+02
        
        utils.writeJSON("%s/results_refit_v1.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=150, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")

             
  
 
    # export and uncertainties
    if True:
    
        exportCfg = {}
        exportCfg["mean"] = ["p4", "p5", "p6", "p7"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p8", "p9", "p10"]
        
        rls.export(exportCfg, "%s/recoil_zmumu_para.json" % outCfgDir, "%s/results_refit.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        outDir_unc = "%s/uncertainties/" % baseDir
        utils.mkdir(outDir_unc, False)
        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)


    
def singlemuon_perp():

    ## TODO
    ## do all steps including v0
    
    procs = ["Data"]
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "Data #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    
    binning_qT = list(utils.drange(0, 100, 1)) + [100]
    
    recoilMin, recoilMax = -150, 150
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-5, yMax=5)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "data_perp"
    fitCfg['func_parms_vals'] = [30, 20, 15, 10, 5]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    fitCfg['func_parms_nparms'] = [5, 0, 0] # sigma, mean, norm 
    
    
    bkgCfg = {}
    bkgCfg['procs'] = ["ttbar", "ewk"]
    bkgCfg['parms'] = [utils.loadJSON("%s/ttbar_perp/results_refit.json"%outDir), utils.loadJSON("%s/ewk_perp/results_refit.json"%outDir)]
    bkgCfg['yields'] = [utils.loadJSON("%s/ttbar_perp/yields.json"%outDir), utils.loadJSON("%s/ewk_perp/yields.json"%outDir)]
    bkgCfg['norms'] = [1.0, 1.0]
    bkgCfg['data_yields'] = utils.loadJSON("%s/yields.json"%baseDir)
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        jsIn = utils.loadJSON("%s/results.json" % outDir_fits_v0)
        
        #fitF, params, cParams = "cpol4", [7.64598e+00, 1.05385e-01, 9.46922e-04, -9.91530e-06, 2.21708e-08], [False, False, False, False, False]
        fitF, params, cParams = "linear", [1, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", cutOffMin=30)

        fitF, params, cParams = "cpol5", [17, 0, 0, 0, 0, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [10, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [6, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)", transform=True)
        
        fitF, params, cParams = "cpol5", [4, 1, 0, 0, 0, 0], [False, False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{5} (GeV)", transform=True)
       
        # add norm parameters
        rls.addParam(jsOut, "p5", "[0]", [0.0005])
        rls.addParam(jsOut, "p6", "[0]", [0.05])
        rls.addParam(jsOut, "p7", "[0]", [0.40])
        rls.addParam(jsOut, "p8", "[0]", [0.45])

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return

    
    
    
    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        #fitCfg = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = utils.loadJSON("%s/zmumu_perp/results_refit.json" % outDir)
        
        fitCfg['func_name'] = "data_perp_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
       
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=150, yTitle = "#sigma_{1} (GeV)", transform=False)
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)", transform=True)
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)", transform=True)
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)", transform=True)
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{5} (GeV)", transform=True)
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", transform=False)
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", transform=False)
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", transform=False)
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", transform=False)
   
    # do a refit absorbing correlations
    if False:

        outDir_refit = "%s/params_refit_v1" % baseDir
        outDir_refit_fits = "%s/fits_refit_v1" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.printVars(fitCfg)
        #quit()
        fitCfg['func_name'] = "data_perp_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-1
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        fitCfg['p0']['p1'] = fitCfg['p0']['p2']
        #fitCfg['p1']['p1'] = fitCfg['p1']['p2']
        #fitCfg['p2']['p1'] = fitCfg['p2']['p2']
        #fitCfg['p3']['p1'] = fitCfg['p3']['p2']

        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        jsOut['func_name'] = "data_perp_cond"
        jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
        jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
        jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
        jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

        jsOut['p0']['p2'] = jsOut['p0']['p1']
        jsOut['p0']['p1'] = jsOut['p0']['p0']*1.3862E+02
        jsOut['p1']['p1'] = jsOut['p1']['p0']*(-1.8742E+01)
        jsOut['p1']['p2'] = jsOut['p7']['p0']*5.2746E+02
        jsOut['p2']['p1'] = jsOut['p2']['p0']*8.3859E+00
        jsOut['p2']['p2'] = jsOut['p8']['p0']*2.8648E+01
        jsOut['p3']['p1'] = jsOut['p3']['p0']*3.6196E+01
        jsOut['p3']['p2'] = jsOut['p8']['p0']*4.0329E+01
        

        for i in range(jsOut['nStatVars']):
        
            jsOut['stat%d_p'%i]['p0']['p2'] = jsOut['stat%d_p'%i]['p0']['p1']
            jsOut['stat%d_p'%i]['p0']['p1'] = jsOut['stat%d_p'%i]['p0']['p0']*1.3862E+02
            jsOut['stat%d_p'%i]['p1']['p1'] = jsOut['stat%d_p'%i]['p1']['p0']*(-1.8742E+01)
            jsOut['stat%d_p'%i]['p1']['p2'] = jsOut['stat%d_p'%i]['p7']['p0']*5.2746E+02
            jsOut['stat%d_p'%i]['p2']['p1'] = jsOut['stat%d_p'%i]['p2']['p0']*8.3859E+00
            jsOut['stat%d_p'%i]['p2']['p2'] = jsOut['stat%d_p'%i]['p8']['p0']*2.8648E+01
            jsOut['stat%d_p'%i]['p3']['p1'] = jsOut['stat%d_p'%i]['p3']['p0']*3.6196E+01
            jsOut['stat%d_p'%i]['p3']['p2'] = jsOut['stat%d_p'%i]['p8']['p0']*4.0329E+01
            
            jsOut['stat%d_m'%i]['p0']['p2'] = jsOut['stat%d_m'%i]['p0']['p1']
            jsOut['stat%d_m'%i]['p0']['p1'] = jsOut['stat%d_m'%i]['p0']['p0']*1.3862E+02
            jsOut['stat%d_m'%i]['p1']['p1'] = jsOut['stat%d_m'%i]['p1']['p0']*(-1.8742E+01)
            jsOut['stat%d_m'%i]['p1']['p2'] = jsOut['stat%d_m'%i]['p7']['p0']*5.2746E+02
            jsOut['stat%d_m'%i]['p2']['p1'] = jsOut['stat%d_m'%i]['p2']['p0']*8.3859E+00
            jsOut['stat%d_m'%i]['p2']['p2'] = jsOut['stat%d_m'%i]['p8']['p0']*2.8648E+01
            jsOut['stat%d_m'%i]['p3']['p1'] = jsOut['stat%d_m'%i]['p3']['p0']*3.6196E+01
            jsOut['stat%d_m'%i]['p3']['p2'] = jsOut['stat%d_m'%i]['p8']['p0']*4.0329E+01
        
            

        utils.writeJSON("%s/results_refit_v1.json" % baseDir, jsOut)
  
        jsIn = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")


    exportCfg = {}
    exportCfg["mean"] = [0, 0, 0, 0, 0]
    exportCfg["sigma"] = ["p0", "p1", "p2", "p3", "p4"]
    exportCfg["norm"] = ["p5", "p6", "p7", "p8"]
    exportCfg["transform"] = ["p1", "p2", "p3", "p4"]
 
    if False: # bkg syst variations
        
        fitCfg = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        fitCfg['func_name'] = "data_perp_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-1
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        fitCfg['p0']['p1'] = fitCfg['p0']['p2']

        if True:
            bkgCfg['norms'] = [1.2, 1.0]
            outDir_refit_fits = "%s/bkg/ttbarUp/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_perp_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*1.3862E+02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*(-1.8742E+01)
            jsOut['p1']['p2'] = jsOut['p7']['p0']*5.2746E+02
            jsOut['p2']['p1'] = jsOut['p2']['p0']*8.3859E+00
            jsOut['p2']['p2'] = jsOut['p8']['p0']*2.8648E+01
            jsOut['p3']['p1'] = jsOut['p3']['p0']*3.6196E+01
            jsOut['p3']['p2'] = jsOut['p8']['p0']*4.0329E+01
        
            utils.writeJSON("%s/results_refit_ttbarUp.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_perp_ttbarUp.json" % outCfgDir, "%s/results_refit_ttbarUp.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [0.8, 1.0]
            outDir_refit_fits = "%s/bkg/ttbarDown/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_perp_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*1.3862E+02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*(-1.8742E+01)
            jsOut['p1']['p2'] = jsOut['p7']['p0']*5.2746E+02
            jsOut['p2']['p1'] = jsOut['p2']['p0']*8.3859E+00
            jsOut['p2']['p2'] = jsOut['p8']['p0']*2.8648E+01
            jsOut['p3']['p1'] = jsOut['p3']['p0']*3.6196E+01
            jsOut['p3']['p2'] = jsOut['p8']['p0']*4.0329E+01
            
            utils.writeJSON("%s/results_refit_ttbarDown.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_perp_ttbarDown.json" % outCfgDir, "%s/results_refit_ttbarDown.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [1.0, 1.2]
            outDir_refit_fits = "%s/bkg/ewkUp/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_perp_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*1.3862E+02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*(-1.8742E+01)
            jsOut['p1']['p2'] = jsOut['p7']['p0']*5.2746E+02
            jsOut['p2']['p1'] = jsOut['p2']['p0']*8.3859E+00
            jsOut['p2']['p2'] = jsOut['p8']['p0']*2.8648E+01
            jsOut['p3']['p1'] = jsOut['p3']['p0']*3.6196E+01
            jsOut['p3']['p2'] = jsOut['p8']['p0']*4.0329E+01
            
            utils.writeJSON("%s/results_refit_ewkUp.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_perp_ewkUp.json" % outCfgDir, "%s/results_refit_ewkUp.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [1.0, 0.8]
            outDir_refit_fits = "%s/bkg/ewkDown/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_perp_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*1.3862E+02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*(-1.8742E+01)
            jsOut['p1']['p2'] = jsOut['p7']['p0']*5.2746E+02
            jsOut['p2']['p1'] = jsOut['p2']['p0']*8.3859E+00
            jsOut['p2']['p2'] = jsOut['p8']['p0']*2.8648E+01
            jsOut['p3']['p1'] = jsOut['p3']['p0']*3.6196E+01
            jsOut['p3']['p2'] = jsOut['p8']['p0']*4.0329E+01
            
            utils.writeJSON("%s/results_refit_ewkDown.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_perp_ewkDown.json" % outCfgDir, "%s/results_refit_ewkDown.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
             

    # export and uncertainties
    if True:
    
        rls.export(exportCfg, "%s/recoil_data_perp.json" % outCfgDir, "%s/results_refit.json" % baseDir) ## refit
        #rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)

        #outDir_unc = "%s/uncertainties/" % baseDir
        #utils.mkdir(outDir_unc, False)
        #fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        #yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        #rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)

        #up = utils.loadJSON("%s/results_refit_ttbarUp.json" % baseDir)
        #dw = utils.loadJSON("%s/results_refit_ttbarDown.json" % baseDir)
        #rls.plotSystUncertainty(fitCfg, up, dw, "syst_ttbar", "TTbar up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
        
        #up = utils.loadJSON("%s/results_refit_ewkUp.json" % baseDir)
        #dw = utils.loadJSON("%s/results_refit_ewkDown.json" % baseDir)
        #rls.plotSystUncertainty(fitCfg, up, dw, "syst_ewk", "EWK up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
   
  
def singlemuon_para():

    procs = ["Data"]
    tag = "singlemuon_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "Data #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    rebin = 1

    outDir_fits_v0 = "%s/fits_v0" % baseDir
    recoilMin, recoilMax = -150, 150
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]
    
    bhist = readProc(groups, "recoil_corr_xy_para_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pw_poly7_poly1"
        fitCfg['parms_init'] = [8.87682e+01, 5.33854e-02 , 4.43228e-01 ,-1.42881e-02 , 2.77539e-04 ,-1.69169e-06 ,-2.38482e-08 , 4.20338e-10 ,-1.79513e-12]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=50)
        return
        

    fitCfg = {} 
    fitCfg['func_name'] = "data_para"
    fitCfg['func_parms_vals'] = [30, 20, 15, 10, 5, 0, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [1, 1, 1, 1, 1, 1, 1, 1, 1]  # 0=float, 1=propagate, 2=TF1
    fitCfg['func_parms_nparms'] = [5, 3, 0] # sigma, mean, norm 


    bkgCfg = {}
    bkgCfg['procs'] = ["ttbar", "ewk"]
    bkgCfg['parms'] = [utils.loadJSON("%s/ttbar_para/results_refit.json"%outDir), utils.loadJSON("%s/ewk_para/results_refit.json"%outDir)]
    bkgCfg['yields'] = [utils.loadJSON("%s/ttbar_para/yields.json"%outDir), utils.loadJSON("%s/ewk_para/yields.json"%outDir)]
    bkgCfg['norms'] = [1.0, 1.0]
    bkgCfg['data_yields'] = utils.loadJSON("%s/singlemuon_para/yields.json"%outDir) # needed for the background fractions

    funcJs = utils.loadJSON("%s/results_refit.json" % baseDir.replace("para", "perp"))

    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, funcJs=funcJs)
        #return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        excludeRegions=[[15, 35]]
        excludeRegions=[]
        fitF, params, cParams = "linear", [40, 1], [False, False] 
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=40, yMax=200, xMax=100, yTitle = "#sigma_{1} (GeV)", cutOffMin=30)

        fitF, params, cParams = "cpol5", [17, 9.09773e-02, 2.68088e-03, -6.23169e-05, 4.45630e-07, -1.02989e-09], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)", transform=True, cutOffMin=17, cutOffMax=35, fitOpts="NS W")
        
        fitF, params, cParams = "cpol5", [13, 9.09773e-02, 2.68088e-03, -6.23169e-05, 4.45630e-07, -1.02989e-09], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", transform=True, fitOpts="NS W")
        
        fitF, params, cParams = "cpol5", [8, 1, 0, 0, 0, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)", transform=True, fitOpts="NS W")
        
        fitF, params, cParams = "cpol5", [4, 9.09773e-02, 2.68088e-03, -6.23169e-05, 4.45630e-07, -1.02989e-09], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=[], fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{5} (GeV)", transform=True, cutOffMax=15, fitOpts="NS W")

        fitF, params, cParams = "cpol5", [3.07948e-01, 4.18206e-01, -9.89649e-03, 8.87369e-05, -2.46014e-07, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)", transform=True, cutOffMin=0, cutOffMax=100, fitOpts="NS W")
        
        fitF, params, cParams = "cpol5", [3.07948e-01, 4.18206e-01, -9.89649e-03, 8.87369e-05, -2.46014e-07, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)", transform=True, cutOffMin=00, cutOffMax=20, fitOpts="NS W")
        
        fitF, params, cParams = "cpol5", [3.07948e-01, 4.18206e-01, -9.89649e-03, 8.87369e-05, -2.46014e-07, 0], [False, False, False, False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)", transform=True, cutOffMin=0, cutOffMax=10, fitOpts="NS W")
      
        fitF, params, cParams = "cpol5", [3.07948e-01, 4.18206e-01, -9.89649e-03, 8.87369e-05, -2.46014e-07, 0], [False, False, False, False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p8", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, excludeRegions=excludeRegions, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{4} (GeV)", transform=True, cutOffMin=-2, cutOffMax=100, fitOpts="NS W")
        
        # add norm parameters
        rls.addParam(jsOut, "p9", "[0]", [0.0005]) 
        rls.addParam(jsOut, "p10", "[0]", [0.05])  
        rls.addParam(jsOut, "p11", "[0]", [0.50])  
        rls.addParam(jsOut, "p12", "[0]", [0.35])  
        
        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return


    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        fitCfg = utils.loadJSON("%s/results.json" % outDir_param_v0)
        #fitCfg = utils.loadJSON("%s/zmumu_para/results_refit.json" % outDir)
        
        fitCfg['func_name'] = "data_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

         
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
         
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)", transform=False)
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)", transform=True)
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)", transform=True)
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)", transform=True)
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{5} (GeV)", transform=True)
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{1} (GeV)", transform=True)
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{2} (GeV)", transform=True)
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{3} (GeV)", transform=True)
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=40, yTitle = "#mu_{4} (GeV)", transform=True)
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", transform=False)
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", transform=False)
        rls.plotParameter("p11", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", transform=False)
        rls.plotParameter("p12", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", transform=False)


    # do a refit absorbing correlations
    if False:

        outDir_refit = "%s/params_refit_v1" % baseDir
        outDir_refit_fits = "%s/fits_refit_v1" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)

        fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.printVars(fitCfg)
        #quit()
        fitCfg['func_name'] = "data_para_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-1
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        fitCfg['p0']['p1'] = fitCfg['p0']['p2']

        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        jsOut['func_name'] = "data_para_cond"
        jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
        jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
        jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
        jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

        jsOut['p0']['p2'] = jsOut['p0']['p1']
        jsOut['p0']['p1'] = jsOut['p0']['p0']*5.0117E-02
        jsOut['p1']['p1'] = jsOut['p1']['p0']*3.0170E+01
        jsOut['p1']['p2'] = jsOut['p9']['p0']*7.7316E+01
        jsOut['p2']['p1'] = jsOut['p2']['p0']*2.6298E+01
        jsOut['p2']['p2'] = jsOut['p9']['p0']*1.1472E+02
        jsOut['p3']['p1'] = jsOut['p3']['p0']*1.0389E+02
        jsOut['p3']['p2'] = jsOut['p9']['p0']*1.5667E+02
        
        
        


        for i in range(jsOut['nStatVars']):
        
            jsOut['stat%d_p'%i]['p0']['p2'] = jsOut['stat%d_p'%i]['p0']['p1'] 
            jsOut['stat%d_p'%i]['p0']['p1'] = jsOut['stat%d_p'%i]['p0']['p0']*5.0117E-02
            jsOut['stat%d_p'%i]['p1']['p1'] = jsOut['stat%d_p'%i]['p1']['p0']*3.0170E+01
            jsOut['stat%d_p'%i]['p1']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*7.7316E+01
            jsOut['stat%d_p'%i]['p2']['p1'] = jsOut['stat%d_p'%i]['p2']['p0']*2.6298E+01
            jsOut['stat%d_p'%i]['p2']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*1.1472E+02
            jsOut['stat%d_p'%i]['p3']['p1'] = jsOut['stat%d_p'%i]['p3']['p0']*1.0389E+02
            jsOut['stat%d_p'%i]['p3']['p2'] = jsOut['stat%d_p'%i]['p9']['p0']*1.5667E+02
            
            jsOut['stat%d_m'%i]['p0']['p2'] = jsOut['stat%d_m'%i]['p0']['p1'] 
            jsOut['stat%d_m'%i]['p0']['p1'] = jsOut['stat%d_m'%i]['p0']['p0']*5.0117E-02
            jsOut['stat%d_m'%i]['p1']['p1'] = jsOut['stat%d_m'%i]['p1']['p0']*3.0170E+01
            jsOut['stat%d_m'%i]['p1']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*7.7316E+01
            jsOut['stat%d_m'%i]['p2']['p1'] = jsOut['stat%d_m'%i]['p2']['p0']*2.6298E+01
            jsOut['stat%d_m'%i]['p2']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*1.1472E+02
            jsOut['stat%d_m'%i]['p3']['p1'] = jsOut['stat%d_m'%i]['p3']['p0']*1.0389E+02
            jsOut['stat%d_m'%i]['p3']['p2'] = jsOut['stat%d_m'%i]['p9']['p0']*1.5667E+02
        
            

        utils.writeJSON("%s/results_refit_v1.json" % baseDir, jsOut)
  
        jsIn = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-50, yMax=50, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p9", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p10", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")


    exportCfg = {}
    exportCfg["mean"] = ["p5", "p5","p6", "p7", "p8"]
    exportCfg["sigma"] = ["p0", "p1", "p2", "p3", "p4"]
    exportCfg["norm"] = ["p9", "p10", "p11", "p12"]
    exportCfg["transform"] = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]
 
    if False: # bkg syst variations
        
        fitCfg = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        fitCfg['func_name'] = "data_para_cond_v1"
        fitCfg['p0']['nParams'] = fitCfg['p0']['nParams']-1
        fitCfg['p1']['nParams'] = fitCfg['p1']['nParams']-2
        fitCfg['p2']['nParams'] = fitCfg['p2']['nParams']-2
        fitCfg['p3']['nParams'] = fitCfg['p3']['nParams']-2
        fitCfg['p0']['p1'] = fitCfg['p0']['p2']

        if True:
            bkgCfg['norms'] = [1.2, 1.0]
            outDir_refit_fits = "%s/bkg/ttbarUp/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_para_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*5.0117E-02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*3.0170E+01
            jsOut['p1']['p2'] = jsOut['p9']['p0']*7.7316E+01
            jsOut['p2']['p1'] = jsOut['p2']['p0']*2.6298E+01
            jsOut['p2']['p2'] = jsOut['p9']['p0']*1.1472E+02
            jsOut['p3']['p1'] = jsOut['p3']['p0']*1.0389E+02
            jsOut['p3']['p2'] = jsOut['p9']['p0']*1.5667E+02
            
            utils.writeJSON("%s/results_refit_ttbarUp.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_para_ttbarUp.json" % outCfgDir, "%s/results_refit_ttbarUp.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [0.8, 1.0]
            outDir_refit_fits = "%s/bkg/ttbarDown/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_para_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*5.0117E-02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*3.0170E+01
            jsOut['p1']['p2'] = jsOut['p9']['p0']*7.7316E+01
            jsOut['p2']['p1'] = jsOut['p2']['p0']*2.6298E+01
            jsOut['p2']['p2'] = jsOut['p9']['p0']*1.1472E+02
            jsOut['p3']['p1'] = jsOut['p3']['p0']*1.0389E+02
            jsOut['p3']['p2'] = jsOut['p9']['p0']*1.5667E+02
            
            utils.writeJSON("%s/results_refit_ttbarDown.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_para_ttbarDown.json" % outCfgDir, "%s/results_refit_ttbarDown.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [1.0, 1.2]
            outDir_refit_fits = "%s/bkg/ewkUp/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_para_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*5.0117E-02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*3.0170E+01
            jsOut['p1']['p2'] = jsOut['p9']['p0']*7.7316E+01
            jsOut['p2']['p1'] = jsOut['p2']['p0']*2.6298E+01
            jsOut['p2']['p2'] = jsOut['p9']['p0']*1.1472E+02
            jsOut['p3']['p1'] = jsOut['p3']['p0']*1.0389E+02
            jsOut['p3']['p2'] = jsOut['p9']['p0']*1.5667E+02
            
            utils.writeJSON("%s/results_refit_ewkUp.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_para_ewkUp.json" % outCfgDir, "%s/results_refit_ewkUp.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
         
        if True:
            bkgCfg['norms'] = [1.0, 0.8]
            outDir_refit_fits = "%s/bkg/ewkDown/fits_refit" % baseDir
            utils.mkdir(outDir_refit_fits, False)
            jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax)
            rls.doFitMultiGauss_plot(bhist, comp, jsOut, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
            
            jsOut['func_name'] = "data_para_cond"
            jsOut['p0']['nParams'] = jsOut['p0']['nParams']+1
            jsOut['p1']['nParams'] = jsOut['p1']['nParams']+2
            jsOut['p2']['nParams'] = jsOut['p2']['nParams']+2
            jsOut['p3']['nParams'] = jsOut['p3']['nParams']+2

            jsOut['p0']['p2'] = jsOut['p0']['p1']
            jsOut['p0']['p1'] = jsOut['p0']['p0']*5.0117E-02
            jsOut['p1']['p1'] = jsOut['p1']['p0']*3.0170E+01
            jsOut['p1']['p2'] = jsOut['p9']['p0']*7.7316E+01
            jsOut['p2']['p1'] = jsOut['p2']['p0']*2.6298E+01
            jsOut['p2']['p2'] = jsOut['p9']['p0']*1.1472E+02
            jsOut['p3']['p1'] = jsOut['p3']['p0']*1.0389E+02
            jsOut['p3']['p2'] = jsOut['p9']['p0']*1.5667E+02
            
            utils.writeJSON("%s/results_refit_ewkDown.json" % baseDir, jsOut)
            rls.export(exportCfg, "%s/recoil_data_para_ewkDown.json" % outCfgDir, "%s/results_refit_ewkDown.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
        
 
    # export and uncertainties
    if True:
    
        rls.export(exportCfg, "%s/recoil_data_para.json" % outCfgDir, "%s/results_refit.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        #outDir_unc = "%s/uncertainties/" % baseDir
        #utils.mkdir(outDir_unc, False)
        #fitCfg = utils.loadJSON("%s/results_refit.json" % baseDir)
        #yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        #rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)

        #up = utils.loadJSON("%s/results_refit_ttbarUp.json" % baseDir)
        #dw = utils.loadJSON("%s/results_refit_ttbarDown.json" % baseDir)
        #rls.plotSystUncertainty(fitCfg, up, dw, "syst_ttbar", "TTbar up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
        
        #up = utils.loadJSON("%s/results_refit_ewkUp.json" % baseDir)
        #dw = utils.loadJSON("%s/results_refit_ewkDown.json" % baseDir)
        #rls.plotSystUncertainty(fitCfg, up, dw, "syst_ewk", "EWK up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
  
  

def readProc(groups, hName, procs):

    label = "%s_%s" % (hName, procs[0]) 
    groups.setHists(hName, "", label=label, procsToRead=procs)
    bhist = groups.groups[procs[0]][label]
    return bhist 



    

    
if __name__ == "__main__":

    met = "DeepMETReso"
    flavor = "mumu"
    rls.__topRight__ = "16.8 fb^{#minus1} (13 TeV)"
    
    groups = datagroups2016("mz_wlike_with_mu_eta_pt_%s.hdf5" % met)
    groups.groups.update({
        "TTbar" : dict(
                members = [groups.datasets[x] for x in ["TTLeptonicPostVFP", "TTSemileptonicPostVFP", "SingleTschanLepDecaysPostVFP", "SingleTtWAntitopPostVFP", "SingleTtchanAntitopPostVFP", "SingleTtchanTopPostVFP", "SingleTtWTopPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
        ),
        "EWK" : dict(
                members = [groups.datasets[x] for x in ["WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "WWTo2L2Nu", "WWTo1L1Nu2QPostVFP", "WZTo3LNuPostVFP", "WZTo2Q2LPostVFP", "WZTo1L1Nu2QPostVFP", "ZZTo2L2NuPostVFP", "ZZTo2Q2LPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
        )
    })
    
    # central
    # "WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "WWTo1L1Nu2QPostVFP",
    
    # ["WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "WWTo2L2Nu", "WWTo1L1Nu2QPostVFP", "WZTo3LNuPostVFP", "WZTo2Q2LPostVFP", "WZTo1L1Nu2QPostVFP", "ZZTo2L2NuPostVFP", "ZZTo2Q2LPostVFP"]

    #   , "WZTo3LNuPostVFP" "WZTo2Q2LPostVFP", "WZTo1L1Nu2QPostVFP",  "ZZTo2L2NuPostVFP", "ZZTo2Q2LPostVFP"
    outCfgDir = f"wremnants/data/recoil/highPU/{flavor}_{met}/"
    outDir = "/eos/user/j/jaeyserm/www/recoil/highPU/%s_%s/" % (flavor, met)
    utils.mkdir(outDir, False)
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]

    import highPU_DeepMETReso_functions as rf
    rls.setFunctionLibs(rf)


    #zmumu_sumet()
    zmumu_para()
    #zmumu_perp()
        
    #ttbar_para()
    #ttbar_perp()
    
    #ewk_para()
    #ewk_perp()

    #singlemuon_para()
    #singlemuon_perp()
 