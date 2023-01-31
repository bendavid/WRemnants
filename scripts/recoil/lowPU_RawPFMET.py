
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
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU

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
    recoilMin, recoilMax = -100, 100
    procNames = groups.getProcNames(to_expand=procs)
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    rebin=2
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-2, yMax=2)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "ewk_perp"
    fitCfg['func_parms_vals'] = [10, 20, 70]
    fitCfg['func_parms_cfg'] = [0, 0, 0]
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-5, yMax=1e3, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [9.90024e-02, 9.34622e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "power", [4.66409e-01, 3.46158e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "linear", [4.66409e-01, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=200, xMax=100, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")

        #rls.addParam(jsOut, "p3", "[0]", [0.75])
        #rls.addParam(jsOut, "p4", "[0]", [0.20])
     
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
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=baseDir, chisq_refit=False, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-5, yMax=1e3, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{3} (GeV)")


def ewk_para_qT():

    procs = ["EWK"]
    tag = "ewk_para_qT"
    comp = "para_qT"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    bhist = readProc(groups, "recoil_corr_xy_para_qT_qTbinned", procs)
    
    recoilMin, recoilMax = -100, 100
    rebin=2
    if False:
        
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.19057e+00, -2.37865e-01, 4.02425e-02, -1.34963e-03, 2.19717e-05, -1.73520e-07, 5.31624e-10]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=0, yMax=50)
        return
    
    fitCfg = {} 
    fitCfg['func_name'] = "ewk_para_qT"
    fitCfg['func_parms_vals'] = [10, 25, 50, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0]
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-5, yMax=1e3, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [1, 1, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "power", [1, 1, 0], [False, True, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "power", [1, 1, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")


        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{1} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{2} (GeV)", fitOpts="NS W", cutOffMin=0)

        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yMin=-50, yTitle = "#mu_{3} (GeV)", fitOpts="NS W", cutOffMin=0)


        rls.addParam(jsOut, "p6", "[0]", [0.33])
        rls.addParam(jsOut, "p7", "[0]", [0.33])

     
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
        fitCfg['func_name'] = "ewk_para_qT_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=baseDir, chisq_refit=False, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        with open("%s/results_refit.json" % baseDir, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        
      
        
        with open("%s/results_refit.json" % baseDir) as f: jsIn = json.load(f)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-5, yMax=1e3, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")  
          
   
def ttbar_perp():

    procs = ["Top"]
    tag = "ttbar_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    recoilMin, recoilMax = -100, 100
    rebin = 2
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-2, yMax=2)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "ttbar_perp"
    fitCfg['func_parms_vals'] = [20, 50, 100]
    fitCfg['func_parms_cfg'] = [0, 0, 0]
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-5, yMax=1e3, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "quadratic", [9.90024e-02, 9.34622e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", fitOpts="NS w")

        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{2} (GeV)", fitOpts="NS w")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=200, xMax=100, yTitle = "#sigma_{3} (GeV)", fitOpts="NS w")

        rls.addParam(jsOut, "p3", "[0]", [0.5])
        rls.addParam(jsOut, "p4", "[0]", [0.4])
     
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
        fitCfg['func_name'] = "ttbar_perp_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=baseDir, chisq_refit=False, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-5, yMax=1e3, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")

    
def ttbar_para_qT():

    procs = ["Top"]
    tag = "ttbar_para_qT"
    comp = "para_qT"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    rebin = 2 
    recoilMin, recoilMax = -100, 100
    bhist = readProc(groups, "recoil_corr_xy_para_qT_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)

    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.99723e+00, -1.47337e-01, 1.61576e-02, -1.42718e-04, -9.91942e-07, 1.47306e-08, -3.50397e-11]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=50, recoilLow=recoilMin, recoilHigh=recoilMax)
        return

    
    fitCfg = {} 
    fitCfg['func_name'] = "ttbar_para_qT"
    fitCfg['func_parms_vals'] = [20, 50, 100, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0]
    
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-5, yMax=1e3, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}

        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "quadratic", [9.90024e-02, 9.34622e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=100, xMax=100, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=200, xMax=100, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")

        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=100, xMax=100, yMin=-100, yTitle = "#mu_{1} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=100, xMax=100, yMin=-100, yTitle = "#mu_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=10, fitMax=100, yMax=100, xMax=100, yMin=-100, yTitle = "#mu_{3} (GeV)", fitOpts="NS W")
        

        rls.addParam(jsOut, "p6", "[0]", [0.65])
        rls.addParam(jsOut, "p7", "[0]", [0.25])
     
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
        fitCfg['func_name'] = "ttbar_para_qT_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=baseDir, chisq_refit=False, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
      
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-5, yMax=1e3, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")

 
def zmumu_perp():

    procs = ["Zmumu"]
    tag = "zmumu_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "DY #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    recoilMin, recoilMax = -100, 100
    doExport = True
    
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-2, yMax=2)
        return


    fitCfg = {} 
    fitCfg['func_name'] = "dy_perp"
    fitCfg['func_parms_vals'] = [15, 5, 8, 10, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        jsIn = utils.loadJSON("%s/results.json" % outDir_fits_v0)
        
        fitF, params, cParams = "power", [9.90024e-02, 9.34622e-01, 1.26746e+01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [4.66409e-01, 3.46158e-01, 3.44635e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [1, 1, 7.5], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", cutOffMin=7.5)
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)", cutOffMax=13)
        
        
        fitF, params, cParams = "linear", [0, 0], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-5, cutOffMax=3, cutOffMin=-3, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "linear", [0, 0], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-5, cutOffMax=3, cutOffMin=-3, yTitle = "#mu_{2} (GeV)")
         
        # add norm parameters
        rls.addParam(jsOut, "p6", "[0]", [0.1])
        rls.addParam(jsOut, "p7", "[0]", [0.20])
        rls.addParam(jsOut, "p8", "[0]", [0.25])

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return


    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "dy_perp_cond"
        #jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, chisq_refit=False, outDir=baseDir, recoilLow=recoilMin, recoilHigh=recoilMax)
        #utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
      
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        #rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
    if doExport:
        exportCfg = {}
        exportCfg["mean"] = ["p4", "p5", "p4", "p5"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p6", "p7", "p8"]
        rls.export(exportCfg, "%s/recoil_zmumu_perp.json" % outCfgDir, "%s/results_refit.json" % baseDir)
         
    
    if False:
    
        bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned_qTrw", procs)
        outDir_refit_fits = "%s/apply_fit" % baseDir
        utils.mkdir(outDir_refit_fits, False)
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)


def zmumu_para_qT():

    procs = ["Zmumu"]
    tag = "zmumu_para_qT"
    comp = "para_qT"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "DY #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    doExport = True
    
    bhist = readProc(groups, "recoil_corr_xy_para_qT_qTbinned", procs)
    recoilMin, recoilMax = -100, 100
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pw_poly6_poly1"
        fitCfg['parms_init'] = [40, -1.64348e-01, 5.56260e-01, -4.26692e-02, 2.10091e-03, -5.89702e-05, 8.59665e-07, -5.01217e-09]
        fitCfg['parms_cte'] = [True, False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-2, yMax=16)
        return

    

    fitCfg = {} 
    fitCfg['func_name'] = "dy_para_qT"
    fitCfg['func_parms_vals'] = [15, 5, 8, 10, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 1, 1, 1] # 0=float, 1=propagate, 2=TF1
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        #fitF, params, cParams = "linear", [9.90024e-02, 9.34622e-01], [False, False, False]
        fitF, params, cParams = "power", [4.54931e+00, 1, -4.25497e+00], [False, True, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=50, yMax=50, xMax=100, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [4.66409e-01, 3.46158e-01, 3.44635e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=50, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [5.18556e-01, 3.95983e-01, 5.85e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=45, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", cutOffMin=6)
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=50, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)")
        
        
        #fitF, params, cParams = "linear", [9.90024e-02, 9.34622e-01], [False, False, False]
        fitF, params, cParams = "power", [4.54931e+00, 2.05145e-01, -1], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)")
        
        #fitF, params, cParams = "pw_poly4_poly1", [15, -4.26860e-01, 8.23660e-01, -8.92751e-02, 5.60430e-03, -1.38168e-04], [True, False, False, False, False, False]
        fitF, params, cParams = "power", [4.54931e+00, 2.05145e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)")
        
        fitF, params, cParams = "power", [4.54931e+00, 2.05145e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)")
        
        #fitF, params, cParams = "pw_poly4_poly1", [15, -4.26860e-01, 8.23660e-01, -8.92751e-02, 5.60430e-03, -1.38168e-04], [True, False, False, False, False, False]
        #fitF, params, cParams = "power", [4.54931e+00, 2.05145e-01, -4.25497e+00], [False, False, False]
        #rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=50, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
        #fitF, params, cParams = "pw_poly4_poly1", [15, -4.26860e-01, 8.23660e-01, -8.92751e-02, 5.60430e-03, -1.38168e-04], [True, False, False, False, False, False]
        #fitF, params, cParams = "power", [4.54931e+00, 2.05145e-01, -4.25497e+00], [False, False, False]
        #rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=50, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{4} (GeV)")
   
        # add norm parameters
        rls.addParam(jsOut, "p7", "[0]", [0.1])
        rls.addParam(jsOut, "p8", "[0]", [0.2])
        rls.addParam(jsOut, "p9", "[0]", [0.25])
        

        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return



    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = jsIn
        fitCfg['func_name'] = "dy_para_qT_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, chisq_refit=False, outDir=baseDir, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=50, yTitle = "#mu_{2} (GeV)")
        #rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        #rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
    if doExport:
        exportCfg = {}
        exportCfg["mean"] = ["p4", "p5", "p6", "p6"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p7", "p8", "p9"]
        rls.export(exportCfg, "%s/recoil_zmumu_para.json" % outCfgDir, "%s/results_refit.json" % baseDir) # , "%s/parametric_mean.json" % baseDir
   
    if False:
    
        bhist = readProc(groups, "recoil_corr_xy_para_qT_qTbinned_qTrw", procs)
        outDir_refit_fits = "%s/apply_fit" % baseDir
        utils.mkdir(outDir_refit_fits, False)
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax)
     
 

   
  
def singlemuon_perp():
    
    procs = ["SingleMuon"]
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "Data #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    doExport = True
    
    recoilMin, recoilMax = -100, 100
    bhist = readProc(groups, "recoil_corr_xy_perp_qTbinned", procs)
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames=procNames, yMin=-2, yMax=2)
        return

    fitCfg = {} 
    fitCfg['func_name'] = "data_perp"
    fitCfg['func_parms_vals'] = [40, 5, 10, 15, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1
    
    bkgCfg = {}
    bkgCfg['procs'] = ["ttbar", "ewk"]
    bkgCfg['parms'] = [utils.loadJSON("%s/ttbar_perp/results_refit.json"%outDir), utils.loadJSON("%s/ewk_perp/results_refit.json"%outDir)]
    bkgCfg['yields'] = [utils.loadJSON("%s/ttbar_perp/yields.json"%outDir), utils.loadJSON("%s/ewk_perp/yields.json"%outDir)]
    bkgCfg['norms'] = [1.0, 1.0]
    bkgCfg['data_yields'] = utils.loadJSON("%s/yields.json"%baseDir)
    
    
    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
        

        #jsIn = utils.loadJSON("%s/results.json" % outDir_param_v0)
        jsIn = utils.loadJSON("%s/zmumu_perp/results_refit.json" % outDir) # take the MC as starting values
        #jsIn['p7']['p0'] = 0.001 # perturb the value
        # single mean for all gauss
        jsIn['p4']['p0'] = 0.0
        jsIn['p4']['p1'] = 0.0
        jsIn['p5'] = copy.deepcopy(jsIn['p6'])
        jsIn['p6'] = copy.deepcopy(jsIn['p7'])
        jsIn['p7'] = copy.deepcopy(jsIn['p8'])
        
        jsIn['nParams'] = jsIn['nParams']-1
        del jsIn['p8']
        
        jsIn['p6'] = copy.deepcopy(jsIn['p7'])
        jsIn['nParams'] = jsIn['nParams']-1
        del jsIn['p7']
        
        jsIn['p5']['p0_bnds'] = (0, 1)
        jsIn['p6']['p0_bnds'] = (0, 1)
        
        
        # remove p0
        #jsIn['p0'] = copy.deepcopy(jsIn['p1'])
        #jsIn['p1'] = copy.deepcopy(jsIn['p2'])
        #jsIn['p2'] = copy.deepcopy(jsIn['p3'])
        #jsIn['p3'] = copy.deepcopy(jsIn['p4'])
        #jsIn['p4'] = copy.deepcopy(jsIn['p5'])
        #jsIn['p5'] = copy.deepcopy(jsIn['p6'])
        #jsIn['nParams'] = jsIn['nParams']-1
        #del jsIn['p6']
        
        fitCfg = jsIn
        fitCfg['func_name'] = "data_perp_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, chisq_refit=False, outDir=baseDir, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        #rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        #rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        #rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
        
    if doExport:
        exportCfg = {}
        exportCfg["mean"] = ["p4", "p4", "p4", "p4"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p5", 0.21, "p6"]
        rls.export(exportCfg, "%s/recoil_data_perp.json" % outCfgDir, "%s/results_refit.json" % baseDir)
         
    
def singlemuon_para_qT():

    procs = ["SingleMuon"]
    tag = "singlemuon_para_qT"
    comp = "para_qT"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "Data #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#parallel}" % (met)
    utils.mkdir(baseDir, False)
    rebin = 1
    doExport = True

    outDir_fits_v0 = "%s/fits_v0" % baseDir
    recoilMin, recoilMax = -100, 100
    
    bhist = readProc(groups, "recoil_corr_xy_para_qT_qTbinned", procs)
    procNames = groups.getProcNames(to_expand=procs)
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pw_poly4_power"
        fitCfg['parms_init'] = [2.37500e+01, 1.41561e-01, 4.30537e-01, -1.62911e-02, 4.85259e-04, -6.32757e-06, 9.75226e-01]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=0, yMax=20)
        return
        
    fitCfg = {} 
    fitCfg['func_name'] = "data_para_qT"
    fitCfg['func_parms_vals'] = [15, 5, 8, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1

    bkgCfg = {}
    bkgCfg['procs'] = ["ttbar", "ewk"]
    bkgCfg['parms'] = [utils.loadJSON("%s/ttbar_para_qT/results_refit.json"%outDir), utils.loadJSON("%s/ewk_para_qT/results_refit.json"%outDir)]
    bkgCfg['yields'] = [utils.loadJSON("%s/ttbar_para_qT/yields.json"%outDir), utils.loadJSON("%s/ewk_para_qT/yields.json"%outDir)]
    bkgCfg['norms'] = [1.0, 1.0]
    bkgCfg['data_yields'] = utils.loadJSON("%s/singlemuon_para_qT/yields.json"%outDir) # needed for the background fractions


    if True:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        #with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        with open("%s/zmumu_para_qT/results_refit.json" % outDir) as f: jsIn = json.load(f) # take the MC as starting values

        # collapse mean1 with mean2
        del jsIn['p4']
        jsIn['p4'] = copy.deepcopy(jsIn['p5'])
        jsIn['p5'] = copy.deepcopy(jsIn['p6'])
        jsIn['p6'] = copy.deepcopy(jsIn['p7'])
        jsIn['p7'] = copy.deepcopy(jsIn['p8'])
        jsIn['p8'] = copy.deepcopy(jsIn['p9'])

        # freeze n3
        #jsIn['p7'] = copy.deepcopy(jsIn['p8'])
        #jsIn['p8'] = copy.deepcopy(jsIn['p9'])        
        jsIn['nParams'] = jsIn['nParams']-2
        del jsIn['p8']
        

        
        fitCfg = jsIn
        fitCfg['func_name'] = "data_para_qT_cond"
        #jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, chisq_refit=False, outDir=baseDir, recoilLow=recoilMin, recoilHigh=recoilMax)
        #with open("%s/results_refit.json" % baseDir, "w") as outfile: json.dump(jsOut, outfile, indent=4)
         
        #with open("%s/results_refit.json" % baseDir) as f: jsIn = json.load(f)
        #rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, bkgCfg=bkgCfg, yMin=1e-3, yMax=1e6, recoilLow=recoilMin, recoilHigh=recoilMax, plotSignal=True)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=50, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=50, yTitle = "#mu_{2} (GeV)")
        #rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-20, yMax=50, yTitle = "#mu_{3} (GeV)")
        #rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{4} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
        #rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
    if doExport: 
        exportCfg = {}
        exportCfg["mean"] = ["p4", "p4", "p5", "p5"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p6", "p7", 0.5886822304617599]
        rls.export(exportCfg, "%s/recoil_data_para.json" % outCfgDir, "%s/results_refit.json" % baseDir)
       
 


def readProc(groups, hName, procs):

    label = "%s_%s" % (hName, procs[0]) 
    groups.setHists(hName, "", label=label, procsToRead=procs)
    bhist = groups.groups[procs[0]][label]
    return bhist 



    

    
if __name__ == "__main__":

    met = "RawPFMET" # DeepMETReso RawPFMET
    flavor = "mumu" # mu, e, mumu, ee

    groups = datagroupsLowPU("lowPU_%s_%s_nnpdf31.pkl.lz4" % (flavor, met), flavor=flavor)

    outCfgDir = f"wremnants/data/recoil/lowPU/{flavor}_{met}/"
    outDir = "/eos/user/j/jaeyserm/www/recoil/lowPU/%s_%s/" % (flavor, met)
    utils.mkdir(outDir, False)
    
    binning_qT = list(utils.drange(0, 100, 0.5)) + [100]

    import lowPU_RawPFMET_functions as rf
    rls.setFunctionLibs(rf)


    #zmumu_para_qT()
    #zmumu_perp()
        
    #ttbar_para_qT()
    #ttbar_perp()
    
    #ewk_para_qT()
    #ewk_perp()
        
    #zz_para_qT_RawPFMET()
    #zz_perp_RawPFMET()

    #singlemuon_para_qT()
    singlemuon_perp()
 