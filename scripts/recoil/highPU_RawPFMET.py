
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
from wremnants.datasets.datagroups2016 import make_datagroups_2016

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
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin, yRatio=1.25, sumw2=False) 
        return 
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [9.90024e-02, 1, 9.34622e-01], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "power", [4.66409e-01, 1, 9.34622e-01], [False, True, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "linear", [4.66409e-01, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")

        rls.addParam(jsOut, "p3", "[0]", [0.76])
        rls.addParam(jsOut, "p4", "[0]", [0.22])
     
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
    
    recoilMin, recoilMax = -150, 150
    rebin=2
    if False:
        
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.19057e+00, -2.37865e-01, 4.02425e-02, -1.34963e-03, 2.19717e-05, -1.73520e-07, 5.31624e-10]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=0, yMax=50)
        return
    
    fitCfg = {} 
    fitCfg['func_name'] = "ewk_para"
    fitCfg['func_parms_vals'] = [10, 25, 50, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0]
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25, sumw2=False)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [1, 1, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "power", [1, 1, 0], [False, True, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "power", [1, 1, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")


        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-50, yTitle = "#mu_{1} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-50, yTitle = "#mu_{2} (GeV)", fitOpts="NS W", cutOffMin=0)

        fitF, params, cParams = "linear", [0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=100, xMax=200, yMin=-50, yTitle = "#mu_{3} (GeV)", fitOpts="NS W", cutOffMin=0)


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
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25)
        return
        
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "quadratic", [9.90024e-02, 9.34622e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "quadratic", [60, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "quadratic", [60, 0, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=200, xMax=200, yTitle = "#sigma_{3} (GeV)")

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
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
        rls.doFitMultiGauss_plot(bhist, comp, jsIn, procLabel, metLabel, outDir_refit_fits, binning_qT, yMin=1e-4, yMax=1e4, rebin=rebin, recoilLow=recoilMin, recoilHigh=recoilMax, statUnc=False, yRatio=1.25)
        
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=200, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
   
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

    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pol6"
        fitCfg['parms_init'] = [2.99723e+00, -1.47337e-01, 1.61576e-02, -1.42718e-04, -9.91942e-07, 1.47306e-08, -3.50397e-11]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=50, recoilLow=recoilMin, recoilHigh=recoilMax, yRatio=1.25)
        return

    
    fitCfg = {} 
    fitCfg['func_name'] = "ttbar_para"
    fitCfg['func_parms_vals'] = [20, 50, 100, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0]
    
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-4, yMax=1e4, recoilLow=recoilMin, recoilHigh=recoilMax, rebin=rebin)
        return
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}

        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "quadratic", [9.90024e-02, 9.34622e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=40, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)", fitOpts="NS W")

        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=40, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{2} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=40, fitMax=100, yMax=200, xMax=200, yTitle = "#sigma_{3} (GeV)", fitOpts="NS W")

        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=40, fitMax=200, yMax=100, xMax=200, yMin=-100, yTitle = "#mu_{1} (GeV)", fitOpts="NS W")
        
        fitF, params, cParams = "quadratic", [4.66409e-01, 3.46158e-01, 0], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=40, fitMax=200, yMax=100, xMax=200, yMin=-100, yTitle = "#mu_{2} (GeV)", fitOpts="NS W")     

        rls.addParam(jsOut, "p5", "[0]", [0.5])
        rls.addParam(jsOut, "p6", "[0]", [0.4])
     
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
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "linear"
        fitCfg['parms_init'] = [0, 0]
        fitCfg['parms_cte'] = [False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-2, yMax=2)
        return


    fitCfg = {} 
    fitCfg['func_name'] = "dy_perp"
    fitCfg['func_parms_vals'] = [25, 5, 8, 10, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0]  # 0=float, 1=propagate, 2=TF1
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        jsIn = utils.loadJSON("%s/results.json" % outDir_fits_v0)
        
        fitF, params, cParams = "power", [1, 1, 1], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=100, xMax=200, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [0.5, 1, 11], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=50, xMax=200, yTitle = "#sigma_{2} (GeV)", cutOffMin=12.5)
        
        fitF, params, cParams = "power", [1, 1, 1], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=50, xMax=200, yTitle = "#sigma_{3} (GeV)")
        
        fitF, params, cParams = "power", [1, 1, 1], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=50, xMax=200, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "linear", [0, 0], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=50, xMax=200, yMin=-50, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "linear", [0, 0], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=200, yMax=50, xMax=200, yMin=-50, yTitle = "#mu_{2} (GeV)")

         
        # add norm parameters
        rls.addParam(jsOut, "p6", "[0]", [0.005])
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
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
      
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
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
 

    # do a refit absorbing correlations
    if True:

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
        exportCfg["mean"] = ["p4", "p5", "p4", "p5"]
        exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
        exportCfg["norm"] = ["p6", "p7", "p8"]
        
        rls.export(exportCfg, "%s/recoil_zmumu_perp.json" % outCfgDir, "%s/results_refit_v1.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit_v1/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        outDir_unc = "%s/uncertainties/" % baseDir
        utils.mkdir(outDir_unc, False)
        fitCfg = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
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
    
    if False:
        fitCfg = {}
        fitCfg['func_name'] = "pw_poly7_poly1"
        fitCfg['parms_init'] = [2.98078e+01, -2.50056e-02, 5.03451e-01, -2.86102e-02, 7.96316e-04, -8.18177e-07, -2.47077e-07, -2.73957e-09,  1.20691e-10]
        fitCfg['parms_cte'] = [False, False, False, False, False, False, False, False, False]
        rls.parameterizeMean(bhist, comp, baseDir, fitCfg, binning_qT, procLabel, metLabel, procNames, yMin=-5, yMax=25)
        return
        

    fitCfg = {} 
    fitCfg['func_name'] = "dy_para"
    fitCfg['func_parms_vals'] = [45, 10, 15, 25, 0, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1
    
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    if False:
        rls.doFitMultiGauss_fit(bhist, comp, fitCfg, procLabel, metLabel, outDir_fits_v0, binning_qT, yMin=1e-2, yMax=1e8, recoilLow=recoilMin, recoilHigh=recoilMax)
        return
        
        
    outDir_param_v0 = "%s/params_v0" % baseDir
    if False:
        jsOut = {}
        outDir_param = outDir_param_v0
        utils.mkdir(outDir_param, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "power", [4.54931e+00, 1, -4.25497e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [4.66409e-01, 1, 3.44635e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [5.18556e-01, 1, 15], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)")
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
    
        # add norm parameters
        rls.addParam(jsOut, "p8", "[0]", [0.05])
        rls.addParam(jsOut, "p9", "[0]", [0.15])
        rls.addParam(jsOut, "p10", "[0]", [0.25])

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
        fitCfg['func_name'] = "dy_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)
        
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
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
    
    procs = ["Data"]
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    procLabel = "Data #rightarrow #mu^{+}#mu^{#minus}"
    metLabel = "%s, U_{#perp}" % (met)
    utils.mkdir(baseDir, False)
    procNames = groups.getProcNames(to_expand=procs)
    doExport = True
    
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
    fitCfg['func_parms_vals'] = [40, 5, 10, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1
    
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
        
        fitF, params, cParams = "power", [1, 1, 13], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{1} (GeV)")

        fitF, params, cParams = "power", [1, 1, 5.5], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [1, 1, 7.5], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", cutOffMin=7.5)

        
        fitF, params, cParams = "linear", [0, 0], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=100, yMax=50, xMax=100, yMin=-5, cutOffMax=3, cutOffMin=-3, yTitle = "#mu_{1} (GeV)")
        

        # add norm parameters
        rls.addParam(jsOut, "p4", "[0]", [0.1])
        rls.addParam(jsOut, "p5", "[0]", [0.4])

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
       
        rls.plotParameter("p0", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        rls.plotParameter("p1", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{2} (GeV)")
        rls.plotParameter("p2", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{3} (GeV)")
        rls.plotParameter("p3", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=50, yTitle = "#sigma_{4} (GeV)")
        rls.plotParameter("p4", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)")
        rls.plotParameter("p5", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)")
        rls.plotParameter("p6", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        rls.plotParameter("p7", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{2} (GeV)")
        rls.plotParameter("p8", jsIn, outDir_refit, binning_qT, procLabel, metLabel, yMin=0, yMax=1, yTitle = "n_{3} (GeV)")
   
    # do a refit absorbing correlations
    if True:

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
    exportCfg["mean"] = ["p4", "p5", "p4", "p5"]
    exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
    exportCfg["norm"] = ["p6", "p7", "p8"]
 
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
    
        rls.export(exportCfg, "%s/recoil_data_perp.json" % outCfgDir, "%s/results_refit_v1.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit_v1/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        outDir_unc = "%s/uncertainties/" % baseDir
        utils.mkdir(outDir_unc, False)
        fitCfg = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)

        up = utils.loadJSON("%s/results_refit_ttbarUp.json" % baseDir)
        dw = utils.loadJSON("%s/results_refit_ttbarDown.json" % baseDir)
        rls.plotSystUncertainty(fitCfg, up, dw, "syst_ttbar", "TTbar up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
        
        up = utils.loadJSON("%s/results_refit_ewkUp.json" % baseDir)
        dw = utils.loadJSON("%s/results_refit_ewkDown.json" % baseDir)
        rls.plotSystUncertainty(fitCfg, up, dw, "syst_ewk", "EWK up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
   
  
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
    fitCfg['func_parms_vals'] = [45, 10, 15, 25, 0, 0, 0, 0]
    fitCfg['func_parms_cfg'] = [0, 0, 0, 0, 0, 0, 0, 0] # 0=float, 1=propagate, 2=TF1


    bkgCfg = {}
    bkgCfg['procs'] = ["ttbar", "ewk"]
    bkgCfg['parms'] = [utils.loadJSON("%s/ttbar_para/results_refit.json"%outDir), utils.loadJSON("%s/ewk_para/results_refit.json"%outDir)]
    bkgCfg['yields'] = [utils.loadJSON("%s/ttbar_para/yields.json"%outDir), utils.loadJSON("%s/ewk_para/yields.json"%outDir)]
    bkgCfg['norms'] = [1.0, 1.0]
    bkgCfg['data_yields'] = utils.loadJSON("%s/singlemuon_para/yields.json"%outDir) # needed for the background fractions


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
        
        fitF, params, cParams = "power", [4.54931e+00, 1, -4.25497e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p0", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=150, yMax=100, xMax=100, yTitle = "#sigma_{1} (GeV)", cutOffMin=25)

        fitF, params, cParams = "power", [4.66409e-01, 1, 3.44635e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p1", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=150, yMax=50, xMax=100, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "power", [5.18556e-01, 1, 16], [False, False, True]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p2", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=150, yMax=50, xMax=100, yTitle = "#sigma_{3} (GeV)", cutOffMin=15)
        
        fitF, params, cParams = "power", [3.01503e-01, 5.26286e-01, 8.13076e+00], [False, False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p3", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=150, yMax=50, xMax=100, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p4", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p5", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{2} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p6", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
        fitF, params, cParams = "linear", [0, 1], [False, False]
        rls.parameterizeGauss(jsIn, jsOut, comp, "p7", fitF, params, outDir_param, binning_qT, procLabel, metLabel, cParams=cParams, fitMin=0, fitMax=80, yMax=50, xMax=100, yMin=-20, yTitle = "#mu_{3} (GeV)")
        
    
        # add norm parameters
        rls.addParam(jsOut, "p8", "[0]", [0.05])
        rls.addParam(jsOut, "p9", "[0]", [0.15])
        rls.addParam(jsOut, "p10", "[0]", [0.25])


        jsOut['nParams'] = len(jsOut)
        utils.writeJSON("%s/results.json" % outDir_param_v0, jsOut)
        return


    if False:
        outDir_refit = "%s/params_refit" % baseDir
        outDir_refit_fits = "%s/fits_refit" % baseDir
        utils.mkdir(outDir_refit, False)
        utils.mkdir(outDir_refit_fits, False)
    
        #fitCfg = utils.loadJSON("%s/results.json" % outDir_param_v0)
        fitCfg = utils.loadJSON("%s/zmumu_para/results_refit.json" % outDir)
        

        fitCfg['func_name'] = "data_para_cond"
        jsOut = rls.combinedFit_scipy(bhist, comp, fitCfg, binning_qT, bkgCfg=bkgCfg, outDir=outDir_refit, recoilLow=recoilMin, recoilHigh=recoilMax)
        utils.writeJSON("%s/results_refit.json" % baseDir, jsOut)

         
        jsIn = utils.loadJSON("%s/results_refit.json" % baseDir)
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


    # do a refit absorbing correlations
    if True:

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
    exportCfg["mean"] = ["p4", "p5", "p6", "p7"]
    exportCfg["sigma"] = ["p0", "p1", "p2", "p3"]
    exportCfg["norm"] = ["p8", "p9", "p10"]
 
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
    
        rls.export(exportCfg, "%s/recoil_data_para.json" % outCfgDir, "%s/results_refit_v1.json" % baseDir, jsInF_mean="%s/parametric_mean.json" % baseDir)
        rls.compareChi2("%s/results.json" % outDir_fits_v0, "%s/fits_refit_v1/results.json" % baseDir, "%s/chi2_comparison" % baseDir, procLabel, metLabel)
    
        outDir_unc = "%s/uncertainties/" % baseDir
        utils.mkdir(outDir_unc, False)
        fitCfg = utils.loadJSON("%s/results_refit_v1.json" % baseDir)
        yieldsIn = utils.loadJSON("%s/yields.json" % baseDir)
        #rls.plotStatUncertainties(fitCfg, yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)

        up = utils.loadJSON("%s/results_refit_ttbarUp.json" % baseDir)
        dw = utils.loadJSON("%s/results_refit_ttbarDown.json" % baseDir)
        rls.plotSystUncertainty(fitCfg, up, dw, "syst_ttbar", "TTbar up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
        
        up = utils.loadJSON("%s/results_refit_ewkUp.json" % baseDir)
        dw = utils.loadJSON("%s/results_refit_ewkDown.json" % baseDir)
        rls.plotSystUncertainty(fitCfg, up, dw, "syst_ewk", "EWK up/down", yieldsIn, exportCfg, comp, procLabel, metLabel, outDir_unc, binning_qT, bkgCfg=bkgCfg, recoilLow=recoilMin, recoilHigh=recoilMax, yMin=1e-2, yMax=1e8)
  
  

def readProc(groups, hName, procs):

    label = "%s_%s" % (hName, procs[0]) 
    groups.setHists(hName, "", label=label, procsToRead=procs)
    bhist = groups.groups[procs[0]][label]
    return bhist 



    

    
if __name__ == "__main__":

    met = "RawPFMET"
    flavor = "mumu"
    rls.__topRight__ = "16.8 fb^{#minus1} (13 TeV)"
    
    groups = make_datagroups_2016("mz_wlike_with_mu_eta_pt_%s.pkl.lz4" % met)
    groups.groups.update({
        "TTbar" : dict(
                members = [groups.datasets[x] for x in ["TTLeptonicPostVFP", "TTSemileptonicPostVFP", "SingleTschanLepDecaysPostVFP", "SingleTtWAntitopPostVFP", "SingleTtchanAntitopPostVFP", "SingleTtchanTopPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
        ),
        "EWK" : dict(
                members = [groups.datasets[x] for x in ["WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "ZZ2l2nuPostVFP", "WZPostVFP", "WWPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
        )
    })

    outCfgDir = f"wremnants/data/recoil/highPU/{flavor}_{met}/"
    outDir = "/eos/user/j/jaeyserm/www/recoil/highPU/%s_%s/" % (flavor, met)
    utils.mkdir(outDir, False)
    
    binning_qT = list(utils.drange(0, 200, 1)) + [200]

    import highPU_RawPFMET_functions as rf
    rls.setFunctionLibs(rf)


    #zmumu_para()
    zmumu_perp()
        
    #ttbar_para()
    #ttbar_perp()
    
    #ewk_para()
    #ewk_perp()

    #singlemuon_para()
    #singlemuon_perp()
 