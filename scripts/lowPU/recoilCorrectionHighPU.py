
import sys,array,math,os
import numpy as np
import ctypes
import json

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import functions
import plotter
import recoilLibs

from wremnants.datasets.datagroups import datagroups2016

import lz4.frame
import narf


def do_zz_perp_RawPFMET():

    proc = "ZZ"
    tag = "zz_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "ZZ #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 30, 50], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.50, 0.40], [0, 0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
        
        
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_zz_para_RawPFMET():

    proc = "ZZ"
    tag = "zz_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "ZZ #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 30, 50], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.50, 0.40], [0, 0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=False, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[1]*x + [0]", [0, 0], [False, False]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)
 
def do_ewk_perp_RawPFMET():

    proc = "EWK"
    tag = "ewk_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [1, 1]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [20, 50], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.9], [1]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
        
        
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_ewk_para_RawPFMET():

    proc = "EWK"
    tag = "ewk_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 50], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.95], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[1]*x + [0]", [0, 0], [False, False]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)

def do_ttbar_perp_RawPFMET():

    proc = "TTbar"
    tag = "ttbar_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [40, 70], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*x*x + [1]*x + [2]", [0, 0, 0], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*x*x + [1]*x + [2]", [0, 0, 0], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=150)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if False:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_ttbar_para_RawPFMET():

    proc = "TTbar"
    tag = "ttbar_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [40, 70], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return

    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
                
        # mean
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=-20, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=-20, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=150, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if False:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)

def do_dymumu_perp_RawPFMET():

    proc = "Zmumu"
    tag = "dymumu_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "DY #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}   
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 15, 20, 35], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.35, 0.25, 0.2], [0, 0, 0]

    #fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    #fitCfg['sigma'], fitCfg['sigma_cfg'] = [8, 15, 30], [0, 0, 0]
    #fitCfg['norm'], fitCfg['norm_cfg'] = [0.2, 0.3], [0, 0]
    
    nGauss = len(fitCfg['mean'])
    
    if False: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
    ###### STEP 2: parameterize norm1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)

        fitF, params, cParams = "[0]*x + [1]", [-.0001, 0.6], [False, False]
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.3, 0.01], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")

        fitF, params, cParams = "[0]*x + [1]", [0, 0.6], [False, False]
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.25, -0.01], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=90, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [0, 0.6], [False, False]
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.25, -0.01], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm3", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        return
        
        
    ###### STEP 3: fit with parameterized norm1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0, 0, 0]
        fitCfg['norm_cfg'] = [1, 1, 1]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v1, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 4: parameterize norm2
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)

        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [2.00274e-02, 9.14330e-01, 1.55842e+01, -7.11068e-05], [False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1, 1, 15], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma4", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.3, 0.01], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")

        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.25, -0.0002], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=90, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.25, -0.0002], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm3", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        
        
        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        return
        
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 0, 2, 0]
        fitCfg['norm_cfg'] = [2, 2, 2]
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v2, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v2 = "%s/param_v2" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "[0]", [0], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#mu_{3} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean4", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#mu_{4} (GeV)", doFit=False)

        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [5.69585e-01, 2.36703e-01, 1.26615e+01, 1.03245e-05], [True, False, True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.27155e-01, 6.58167e-01, 1.72465e+01], [False, True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [-9.06199e-05, 1.69713e+00, 2.04148e+01], [True, False, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [4.29187e-02, 9.24539e-01, 2.40392e+01], [False, True, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma4", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [3.15032e-01, 9.16645e-04], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")

        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.34331e-01, -3.01182e-04], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=90, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.82460e-01, -2.75679e-04], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm3", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{3} (GeV)")
        
        

        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
        
        
    
        
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, qTmax=150)
        with open("%s/results_refit.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2, 2]
        fitCfg['norm_cfg'] = [2, 2, 2]

        outDir_refit = "%s/refit_v2" % baseDir
        with open("%s/results_refit.json" % outDir_param_v2) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v2_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "n_{1} (GeV)")
        recoilLibs.doPlot(comp, "norm2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "n_{2} (GeV)")
        recoilLibs.doPlot(comp, "norm3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "n_{3} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        recoilLibs.doPlot(comp, "sigma3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        recoilLibs.doPlot(comp, "sigma4", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{4} (GeV)")

def do_dymumu_para_RawPFMET():

    proc = "Zmumu"
    tag = "dymumu_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "DY #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [0, 0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 15, 20, 35], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.35, 0.25, 0.2], [0, 0, 0]

    '''
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [8, 15, 30], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.2, 0.3], [0, 0]
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 20, 30], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.4, 0.3], [0, 0]
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 20, 30], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.35, 0.3], [0, 0]
    
    
    '''
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150, propagate=False)
        return

    
    ###### STEP 2: parameterize norm1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        # means
        #fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 - [4]*20) + [4]*x   )", [2.22102e-01, 3.76322e-01, -1.31296e-02, 1.89339e-04, 6.49759e-02], [False, False, False, False, False]
        fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 + [4]*20*20*20*20 - [5]*20) + [5]*x   )", [7.45825e-02, 4.74556e-01, -2.82073e-02, 1.05968e-03, -1.61073e-05, 6.47372e-02], [False, False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{3} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean4", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{4} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [7.26314e-01, 2.49241e-01, 4.41560e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma4", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [3.20601e-01, 7.92963e-04], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{1} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.50422e-01, 1.05541e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [3.34560e-01, 5.32842e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm3", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{3} (GeV)")
        
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        return
        
        
    ###### STEP 3: fit with parameterized norm1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0, 0, 0]
        fitCfg['norm_cfg'] = [2, 2, 2]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v1, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 4: parameterize sigma2
    outDir_param_v1 = "%s/param_v1" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        # means
        #fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 - [4]*20) + [4]*x   )", [2.22102e-01, 3.76322e-01, -1.31296e-02, 1.89339e-04, 6.49759e-02], [False, False, False, False, False]
        fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 + [4]*20*20*20*20 - [5]*20) + [5]*x   )", [7.45825e-02, 4.74556e-01, -2.82073e-02, 1.05968e-03, -1.61073e-05, 6.47372e-02], [False, False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{3} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean4", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{4} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.20049e-01, 7.81115e-01, 1.31185e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [3.62778e-03, 1.33033e+00, 1.80832e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.66941e-01, 5.78182e-01, 1.99598e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [5.04482e-03, 1.38706e+00, 2.52433e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma4", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{4} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [3.20601e-01, 7.92963e-04], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{1} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.50422e-01, 1.05541e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [3.34560e-01, 5.32842e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm3", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "n_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
       
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, qTmax=450, singleMean=True)
        with open("%s/results_refit.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #return
       
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2, 2]
        fitCfg['norm_cfg'] = [2, 2, 2]

        outDir_refit = "%s/refit_v1" % baseDir
        with open("%s/results_refit.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v1_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "mean1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)")
        recoilLibs.doPlot(comp, "mean2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)")
        recoilLibs.doPlot(comp, "mean3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{3} (GeV)")
        recoilLibs.doPlot(comp, "mean4", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{3} (GeV)")
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#n_{1} (GeV)")
        recoilLibs.doPlot(comp, "norm2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#n_{2} (GeV)")
        recoilLibs.doPlot(comp, "norm3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#n_{3} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        recoilLibs.doPlot(comp, "sigma3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        recoilLibs.doPlot(comp, "sigma4", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        return
        
        
        
       
       
       
       
       
       
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 2, 0]
        fitCfg['norm_cfg'] = [2, 2]
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v2, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v2 = "%s/param_v2" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)
        
        # means
        #fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 - [4]*20) + [4]*x   )", [2.22102e-01, 3.76322e-01, -1.31296e-02, 1.89339e-04, 6.49759e-02], [False, False, False, False, False]
        fitF, params, cParams = "(x<25)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x) +0+ (x>25)*1*(  ([0] + [1]*25 + [2]*25*25 + [3]*25*25*25 + + [4]*25*25*25*25 - [5]*25) + [5]*x   )", [-7.04087e-03, 5.05113e-01, -3.14911e-02, 1.13700e-03, -1.60587e-05, 6.65909e-02], [False, False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [6.25579e-02, 8.54263e-01, 1.21481e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.05814e-02, 1.06093e+00, 1.93983e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.42627e-01, 2.77646e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [7.22781e-01, 9.43120e-04], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{2} (GeV)")
        
        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        return
        
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v3 = "%s/fits_v3" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 0]
        fitCfg['norm_cfg'] = [2, 2]
        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v3, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v3 = "%s/param_v3" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v3, True)
        with open("%s/results.json" % outDir_fits_v3) as f: jsIn = json.load(f)
        
        # means
        #fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*(  ([0] + [1]*20 + [2]*20*20 + [3]*20*20*20 - [4]*20) + [4]*x   )", [2.22102e-01, 3.76322e-01, -1.31296e-02, 1.89339e-04, 6.49759e-02], [False, False, False, False, False]
        fitF, params, cParams = "(x<25)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x) +0+ (x>25)*1*(  ([0] + [1]*25 + [2]*25*25 + [3]*25*25*25 + + [4]*25*25*25*25 - [5]*25) + [5]*x   )", [-7.04087e-03, 5.05113e-01, -3.14911e-02, 1.13700e-03, -1.60587e-05, 6.65909e-02], [False, False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [6.25579e-02, 8.54263e-01, 1.21481e+01], [True, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.05814e-02, 1.06093e+00, 1.93983e+01], [True, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        fitF, params, cParams = "[0] + [1]*x + [2]*x*x + [3]*x*x*x", [0, 0, 0], [False, False, False, False]
        fitF, params, cParams = "(x<60)*1*([0] + [1]*x + [2]*x*x) +0+ (x>60)*1*( [0] + [1]*60 + [2]*60*60 -[3]*60 + [3]*x  )", [2.99481e+01, -9.80642e-02, 1.17718e-03, 4.73951e-02], [False, False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [2.42627e-01, 2.77646e-03], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{1} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [7.22781e-01, 9.43120e-04], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v3, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=1, xMax=500, yTitle = "#mu_{2} (GeV)")
        
        with open("%s/results.json" % outDir_param_v3, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
        
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v3) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, qTmax=450, singleMean=True)
        with open("%s/results_refit.json" % outDir_param_v3, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2]
        fitCfg['norm_cfg'] = [2, 2]

        outDir_refit = "%s/refit_v3" % baseDir
        with open("%s/results_refit.json" % outDir_param_v3) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v3_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "mean1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)")
        recoilLibs.doPlot(comp, "mean2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)")
        recoilLibs.doPlot(comp, "mean3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{3} (GeV)")
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "norm2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{2} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        recoilLibs.doPlot(comp, "sigma3", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{3} (GeV)")
        

def do_singlemuon_perp_RawPFMET_2Gauss():

    proc = "Data"
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "SingleMuon #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [15, 30], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    

    bkgCfg = {}
    bkgCfg['TTbar'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ttbar_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['EWK'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ewk_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['ZZ'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/zz_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    nGauss = len(fitCfg['mean'])
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
    
    ###### STEP 2: parameterize norm1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [7.26314e-01, 2.49241e-01, 4.41560e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [0, 0.6], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        return
        
        
    ###### STEP 3: fit with parameterized norm1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v1, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 4: parameterize sigma2
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [7.26314e-01, 2.49241e-01, 4.41560e+00, 0], [False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=200, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [-0.001, 0.6], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        return
        
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 2]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v2, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v2 = "%s/param_v2" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [2.00274e-02, 9.14330e-01, 1.55842e+01, -7.11068e-05], [False, True, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=450, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, True, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=300, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [-1.21169e-03, 5.89889e-01], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
        
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, bkgCfg=bkgCfg, qTmax=450)
        with open("%s/results_refit.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [2]

        outDir_refit = "%s/refit_v2" % baseDir
        with open("%s/results_refit.json" % outDir_param_v2) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v2_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
 

def do_singlemuon_perp_RawPFMET():

    proc = "Data"
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "SingleMuon #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [15, 30], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    

    bkgCfg = {}
    bkgCfg['TTbar'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ttbar_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['EWK'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ewk_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['ZZ'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/zz_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    nGauss = len(fitCfg['mean'])
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
    
    ###### STEP 2: parameterize norm1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [7.26314e-01, 2.49241e-01, 4.41560e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [0, 0.6], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        return
        
        
    ###### STEP 3: fit with parameterized norm1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v1, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 4: parameterize sigma2
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [7.26314e-01, 2.49241e-01, 4.41560e+00, 0], [False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=200, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [-0.001, 0.6], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        return
        
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 2]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v2, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v2 = "%s/param_v2" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=False)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=False)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [2.00274e-02, 9.14330e-01, 1.55842e+01, -7.11068e-05], [False, True, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=450, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, True, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=300, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*x + [1]", [-1.21169e-03, 5.89889e-01], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
        
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, bkgCfg=bkgCfg, qTmax=450)
        with open("%s/results_refit.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [2]

        outDir_refit = "%s/refit_v2" % baseDir
        with open("%s/results_refit.json" % outDir_param_v2) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v2_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
 
 
 
def do_singlemuon_para_RawPFMET():

    proc = "Data"
    tag = "singlemuon_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "SingleMuon #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1
    #rebin = [-300, -100, -80, -70, -60, -50, -46, -42, -38, -34] + list(range(-30, 30, 1)) + [30, 34, 38, 42, 46, 50, 60, 70, 80, 100, 300]


    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [6, 10, 15], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.3, 0.6], [1, 1]
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [5, 15], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.8], [0]
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [15, 30], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]

    bkgCfg = {}
    bkgCfg['TTbar'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ttbar_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['EWK'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/ewk_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['ZZ'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_RawPFMET/zz_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    nGauss = len(fitCfg['mean'])
    if False: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        sys.exit()

    
    ###### STEP 2: parameterize norm1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*([4] + [5]*x)", [0, 0, 0, 0, 0, 0], [False, False, False, False, False, False]
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(  ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x   )", [7.22670e-02, 4.21357e-01, -1.14848e-02, 1.47450e-04, 9.26982e-02], [False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [7.26314e-01, 2.49241e-01, 4.41560e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [5.93495e-01, 3.62521e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v0, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        return
        
        
    ###### STEP 3: fit with parameterized norm1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v1, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 4: parameterize sigma2
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*([4] + [5]*x)", [0, 0, 0, 0, 0, 0], [False, False, False, False, False, False]
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(  ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x   )", [7.22670e-02, 4.21357e-01, -1.14848e-02, 1.47450e-04, 9.26982e-02], [False, False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [4.06103e-02, 9.18130e-01, 1.54776e+01, -2.53522e-05], [False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.15048e-02, 9.60879e-01, 2.26090e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=200, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [5.93495e-01, 3.62521e-03], [False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        return
        
    ###### STEP 5: fit with parameterized sigma1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 0]
        fitCfg['norm_cfg'] = [2]
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v2, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
        
        
    ###### STEP 6: parameterize other parameters
    outDir_param_v2 = "%s/param_v2" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(  ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x   )", [7.22670e-02, 4.21357e-01, -1.14848e-02, 1.47450e-04, 9.26982e-02], [False, True, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)", doFit=True)
        recoilLibs.doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2] + [3]*x*x", [4.06103e-02, 9.18130e-01, 1.54776e+01, -2.53522e-05], [False, False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
 
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.02607e-02, 1.24859e+00, 2.25236e+01], [False, False, False]
        recoilLibs.doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=500, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [5.93495e-01, 3.62521e-03], [True, True]
        recoilLibs.doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, recoil_qTbins, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        #return
        
    ###### STEP 7: refit
    if True:
        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        jsOut = recoilLibs.combinedFit(rFile, proc, comp, jsIn, recoil_qTbins, rebin=rebin, bkgCfg=bkgCfg, qTmax=450, singleMean=True)
        with open("%s/results_refit.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [2]

        outDir_refit = "%s/refit_v2" % baseDir
        with open("%s/results_refit.json" % outDir_param_v2) as f: jsIn = json.load(f)
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_refit, recoil_qTbins, funcJs=jsIn, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, xMin=-150, xMax=150)
 
        outDir_refit = "%s/param_v2_refit" % baseDir
        functions.prepareDir(outDir_refit)
        recoilLibs.doPlot(comp, "mean1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{1} (GeV)")
        recoilLibs.doPlot(comp, "mean2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#mu_{2} (GeV)")
        recoilLibs.doPlot(comp, "norm1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=1, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma1", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=30, xMax=500, yTitle = "#sigma_{1} (GeV)")
        recoilLibs.doPlot(comp, "sigma2", jsIn, outDir_refit, recoil_qTbins, label, yMin=0, yMax=50, xMax=500, yTitle = "#sigma_{2} (GeV)")
        
        
        
        
 
def do_zz_perp_DeepMETReso():

    proc = "ZZ"
    tag = "zz_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "ZZ #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 30, 50], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.50, 0.40], [0, 0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
        
        
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_zz_para_DeepMETReso():

    proc = "ZZ"
    tag = "zz_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "ZZ #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0], [0, 0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 30, 50], [0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.50, 0.40], [0, 0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=False, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[1]*x + [0]", [0, 0], [False, False]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)
 
def do_ewk_perp_DeepMETReso():

    proc = "EWK"
    tag = "ewk_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [1, 1]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [20, 50], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.9], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
        
        
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_ewk_para_DeepMETReso():

    proc = "EWK"
    tag = "ewk_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "EWK #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [10, 50], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.95], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=False, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.69746e-01, 5.59970e-01, 8], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[1]*x + [0]", [0, 0], [False, False]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=80, yMin=-5, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=50, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)

def do_ttbar_perp_DeepMETReso():

    proc = "TTbar"
    tag = "ttbar_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [40, 70], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-200, xMax=200)
        return
    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*x*x + [1]*x + [2]", [0, 0, 0], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*x*x + [1]*x + [2]", [0, 0, 0], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
       
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=150)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if False:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200)

def do_ttbar_para_DeepMETReso():

    proc = "TTbar"
    tag = "ttbar_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "TTbar #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 2

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0], [0, 0]
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [40, 70], [0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.5], [0]
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=False, xMin=-200, xMax=200)
        return

    
    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=150, yMin=0, yMax=100, yTitle = "#sigma_{2} (GeV)")
                
        # mean
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=-20, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        fitF, params, cParams = "[0]*x*x*x + [1]*x*x + [2]*x + [3]", [0, 0, 0, 0], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=150, yMin=-20, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if False:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(bkg, comp, jsIn, rebin=rebin, qTmax=150, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if False:
        fitCfg['mean_cfg'] = [1, 1]
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(bkg, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, xMin=-200, xMax=200, singleMean=True)

 
def do_dymumu_perp_DeepMETReso():

    proc = "Zmumu"
    tag = "dymumu_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "DY #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}

    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [25, 5, 9, 14], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.007, 0.30, 0.55], [1, 0, 0]
    
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return
        
        
  

    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.98628e-01, 5.72942e-01, 3.83471e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.45588e-01, 6.81281e-01, 7.02926e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{2} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [4.61160e-01, 4.48764e-01, 9.41019e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{3} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.24638e-01, 1, 1.36732e+01], [False, False, False]
        doFit(jsIn, jsOut, "sigma4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{4} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.24638e-01, 1, 1.36732e+01], [False, False, False]
        doFit(jsIn, jsOut, "sigma5", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=10, yMax=50, yTitle = "#sigma_{5} (GeV)")

        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][2]], [True]
        doFit(jsIn, jsOut, "mean3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][3]], [True]
        doFit(jsIn, jsOut, "mean4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{4} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][4]], [True]
        doFit(jsIn, jsOut, "mean5", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{5} (GeV)", doFit=False)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][1]], [True]
        doFit(jsIn, jsOut, "norm2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][2]], [True]
        doFit(jsIn, jsOut, "norm3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][3]], [True]
        doFit(jsIn, jsOut, "norm4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", doFit=False)

        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(proc, comp, jsIn, rebin=rebin)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2, 2]
        fitCfg['norm_cfg'] = [1, 1, 1]
        
        fitCfg['mean_cfg'] = [2, 2, 2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2, 2, 2]
        fitCfg['norm_cfg'] = [1, 1, 1, 1]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin)


def do_dymumu_para_DeepMETReso_old():

    proc = "Zmumu"
    tag = "dymumu_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "DY #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    outDir_fits = "%s/fits_v0" % baseDir
    fitCfg = {}

    ## 1 GeV; floating means
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [0, 0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [30, 5, 9, 14], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.005, 0.30, 0.55], [1, 0, 0]
    
    # 0.5 GeV
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [0, 0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [30, 5, 8, 10], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.005, 0.25, 0.25], [0, 0, 0]
    
    nGauss = len(fitCfg['mean'])
    
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits, recoil_qTbins, qTmax=500, rebin=rebin, singleMean=True, xMin=-150, xMax=150)
        return

    ###### STEP 2: parameterize
    outDir_param = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param, False)
        with open("%s/results.json" % outDir_fits) as f: jsIn = json.load(f)
        
        # sigma
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.98628e-01, 5.72942e-01, 3.83471e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.45588e-01, 6.81281e-01, 7.02926e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma2", fitF, params, outDir_param, label,cParams=cParams,  fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{2} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [4.61160e-01, 4.48764e-01, 9.41019e+00], [False, False, False]
        doFit(jsIn, jsOut, "sigma3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{3} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.24638e-01, 1, 1.36732e+01], [False, False, False]
        doFit(jsIn, jsOut, "sigma4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#sigma_{4} (GeV)")
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.24638e-01, 1, 1.36732e+01], [False, False, False]
        doFit(jsIn, jsOut, "sigma5", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=10, yMax=50, yTitle = "#sigma_{5} (GeV)")

        # mean
        fitF, params, cParams = "(x<20)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>20)*1*([4] + [5]*x)", [0, 0, 0, 0, 0, 0], [False, False, False, False, False, False]
        doFit(jsIn, jsOut, "mean1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=20, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#mu_{2} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#mu_{3} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#mu_{4} (GeV)", doFit=True)
        doFit(jsIn, jsOut, "mean5", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=20, yTitle = "#mu_{5} (GeV)", doFit=True)

        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, "norm1", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][1]], [True]
        doFit(jsIn, jsOut, "norm2", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][2]], [True]
        doFit(jsIn, jsOut, "norm3", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{3} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][3]], [True]
        doFit(jsIn, jsOut, "norm4", fitF, params, outDir_param, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{4} (GeV)", doFit=False)
        
        with open("%s/results.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)
    

    ###### STEP 3: refit
    if True:
        with open("%s/results.json" % outDir_param) as f: jsIn = json.load(f)
        jsOut = combinedFit(proc, comp, jsIn, rebin=rebin, singleMean=True)
        with open("%s/results_refit.json" % outDir_param, "w") as outfile: json.dump(jsOut, outfile, indent=4)


    ###### STEP 4: validate
    if True:
        fitCfg['mean_cfg'] = [2, 2, 2, 2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2, 2, 2]
        fitCfg['norm_cfg'] = [2, 2, 2, 2]
        outDir_refit = "%s/refit_v0" % baseDir
        with open("%s/results_refit.json" % outDir_param) as f: jsIn = json.load(f)
        #doFitMultiGauss(proc, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin)

        outDir_refit = "%s/param_v0_refit" % baseDir
        functions.prepareDir(outDir_refit)
        doPlot(comp, "sigma1", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        doPlot(comp, "sigma2", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{2} (GeV)")
        doPlot(comp, "sigma3", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{3} (GeV)")
        doPlot(comp, "sigma4", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{4} (GeV)")
        doPlot(comp, "sigma5", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{5} (GeV)")
        doPlot(comp, "mean1", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#mu_{1} (GeV)")
        doPlot(comp, "mean2", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#mu_{2} (GeV)")
        doPlot(comp, "mean3", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#mu_{3} (GeV)")

 
def do_singlemuon_perp_DeepMETReso():

    proc = "Data"
    tag = "singlemuon_perp"
    comp = "perp"
    baseDir = "%s/%s" % (outDir, tag)
    name = "SingleMuon #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1

    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0,0], [0, 0, 0,0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [35, 5, 8, 20], [0, 0, 0,0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.005, 0.1, 0.6], [1, 0,0]
    
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [1, 1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [25, 5, 9, 14], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.007, 0.30, 0.55], [1, 0, 0]

    bkgCfg = {}
    bkgCfg['TTbar'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/ttbar_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['EWK'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/ewk_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['ZZ'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/zz_perp/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    nGauss = len(fitCfg['mean'])
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, xMin=-150, xMax=150)
        return
    
    ###### STEP 2: parameterize sigma1
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [7.26314e-01, 2.49241e-01, 4.41560e+00], [False, False, False]
        ###doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
 
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=20, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1, 0, 10], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=30, yTitle = "#sigma_{3} (GeV)")
        
        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        sys.exit()
        
        
    ###### STEP 3: fit with parameterized sigma1
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 0, 0]
        fitCfg['norm_cfg'] = [1, 1]
        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_fits_v1, funcJs=jsIn, qTmax=200, rebin=rebin, bkgCfg=bkgCfg)    
        sys.exit()
        
        
    ###### STEP 4: parameterize other parameters
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [9.62059e+00, 2.88035e-02, -4.51838e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=60, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=20, yTitle = "#sigma_{2} (GeV)")
        
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.33181e-02, 1.33470e+00, 1.35768e+01], [False, True, False]
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [1.33181e-02, 1., 1.35768e+01], [False, True, False]
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [8.121e-02, 1., 1.306e+01], [False, True, False]
        doFit(jsIn, jsOut, comp, "sigma3", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=30, yTitle = "#sigma_{3} (GeV)")
        #sys.exit()
        # mean
        fitF, params, cParams = "[0]", [fitCfg['mean'][0]], [True]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][1]], [True]
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{2} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['mean'][2]], [True]
        doFit(jsIn, jsOut, comp, "mean3", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=50, yMin=-5, yMax=5, yTitle = "#mu_{3} (GeV)", doFit=False)


        # norm
        fitF, params, cParams = "[0]", [fitCfg['norm'][0]], [True]
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)
        fitF, params, cParams = "[0]", [fitCfg['norm'][1]], [True]
        doFit(jsIn, jsOut, comp, "norm2", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=50, yMin=0, yMax=1, yTitle = "n_{2} (GeV)", doFit=False)

        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)    
        sys.exit()
        
    ###### STEP 5: refit
    if True:
        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        jsOut = combinedFit(proc, comp, jsIn, rebin=rebin, bkgCfg=bkgCfg, qTmax=50)
        with open("%s/results_refit.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
        
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [1, 1, 1] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2, 2]
        fitCfg['norm_cfg'] = [1, 1]

        outDir_refit = "%s/refit_v1" % baseDir
        with open("%s/results_refit.json" % outDir_param_v1) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=200, rebin=rebin, bkgCfg=bkgCfg)
 
        outDir_refit = "%s/param_v1_refit" % baseDir
        functions.prepareDir(outDir_refit)
        doPlot(comp, "sigma1", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#sigma_{1} (GeV)")
        doPlot(comp, "sigma2", jsIn, outDir_refit, label, yMin=0, yMax=30, yTitle = "#sigma_{2} (GeV)")
        doPlot(comp, "sigma3", jsIn, outDir_refit, label, yMin=0, yMax=30, yTitle = "#sigma_{3} (GeV)")
 
def do_singlemuon_para_DeepMETReso_old():

    proc = "Data"
    tag = "singlemuon_para"
    comp = "para"
    baseDir = "%s/%s" % (outDir, tag)
    name = "SingleMuon #rightarrow #mu^{+}#mu^{#minus}"
    label = "%s, %s" % (name, met)
    functions.prepareDir(baseDir, False)
    rebin = 1
    #rebin = [-300, -100, -80, -70, -60, -50, -46, -42, -38, -34] + list(range(-30, 30, 1)) + [30, 34, 38, 42, 46, 50, 60, 70, 80, 100, 300]


    ###### STEP 1: floating fits
    iteration = 0
    outDir_fits_v0 = "%s/fits_v0" % baseDir
    fitCfg = {}
    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0,0], [1, 0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [35, 5, 8, 20], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.005, 0.1, 0.6], [1, 0, 0]

    fitCfg['mean'], fitCfg['mean_cfg'] = [0, 0, 0, 0], [0, 0, 0, 0] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
    fitCfg['sigma'], fitCfg['sigma_cfg'] = [30, 5, 9, 14], [0, 0, 0, 0]
    fitCfg['norm'], fitCfg['norm_cfg'] = [0.005, 0.30, 0.55], [1, 0, 0]

    bkgCfg = {}
    bkgCfg['TTbar'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/ttbar_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['EWK'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/ewk_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    bkgCfg['ZZ'] = { "filename": "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_mumu_DeepMETReso/zz_para/fits_v0/results.json", "norm": 1.0, "float": False, "subtract": False }
    nGauss = len(fitCfg['mean'])
    if True: 
        recoilLibs.doFitMultiGauss(rFile, proc, comp, fitCfg, label, outDir_fits_v0, recoil_qTbins, qTmax=500, rebin=rebin, bkgCfg=bkgCfg, singleMean=False, xMin=-150, xMax=150)
        return

    
    ###### STEP 2: parameterize mean
    outDir_param_v0 = "%s/param_v0" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v0, True)
        with open("%s/results.json" % outDir_fits_v0) as f: jsIn = json.load(f)

        # mean
        #fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(    ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x)", [1.54013e-01, 4.20263e-01, -1.38828e-02, 2.80192e-04, 9.44310e-02], [False, False, False, False, False, False]
        #fitF, params, cParams = "(x<[4])*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>[4])*1*(    ([0] + [1]*[4] + [2]*[4]*[4] + [3]*[4]*[4]*[4] - ([1] + 2*[2]*[4] + 3*[3]*[4]*[4])*[4]) + ([1] + 2*[2]*[4] + 3*[3]*[4]*[4])*x)", [1.54013e-01, 4.20263e-01, -1.38828e-02, 2.80192e-04, 30], [False, False, False, False, False, False]
        #fitF, params, cParams = "(x<36)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>36)*1*(    ([0] + [1]*36 + [2]*36*36 + [3]*36*36*36 - ([1] + 2*[2]*36 + 3*[3]*36*36)*36) + ([1] + 2*[2]*36 + 3*[3]*36*36)*x  )", [0, 4.20263e-01, -1.38828e-02, 2.80192e-04], [True, False, False, False]
        #fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(    ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - ([1] + 2*[2]*30 + 3*[3]*30*30 - 2*[4]*30)*30 - [4]*30*30) + ([1] + 2*[2]*30 + 3*[3]*30*30 - 2*[4]*30)*x + [4]*x*x  )", [1.58212e-01, 4.03838e-01, -1.03807e-02, 1.21562e-04, -2.26096e-04], [False, False, False, False, False]
        #fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(    ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - ([1] + 2*[2]*30 + 3*[3]*30*30)*30) + ([1] + 2*[2]*30 + 3*[3]*30*30)*x )", [1.58212e-01, 4.03838e-01, -1.03807e-02, 1.21562e-04], [False, False, False, False]

        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(  ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x   )", [2.45501e-01, 3.67274e-01, -7.26948e-03, 5.62285e-05, 8.95248e-02], [False, False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v0, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)

        with open("%s/results.json" % outDir_param_v0, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        sys.exit()
        
    ###### STEP 3: refit with fixed means
    outDir_fits_v1 = "%s/fits_v1" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0]
        fitCfg['norm_cfg'] = [0]

        with open("%s/results.json" % outDir_param_v0) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_fits_v1, funcJs=jsIn, qTmax=150, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, orderSigmas=True)
        sys.exit()

    ###### STEP 4: parameterize norm1
    outDir_param_v1 = "%s/param_v1" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v1, True)
        with open("%s/results.json" % outDir_fits_v1) as f: jsIn = json.load(f)
        
        # means
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(  ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - [4]*30) + [4]*x   )", [2.45501e-01, 3.67274e-01, -7.26948e-03, 5.62285e-05, 8.95248e-02], [False, False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        # sigma1
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 4.91301e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=30, yTitle = "#sigma_{1} (GeV)", doFit=True)
        
        # sigma2
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [2.44857e-01, 5.19971e-01, 11], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=30, yTitle = "#sigma_{2} (GeV)", doFit=True)
        
        # norm
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.72, 0.015], [True, True] # ad-hoc exp
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v1, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)

        with open("%s/results.json" % outDir_param_v1, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        sys.exit()
    
    ###### STEP 5: refit with fixed means and fixed norm1
    outDir_fits_v2 = "%s/fits_v2" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 0]
        fitCfg['norm_cfg'] = [2]

        with open("%s/results.json" % outDir_param_v1) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_fits_v2, funcJs=jsIn, qTmax=150, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, orderSigmas=False)
        sys.exit()
    
    ###### STEP 5: parameterize sigma2
    outDir_param_v2 = "%s/param_v2" % baseDir
    if False:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v2, True)
        with open("%s/results.json" % outDir_fits_v2) as f: jsIn = json.load(f)

        # means
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(    ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - ([1] + 2*[2]*30 + 3*[3]*30*30)*30) + ([1] + 2*[2]*30 + 3*[3]*30*30)*x )", [1.58212e-01, 4.03838e-01, -1.03807e-02, 1.21562e-04], [False, False, False, False]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v2, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v2, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        # sigma1
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [3.95805e-01, 4.45398e-01, 5.94177e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v2, label, cParams=cParams, fitMin=0, fitMax=80, yMin=0, yMax=30, yTitle = "#sigma_{1} (GeV)")
        
        # sigma2
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [3.06666e-02, 1.07437e+00, 1.11954e+01], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v2, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=30, yTitle = "#sigma_{2} (GeV)")
        
        # norm
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.72, 0.015], [True, True] # ad-hoc exp
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v2, label, cParams=cParams, fitMin=0, fitMax=100, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=True)

        with open("%s/results.json" % outDir_param_v2, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        sys.exit()



    ###### STEP 5: refit with fixed means, sigma2 and norm1
    outDir_fits_v3 = "%s/fits_v3" % baseDir
    if False:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [0, 2]
        fitCfg['norm_cfg'] = [2]

        with open("%s/results.json" % outDir_param_v2) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_fits_v3, funcJs=jsIn, qTmax=150, rebin=rebin, bkgCfg=bkgCfg, singleMean=True)
        sys.exit()
        
    ###### STEP 5: parameterize sigma2 and sigma3
    outDir_param_v3 = "%s/param_v3" % baseDir
    if True:
        jsOut = {}
        jsOut['nGauss'] = nGauss
        functions.prepareDir(outDir_param_v3, True)
        with open("%s/results.json" % outDir_fits_v3) as f: jsIn = json.load(f)

        # means
        fitF, params, cParams = "(x<30)*1*([0] + [1]*x + [2]*x*x + [3]*x*x*x) +0+ (x>30)*1*(    ([0] + [1]*30 + [2]*30*30 + [3]*30*30*30 - ([1] + 2*[2]*30 + 3*[3]*30*30)*30) + ([1] + 2*[2]*30 + 3*[3]*30*30)*x )", [3.08663e-01, 3.48668e-01, -6.06558e-03, 3.62872e-05], [False, False, True, True]
        doFit(jsIn, jsOut, comp, "mean1", fitF, params, outDir_param_v3, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, yTitle = "#mu_{1} (GeV)", doFit=True)
        doFit(jsIn, jsOut, comp, "mean2", fitF, params, outDir_param_v3, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=50, yTitle = "#mu_{2} (GeV)", doFit=True)
       
        # sigma1
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [3.15325e-01, 5.07352e-01, 6.04702e+00], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma1", fitF, params, outDir_param_v3, label, cParams=cParams, fitMin=0, fitMax=70, yMin=0, yMax=30, yTitle = "#sigma_{1} (GeV)", doFit=True)
        
        # sigma2
        fitF, params, cParams = "[0]*TMath::Power(x, [1]) + [2]", [3.04763e-02, 1.07898e+00, 1.11939e+01], [False, False, False]
        doFit(jsIn, jsOut, comp, "sigma2", fitF, params, outDir_param_v3, label, cParams=cParams, fitMin=0, fitMax=140, yMin=0, yMax=30, yTitle = "#sigma_{2} (GeV)")
        
        # norm
        fitF, params, cParams = "[0]*TMath::Exp(-[1]*x)", [0.72, 0.015], [True, True] # ad-hoc exp
        doFit(jsIn, jsOut, comp, "norm1", fitF, params, outDir_param_v3, label, cParams=cParams, fitMin=0, fitMax=120, yMin=0, yMax=1, yTitle = "n_{1} (GeV)", doFit=False)

        with open("%s/results.json" % outDir_param_v3, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
    
    ###### STEP 5: refit
    if True:
        with open("%s/results.json" % outDir_param_v3) as f: jsIn = json.load(f)
        jsOut = combinedFit(proc, comp, jsIn, rebin=rebin, bkgCfg=bkgCfg, singleMean=True, qTmax=70) # 95 was good
        with open("%s/results_refit.json" % outDir_param_v3, "w") as outfile: json.dump(jsOut, outfile, indent=4)
        #sys.exit()
   
    ###### STEP 6: final validation
    if True:
        fitCfg['mean_cfg'] = [2, 2] # 0 = floating, 1 = fixed, 2 = parametric (fixed), 3 = parametric (floating)
        fitCfg['sigma_cfg'] = [2, 2]
        fitCfg['norm_cfg'] = [2]

        outDir_refit = "%s/refit_v3" % baseDir
        with open("%s/results_refit.json" % outDir_param_v3) as f: jsIn = json.load(f)
        doFitMultiGauss(proc, comp, fitCfg, label, outDir_refit, funcJs=jsIn, qTmax=150, rebin=rebin, bkgCfg=bkgCfg, singleMean=True)
 
        outDir_refit = "%s/param_v3_refit" % baseDir
        functions.prepareDir(outDir_refit)
        doPlot(comp, "norm1", jsIn, outDir_refit, label, yMin=0, yMax=1, yTitle = "n_{1} (GeV)")
        doPlot(comp, "sigma1", jsIn, outDir_refit, label, yMin=0, yMax=30, yTitle = "#sigma_{1} (GeV)")
        doPlot(comp, "sigma2", jsIn, outDir_refit, label, yMin=0, yMax=30, yTitle = "#sigma_{2} (GeV)")
        doPlot(comp, "mean1", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#mu_{1} (GeV)")
        doPlot(comp, "mean2", jsIn, outDir_refit, label, yMin=0, yMax=20, yTitle = "#mu_{2} (GeV)")



def prepareFile(fInName, fOutName):

    def readProc(groups_mumu, hName, procName):

        if isinstance(procName, str):
            label = "%s_%s" % (hName, procName) 
            groups_mumu.setHists(hName, "", label=label, procsToRead=[procName], selectSignal=False)
            bhist = groups_mumu.groups[procName][label]
        else:
            label = "%s_%s" % (hName, procName[0]) 
            groups_mumu.setHists(hName, "", label=label, procsToRead=procName, selectSignal=False)
            bhist = groups_mumu.groups[procName[0]][label]
        return bhist 


    datagroups = datagroups2016(fInName)
    datagroups.groups.update({
        "TTbar" : dict(
                members = [datagroups.datasets[x] for x in ["TTLeptonicPostVFP", "TTSemileptonicPostVFP", "SingleTschanLepDecaysPostVFP", "SingleTtWAntitopPostVFP", "SingleTtchanAntitopPostVFP", "SingleTtchanTopPostVFP"]],
                label = r"Z$\to\tau\tau$",
                color = "darkblue",
        ),
        "EWK" : dict(
                members = [datagroups.datasets[x] for x in ["WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "WZPostVFP", "WWPostVFP"]],
                #members = [datagroups.datasets[x] for x in ["WplusmunuPostVFP", "WminusmunuPostVFP", "WplustaunuPostVFP", "WminustaunuPostVFP", "ZtautauPostVFP", "ZZ2l2nuPostVFP", ]],
                label = "EWK", 
                color = "darkblue",
        ),
        "ZZ" : dict(
                members = [datagroups.datasets[x] for x in ["ZZ2l2nuPostVFP"]],
                label = "ZZ",
                color = "darkblue",
        )
    })
    
    fOut = ROOT.TFile(fOutName, "RECREATE")


    procs = ["Data", "Zmumu", "EWK", "TTbar", "ZZ"]
    for proc in procs:
        print(proc)
        bhist_para_qT = readProc(datagroups, "recoil_corr_xy_para_qT_qTbinned", proc) # recoil_corr_xy_para_qT_qTbinned recoil_uncorr_para_qT_qTbinned
        bhist_perp = readProc(datagroups, "recoil_corr_xy_perp_qTbinned", proc) # recoil_corr_xy_perp_qTbinned recoil_uncorr_perp_qTbinned
        rhist_para_qT = narf.hist_to_root(bhist_para_qT)
        rhist_perp = narf.hist_to_root(bhist_perp)
        
        rhist_para_qT.SetName("%s_para" % (proc))
        rhist_perp.SetName("%s_perp" % (proc))
        
        rhist_para_qT.Write()
        rhist_perp.Write()
        
        '''
        for iBin in range(1, rhist_para_qT.GetNbinsX()+2):

            hist_para = rhist_para_qT.ProjectionY("%s_para_bin%d" % (proc, iBin), iBin, iBin)
            hist_para.Write()
            print(iBin, hist_para.Integral())
            
            hist_perp = rhist_perp.ProjectionY("%s_perp_bin%d" % (proc, iBin), iBin, iBin)
            hist_perp.Write()
        '''
        
        
    '''
    # GEN level qT 
    if flavor == "mumu": proc = "DYmumu"
    if flavor == "ee": proc = "DYee"
    if flavor == "mu": proc = ["WplusJetsToMuNu", "WminusJetsToMuNu"] # merge plus and minus
    if flavor == "e": proc = ["WplusJetsToENu", "WminusJetsToENu"] # merge plus and minus
    bhist_para_qT = readProc(datagroups, "recoil_corr_xy_para_qT_qTbinned_gen", proc)
    bhist_perp = readProc(datagroups, "recoil_corr_xy_perp_qTbinned_gen", proc)
    rhist_para_qT = narf.hist_to_root(bhist_para_qT)
    rhist_perp = narf.hist_to_root(bhist_perp)
    for iBin in range(1, rhist_para_qT.GetNbinsX()+1):
            
        hist_para = rhist_para_qT.ProjectionY("%s_gen_para_bin%d" % (sig, iBin), iBin, iBin)
        print(hist_para)
        hist_para.Write()
            
        hist_perp = rhist_perp.ProjectionY("%s_gen_perp_bin%d" % (sig, iBin), iBin, iBin)
        hist_perp.Write()
    '''    
    fOut.ls()
    fOut.Close()
     
if __name__ == "__main__":

    met = "DeepMETReso" # DeepMETReso RawPFMET
    flavor = "mumu" # mu, e, mumu, ee
    
  
    rFile = "wremnants/data/lowPU/recoil/highPU_%s_%s.root" % (flavor, met)
    #prepareFile("mz_wlike_with_mu_eta_pt_%s.pkl.lz4" % met, rFile)
    #sys.exit()
    
    #recoil_qTbins = list(range(0, 300, 1)) + list(range(300, 350, 2))  + list(range(350, 400, 5)) + list(range(400, 520, 20))
    #recoil_qTbins = list(range(0, 150, 1)) + list(range(150, 200, 2)) + list(range(200, 250, 5)) + list(range(250, 300, 10)) + list(range(300, 500, 50)) + [500]
    recoil_qTbins = list(functions.drange(0, 30, 0.5)) + list(range(30, 150, 1)) + list(range(150, 200, 2)) + list(range(200, 250, 5)) + list(range(250, 300, 10)) + list(range(300, 500, 50)) + [500]
    #recoil_qTbins = list(range(300, 350, 2))  + list(range(350, 400, 5)) + list(range(400, 520, 20))
    #recoil_qTbins = list(functions.drange(0, 30, 0.5))
    outDir = "/eos/user/j/jaeyserm/www/wmass/highPU/recoilCorrection/recoil_%s_%s/" % (flavor, met)
    functions.prepareDir(outDir, False)

    if met == "RawPFMET":
    
        #do_dymumu_para_RawPFMET()
        do_dymumu_perp_RawPFMET()
        
        #do_ttbar_para_RawPFMET()
        #do_ttbar_perp_RawPFMET()
        
        #do_ewk_para_RawPFMET()
        #do_ewk_perp_RawPFMET()
        
        #do_zz_para_RawPFMET()
        #do_zz_perp_RawPFMET()
        
        #do_singlemuon_para_RawPFMET()
        #do_singlemuon_perp_RawPFMET()
        pass

    if met == "DeepMETReso":
        
        do_dymumu_para_DeepMETReso()
        #do_dymumu_perp_DeepMETReso()
        
        #do_ttbar_para_DeepMETReso()
        #do_ttbar_perp_DeepMETReso()
        
        #do_ewk_para_DeepMETReso()
        #do_ewk_perp_DeepMETReso()
        
        #do_zz_para_DeepMETReso()
        #do_zz_perp_DeepMETReso()
        
        #do_singlemuon_para_DeepMETReso()
        #do_singlemuon_perp_DeepMETReso()