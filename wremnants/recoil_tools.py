import ROOT
import hist
import numpy as np
import copy
import logging
import sys
import decimal
import json
import os
import array
from utilities import common as common

ROOT.gInterpreter.Declare('#include "lowpu_recoil.h"')
logger = logging.getLogger("wremnants").getChild(__name__.split(".")[-1])


def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)
        
        
class Recoil:

    def __init__(self, type_, args, flavor="mu"):

        self.met = args.met
        self.flavor = flavor 
        self.parametric = True ## remove
        self.highPU = False
        self.args = args
        self.smearWeights = True
        self.storeHists = args.recoilHists
        setattr(ROOT.wrem, "recoil_verbose", True if args.verbose > 3 else False)

        if type_ == "highPU":
            self.highPU = True
            setattr(ROOT.wrem, "recoil_correction_qTmax", 1000)
            self.recoil_qTbins = list(drange(0, 500, 1)) + [500]
            
            self.met_xycorr_setup(f"{common.data_dir}/recoil/highPU/{self.flavor}_{self.met}/met_xy_correction.json") # MET XY correction
            if "mu" in self.flavor:
                self.set_qT_weights(f"{common.data_dir}/recoil/highPU/mumu_{self.met}/qT_reweighting.json") # qT reweigthing
            
            # recoil calibrations
            if self.met == "RawPFMET":
                logger.info(f"Apply recoil corrections for {self.met}")
                flavor_ = "mumu" # both mu and mumu
                self.addParametric("target_para", f"{common.data_dir}/recoil/highPU/{flavor_}_{self.met}/recoil_data_para.json", doUnc=False)
                self.addParametric("source_para", f"{common.data_dir}/recoil/highPU/{flavor_}_{self.met}/recoil_zmumu_para.json", doUnc=False)
                self.addParametric("target_perp", f"{common.data_dir}/recoil/highPU/{flavor_}_{self.met}/recoil_data_perp.json", doUnc=False)
                self.addParametric("source_perp", f"{common.data_dir}/recoil/highPU/{flavor_}_{self.met}/recoil_zmumu_perp.json", doUnc=False)
                setattr(ROOT.wrem, "applyRecoilCorrection", True)
                
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ttbarUp.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ttbarDown.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ewkUp.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ewkDown.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ttbarUp.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ttbarDown.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ewkUp.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ewkDown.json")
                
            else:
                logger.warning(f"Recoil corrections for {self.met} not available, use default XY-corrected MET")
                setattr(ROOT.wrem, "applyRecoilCorrection", False)
    
        elif type_ == "lowPU":
            self.highPU = False
            setattr(ROOT.wrem, "recoil_correction_qTmax", 9e99)
            self.recoil_qTbins = list(drange(0, 150, 0.5)) + [150]
            
            self.met_xycorr_setup(f"{common.data_dir}/recoil/lowPU/{self.flavor}_{self.met}/met_xy_correction.json") # MET XY correction
            if "mu" in self.flavor: # todo ee
                self.set_qT_weights(f"{common.data_dir}/recoil/lowPU/mumu_{self.met}/qT_reweighting.json") # qT reweigthing

            # recoil calibrations
            if self.met == "RawPFMET":
                logger.info(f"Apply recoil corrections for {self.met}")
                flavor_ = "mumu" # both mu and mumu
                self.addParametric("target_para", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para.json")
                self.addParametric("source_para", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_zmumu_para.json")
                self.addParametric("target_perp", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp.json")
                self.addParametric("source_perp", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_zmumu_perp.json")
                setattr(ROOT.wrem, "applyRecoilCorrection", True)
                
                # syst variations should contain target/source and para/perp, Needed for weights
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ttbarUp.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ttbarDown.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ewkUp.json")
                self.addParametricUnc("target_para_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_para_ewkDown.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ttbarUp.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ttbarDown.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ewkUp.json")
                self.addParametricUnc("target_perp_bkg", f"{common.data_dir}/recoil/lowPU/{flavor_}_{self.met}/recoil_data_perp_ewkDown.json")
                
            else:
                logger.warning(f"Recoil corrections for {self.met} not available, use default XY-corrected MET")
                setattr(ROOT.wrem, "applyRecoilCorrection", False)
                
        else: sys.exit("Recoil highPU or lowPU")
        
        
        
        # define axes
        if self.flavor == "mu": self.axis_MET_pt = hist.axis.Regular(200, 0, 200, name = "recoil_MET_pt", underflow=False)
        else: self.axis_MET_pt = hist.axis.Regular(400, 0, 200, name = "recoil_MET_pt", underflow=False)
        self.axis_MET_phi = hist.axis.Regular(50, -4, 4, name = "recoil_MET_phi")
        self.axis_MET_dphi = hist.axis.Regular(41, -4.1, 4.1, name = "recoil_MET_phi")
        
       
        self.axis_recoil_magn = hist.axis.Regular(300, 0, 300, name = "recoil_magn", underflow=False)
        self.axis_recoil_para_qT = hist.axis.Regular(400, -300, 100, name = "recoil_para_qT")
        self.axis_recoil_para = hist.axis.Variable(list(range(-150, 150, 2)) + [150], name = "recoil_para")
        self.axis_recoil_perp = hist.axis.Variable(list(range(-150, 150, 2)) + [150], name = "recoil_perp")
        
        self.axis_recoil_magn_fine = hist.axis.Regular(300, 0, 300, name = "recoil_magn", underflow=False)
        self.axis_recoil_para_fine = hist.axis.Regular(1000, -500, 500, name = "recoil_para")
        self.axis_recoil_perp_fine = hist.axis.Regular(1000, -500, 500, name = "recoil_perp")
        self.axis_MET_xy = hist.axis.Regular(200, -100, 100, name = "MET_xy")

        self.axis_mt = hist.axis.Regular(200, 0., 200., name = "mt", underflow=False)
        self.axis_njets = hist.axis.Regular(30, 0.5, 30.5, name = "recoil_njets")
        self.axis_npv = hist.axis.Regular(100, 0.5, 100.5, name = "recoil_npv")
        self.axis_sumEt = hist.axis.Regular(400, 0, 4000, name = "recoil_sumEt")
        self.axis_rapidity = hist.axis.Regular(24, -2.4, 2.4, name = "recoil_rapidity")

        self.axis_qT = hist.axis.Variable(self.recoil_qTbins, name = "qTbinned", underflow=False)
        
        self.axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
        self.axis_passIso = hist.axis.Boolean(name = "passIso")
        self.axis_passMT = hist.axis.Boolean(name = "passMT")
        
        
        
        self.axis_recoil_2d_para = hist.axis.Regular(40, -10, 10, name = "axis_recoil_2d_para", underflow=False, overflow=False)
        self.axis_recoil_2d_perp = hist.axis.Regular(40, -10, 10, name = "axis_recoil_2d_perp", underflow=False, overflow=False)


    def set_qT_weights(self, fIn):
        if not os.path.exists(fIn):
            logger.warning(f"qT weights file {fIn} not found")
            return
        else: 
            logger.info(f"Add qT weights {fIn}")
        with open(fIn) as f: jsIn = json.load(f)
        qTrw_bins = getattr(ROOT.wrem, "qTrw_bins")
        qTrw_weights = getattr(ROOT.wrem, "qTrw_weights")
        for qTbin in jsIn['qT_bins']:
            qTrw_bins.push_back(float(qTbin))
        for w in jsIn['weights']:
            qTrw_weights.push_back(w)
        setattr(ROOT.wrem, "applyqTReweighting", True)


    def addParametric(self, tag, fIn, doUnc=True):

        if not os.path.exists(fIn):
            logger.warning(f"Recoil parametric file {fIn} not found")
            return
        else: 
            logger.info(f"Add recoil parametric set for {tag}")
        with open(fIn) as f: jsIn = json.load(f)
        nGauss = jsIn['nGauss']
        recoil_param_nGauss = getattr(ROOT.wrem, "recoil_param_nGauss")
        recoil_param_funcs = getattr(ROOT.wrem, "recoil_param_funcs")
        recoil_param_nGauss.insert((tag, nGauss))
        
        for iGauss in range(1, nGauss+1):
            
            label = "mean%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s" % (tag, label)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)            
            
            label = "sigma%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s" % (tag, label)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)
            
            if iGauss == nGauss: continue
            label = "norm%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s" % (tag, label)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)
            
        
        # do statistical variations
        if not doUnc: return
        nStatVars = jsIn['nStatVars']
        recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
        recoil_param_nStat.insert((tag, nStatVars*2))
        iiStat = 0
        test = {}
        for nStat in range(0, nStatVars):
            for statVar in ["p", "m"]:
                systLabel = "syst%d" % iiStat
                statLabel = "stat%d_%s" % (nStat, statVar)
                for iGauss in range(1, nGauss+1):
                    
                    label = "mean%d" % iGauss
                    pars = [jsIn[statLabel][label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
                    funcName = "%s_%s_%s" % (tag, label, systLabel)
                    funcExpr = jsIn[label]['func']
                    ROOT.wrem.insertFunction(funcName, funcExpr, pars)

                    label = "sigma%d" % iGauss
                    pars = [jsIn[statLabel][label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
                    funcName = "%s_%s_%s" % (tag, label, systLabel)
                    funcExpr = jsIn[label]['func']
                    ROOT.wrem.insertFunction(funcName, funcExpr, pars)

                    if iGauss == nGauss: continue
                    label = "norm%d" % iGauss
                    pars = [jsIn[statLabel][label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
                    funcName = "%s_%s_%s" % (tag, label, systLabel)
                    funcExpr = jsIn[label]['func']
                    ROOT.wrem.insertFunction(funcName, funcExpr, pars)
                iiStat += 1

    def addParametricUnc(self, tag, fIn):

        logger.info(f"Add recoil parametric uncertainty set for {tag}")
        with open(fIn) as f: jsIn = json.load(f)
        nGauss = jsIn['nGauss']
        recoil_param_nGauss = getattr(ROOT.wrem, "recoil_param_nGauss")
        recoil_param_funcs = getattr(ROOT.wrem, "recoil_param_funcs")
        recoil_param_nGauss.insert((tag, nGauss))
        
        recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
        if not tag in recoil_param_nStat:
            recoil_param_nStat.insert((tag, 0))
  
        nStatVars = recoil_param_nStat[tag] + 1
        recoil_param_nStat[tag] = nStatVars
        nStat = nStatVars-1
        systLabel = "syst%d" % nStat
        for iGauss in range(1, nGauss+1):
            
            label = "mean%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s_%s" % (tag, label, systLabel)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)
            
            label = "sigma%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s_%s" % (tag, label, systLabel)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)
            
            if iGauss == nGauss: continue
            label = "norm%d" % iGauss
            pars = [jsIn[label]['p%d' % iParam] for iParam in range(0, jsIn[label]['nParams'])]
            funcName = "%s_%s_%s" % (tag, label, systLabel)
            funcExpr = jsIn[label]['func']
            ROOT.wrem.insertFunction(funcName, funcExpr, pars)
            
        
        
            

    def addBinned(self, tag, fIn):
    
        print("Add recoil binned set for %s" % tag)
        with open(fIn) as f: jsIn = json.load(f)
        nGauss = jsIn['nGauss']
        qTbins = jsIn['qTbins']
        recoil_binned_nGauss = getattr(ROOT.wrem, "recoil_binned_nGauss")
        recoil_binned_qTbins = getattr(ROOT.wrem, "recoil_binned_qTbins")
        recoil_binned_components = getattr(ROOT.wrem, "recoil_binned_components")
        recoil_binned_nGauss.insert((tag, nGauss))
        
        vec = ROOT.std.vector("float")()
        for qTbin in qTbins: 
            vec.push_back(1.0*qTbin)
        recoil_binned_qTbins.insert((tag, vec))

        vec_qt = ROOT.std.vector("std::vector<std::vector<float>>")()
        for iBin in range(1, len(qTbins)):
            dd = jsIn[str(iBin)]
            n = 0
            vec_gauss = ROOT.std.vector("std::vector<float>")()
            for iGauss in range(1, nGauss+1):
                vec_comp = ROOT.std.vector("float")()
                vec_comp.push_back(dd["mean%d" % iGauss])
                vec_comp.push_back(dd["sigma%d" % iGauss])
                if iGauss == nGauss: vec_comp.push_back(1.-n)
                else:                
                    vec_comp.push_back(dd["norm%d" % iGauss])
                    n += dd["norm%d" % iGauss]
                vec_gauss.push_back(vec_comp)
            vec_qt.push_back(vec_gauss)
        recoil_binned_components.insert((tag, vec_qt))
        
        


    def met_xycorr_setup(self, fIn):
    
        def loadCoeff(d, out):
            npol = d['polyOrder']
            coeff = ROOT.std.vector["float"]()
            for v in range(0, npol+1): coeff.push_back(d['p%d' % v])
            setattr(ROOT.wrem, out, coeff)

        # load the polynomial coefficients
        if not os.path.exists(fIn):
            logger.warning(f"XY correction file {fIn} not found")
            return
        else: 
            logger.info(f"Add MET XY correction {fIn}")
        with open(fIn) as f: jsdata = json.load(f)
        loadCoeff(jsdata['x']['data']['nominal'], 'met_xy_corr_x_data_nom')
        loadCoeff(jsdata['y']['data']['nominal'], 'met_xy_corr_y_data_nom')
        loadCoeff(jsdata['x']['mc']['nominal'], 'met_xy_corr_x_mc_nom')
        loadCoeff(jsdata['y']['mc']['nominal'], 'met_xy_corr_y_mc_nom')
        setattr(ROOT.wrem, "applyMETxyCorrection", True)


    def recoil_Z(self, df, results, dataset, datasets_to_apply, lep_cols, trg_cols):
    
        self.df = df
        self.results = results
        self.dataset = dataset
        self.datasets_to_apply = datasets_to_apply
        
        self.leptons_pt = lep_cols[0]
        self.leptons_phi = lep_cols[1]
        self.leptons_uncorr_pt = lep_cols[2]
        
        self.trgLep_pt = trg_cols[0]
        self.trgLep_phi = trg_cols[1]
        self.nonTrgLep_pt = trg_cols[2]
        self.nonTrgLep_phi = trg_cols[3]
        
        self.setup_MET()
        self.setup_recoil_Z()
        self.setup_recoil_gen()
        self.apply_recoil_Z()
        if self.storeHists:
            self.auxHists()
        
        
        return self.df
   
    def recoil_W(self, df, results, dataset, datasets_to_apply, lep_cols):
    
        self.df = df
        self.results = results
        self.dataset = dataset
        self.datasets_to_apply = datasets_to_apply
        
        self.leptons_pt = lep_cols[0]
        self.leptons_phi = lep_cols[1]
        self.leptons_charge = lep_cols[2]
        self.leptons_uncorr_pt = lep_cols[3]

        self.setup_MET()
        self.setup_recoil_gen()
        self.apply_recoil_W()

        return self.df
   
        
    def setup_MET(self):
        '''
        Define MET variables, apply lepton and XY corrections
        '''
        
        # for the Z, leptons_pt, leptons_phi, leptons_uncorr_pt should be Vec_f with size 2
        # for the W, it should contain only a float/double

    
        if self.met == "PFMET": met_pt, met_phi = "MET_pt", "MET_phi"
        elif self.met == "RawPFMET": met_pt, met_phi = "RawMET_pt", "RawMET_phi"
        elif self.met == "DeepMETReso": met_pt, met_phi = "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi"
        else: sys.exit("MET type %s not supported" % self.met)
     
        # uncorrected MET
        self.df = self.df.Alias("MET_uncorr_pt", met_pt)
        self.df = self.df.Alias("MET_uncorr_phi", met_phi)
        self.df = self.df.Define("METx_uncorr", "MET_uncorr_pt*cos(MET_uncorr_phi)")
        self.df = self.df.Define("METy_uncorr", "MET_uncorr_pt*sin(MET_uncorr_phi)")

        # lepton corrected MET
        self.df = self.df.Define("MET_corr_lep", "wrem::METLeptonCorrection(MET_uncorr_pt, MET_uncorr_phi, %s, %s, %s)" % (self.leptons_pt, self.leptons_uncorr_pt, self.leptons_phi))
        self.df = self.df.Define("MET_corr_lep_pt", "MET_corr_lep[0]")
        self.df = self.df.Define("MET_corr_lep_phi", "MET_corr_lep[1]")
        self.df = self.df.Define("METx_corr_lep", "MET_corr_lep_pt*cos(MET_corr_lep_phi)")
        self.df = self.df.Define("METy_corr_lep", "MET_corr_lep_pt*sin(MET_corr_lep_phi)")

        # phi corrected MET (XY corrections)
        self.df = self.df.Define("MET_corr_xy", "wrem::METXYCorrection(MET_corr_lep_pt, MET_corr_lep_phi, PV_npvs, %d)" % self.dataset.is_data)
        self.df = self.df.Define("MET_corr_xy_pt", "MET_corr_xy[0]")
        self.df = self.df.Define("MET_corr_xy_phi", "MET_corr_xy[1]")
        self.df = self.df.Define("METx_corr_xy", "MET_corr_xy_pt*cos(MET_corr_xy_phi)")
        self.df = self.df.Define("METy_corr_xy", "MET_corr_xy_pt*sin(MET_corr_xy_phi)")

        if not self.storeHists: 
            return
            
        self.results.append(self.df.HistoBoost("MET_uncorr_pt", [self.axis_MET_pt], ["MET_uncorr_pt", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_uncorr_phi", [self.axis_MET_phi], ["MET_uncorr_phi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METx_uncorr", [self.axis_MET_xy], ["METx_uncorr", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_uncorr", [self.axis_MET_xy], ["METy_uncorr", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("MET_corr_lep_pt", [self.axis_MET_pt], ["MET_corr_lep_pt", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_lep_phi", [self.axis_MET_phi], ["MET_corr_lep_phi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METx_corr_lep", [self.axis_MET_xy], ["METx_corr_lep", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_corr_lep", [self.axis_MET_xy], ["METy_corr_lep", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("MET_corr_xy_pt", [self.axis_MET_pt], ["MET_corr_xy_pt", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_xy_phi", [self.axis_MET_phi], ["MET_corr_xy_phi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METx_corr_xy", [self.axis_MET_xy], ["METx_corr_xy", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_corr_xy", [self.axis_MET_xy], ["METy_corr_xy", "nominal_weight"]))
        

        # histograms as function of npv, to derive/closure test the XY correction
        self.results.append(self.df.HistoBoost("METx_corr_lep_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METx_corr_lep", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_corr_lep_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METy_corr_lep", "nominal_weight"]))
            
        self.results.append(self.df.HistoBoost("METx_corr_xy_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METx_corr_xy", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_corr_xy_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METy_corr_xy", "nominal_weight"]))
       
    
   

    def setup_recoil_Z(self):
        '''
        Setup the uncorrected recoil components
        '''
        
        # construct the 2D vector sum of the two leptons
        self.df = self.df.Define("Lep1_mom2", "ROOT::Math::Polar2DVectorD(%s[0], %s[0])" % (self.leptons_pt, self.leptons_phi))
        self.df = self.df.Define("Lep2_mom2", "ROOT::Math::Polar2DVectorD(%s[1], %s[1])" % (self.leptons_pt, self.leptons_phi))
        self.df = self.df.Define("Z_mom2", "Lep1_mom2 + Lep2_mom2")
        self.df = self.df.Define("qT", "Z_mom2.R()")
        self.df = self.df.Define("qTbin", "wrem::getqTbin(qT)")
              
        if self.dataset.name in self.datasets_to_apply:
            self.df = self.df.Define("qTweight_", "wrem::qTweight(qT)")
            if not self.args.theoryCorrAltOnly and self.args.theoryCorr:
                self.df = self.df.Define("nominal_weight_qTrw", "nominal_weight_uncorr*qTweight_")
            else:
                self.df = self.df.Define("nominal_weight_qTrw", "nominal_weight*qTweight_")

        else:
            self.df = self.df.Alias("nominal_weight_qTrw", "nominal_weight")
        
        # uncorrected recoil
        self.df = self.df.Define("recoil_uncorr", "wrem::recoilComponents(MET_uncorr_pt, MET_uncorr_phi, qT, Z_mom2.Phi())")
        self.df = self.df.Define("recoil_uncorr_magn", "recoil_uncorr[0]")
        self.df = self.df.Define("recoil_uncorr_para_qT", "recoil_uncorr[1]")
        self.df = self.df.Define("recoil_uncorr_para", "recoil_uncorr[1] + qT")
        self.df = self.df.Define("recoil_uncorr_perp", "recoil_uncorr[2]")
        
        # lep corrected recoil
        self.df = self.df.Define("recoil_corr_lep", "wrem::recoilComponents(MET_corr_lep_pt, MET_corr_lep_phi, qT, Z_mom2.Phi())")
        self.df = self.df.Define("recoil_corr_lep_magn", "recoil_corr_lep[0]")
        self.df = self.df.Define("recoil_corr_lep_para_qT", "recoil_corr_lep[1]")
        self.df = self.df.Define("recoil_corr_lep_para", "recoil_corr_lep[1] + qT")
        self.df = self.df.Define("recoil_corr_lep_perp", "recoil_corr_lep[2]")
        
        # MET XY corrected recoil
        self.df = self.df.Define("recoil_corr_xy", "wrem::recoilComponents(MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())")
        self.df = self.df.Define("recoil_corr_xy_magn", "recoil_corr_xy[0]")
        self.df = self.df.Define("recoil_corr_xy_para_qT", "recoil_corr_xy[1]")
        self.df = self.df.Define("recoil_corr_xy_para", "recoil_corr_xy[1] + qT")
        self.df = self.df.Define("recoil_corr_xy_perp", "recoil_corr_xy[2]")
        
        # transverse mass
        self.df = self.df.Define("mT_uncorr", f"wrem::mt_wlike_nano({self.trgLep_pt}, {self.trgLep_phi}, {self.nonTrgLep_pt}, {self.nonTrgLep_phi}, MET_uncorr_pt, MET_uncorr_phi)")
        self.df = self.df.Define("mT_corr_lep", f"wrem::mt_wlike_nano({self.trgLep_pt}, {self.trgLep_phi}, {self.nonTrgLep_pt}, {self.nonTrgLep_phi}, MET_corr_lep_pt, MET_corr_lep_phi)")
        self.df = self.df.Define("mT_corr_xy", f"wrem::mt_wlike_nano({self.trgLep_pt}, {self.trgLep_phi}, {self.nonTrgLep_pt}, {self.nonTrgLep_phi}, MET_corr_xy_pt, MET_corr_xy_phi)")
        

   
        if not self.storeHists: 
            return
            
        self.results.append(self.df.HistoBoost("recoil_uncorr_magn", [self.axis_recoil_magn], ["recoil_uncorr_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_uncorr_para_qT", [self.axis_recoil_para_qT], ["recoil_uncorr_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_uncorr_para", [self.axis_recoil_para], ["recoil_uncorr_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_uncorr_perp", [self.axis_recoil_perp], ["recoil_uncorr_perp", "nominal_weight"]))
            
        self.results.append(self.df.HistoBoost("recoil_corr_lep_magn", [self.axis_recoil_magn], ["recoil_corr_lep_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_lep_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_lep_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_lep_para", [self.axis_recoil_para], ["recoil_corr_lep_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_lep_perp", [self.axis_recoil_perp], ["recoil_corr_lep_perp", "nominal_weight"]))
            
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn", [self.axis_recoil_magn], ["recoil_corr_xy_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para", [self.axis_recoil_para], ["recoil_corr_xy_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp", [self.axis_recoil_perp], ["recoil_corr_xy_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn_qTrw", [self.axis_recoil_magn], ["recoil_corr_xy_magn", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_qTrw", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qTrw", [self.axis_recoil_para], ["recoil_corr_xy_para", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_qTrw", [self.axis_recoil_perp], ["recoil_corr_xy_perp", "nominal_weight_qTrw"]))
        
        self.results.append(self.df.HistoBoost("MET_corr_xy_pt_qTrw", [self.axis_MET_pt], ["MET_corr_xy_pt", "nominal_weight_qTrw"]))
        
        
        self.results.append(self.df.HistoBoost("mT_uncorr", [self.axis_mt], ["mT_uncorr", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("mT_corr_lep", [self.axis_mt], ["mT_corr_lep", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("mT_corr_xy", [self.axis_mt], ["mT_corr_xy", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("mT_corr_xy_qTrw", [self.axis_mt], ["mT_corr_xy", "nominal_weight_qTrw"]))
        
        

        # for the correction, binned in qT
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn_qTbinned", [self.axis_qT, self.axis_recoil_magn_fine], ["qT", "recoil_corr_xy_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qTbinned", [self.axis_qT, self.axis_recoil_para_fine], ["qT", "recoil_corr_xy_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_qTbinned", [self.axis_qT, self.axis_recoil_para_fine], ["qT", "recoil_corr_xy_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_qTbinned", [self.axis_qT, self.axis_recoil_perp_fine], ["qT", "recoil_corr_xy_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn_qTbinned_qTrw", [self.axis_qT, self.axis_recoil_magn_fine], ["qT", "recoil_corr_xy_magn", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qTbinned_qTrw", [self.axis_qT, self.axis_recoil_para_fine], ["qT", "recoil_corr_xy_para", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_qTbinned_qTrw", [self.axis_qT, self.axis_recoil_para_fine], ["qT", "recoil_corr_xy_para_qT", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_qTbinned_qTrw", [self.axis_qT, self.axis_recoil_perp_fine], ["qT", "recoil_corr_xy_perp", "nominal_weight_qTrw"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_npv", [self.axis_npv, self.axis_recoil_para_fine], ["PV_npvs", "recoil_corr_xy_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_npv", [self.axis_npv, self.axis_recoil_perp_fine], ["PV_npvs", "recoil_corr_xy_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_y", [self.axis_rapidity, self.axis_recoil_para_fine], ["yZ", "recoil_corr_xy_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_y", [self.axis_rapidity, self.axis_recoil_perp_fine], ["yZ", "recoil_corr_xy_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_sumEt", [self.axis_sumEt, self.axis_recoil_para_fine], ["RawMET_sumEt", "recoil_corr_xy_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_sumEt", [self.axis_sumEt, self.axis_recoil_perp_fine], ["RawMET_sumEt", "recoil_corr_xy_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_perp", [self.axis_recoil_para, self.axis_recoil_perp], ["recoil_corr_xy_para", "recoil_corr_xy_perp", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("qT", [self.axis_qT], ["qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("qT_qTrw", [self.axis_qT], ["qT", "nominal_weight_qTrw"]))


    def auxHists(self):
        
        self.df = self.df.Define("njets", "Jet_pt.size()")
        self.results.append(self.df.HistoBoost("njets", [self.axis_njets], ["njets", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("npv", [self.axis_npv], ["PV_npvs", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("RawMET_sumEt", [self.axis_sumEt], ["RawMET_sumEt", "nominal_weight"]))
        
        # recoil resolutions vs any
        

        

        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_njets", [self.axis_njets, self.axis_recoil_para_fine], ["njets", "recoil_corr_xy_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_njets", [self.axis_njets, self.axis_recoil_perp_fine], ["njets", "recoil_corr_xy_perp", "nominal_weight"]))
        
        # correlation sumET vs qT
        self.results.append(self.df.HistoBoost("qT_sumEt", [self.axis_qT, self.axis_sumEt], ["qT", "RawMET_sumEt", "nominal_weight"]))
 

    def setup_recoil_gen(self):
        
        if not self.dataset.name in self.datasets_to_apply: 
            return 
        
        # decompose MET and dilepton (for Z) or lepton (for W) along the generator boson direction
        if self.flavor == "mumu" or self.flavor == "ee": # dilepton
            self.df = self.df.Define("recoil_corr_xy_gen", "wrem::recoilComponentsGen(MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi(), phiVgen)")
        else: # single lepton
            self.df = self.df.Define("recoil_corr_xy_gen", f"wrem::recoilComponentsGen(MET_corr_xy_pt, MET_corr_xy_phi, {self.leptons_pt}, {self.leptons_phi}, phiVgen)")


        self.df = self.df.Define("qT_gen", "ptVgen") # pre-fsr defines should be loaded
        
        self.df = self.df.Define("recoil_corr_xy_magn_gen", "recoil_corr_xy_gen[0]")
        self.df = self.df.Define("recoil_corr_xy_para_gen", "recoil_corr_xy_gen[1] + qT_gen")
        self.df = self.df.Define("recoil_corr_xy_para_qT_gen", "recoil_corr_xy_gen[1]")
        self.df = self.df.Define("recoil_corr_xy_perp_gen", "recoil_corr_xy_gen[2]")
        
        if not self.storeHists: 
            return
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn_gen", [self.axis_recoil_magn], ["recoil_corr_xy_magn_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_gen", [self.axis_recoil_para], ["recoil_corr_xy_para_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_gen", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_gen", [self.axis_recoil_perp], ["recoil_corr_xy_perp_gen", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_xy_magn_gen_qTbinned", [self.axis_qT, self.axis_recoil_magn_fine], ["qT_gen", "recoil_corr_xy_magn_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_gen_qTbinned", [self.axis_qT, self.axis_recoil_para_fine], ["qT_gen", "recoil_corr_xy_para_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_para_qT_gen_qTbinned", [self.axis_qT, self.axis_recoil_para_fine], ["qT_gen", "recoil_corr_xy_para_qT_gen", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_xy_perp_gen_qTbinned", [self.axis_qT, self.axis_recoil_perp_fine], ["qT_gen", "recoil_corr_xy_perp_gen", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("qT_gen", [self.axis_qT], ["qT_gen", "nominal_weight"]))
        
     

    def apply_recoil_Z(self): 

        if self.dataset.name in self.datasets_to_apply:
            
            if self.parametric: self.df = self.df.Define("recoil_corr_rec", "wrem::recoilCorrectionParametric(recoil_corr_xy_para, recoil_corr_xy_perp, qT)")
            else: self.df = self.df.Define("recoil_corr_rec", "wrem::recoilCorrectionBinned(recoil_corr_xy_para_qT, recoil_corr_xy_perp, qTbin, qT)") # 
            self.df = self.df.Define("recoil_corr_rec_magn", "recoil_corr_rec[0]")
            self.df = self.df.Define("recoil_corr_rec_para_qT", "recoil_corr_rec[1] - qT")
            self.df = self.df.Define("recoil_corr_rec_para", "recoil_corr_rec[1]")
            self.df = self.df.Define("recoil_corr_rec_perp", "recoil_corr_rec[2]")
            
            self.df = self.df.Define("MET_corr_rec", "wrem::METCorrection(MET_corr_xy_pt, MET_corr_xy_phi, recoil_corr_rec_para_qT, recoil_corr_rec_perp, qT, Z_mom2.Phi()) ")
            self.df = self.df.Define("MET_corr_rec_pt", "MET_corr_rec[0]")
            self.df = self.df.Define("MET_corr_rec_phi", "MET_corr_rec[1]")
            
        else:
            self.df = self.df.Alias("recoil_corr_rec_magn", "recoil_corr_xy_magn")
            self.df = self.df.Alias("recoil_corr_rec_para_qT", "recoil_corr_xy_para_qT")
            self.df = self.df.Define("recoil_corr_rec_para", "recoil_corr_xy_para")
            self.df = self.df.Alias("recoil_corr_rec_perp", "recoil_corr_xy_perp")
            self.df = self.df.Alias("MET_corr_rec_pt", "MET_corr_xy_pt")
            self.df = self.df.Alias("MET_corr_rec_phi", "MET_corr_xy_phi")

       
        self.df = self.df.Define("MET_corr_rec_xy_dPhi", "wrem::deltaPhi(MET_corr_rec_phi,MET_corr_xy_phi)")
        self.df = self.df.Define("METx_corr_rec", "MET_corr_rec_pt*cos(MET_corr_rec_phi)")
        self.df = self.df.Define("METy_corr_rec", "MET_corr_rec_pt*sin(MET_corr_rec_phi)")
        self.df = self.df.Define("mT_corr_rec", f"wrem::mt_wlike_nano({self.trgLep_pt}, {self.trgLep_phi}, {self.nonTrgLep_pt}, {self.nonTrgLep_phi}, MET_corr_rec_pt, MET_corr_rec_phi)")
        
        if not self.storeHists: 
            return
        
        self.results.append(self.df.HistoBoost("recoil_corr_rec_magn", [self.axis_recoil_magn], ["recoil_corr_rec_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_rec_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para", [self.axis_recoil_para], ["recoil_corr_rec_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_perp", [self.axis_recoil_perp], ["recoil_corr_rec_perp", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("MET_corr_rec_pt", [self.axis_MET_pt], ["MET_corr_rec_pt", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_phi", [self.axis_MET_phi], ["MET_corr_rec_phi", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("MET_corr_rec_pt_qT", [self.axis_qT, self.axis_MET_pt], ["qT", "MET_corr_rec_pt", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("recoil_corr_rec_magn_qTbinned", [self.axis_qT, self.axis_recoil_magn], ["qT", "recoil_corr_rec_magn", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_qT_qTbinned", [self.axis_qT, self.axis_recoil_para_qT], ["qT", "recoil_corr_rec_para_qT", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_qTbinned", [self.axis_qT, self.axis_recoil_para], ["qT", "recoil_corr_rec_para", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_perp_qTbinned", [self.axis_qT, self.axis_recoil_perp], ["qT", "recoil_corr_rec_perp", "nominal_weight"]))
        
        
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_perp_2d", [self.axis_recoil_2d_para, self.axis_recoil_2d_perp], ["recoil_corr_rec_para", "recoil_corr_rec_perp", "nominal_weight"]))

        self.results.append(self.df.HistoBoost("METx_corr_rec", [self.axis_MET_xy], ["METx_corr_rec", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("METy_corr_rec", [self.axis_MET_xy], ["METy_corr_rec", "nominal_weight"]))
        
        self.results.append(self.df.HistoBoost("mT_corr_rec", [self.axis_mt], ["mT_corr_rec", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("mT_corr_rec_qTrw", [self.axis_mt], ["mT_corr_rec", "nominal_weight_qTrw"]))
        
        
        self.df = self.df.Define("MET_corr_xy_ll_dPhi", "wrem::deltaPhi(MET_corr_xy_phi,Z_mom2.Phi())")
        self.df = self.df.Define("MET_corr_rec_ll_dPhi", "wrem::deltaPhi(MET_corr_rec_phi,Z_mom2.Phi())")
        self.df = self.df.Define("ll_phi", "Z_mom2.Phi()")
        self.results.append(self.df.HistoBoost("MET_corr_xy_ll_dPhi", [self.axis_MET_phi], ["MET_corr_xy_ll_dPhi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_ll_dPhi", [self.axis_MET_phi], ["MET_corr_rec_ll_dPhi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_xy_dPhi", [self.axis_MET_dphi], ["MET_corr_rec_xy_dPhi", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("ll_phi", [self.axis_MET_phi], ["ll_phi", "nominal_weight"]))
        

        # for validation, reoil plots with weighted qT spectrum
        self.results.append(self.df.HistoBoost("recoil_corr_rec_magn_qTrw", [self.axis_recoil_magn], ["recoil_corr_rec_magn", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_qT_qTrw", [self.axis_recoil_para_qT], ["recoil_corr_rec_para_qT", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_para_qTrw", [self.axis_recoil_para], ["recoil_corr_rec_para", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_perp_qTrw", [self.axis_recoil_perp], ["recoil_corr_rec_perp", "nominal_weight_qTrw"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_pt_qTrw", [self.axis_MET_pt], ["MET_corr_rec_pt", "nominal_weight_qTrw"]))

    
     
    def apply_recoil_W(self): 
    

        if self.dataset.name in self.datasets_to_apply:
        
            self.df = self.df.Define("qTweight_", "wrem::qTweight(qT_gen)")
            if not self.args.theoryCorrAltOnly and self.args.theoryCorr:
                self.df = self.df.Define("nominal_weight_qTrw", "nominal_weight_uncorr*qTweight_")
            else:
                self.df = self.df.Define("nominal_weight_qTrw", "nominal_weight*qTweight_")
            
            ''' 
            self.df = self.df.Define("recoil_corr_wz", "wrem::recoilCorrectionBinnedWtoZ(Lep_charge, recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen)")
            self.df = self.df.Define("recoil_corr_wz_magn", "recoil_corr_wz[0]")
            self.df = self.df.Define("recoil_corr_wz_para", "recoil_corr_wz[1]")
            self.df = self.df.Define("recoil_corr_wz_para_qT", "recoil_corr_wz[1] + qT_gen")
            self.df = self.df.Define("recoil_corr_wz_perp", "recoil_corr_wz[2]")
            self.df = self.df.Define("recoil_corr_wz_para_qT_perp", "recoil_corr_wz[1] + qT_gen + recoil_corr_wz[2]")

            self.df = self.df.Define("MET_corr_wz", "wrem::METCorrectionGen(recoil_corr_wz_para, recoil_corr_wz_perp, Lep_pt, Lep_phi, phiVgen) ") # was qTgen
            self.df = self.df.Define("MET_corr_wz_pt", "MET_corr_wz[0]")
            self.df = self.df.Define("MET_corr_wz_phi", "MET_corr_wz[1]")
            
            self.results.append(self.df.HistoBoost("recoil_corr_wz_magn", [self.axis_recoil_magn], ["recoil_corr_wz_magn", "nominal_weight"]))
            self.results.append(self.df.HistoBoost("recoil_corr_wz_para", [self.axis_recoil_para], ["recoil_corr_wz_para", "nominal_weight"]))
            self.results.append(self.df.HistoBoost("recoil_corr_wz_para_qT", [self.axis_recoil_perp], ["recoil_corr_wz_para_qT", "nominal_weight"]))
            self.results.append(self.df.HistoBoost("recoil_corr_wz_perp", [self.axis_recoil_perp], ["recoil_corr_wz_perp", "nominal_weight"]))
            '''     

            
            if self.parametric: self.df = self.df.Define("recoil_corr_rec", "wrem::recoilCorrectionParametric(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qT_gen)")
            else: self.df = self.df.Define("recoil_corr_rec", "wrem::recoilCorrectionBinned(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen)")
            self.df = self.df.Define("recoil_corr_rec_magn_gen", "recoil_corr_rec[0]")
            self.df = self.df.Define("recoil_corr_rec_para_gen", "recoil_corr_rec[1]")
            self.df = self.df.Define("recoil_corr_rec_para_qT_gen", "recoil_corr_rec[1] - qT_gen")
            self.df = self.df.Define("recoil_corr_rec_perp_gen", "recoil_corr_rec[2]")
            
            self.df = self.df.Define("MET_corr_rec", f"wrem::METCorrectionGen(recoil_corr_rec_para_qT_gen, recoil_corr_rec_perp_gen, {self.leptons_pt}, {self.leptons_phi}, phiVgen) ") 
            self.df = self.df.Define("MET_corr_rec_pt", "MET_corr_rec[0]")
            self.df = self.df.Define("MET_corr_rec_phi", "MET_corr_rec[1]")
            
            # compute recoil
            self.df = self.df.Define("recoil_corr_rec_magn", f"wrem::recoilComponents(MET_corr_rec_pt, MET_corr_rec_phi, {self.leptons_pt}, {self.leptons_phi})")
            
            
        else:
        
            self.df = self.df.Alias("nominal_weight_qTrw", "nominal_weight")
            # recoil not defined for backgrounds, only MET
            #self.df = self.df.Alias("recoil_corr_wz_magn", "recoil_corr_xy_magn")
            #self.df = self.df.Alias("recoil_corr_wz_para", "recoil_corr_xy_para")
            #self.df = self.df.Define("recoil_corr_wz_para_qT", "recoil_corr_xy_para + qT")
            #self.df = self.df.Alias("recoil_corr_wz_perp", "recoil_corr_xy_perp")
            #self.df = self.df.Alias("recoil_corr_wz_para_qT_perp", "recoil_corr_xy_para_qT_perp")
            
            self.df = self.df.Alias("MET_corr_wz_pt", "MET_corr_xy_pt")
            self.df = self.df.Alias("MET_corr_wz_phi", "MET_corr_xy_phi")
            
            self.df = self.df.Define("MET_corr_rec_pt", "MET_corr_xy_pt")
            self.df = self.df.Define("MET_corr_rec_phi", "MET_corr_xy_phi")
            
            self.df = self.df.Define("recoil_corr_rec_magn", f"wrem::recoilComponents(MET_corr_xy_pt, MET_corr_xy_phi, {self.leptons_pt}, {self.leptons_phi})")
        
        self.df = self.df.Define("mT_corr_lep", f"wrem::mt_2({self.leptons_pt}, {self.leptons_phi}, MET_corr_lep_pt, MET_corr_lep_phi)")
        self.df = self.df.Define("mT_corr_xy", f"wrem::mt_2({self.leptons_pt}, {self.leptons_phi}, MET_corr_xy_pt, MET_corr_xy_phi)")
        self.df = self.df.Define("mT_corr_rec", f"wrem::mt_2({self.leptons_pt}, {self.leptons_phi}, MET_corr_rec_pt, MET_corr_rec_phi)")
        self.df = self.df.Define("passMT_lep", "mT_corr_lep > 40")
        self.df = self.df.Define("passMT_xy", "mT_corr_xy > 40")
        self.df = self.df.Define("passMT_rec", "mT_corr_rec > 40")
        
        
        if not self.storeHists: 
            return
        

        cols_mT = ["mT_corr_rec", self.leptons_charge, "passMT_rec", "passIso"]
        cols_MET_pt = ["MET_corr_rec_pt", self.leptons_charge, "passMT_rec", "passIso"]
        cols_MET_phi = ["MET_corr_rec_phi", self.leptons_charge, "passMT_rec", "passIso"]
        
        cols_mT_xy = ["mT_corr_xy", self.leptons_charge, "passMT_xy", "passIso"]
        cols_MET_pt_xy = ["MET_corr_xy_pt", self.leptons_charge, "passMT_xy", "passIso"]
        cols_MET_phi_xy = ["MET_corr_xy_phi", self.leptons_charge, "passMT_xy", "passIso"]
        
        cols_mT_lep = ["mT_corr_lep", self.leptons_charge, "passMT_lep", "passIso"]
        cols_MET_pt_lep = ["MET_corr_lep_pt", self.leptons_charge, "passMT_lep", "passIso"]
        cols_MET_phi_lep = ["MET_corr_lep_phi", self.leptons_charge, "passMT_lep", "passIso"]
        
        axes_mT = [self.axis_mt, self.axis_charge, self.axis_passMT, self.axis_passIso]
        axes_MET_pt = [self.axis_MET_pt, self.axis_charge, self.axis_passMT, self.axis_passIso]
        axes_MET_phi = [self.axis_MET_phi, self.axis_charge, self.axis_passMT, self.axis_passIso]
        
        self.results.append(self.df.HistoBoost("mT_corr_lep", axes_mT, cols_mT_lep + ["nominal_weight"]))     
        self.results.append(self.df.HistoBoost("mT_corr_xy", axes_mT, cols_mT_xy + ["nominal_weight"]))      
        self.results.append(self.df.HistoBoost("mT_corr_rec", axes_mT, cols_mT + ["nominal_weight"]))  
        self.results.append(self.df.HistoBoost("mT_corr_rec_qTrw", axes_mT, cols_mT + ["nominal_weight_qTrw"]))  

        self.results.append(self.df.HistoBoost("MET_corr_lep_pt", axes_MET_pt, cols_MET_pt_lep + ["nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_xy_pt", axes_MET_pt, cols_MET_pt_xy + ["nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_pt", axes_MET_pt, cols_MET_pt + ["nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_pt_qTrw", axes_MET_pt, cols_MET_pt + ["nominal_weight_qTrw"]))
            
        self.results.append(self.df.HistoBoost("MET_corr_lep_phi", axes_MET_phi, cols_MET_phi_lep + ["nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_xy_phi", axes_MET_phi, cols_MET_phi_xy + ["nominal_weight"]))
        self.results.append(self.df.HistoBoost("MET_corr_rec_phi", axes_MET_phi, cols_MET_phi + ["nominal_weight"]))
            
        self.results.append(self.df.HistoBoost("recoil_corr_rec_magn", [self.axis_recoil_magn, self.axis_charge, self.axis_passMT, self.axis_passIso], ["recoil_corr_rec_magn", self.leptons_charge, "passMT_rec", "passIso", "nominal_weight"]))
        self.results.append(self.df.HistoBoost("recoil_corr_rec_magn_qTrw", [self.axis_recoil_magn, self.axis_charge, self.axis_passMT, self.axis_passIso], ["recoil_corr_rec_magn", self.leptons_charge, "passMT_rec", "passIso", "nominal_weight_qTrw"]))
        
        
    

    def recoil_Z_statUnc_lowPU(self, df, results, axis_recoil_gen, axis_recoil_reco, axis_mt, axis_mll):
        
        # lowPU: propagate uncertainties to recoil, MET and mT
        
        axis_recoil_stat_unc = hist.axis.Regular(len(self.recoil_qTbins), 0, len(self.recoil_qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")
        
        # make copies of axes for perturbation
        axis_recoil_reco_pert = hist.axis.Variable(axis_recoil_reco.edges, name = "recoil_reco_pert", underflow=False, overflow=True)
        axis_mt_pert = hist.axis.Variable(axis_mt.edges, name = "mt_pert",underflow=False, overflow=True)
        #axis_recoil_reco_pert = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 150], name = "recoil_reco_pert", underflow=False, overflow=True)
        axis_MET_pert = hist.axis.Variable(self.axis_MET_pt.edges, name = "recoil_MET_pt_pert", underflow=False, overflow=True)
        axis_recoil_magn_pert = hist.axis.Variable(self.axis_recoil_magn.edges, name = "recoil_magn_pert", underflow=False, overflow=True)
        
        
        axis_recoil_para_pert = hist.axis.Variable(self.axis_recoil_para.edges, name = "recoil_para_pert")
        axis_recoil_para_qT_pert = hist.axis.Variable(self.axis_recoil_para_qT.edges, name = "recoil_para_qT_pert", underflow=False, overflow=False)
        axis_recoil_perp_pert = hist.axis.Variable(self.axis_recoil_perp.edges, name = "recoil_perp_pert", underflow=False, overflow=False)
    
        # recoil stat uncertainty  (para_data/perp_data/para_mc/perp_mc, number of variations)
        self.df = self.df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(self.recoil_qTbins)))
        recoil_vars = [(1,2), (1,3), (1,4),   (2,2), (2,3),   (3,2), (3,3), (3,4),    (4,2), (4,3)]
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),(3,6),(3,7),  (4,2),(4,3),(4,4),(4,5)] # 
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8),  (4,2),(4,3),(4,4),(4,5)] # 
        
        
        # 2 gauss for data = 4, 3 gauss for MC
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
        #recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
        for k in recoil_vars:

            # perturbation for current qTbin
            self.df = self.df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT, %d, %d)" % (k[0], k[1]))
           
            self.df = self.df.Define("recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_magn_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            self.results.append(self.df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "nominal_weight"]))
            
            
            #self.df = self.df.Define("recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            #self.results.append(self.df.HistoBoost("gen_reco_para_qT_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
            self.results.append(self.df.HistoBoost("recoil_magn_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_magn, axis_recoil_magn_pert, axis_recoil_stat_unc], ["recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
            
            self.df = self.df.Define("recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_qT_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, qT)" % (k[0], k[1]))
            self.results.append(self.df.HistoBoost("recoil_para_qT_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para_qT, axis_recoil_para_qT_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para_qT", "recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
            self.df = self.df.Define("recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_perp_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            self.results.append(self.df.HistoBoost("recoil_perp_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_perp, axis_recoil_perp_pert, axis_recoil_stat_unc], ["recoil_corr_rec_perp", "recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
            
            self.df = self.df.Define("recoil_corr_MET_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
            self.results.append(self.df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc], ["MET_corr_rec_pt", "recoil_corr_MET_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
            self.df = self.df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_lep_pt, MET_corr_lep_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
            
            
            self.df = self.df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_StatUnc(recoil_corr_MET_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (k[0], k[1], k[0], k[1]))
            self.results.append(self.df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
        
        return df
        

    def recoil_Z_unc(self, df, results, dataset, datasets_to_apply, hNames=[], cols=[], axes=[]):
        
        if not dataset.name in datasets_to_apply:
            return df
        
        cols_base = ["recoil_corr_rec_para_qT", "recoil_corr_rec_para", "recoil_corr_rec_perp", "recoil_corr_rec_magn", "MET_corr_rec_pt", "mT_corr_rec"]
        if self.storeHists:
            cols.extend(cols_base)
            hNames.extend(cols_base)
            axes.extend([self.axis_recoil_para_qT, self.axis_recoil_para, self.axis_recoil_perp, self.axis_recoil_magn, self.axis_MET_pt, self.axis_mt])
        
        if self.parametric and self.smearWeights:
            
            recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
            for item in recoil_param_nStat:
                tag = item.first
                nVars = recoil_param_nStat[tag]
              
                val = "recoil_corr_xy_para" if "para" in tag else "recoil_corr_xy_perp"
                tag_nom = "%s_%s" % ("target" if "target" in tag else "source", "para" if "para" in tag else "perp")
                recoilWeights = f"recoilWeights_{tag}"
                recoilTensorWeights = f"recoilTensorWeights_{tag}"
                df = df.Define(recoilWeights, "wrem::recoilCorrectionParametricUncWeights(%s, qT,  \"%s\", \"%s\")" % (val, tag_nom, tag))
                df = df.Define(recoilTensorWeights, "Eigen::TensorFixedSize<double, Eigen::Sizes<%d>> res; auto w = nominal_weight*%s; std::copy(std::begin(w), std::end(w), res.data()); return res;" % (nVars, recoilWeights))

                for hName, col, ax in zip(hNames, cols, axes):
                    results.append(df.HistoBoost(f"{hName}_recoilSyst_{tag}", ax if isinstance(col, list) else [ax], (col if isinstance(col, list) else [col]) + [recoilTensorWeights]))
                    
                    
                # qT reweighted uncertainties, only for base_cols
                recoilTensorWeights_qTrw = f"recoilTensorWeights_qTrw_{tag}"
                df = df.Define(recoilTensorWeights_qTrw, "Eigen::TensorFixedSize<double, Eigen::Sizes<%d>> res; auto w = nominal_weight_qTrw*%s; std::copy(std::begin(w), std::end(w), res.data()); return res;" % (nVars, recoilWeights))
                for hName, col, ax in zip(hNames, cols, axes):
                    if not hName in cols_base and not self.storeHists: continue
                    results.append(df.HistoBoost(f"{hName}_qTrw_recoilSyst_{tag}", ax if isinstance(col, list) else [ax], (col if isinstance(col, list) else [col]) + [recoilTensorWeights_qTrw]))

               
        elif self.parametric and not self.smearWeights:
        
            recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
            for item in recoil_param_nStat:
                tag = item.first
                nVars = recoil_param_nStat[tag]
                


                axis_recoil_unc = hist.axis.Regular(nVars, 0, nVars, underflow=False, overflow=False, name = "axis_recoil_unc_%s" % tag)
                df = df.Define("recoil_corr_rec_syst_%s" % tag, "wrem::recoilCorrectionParametricUnc(recoil_corr_xy_para, recoil_corr_xy_perp, qT, \"%s\")" % tag)
                df = df.Define("recoil_corr_rec_systIdx_%s" % tag, "wrem::indices_(%d, 0)" % nVars)
                    
                # MET
                #df = df.Define("recoil_corr_MET_pt_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_pt_unc(recoil_corr_rec_syst_%s, qT, Z_mom2.Phi())" % tag)
                #df = df.Define("recoil_corr_MET_phi_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_phi_unc(recoil_corr_rec_syst_%s, qT, Z_mom2.Phi())" % tag)
                #results.append(df.HistoBoost("MET_recoilSyst_%s" % tag, [axis_MET_pt, axis_recoil_unc], ["recoil_corr_MET_pt_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight"]))
                
                
                # recoil
                df = df.Define("recoil_corr_para_qT_syst_%s" % tag, "wrem::recoilCorrectionParametric_para_qT_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_para_syst_%s" % tag, "wrem::recoilCorrectionParametric_para_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_perp_syst_%s" % tag, "wrem::recoilCorrectionParametric_perp_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_magn_syst_%s" % tag, "wrem::recoilCorrectionParametric_magn_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                results.append(df.HistoBoost("recoil_corr_rec_para_qT_recoilSyst_%s" % tag, [self.axis_recoil_para_qT, axis_recoil_unc], ["recoil_corr_para_qT_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight_qTrw"]))
                results.append(df.HistoBoost("recoil_corr_rec_para_recoilSyst_%s" % tag, [self.axis_recoil_para, axis_recoil_unc], ["recoil_corr_para_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight_qTrw"]))
                results.append(df.HistoBoost("recoil_corr_rec_perp_recoilSyst_%s" % tag, [self.axis_recoil_perp, axis_recoil_unc], ["recoil_corr_perp_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight_qTrw"]))
                results.append(df.HistoBoost("recoil_corr_rec_magn_recoilSyst_%s" % tag, [self.axis_recoil_magn, axis_recoil_unc], ["recoil_corr_magn_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight_qTrw"]))
                
                # mT
                #df = df.Define("recoil_corr_mT_recoilSyst_%s" % tag, "wrem::recoilCorrectionParametric_mT_unc(recoil_corr_MET_pt_syst_%s, recoil_corr_MET_phi_syst_%s, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (tag, tag))
                #results.append(df.HistoBoost("mT_corr_rec_recoilSyst_%s" % tag, [axis_mt, axis_recoil_unc], ["recoil_corr_mT_recoilSyst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight"]))
            
                # Z peak
                #results.append(df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "nominal_weight"]))
              

        else:

            axis_recoil_stat_unc = hist.axis.Regular(len(recoil_qTbins), 0, len(recoil_qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")
            
            # make copies of axes for perturbation
            axis_recoil_reco_pert = hist.axis.Variable(axis_recoil_reco.edges, name = "recoil_reco_pert", underflow=False, overflow=True)
            axis_mt_pert = hist.axis.Variable(axis_mt.edges, name = "mt_pert",underflow=False, overflow=True)
            axis_MET_pert = hist.axis.Variable(axis_MET_pt.edges, name = "recoil_MET_pt_pert", underflow=False, overflow=True)
            
            axis_recoil_para_qT_pert = hist.axis.Variable(axis_recoil_para_qT.edges, name = "recoil_para_qT_pert", underflow=False, overflow=False)
            axis_recoil_para_pert = hist.axis.Variable(axis_recoil_para.edges, name = "recoil_para_pert", underflow=False, overflow=False)
            axis_recoil_magn_pert = hist.axis.Variable(axis_recoil_magn.edges, name = "recoil_magn_pert", underflow=False)
            axis_recoil_perp_pert = hist.axis.Variable(axis_recoil_perp.edges, name = "recoil_perp_pert", underflow=False, overflow=False)
            
        
            # recoil stat uncertainty  (para_data/perp_data/para_mc/perp_mc, number of variations)
            df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(recoil_qTbins)))
            recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)] # 2 gauss for data = 4, 3 gauss for MC
            for k in recoil_vars:

                # perturbation for current qTbin
                df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT, %d, %d)" % (k[0], k[1]))
               
                # MET
                df = df.Define("recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
                df = df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
                results.append(df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc], ["MET_corr_rec_pt", "recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))

                # recoil
                df = df.Define("recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_qT_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, qT)" % (k[0], k[1]))
                df = df.Define("recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                df = df.Define("recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_perp_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                df = df.Define("recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_magn_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                results.append(df.HistoBoost("recoil_para_qT_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_para_qT, axis_recoil_para_qT_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para_qT", "recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
                results.append(df.HistoBoost("recoil_para_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
                results.append(df.HistoBoost("recoil_perp_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_perp, axis_recoil_perp_pert, axis_recoil_stat_unc], ["recoil_corr_rec_perp", "recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
                results.append(df.HistoBoost("recoil_magn_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_magn, axis_recoil_magn_pert, axis_recoil_stat_unc], ["recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))

                # mT
                df = df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_StatUnc(recoil_corr_MET_pt_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (k[0], k[1], k[0], k[1]))
                results.append(df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
                
                # Z peak
                #results.append(df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "nominal_weight"]))
                #results.append(df.HistoBoost("gen_reco_para_qT_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))

  
        return df
        
        
    def recoil_W_unc(self, df, results, dataset, datasets_to_apply, hNames=[], cols=[], axes=[]):
    
        if not dataset.name in datasets_to_apply:
            return df
    
        cols_base = ["MET_corr_rec_pt", "mT_corr_rec"]
        if self.storeHists:
            cols.extend(cols_base)
            hNames.extend(cols_base)
            axes.extend([self.axis_MET_pt, self.axis_mt])
    
        if self.parametric and self.smearWeights:
            
            recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
            for item in recoil_param_nStat:
                tag = item.first
                nVars = recoil_param_nStat[tag]
              
                val = "recoil_corr_xy_para_gen" if "para" in tag else "recoil_corr_xy_perp_gen"
                tag_nom = "%s_%s" % ("target" if "target" in tag else "source", "para" if "para" in tag else "perp")
                recoilWeights = f"recoilWeights_{tag}"
                recoilTensorWeights = f"recoilTensorWeights_{tag}"
                df = df.Define(recoilWeights, "wrem::recoilCorrectionParametricUncWeights(%s, qT_gen,  \"%s\", \"%s\")" % (val, tag_nom, tag))
                df = df.Define(recoilTensorWeights, "Eigen::TensorFixedSize<double, Eigen::Sizes<%d>> res; auto w = nominal_weight*%s; std::copy(std::begin(w), std::end(w), res.data()); return res;" % (nVars, recoilWeights))

                for hName, col, ax in zip(hNames, cols, axes):
                    results.append(df.HistoBoost(f"{hName}_recoilSyst_{tag}", ax if isinstance(ax, list) else [ax], (col if isinstance(col, list) else [col]) + [recoilTensorWeights]))
                    
                    
                # qT reweighted uncertainties
                recoilTensorWeights_qTrw = f"recoilTensorWeights_qTrw_{tag}"
                df = df.Define(recoilTensorWeights_qTrw, "Eigen::TensorFixedSize<double, Eigen::Sizes<%d>> res; auto w = nominal_weight_qTrw*%s; std::copy(std::begin(w), std::end(w), res.data()); return res;" % (nVars, recoilWeights))
                for hName, col, ax in zip(hNames, cols, axes):
                    if not hName in cols_base and not self.storeHists: continue
                    results.append(df.HistoBoost(f"{hName}_qTrw_recoilSyst_{tag}", ax if isinstance(ax, list) else [ax], (col if isinstance(col, list) else [col]) + [recoilTensorWeights_qTrw]))

        
        elif self.parametric and not self.smearWeights:
        
            recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
            for item in recoil_param_nStat:
                tag = item.first
                nVars = recoil_param_nStat[tag]

                axis_recoil_unc = hist.axis.Regular(nVars, 0, nVars, underflow=False, overflow=False, name = "axis_recoil_unc_%s" % tag)
                df = df.Define("recoil_corr_rec_syst_%s" % tag, "wrem::recoilCorrectionParametricUnc(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qT_gen, \"%s\")" % tag)
                df = df.Define("recoil_corr_rec_systIdx_%s" % tag, "wrem::indices_(%d, 0)" % nVars)
                    
                # MET 
                df = df.Define("recoil_corr_MET_pt_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_pt_gen_unc(recoil_corr_rec_syst_%s, Lep_pt, Lep_phi, phiVgen)" % tag)
                df = df.Define("recoil_corr_MET_phi_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_phi_gen_unc(recoil_corr_rec_syst_%s, Lep_pt, Lep_phi, phiVgen)" % tag)
                results.append(df.HistoBoost("MET_recoilSyst_%s" % tag, [axis_MET_pt, axis_recoil_unc, axis_charge, axis_passIso], ["recoil_corr_MET_pt_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "Lep_charge", "passIso", "nominal_weight"]))
                
                # mT
                df = df.Define("recoil_corr_mT_recoilSyst_%s" % tag, "wrem::recoilCorrectionParametric_mT_2_unc(recoil_corr_MET_pt_syst_%s, recoil_corr_MET_phi_syst_%s, Lep_pt, Lep_phi)" % (tag, tag))
                results.append(df.HistoBoost("mT_corr_rec_recoilSyst_%s" % tag, [axis_mt, axis_charge, axis_passIso, axis_recoil_unc], ["recoil_corr_mT_recoilSyst_%s" % tag, "Lep_charge", "passIso", "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight"]))
           
           
                val = "recoil_corr_xy_para_qT_gen" if "para" in tag else "recoil_corr_xy_perp_gen"
                tag_nom = "%s_%s" % ("target" if "target" in tag else "source", "para" if "para" in tag else "perp")
                tmp = "recoil_corr_rec_syst_weight_%s" % tag
                df = df.Define(tmp, "wrem::recoilCorrectionParametricUncWeights(%s, qT_gen,  \"%s\", \"%s\")" % (val, tag_nom, tag))
                df = df.Define(tmp+"_tensor", "Eigen::TensorFixedSize<double, Eigen::Sizes<%d>> res; auto w = nominal_weight*%s; std::copy(std::begin(w), std::end(w), res.data()); return res;" % (nVars, tmp))

                results.append(df.HistoBoost("mT_corr_rec_recoilSystWeight_%s" % tag, [axis_mt, axis_charge, axis_passMT, axis_passIso], ["mT_corr_rec", "Lep_charge", "passMT", "passIso", tmp+"_tensor"]))
                ###results.append(df.HistoBoost("mT_corr_rec_recoilSyst_weight_%s" % tag, [axis_mt, axis_charge, axis_passIso], ["mT_corr_rec", "Lep_charge", "passIso", tmp+"_tensor"]))
                results.append(df.HistoBoost("recoil_magn_recoilSystWeight_%s" % tag, [axis_recoil_magn, axis_charge, axis_passMT, axis_passIso], ["recoil_corr_rec_magn", "Lep_charge", "passMT", "passIso", tmp+"_tensor"]))
                results.append(df.HistoBoost("MET_recoilSystWeight_%s" % tag, [axis_MET_pt, axis_charge, axis_passMT, axis_passIso], ["MET_corr_rec_pt", "Lep_charge", "passMT", "passIso", tmp+"_tensor"]))

           
        else:

            axis_recoil_stat_unc = hist.axis.Regular(len(recoil_qTbins), 0, len(recoil_qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")
            axis_mt_pert = hist.axis.Variable(axis_mt.edges, name = "mt_pert", underflow=False, overflow=True)
            axis_MET_pert = hist.axis.Variable(axis_MET_pt.edges, name = "recoil_MET_pt_pert", underflow=False, overflow=True)
            
            # recoil stat uncertainty  (para_data/perp_data/para_mc/perp_mc, number of variations)
            df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(recoil_qTbins)))

            # 2 gauss for data = 4, 3 gauss for MC
            recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
            for k in recoil_vars:

                # perturbation for current qTbin
                df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen, %d, %d)" % (k[0], k[1]))
               
                # MET
                df = df.Define("recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_gen_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi, phiVgen)" % (k[0], k[1]))
                df = df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_gen_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi, phiVgen)" % (k[0], k[1]))
                results.append(df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc, axis_eta, axis_charge], ["MET_corr_rec_pt", "recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "Lep_abs_eta", "Lep_charge", "nominal_weight"]))

                # mT
                df = df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_2_StatUnc(recoil_corr_MET_pt_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi)" % (k[0], k[1], k[0], k[1]))
                results.append(df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc, axis_eta, axis_charge], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "Lep_abs_eta", "Lep_charge", "nominal_weight"]))

        return df