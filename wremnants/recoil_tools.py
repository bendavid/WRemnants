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

def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)
        

class Recoil:

    def __init__(self, type_, flavor_, met_):

        self.met = met_
        self.flavor = flavor_   
        self.parametric = True
   
        if type_ == "highPU":
        
            #self.recoil_qTbins = list(drange(0, 30, 0.5)) + list(range(30, 60, 2)) + list(range(60, 100, 5)) + list(range(100, 210, 10)) + [10000]
            self.recoil_qTbins = list(range(0, 50, 1)) + list(range(50, 100, 5)) + list(range(100, 200, 10)) + list(range(200, 400, 20)) + [10000]
            rInput_binned = "wremnants/data/recoil_fits_Z.root"
            rInput_parametric = ""
    
        elif type_ == "lowPU":
        
            self.recoil_qTbins = list(range(0, 30, 1)) + list(range(30, 50, 2)) + list(range(50, 100, 5)) + list(range(100, 150, 10)) + [150, 175, 200, 10000]
            self.recoil_qTbins = list(range(0, 30, 1)) + list(range(30, 50, 2)) + list(range(50, 70, 5)) + list(range(70, 150, 10)) + [150, 200, 300, 10000]
            self.recoil_qTbins = list(range(0, 50, 1)) + list(range(50, 70, 2)) + list(range(70, 100, 5)) + list(range(100, 150, 10)) + [150, 200, 300, 10000]
            #self.recoil_qTbins = [0, 0.5, 1, 1.5] + list(range(2, 50, 1)) + list(range(50, 70, 2)) + list(range(70, 100, 5)) + list(range(100, 150, 10)) + [150, 200, 300, 10000]
            #self.recoil_qTbins = list(drange(0, 30, 0.5)) + list(range(30, 50, 1)) + list(range(50, 70, 2)) + list(range(70, 100, 5)) + list(range(100, 150, 10)) + [150, 200, 300, 10000]
            #self.recoil_qTbins = list(range(0, 300, 1))
            self.recoil_qTbins__ = list(drange(0, 300, 0.5))
            self.recoil_qTbins__ = list(range(0, 300, 1))
            #self.recoil_qTbins = list(range(0, 20, 1)) + list(range(20, 40, 2)) + list(range(40, 55, 3)) + list(range(55, 80, 5)) + list(range(80, 100, 10)) + [100, 125, 150, 10000]
            #print(self.recoil_qTbins)
            
            #rInput_binned = "wremnants/data/lowPU/recoil_fits_Z.root"
            #rInput_binned = "/eos/user/j/jaeyserm/analysis/lowPU/recoil_mumu_DeepMET_lowPU_RoccoR_lowPU_v0_param.root"
            rInput_binned = "/eos/user/j/jaeyserm/analysis/lowPU/recoil_mumu_PFMET_lowPU_RoccoR_lowPU_v0_param.root"
            #if self.flavor == "mu": self.flavor = "mumu"
            #if self.flavor == "ee": self.flavor = "mumu"
            rInput_binned = "wremnants/data/lowPU/recoil_param_%s_%s.root" % ("mumu" if self.flavor == "mu" else self.flavor, self.met)
            rInput_parametric = ""
            
            self.met_xycorr_setup("wremnants/data/lowPU/MET_xy_corr_coeff_%s_%s.json" % (self.flavor, self.met))
            
            if self.flavor == "mu":
                rInput_binned_gen_W = "wremnants/data/lowPU/recoil/param_mu_%s_gen.root" % self.met
                rInput_binned_gen_Z = "wremnants/data/lowPU/recoil/param_mumu_%s_gen.root" % self.met
                ROOT.wrem.recoil_init_gen(rInput_binned_gen_W, rInput_binned_gen_Z)
                
                
            # set the parametric forms 
            #self.addParametric("target_para", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/singlemuon_para/param_v3/results_refit.json")
            self.addBinned("target_para", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/singlemuon_para/fits_v0/results.json")
            self.addParametric("target_perp", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/singlemuon_perp/param_v1/results_refit.json")
            #self.addParametric("source_para", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/dymumu_para/param_v0/results_refit.json")
            self.addBinned("source_para", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/dymumu_para/fits_v0/results.json")
            self.addParametric("source_perp", "/eos/user/j/jaeyserm/www/wmass/lowPU/recoilCorrectionParametric/recoil_mumu_RawPFMET/dymumu_perp/param_v0/results_refit.json")
            
        else: sys.exit("Recoil highPU or lowPU")
        


        # set the recoil correction bins in the analyzer
        qTbins_vec = ROOT.std.vector["float"]()
        for v in self.recoil_qTbins: qTbins_vec.push_back(v)
        setattr(ROOT.wrem, "qTbins", qTbins_vec)
        
        # load recoil hists
        ROOT.wrem.recoil_init(rInput_binned)
        
        
        # define axes
        if self.flavor == "mu": self.axis_MET_pt = hist.axis.Regular(200, 0, 200, name = "recoil_MET_pt", underflow=False)
        else: self.axis_MET_pt = hist.axis.Regular(400, 0, 200, name = "recoil_MET_pt", underflow=False)
        self.axis_MET_phi = hist.axis.Regular(50, -4, 4, name = "recoil_MET_phi")
        self.axis_MET_dphi = hist.axis.Regular(41, -4.1, 4.1, name = "recoil_MET_phi")
        
        
        # Check overflow and binned parametric recoil unc! Seems not to work when there is an under/overflow??
        # Likely there is an interference with recoil_corr_stat_idx, which is defined at 0?
        self.axis_recoil_magn = hist.axis.Regular(150, 0, 300, name = "recoil_magn", underflow=False)
        self.axis_recoil_para = hist.axis.Regular(400, -300, 100, name = "recoil_para", underflow=False, overflow=False)
        #self.axis_recoil_para_qT = hist.axis.Variable([-10000,-100,-80,-70,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,100,10000], name = "recoil_para_qT", underflow=False, overflow=False)
        #self.axis_recoil_perp = hist.axis.Variable([-10000,-100,-80,-70,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,100,10000], name = "recoil_perp", underflow=False, overflow=False)
        self.axis_recoil_para_qT = hist.axis.Variable([-10000] + list(range(-100, 100, 2)) + [100,10000], name = "recoil_para_qT", underflow=False, overflow=False)
        self.axis_recoil_perp = hist.axis.Variable([-10000] + list(range(-100, 100, 2)) + [100,10000], name = "recoil_perp", underflow=False, overflow=False)
        
        self.axis_recoil_para_perp_abs = hist.axis.Variable(list(range(0, 10, 1)) + list(range(10, 30, 2)) + [30, 34, 38, 42, 46, 50, 60, 70, 80, 100, 10000], name = "axis_recoil_para_perp_abs", underflow=False, overflow=False)
        self.axis_recoil_para_perp_squared = hist.axis.Regular(1000, 0, 1000, name = "axis_recoil_para_perp_squared", underflow=False, overflow=False)

        
        self.axis_npv = hist.axis.Regular(50, 0.5, 50.5, name = "recoil_npv")
        self.axis_RawMET_sumEt = hist.axis.Regular(1000, 0, 1000, name = "RawMET_sumEt")
        
        
        
        self.axis_qTbinned = hist.axis.Variable(self.recoil_qTbins__, name = "qTbinned")
        self.axis_qT = hist.axis.Regular(600, 0, 300, name = "qT")
        
        self.axis_recoil_magn_fine = hist.axis.Regular(300, 0, 300, name = "recoil_magn")
        self.axis_recoil_para_fine = hist.axis.Regular(800, -500, 300, name = "recoil_para")
        self.axis_recoil_para_qT_fine = hist.axis.Regular(600, -300, 300, name = "recoil_para_qT")
        self.axis_recoil_perp_fine = hist.axis.Regular(600, -300, 300, name = "recoil_perp")
        self.axis_MET_xy = hist.axis.Regular(200, -100, 100, name = "MET_xy")
        
        self.axis_recoil_2d_para = hist.axis.Regular(40, -10, 10, name = "axis_recoil_2d_para", underflow=False, overflow=False)
        self.axis_recoil_2d_perp = hist.axis.Regular(40, -10, 10, name = "axis_recoil_2d_perp", underflow=False, overflow=False)

    def addParametric(self, tag, fIn):
    
        print("Add recoil parametric set for %s" % tag)
        with open(fIn) as f: jsIn = json.load(f)
        nGauss = jsIn['nGauss']
        recoil_param_nGauss = getattr(ROOT.wrem, "recoil_param_nGauss")
        recoil_param_funcs = getattr(ROOT.wrem, "recoil_param_funcs")
        recoil_param_nGauss.insert((tag, nGauss))
        
        for iGauss in range(1, nGauss+1):
            
            label = "mean%d" % iGauss
            func = ROOT.TF1("%s_%s" % (tag, label), jsIn[label]['func'], 0, 300)
            for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[label]['p%d' % iParam])
            recoil_param_funcs.insert((func.GetName(), func))
            
            label = "sigma%d" % iGauss
            func = ROOT.TF1("%s_%s" % (tag, label), jsIn[label]['func'], 0, 300)
            for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[label]['p%d' % iParam])
            recoil_param_funcs.insert((func.GetName(), func))
            
            if iGauss == nGauss: continue
            label = "norm%d" % iGauss
            func = ROOT.TF1("%s_%s" % (tag, label), jsIn[label]['func'], 0, 300)
            for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[label]['p%d' % iParam])
            recoil_param_funcs.insert((func.GetName(), func))
        
        # do statistical variations
        nStatVars = jsIn['nStatVars']
        recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
        recoil_param_nStat.insert((tag, nStatVars))
        for nStat in range(0, nStatVars):
            statLabel = "stat%d" % nStat
            systLabel = "syst%d" % nStat
            for iGauss in range(1, nGauss+1):
                
                label = "mean%d" % iGauss
                func = ROOT.TF1("%s_%s_%s" % (tag, label, systLabel), jsIn[label]['func'], 0, 300)
                for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[statLabel][label]['p%d' % iParam])
                recoil_param_funcs.insert((func.GetName(), func))
                
                label = "sigma%d" % iGauss
                func = ROOT.TF1("%s_%s_%s" % (tag, label, systLabel), jsIn[label]['func'], 0, 300)
                for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[statLabel][label]['p%d' % iParam])
                recoil_param_funcs.insert((func.GetName(), func))
                
                if iGauss == nGauss: continue
                label = "norm%d" % iGauss
                func = ROOT.TF1("%s_%s_%s" % (tag, label, systLabel), jsIn[label]['func'], 0, 300)
                for iParam in range(0, jsIn[label]['nParams']): func.SetParameter(iParam, jsIn[statLabel][label]['p%d' % iParam])
                recoil_param_funcs.insert((func.GetName(), func))        
        

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
        for qTbin in qTbins: vec.push_back(1.0*qTbin)
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
            print("XY correction file %s not found" % fIn)
            return
        f = open(fIn)
        jsdata = json.load(f)
        loadCoeff(jsdata['x']['data']['nominal'], 'met_xy_corr_x_data_nom')
        loadCoeff(jsdata['y']['data']['nominal'], 'met_xy_corr_y_data_nom')
        loadCoeff(jsdata['x']['mc']['nominal'], 'met_xy_corr_x_mc_nom')
        loadCoeff(jsdata['y']['mc']['nominal'], 'met_xy_corr_y_mc_nom')
        
    def recoil_setup_MET(self, df, results, dataset, leptons_pt, leptons_phi, leptons_uncorr_pt):
        '''
        Define MET variables, apply lepton and XY corrections
        '''
        
        self.leptons_pt = leptons_pt
        self.leptons_phi = leptons_phi
        self.leptons_uncorr_pt = leptons_uncorr_pt
        
        if self.met == "RawPFMET": met_pt, met_phi = "RawMET_pt", "RawMET_phi"
        elif self.met == "DeepMETReso": met_pt, met_phi = "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi"
        else: sys.exit("MET type %s not supported" % self.met)
     
        # uncorrected MET
        df = df.Alias("MET_uncorr_pt", met_pt)
        df = df.Alias("MET_uncorr_phi", met_phi)
        df = df.Define("METx_uncorr", "MET_uncorr_pt*cos(MET_uncorr_phi)")
        df = df.Define("METy_uncorr", "MET_uncorr_pt*sin(MET_uncorr_phi)")

        # lepton corrected MET
        df = df.Define("MET_corr_lep", "wrem::METLeptonCorrection(MET_uncorr_pt, MET_uncorr_phi, %s, %s, %s)" % (leptons_pt, leptons_uncorr_pt, leptons_phi))
        df = df.Define("MET_corr_lep_pt", "MET_corr_lep[0]")
        df = df.Define("MET_corr_lep_phi", "MET_corr_lep[1]")
        df = df.Define("METx_corr_lep", "MET_corr_lep_pt*cos(MET_corr_lep_phi)")
        df = df.Define("METy_corr_lep", "MET_corr_lep_pt*sin(MET_corr_lep_phi)")

        # phi corrected MET (XY corrections)
        df = df.Define("MET_corr_xy", "wrem::METXYCorrection(MET_corr_lep_pt, MET_corr_lep_phi, PV_npvs, %d)" % dataset.is_data)
        df = df.Define("MET_corr_xy_pt", "MET_corr_xy[0]")
        df = df.Define("MET_corr_xy_phi", "MET_corr_xy[1]")
        df = df.Define("METx_corr_xy", "MET_corr_xy_pt*cos(MET_corr_xy_phi)")
        df = df.Define("METy_corr_xy", "MET_corr_xy_pt*sin(MET_corr_xy_phi)")
        
        # MET unclustered correction
        df = df.Define("METx_corr_xy_unclUp", "METx_corr_xy+MET_MetUnclustEnUpDeltaX")
        df = df.Define("METy_corr_xy_unclUp", "METy_corr_xy+MET_MetUnclustEnUpDeltaY")
        df = df.Define("METx_corr_xy_unclDw", "METx_corr_xy-MET_MetUnclustEnUpDeltaX")
        df = df.Define("METy_corr_xy_unclDw", "METy_corr_xy-MET_MetUnclustEnUpDeltaY")
        df = df.Define("MET_corr_xy_pt_unclUp", "sqrt(METx_corr_xy_unclUp*METx_corr_xy_unclUp + METy_corr_xy_unclUp*METy_corr_xy_unclUp)")
        df = df.Define("MET_corr_xy_phi_unclUp", "atan2(METy_corr_xy_unclUp, METx_corr_xy_unclUp)")
        df = df.Define("MET_corr_xy_pt_unclDw", "sqrt(METx_corr_xy_unclDw*METx_corr_xy_unclDw + METy_corr_xy_unclDw*METy_corr_xy_unclDw)")
        df = df.Define("MET_corr_xy_phi_unclDw", "atan2(METy_corr_xy_unclDw, METx_corr_xy_unclDw)")
        
        
        results.append(df.HistoBoost("MET_uncorr_pt", [self.axis_MET_pt], ["MET_uncorr_pt", "nominal_weight"]))
        results.append(df.HistoBoost("MET_uncorr_phi", [self.axis_MET_phi], ["MET_uncorr_phi", "nominal_weight"]))
        results.append(df.HistoBoost("METx_uncorr", [self.axis_MET_xy], ["METx_uncorr", "nominal_weight"]))
        results.append(df.HistoBoost("METy_uncorr", [self.axis_MET_xy], ["METy_uncorr", "nominal_weight"]))
        
        results.append(df.HistoBoost("MET_corr_lep_pt", [self.axis_MET_pt], ["MET_corr_lep_pt", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_lep_phi", [self.axis_MET_phi], ["MET_corr_lep_phi", "nominal_weight"]))
        results.append(df.HistoBoost("METx_corr_lep", [self.axis_MET_xy], ["METx_corr_lep", "nominal_weight"]))
        results.append(df.HistoBoost("METy_corr_lep", [self.axis_MET_xy], ["METy_corr_lep", "nominal_weight"]))
        
        results.append(df.HistoBoost("MET_corr_xy_pt", [self.axis_MET_pt], ["MET_corr_xy_pt", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_xy_phi", [self.axis_MET_phi], ["MET_corr_xy_phi", "nominal_weight"]))
        results.append(df.HistoBoost("METx_corr_xy", [self.axis_MET_xy], ["METx_corr_xy", "nominal_weight"]))
        results.append(df.HistoBoost("METy_corr_xy", [self.axis_MET_xy], ["METy_corr_xy", "nominal_weight"]))
        
        results.append(df.HistoBoost("MET_corr_xy_pt_unclUp", [self.axis_MET_pt], ["MET_corr_xy_pt_unclUp", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_xy_phi_unclUp", [self.axis_MET_phi], ["MET_corr_xy_phi_unclUp", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_xy_pt_unclDw", [self.axis_MET_pt], ["MET_corr_xy_pt_unclDw", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_xy_phi_unclDw", [self.axis_MET_phi], ["MET_corr_xy_phi_unclDw", "nominal_weight"]))
        
        results.append(df.HistoBoost("npv", [self.axis_npv], ["PV_npvs", "nominal_weight"]))
        results.append(df.HistoBoost("RawMET_sumEt", [self.axis_RawMET_sumEt], ["RawMET_sumEt", "nominal_weight"]))
        
        # MET unclustered correction
        
        
            
        # histograms as function of npv, to derive/closure test the XY correction
        results.append(df.HistoBoost("METx_corr_lep_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METx_corr_lep", "nominal_weight"]))
        results.append(df.HistoBoost("METy_corr_lep_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METy_corr_lep", "nominal_weight"]))
            
        results.append(df.HistoBoost("METx_corr_xy_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METx_corr_xy", "nominal_weight"]))
        results.append(df.HistoBoost("METy_corr_xy_npv", [self.axis_npv, self.axis_MET_xy], ["PV_npvs", "METy_corr_xy", "nominal_weight"]))
        
        return df    
    
   
    def recoil_setup_Z(self, df, results, dataset_name):
        '''
        Setup the uncorrected recoil components
        '''
        
        # construct the 2D vector sum of the two leptons
        df = df.Define("Lep1_mom2", "ROOT::Math::Polar2DVectorD(%s[0], %s[0])" % (self.leptons_pt, self.leptons_phi))
        df = df.Define("Lep2_mom2", "ROOT::Math::Polar2DVectorD(%s[1], %s[1])" % (self.leptons_pt, self.leptons_phi))
        df = df.Define("Z_mom2", "Lep1_mom2 + Lep2_mom2")
        df = df.Define("qT", "Z_mom2.R()")
        df = df.Define("qTbin", "wrem::getqTbin(qT)")
        
        # uncorrected recoil
        df = df.Define("recoil_uncorr", "wrem::recoilComponents(MET_uncorr_pt, MET_uncorr_phi, qT, Z_mom2.Phi())")
        df = df.Define("recoil_uncorr_magn", "recoil_uncorr[0]")
        df = df.Define("recoil_uncorr_para", "recoil_uncorr[1]")
        df = df.Define("recoil_uncorr_para_qT", "recoil_uncorr[1] + qT")
        df = df.Define("recoil_uncorr_perp", "recoil_uncorr[2]")
        df = df.Define("recoil_uncorr_para_qT_perp", "recoil_uncorr[1] + qT + recoil_uncorr[2]")
        
        # lep corrected recoil
        df = df.Define("recoil_corr_lep", "wrem::recoilComponents(MET_corr_lep_pt, MET_corr_lep_phi, qT, Z_mom2.Phi())")
        df = df.Define("recoil_corr_lep_magn", "recoil_corr_lep[0]")
        df = df.Define("recoil_corr_lep_para", "recoil_corr_lep[1]")
        df = df.Define("recoil_corr_lep_para_qT", "recoil_corr_lep[1] + qT")
        df = df.Define("recoil_corr_lep_perp", "recoil_corr_lep[2]")
        df = df.Define("recoil_corr_lep_para_qT_perp", "recoil_corr_lep[1] + qT + recoil_corr_lep[2]")
        
        # MET XY corrected recoil
        df = df.Define("recoil_corr_xy", "wrem::recoilComponents(MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())")
        df = df.Define("recoil_corr_xy_magn", "recoil_corr_xy[0]")
        df = df.Define("recoil_corr_xy_para", "recoil_corr_xy[1]")
        df = df.Define("recoil_corr_xy_para_qT", "recoil_corr_xy[1] + qT")
        df = df.Define("recoil_corr_xy_perp", "recoil_corr_xy[2]")
        df = df.Define("recoil_corr_xy_para_qT_perp", "recoil_corr_xy[1] + qT + recoil_corr_xy[2]")
        
   
        results.append(df.HistoBoost("recoil_uncorr_magn", [self.axis_recoil_magn], ["recoil_uncorr_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_para", [self.axis_recoil_para], ["recoil_uncorr_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_para_qT", [self.axis_recoil_para_qT], ["recoil_uncorr_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_perp", [self.axis_recoil_perp], ["recoil_uncorr_perp", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_para_qT_perp", [self.axis_recoil_para_qT], ["recoil_uncorr_para_qT_perp", "nominal_weight"]))
            
        results.append(df.HistoBoost("recoil_corr_lep_magn", [self.axis_recoil_magn], ["recoil_corr_lep_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_para", [self.axis_recoil_para], ["recoil_corr_lep_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_lep_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_perp", [self.axis_recoil_perp], ["recoil_corr_lep_perp", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_para_qT_perp", [self.axis_recoil_para_qT], ["recoil_corr_lep_para_qT_perp", "nominal_weight"]))
            
        results.append(df.HistoBoost("recoil_corr_xy_magn", [self.axis_recoil_magn], ["recoil_corr_xy_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para", [self.axis_recoil_para], ["recoil_corr_xy_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_perp", [self.axis_recoil_perp], ["recoil_corr_xy_perp", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_perp", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT_perp", "nominal_weight"]))
        
        # for the correction, binned in qT
        results.append(df.HistoBoost("recoil_uncorr_magn_qTbinned", [self.axis_qTbinned, self.axis_recoil_magn_fine], ["qT", "recoil_uncorr_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_para_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_fine], ["qT", "recoil_uncorr_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_para_qT_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_qT_fine], ["qT", "recoil_uncorr_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_uncorr_perp_qTbinned", [self.axis_qTbinned, self.axis_recoil_perp_fine], ["qT", "recoil_uncorr_perp", "nominal_weight"]))
        
        results.append(df.HistoBoost("recoil_corr_lep_magn_qTbinned", [self.axis_qTbinned, self.axis_recoil_magn_fine], ["qT", "recoil_corr_lep_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_para_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_fine], ["qT", "recoil_corr_lep_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_para_qT_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_qT_fine], ["qT", "recoil_corr_lep_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_lep_perp_qTbinned", [self.axis_qTbinned, self.axis_recoil_perp_fine], ["qT", "recoil_corr_lep_perp", "nominal_weight"]))
        
        results.append(df.HistoBoost("recoil_corr_xy_magn_qTbinned", [self.axis_qTbinned, self.axis_recoil_magn_fine], ["qT", "recoil_corr_xy_magn", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_fine], ["qT", "recoil_corr_xy_para", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_qT_fine], ["qT", "recoil_corr_xy_para_qT", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_perp_qTbinned", [self.axis_qTbinned, self.axis_recoil_perp_fine], ["qT", "recoil_corr_xy_perp", "nominal_weight"]))
            
        results.append(df.HistoBoost("qT", [self.axis_qT], ["qT", "nominal_weight"]))

        return df


    def recoil_setup_gen(self, df, results, dataset):
    
        

        if (dataset.name == "DYmumu" and self.flavor == "mumu") or (dataset.name == "DYee" and self.flavor == "ee"):
            df = df.Define("recoil_corr_xy_gen", "wrem::recoilComponentsGen(MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi(), phiVgen)")
        elif (dataset.name in ["WplusJetsToMuNu", "WminusJetsToMuNu"] and self.flavor == "mu") or (dataset.name in ["WplusJetsToENu", "WminusJetsToENu"] and self.flavor == "e"):
            df = df.Define("recoil_corr_xy_gen", "wrem::recoilComponentsGen(MET_corr_xy_pt, MET_corr_xy_phi, Lep_pt, Lep_phi, phiVgen)")
        else: return df
        
        # decompose MET and dilepton (for Z) or lepton (for W) along the generator boson direction
        df = df.Define("qT_gen", "ptVgen") # pre-fsr defines should be loaded
        df = df.Define("qTbin_gen", "wrem::getqTbin(qT_gen)") 
        
        df = df.Define("recoil_corr_xy_magn_gen", "recoil_corr_xy_gen[0]")
        df = df.Define("recoil_corr_xy_para_gen", "recoil_corr_xy_gen[1]")
        df = df.Define("recoil_corr_xy_para_qT_gen", "recoil_corr_xy_gen[1] + qT_gen")
        df = df.Define("recoil_corr_xy_perp_gen", "recoil_corr_xy_gen[2]")
        df = df.Define("recoil_corr_xy_para_qT_perp_gen", "recoil_corr_xy_gen[1] + qT_gen + recoil_corr_xy_gen[2]")
        
        results.append(df.HistoBoost("recoil_corr_xy_magn_gen", [self.axis_recoil_magn], ["recoil_corr_xy_magn_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_gen", [self.axis_recoil_para], ["recoil_corr_xy_para_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_gen", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_perp_gen", [self.axis_recoil_perp], ["recoil_corr_xy_perp_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_perp_gen", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT_perp_gen", "nominal_weight"]))
        
        results.append(df.HistoBoost("recoil_corr_xy_magn_qTbinned_gen", [self.axis_qTbinned, self.axis_recoil_magn_fine], ["qT_gen", "recoil_corr_xy_magn_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qTbinned_gen", [self.axis_qTbinned, self.axis_recoil_para_fine], ["qT_gen", "recoil_corr_xy_para_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_qTbinned_gen", [self.axis_qTbinned, self.axis_recoil_para_qT_fine], ["qT_gen", "recoil_corr_xy_para_qT_gen", "nominal_weight"]))
        results.append(df.HistoBoost("recoil_corr_xy_perp_qTbinned_gen", [self.axis_qTbinned, self.axis_recoil_perp_fine], ["qT_gen", "recoil_corr_xy_perp_gen", "nominal_weight"]))
        
        results.append(df.HistoBoost("qT_gen", [self.axis_qT], ["qT_gen", "nominal_weight"]))
        
        return df

    def recoil_apply_Z(self, df, results, dataset_name, datasets_to_apply): 

        if dataset_name in datasets_to_apply:
            
            if self.parametric: df = df.Define("recoil_corr_rec", "wrem::recoilCorrectionParametric(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT)")
            else: df = df.Define("recoil_corr_rec", "wrem::recoilCorrectionBinned(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT)") # 
            df = df.Define("recoil_corr_rec_magn", "recoil_corr_rec[0]")
            df = df.Define("recoil_corr_rec_para", "recoil_corr_rec[1]")
            df = df.Define("recoil_corr_rec_para_qT", "recoil_corr_rec[1] + qT")
            df = df.Define("recoil_corr_rec_perp", "recoil_corr_rec[2]")
            df = df.Define("recoil_corr_rec_para_qT_perp", "recoil_corr_rec[1] + qT + recoil_corr_rec[2]")
            
            df = df.Define("MET_corr_rec", "wrem::METCorrection(MET_corr_xy_pt, MET_corr_xy_phi, recoil_corr_rec_para, recoil_corr_rec_perp, qT, Z_mom2.Phi()) ")
            df = df.Define("MET_corr_rec_pt", "MET_corr_rec[0]")
            df = df.Define("MET_corr_rec_phi", "MET_corr_rec[1]")
            
        else:
            df = df.Alias("recoil_corr_rec_magn", "recoil_corr_xy_magn")
            df = df.Alias("recoil_corr_rec_para", "recoil_corr_xy_para")
            df = df.Define("recoil_corr_rec_para_qT", "recoil_corr_xy_para + qT")
            df = df.Alias("recoil_corr_rec_perp", "recoil_corr_xy_perp")
            df = df.Alias("MET_corr_rec_pt", "MET_corr_xy_pt")
            df = df.Alias("MET_corr_rec_phi", "MET_corr_xy_phi")
            df = df.Alias("recoil_corr_rec_para_qT_perp", "recoil_corr_xy_para_qT_perp")
        
        #df = df.Filter("recoil_corr_rec_para_qT < 40 && recoil_corr_rec_para_qT > -40")
        #df = df.Filter("recoil_corr_rec_perp < 40 && recoil_corr_rec_perp > -40")
        
        #df = df.Filter("recoil_corr_rec_para_qT < -3 || recoil_corr_rec_para_qT > 3")
        #df = df.Filter("recoil_corr_rec_perp < -1 || recoil_corr_rec_perp > 1")
        #df = df.Filter("recoil_corr_rec_perp < 1 && recoil_corr_rec_perp > -1")
        df = df.Define("MET_corr_rec_xy_dPhi", "wrem::deltaPhi(MET_corr_rec_phi,MET_corr_xy_phi)")
        #df = df.Filter("MET_corr_rec_xy_dPhi < 3.14/2 && MET_corr_rec_xy_dPhi > -3.14/2")
        #df = df.Filter("MET_corr_rec_xy_dPhi > 3.14/2 || MET_corr_rec_xy_dPhi < -3.14/2")
        
        #df = df.Filter("qT < 50")
        #df = df.Filter("qT < 8")
        #df = df.Filter("qT > 7")
        
        if dataset_name in datasets_to_apply:
            df = df.Define("qTweight_", "wrem::qTweight(qTbin)")
            #df = df.Define("qTweight_", "1.0")
        else: 
            df = df.Define("qTweight_", "1.0")
        df = df.Define("weight_qTreweighted", "nominal_weight*qTweight_")
        
        df = df.Define("METx_corr_rec", "MET_corr_rec_pt*cos(MET_corr_rec_phi)")
        df = df.Define("METy_corr_rec", "MET_corr_rec_pt*sin(MET_corr_rec_phi)")
        
        results.append(df.HistoBoost("qT_qTreweighted", [self.axis_qT], ["qT", "weight_qTreweighted"]))
        results.append(df.HistoBoost("qT_qTreweighted_cut", [self.axis_qT], ["qT", "weight_qTreweighted"]))

        results.append(df.HistoBoost("recoil_corr_rec_magn", [self.axis_recoil_magn], ["recoil_corr_rec_magn", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para", [self.axis_recoil_para], ["recoil_corr_rec_para", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para_qT", [self.axis_recoil_para_qT], ["recoil_corr_rec_para_qT", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_perp", [self.axis_recoil_perp], ["recoil_corr_rec_perp", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_rec_pt", [self.axis_MET_pt], ["MET_corr_rec_pt", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_rec_pt_qT", [self.axis_qT, self.axis_MET_pt], ["qT", "MET_corr_rec_pt", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_rec_phi", [self.axis_MET_phi], ["MET_corr_rec_phi", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para_qT_perp", [self.axis_recoil_para_qT], ["recoil_corr_rec_para_qT_perp", "weight_qTreweighted"]))
        
        df = df.Define("recoil_corr_rec_para_qT_abs", "abs(recoil_corr_rec_para_qT)")
        df = df.Define("recoil_corr_rec_perp_abs", "abs(recoil_corr_rec_perp)")
        results.append(df.HistoBoost("recoil_corr_rec_para_qT_abs", [self.axis_recoil_para_perp_abs], ["recoil_corr_rec_para_qT_abs", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_perp_abs", [self.axis_recoil_para_perp_abs], ["recoil_corr_rec_perp_abs", "weight_qTreweighted"]))
        
        df = df.Define("recoil_corr_rec_para_plus_qt_squared", "(recoil_corr_rec_para + qT)*(recoil_corr_rec_para + qT)")
        df = df.Define("recoil_corr_rec_para_qT_squared", "recoil_corr_rec_para_qT*recoil_corr_rec_para_qT")
        df = df.Define("recoil_corr_rec_perp_squared", "recoil_corr_rec_perp*recoil_corr_rec_perp")
        results.append(df.HistoBoost("recoil_corr_rec_para_plus_qt_squared", [self.axis_recoil_para_perp_squared], ["recoil_corr_rec_para_plus_qt_squared", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para_qT_squared", [self.axis_recoil_para_perp_squared], ["recoil_corr_rec_para_qT_squared", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_perp_squared", [self.axis_recoil_para_perp_squared], ["recoil_corr_rec_perp_squared", "weight_qTreweighted"]))
        
        results.append(df.HistoBoost("recoil_corr_rec_magn_qTbinned", [self.axis_qTbinned, self.axis_recoil_magn], ["qT", "recoil_corr_rec_magn", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para_qTbinned", [self.axis_qTbinned, self.axis_recoil_para], ["qT", "recoil_corr_rec_para", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_para_qT_qTbinned", [self.axis_qTbinned, self.axis_recoil_para_qT], ["qT", "recoil_corr_rec_para_qT", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_rec_perp_qTbinned", [self.axis_qTbinned, self.axis_recoil_perp], ["qT", "recoil_corr_rec_perp", "weight_qTreweighted"]))
        
        
        results.append(df.HistoBoost("recoil_corr_rec_para_perp_2d", [self.axis_recoil_2d_para, self.axis_recoil_2d_perp], ["recoil_corr_rec_para_qT", "recoil_corr_rec_perp", "weight_qTreweighted"]))
        
        results.append(df.HistoBoost("met_recoil_uncorr", [self.axis_MET_pt, self.axis_recoil_magn], ["MET_corr_xy_pt", "recoil_corr_xy_magn", "weight_qTreweighted"]))
        results.append(df.HistoBoost("met_recoil_corr_rec", [self.axis_MET_pt, self.axis_recoil_magn], ["MET_corr_rec_pt", "recoil_corr_rec_magn", "weight_qTreweighted"]))
        
        results.append(df.HistoBoost("METx_corr_rec", [self.axis_MET_xy], ["METx_corr_rec", "weight_qTreweighted"]))
        results.append(df.HistoBoost("METy_corr_rec", [self.axis_MET_xy], ["METy_corr_rec", "weight_qTreweighted"]))
        
        
        df = df.Define("MET_corr_xy_ll_dPhi", "wrem::deltaPhi(MET_corr_xy_phi,Z_mom2.Phi())")
        df = df.Define("MET_corr_rec_ll_dPhi", "wrem::deltaPhi(MET_corr_rec_phi,Z_mom2.Phi())")
        #df = df.Define("MET_corr_rec_xy_dPhi", "wrem::deltaPhi(MET_corr_rec_phi,MET_corr_xy_phi)")
        df = df.Define("ll_phi", "Z_mom2.Phi()")
        results.append(df.HistoBoost("MET_corr_xy_ll_dPhi", [self.axis_MET_phi], ["MET_corr_xy_ll_dPhi", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_rec_ll_dPhi", [self.axis_MET_phi], ["MET_corr_rec_ll_dPhi", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_rec_xy_dPhi", [self.axis_MET_dphi], ["MET_corr_rec_xy_dPhi", "weight_qTreweighted"]))
        results.append(df.HistoBoost("ll_phi", [self.axis_MET_phi], ["ll_phi", "weight_qTreweighted"]))
        
        results.append(df.HistoBoost("recoil_corr_xy_magn_cut", [self.axis_recoil_magn], ["recoil_corr_xy_magn", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_cut", [self.axis_recoil_para], ["recoil_corr_xy_para", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_cut", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_xy_perp_cut", [self.axis_recoil_perp], ["recoil_corr_xy_perp", "weight_qTreweighted"]))
        results.append(df.HistoBoost("recoil_corr_xy_para_qT_perp_cut", [self.axis_recoil_para_qT], ["recoil_corr_xy_para_qT_perp", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_xy_pt_cut", [self.axis_MET_pt], ["MET_corr_xy_pt", "weight_qTreweighted"]))
        results.append(df.HistoBoost("MET_corr_xy_phi_cut", [self.axis_MET_phi], ["MET_corr_xy_phi", "weight_qTreweighted"]))

        '''
        if dataset_name in datasets_to_apply:
            
            df = df.Define("recoil_corr_rec_cut", "wrem::recoilCorrectionBinnedCut(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT, Lep_pt, Lep_eta, Lep_phi, MET_corr_xy_phi, Z_mom2.Phi(), MET_corr_xy_pt)")
            df = df.Define("recoil_corr_rec_magn_cut", "recoil_corr_rec_cut[0]")
            df = df.Define("recoil_corr_rec_para_cut", "recoil_corr_rec_cut[1]")
            df = df.Define("recoil_corr_rec_para_qT_cut", "recoil_corr_rec_cut[1] + qT")
            df = df.Define("recoil_corr_rec_perp_cut", "recoil_corr_rec_cut[2]")
            df = df.Define("recoil_corr_rec_para_qT_perp_cut", "recoil_corr_rec_cut[1] + qT + recoil_corr_rec_cut[2]")
            results.append(df.HistoBoost("test", [self.axis_MET_xy], ["recoil_corr_rec_magn_cut", "nominal_weight"]))

        df = df.Define("counts", "1")
        results.append(df.HistoBoost("counts", [self.axis_MET_xy], ["counts"]))
        results.append(df.HistoBoost("counts_weighted", [self.axis_MET_xy], ["counts", "nominal_weight"]))
        '''
        
        
        
        '''
        df = df.Define("METx_rec_corr", "MET_corr_rec_pt*cos(MET_corr_rec_phi)")
        df = df.Define("METy_rec_corr", "MET_corr_rec_pt*sin(MET_corr_rec_phi)")
        
        results.append(df.HistoBoost("METx_rec_corr", [self.axis_MET_xy], ["METx_rec_corr", "nominal_weight"]))
        results.append(df.HistoBoost("METy_rec_corr", [self.axis_MET_xy], ["METy_rec_corr", "nominal_weight"]))
        
        results.append(df.HistoBoost("METx_rec_corr_qTbinned", [self.axis_qTbinned, self.axis_MET_xy], ["qT", "METx_rec_corr", "nominal_weight"]))
        results.append(df.HistoBoost("METy_rec_corr_qTbinned", [self.axis_qTbinned, self.axis_MET_xy], ["qT", "METy_rec_corr", "nominal_weight"]))
        
        
        if dataset_name == "singlemuon": df = df.Define("met_xy_corr", "wrem::recoilXYCorrection(MET_corr_rec_pt, MET_corr_rec_phi, qTbin, 1)")
        else: df = df.Define("met_xy_corr", "wrem::recoilXYCorrection(MET_corr_rec_pt, MET_corr_rec_phi, qTbin, 0)")
        df = df.Define("met_xy_corr_pt", "met_xy_corr[0]")
        df = df.Define("met_xy_corr_phi", "met_xy_corr[1]")
        
        results.append(df.HistoBoost("MET_corr_rec_xy_pt", [self.axis_MET_pt], ["met_xy_corr_pt", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_rec_xy_phi", [self.axis_MET_phi], ["met_xy_corr_phi", "nominal_weight"]))
        
        df = df.Define("METx_rec_xy_corr", "met_xy_corr_pt*cos(met_xy_corr_phi)")
        df = df.Define("METy_rec_xy_corr", "met_xy_corr_pt*sin(met_xy_corr_phi)")
        
        results.append(df.HistoBoost("METx_rec_xy_corr", [self.axis_MET_xy], ["METx_rec_xy_corr", "nominal_weight"]))
        results.append(df.HistoBoost("METy_rec_xy_corr", [self.axis_MET_xy], ["METy_rec_xy_corr", "nominal_weight"]))
        '''
    
        return df
        
    def recoil_W(self, df, results, dataset, datasets_to_apply, auxPlots=False): 
    
        if dataset.name in datasets_to_apply:

            df = df.Define("recoil_corr_wz", "wrem::recoilCorrectionBinnedWtoZ(Lep_charge, recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen)")
            df = df.Define("recoil_corr_wz_magn", "recoil_corr_wz[0]")
            df = df.Define("recoil_corr_wz_para", "recoil_corr_wz[1]")
            df = df.Define("recoil_corr_wz_para_qT", "recoil_corr_wz[1] + qT_gen")
            df = df.Define("recoil_corr_wz_perp", "recoil_corr_wz[2]")
            df = df.Define("recoil_corr_wz_para_qT_perp", "recoil_corr_wz[1] + qT_gen + recoil_corr_wz[2]")

            df = df.Define("MET_corr_wz", "wrem::METCorrectionGen(recoil_corr_wz_para, recoil_corr_wz_perp, Lep_pt, Lep_phi, phiVgen) ") # was qTgen
            df = df.Define("MET_corr_wz_pt", "MET_corr_wz[0]")
            df = df.Define("MET_corr_wz_phi", "MET_corr_wz[1]")
            
            results.append(df.HistoBoost("recoil_corr_wz_magn", [self.axis_recoil_magn], ["recoil_corr_wz_magn", "nominal_weight"]))
            results.append(df.HistoBoost("recoil_corr_wz_para", [self.axis_recoil_para], ["recoil_corr_wz_para", "nominal_weight"]))
            results.append(df.HistoBoost("recoil_corr_wz_para_qT", [self.axis_recoil_perp], ["recoil_corr_wz_para_qT", "nominal_weight"]))
            results.append(df.HistoBoost("recoil_corr_wz_perp", [self.axis_recoil_perp], ["recoil_corr_wz_perp", "nominal_weight"]))
            
            
            if self.parametric: df = df.Define("recoil_corr_rec", "wrem::recoilCorrectionParametric(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen)")
            else: df = df.Define("recoil_corr_rec", "wrem::recoilCorrectionBinned(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen)")
            df = df.Define("recoil_corr_rec_magn", "recoil_corr_rec[0]")
            df = df.Define("recoil_corr_rec_para", "recoil_corr_rec[1]")
            df = df.Define("recoil_corr_rec_para_qT", "recoil_corr_rec[1] + qT_gen")
            df = df.Define("recoil_corr_rec_perp", "recoil_corr_rec[2]")
            df = df.Define("recoil_corr_rec_para_qT_perp", "recoil_corr_rec[1] + qT_gen + recoil_corr_rec[2]")
            
            df = df.Define("MET_corr_rec", "wrem::METCorrectionGen(recoil_corr_rec_para, recoil_corr_rec_perp, Lep_pt, Lep_phi, phiVgen) ") 
            df = df.Define("MET_corr_rec_pt", "MET_corr_rec[0]")
            df = df.Define("MET_corr_rec_phi", "MET_corr_rec[1]")
            
        else:
            # recoil not defined for backgrounds, only MET
            #df = df.Alias("recoil_corr_wz_magn", "recoil_corr_xy_magn")
            #df = df.Alias("recoil_corr_wz_para", "recoil_corr_xy_para")
            #df = df.Define("recoil_corr_wz_para_qT", "recoil_corr_xy_para + qT")
            #df = df.Alias("recoil_corr_wz_perp", "recoil_corr_xy_perp")
            #df = df.Alias("recoil_corr_wz_para_qT_perp", "recoil_corr_xy_para_qT_perp")
            
            df = df.Alias("MET_corr_wz_pt", "MET_corr_xy_pt")
            df = df.Alias("MET_corr_wz_phi", "MET_corr_xy_phi")
            
            df = df.Define("MET_corr_rec_pt", "MET_corr_xy_pt")
            df = df.Define("MET_corr_rec_phi", "MET_corr_xy_phi")
        
        
        results.append(df.HistoBoost("MET_corr_wz_pt", [self.axis_MET_pt], ["MET_corr_wz_pt", "nominal_weight"]))
        results.append(df.HistoBoost("MET_corr_wz_phi", [self.axis_MET_phi], ["MET_corr_wz_phi", "nominal_weight"]))
        
        return df        
    

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
        df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(self.recoil_qTbins)))
        recoil_vars = [(1,2), (1,3), (1,4),   (2,2), (2,3),   (3,2), (3,3), (3,4),    (4,2), (4,3)]
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),(3,6),(3,7),  (4,2),(4,3),(4,4),(4,5)] # 
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8),  (4,2),(4,3),(4,4),(4,5)] # 
        
        
        # 2 gauss for data = 4, 3 gauss for MC
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
        #recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
        for k in recoil_vars:

            # perturbation for current qTbin
            df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT, %d, %d)" % (k[0], k[1]))
           
            df = df.Define("recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_magn_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            results.append(df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            
            #df = df.Define("recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            #results.append(df.HistoBoost("gen_reco_para_qT_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            results.append(df.HistoBoost("recoil_magn_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_magn, axis_recoil_magn_pert, axis_recoil_stat_unc], ["recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            
            df = df.Define("recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_qT_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, qT)" % (k[0], k[1]))
            results.append(df.HistoBoost("recoil_para_qT_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para_qT, axis_recoil_para_qT_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para_qT", "recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            df = df.Define("recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_perp_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
            results.append(df.HistoBoost("recoil_perp_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_perp, axis_recoil_perp_pert, axis_recoil_stat_unc], ["recoil_corr_rec_perp", "recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            
            df = df.Define("recoil_corr_MET_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
            results.append(df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc], ["MET_corr_rec_pt", "recoil_corr_MET_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
            
            df = df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_lep_pt, MET_corr_lep_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
            
            
            df = df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_StatUnc(recoil_corr_MET_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (k[0], k[1], k[0], k[1]))
            results.append(df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
        
        return df
        

    def recoil_Z_unc_lowPU(self, df, results, axis_recoil_gen, axis_recoil_reco, axis_mt, axis_mll):
    
        return df
        if self.parametric:
        
            recoil_param_nStat = getattr(ROOT.wrem, "recoil_param_nStat")
            for item in recoil_param_nStat:
                tag = item.first
                nVars = recoil_param_nStat[tag]

                axis_recoil_unc = hist.axis.Regular(nVars, 0, nVars, underflow=False, overflow=False, name = "axis_recoil_unc_%s" % tag)
                df = df.Define("recoil_corr_rec_syst_%s" % tag, "wrem::recoilCorrectionParametricUnc(recoil_corr_xy_para, recoil_corr_xy_perp, qT, \"%s\")" % tag)
                df = df.Define("recoil_corr_rec_systIdx_%s" % tag, "wrem::indices_(%d, 0)" % nVars)
                    
                # MET
                df = df.Define("recoil_corr_MET_pt_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_pt_unc(recoil_corr_rec_syst_%s, qT, Z_mom2.Phi())" % tag)
                df = df.Define("recoil_corr_MET_phi_syst_%s" % tag, "wrem::recoilCorrectionParametric_MET_phi_unc(recoil_corr_rec_syst_%s, qT, Z_mom2.Phi())" % tag)
                results.append(df.HistoBoost("MET_recoilSyst_%s" % tag, [self.axis_MET_pt, axis_recoil_unc], ["recoil_corr_MET_pt_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
                
                
                # recoil
                df = df.Define("recoil_corr_para_qT_syst_%s" % tag, "wrem::recoilCorrectionParametric_para_qT_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_para_syst_%s" % tag, "wrem::recoilCorrectionParametric_para_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_perp_syst_%s" % tag, "wrem::recoilCorrectionParametric_perp_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                df = df.Define("recoil_corr_magn_syst_%s" % tag, "wrem::recoilCorrectionParametric_magn_unc(recoil_corr_rec_syst_%s, qT)"% tag)
                results.append(df.HistoBoost("recoil_para_qT_recoilSyst_%s" % tag, [self.axis_recoil_para_qT, axis_recoil_unc], ["recoil_corr_para_qT_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_para_recoilSyst_%s" % tag, [self.axis_recoil_para, axis_recoil_unc], ["recoil_corr_para_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_perp_recoilSyst_%s" % tag, [self.axis_recoil_perp, axis_recoil_unc], ["recoil_corr_perp_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_magn_recoilSyst_%s" % tag, [self.axis_recoil_magn, axis_recoil_unc], ["recoil_corr_magn_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
                
                # mT
                df = df.Define("recoil_corr_mT_recoilSyst_%s" % tag, "wrem::recoilCorrectionParametric_mT_unc(recoil_corr_MET_pt_syst_%s, recoil_corr_MET_phi_syst_%s, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (tag, tag))
                results.append(df.HistoBoost("mT_corr_rec_recoilSyst_%s" % tag, [axis_mt, axis_recoil_unc], ["recoil_corr_mT_recoilSyst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "weight_qTreweighted"]))
            
                # Z peak
                #results.append(df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "weight_qTreweighted"]))
                
        else:

            axis_recoil_stat_unc = hist.axis.Regular(len(self.recoil_qTbins), 0, len(self.recoil_qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")
            
            # make copies of axes for perturbation
            axis_recoil_reco_pert = hist.axis.Variable(axis_recoil_reco.edges, name = "recoil_reco_pert", underflow=False, overflow=True)
            axis_mt_pert = hist.axis.Variable(axis_mt.edges, name = "mt_pert",underflow=False, overflow=True)
            axis_MET_pert = hist.axis.Variable(self.axis_MET_pt.edges, name = "recoil_MET_pt_pert", underflow=False, overflow=True)
            
            axis_recoil_para_qT_pert = hist.axis.Variable(self.axis_recoil_para_qT.edges, name = "recoil_para_qT_pert", underflow=False, overflow=False)
            axis_recoil_para_pert = hist.axis.Variable(self.axis_recoil_para.edges, name = "recoil_para_pert", underflow=False, overflow=False)
            axis_recoil_magn_pert = hist.axis.Variable(self.axis_recoil_magn.edges, name = "recoil_magn_pert", underflow=False)
            axis_recoil_perp_pert = hist.axis.Variable(self.axis_recoil_perp.edges, name = "recoil_perp_pert", underflow=False, overflow=False)
            
        
            # recoil stat uncertainty  (para_data/perp_data/para_mc/perp_mc, number of variations)
            df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(self.recoil_qTbins)))
            recoil_vars = [(1,2),(1,3),(1,4),(1,5),(1,6),  (2,2),(2,3),(2,4),(2,5),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)] # 2 gauss for data = 4, 3 gauss for MC
            for k in recoil_vars:

                # perturbation for current qTbin
                df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para, recoil_corr_xy_perp, qTbin, qT, %d, %d)" % (k[0], k[1]))
               
                # MET
                df = df.Define("recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
                df = df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, MET_corr_xy_pt, MET_corr_xy_phi, qT, Z_mom2.Phi())" % (k[0], k[1]))
                results.append(df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc], ["MET_corr_rec_pt", "recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))

                # recoil
                df = df.Define("recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_qT_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin, qT)" % (k[0], k[1]))
                df = df.Define("recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_para_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                df = df.Define("recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_perp_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                df = df.Define("recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_magn_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin)" % (k[0], k[1]))
                results.append(df.HistoBoost("recoil_para_qT_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para_qT, axis_recoil_para_qT_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para_qT", "recoil_corr_para_qT_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_para_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_rec_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_perp_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_perp, axis_recoil_perp_pert, axis_recoil_stat_unc], ["recoil_corr_rec_perp", "recoil_corr_perp_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))
                results.append(df.HistoBoost("recoil_magn_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_magn, axis_recoil_magn_pert, axis_recoil_stat_unc], ["recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))

                # mT
                df = df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_StatUnc(recoil_corr_MET_pt_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin, TrigMuon_pt, TrigMuon_phi, NonTrigMuon_pt, NonTrigMuon_phi)" % (k[0], k[1], k[0], k[1]))
                results.append(df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "nominal_weight"]))
            
                
                # Z peak
                #results.append(df.HistoBoost("gen_reco_magn_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_recoil_gen, axis_recoil_reco, axis_recoil_reco_pert, axis_mll, axis_recoil_stat_unc], ["ptVgen", "recoil_corr_rec_magn", "recoil_corr_magn_stat_unc_%d_%d" % (k[0], k[1]), "massZ", "recoil_corr_stat_idx", "weight_qTreweighted"]))
                #results.append(df.HistoBoost("gen_reco_para_qT_mll_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc], ["recoil_corr_para", "recoil_corr_para_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "weight_qTreweighted"]))

        return df

        

        
        
    def recoil_W_unc_lowPU(self, df, results, axis_charge, axis_mt, axis_eta, axis_passIso):
    
        if self.parametric:
        
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
                results.append(df.HistoBoost("MET_recoilSyst_%s" % tag, [self.axis_MET_pt, axis_recoil_unc, axis_eta, axis_charge, axis_passIso], ["recoil_corr_MET_pt_syst_%s" % tag, "recoil_corr_rec_systIdx_%s" % tag, "Lep_abs_eta", "Lep_charge", "passIso", "nominal_weight"]))
                
                # mT
                df = df.Define("recoil_corr_mT_recoilSyst_%s" % tag, "wrem::recoilCorrectionParametric_mT_2_unc(recoil_corr_MET_pt_syst_%s, recoil_corr_MET_phi_syst_%s, Lep_pt, Lep_phi)" % (tag, tag))
                results.append(df.HistoBoost("mT_corr_rec_recoilSyst_%s" % tag, [axis_mt, axis_eta, axis_charge, axis_passIso, axis_recoil_unc], ["recoil_corr_mT_recoilSyst_%s" % tag, "Lep_abs_eta", "Lep_charge", "passIso", "recoil_corr_rec_systIdx_%s" % tag, "nominal_weight"]))
           
        else:

            axis_recoil_stat_unc = hist.axis.Regular(len(self.recoil_qTbins), 0, len(self.recoil_qTbins), underflow=False, overflow=False, name = "recoil_stat_unc_var")
            axis_mt_pert = hist.axis.Variable(axis_mt.edges, name = "mt_pert", underflow=False, overflow=True)
            axis_MET_pert = hist.axis.Variable(self.axis_MET_pt.edges, name = "recoil_MET_pt_pert", underflow=False, overflow=True)
            
            # recoil stat uncertainty  (para_data/perp_data/para_mc/perp_mc, number of variations)
            df = df.Define("recoil_corr_stat_idx", "wrem::indices_(%d, 0)" % (len(self.recoil_qTbins)))

            # 2 gauss for data = 4, 3 gauss for MC
            recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
            for k in recoil_vars:

                # perturbation for current qTbin
                df = df.Define("recoil_corr_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned(recoil_corr_xy_para_gen, recoil_corr_xy_perp_gen, qTbin_gen, qT_gen, %d, %d)" % (k[0], k[1]))
               
                # MET
                df = df.Define("recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_pt_gen_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi, phiVgen)" % (k[0], k[1]))
                df = df.Define("recoil_corr_MET_phi_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_MET_phi_gen_StatUnc(recoil_corr_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi, phiVgen)" % (k[0], k[1]))
                results.append(df.HistoBoost("MET_recoilStatUnc_%d_%d" % (k[0], k[1]), [self.axis_MET_pt, axis_MET_pert, axis_recoil_stat_unc, axis_eta, axis_charge], ["MET_corr_rec_pt", "recoil_corr_MET_pt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "Lep_abs_eta", "Lep_charge", "nominal_weight"]))

                # mT
                df = df.Define("recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "wrem::recoilCorrectionBinned_mt_2_StatUnc(recoil_corr_MET_pt_stat_unc_%d_%d, recoil_corr_MET_phi_stat_unc_%d_%d, qTbin_gen, Lep_pt, Lep_phi)" % (k[0], k[1], k[0], k[1]))
                results.append(df.HistoBoost("mT_corr_rec_recoilStatUnc_%d_%d" % (k[0], k[1]), [axis_mt, axis_mt_pert, axis_recoil_stat_unc, axis_eta, axis_charge], ["mT_corr_rec", "recoil_corr_mt_stat_unc_%d_%d" % (k[0], k[1]), "recoil_corr_stat_idx", "Lep_abs_eta", "Lep_charge", "nominal_weight"]))
            
        
        return df
    
