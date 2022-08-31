#ifndef WREMNANTS_LOWPU_RECOIL_H
#define WREMNANTS_LOWPU_RECOIL_H


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include "defines.h"

namespace wrem {

struct recoil {
  float para;
  float perp;
};    
    
using Vec_recoil = ROOT::VecOps::RVec<recoil>;
    
TFile *recoil_fits_Z_param = new TFile("wremnants/data/lowPU/recoil_fits_Z_param_refit.root", "READ");
//TFile *recoil_fits_Z_param_refit = new TFile("wremnants/data/lowPU/recoil_fits_Z_param_refit.root", "READ");
//TFile *recoil_fits_Z = new TFile("wremnants/data/lowPU/recoil_fits_Z.root", "READ");
//TFile *recoil_fits_Z = new TFile("wremnants/data/lowPU/recoil_fits_Z.root", "READ");
std::map<std::string, TH3D*> recoil_hists;
std::map<std::string, TH3D*> recoil_hists_gen_W;
std::map<std::string, TH3D*> recoil_hists_gen_Z;
std::map<std::string, TH1D*> recoil_hists_param;

// recoil parametric maps
std::map<std::string, int> recoil_param_nGauss;


std::map<std::string, TF1*> recoil_param_funcs;
std::map<std::string, int> recoil_param_nStat;

std::map<std::string, int> recoil_binned_nGauss;
std::map<std::string, std::vector<float>> recoil_binned_qTbins;
std::map<std::string, std::vector<  std::vector<  std::vector< float > >  >> recoil_binned_components; // tag qTbin nGauss, component

//std::map<std::string, TF1*> recoil_param_funcs;
//std::map<std::string, std::map<std::string, TF1*>> recoil_param_funcs;
    
std::vector<float> qTbins;

//std::vector<float> qTweights{ 1.042, 0.977, 0.926, 0.914, 0.921, 0.949, 0.983, 0.973, 0.999, 1.025, 1.049, 1.005, 1.087, 1.048, 1.057, 1.068, 1.068, 1.098, 1.036, 1.049, 1.081, 1.042, 1.064, 1.052, 1.123, 1.049, 1.000, 1.028, 1.047, 1.005, 0.972, 1.032, 1.067, 1.033, 1.049, 1.004, 0.994, 1.056, 1.056, 1.013, 0.988, 0.985, 0.956, 0.986, 0.999, 1.053, 1.008, 1.032, 1.044, 0.920, 1.000, 0.921, 1.016, 0.905, 0.845 };
//std::vector<float> qTweights{   1.041, 0.977, 0.925, 0.913, 0.920, 0.948, 0.981, 0.972, 0.997, 1.024, 1.047, 1.004, 1.085, 1.046, 1.055, 1.066, 1.066, 1.095, 1.033, 1.047, 1.078, 1.039, 1.061, 1.049, 1.120, 1.045, 0.996, 1.024, 1.042, 1.000, 0.967, 1.026, 1.060, 1.025, 1.041, 0.995, 0.984, 1.046, 1.044, 1.000, 0.974, 0.968, 0.936, 0.962, 0.971, 1.021, 0.974, 1.000, 1.011, 0.889, 0.974, 0.897, 0.998, 0.891, 1 };
//std::vector<float> qTweights{  1.078, 1.029, 0.981, 0.973, 0.921, 0.929, 0.921, 0.905, 0.915, 0.925, 0.952, 0.944, 1.005, 0.956, 0.954, 0.991, 0.995, 1.000, 1.016, 1.032, 1.030, 1.065, 1.020, 0.986, 1.073, 1.097, 1.046, 1.046, 1.048, 1.062, 1.050, 1.082, 1.079, 1.052, 1.075, 1.117, 1.027, 1.040, 1.039, 1.055, 1.080, 1.075, 1.064, 1.012, 1.061, 1.061, 1.034, 1.064, 1.143, 1.095, 1.057, 1.034, 0.992, 1.001, 1.044, 1.003, 1.037, 1.048, 1.000, 1.000, 0.973, 0.961, 1.025, 1.028, 1.113, 1.004, 1.039, 1.011, 1.070, 1.010, 0.979, 1.011, 0.980, 0.989, 1.052, 1.039, 0.985, 1.108, 1.023, 0.975, 0.956, 0.952, 1.041, 0.916, 0.994, 0.926, 0.946, 0.968, 0.972, 0.932, 1.004, 0.932, 0.999, 1.047, 0.958, 0.992, 1.000, 1.011, 0.889, 0.973, 0.898, 0.997, 0.891, 1.000, 1.000 };

std::vector<float> qTweights { 1.041, 0.977, 0.925, 0.913, 0.920, 0.948, 0.981, 0.972, 0.997, 1.024, 1.047, 1.004, 1.085, 1.046, 1.055, 1.066, 1.066, 1.095, 1.033, 1.047, 1.078, 1.039, 1.061, 1.049, 1.120, 1.045, 0.996, 1.024, 1.042, 1.000, 0.973, 0.961, 1.025, 1.028, 1.113, 1.004, 1.039, 1.011, 1.070, 1.010, 0.979, 1.011, 0.980, 0.989, 1.051, 1.039, 0.985, 1.108, 1.023, 0.975, 0.956, 0.952, 1.041, 0.916, 0.994, 0.926, 0.946, 0.968, 0.972, 0.932, 1.004, 0.932, 0.999, 1.047, 0.958, 0.992, 1.000, 1.011, 0.889, 0.973, 0.898, 0.997, 0.891, 1.000, 1.000, 1, 1, 1, 1, 1  };

//std::vector<float> qTweights { 1.041, 0.977, 0.925, 0.913, 0.920, 0.948, 0.981, 0.972, 0.997, 1.024, 1.047, 1.004, 1.085, 1.046, 1.055, 1.066, 1.066, 1.095, 1.033, 1.047, 1.059, 1.055, 1.083, 1.010, 1.022, 0.967, 1.026, 1.060, 1.025, 1.041, 0.990, 1.026, 1.038, 0.963, 0.986, 0.968, 0.936, 0.962, 1.004, 0.932, 1.021, 0.974, 0.988, 0.926, 0.960, 1, 1, 1 };




std::vector<float> corr_x { -0.595, -0.660, -0.548, -0.732, -0.678, -0.643, -0.390, -0.725, -0.540, -0.461, -0.677, -0.755, -0.730, -0.317, -0.340, -0.782, -0.448, -0.261, -0.643, -0.649, -0.651, -0.598, -0.260, -0.975, -0.632, -0.655, -0.760, -0.664, -0.911, -0.806, -0.858, -0.980, -0.604, -0.298, -1.037, -0.924, -0.695, -0.568, -0.633, -0.992, -1.186, -0.389, -0.987, -0.593, -0.986, -0.008, -0.155, -0.177, -0.862, -1.745, -0.844, -0.607, -0.191, -1.213, -0.317, -0.787, -1.178, -0.900, -2.065, -1.141, -1.223, -1.454, 0.331, -0.534, -0.834, 0.763, -1.153, -0.515, -0.659, -0.055, -2.063, -1.534, -0.426, 0 };
std::vector<float> corr_y { 0.507, 0.243, 0.356, 0.229, 0.339, 0.327, 0.307, 0.382, 0.339, 0.539, 0.287, 0.281, 0.369, 0.551, 0.754, 0.009, 0.777, 0.610, 0.304, 0.336, 0.042, 0.542, 0.220, 0.550, 0.710, 0.409, 0.170, 0.138, 0.433, 0.263, 0.841, -0.268, 0.614, -0.086, 0.663, 0.293, -0.154, 1.577, 0.705, 1.028, 0.432, 0.418, 0.363, -0.001, 0.948, -0.070, 1.130, 1.444, 1.592, -0.414, 0.133, 1.050, 0.716, 1.217, 0.690, -0.025, -0.360, 0.799, -0.531, -0.327, 1.601, -0.434, -0.259, -0.576, 2.648, 0.039, 0.807, 1.053, -1.098, -0.495, 2.614, 0.960, 4.728, 0.449};


std::vector<float> met_xy_corr_x_data_nom, met_xy_corr_y_data_nom, met_xy_corr_x_mc_nom, met_xy_corr_y_mc_nom;


void recoil_init(char* name) {
    
    TH3D *h3;
    TH1D *h1;
    
    
    TFile *recoil_fits_Z = new TFile(name, "READ");
    
    /*
    for(TKey *key: ROOT::RangeStaticCast<TKey*>(*recoil_fits_Z->GetListOfKeys())) {
        
        auto *h = recoil_fits_Z->Get(key->GetName());
        if(strncmp(key->GetClassName(), "TH3D", 4) == 0) {
            cout << "Load recoil histo " << key->GetName() << endl;
            h3 = (TH3D*)(recoil_fits_Z->Get(key->GetName()));
            recoil_hists.insert({key->GetName(), h3});
        }
    }
    
    */
    

 
    h3 = (TH3D*)recoil_fits_Z->Get("mc_para");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("mc_perp");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("data_para");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("data_perp");
    recoil_hists.insert({h3->GetName(), h3});
/*
    h1 = (TH1D*)recoil_fits_Z_param->Get("para_mc");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("perp_mc");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("para_data");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("perp_data");
    recoil_hists_param.insert({h1->GetName(), h1});
*/
}


void recoil_init_gen(char* name_W, char* name_Z) {
    
    TH3D *h3;

    TFile *recoil_fits_W = new TFile(name_W, "READ");
    TH3D *gen_plus_para = (TH3D*)recoil_fits_W->Get("gen_plus_para");
    TH3D *gen_minus_para = (TH3D*)recoil_fits_W->Get("gen_minus_para");
    TH3D *gen_plus_perp = (TH3D*)recoil_fits_W->Get("gen_plus_perp");
    TH3D *gen_minus_perp = (TH3D*)recoil_fits_W->Get("gen_minus_perp");
    gen_plus_para->SetDirectory(0);
    gen_minus_para->SetDirectory(0);
    gen_plus_perp->SetDirectory(0);
    gen_minus_perp->SetDirectory(0);
    recoil_hists_gen_W.insert({gen_plus_para->GetName(), gen_plus_para});
    recoil_hists_gen_W.insert({gen_minus_para->GetName(), gen_minus_para});
    recoil_hists_gen_W.insert({gen_plus_perp->GetName(), gen_plus_perp});
    recoil_hists_gen_W.insert({gen_minus_perp->GetName(), gen_minus_perp});
    recoil_fits_W->Close();
    
    TFile *recoil_fits_Z = new TFile(name_Z, "READ");
    TH3D *gen_para_Z = (TH3D*)recoil_fits_Z->Get("gen_para");
    TH3D *gen_perp_Z = (TH3D*)recoil_fits_Z->Get("gen_perp");
    gen_para_Z->SetDirectory(0);
    gen_perp_Z->SetDirectory(0);
    recoil_hists_gen_Z.insert({gen_para_Z->GetName(), gen_para_Z});
    recoil_hists_gen_Z.insert({gen_perp_Z->GetName(), gen_perp_Z});
    recoil_fits_Z->Close();
}


int getqTbin(float qT) {

    for(unsigned int i=0; i < qTbins.size()-1; i++) {
        
        if(qT >= qTbins.at(i) and qT < qTbins.at(i+1)) return i;
    }
    
    //cout << "[getQtBin] Bin not found qT=" << qT << endl;
    return -1;
}

double qTweight(int qTbinIdx) {
    
    //return 1.;
    return qTweights.at(qTbinIdx);
}





Vec_f METLeptonCorrection(double MET_pt, double MET_phi, Vec_f lep_pt_uncorr, Vec_f lep_pt, Vec_f lep_phi) {

	// correct MET for muon scale corrections (in lab frame pT-phi)
	TVector2 lep1_raw(lep_pt_uncorr[0]*cos(lep_phi[0]), lep_pt_uncorr[0]*sin(lep_phi[0]));
	TVector2 lep2_raw(lep_pt_uncorr[1]*cos(lep_phi[1]), lep_pt_uncorr[1]*sin(lep_phi[1]));
    TVector2 lep1_corr(lep_pt[0]*cos(lep_phi[0]), lep_pt[0]*sin(lep_phi[0]));
    TVector2 lep2_corr(lep_pt[1]*cos(lep_phi[1]), lep_pt[1]*sin(lep_phi[1]));

	TVector2 MET(MET_pt*cos(MET_phi), MET_pt*sin(MET_phi));
	TVector2 MET_corr = MET + lep1_raw + lep2_raw - lep1_corr - lep2_corr;
	double MET_pt_corr = MET_corr.Mod();
    double MET_phi_corr = MET_corr.Phi_mpi_pi(MET_corr.Phi());
	
    Vec_f res(2, 0);
    res[0] = MET_pt_corr;
    res[1] = MET_phi_corr;
    
    //cout << "** Lepton angles ** " << MET_corr.Phi() << " " << MET_corr.Phi_mpi_pi(MET_corr.Phi()) << endl;
     

	
	return res;
	
}


Vec_f METLeptonCorrection(double MET_pt, double MET_phi, double lep_pt_uncorr, double lep_pt, double lep_phi) {

    Vec_f lep_pt_uncorr_(1, lep_pt_uncorr);
    Vec_f lep_pt_(1, lep_pt);
    Vec_f lep_phi_(1, lep_phi);
	return METLeptonCorrection(MET_pt, MET_phi, lep_pt_uncorr_, lep_pt_, lep_phi_);
}

Vec_f recoilComponents(double MET_pt, double MET_phi, double lep_pt, double lep_phi) {
    
    // computes the hadronic recoil, defined as -(MET + V), V being the vector sum of the (di)lepton system in the lab frame
    
    
    
    double pUx  = MET_pt*cos(MET_phi) + lep_pt*cos(lep_phi); // recoil x in lab frame
	double pUy  = MET_pt*sin(MET_phi) + lep_pt*sin(lep_phi); // recoil y in lab frame
	double pU   = std::hypot(pUx, pUy);
	double Ux	= - (pUx*cos(lep_phi) + pUy*sin(lep_phi));
	double Uy	= - (- pUx*sin(lep_phi) + pUy*cos(lep_phi));    

    Vec_f res(3, 0);
    res[0] = pU;
    res[1] = Ux;
    res[2] = Uy;
    
    //cout << " MET_phi=" << MET_phi << "  lep_phi=" << lep_phi<< "  pU=" << pU<< "  Ux=" << Ux << "  Uy=" << Uy << endl;
    
	return res;
}


Vec_f recoilComponentsGen(double MET_pt, double MET_phi, double lep_pt, double lep_phi, double gen_lep_phi) {
    
    // computes the hadronic recoil, defined as -(MET + V), V being the vector sum of the (di)lepton system in the lab frame
    double pUx  = MET_pt*cos(MET_phi) + lep_pt*cos(lep_phi); // recoil x in lab frame
	double pUy  = MET_pt*sin(MET_phi) + lep_pt*sin(lep_phi); // recoil y in lab frame
	double pU   = std::hypot(pUx, pUy);
	double Ux	= - (pUx*cos(gen_lep_phi) + pUy*sin(gen_lep_phi));
	double Uy	= - (- pUx*sin(gen_lep_phi) + pUy*cos(gen_lep_phi));

    Vec_f res(3, 0);
    res[0] = pU;
    res[1] = Ux;
    res[2] = Uy;
        
	return res;
}

Vec_f METXYCorrection(double MET_pt, double MET_phi, int npv, int isData) {
    
    // preserve the MET magnitude??
    
    double delta;
    double pUx = MET_pt*cos(MET_phi);
    double pUy = MET_pt*sin(MET_phi);

    if(isData == 1) {
    
        // x
        delta = 0;
        for(int k=0; k<met_xy_corr_x_data_nom.size(); k++) delta += met_xy_corr_x_data_nom.at(k)*std::pow(npv, k);
        pUx -= delta;
        
        // y
        delta = 0;
        for(int k=0; k<met_xy_corr_y_data_nom.size(); k++) delta += met_xy_corr_y_data_nom.at(k)*std::pow(npv, k);
        pUy -= delta;    
    }
    else {
        // x
        delta = 0;
        for(int k=0; k<met_xy_corr_x_mc_nom.size(); k++) delta += met_xy_corr_x_mc_nom.at(k)*std::pow(npv, k);
        pUx -= delta;
        
        // y
        delta = 0;
        for(int k=0; k<met_xy_corr_y_mc_nom.size(); k++) delta += met_xy_corr_y_mc_nom.at(k)*std::pow(npv, k);
        pUy -= delta;     
    }

    Vec_f res(2, 0);
    res[0] = hypot(pUx, pUy);
    res[1] = atan2(pUy, pUx);
    
    //cout << "** XY angles ** " << atan(pUy/pUx) << " " << res[1] << " (pUy=" << pUy << " pUx=" << pUx << ")" << endl;
    
    return res;
}


double recoilParameterization_sigma(double a, double b, double c, double qT) {
    
    return a*std::pow(qT + b, c);
}



class GaussianSum {
    
    private:
    
        std::vector<double> mean;
        std::vector<double> sigma;
        std::vector<double> norm;
        double totNorm;
        double extraNorm =1;
        
        //int MaxIterations = 512; 
        //double _tol = 2.2204460492503131e-16;
        
        int MaxIterations = 1024; 
        double _tol = 2.2204460492503131e-16;
        
        bool verbose = true;
    
    
    public:
    
        GaussianSum() {};
        
        void addTerm(double mean_, double sigma_, double norm_) {
            
            mean.push_back(mean_);
            sigma.push_back(sigma_);
            norm.push_back(norm_*extraNorm);
            totNorm = std::accumulate(norm.begin(), norm.end(), 0.0f); // recompute norm
        }
        
        
        void setNorm(double norm) {
            
            extraNorm = norm;
        }
        
        double evalCDF(double pval) {

            double ret = 0;
            double tmp;
            for(unsigned int i = 0; i < mean.size(); i++) {
                
                tmp = (pval-mean.at(i))/sigma.at(i);
                ret += norm.at(i)*std::erfc(-tmp/TMath::Sqrt2());
            }
            return 0.5*ret/totNorm;
        }
        
        double findRoot(double value, double xMin, double xMax) {
            
            double a(xMin), b(xMax);
            
            if(totNorm < 0.999 or totNorm > 1.001) cout << "[nGaussian::findRoot]: Total norm not equal to 1: " << totNorm  << endl;
            
            double fa = this->evalCDF(a) - value;
            double fb = this->evalCDF(b) - value;
            
     
            if(fb*fa > 0) {
                if(verbose) cout << "[nGaussian::findRoot]: initial interval does not bracket a root (fa=" << fa << ", fb=" << fb << ", pval=" << value << ", xMin=" << a << ", xMax=" << b <<  ")" << endl;
                return 999999;
            }

            Bool_t ac_equal(kFALSE);
            double fc = fb;
            double c(0),d(0),e(0);
            for(Int_t iter= 0; iter <= MaxIterations; iter++) {
         
                if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
                  
                    // Rename a,b,c and adjust bounding interval d
                    ac_equal = kTRUE;
                    c = a;
                    fc = fa;
                    d = b - a;
                    e = b - a;
                }
              
                if (fabs (fc) < fabs (fb)) {
                  
                    ac_equal = kTRUE;
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
             
                double tol = 0.5 * _tol * fabs(b);
                double m = 0.5 * (c - b);
             
             
                if (fb == 0 || fabs(m) <= tol) {
                  
                    //cout << "RooBrentRootFinder: iter = " << iter << " m = " << m << " tol = " << tol << endl ;
                    return b;
                }
              
                if (fabs (e) < tol || fabs (fa) <= fabs (fb)) {
                  
                    // Bounds decreasing too slowly: use bisection
                    d = m;
                    e = m;
                }
                else {
                  
                    // Attempt inverse cubic interpolation
                    double p, q, r;
                    double s = fb / fa;
                  
                    if (ac_equal) {
                    
                        p = 2 * m * s;
                        q = 1 - s;
                    }
                    else {
                        
                        q = fa / fc;
                        r = fb / fc;
                        p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                        q = (q - 1) * (r - 1) * (s - 1);
                    }
                  
                    // Check whether we are in bounds
                    if (p > 0) {
                        
                        q = -q;
                    }
                    else {
                        
                        p = -p;
                    }
                  
                    double min1= 3 * m * q - fabs (tol * q);
                    double min2= fabs (e * q);
                    
                    if (2 * p < (min1 < min2 ? min1 : min2)) {
               
                        // Accept the interpolation
                        e = d;
                        d = p / q;
                    }
                    else {
               
                        // Interpolation failed: use bisection.
                        d = m;
                        e = m;
                    }
                }
                
                // Move last best guess to a
                a = b;
                fa = fb;
                // Evaluate new trial root
                if (fabs (d) > tol) {
                  
                    b += d;
                }
                else {
                  
                    b += (m > 0 ? +tol : -tol);
                }
                
                fb = this->evalCDF(b) - value;
            }
            
            
            if(verbose) cout << "[nGaussian::findRoot]: maximum iterations exceeded" << endl;
            //return b; // return best-estimate of root
            return 999999; // return default value, reject this 
        }
    
};



// recoil correction
Vec_f recoilCorrectionBinned(double pU1, double pU2, double qTbinIdx, double qT, int corrType=0, int corrParam=1) {
    
    Vec_f res(3, 0);
    //cout << pU1 << " " << pU2 << " " << qTbinIdx << endl; 
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    //return res;
    
    TH3D *h3;
    int corrIdx = 1;
    
    double norm_para_data = 1;
    double norm_perp_data = 1;
    double norm_para_mc = 1;
    double norm_perp_mc = 1;
    
    //cout << "************* " << qT << endl;
    
    GaussianSum u1_data;
    h3 = recoil_hists["data_para"];
    corrIdx = 1;
    if(corrType == 1) corrIdx = corrParam;
    norm_para_data =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u1_data.setNorm(norm_para_data);
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
/*
    if(corrType == 1) {
        if(qTbinIdx == 40) {
            cout << corrIdx << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 4, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 5, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 6, corrIdx) << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 7, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 8, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 9, corrIdx) << endl;
        }
    }*/
    /*
    if(qT > 280) {
        cout << "************* " << qT << " " << qTbinIdx << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 4, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 5, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 6, corrIdx) << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 7, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 8, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 9, corrIdx) << endl;
    }
    */
    //if(qTbinIdx == 3) {
    //    cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << " " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << " " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
    //}
    
    
    
    GaussianSum u2_data;
    h3 = recoil_hists["data_perp"];
    corrIdx = 1;
    if(corrType == 2) corrIdx = corrParam;
    norm_perp_data =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u2_data.setNorm(norm_perp_data);
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
    /*
    // parametric
    
    GaussianSum u2_data;
    double mean = 0;
    double sigma1 = recoilParameterization_sigma(5.02232e+00, 1.39195e+00, 6.33123e-02, qT);
    double sigma2 = recoilParameterization_sigma(3.71235e-01, 1.07052e+02, 6.74092e-01, qT);
    double sigma3 = recoilParameterization_sigma(5.69628e-01, 1.19781e+02, 6.57171e-01, qT);
    double n1 = 0.3;
    double n2 = 0.6;
    double n3 = 1. - n1 - n2;
    u2_data.addTerm(mean, sigma1, n1);
    u2_data.addTerm(mean, sigma2, n2);
    u2_data.addTerm(mean, sigma3, n3);
    
    */
    
    ///////////////////

    
    GaussianSum u1_mc;
    h3 = recoil_hists["mc_para"];
    corrIdx = 1;
    if(corrType == 3) corrIdx = corrParam;
    norm_para_mc =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u1_mc.setNorm(norm_para_mc);
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    

    
    GaussianSum u2_mc;
    h3 = recoil_hists["mc_perp"];
    corrIdx = 1;
    if(corrType == 4) corrIdx = corrParam;
    norm_perp_mc =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u2_mc.setNorm(norm_perp_mc);
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
    /*
    GaussianSum u2_mc;
    double mean_u2_mc = 0;
    double sigma1_u2_mc = recoilParameterization_sigma(3.01405e+00, 4.94889e+00, 1.54923e-01, qT);
    double sigma2_u2_mc = recoilParameterization_sigma(5.49714e+00, 8.45153e+00, 1.21679e-01, qT);
    double sigma3_u2_mc = recoilParameterization_sigma(8.08396e+00, 7.06204e+00, 9.56741e-02, qT);
    double sigma4_u2_mc = recoilParameterization_sigma(9.64523e-02, 1.20944e+02, 1.03714e+00, qT);
    double n1_u2_mc = 0.22;
    double n2_u2_mc = 0.5;
    double n3_u2_mc = 0.25;
    double n4_u2_mc = 1. - n1_u2_mc - n2_u2_mc - n3_u2_mc;
    u2_mc.addTerm(mean_u2_mc, sigma1_u2_mc, n1_u2_mc);
    u2_mc.addTerm(mean_u2_mc, sigma2_u2_mc, n2_u2_mc);
    u2_mc.addTerm(mean_u2_mc, sigma3_u2_mc, n3_u2_mc);
    u2_mc.addTerm(mean_u2_mc, sigma4_u2_mc, n4_u2_mc);
    */



    double pVal_u1_mc = u1_mc.evalCDF(pU1 + qT);
    double pVal_u2_mc = u2_mc.evalCDF(pU2);
    /*
    double pU1_mc = u1_mc.findRoot(pVal_u1_mc, -1000, 1000);
    double pU2_mc = u2_mc.findRoot(pVal_u2_mc, -500, 500);
    
    if(abs(pU2_mc-pU2) > 0.1) {
        
        cout << pU2_mc << " " << pU2 << endl;
    }*/
    
    /*
    if(qTbinIdx == 30) {
             cout << norm_para_data << " " << norm_perp_data << " " << norm_para_mc << " " << norm_perp_mc << " " << pVal_u1_mc << " " << pVal_u2_mc << endl;
    }
    else return res;
    */
    
    double pU1_data = u1_data.findRoot(pVal_u1_mc, -1000, 1000);
    double pU2_data = u2_data.findRoot(pVal_u2_mc, -500, 500);
    //double pU2_data=  999999;
    //double pU1_data=  999999;

  
    if(pU1_data == 999999) cout << "************* " << qT << endl;
    if(pU2_data == 999999) cout << "************* " << qT << endl;
    
    
    
    //if(qTbinIdx == 1) cout << norm_para_data << " " << norm_perp_data << " " << norm_para_mc << " " << norm_perp_mc << endl;

    if(pU1_data == 999999) pU1_data = pU1;
    else pU1_data -= qT;
    
    
    //else pU1_data -= qT;
    //else {
    //   if((pU1+qT)*(pU1_data+qT) > 0) pU1_data -= qT; // SS
    //   else pU1_data += qT; // OS
    //}
    //else pU1_data -= qT;
    
    if((pU1+qT)*(pU1_data+qT) < 0) {
        
        //cout << "OS " << qT << " " << qTbinIdx << " " << pU1+qT << " " << pU1_data+qT << " " << pU2 << " " << pU2_data << endl;
    }
    
    if(pU2_data == 999999) pU1_data = pU2;    
        
    // correct MC to DATA
    res[1] = pU1_data;
    res[2] = pU2_data;
    res[0] = std::hypot(res[1], res[2]);

	return res;
}


Vec_f recoilCorrectionBinnedWtoZ(int q, double pU1, double pU2, double qTbinIdx, double qT, int corrType=0, int corrParam=1) {
    
    Vec_f res(3, 0);
    //cout << pU1 << " " << pU2 << " " << qTbinIdx << endl; 
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    return res;
    
    TH3D *h3;
    int corrIdx = 1;
    
    
    //cout << "************* " << qT << endl;
    
    GaussianSum u1_Z;
    h3 = recoil_hists_gen_Z["gen_para"];
    corrIdx = 1;
    if(corrType == 1) corrIdx = corrParam;
    u1_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
    
    GaussianSum u2_Z;
    h3 = recoil_hists_gen_Z["gen_perp"];
    corrIdx = 1;
    if(corrType == 2) corrIdx = corrParam;
    u2_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_Z.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
    
    
    GaussianSum u1_W;
    if(q > 0) h3 = recoil_hists_gen_W["gen_plus_para"];
    else h3 = recoil_hists_gen_W["gen_minus_para"];
    corrIdx = 1;
    if(corrType == 3) corrIdx = corrParam;
    u1_W.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_W.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_W.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    

    
    GaussianSum u2_W;
    if(q > 0) h3 = recoil_hists_gen_W["gen_plus_perp"];
    else h3 = recoil_hists_gen_W["gen_minus_perp"];
    corrIdx = 1;
    if(corrType == 4) corrIdx = corrParam;
    u2_W.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_W.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_W.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));


    
    double pVal_u1_W = u1_W.evalCDF(pU1 + qT);
    double pVal_u2_W = u2_W.evalCDF(pU2);

    
    double pU1_Z = u1_Z.findRoot(pVal_u1_W, -1000, 1000);
    double pU2_Z = u2_Z.findRoot(pVal_u2_W, -500, 500);
  
    if(pU1_Z == 999999) cout << "************* " << qT << endl;
    if(pU2_Z == 999999) cout << "************* " << qT << endl;
    
    if(pU1_Z == 999999) pU1_Z = pU1;
    else pU1_Z -= qT;
    

    if(pU2_Z == 999999) pU2_Z = pU2;    
        
    // correct W to Z
    res[1] = pU1_Z;
    res[2] = pU2_Z;
    res[0] = std::hypot(res[1], res[2]);

	return res;
}




// recoil correction
Vec_f recoilCorrectionBinnedCut(double pU1, double pU2, double qTbinIdx, double qT, Vec_f lep_pt, Vec_f lep_eta, Vec_f lep_phi, double MET_phi, double lep_phi_, double MET, int corrType=0, int corrParam=1) {
    
    Vec_f res(3, 0);
    //cout << pU1 << " " << pU2 << " " << qTbinIdx << endl; 
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    
    TH3D *h3;
    int corrIdx = 1;
    
    double norm_para_data = 1;
    double norm_perp_data = 1;
    double norm_para_mc = 1;
    double norm_perp_mc = 1;
    
    //cout << "************* " << qT << endl;
    
    GaussianSum u1_data;
    h3 = recoil_hists["data_para"];
    corrIdx = 1;
    if(corrType == 1) corrIdx = corrParam;
    norm_para_data =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u1_data.setNorm(norm_para_data);
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    
/*
    if(corrType == 1) {
        if(qTbinIdx == 40) {
            cout << corrIdx << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 4, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 5, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 6, corrIdx) << endl;
            cout << h3->GetBinContent(qTbinIdx+1, 7, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 8, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 9, corrIdx) << endl;
        }
    }*/
    /*
    if(qT > 280) {
        cout << "************* " << qT << " " << qTbinIdx << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 4, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 5, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 6, corrIdx) << endl;
    cout << h3->GetBinContent(qTbinIdx+1, 7, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 8, corrIdx) << "  " << h3->GetBinContent(qTbinIdx+1, 9, corrIdx) << endl;
    }
    */
    //if(qTbinIdx == 3) {
    //    cout << h3->GetBinContent(qTbinIdx+1, 1, corrIdx) << " " << h3->GetBinContent(qTbinIdx+1, 2, corrIdx) << " " << h3->GetBinContent(qTbinIdx+1, 3, corrIdx) << endl;
    //}
    
    GaussianSum u2_data;
    h3 = recoil_hists["data_perp"];
    corrIdx = 1;
    if(corrType == 2) corrIdx = corrParam;
    norm_perp_data =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u2_data.setNorm(norm_perp_data);
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    

    
    GaussianSum u1_mc;
    h3 = recoil_hists["mc_para"];
    corrIdx = 1;
    if(corrType == 3) corrIdx = corrParam;
    norm_para_mc =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u1_mc.setNorm(norm_para_mc);
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    

    
    GaussianSum u2_mc;
    h3 = recoil_hists["mc_perp"];
    corrIdx = 1;
    if(corrType == 4) corrIdx = corrParam;
    norm_perp_mc =  h3->GetBinContent(qTbinIdx+1, 10, 1);
    //u2_mc.setNorm(norm_perp_mc);
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 4, corrIdx), h3->GetBinContent(qTbinIdx+1, 5, corrIdx), h3->GetBinContent(qTbinIdx+1, 6, corrIdx));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 7, corrIdx), h3->GetBinContent(qTbinIdx+1, 8, corrIdx), h3->GetBinContent(qTbinIdx+1, 9, corrIdx));
    

    
    
    
    double pVal_u1_mc = u1_mc.evalCDF(pU1 + qT);
    double pVal_u2_mc = u2_mc.evalCDF(pU2);
    /*
    if(qTbinIdx == 30) {
             cout << norm_para_data << " " << norm_perp_data << " " << norm_para_mc << " " << norm_perp_mc << " " << pVal_u1_mc << " " << pVal_u2_mc << endl;
    }
    else return res;
    */
    
    double pU1_data = u1_data.findRoot(pVal_u1_mc, -1000, 1000);
    double pU2_data = u2_data.findRoot(pVal_u2_mc, -500, 500);
    //double pU2_data=  999999;
    //double pU1_data=  999999;
    
    cout << "***********************************************************" << endl;
    cout << qT << " " << qTbinIdx << endl;
    
    cout << hypot(pU1 + qT, pU2) << " " << hypot(pU1_data, pU2_data) << endl;
    cout << "pU1_data=" << pU1_data << " pU2_data=" << pU2_data << endl; 
  
    if(pU1_data == 999999) cout << "************* " << qT << endl;
    if(pU2_data == 999999) cout << "************* " << qT << endl;
    
    
    
    //if(qTbinIdx == 1) cout << norm_para_data << " " << norm_perp_data << " " << norm_para_mc << " " << norm_perp_mc << endl;

    if(pU1_data == 999999) pU1_data = pU1;
    //else {
    //   if(pU1_data > 0) pU1_data -= qT;
    //   else pU1_data += qT;
    //}
    else pU1_data -= qT;
    
    if(pU2_data == 999999) pU1_data = pU2;
    
    
    cout << "PARA+qT " << (pU1+qT) << " " << (pU1_data+qT) << endl; 
    cout << "PARA " << pU1 << " " << pU1_data << endl; 
    cout << "PERP " << pU2 << " " << pU2_data << endl; 
    
    cout << "lep_pt1=" << lep_pt[0] << " lep_pt2=" << lep_pt[1] << endl;  
    cout << "lep_eta1=" << lep_eta[0] << " lep_eta2=" << lep_eta[1] << endl;  
    cout << "lep_phi1=" << lep_phi[0] << " lep_phi2=" << lep_phi[1] << endl;  
    cout << "MET=" << MET <<" MET_phi=" << MET_phi << " lep_phi=" << lep_phi_ << endl;  
    cout << "***********************************************************" << endl;
    
    if(qT < 36.3028 and qT > 36.3024) {
        
        exit(0);
    }
    
    // correct MC to DATA
    res[1] = pU1_data;
    res[2] = pU2_data;
    res[0] = std::hypot(res[1], res[2]);

	return res;
}



// recoil correction
Vec_f recoilCorrectionBinned_old(double pU1, double pU2, double qTbinIdx, int corrType=0, int corrParam=1) {
    
    Vec_f res(3, 0);
    //cout << pU1 << " " << pU2 << " " << qTbinIdx << endl; 
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    
    TH3D *h3;
    int corrIdx = 1;
    
    GaussianSum u1_data;
    h3 = recoil_hists["para_data"];
    corrIdx = 1;
    if(corrType == 1) corrIdx = corrParam;
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx), 1.-h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    
    GaussianSum u2_data;
    h3 = recoil_hists["perp_data"];
    corrIdx = 1;
    if(corrType == 2) corrIdx = corrParam;
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx), 1.-h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    
    GaussianSum u1_mc;
    h3 = recoil_hists["para_mc"];
    corrIdx = 1;
    if(corrType == 3) corrIdx = corrParam;
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx), 1.-h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    
    GaussianSum u2_mc;
    h3 = recoil_hists["perp_mc"];
    corrIdx = 1;
    if(corrType == 4) corrIdx = corrParam;
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 2, corrIdx), h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, corrIdx), h3->GetBinContent(qTbinIdx+1, 3, corrIdx), 1.-h3->GetBinContent(qTbinIdx+1, 4, corrIdx));
    
    double pVal_u1_mc = u1_mc.evalCDF(pU1);
    double pVal_u2_mc = u2_mc.evalCDF(pU2);
    
    double pU1_data = u1_data.findRoot(pVal_u1_mc, -500, 500);
    double pU2_data = u2_data.findRoot(pVal_u2_mc, -500, 500);

    if(pU1_data == 999999) pU1_data = pU1;
    if(pU2_data == 999999) pU1_data = pU2;    
        
    // correct MC to DATA
    res[1] = pU1_data;
    res[2] = pU2_data;
    res[0] = std::hypot(res[1], res[2]);

	return res;
}

// recoil correction



double recoilParameterization_mean(double a, double b, double c, double d, double e, double f, double g, double qT) {
    
    return (a*qT*qT*qT*qT + b*qT*qT*qT + c*qT*qT + d*qT)/(e*qT*qT*qT + f*qT*qT + g*qT + 1);
}


GaussianSum constructParametricGauss(string tag, double qT, int syst=-1) {
    //cout << "tag=" << tag << "qT=" << qT << endl;
    GaussianSum gauss;
    double totNorm = 0.;
    double mean, sigma, norm;
    string systLabel = "";
    
    //if(tag == "target_perp" and syst==1) syst = -1;
    if(syst != -1) systLabel = "_syst" + to_string(syst);
     
    if(recoil_param_nGauss.find(tag) != recoil_param_nGauss.end()) { // parametric
        //cout << " PARAMETIC " << endl; 
        int nGauss = recoil_param_nGauss[tag];
        for(int i=1; i<=nGauss; i++) {
            if(i == nGauss) norm = 1.-totNorm;
            else {
                norm = recoil_param_funcs[tag + "_norm" + to_string(i) + systLabel]->Eval(qT);
                totNorm += norm;
            }
            mean = recoil_param_funcs[tag + "_mean" + to_string(i) + systLabel]->Eval(qT);
            sigma = recoil_param_funcs[tag + "_sigma" + to_string(i) + systLabel]->Eval(qT);
            gauss.addTerm(mean, sigma, norm);
            //cout << "  iGauss=" << i << " norm=" << norm << " mean=" << mean << " sigma=" << sigma << endl;
        }
    }
    else if (recoil_binned_nGauss.find(tag) != recoil_binned_nGauss.end()) {
        //cout << " BINNED " << endl; 
        // find the bin
        int qTbin = -1;
        int nGauss = recoil_binned_nGauss[tag];
        for(unsigned int i=0; i < recoil_binned_qTbins[tag].size()-1; i++) {
            if(qT >= recoil_binned_qTbins[tag].at(i) and qT < recoil_binned_qTbins[tag].at(i+1)) {
                qTbin = i;
                break;
            }
        }
        if(qTbin == -1) {
            qTbin = recoil_binned_qTbins[tag].size()-2;
            cout << "Bin outside range for tag=" << tag << " and qT=" << qT << endl; 
        }
        
        // construnct the Gauss
        for(int i=0; i<nGauss; i++) {
            mean = recoil_binned_components[tag].at(qTbin).at(i).at(0);
            sigma = recoil_binned_components[tag].at(qTbin).at(i).at(1);
            norm = recoil_binned_components[tag].at(qTbin).at(i).at(2);
            gauss.addTerm(mean, sigma, norm);
            //cout << " ---> qT=" << qT << " qTbin=" << qTbin << " nGauss=" << i << " mean=" << mean << " sigma=" << sigma << " norm=" << norm << endl;
        }  
    }
    else {
        exit(1);
    }
    return gauss;    
}


// recoil correction
Vec_f recoilCorrectionParametric(double para, double perp, double qTbinIdx, double qT, string systTag="", int systIdx=-1) {
    
    Vec_f res(3, 0);
    //cout << pU1 << " " << pU2 << " " << qTbinIdx << endl; 
    
    res[1] = para;
    res[2] = perp;
    res[0] = std::hypot(res[1], res[2]);
    //return res;
    
    int corrIdx = 1;
    
    double norm_para_data = 1;
    double norm_perp_data = 1;
    double norm_para_mc = 1;
    double norm_perp_mc = 1;
    
    GaussianSum data_para;
    GaussianSum data_perp;
    GaussianSum dy_para;
    GaussianSum dy_perp;
    
    if(systIdx != -1 and systTag == "target_para") data_para = constructParametricGauss("target_para", qT, systIdx);
    else data_para = constructParametricGauss("target_para", qT);
    if(systIdx != -1 and systTag == "target_perp") data_perp = constructParametricGauss("target_perp", qT, systIdx);
    else data_perp = constructParametricGauss("target_perp", qT);
    if(systIdx != -1 and systTag == "source_para") dy_para = constructParametricGauss("source_para", qT, systIdx);
    else dy_para = constructParametricGauss("source_para", qT);
    if(systIdx != -1 and systTag == "source_perp") dy_perp = constructParametricGauss("source_perp", qT, systIdx);
    else dy_perp = constructParametricGauss("source_perp", qT);
    
    //GaussianSum data_perp = constructParametricGauss("target_perp", qT);
    //GaussianSum dy_para = constructParametricGauss("source_para", qT);
    //GaussianSum dy_perp = constructParametricGauss("source_perp", qT);
    
    
    double pVal_para_dy = dy_para.evalCDF(para + qT);
    double pVal_perp_dy = dy_perp.evalCDF(perp);
    
    
    double para_corr = data_para.findRoot(pVal_para_dy, -1000, 1000);
    double perp_corr = data_perp.findRoot(pVal_perp_dy, -500, 500);

  
    if(para_corr == 999999) cout << "************* " << qT << endl;
    if(perp_corr == 999999) cout << "************* " << qT << endl;
    
    
  
    if(para_corr == 999999) para_corr = para;
    else para_corr -= qT;
    
    
    //else pU1_data -= qT;
    //else {
    //   if((pU1+qT)*(pU1_data+qT) > 0) pU1_data -= qT; // SS
    //   else pU1_data += qT; // OS
    //}
    //else pU1_data -= qT;
    
    //if((pU1+qT)*(pU1_data+qT) < 0) {
        
        //cout << "OS " << qT << " " << qTbinIdx << " " << pU1+qT << " " << pU1_data+qT << " " << pU2 << " " << pU2_data << endl;
    //}
    
    if(perp_corr == 999999) perp_corr = perp;    
        
    // correct MC to DATA
    res[1] = para_corr;
    res[2] = perp_corr;
    res[0] = std::hypot(res[1], res[2]);

	return res;    
}


// recoil correction
Vec_recoil recoilCorrectionParametricUnc(double para, double perp, double qT, string tag) {
    
    int nSysts = recoil_param_nStat[tag];
    Vec_recoil res(nSysts);
    for(int iSyst = 0; iSyst<nSysts; iSyst++) {
        Vec_f ret = recoilCorrectionParametric(para, perp, 0, qT, tag, iSyst);
        res[iSyst].para = ret[1];
        res[iSyst].perp = ret[2];
    }
    return res;
}

Vec_f recoilCorrectionParametric_Old(double pU1, double pU2, double qT, int corrType=0, int corrIdx=1) {
    
    /*
    corrType:
        - 0: none
        - 1: para data
        - 1: perp data
        - 1: para mc
        - 1: perp mc
        
    corrIdx: according to the variation/diagonalization of covariance matrix
        
    */
    
    Vec_f res(3, 0);
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    
    TH1D *h1;
    double m, s1, s2, n;
    int sIdx = 1;
    
    GaussianSum u1_data;
    h1 = recoil_hists_param["para_data"];
    if(corrType == 1) sIdx = corrIdx;
    m = recoilParameterization_mean(h1->GetBinContent(1, 1, sIdx), h1->GetBinContent(1, 2, sIdx), h1->GetBinContent(1, 3, sIdx), h1->GetBinContent(1, 4, sIdx), h1->GetBinContent(1, 5, sIdx), h1->GetBinContent(1, 6, sIdx), h1->GetBinContent(1, 7, sIdx), qT);
    s1 = recoilParameterization_sigma(h1->GetBinContent(2, 1, sIdx), h1->GetBinContent(2, 2, sIdx), h1->GetBinContent(2, 3, sIdx), qT);
    s2 = recoilParameterization_sigma(h1->GetBinContent(3, 1, sIdx), h1->GetBinContent(3, 2, sIdx), h1->GetBinContent(3, 3, sIdx), qT);
    n = h1->GetBinContent(4, 1, sIdx);
    u1_data.addTerm(m, s1, n);
    u1_data.addTerm(m, s2, 1.-n);
    
    GaussianSum u2_data;
    h1 = recoil_hists_param["perp_data"];
    if(corrType == 2) sIdx = corrIdx;
    s1 = recoilParameterization_sigma(h1->GetBinContent(2, 1, sIdx), h1->GetBinContent(2, 2, sIdx), h1->GetBinContent(2, 3, sIdx), qT);
    s2 = recoilParameterization_sigma(h1->GetBinContent(3, 1, sIdx), h1->GetBinContent(3, 2, sIdx), h1->GetBinContent(3, 3, sIdx), qT);
    n = h1->GetBinContent(4, 1, sIdx);
    u2_data.addTerm(0, s1, n);
    u2_data.addTerm(0, s2, 1.-n);
    
    GaussianSum u1_mc;
    h1 = recoil_hists_param["para_mc"];
    if(corrType == 3) sIdx = corrIdx;
    m = recoilParameterization_mean(h1->GetBinContent(1, 1, sIdx), h1->GetBinContent(1, 2, sIdx), h1->GetBinContent(1, 3, sIdx), h1->GetBinContent(1, 4, sIdx), h1->GetBinContent(1, 5, sIdx), h1->GetBinContent(1, 6, sIdx), h1->GetBinContent(1, 7, sIdx), qT);
    s1 = recoilParameterization_sigma(h1->GetBinContent(2, 1, sIdx), h1->GetBinContent(2, 2, sIdx), h1->GetBinContent(2, 3, sIdx), qT);
    s2 = recoilParameterization_sigma(h1->GetBinContent(3, 1, sIdx), h1->GetBinContent(3, 2, sIdx), h1->GetBinContent(3, 3, sIdx), qT);
    n = h1->GetBinContent(4, 1, sIdx);
    u1_mc.addTerm(m, s1, n);
    u1_mc.addTerm(m, s2, 1.-n);
    
    GaussianSum u2_mc;
    h1 = recoil_hists_param["perp_mc"];
    if(corrType == 4) sIdx = corrIdx;
    s1 = recoilParameterization_sigma(h1->GetBinContent(2, 1, sIdx), h1->GetBinContent(2, 2, sIdx), h1->GetBinContent(2, 3, sIdx), qT);
    s2 = recoilParameterization_sigma(h1->GetBinContent(3, 1, sIdx), h1->GetBinContent(3, 2, sIdx), h1->GetBinContent(3, 3, sIdx), qT);
    n = h1->GetBinContent(4, 1, sIdx);
    u2_mc.addTerm(0, s1, n);
    u2_mc.addTerm(0, s2, 1.-n);
    
    double pVal_u1_mc = u1_mc.evalCDF(pU1);
    double pVal_u2_mc = u2_mc.evalCDF(pU2);
    
    double pU1_data = u1_data.findRoot(pVal_u1_mc, -1000, 1000);
    double pU2_data = u2_data.findRoot(pVal_u2_mc, -500, 500);

    if(pU1_data == 999999) pU1_data = pU1;
    if(pU2_data == 999999) pU1_data = pU2;    
        
    // correct MC to DATA
    res[1] = pU1_data;
    res[2] = pU2_data;
    res[0] = std::hypot(res[1], res[2]);


	return res;	
}




Vec_f METCorrection(double MET_pt, double MET_phi, double rec_para, double rec_perp, double V_pt, double V_phi) {

    Vec_f res(2, 0);

        
    double lMX = -V_pt*cos(V_phi) - rec_para*cos(V_phi) + rec_perp*sin(V_phi);
    double lMY = -V_pt*sin(V_phi) - rec_para*sin(V_phi) - rec_perp*cos(V_phi);

        
    //res[0] = sqrt(lMX*lMX + lMY*lMY);
    res[0] = sqrt((rec_para+V_pt)*(rec_para+V_pt) + rec_perp*rec_perp);
    res[1] = atan2(lMY, lMX);
    //if(lMX > 0) res[1] = atan(lMY/lMX);
    //else res[1] = (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX);
  
	return res;
}

Vec_f METCorrectionGen(double rec_para, double rec_perp, double lep_pt, double lep_phi, double V_phi) {

    Vec_f res(2, 0);

        
    double lMX = -lep_pt*cos(lep_phi) - rec_para*cos(V_phi) + rec_perp*sin(V_phi);
    double lMY = -lep_pt*sin(lep_phi) - rec_para*sin(V_phi) - rec_perp*cos(V_phi);

        
    res[0] = sqrt(lMX*lMX + lMY*lMY);
    //res[0] = sqrt((rec_para+lep_pt)*(rec_para+lep_pt) + rec_perp*rec_perp);
    res[1] = atan2(lMY, lMX);
    //if(lMX > 0) res[1] = atan(lMY/lMX);
    //else res[1] = (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX);
  
	return res;
}

/*
Vec_f METCorrection_StatUnc(double MET_pt, double MET_phi, double rec_para, double rec_perp, double V_pt, double V_phi) {

    Vec_f res(2, 0);

        
    double lMX = -V_pt*cos(V_phi) - rec_para*cos(V_phi) + rec_perp*sin(V_phi);
    double lMY = -V_pt*sin(V_phi) - rec_para*sin(V_phi) - rec_perp*cos(V_phi);

        
    res[0] = sqrt(lMX*lMX + lMY*lMY);
    if(lMX > 0) res[1] = atan(lMY/lMX);
    else res[1] = (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX);
  
	return res;
}*/




Vec_f recoilCorrectionParametric_MET_pt_unc(Vec_recoil recoil, double V_pt, double V_phi) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) {
       
        double lMX = -V_pt*cos(V_phi) - recoil[iSyst].para*cos(V_phi) + recoil[iSyst].perp*sin(V_phi);
        double lMY = -V_pt*sin(V_phi) - recoil[iSyst].para*sin(V_phi) - recoil[iSyst].perp*cos(V_phi);     
        res[iSyst] = sqrt(lMX*lMX + lMY*lMY);
    }
	return res;
}

Vec_f recoilCorrectionParametric_MET_phi_unc(Vec_recoil recoil, double V_pt, double V_phi) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) {
       
        double lMX = -V_pt*cos(V_phi) - recoil[iSyst].para*cos(V_phi) + recoil[iSyst].perp*sin(V_phi);
        double lMY = -V_pt*sin(V_phi) - recoil[iSyst].para*sin(V_phi) - recoil[iSyst].perp*cos(V_phi);     
        res[iSyst] = atan2(lMY, lMX);
    }
	return res;
}

Vec_f recoilCorrectionParametric_para_qT_unc(Vec_recoil recoil, double qT) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) res[iSyst] = recoil[iSyst].para + qT;
	return res;
}

Vec_f recoilCorrectionParametric_para_unc(Vec_recoil recoil, double qT) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) res[iSyst] = recoil[iSyst].para;
	return res;
}

Vec_f recoilCorrectionParametric_perp_unc(Vec_recoil recoil, double qT) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) res[iSyst] = recoil[iSyst].perp;
	return res;
}

Vec_f recoilCorrectionParametric_magn_unc(Vec_recoil recoil, double qT) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) res[iSyst] = sqrt(recoil[iSyst].para*recoil[iSyst].para + recoil[iSyst].perp*recoil[iSyst].perp);
	return res;
}

Vec_f recoilCorrectionParametric_mT_unc(Vec_f met_pt, Vec_f met_phi, float pt, float phi, float ptOther, float phiOther) {
    
    int nBins = met_pt.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) {
        
        TVector2 pl = TVector2();
        pl.SetMagPhi(ptOther, phiOther);

        TVector2 met_wlike = TVector2();
        met_wlike.SetMagPhi(met_pt[iSyst], met_phi[iSyst]);
        met_wlike = pl + met_wlike;

        res[iSyst] = std::sqrt(2*pt*met_wlike.Mod()*(1-std::cos(phi-met_wlike.Phi())));  
    }
	return res;
}

Vec_f recoilCorrectionParametric_MET_pt_gen_unc(Vec_recoil recoil, double lep_pt, double lep_phi, double V_phi) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) {
       
        double lMX = -lep_pt*cos(lep_phi) - recoil[iSyst].para*cos(V_phi) + recoil[iSyst].perp*sin(V_phi);
        double lMY = -lep_pt*sin(lep_phi) - recoil[iSyst].para*sin(V_phi) - recoil[iSyst].perp*cos(V_phi);   
        res[iSyst] = sqrt(lMX*lMX + lMY*lMY);
    }
	return res;
}

Vec_f recoilCorrectionParametric_MET_phi_gen_unc(Vec_recoil recoil, double lep_pt, double lep_phi, double V_phi) {
    
    int nBins = recoil.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) {
       
        double lMX = -lep_pt*cos(lep_phi) - recoil[iSyst].para*cos(V_phi) + recoil[iSyst].perp*sin(V_phi);
        double lMY = -lep_pt*sin(lep_phi) - recoil[iSyst].para*sin(V_phi) - recoil[iSyst].perp*cos(V_phi);   
        res[iSyst] = atan2(lMY, lMX);
    }
	return res;
}

Vec_f recoilCorrectionParametric_mT_2_unc(Vec_f met_pt, Vec_f met_phi, float pt, float phi) {
    
    int nBins = met_pt.size();
    Vec_f res(nBins, -1e5);
    for(int iSyst=0; iSyst<nBins; iSyst++) res[iSyst] = std::sqrt(2*pt*met_pt[iSyst]*(1-std::cos(phi-met_phi[iSyst])));
	return res;
}


/*
Binned recoil correction functions
*/
Vec_f recoilCorrectionBinned_magn_StatUnc(Vec_f pU, int qTbinIdx) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    res[qTbinIdx] = pU[0];
	return res;
}

Vec_f recoilCorrectionBinned_para_StatUnc(Vec_f pU, int qTbinIdx) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    res[qTbinIdx] = pU[1];
	return res;
}

Vec_f recoilCorrectionBinned_para_qT_StatUnc(Vec_f pU, int qTbinIdx, double qT) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5); // BE CAREFUL, SHOULD BE UNDERFOW
    res[qTbinIdx] = pU[1] + qT;
	return res;
}

Vec_f recoilCorrectionBinned_perp_StatUnc(Vec_f pU, int qTbinIdx) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5); // // BE CAREFUL, SHOULD BE UNDERFOW
    res[qTbinIdx] = pU[2];
	return res;
}


Vec_f recoilCorrectionBinned_MET_pt_StatUnc(Vec_f pU, int qTbinIdx, double MET_pt, double MET_phi, double V_pt, double V_phi) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    
            
    double lMX = -V_pt*cos(V_phi) - pU[1]*cos(V_phi) + pU[2]*sin(V_phi);
    double lMY = -V_pt*sin(V_phi) - pU[1]*sin(V_phi) - pU[2]*cos(V_phi);

    res[qTbinIdx] = sqrt(lMX*lMX + lMY*lMY);
	return res;
}

Vec_f recoilCorrectionBinned_MET_phi_StatUnc(Vec_f pU, int qTbinIdx, double MET_pt, double MET_phi, double V_pt, double V_phi) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    
            
    double lMX = -V_pt*cos(V_phi) - pU[1]*cos(V_phi) + pU[2]*sin(V_phi);
    double lMY = -V_pt*sin(V_phi) - pU[1]*sin(V_phi) - pU[2]*cos(V_phi);

    res[qTbinIdx] = atan2(lMY, lMX);
	return res;
}






Vec_f recoilCorrectionBinned_MET_pt_gen_StatUnc(Vec_f pU, int qTbinIdx, double lep_pt, double lep_phi, double V_phi) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    
            
    double lMX = -lep_pt*cos(lep_phi) - pU[1]*cos(V_phi) + pU[2]*sin(V_phi);
    double lMY = -lep_pt*sin(lep_phi) - pU[1]*sin(V_phi) - pU[2]*cos(V_phi);

    res[qTbinIdx] = sqrt(lMX*lMX + lMY*lMY);
	return res;
}

Vec_f recoilCorrectionBinned_MET_phi_gen_StatUnc(Vec_f pU, int qTbinIdx, double lep_pt, double lep_phi, double V_phi) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    
            
    double lMX = -lep_pt*cos(lep_phi) - pU[1]*cos(V_phi) + pU[2]*sin(V_phi);
    double lMY = -lep_pt*sin(lep_phi) - pU[1]*sin(V_phi) - pU[2]*cos(V_phi);

    res[qTbinIdx] = atan2(lMY, lMX);
	return res;
}




Vec_f recoilCorrectionBinned_mt_StatUnc(Vec_f met_pt, Vec_f met_phi, int qTbinIdx, float pt, float phi, float ptOther, float phiOther) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
    
    TVector2 pl = TVector2();
    pl.SetMagPhi(ptOther,phiOther);

    TVector2 met_wlike = TVector2();
    met_wlike.SetMagPhi(met_pt[qTbinIdx], met_phi[qTbinIdx]);
    met_wlike = pl + met_wlike;

    res[qTbinIdx] = std::sqrt(2*pt*met_wlike.Mod()*(1-std::cos(phi-met_wlike.Phi())));  
	return res;
}


Vec_f recoilCorrectionBinned_mt_2_StatUnc(Vec_f met_pt, Vec_f met_phi, int qTbinIdx, float pt, float phi) {
    
    int nqTbins = qTbins.size();
    Vec_f res(nqTbins, -1e5);
  
    res[qTbinIdx] = std::sqrt(2*pt*met_pt[qTbinIdx]*(1-std::cos(phi-met_phi[qTbinIdx])));;
	return res;
}






// recoil correction
/*
Vec_f recoilCorrectionBinned_StatUnc(double pU1, double pU2, double qTbinIdx, int corrType, int corrIdx) {
    
    int nqTbins = qTbins.size();

    Vec_f t(3, -1);
    Vec_f res(nqTbins, -1); // per qT bin
    
    t = recoilCorrectionBinned(pU1, pU2, qTbinIdx, corrType, corrIdx);
    res[qTbinIdx] = t[0];
    
	return res;
}


Vec_recoilType recoilCorrectionBinned_StatUnc(double pU1, double pU2, double qTbinIdx, int corrType, int corrIdx) {
    
    int nqTbins = qTbins.size();

    Vec_f t(3, -1);
    
    RecoilType tmp;
    tmp.pu1 = 0;
    tmp.pu2 = 0;
    
    Vec_recoilType res(nqTbins, tmp);
    
    t = recoilCorrectionBinned(pU1, pU2, qTbinIdx, corrType, corrIdx);
    res[qTbinIdx].pu1 = t[1];
    res[qTbinIdx].pu1 = t[2];
    
	return res;
}

Vec_f recoilCorrectionBinned_magn_StatUnc(Vec_recoilType pU) {
    
    unsigned size = pU.size();
    Vec_f res(size, -1);
    
    for(unsigned int i=0; i < size; i++) {
        res[i] = sqrt(pU[i].pu1*pU[i].pu1 + pU[i].pu2*pU[i].pu2);
    }

	return res;
}

*/
Vec_f recoilCorrectionBinned_StatUnc(double pU1, double pU2, int qTbinIdx, int corrType, int corrIdx) {
    
    Vec_f res = recoilCorrectionBinned(pU1, pU2, qTbinIdx, corrType, corrIdx);
	return res;
}






}

#endif
