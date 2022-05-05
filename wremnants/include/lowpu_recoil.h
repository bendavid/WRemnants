#ifndef WREMNANTS_LOWPU_RECOIL_H
#define WREMNANTS_LOWPU_RECOIL_H


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include "defines.h"

namespace wrem {
    
TFile *recoil_fits_Z_param = new TFile("wremnants/data/lowPU/recoil_fits_Z_param_refit.root", "READ");
//TFile *recoil_fits_Z_param_refit = new TFile("wremnants/data/lowPU/recoil_fits_Z_param_refit.root", "READ");
TFile *recoil_fits_Z = new TFile("wremnants/data/lowPU/recoil_fits_Z.root", "READ");
std::map<std::string, TH3D*> recoil_hists;
std::map<std::string, TH1D*> recoil_hists_param;

std::vector<float> qTbins;

void recoil_init() {
    
    TH3D *h3;
    TH1D *h1;
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
    

 
    h3 = (TH3D*)recoil_fits_Z->Get("para_mc");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("perp_mc");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("para_data");
    recoil_hists.insert({h3->GetName(), h3});
    
    h3 = (TH3D*)recoil_fits_Z->Get("perp_data");
    recoil_hists.insert({h3->GetName(), h3});

    h1 = (TH1D*)recoil_fits_Z_param->Get("para_mc");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("perp_mc");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("para_data");
    recoil_hists_param.insert({h1->GetName(), h1});
    
    h1 = (TH1D*)recoil_fits_Z_param->Get("perp_data");
    recoil_hists_param.insert({h1->GetName(), h1});

}

int getqTbin(float qT) {

    for(unsigned int i=0; i < qTbins.size()-1; i++) {
        
        if(qT >= qTbins.at(i) and qT < qTbins.at(i+1)) return i;
    }
    
    cout << "[getQtBin] Bin not found qT=" << qT << endl;
    return -1;
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
    double MET_phi_corr = MET_corr.Phi();
	
    Vec_f res(2, 0);
    res[0] = MET_pt_corr;
    res[1] = MET_phi_corr;
	
	return res;
	
}

Vec_f recoilComponents(double MET_pt, double MET_phi, double lep_pt, double lep_phi) {
    
    // computes the hadronic recoil, defined as -(MET + V), V being the vector sum of the (di)lepton system in the lab frame

    
    double pUx  = MET_pt*cos(MET_phi) + lep_pt*cos(lep_phi); // recoil x in lab frame
	double pUy  = MET_pt*sin(MET_phi) + lep_pt*sin(lep_phi); // recoil y in lab frame
	double pU   = std::hypot(pUx, pUy);
	double pCos = - (pUx*cos(lep_phi) + pUy*sin(lep_phi))/pU;
	double pSin =   (pUx*sin(lep_phi) - pUy*cos(lep_phi))/pU;
	double Ux	= pU*pCos; // recoil x along boson pT
	double Uy	= pU*pSin; // recoil y along boson pT
    
    Vec_f res(3, 0);
    res[0] = pU;
    res[1] = Ux;
    res[2] = Uy;
    
	return res;
}



class GaussianSum {
    
    private:
    
        std::vector<double> mean;
        std::vector<double> sigma;
        std::vector<double> norm;
        double totNorm;
        
        int MaxIterations = 512;
        double _tol = 2.2204460492503131e-16;
        
        bool verbose;
    
    
    public:
    
        GaussianSum() {};
        
        void addTerm(double mean_, double sigma_, double norm_) {
            
            mean.push_back(mean_);
            sigma.push_back(sigma_);
            norm.push_back(norm_);
            totNorm = std::accumulate(norm.begin(), norm.end(), 0.0f); // recompute norm
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
            return b;
        }
    
};


// recoil correction
Vec_f recoilCorrectionBinned(double pU1, double pU2, double qTbinIdx) {
    
    Vec_f res(3, 0);
    
    res[1] = pU1;
    res[2] = pU2;
    res[0] = std::hypot(res[1], res[2]);
    
    TH3D *h3;
    
    GaussianSum u1_data;
    h3 = recoil_hists["para_data"];
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 2, 1), h3->GetBinContent(qTbinIdx+1, 4, 1));
    u1_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 3, 1), 1.-h3->GetBinContent(qTbinIdx+1, 4, 1));
    
    GaussianSum u2_data;
    h3 = recoil_hists["perp_data"];
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 2, 1), h3->GetBinContent(qTbinIdx+1, 4, 1));
    u2_data.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 3, 1), 1.-h3->GetBinContent(qTbinIdx+1, 4, 1));
    
    GaussianSum u1_mc;
    h3 = recoil_hists["para_mc"];
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 2, 1), h3->GetBinContent(qTbinIdx+1, 4, 1));
    u1_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 3, 1), 1.-h3->GetBinContent(qTbinIdx+1, 4, 1));
    
    GaussianSum u2_mc;
    h3 = recoil_hists["perp_mc"];
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 2, 1), h3->GetBinContent(qTbinIdx+1, 4, 1));
    u2_mc.addTerm(h3->GetBinContent(qTbinIdx+1, 1, 1), h3->GetBinContent(qTbinIdx+1, 3, 1), 1.-h3->GetBinContent(qTbinIdx+1, 4, 1));
    
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

double recoilParameterization_sigma(double a, double b, double c, double qT) {
    
    return a*std::pow(qT + b, c);
}

Vec_f recoilCorrectionParametric(double pU1, double pU2, double qT, int corrType=0, int corrIdx=1) {
    
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

        
    res[0] = sqrt(lMX*lMX + lMY*lMY);
    if(lMX > 0) res[1] = atan(lMY/lMX);
    else res[1] = (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX);
  
	return res;
}


}

#endif