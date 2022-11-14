#ifndef WREMNANTS_EtaPtCorrelatedEfficiency_h
#define WREMNANTS_EtaPtCorrelatedEfficiency_h

#include "TROOT.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TMath.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
//#include <stdlib.h>
//#include <stdio.h>
#include <cstdlib> 
#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>

namespace wrem {

    class EtaPtCorrelatedEfficiency {
  
    public:
  
        EtaPtCorrelatedEfficiency(TH3D* histocov, TH2D* histoerf, double ptmin, double ptmax);
        double DoEffSyst(int etabin, double pt, double *variations, bool getDiff=false);
        std::vector<double> DoEffSyst(int etabin, int ipar);
        void setPtRange(double ptmin, double ptmax) { ptmin_ = ptmin; ptmax_ = ptmax; }
        void setSmoothingFunction(const std::string& name);
    
    protected:

        Eigen::MatrixXd covariance(int etabin);
        void DoHessianShifts(int etabin, int ipar, double *inpars, double *outpars);

        std::string smoothFunctionName_ = "cheb3";
        TH3D *covhist_;
        TH2D *parhist_;
        int ndim_ = 4;
        double ptmin_ = 0.0;
        double ptmax_ = 100.0;
        TF1* function_ = nullptr;
        // list of predefined functions
        TF1* tf1_cheb3_ = new TF1("tf1_cheb3", "cheb3", ptmin_, ptmax_);
        TF1* tf1_cheb2_ = new TF1("tf1_cheb2", "cheb2", ptmin_, ptmax_);
    };


    EtaPtCorrelatedEfficiency::EtaPtCorrelatedEfficiency(TH3D* histocov, TH2D* histoerf, double ptmin, double ptmax) {
        covhist_ = histocov;
        int ny = covhist_->GetNbinsY();
        int nz = covhist_->GetNbinsZ();
        assert(ny==nz);
        ndim_ = ny;
        parhist_ = histoerf;
        ptmin_ = ptmin;
        ptmax_ = ptmax;
        setSmoothingFunction("cheb3");
    
    }

    void EtaPtCorrelatedEfficiency::setSmoothingFunction(const std::string& name) {
        smoothFunctionName_ = name;
        if (name.find("cheb3") != std::string::npos) {
            function_ = tf1_cheb3_;
            ndim_ = 4;   
        } else if (name.find("cheb2") != std::string::npos) {
            function_ = tf1_cheb2_;
            ndim_ = 3;
        } else {
            std::cout << "Smoothing function " << name << " not implemented. Abort" << std::endl;
            exit(EXIT_FAILURE);
        }
        return;
    }


    Eigen::MatrixXd EtaPtCorrelatedEfficiency::covariance(int etabin) {
        Eigen::MatrixXd covMat(ndim_, ndim_);
        for (int i = 0; i < ndim_; ++i) {
            for (int j = 0; j < ndim_; ++j) {
                covMat(i,j) = covhist_->GetBinContent(etabin, i+1, j+1);
            }
        }
        // std::cout << "covariance matrix = " << std::endl << covMat << std::endl;
        return covMat;
    }

    void EtaPtCorrelatedEfficiency::DoHessianShifts(int etabin, int ipar, double *inpars, double *outpars) {
    
        // diagonalize the covariance matrix
        Eigen::MatrixXd covMat = covariance(etabin);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(covMat);
        Eigen::VectorXd eigenv = es.eigenvalues();
        Eigen::MatrixXd transformation = es.eigenvectors();
        // std::cout << "Transformation = " << std::endl << transformation << std::endl;

        // transform the pars in the diagonal basis
        const unsigned int npars = transformation.rows();
        const unsigned int neigenvectors = transformation.cols(); 
        Eigen::VectorXd inparv(npars);
        for (unsigned int jpar = 0; jpar < npars; ++jpar) {
            inparv[ipar] = inpars[jpar];
        }
        // std::cout << "inparv = " << std::endl << inparv << std::endl;
        Eigen::VectorXd diagbasisv = transformation.transpose()*inparv;
        // std::cout << "diagbasisv = " << std::endl << diagbasisv << std::endl;

        // shift one of them by the diagonal uncertainty (uncorrelated in this basis)
        diagbasisv[ipar] += sqrt(eigenv[ipar]);

        // transform the pars back in the original basis
        Eigen::VectorXd outparv = transformation*diagbasisv;
        for (unsigned int ieig = 0; ieig < neigenvectors; ++ieig) {
            outpars[ieig] = outparv[ieig];
        }
        // std::cout << "outparv = " << std::endl << outparv << std::endl;
        return;
    }

    double EtaPtCorrelatedEfficiency::DoEffSyst(int etabin, double pt, double *variations, bool getDiff=false) {

        double inpars[ndim_], outpars[ndim_];
    
        for (int ipar = 0; ipar < ndim_; ++ipar) {
            inpars[ipar] = parhist_->GetBinContent(etabin, ipar+1);
        }
    
        double nominal = function_->EvalPar(&pt, inpars);
        
        for (int ipar = 0; ipar < ndim_; ++ipar) {
            DoHessianShifts(etabin, ipar, inpars, outpars);
            variations[ipar] = function_->EvalPar(&pt, outpars);
            if (getDiff) variations[ipar] -= nominal;
        }
        return nominal;
    }

    std::vector<double> EtaPtCorrelatedEfficiency::DoEffSyst(int etabin, int ipar) {

        if (ipar >= ndim_) {
            std::cout << "EtaPtCorrelatedEfficiency::DoEffSyst(int etabin, int ipar):  ipar >= " << ndim_ << " (" << ipar << ")" << std::endl;
            exit(EXIT_FAILURE);
        }
        double inpars[ndim_], outpars[ndim_];
    
        for (int jpar = 0; jpar < ndim_; ++jpar) {
            inpars[jpar] = parhist_->GetBinContent(etabin, jpar+1);
        }
        DoHessianShifts(etabin, ipar, inpars, outpars);
        std::vector<double> ret(ndim_, 0.0);
        for (int jpar = 0; jpar < ndim_; ++jpar) {
            ret[jpar] = outpars[jpar];
        }
        return ret;
    }

    
}
    
#endif
