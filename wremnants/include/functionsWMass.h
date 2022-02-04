#ifndef WREMNANTS_FUNCTIONS_WMASS_H
#define WREMNANTS_FUNCTIONS_WMASS_H

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include "defines.h"

namespace wrem {

TF1 * helicityFractionSimple_0 = new TF1("helicityFraction_0", "3./4*(TMath::Sqrt(1-x*x))^2", -1., 1.);
TF1 * helicityFractionSimple_L = new TF1("helicityFraction_L", "3./8.*(1-x)^2"              , -1., 1.);
TF1 * helicityFractionSimple_R = new TF1("helicityFraction_R", "3./8.*(1+x)^2"              , -1., 1.);

TFile *_file_helicityFractionsSimple = NULL;
TH2 * helicityFractionsSimple_0 = NULL;
TH2 * helicityFractionsSimple_L = NULL;
TH2 * helicityFractionsSimple_R = NULL;

string getEnvironmentVariable(const string& env_var_name = "CMSSW_BASE") {

  char* _env_var_ptr = getenv(env_var_name.c_str());
  if (_env_var_ptr == nullptr) {
    cout << "Error: environment variable " << env_var_name << " not found. Exit" << endl;
    exit(EXIT_FAILURE);
  } else {
    string str = string(_env_var_ptr);
    return str;
  }

}

bool etaptIsInBin(const TH2& h,
		  const int& ietabin, const int& iptbin,
		  const float& eta, const float& pt) {
  if (pt >= h.GetYaxis()->GetBinLowEdge(iptbin) and pt < h.GetYaxis()->GetBinUpEdge(iptbin) and \
      eta >= h.GetXaxis()->GetBinLowEdge(ietabin) and eta < h.GetXaxis()->GetBinUpEdge(ietabin)
      )
    return true;
  else
    return false;

}


double getValFromTH2(const TH2& h, const float& x, const float& y, const float& sumError=0.0) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  if (sumError)
    return h.GetBinContent(xbin, ybin) + sumError * h.GetBinError(xbin, ybin);
  else
    return h.GetBinContent(xbin, ybin);
}

double getValFromTH2bin(const TH2& h, const int& xbin, const int& ybin, const float& sumError=0.0) {
  if (sumError)
    return h.GetBinContent(xbin, ybin) + sumError * h.GetBinError(xbin, ybin);
  else
    return h.GetBinContent(xbin, ybin);
}

double getRelUncertaintyFromTH2(const TH2& h, const float& x, const float& y, const float valBadRatio = 1.0) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  if (h.GetBinContent(xbin, ybin) != 0.0)
    return h.GetBinError(xbin, ybin)/h.GetBinContent(xbin, ybin);
  else
    return valBadRatio;
}

double getRelUncertaintyFromTH2bin(const TH2& h, const int& xbin, const int& ybin, const float valBadRatio = 1.0) {
  if (h.GetBinContent(xbin, ybin) != 0.0)
    return h.GetBinError(xbin, ybin)/h.GetBinContent(xbin, ybin);
  else
    return valBadRatio;
}

double getAbsUncertaintyFromTH2(const TH2& h, const float& x, const float& y) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  return h.GetBinError(xbin, ybin);
}

double getAbsUncertaintyFromTH2bin(const TH2& h, const int& xbin, const int& ybin) {
  return h.GetBinError(xbin, ybin);
}


float returnChargeVal(float val1, int ch1, float val2, int ch2, ULong64_t evt){

    float retVal = -999.;

    if (evt%2 ) retVal = ch1 > 0 ? val1 : val2; //odd event numbers
    else        retVal = ch1 > 0 ? val2 : val1; //even event numbers

    return retVal;

}


float returnPlusVal(float val1, int ch1, float val2, int ch2){

  // return value of the lepton with desired charge, without looking at event parity
  // assumes opposite charges
  return (ch1 > 0) ? val1 : val2;

}

float returnMinusVal(float val1, int ch1, float val2, int ch2){

  // return value of the lepton with desired charge, without looking at event parity
  // assumes opposite charges
  return (ch1 < 0) ? val1 : val2;

}

float returnChargeValAllEvt(int desiredCharge, float val1, int ch1, float val2, int ch2, float returnForSameCharge = 0.0){

  if (ch1*ch2 > 0) {
    return returnForSameCharge;
  } else {
    if (desiredCharge > 0) return (ch1 > 0) ? val1 : val2;
    else                   return (ch1 < 0) ? val1 : val2;
  }

}


float mt_wlike(float pt1, float phi1, int ch1, float pt2, float phi2, int ch2, float met, float phimet, ULong64_t evt) {

  //if (ch1 == ch2) return mt_wlike_samesign(pt1,phi1,pt2,phi2,met,phimet);
  float ptL  = 0.0;
  float phiL = 0.0;

  float ptl  = 0.0;
  float phil = 0.0;

  // positive (negative) leptons on odd (even) events
  if (isOddEvent(evt)) {
    if (ch1 > 0) {
      ptL = pt1;
      phiL = phi1;
      ptl = pt2;
      phil = phi2;
    } else {
      ptL = pt2;
      phiL = phi2;
      ptl = pt1;
      phil = phi1;
    }
  } else {
    if (ch1 < 0) {
      ptL = pt1;
      phiL = phi1;
      ptl = pt2;
      phil = phi2;
    } else {
      ptL = pt2;
      phiL = phi2;
      ptl = pt1;
      phil = phi1;
    }
  }

  TVector2 pl = TVector2();
  pl.SetMagPhi(ptl,phil);

  TVector2 metv = TVector2();
  metv.SetMagPhi(met,phimet);
  TVector2 met_wlike = pl+metv;

  return std::sqrt(2*ptL*met_wlike.Mod()*(1-std::cos(phiL-met_wlike.Phi())));

}

float lepPlusMetPt(float pt, float phi, float met, float phimet) {

  TVector2 pl = TVector2();
  pl.SetMagPhi(pt, phi);

  TVector2 met_wlike = TVector2();
  met_wlike.SetMagPhi(met, phimet);
  met_wlike = pl + met_wlike;
  return met_wlike.Mod();

}

float mt_wlike_nano(float pt, float phi, float ptOther, float phiOther, float met, float phimet) {

  TVector2 pl = TVector2();
  pl.SetMagPhi(ptOther,phiOther);

  TVector2 met_wlike = TVector2();
  met_wlike.SetMagPhi(met,phimet);
  met_wlike = pl + met_wlike;

  return std::sqrt(2*pt*met_wlike.Mod()*(1-std::cos(phi-met_wlike.Phi())));

}

Vec_f lepPlusMetPt_allLep(Vec_f pt, Vec_f phi, float met, float phimet) {

  unsigned int size = 2;
  Vec_f res(size, 0.0); // 2 elements initialized to 0

  for (unsigned int i = 0; i < size; ++i) {

    unsigned int indexOther = size - i - 1;
    TVector2 pl = TVector2();
    pl.SetMagPhi(pt[indexOther], phi[indexOther]);
    TVector2 met_wlike = TVector2();
    met_wlike.SetMagPhi(met, phimet);
    met_wlike = pl + met_wlike;
    res[i] = met_wlike.Mod();

  }

  return res;

}

Vec_f mt_wlikeSS_nano(Vec_f pt, Vec_f phi, float met, float phimet) {

  unsigned int size = 2;
  Vec_f res(size, 0.0); // 2 elements initialized to 0

  for (unsigned int i = 0; i < size; ++i) {

    unsigned int indexOther = size - i - 1;
    TVector2 pl = TVector2();
    pl.SetMagPhi(pt[indexOther], phi[indexOther]);
    TVector2 met_wlike = TVector2();
    met_wlike.SetMagPhi(met, phimet);
    met_wlike = pl + met_wlike;
    res[i] = std::sqrt(2 * pt[i] * met_wlike.Mod() * (1-std::cos(phi[i] - met_wlike.Phi())));

  }

  return res;

}

float helicityWeightSimple(float yw, float ptw, float costheta, int pol)
{

  if (!helicityFractionsSimple_0 || !helicityFractionsSimple_L || !helicityFractionsSimple_R) {
    _file_helicityFractionsSimple = new TFile("w-mass-13TeV/fractionReweighting/fractions.root","read");
    helicityFractionsSimple_0 = (TH2F*)(_file_helicityFractionsSimple->Get("fraction0_plus_sym"));
    helicityFractionsSimple_L = (TH2F*)(_file_helicityFractionsSimple->Get("fractionL_plus_sym"));
    helicityFractionsSimple_R = (TH2F*)(_file_helicityFractionsSimple->Get("fractionR_plus_sym"));
  }

  if (std::abs(costheta) > 1.) {
    //std::cout << " found an event with weird cosTheta = " << costheta << std::endl;
    //std::cout << " setting event weight to 0" << std::endl;
    return 0;
  }

  TH2 *hist_f0 = helicityFractionsSimple_0;
  TH2 *hist_fL = helicityFractionsSimple_L;
  TH2 *hist_fR = helicityFractionsSimple_R;

  // float yval  = std::abs(yw) > hist_f0->GetXaxis()->GetXmax() ? hist_f0->GetXaxis()->GetXmax() : yw;
  // float ptval = ptw > hist_f0->GetYaxis()->GetXmax() ? hist_f0->GetYaxis()->GetXmax() : ptw;

  int ywbin = std::max(1, std::min(hist_f0->GetNbinsX(), hist_f0->GetXaxis()->FindBin(yw )));
  int ptbin = std::max(1, std::min(hist_f0->GetNbinsY(), hist_f0->GetYaxis()->FindBin(ptw)));

  float f0 = hist_f0->GetBinContent(ywbin, ptbin);
  float fL = hist_fL->GetBinContent(ywbin, ptbin);
  float fR = hist_fR->GetBinContent(ywbin, ptbin);

  float f0Term = helicityFractionSimple_0->Eval(costheta);
  float fLTerm = helicityFractionSimple_L->Eval(costheta);
  float fRTerm = helicityFractionSimple_R->Eval(costheta);

  if      (pol == 0) return f0*f0Term/(f0*f0Term+fL*fLTerm+fR*fRTerm);
  else if (pol == 1) return fL*fLTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm);
  else if (pol == 2) return fR*fRTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm);

  std::cout << "something went wrong in the helicity reweighting" << std::endl;
  return -99999.;

}

float mydeltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float mydeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1-eta2;
  float dphi = mydeltaPhi(phi1,phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}


//-------------------

TFile *_file_ratio_FSRoverNoFSR_etaPt_mu = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_mu = NULL;
TFile *_file_ratio_FSRoverNoFSR_etaPt_el = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_el = NULL;

// old function, should not be needed this time
float FSRscaleFactorEtaPt(int pdgId, float dresspt, float dresseta) {

  // these histograms are the gen-level xsec before any cut (except the gen_decayId)
  // the yields and ratio of fsr and no-fsr (both on top of W-pt reweighting) are here:
  // muon: http://mciprian.web.cern.ch/mciprian/wmass/13TeV/distribution/FSR_atGenLevel_muon_genBinAnalysis/
  // electron: http://mciprian.web.cern.ch/mciprian/wmass/13TeV/distribution/FSR_atGenLevel_electron_genBinAnalysis/
  // FSR should not change the gen-level xsec, but our fsr weights can artificially do it
  // so, we use the ratio to rescale the fsr weight in bins of gen-level pt and eta
  // the histogram uses the binning for the 2D xsec measurement (as each of them should keep the same xsec as without fsr)


  TH2D* hratio_FSR_noFSR = NULL;
  Double_t outlierWgt = 0.0;  // hardcoded below, taken from ratio of yields outside acceptance

  if (abs(pdgId)==11) {

    if (!ratio_FSRoverNoFSR_etaPt_el) {
      _file_ratio_FSRoverNoFSR_etaPt_el = new TFile("w-mass-13TeV/theoryReweighting/FSR_atGenLevel_electron_genEtaPtAnalysis.root","read");
      ratio_FSRoverNoFSR_etaPt_el = (TH2D*) _file_ratio_FSRoverNoFSR_etaPt_el->Get("ratio__Wzpt_el__Wfsr_el");
    }
    hratio_FSR_noFSR = ratio_FSRoverNoFSR_etaPt_el;
    outlierWgt = 0.964796;

  } else if (abs(pdgId)==13) {

     if (!ratio_FSRoverNoFSR_etaPt_mu) {
      _file_ratio_FSRoverNoFSR_etaPt_mu = new TFile("w-mass-13TeV/theoryReweighting/FSR_atGenLevel_muon_genEtaPtAnalysis.root","read");
      ratio_FSRoverNoFSR_etaPt_mu = (TH2D*) _file_ratio_FSRoverNoFSR_etaPt_mu->Get("ratio__Wzpt_mu__Wfsr_mu");
    }
    hratio_FSR_noFSR = ratio_FSRoverNoFSR_etaPt_mu;
    outlierWgt = 0.914322;

  } else {

    return 1.0;

  }

  int etabin = hratio_FSR_noFSR->GetXaxis()->FindFixBin(fabs(dresseta));
  int ptbin = hratio_FSR_noFSR->GetXaxis()->FindFixBin(fabs(dresspt));
  if (ptbin == 0 or ptbin > hratio_FSR_noFSR->GetNbinsY() or etabin == 0 or etabin> hratio_FSR_noFSR->GetNbinsX()) {
    return outlierWgt;
  } else {
    return 1./hratio_FSR_noFSR->GetBinContent(etabin,ptbin); // ratio is FSR/no-FSR, but need the opposite
  }

}

//-------------------

TFile *_file_fsrWeights_simple = NULL;
TH1F * fsrWeights_el = NULL;
TH1F * fsrWeights_mu = NULL;

float fsrPhotosWeightSimple(int pdgId, float dresspt, float barept, bool normToSameGenArea = false, float dresseta = 0) {

  if (!fsrWeights_el || !fsrWeights_mu) {
    _file_fsrWeights_simple = new TFile("w-mass-13TeV/theoryReweighting/photos_rwgt_integrated.root","read");
    fsrWeights_el  = (TH1F*)(_file_fsrWeights_simple->Get("w_e_h_lptBareOverDressed_ratio"));
    fsrWeights_mu  = (TH1F*)(_file_fsrWeights_simple->Get("w_mu_h_lptBareOverDressed_ratio"));
  }
  TH1F *fsrWeights = 0;
  if      (abs(pdgId)==11) fsrWeights = fsrWeights_el;
  else if (abs(pdgId)==13) fsrWeights = fsrWeights_mu;
  else return 1;

  int ratiobin  = std::max(1, std::min(fsrWeights->GetNbinsX(), fsrWeights->GetXaxis()->FindFixBin(barept/dresspt)));

  // normfactor is used to keep the total xsec unchanged when using the fsr
  // it was obtained from the gen level cross section for e/mu without any cut (except the gen-pdgid), as the ratio of nomi/fsr (nomi/fsr < 1)
  //Double_t normfactor = (abs(pdgId) == 11) ? 0.964792 : 0.914320;
  Double_t normfactor = (normToSameGenArea) ? FSRscaleFactorEtaPt(pdgId, dresspt, dresseta) : 1.0;
  return normfactor * fsrWeights->GetBinContent(ratiobin);

}


TFile *_file_fsrWeights = NULL;
TH3F * fsrWeights_elplus  = NULL;
TH3F * fsrWeights_elminus = NULL;
TH3F * fsrWeights_muplus  = NULL;
TH3F * fsrWeights_muminus = NULL;

float fsrPhotosWeight(int pdgId, float dresseta, float dresspt, float barept) {
  if (!fsrWeights_elplus || !fsrWeights_elminus || !fsrWeights_muplus || !fsrWeights_muminus) {
    _file_fsrWeights = new TFile("w-mass-13TeV/theoryReweighting/photos_rwgt.root","read");
    fsrWeights_elplus  = (TH3F*)(_file_fsrWeights->Get("qed_weights_wp_e"));
    fsrWeights_elminus = (TH3F*)(_file_fsrWeights->Get("qed_weights_wm_e"));
    fsrWeights_muplus  = (TH3F*)(_file_fsrWeights->Get("qed_weights_wp_mu"));
    fsrWeights_muminus = (TH3F*)(_file_fsrWeights->Get("qed_weights_wm_mu"));
  }
  TH3F *fsrWeights = 0;
  if      (abs(pdgId)==11) fsrWeights = ( pdgId>0 ? fsrWeights_elplus : fsrWeights_elminus );
  else if (abs(pdgId)==13) fsrWeights = ( pdgId>0 ? fsrWeights_muplus : fsrWeights_muminus );
  else return 1;

  int etabin = std::max(1, std::min(fsrWeights->GetNbinsX(), fsrWeights->GetXaxis()->FindFixBin(fabs(dresseta))));
  int ptbin  = std::max(1, std::min(fsrWeights->GetNbinsY(), fsrWeights->GetYaxis()->FindFixBin(dresspt)));
  int zbin  = std::max(1, std::min(fsrWeights->GetNbinsZ(), fsrWeights->GetZaxis()->FindFixBin(barept/dresspt)));

  return fsrWeights->GetBinContent(etabin,ptbin,zbin);
}

TFile *_file_dyptWeights = NULL;
TH1F *amcnlody = NULL;

float dyptWeight(float pt2l, int isZ, bool scaleNormWToGenXsecBeforeCuts = false) {
  if (!amcnlody) {
    _file_dyptWeights = new TFile("w-mass-13TeV/theoryReweighting/zpt_weights.root");
    amcnlody = (TH1F*)(_file_dyptWeights->Get("amcnlo"));
  }
  int ptbin = std::max(1, std::min(amcnlody->GetNbinsX(), amcnlody->GetXaxis()->FindFixBin(pt2l)));
  // for the Z, change the pT *and* normalization to the measured one in data
  // for the W, change only the pT, but scale back to the total xsec of MC@NLO (factor comptued for total xsec in acceptance)
  // when using pt range 26-56 instead of 26-45, there is an additional factor 0.987, so the actual number would be 0.946

  // when using only 0.958 I see that total gen-xsec no-Wpt over with-Wpt is 1.0138013 for W+ and 1.0144267 for W-
  // let's take only 1.014 without distinguishing charges
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  float scaleW = scaleNormWToGenXsecBeforeCuts ? 0.9714120 : 0.958;  //1.014 * 0.958
  float scaleToMCaNLO = isZ ? 1. : scaleW;
  // plots are MC/data
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  return scaleToMCaNLO / amcnlody->GetBinContent(ptbin);
}
//=================
TFile *_file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO = NULL;
TH1F *zptRatio_MGaMCAtNLO_PowhegMiNNLO = NULL;

float dyptWeight_PowhegMiNNLO(float pt2l, int isZ, bool scaleNormWToGenXsecBeforeCuts = false, bool usePhotos = false) {
  if (!amcnlody) {
    _file_dyptWeights = new TFile("w-mass-13TeV/theoryReweighting/zpt_weights.root");
    amcnlody = (TH1F*)(_file_dyptWeights->Get("amcnlo"));
  }
  int ptbin = std::max(1, std::min(amcnlody->GetNbinsX(), amcnlody->GetXaxis()->FindFixBin(pt2l)));
  // for the Z, change the pT *and* normalization to the measured one in data
  // for the W, change only the pT, but scale back to the total xsec of MC@NLO (factor comptued for total xsec in acceptance)
  // when using pt range 26-56 instead of 26-45, there is an additional factor 0.987, so the actual number would be 0.946

  // when using only 0.958 I see that total gen-xsec no-Wpt over with-Wpt is 1.0138013 for W+ and 1.0144267 for W-
  // let's take only 1.014 without distinguishing charges
  //float scaleToMCaNLO = isZ ? 1. : 0.958;

  if (!zptRatio_MGaMCAtNLO_PowhegMiNNLO) {
    // 2 GeV granularity from 0 to 100 GeV (last bin is overflow)
    _file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO = new TFile("w-mass-13TeV/theoryReweighting/ZptRatio_MGaMCAtNLO_PowhegMiNNLO_dressed_noCuts.root");
    string hname_MGaMCAtNLO_PowhegMiNNLO = "pythia8";
    if (usePhotos) hname_MGaMCAtNLO_PowhegMiNNLO = "pythia8_photos";
    zptRatio_MGaMCAtNLO_PowhegMiNNLO = (TH1F*)(_file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO->Get(hname_MGaMCAtNLO_PowhegMiNNLO.c_str()));
  }
  int ptbinCorr = std::max(1, std::min(zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetNbinsX(),
				       zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetXaxis()->FindFixBin(pt2l))
			   );

  float scaleW = scaleNormWToGenXsecBeforeCuts ? 0.9714120 : 0.958;  //1.014 * 0.958 // to be revisited with new MC
  float scaleToMCaNLO = isZ ? 1. : scaleW;
  // plots are MC/data
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  return zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetBinContent(ptbinCorr) * scaleToMCaNLO / amcnlody->GetBinContent(ptbin);
}


//=================

TFile *_file_postfitWeights = NULL;
TH1F *hist_scales_0plus  = NULL;
TH1F *hist_scales_Lplus  = NULL;
TH1F *hist_scales_Rplus  = NULL;
TH1F *hist_scales_0minus = NULL;
TH1F *hist_scales_Lminus = NULL;
TH1F *hist_scales_Rminus = NULL;

float postfitQCDWeight(float pt2l, int pol, int charge, int flav=0) {
  TH1F* hist_scales = NULL;
  if (!_file_postfitWeights) {
    if (flav==0)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitlep.root", "read");
    else if (flav==13)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitmu.root", "read");
    else if (flav==11)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitel.root","read");
    else {
      std::cout << "ERROR! Unknown flavor: " << flav << " Returning 0 weight." << std::endl;
      return 0;
    }
    hist_scales_0plus  = (TH1F*)(_file_postfitWeights->Get("weights_longplus"));
    hist_scales_Lplus  = (TH1F*)(_file_postfitWeights->Get("weights_leftplus"));
    hist_scales_Rplus  = (TH1F*)(_file_postfitWeights->Get("weights_rightplus"));
    hist_scales_0minus = (TH1F*)(_file_postfitWeights->Get("weights_longminus"));
    hist_scales_Lminus = (TH1F*)(_file_postfitWeights->Get("weights_leftminus"));
    hist_scales_Rminus = (TH1F*)(_file_postfitWeights->Get("weights_rightminus"));
  }
  if (charge>0) {
    if (pol==0)      hist_scales = hist_scales_0plus;
    else if (pol==1) hist_scales = hist_scales_Lplus;
    else if (pol==2) hist_scales = hist_scales_Rplus;
    else std::cerr << "ERROR: polarization " << pol << " not defined for postfitQCDWeight()" << std::endl;
  } else {
    if (pol==0)      hist_scales = hist_scales_0minus;
    else if (pol==1) hist_scales = hist_scales_Lminus;
    else if (pol==2) hist_scales = hist_scales_Rminus;
    else std::cerr << "ERROR: polarization " << pol << " not defined for postfitQCDWeight()" << std::endl;
  }
  int ptbin = std::max(1, std::min(hist_scales->GetNbinsX(), hist_scales->GetXaxis()->FindFixBin(pt2l)));
  return hist_scales->GetBinContent(ptbin);
}


//==================================================
TRandom3 *rng = NULL;

float getSmearedVar(float var, float smear, ULong64_t eventNumber, int isData, bool smearOnlyMC=false) {

  if (smearOnlyMC && isData) return var;

  if(!rng) rng = new TRandom3();
  // use eventNumber as seed, otherwise each time the function is called for the same event, the smearer produce a different smeared value
  rng->SetSeed(eventNumber); // make it really random across different jobs
  return var * ( 1.0 + smear * rng->Gaus());

}

//==================================================

// Sorry you have to manually keep these consistent
std::unordered_map<DataEra, std::string> eraNames = { {BToF, "BtoF"}, {GToH, "GtoH"} };
std::unordered_map<DataType, std::string> datatypeNames = { {MC, "MC"}, {Data, "Data"} };
std::unordered_map<ScaleFactorType, std::string> scalefactorNames = { {isoTrigPlus, "isoTrigPlus"}, {isoTrigMinus, "isoTrigMinus"}, {isoNotrig, "isoNotrig"}, {noisoTrigPlus, "noisoTrigPlus"}, {noisoTrigMinus, "noisoTrigMinus"}, {noisoNotrig, "noisoNotrig"}, {antiisoTrigPlus, "antiisoTrigPlus"}, {antiisoTrigMinus, "antiisoTrigMinus"}, {antiisoNotrig, "antiisoNotrig"}, {isoOnly, "isoOnly"}, {antiisoOnly, "antiisoOnly"}, {isoNotrigOnly, "isoNotrigOnly"}, {antiisoNotrigOnly, "antiisoNotrigOnly"}, {tracking, "tracking"}, {reco, "reco"}, {trackingReco, "trackingReco"} };
// FOR TESTS WITH EFFICIENCIES (F is preVFP part)
std::unordered_map<DataEra, std::string> runEraNames = { {GToH, "GtoH"}, {BToF, "BtoF"}, {B, "B"}, {C, "C"}, {D, "D"}, {E, "E"}, {F, "F"}, {G, "G"}, {H, "H"} };

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &p) const
    {
      std::size_t seed = 0;
      boost::hash_combine(seed, p.first);
      boost::hash_combine(seed, p.second);
      return seed;
    }
};

std::unordered_map<DataEra, TH2D> hMuonPrefiringNew = {}; // this has the modelling of prefiring versus pt
TH2D* hMuonPrefiringNew_hotspot = nullptr;
std::unordered_map<DataEra, TH1F> hMuonPrefiring = {}; // will store pre and post (only BToF and GToH)
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> scaleFactorHist = {};
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> scaleFactorHist_dataAltSig = {};


void initializeScaleFactors(const string &datadir, const string& _filename_allSF = "./testMuonSF/scaleFactorProduct_31Mar2021.root", const bool oldSFname = false) {

  //TFile _file_allSF_hardcoded = TFile("/testMuonSF/scaleFactorProduct_31May2021_nodz_dxybs.root", "read");
  TFile _file_allSF = TFile(_filename_allSF.c_str(), "read");
  if (!_file_allSF.IsOpen()) {
    std::cerr << "WARNING: Failed to open scaleFactors file " << _filename_allSF << "! No scale factors will be applied\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "INFO >>> Initializing histograms for SF from file " << _filename_allSF << std::endl;

  // nominal scale factors
  for (auto& corr : scalefactorNames) {
    for (auto& era : eraNames) {
      std::string sfNamePrefix = oldSFname ? "fullSF2D" : "fullSF2D_nominal";
      std::vector<std::string> vars = {sfNamePrefix, corr.second, era.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr) {
          std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
                    << _filename_allSF << "! Aborting" << std::endl;
          exit(EXIT_FAILURE);
      }
      histptr->SetDirectory(0);
      DataEra eraVal = era.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, eraVal);
      scaleFactorHist[corrKey] = *dynamic_cast<TH2F*>(histptr);
    }

  }

  if (not oldSFname) {

      // now all alternative scale factors where data was fitted with analytic function rather than MC template
      for (auto& corr : scalefactorNames) {
          for (auto& era : eraNames) {
              std::vector<std::string> vars = {"fullSF2D_dataAltSig", corr.second, era.second};
              std::string corrname = boost::algorithm::join(vars, "_");
              auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
              if (histptr == nullptr) {
                  std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
                            << _filename_allSF << "! Aborting" << std::endl;
                  exit(EXIT_FAILURE);
              }
              histptr->SetDirectory(0);
              DataEra eraVal = era.first;
              ScaleFactorType key = corr.first;
              // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
              auto corrKey = std::make_pair(key, eraVal);
              scaleFactorHist_dataAltSig[corrKey] = *dynamic_cast<TH2F*>(histptr);
          }

      }

  }

  _file_allSF.Close(); // should work since we used TH1F::SetDirectory(0) to detach histogram from file

  // new prefiring from Jan
  std::string _filename_prefiringNew = datadir + "/testMuonSF/L1MuonPrefiringParametriations_histograms.root";
  TFile _file_prefiringNew = TFile(_filename_prefiringNew.c_str(), "read");
  if (!_file_prefiringNew.IsOpen()) {
    std::cerr << "WARNING: Failed to open prefiring file " << _filename_prefiringNew << "\n";
    exit(EXIT_FAILURE);
  }
  std::cout << "INFO >>> Initializing histograms for prefiring from file " << _filename_prefiringNew << std::endl;
  for (auto& era : runEraNames) {
    // std::cout << era.second << std::endl;
    if (era.first == H) {
      hMuonPrefiringNew[era.first] = *(dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016H")));
      hMuonPrefiringNew[era.first].SetDirectory(0);
    } else if (era.first != GToH) {
      // BG should be like preVFP, but more data was used to derive corrections
      //hMuonPrefiringNew[era.first] = *(dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016preVFP")));
      hMuonPrefiringNew[era.first] = *(dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016BG")));
      hMuonPrefiringNew[era.first].SetDirectory(0);
    } else {
      hMuonPrefiringNew[GToH] = *(dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016postVFP")));
      hMuonPrefiringNew[GToH].SetDirectory(0);
    }
  }
  hMuonPrefiringNew_hotspot = (dynamic_cast<TH2D*>(_file_prefiringNew.Get("L1prefiring_muonparam_2016_hotspot")));
  hMuonPrefiringNew_hotspot->SetDirectory(0);
  _file_prefiringNew.Close();

}

////=====================================================================================

std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2D, pair_hash> efficiencyMCtruthPerEra = {};
std::unordered_map<ScaleFactorType, std::string> efficiencyNames = { {isoTrigPlus, "idipANDtrigANDiso"}, {isoTrigMinus, "idipANDtrigANDiso"}, {isoNotrig, "idipANDisonotrig"}, {noisoTrigPlus, "idipANDtrig"}, {noisoTrigMinus, "idipANDtrig"}, {noisoNotrig, "idip"}, {trackingReco, "trackerOrGlobal"}, {isoTrigPlusStepOfTnp, "isoTrigPlusProdStepOfTnP"}, {isoNotrigStepOfTnP, "isoNotrigProdStepOfTnP"}, {global, "global"} };
// std::unordered_map<ScaleFactorType, std::string> efficiencyNamesAsTNP = { {isoTrigPlus, {"idipStepOfTnP", "triggerStepOfTnP", "isoStepOfTnP"}},
// 									  {isoTrigMinus, {"idipStepOfTnP", "triggerStepOfTnP", "isoStepOfTnP"}},
// 									  {isoNotrig, {"idipStepOfTnP", "isonotrigStepOfTnP"}},
// 									  {noisoTrigPlus, {"idipStepOfTnP", "triggerStepOfTnP"}},
// 									  {noisoTrigMinus, {"idipStepOfTnP", "triggerStepOfTnP"}},
// 									  {noisoNotrig, {"idipStepOfTnP"}} };

// for test, the trigger efficiency is for plus only at the moment, and all are made with no dz cut and using dxybs
void initializeEfficiencyMCtruth(const string& _filename_allSF = "./testMuonSF/mcTruthEff.root") {

  TFile _file_allSF = TFile(_filename_allSF.c_str(), "read");
  if (!_file_allSF.IsOpen()) {
    std::cerr << "WARNING: Failed to open scaleFactors file " << _filename_allSF << "! No scale factors will be applied\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "INFO >>> Initializing histograms for MC truth efficiency from file " << _filename_allSF << std::endl;

  for (auto& corr : efficiencyNames) {

    for (auto& dataRunEra : runEraNames) {
      std::string dataRunEraStr = dataRunEra.second;
      // adjust some names where the convention was not always consistent
      if (dataRunEraStr.find("GtoH") != std::string::npos) dataRunEraStr = "GToH";
      if (dataRunEraStr.find("BtoF") != std::string::npos) dataRunEraStr = "BToF";
      std::vector<std::string> vars = {"mcTruthEff", corr.second, dataRunEraStr};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = dynamic_cast<TH2D*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr) {
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! Will continue neglecting them, but be careful" << std::endl;
	// exit(EXIT_FAILURE); // these are test SF and sometimes we miss some of them, but as long as we don't try to use them it is fine
	continue;
      }
      histptr->SetDirectory(0);
      DataEra typeVal = dataRunEra.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, typeVal);
      efficiencyMCtruthPerEra[corrKey] = *dynamic_cast<TH2D*>(histptr);
    }

  }

}

double _get_MCtruthEffPerEra(float pt,      float eta,      int charge,
			     float ptOther, float etaOther,
			     DataEra dtype = C,
			     bool isoSF1 = true, // to use SF for iso or antiiso
			     bool isoSF2 = true,
			     bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
			     ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (neglectIso) {
    sftype = charge > 0 ? noisoTrigPlus : noisoTrigMinus;
  } else if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, dtype);
  const TH2D& hcorr = efficiencyMCtruthPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);
  if (sf <= 0.0) {
    std::cout << "MC truth eff = " << sf << " for era " << runEraNames[dtype] << "   pt = " << pt << "    eta = " << eta << "   charge = " << charge << std::endl;
  }
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    if (neglectIso) {
      sftype = noisoNotrig;
    } else {
      sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    }
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2D& hcorrOther = efficiencyMCtruthPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_MCtruthEffAsTnpPerEra(float pt,      float eta,      int charge,
				  float ptOther, float etaOther,
				  DataEra dtype = C
				  ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlusStepOfTnp : isoTrigPlusStepOfTnp; // only charge plus available for now

  auto const key = std::make_pair(sftype, dtype);
  const TH2D& hcorr = efficiencyMCtruthPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);
  if (sf <= 0.0) {
    std::cout << "MC truth eff = " << sf << " for era " << runEraNames[dtype] << "   pt = " << pt << "    eta = " << eta << "   charge = " << charge << std::endl;
  }
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    sftype = isoNotrig;
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2D& hcorrOther = efficiencyMCtruthPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}


double _get_MCtruthRecoAndTrackingEffPerEra(float pt,      float eta,      int charge,
					    float ptOther, float etaOther,
					    DataEra dtype = C
					    ) {

  ScaleFactorType sftype = trackingReco;

  auto const key = std::make_pair(sftype, dtype);
  const TH2D& hcorr = efficiencyMCtruthPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);

  if (ptOther > 0.0) {
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2D& hcorrOther = efficiencyMCtruthPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_MCtruthGlobalMuonEffPerEra(float pt,      float eta,      int charge,
				       float ptOther, float etaOther,
				       DataEra dtype = C
				       ) {

  ScaleFactorType sftype = global;

  auto const key = std::make_pair(sftype, dtype);
  const TH2D& hcorr = efficiencyMCtruthPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);

  if (ptOther > 0.0) {
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2D& hcorrOther = efficiencyMCtruthPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

//================

//std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> scaleFactorDataPerEra = {};
//std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> scaleFactorMCPerEra = {};
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> efficiencyDataPerEra = {};
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> efficiencyMCPerEra = {};
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2F, pair_hash> scaleFactorPerEra = {};
std::unordered_map<ScaleFactorType, TH2F> lumiAverageScaleFactorPreVFP = {};

std::unordered_map<ScaleFactorType, std::string> scalefactorNamesTest = { {isoTrigPlus, "isoTrigPlus"}, {isoTrigMinus, "isoTrigMinus"}, {isoNotrig, "isoNotrig"}, {noisoTrigPlus, "noisoTrigPlus"}, {noisoTrigMinus, "noisoTrigMinus"}, {noisoNotrig, "noisoNotrig"}, {antiisoTrigPlus, "antiisoTrigPlus"}, {antiisoTrigMinus, "antiisoTrigMinus"}, {antiisoNotrig, "antiisoNotrig"}, {trackingReco, "trackingReco"}, {tracking, "tracking"}, {reco, "reco"} };


void initializeScaleFactorsTest(const string& _filename_allSF = "./testMuonSF/productEffAndSFperEra_nodz_dxybs.root") {

  TFile _file_allSF = TFile(_filename_allSF.c_str(), "read");
  if (!_file_allSF.IsOpen()) {
    std::cerr << "WARNING: Failed to open scaleFactors file " << _filename_allSF << "! No scale factors will be applied\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "INFO >>> Initializing histograms for test SF from file " << _filename_allSF << std::endl;

  // data/data
  for (auto& corr : scalefactorNamesTest) {

    for (auto& dataRunEra : runEraNames) {

      std::vector<std::string> vars = {"fullEffData2D", corr.second, dataRunEra.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr) {
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! Will continue neglecting them, but be careful" << std::endl;
	// exit(EXIT_FAILURE); // these are test SF and sometimes we miss some of them, but as long as we don't try to use them it is fine
	continue;
      }
      histptr->SetDirectory(0);
      DataEra typeVal = dataRunEra.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, typeVal);
      efficiencyDataPerEra[corrKey] = *dynamic_cast<TH2F*>(histptr);
    }

    // now for MC/MC
    for (auto& mcRunEra : runEraNames) {

      std::vector<std::string> vars = {"fullEffMC2D", corr.second, mcRunEra.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr) {
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! Will continue neglecting them, but be careful" << std::endl;
	// exit(EXIT_FAILURE); // these are test SF and sometimes we miss some of them, but as long as we don't try to use them it is fine
	continue;
      }
      histptr->SetDirectory(0);
      DataEra typeVal = mcRunEra.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, typeVal);
      efficiencyMCPerEra[corrKey] = *dynamic_cast<TH2F*>(histptr);
    }


    // repeat for data/MC
    for (auto& dataRunEra : runEraNames) {
      std::vector<std::string> vars = {"fullSF2D", corr.second, dataRunEra.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr) {
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! Will continue neglecting them, but be careful" << std::endl;
	// exit(EXIT_FAILURE); // these are test SF and sometimes we miss some of them, but as long as we don't try to use them it is fine
	continue;
      }
      histptr->SetDirectory(0);
      DataEra typeVal = dataRunEra.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, typeVal);
      scaleFactorPerEra[corrKey] = *dynamic_cast<TH2F*>(histptr);
    }

    // std::string corrname = Form("lumiAveScaleFactor_%s_BtoF", corr.second.c_str());
    // auto* histptr = dynamic_cast<TH2F*>(_file_allSF.Get(corrname.c_str()));
    // if (histptr == nullptr) {
    //   std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
    // 		<< _filename_allSF << "! Will continue neglecting them, but be careful" << std::endl;
    //   // exit(EXIT_FAILURE); // these are test SF and sometimes we miss some of them, but as long as we don't try to use them it is fine
    //   continue;
    // }
    // histptr->SetDirectory(0);
    // ScaleFactorType corrKey = corr.first;
    // lumiAverageScaleFactorPreVFP[corrKey] = *dynamic_cast<TH2F*>(histptr);

  }

}

double _get_TnpRecoAndTrackingEffMCPerEra(float pt,      float eta,      int charge,
					  float ptOther, float etaOther,
					  DataEra dtype = G
					  ) {

  ScaleFactorType sftype = trackingReco;

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = efficiencyMCPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);

  if (ptOther > 0.0) {
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = efficiencyMCPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_TnpRecoAndTrackingEffDataPerEra(float pt,      float eta,      int charge,
					    float ptOther, float etaOther,
					    DataEra dtype = G
					    ) {

  ScaleFactorType sftype = trackingReco;

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = efficiencyDataPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);

  if (ptOther > 0.0) {
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = efficiencyDataPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_TnpRecoAndTrackingSFPerEra(float pt,      float eta,      int charge,
				       float ptOther, float etaOther,
				       DataEra dtype = G
				       ) {

  ScaleFactorType sftype = trackingReco;

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = scaleFactorPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);

  if (ptOther > 0.0) {
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = scaleFactorPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_fullMuonDataEfficiencyEra(float pt,      float eta,      int charge,
				      float ptOther, float etaOther,
				      DataEra dtype = C,
				      bool isoSF1 = true, // to use SF for iso or antiiso
				      bool isoSF2 = true,
				      bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
				      ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (neglectIso) {
    sftype = charge > 0 ? noisoTrigPlus : noisoTrigMinus;
  } else if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = efficiencyDataPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    if (neglectIso) {
      sftype = noisoNotrig;
    } else {
      sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    }
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = efficiencyDataPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_fullMuonMCEfficiencyEra(float pt,      float eta,      int charge,
				    float ptOther, float etaOther,
				    DataEra dtype = C,
				    bool isoSF1 = true, // to use SF for iso or antiiso
				    bool isoSF2 = true,
				    bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
				    ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (neglectIso) {
    sftype = charge > 0 ? noisoTrigPlus : noisoTrigMinus;
  } else if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = efficiencyMCPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    if (neglectIso) {
      sftype = noisoNotrig;
    } else {
      sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    }
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = efficiencyMCPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}


double _get_fullMuonDataEfficiencyEraWithTrackingReco(float pt,      float eta,      int charge,
				    float ptOther, float etaOther,
				    DataEra dtype = C,
				    bool isoSF1 = true, // to use SF for iso or antiiso
				    bool isoSF2 = true,
				    bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
				    ) {

    return _get_fullMuonDataEfficiencyEra(pt, eta, charge, ptOther, etaOther, dtype, isoSF1, isoSF2, neglectIso) *_get_TnpRecoAndTrackingEffDataPerEra(pt, eta, charge, ptOther, etaOther, dtype);

}

double _get_fullMuonSF_perDataEra(float pt,      float eta,      int charge,
				  float ptOther, float etaOther,
				  DataEra dtype = B,
				  bool isoSF1 = true, // to use SF for iso or antiiso
				  bool isoSF2 = true,
				  bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
				  ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (neglectIso) {
    sftype = charge > 0 ? noisoTrigPlus : noisoTrigMinus;
  } else if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, dtype);
  const TH2F& hcorr = scaleFactorPerEra.at(key);
  double sf = getValFromTH2(hcorr, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    if (neglectIso) {
      sftype = noisoNotrig;
    } else {
      sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    }
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2F& hcorrOther = scaleFactorPerEra.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

double _get_lumiAverageScaleFactorPreVFP(float pt,      float eta,      int charge,
					 float ptOther, float etaOther,
					 bool isoSF1 = true, // to use SF for iso or antiiso
					 bool isoSF2 = true,
					 bool neglectIso = false // to neglect iso on both legs, overriding isoSF1 and isoSF2
					 ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (neglectIso) {
    sftype = charge > 0 ? noisoTrigPlus : noisoTrigMinus;
  } else if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  const TH2F& hcorr = lumiAverageScaleFactorPreVFP.at(sftype);
  double sf = getValFromTH2(hcorr, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    if (neglectIso) {
      sftype = noisoNotrig;
    } else {
      sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    }
    const TH2F& hcorrOther = lumiAverageScaleFactorPreVFP.at(sftype);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}

////=====================================================================================

// may just add this in a function, to avoid defining a column
Vec_b prefirableMuon(const Vec_f& pt, const Vec_b& looseId) {

  Vec_b res(pt.size(),false); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    if (pt[i] < 22) continue;
    if (not looseId[i]) continue;
    res[i] = true;
  }
  return res;

}

double _get_newMuonPrefiringSF(const Vec_f& eta, const Vec_f& pt, const Vec_f& phi, const Vec_b& looseId, DataEra era = BToF) {

  double sf = 1.0;
  // std::cout << "PREFIRING FOR: " << eraNames[era] << std::endl;

  const TH2D& hprefire = hMuonPrefiringNew[era];
  int nBins = hprefire.GetNbinsX();
  int prefireBin = 0;
  double prefiringProbability = 0.0;
  for (unsigned int i = 0; i < eta.size(); ++i) {
    if (not looseId[i]) continue;
    if (eta[i] > 1.24 and eta[i] < 1.6 and phi[i] > 2.44346 and phi[i] < 2.79253) {
      prefiringProbability = hMuonPrefiringNew_hotspot->GetBinContent(1, 3)/(TMath::Exp( (pt[i] - hMuonPrefiringNew_hotspot->GetBinContent(1, 1)) / hMuonPrefiringNew_hotspot->GetBinContent(1, 2) ) + 1);
    } else {
      prefireBin = std::max(1, std::min(hprefire.GetXaxis()->FindFixBin(fabs(eta[i])), nBins));
      prefiringProbability = hprefire.GetBinContent(prefireBin, 3)/(TMath::Exp( (pt[i] - hprefire.GetBinContent(prefireBin, 1)) / hprefire.GetBinContent(prefireBin, 2) ) + 1);
    }
    if (prefiringProbability < 0) prefiringProbability = 0.0;
    else if (prefiringProbability > 1) prefiringProbability = 1.0;
    //if (prefiringProbability < 0 || prefiringProbability > 1) {
      //std::cout << "params 0,1,2 = " << hprefire.GetBinContent(prefireBin, 1) << "," << hprefire.GetBinContent(prefireBin, 2) << ", " << hprefire.GetBinContent(prefireBin, 3) << std::endl;
      //std::cout << "eta,pt = " << eta[i] << "," << pt[i] << "   prefiringProbability = " << prefiringProbability << std::endl;
      //prefiringProbability = 0.0;
    //}
    sf *= (1.0 - prefiringProbability);
  }
  // std::cout << "PREFIRING FOR: " << eraNames[era] << "   nMuons = " << eta.size() << "   sf = " << sf << std::endl;
  return sf;

}

// Vec_d _get_newMuonPrefiringSFvariationSyst(int nSystVar, const Vec_f& eta, const Vec_f& pt, const Vec_f& phi, const Vec_b& looseId, DataEra era = BToF) {
Vec_d _get_newMuonPrefiringSFvariationSyst(const Vec_f& eta, const Vec_f& pt, const Vec_f& phi, const Vec_b& looseId, DataEra era = BToF) {

  // should be called with nSystVar = 3 (at least, can codify more variations) representing nominal, Up, Down
  Vec_d res(3, 1.0); // initialize to 1
  // std::cout << "PREFIRING FOR: " << eraNames[era] << std::endl;
  const TH2D& hprefire = hMuonPrefiringNew[era];
  int nBins = hprefire.GetNbinsX();
  int prefireBin = 0;
  double prefiringProbability = 0.0;

  for (unsigned int i = 0; i < eta.size(); ++i) {
    if (not looseId[i]) continue;
    if (eta[i] > 1.24 and eta[i] < 1.6 and phi[i] > 2.44346 and phi[i] < 2.79253) {
      prefiringProbability = hMuonPrefiringNew_hotspot->GetBinContent(1, 3)/(TMath::Exp( (pt[i] - hMuonPrefiringNew_hotspot->GetBinContent(1, 1)) / hMuonPrefiringNew_hotspot->GetBinContent(1, 2) ) + 1);
    } else {
      prefireBin = std::max(1, std::min(hprefire.GetXaxis()->FindFixBin(fabs(eta[i])), nBins));
      prefiringProbability = hprefire.GetBinContent(prefireBin, 3)/(TMath::Exp( (pt[i] - hprefire.GetBinContent(prefireBin, 1)) / hprefire.GetBinContent(prefireBin, 2) ) + 1);
    }
    if (prefiringProbability < 0) prefiringProbability = 0.0;
    res[0] *= (1.0 - std::min(1.0, prefiringProbability));
    res[1] *= (1.0 - std::min(1.0, 1.11*prefiringProbability));
    res[2] *= (1.0 - std::min(1.0, 0.89*prefiringProbability));
    //if (prefiringProbability < 0 || prefiringProbability > 1) {
      //std::cout << "params 0,1,2 = " << hprefire.GetBinContent(prefireBin, 1) << "," << hprefire.GetBinContent(prefireBin, 2) << ", " << hprefire.GetBinContent(prefireBin, 3) << std::endl;
      //std::cout << "eta,pt = " << eta[i] << "," << pt[i] << "   prefiringProbability = " << prefiringProbability << std::endl;
      //prefiringProbability = 0.0;
    //}
  }
  // std::cout << "PREFIRING FOR: " << eraNames[era] << "   nMuons = " << eta.size() << "   sf = " << sf << std::endl;
  return res;

}

Vec_d _get_newMuonPrefiringSFvariationStat(int n_prefireBinNuisance, const Vec_f& eta, const Vec_f& pt, const Vec_f& phi, const Vec_b& looseId, DataEra era = BToF) {

  // this function directly provides the alternative prefiring SF, not the variation on the original one
  // it is supposed to be called instead of the nominal weight
  // it returns a vector used as an event weight to get all variations in the same TH3 (eta-pt-prefireBin)

  Vec_d res(n_prefireBinNuisance, 1.0); // initialize to 1

  double tmpval = 0.0;
  double prefiringProbability = 0.0;
  double prefiringProbabilityErr = 0.0;
  double prefiringProbabilityUp = 0.0;
  const TH2D& hprefire = hMuonPrefiringNew[era];
  int nBins = hprefire.GetNbinsX();
  int prefireBin = 0;

  for (unsigned int i = 0; i < eta.size(); ++i) {

    if (not looseId[i]) continue;
    prefireBin = std::max(1, std::min(hprefire.GetXaxis()->FindFixBin(fabs(eta[i])), nBins));
    if (eta[i] > 1.24 and eta[i] < 1.6 and phi[i] > 2.44346 and phi[i] < 2.79253) {
      tmpval = 1.0 / (TMath::Exp( (pt[i] - hMuonPrefiringNew_hotspot->GetBinContent(1, 1)) / hMuonPrefiringNew_hotspot->GetBinContent(1, 2) ) + 1);
      prefiringProbability    = hMuonPrefiringNew_hotspot->GetBinContent(1, 3);
      prefiringProbabilityErr = hMuonPrefiringNew_hotspot->GetBinError(1, 3);

    } else {
      tmpval = 1.0 / (TMath::Exp( (pt[i] - hprefire.GetBinContent(prefireBin, 1)) / hprefire.GetBinContent(prefireBin, 2) ) + 1);
      prefiringProbability    = hprefire.GetBinContent(prefireBin, 3);
      prefiringProbabilityErr = hprefire.GetBinError(prefireBin, 3);
    }
    // fill the vector with the nominal weight in each bin,
    // while the vector bin corresponding to prefireBin will be filled with prefiring probability moved by its uncertainty
    prefiringProbabilityUp = tmpval * (prefiringProbability + prefiringProbabilityErr);
    prefiringProbability *= tmpval;
    tmpval = res[prefireBin-1];
    res *= std::min(1.0, std::max(0.0, 1.0 - prefiringProbability));
    res[prefireBin-1] = std::min(1.0, std::max(0.0, tmpval * (1.0 - prefiringProbabilityUp)));

  }

  return res;

}

double _get_fullMuonSF_EtaVsZ(float muonZ,      float eta,      int charge,
			      float muonZOther, float etaOther,
			      DataEra era = BToF,
			      bool altSF = false
			      ) {

  // similar to _get_fullMuonSF but for SF vs eta and Z (PV_z + Muon_dz).
  // only the triggering Muon_dz will be used, the other should have the same dz and in any case dz << PV_z
  // muonZOther is only used to assess whether there are 2 leptons (Z boson) or only one (W bobosn)
  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;

  auto const key = std::make_pair(sftype, era);
  // const TH2F& hsf = useDataAltSig ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  double sf = getValFromTH2(hsf, eta, muonZ);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (muonZOther > 0.0) {
    ScaleFactorType sftypeOther = isoNotrig;
    //ScaleFactorType sftypeOther = noisoNotrig;
    auto const keyOther = std::make_pair(sftypeOther, era);
    // const TH2F& hsfOther = useDataAltSig ? scaleFactorHist_dataAltSig.at(keyOther) : scaleFactorHist.at(keyOther);
    const TH2F& hsfOther = altSF ? scaleFactorHist_dataAltSig.at(keyOther) : scaleFactorHist.at(keyOther);
    sf *= getValFromTH2(hsfOther, etaOther, muonZ);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}


double _get_fullMuonSF(float pt,      float eta,      int charge,
		       float ptOther, float etaOther,
		       DataEra era = BToF,
		       bool isoSF1 = true,
		       bool isoSF2 = true
		       //bool useDataAltSig = false
		       ) {

  // function to get full muon scale factor for  analysis (except prefiring, handled elsewhere)
  // first three arguments are for the triggering muon, second two for the non triggering one
  // isoSF1 and isoSF2 are to use SF for isolation or antiisolation, for triggering and non triggering muons respectively
  // if ptOther < 0, the second lepton is ignored, so we can use this function for Wmass as well (etaOther is not used)
  // may actually use another function for Wmass

  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  // std::cout << "scale factors for " << eraNames[era] << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, era);
  // const TH2F& hsf = useDataAltSig ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  const TH2F& hsf = scaleFactorHist.at(key);
  double sf = getValFromTH2(hsf, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    ScaleFactorType sftypeOther = isoSF2 ? isoNotrig : antiisoNotrig;
    //ScaleFactorType sftypeOther = isoSF2 ? noisoNotrig : antiisoNotrig;
    auto const keyOther = std::make_pair(sftypeOther, era);
    // const TH2F& hsfOther = useDataAltSig ? scaleFactorHist_dataAltSig.at(keyOther) : scaleFactorHist.at(keyOther);
    const TH2F& hsfOther = scaleFactorHist.at(keyOther);
    sf *= getValFromTH2(hsfOther, etaOther, ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

double _get_fullMuonSF_dataAltSig(float pt,      float eta,      int charge,
				  float ptOther, float etaOther,
				  DataEra era = BToF,
				  bool isoSF1 = true,
				  bool isoSF2 = true
				  ) {

  // function to get full muon scale factor for  analysis (except prefiring, handled elsewhere)
  // first three arguments are for the triggering muon, second two for the non triggering one
  // isoSF1 and isoSF2 are to use SF for isolation or antiisolation, for triggering and non triggering muons respectively
  // if ptOther < 0, the second lepton is ignored, so we can use this function for Wmass as well (etaOther is not used)
  // may actually use another function for Wmass

  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  // std::cout << "scale factors for " << eraNames[era] << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = scaleFactorHist_dataAltSig.at(key);
  double sf = getValFromTH2(hsf, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    ScaleFactorType sftypeOther = isoSF2 ? isoNotrig : antiisoNotrig;
    auto const keyOther = std::make_pair(sftypeOther, era);
    const TH2F& hsfOther = scaleFactorHist_dataAltSig.at(keyOther);
    sf *= getValFromTH2(hsfOther, etaOther, ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

double _get_fullMuonSF_dataAltSig_splitIso(bool alternateIsoStep,
                                           float pt,      float eta,      int charge,
                                           float ptOther, float etaOther,
                                           DataEra era = BToF,
                                           bool isoSF1 = true,
                                           bool isoSF2 = true
    ) {

    // function to get full muon scale factor for  analysis (except prefiring, handled elsewhere)
    // first three arguments are for the triggering muon, second two for the non triggering one
    // isoSF1 and isoSF2 are to use SF for isolation or antiisolation, for triggering and non triggering muons respectively
    // if ptOther < 0, the second lepton is ignored, so we can use this function for Wmass as well (etaOther is not used)
    // may actually use another function for Wmass

    //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
    //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
    // std::cout << "scale factors for " << eraNames[era] << std::endl;

    ScaleFactorType sftypeWithoutIsoStep = charge > 0 ? noisoTrigPlus : noisoTrigMinus; // reco*tracking not included here
    ScaleFactorType sftypeOnlyIsoStep = isoOnly;
    if (not isoSF1)
        sftypeOnlyIsoStep = antiisoOnly;

    //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
    //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
    //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

    auto const keyWithoutIsoStep = std::make_pair(sftypeWithoutIsoStep, era);
    const TH2F& hsfWithoutIsoStep = alternateIsoStep ? scaleFactorHist.at(keyWithoutIsoStep) : scaleFactorHist_dataAltSig.at(keyWithoutIsoStep); // if alternateIsoStep use nominal SF for the other pieces

    auto const keyOnlyIsoStep = std::make_pair(sftypeOnlyIsoStep, era);
    const TH2F& hsfOnlyIsoStep = alternateIsoStep ? scaleFactorHist_dataAltSig.at(keyOnlyIsoStep) : scaleFactorHist.at(keyOnlyIsoStep); // if alternateIsoStep use alternate SF for the isolation

    double sf = getValFromTH2(hsfWithoutIsoStep, eta, pt) * getValFromTH2(hsfOnlyIsoStep, eta, pt);
    //std::cout << "scale factor main leg -> " << sf << std::endl;

    if (ptOther > 0.0) {

        ScaleFactorType sftypeOtherWithoutIsoStep = noisoNotrig; // basically only idip (reco*tracking are applied elsewhere)
        ScaleFactorType sftypeOtherOnlyIsoStep = isoNotrigOnly;
        if (not isoSF2)
            sftypeOtherOnlyIsoStep = antiisoNotrigOnly;

        auto const keyOtherWithoutIsoStep = std::make_pair(sftypeOtherWithoutIsoStep, era);
        const TH2F& hsfOtherWithoutIsoStep = alternateIsoStep ? scaleFactorHist.at(keyOtherWithoutIsoStep) : scaleFactorHist_dataAltSig.at(keyOtherWithoutIsoStep);

        auto const keyOtherOnlyIsoStep = std::make_pair(sftypeOtherOnlyIsoStep, era);
        const TH2F& hsfOtherOnlyIsoStep = alternateIsoStep ? scaleFactorHist_dataAltSig.at(keyOtherOnlyIsoStep) : scaleFactorHist.at(keyOtherOnlyIsoStep);

        sf *= (getValFromTH2(hsfOtherWithoutIsoStep, etaOther, ptOther) * getValFromTH2(hsfOtherOnlyIsoStep, etaOther, ptOther));

    }
    //std::cout << "final scale factor -> " << sf << std::endl;
    return sf;

}


// tnp reco scale factors, not included in the products used above
double _get_tnpTrackingSF(float pt,      float eta,      int charge,
			  float ptOther, float etaOther,
			  DataEra era = BToF,
			  bool altSF = false,
			  ScaleFactorType sftype = tracking
			  ) {

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  double sf = getValFromTH2(hsf, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    sf *= getValFromTH2(hsf, etaOther, ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

// tnp reco scale factors, not included in the products used above
double _get_tnpTrackingSF_EtaVsZ(float muonZ,      float eta,      int charge,
				 float muonZOther, float etaOther,
				 DataEra era = BToF,
				 bool altSF = false,
				 ScaleFactorType sftype = tracking
				 ) {

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  double sf = getValFromTH2(hsf, eta, muonZ);
  double tmp = 0.0;
  // temporary patch for some unstable bins
  if (sf < 0.99 or sf > 1.01) sf = 1.0;
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (muonZOther > 0.0) {
    tmp = getValFromTH2(hsf, etaOther, muonZ);
    if (tmp > 0.99 and tmp < 1.01) sf *=  tmp;
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}


// tnp reco scale factors, not included in the products used above
double _get_tnpRecoSF(float pt,      float eta,      int charge,
		      float ptOther, float etaOther,
		      DataEra era = BToF,
		      bool altSF = false,
		      ScaleFactorType sftype = reco
		      ) {

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  double sf = getValFromTH2(hsf, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    sf *= getValFromTH2(hsf, etaOther, ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

// tnp reco scale factors vs eta-Z, not included in the products used above
double _get_tnpRecoSF_EtaVsZ(float muonZ,      float eta,      int charge,
			     float muonZOther, float etaOther,
			     DataEra era = BToF,
			     bool altSF = false,
			     ScaleFactorType sftype = reco
			     ) {

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
  // temporary path for 1 bin, use value for alternative fit to data rather than nominal
  double sf = (era == BToF and eta < -2.3 and muonZ < -7.5) ? 0.97740746 : getValFromTH2(hsf, eta, muonZ);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (muonZOther > 0.0) {
    if (era == BToF and etaOther < -2.3 and muonZ < -7.5)
      sf *= 0.97740746;
    else
      sf *= getValFromTH2(hsf, etaOther, muonZ);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}


double _get_tnpTrackingRecoSF(float pt,      float eta,      int charge,
                              float ptOther, float etaOther,
                              DataEra era = BToF,
                              bool altSF = false
    ) {

    auto const key = std::make_pair(trackingReco, era);
    const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
    double sf = getValFromTH2(hsf, eta, pt);
    //std::cout << "scale factor main leg -> " << sf << std::endl;

    if (ptOther > 0.0) {
        sf *= getValFromTH2(hsf, etaOther, ptOther);
    }
    //std::cout << "final scale factor -> " << sf << std::endl;
    return sf;

}

double _get_tnpTrackingRecoSFvariation(int n_tnpBinNuisance,
                                       float pt,      float eta,      int charge,
                                       float ptOther, float etaOther,
                                       DataEra era = BToF,
                                       bool altSF = false
    ) {

    auto const key = std::make_pair(trackingReco, era);
    const TH2F& hsf = altSF ? scaleFactorHist_dataAltSig.at(key) : scaleFactorHist.at(key);
    //std::cout << "scale factor main leg -> " << sf << std::endl;

    int nEtaBins = hsf.GetNbinsX();
    int nPtBins  = hsf.GetNbinsY();
    int ietaTnP = std::min(nEtaBins, std::max(1, hsf.GetXaxis()->FindFixBin(eta)));
    int iptTnP  = std::min(nPtBins,  std::max(1, hsf.GetYaxis()->FindFixBin(pt)));
    int tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
    double sf = hsf.GetBinContent(ietaTnP, iptTnP);
    // initialize to nominal SF
    Vec_d res(n_tnpBinNuisance, sf);
    if (tnpBinNuisance <= n_tnpBinNuisance) {
        res[tnpBinNuisance - 1] += hsf.GetBinError(ietaTnP, iptTnP);
    }

    if (ptOther > 0.0) {
        ietaTnP = std::min(nEtaBins, std::max(1, hsf.GetXaxis()->FindFixBin(etaOther)));
        iptTnP = std::min(nPtBins, std::max(1, hsf.GetYaxis()->FindFixBin(ptOther)));
        tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
        double tmp = res[tnpBinNuisance-1]; // save value of sf for the first lepton only, it is needed below if using variations
        sf = hsf.GetBinContent(ietaTnP, iptTnP);
        res *= sf; // multiply all cells by sf for second lepton
        if (tnpBinNuisance <= n_tnpBinNuisance) {
            res[tnpBinNuisance-1] = tmp * (sf + hsf.GetBinError(ietaTnP, iptTnP));
        }
    }
    //std::cout << "final scale factor -> " << sf << std::endl;
    return sf;

}


Vec_d _get_fullMuonSFvariation(int n_tnpBinNuisance,
                               float pt,      float eta, int charge,
                               float ptOther=-1, float etaOther=-1,
                               DataEra era = BToF,
                               bool isoSF1 = true,
                               bool isoSF2 = true)
{

  // this is an helper function to define the effSystvariations for the Wlike analysis
  // idea is to fill again the alternative template for each effStat nuisance, where the
  // nuisance is defined for each single TnP eta-pt bin, and they will be considered as uncorrelated
  // so not as done in SMP-18-012 from the fit to efficiencies with Error function
  //
  // this should replace the nominal SF weight, and return SF for any bin except the specific one corresponding
  // to the nuisance parameter, and SF+err for that one
  // in order to facilitate the usage with the ReplaceWeight functionality to customize event weight per histogram,
  // the input arguments should be the same as in _get_fullMuonSF, and new ones should appear in the beginning

  // n_tnpBinNuisance is the number of TnP bins, used to set the size of the RVec that will be returned

  // tnpBinNuisance is supposed to start from 1 and be mapped into SF histogram bins as shown below
  //
  // pt | 7 | 8 | 9 |
  //     --- --- ---
  //    | 4 | 5 | 6 |
  //     --- --- ---
  //    | 1 | 2 | 3 |
  //              eta

  // if ptOther < 0 it is assumed only one lepton exists (so this function could also be used for wmass)
  // in that case the values are not used (and etaOther could actually take any value)

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  auto const key = std::make_pair(sftype, era);
  const TH2F& hsf = scaleFactorHist.at(key);
  int nEtaBins = hsf.GetNbinsX();
  int nPtBins  = hsf.GetNbinsY();

  int ietaTnP = std::min(nEtaBins, std::max(1, hsf.GetXaxis()->FindFixBin(eta)));
  int iptTnP  = std::min(nPtBins,  std::max(1, hsf.GetYaxis()->FindFixBin(pt)));
  int tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
  // watch out, tnpBinNuisance  must not get larger than n_tnpBinNuisance. If you want to use less tnp bins than the histogram actually has, you need a pt cut accordingly
  // in case it is larger, we just do not change it, so to avoid reading res[i] at the i-th cell, which would not be within the vector size

   // initialize to nominal SF
  Vec_d res(n_tnpBinNuisance, hsf.GetBinContent(ietaTnP, iptTnP));
  // sum or subtract error in specific bin
  // for isolation, one has to account for anticorrelation between isolation and anti-isolation efficiency
  // here we act on the scale factors, but it should be a reasonable approximation anyway
  if (tnpBinNuisance <= n_tnpBinNuisance) {
    if (isoSF1)
      res[tnpBinNuisance-1] += hsf.GetBinError(ietaTnP, iptTnP);
    else
      res[tnpBinNuisance-1] -= hsf.GetBinError(ietaTnP, iptTnP);
  }

  if (ptOther > 0) {
    ScaleFactorType sftypeOther = isoNotrig;
    if (not isoSF2)
      sftypeOther = antiisoNotrig;
    auto const keyOther = std::make_pair(sftypeOther, era);
    const TH2F& hsfOther = scaleFactorHist.at(keyOther);
    ietaTnP = std::min(nEtaBins, std::max(1, hsfOther.GetXaxis()->FindFixBin(etaOther)));
    iptTnP  = std::min(nPtBins,  std::max(1, hsfOther.GetYaxis()->FindFixBin(ptOther)));
    tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
    double tmp = res[tnpBinNuisance-1]; // save value of sf for the first lepton only, it is needed below if using variations
    double sfOther = hsfOther.GetBinContent(ietaTnP, iptTnP);
    res *= sfOther; // multiply all cells by sf for second lepton
    if (tnpBinNuisance <= n_tnpBinNuisance) {
      // now modify one cell to account for the sf variation for the second lepton
      // see comment above for isolation part
      if (isoSF2)
          res[tnpBinNuisance-1] = tmp * (sfOther + hsfOther.GetBinError(ietaTnP, iptTnP));
      else
          res[tnpBinNuisance-1] = tmp * (sfOther - hsfOther.GetBinError(ietaTnP, iptTnP));
    }
  }

  return res;
}


Vec_d _get_fullMuonSFvariation_splitIso(int n_tnpBinNuisance,
                                        float pt,      float eta, int charge,
                                        float ptOther=-1, float etaOther=-1,
                                        DataEra era = BToF,
                                        bool isoSF1 = true,
                                        bool isoSF2 = true)
{

    // this is an helper function to define the effSystvariations for the Wlike analysis
    // idea is to fill again the alternative template for each effStat nuisance, where the
    // nuisance is defined for each single TnP eta-pt bin, and they will be considered as uncorrelated
    // so not as done in SMP-18-012 from the fit to efficiencies with Error function
    //
    // this should replace the nominal SF weight, and return SF for any bin except the specific one corresponding
    // to the nuisance parameter, and SF+err for that one
    // in order to facilitate the usage with the ReplaceWeight functionality to customize event weight per histogram,
    // the input arguments should be the same as in _get_fullMuonSF, and new ones should appear in the beginning

    // n_tnpBinNuisance is the number of TnP bins, used to set the size of the RVec that will be returned

    // tnpBinNuisance is supposed to start from 1 and be mapped into SF histogram bins as shown below
    //
    // pt | 7 | 8 | 9 |
    //     --- --- ---
    //    | 4 | 5 | 6 |
    //     --- --- ---
    //    | 1 | 2 | 3 |
    //              eta

    // if ptOther < 0 it is assumed only one lepton exists (so this function could also be used for wmass)
    // in that case the values are not used (and etaOther could actually take any value)

    // isolation is considered separately wrt other sf, because that uncertainty will be anticorrelated between isoalte and antiisolated regions (latter are for the fakes)
    // so effectively we return M variations where M = 2 * N, N being the number of tnp bins

    ScaleFactorType sftypeWithoutIsoStep = charge > 0 ? noisoTrigPlus : noisoTrigMinus; // reco*tracking not included here
    ScaleFactorType sftypeOnlyIsoStep = isoOnly;
    if (not isoSF1)
        sftypeOnlyIsoStep = antiisoOnly;

    //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
    //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
    //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

    auto const keyWithoutIsoStep = std::make_pair(sftypeWithoutIsoStep, era);
    const TH2F& hsfWithoutIsoStep = scaleFactorHist.at(keyWithoutIsoStep);

    auto const keyOnlyIsoStep = std::make_pair(sftypeOnlyIsoStep, era);
    const TH2F& hsfOnlyIsoStep = scaleFactorHist.at(keyOnlyIsoStep);

    // binning is the same for all SF histograms, can get number of bins with this one
    int nEtaBins = hsfWithoutIsoStep.GetNbinsX();
    int nPtBins  = hsfWithoutIsoStep.GetNbinsY();
    int ietaTnP = std::min(nEtaBins, std::max(1, hsfWithoutIsoStep.GetXaxis()->FindFixBin(eta)));
    int iptTnP  = std::min(nPtBins,  std::max(1, hsfWithoutIsoStep.GetYaxis()->FindFixBin(pt)));
    int tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
    // watch out, tnpBinNuisance  must not get larger than n_tnpBinNuisance. If you want to use less tnp bins than the histogram actually has, you need a pt cut accordingly
    // in case it is larger, we just do not change it, so to avoid reading res[i] at the i-th cell, which would not be within the vector size

    // initialize to nominal SF, using 2 * N_tnpBins
    double sfWithoutIsoStep = hsfWithoutIsoStep.GetBinContent(ietaTnP, iptTnP);
    double sfOnlyIsoStep    = hsfOnlyIsoStep.GetBinContent(ietaTnP, iptTnP);
    Vec_d res(2 * n_tnpBinNuisance, sfWithoutIsoStep * sfOnlyIsoStep);
    // sum or subtract error in specific bin
    // for isolation, one has to account for anticorrelation between isolation and anti-isolation efficiency
    // here we act on the scale factors rather than the efficiency, but it should be a reasonable approximation anyway, with our numbers we verified also the SF are almost 100% anticorrelated
    if (tnpBinNuisance <= n_tnpBinNuisance) {
        // first fill element for part without isolation
        res[tnpBinNuisance-1] = (sfWithoutIsoStep + hsfWithoutIsoStep.GetBinError(ietaTnP, iptTnP)) * sfOnlyIsoStep;
        // and then element in the second set of N elements, with only isolation part
        if (isoSF1)
            res[n_tnpBinNuisance + tnpBinNuisance-1] = sfWithoutIsoStep * (sfOnlyIsoStep + hsfOnlyIsoStep.GetBinError(ietaTnP, iptTnP));
        else
            res[n_tnpBinNuisance + tnpBinNuisance-1] = sfWithoutIsoStep * std::max(0.0, sfOnlyIsoStep - hsfOnlyIsoStep.GetBinError(ietaTnP, iptTnP));
    }

    if (ptOther > 0) {

        ScaleFactorType sftypeOtherWithoutIsoStep = noisoNotrig; // basically only idip (reco*tracking are applied elsewhere)
        ScaleFactorType sftypeOtherOnlyIsoStep = isoNotrigOnly;
        if (not isoSF2)
            sftypeOtherOnlyIsoStep = antiisoNotrigOnly;

        auto const keyOtherWithoutIsoStep = std::make_pair(sftypeOtherWithoutIsoStep, era);
        const TH2F& hsfOtherWithoutIsoStep = scaleFactorHist.at(keyOtherWithoutIsoStep);

        auto const keyOtherOnlyIsoStep = std::make_pair(sftypeOtherOnlyIsoStep, era);
        const TH2F& hsfOtherOnlyIsoStep = scaleFactorHist.at(keyOtherOnlyIsoStep);

        ietaTnP = std::min(nEtaBins, std::max(1, hsfOtherWithoutIsoStep.GetXaxis()->FindFixBin(etaOther)));
        iptTnP  = std::min(nPtBins,  std::max(1, hsfOtherWithoutIsoStep.GetYaxis()->FindFixBin(ptOther)));

        tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
        double tmpWithoutIsoStep = res[tnpBinNuisance-1]; // save value of sf for the first lepton only, it is needed below if using variations
        double tmpOnlyIsoStep    = res[n_tnpBinNuisance + tnpBinNuisance-1]; // save value of sf for the first lepton only, it is needed below if using variations

        sfWithoutIsoStep = hsfOtherWithoutIsoStep.GetBinContent(ietaTnP, iptTnP);
        sfOnlyIsoStep    = hsfOtherOnlyIsoStep.GetBinContent(ietaTnP, iptTnP);
        res *= (sfWithoutIsoStep * sfOnlyIsoStep); // multiply all cells by sf for second lepton

        if (tnpBinNuisance <= n_tnpBinNuisance) {
            // first fill element for part without isolation
            res[tnpBinNuisance-1] = tmpWithoutIsoStep * (sfWithoutIsoStep + hsfOtherWithoutIsoStep.GetBinError(ietaTnP, iptTnP)) * sfOnlyIsoStep;
            // and then element in the second set of N elements, with only isolation part
            // see comment above for isolation part
            if (isoSF2)
                res[n_tnpBinNuisance + tnpBinNuisance-1] = tmpOnlyIsoStep * sfWithoutIsoStep * (sfOnlyIsoStep + hsfOtherOnlyIsoStep.GetBinError(ietaTnP, iptTnP));
            else
                res[n_tnpBinNuisance + tnpBinNuisance-1] = tmpOnlyIsoStep * sfWithoutIsoStep * std::max(0.0, sfOnlyIsoStep - hsfOtherOnlyIsoStep.GetBinError(ietaTnP, iptTnP));
        }
    }

    return res;

}


double qcdScaleWeight_VptBinned(const double& qcdscale, const double& vpt, const double& ptlow, const double& pthigh) {

  if (vpt >= ptlow and vpt < pthigh)
    return qcdscale;
  else
    return 1.0;

}

Vec_f qcdScaleWeight_VptBinned(const Vec_f& qcdscale, const double& vpt, const double& ptlow, const double& pthigh) {

  if (vpt >= ptlow and vpt < pthigh) {
    return qcdscale;
  } else {
    Vec_f res(qcdscale.size(),1.0); // initialize to 1
    return res;
  }

}


bool isInAccEtaPt(Float_t eta, Float_t pt,
		  Float_t etalow, Float_t etahigh,
		  Float_t ptlow, Float_t pthigh) {

  if (eta > etalow and eta < etahigh and pt > ptlow and pt < pthigh)
    return 1;
  else
    return 0;

}

bool isOutAccEtaPt(Float_t eta, Float_t pt,
		   Float_t etalow, Float_t etahigh,
		   Float_t ptlow, Float_t pthigh) {

  if (isInAccEtaPt(eta,pt,etalow,etahigh,ptlow,pthigh))
    return 0;
  else
    return 1;

}


float safeRatio(float num, float den, float safe = 1.0) {
  return (den != 0.0) ? (num/den) : safe;
}

int regionIsoMt(bool lowIso, bool lowMt) {

  if      (not lowIso and     lowMt) return 0; // fakes region (failing isolation)
  else if (    lowIso and     lowMt) return 1; // fakes region (passing isolation)
  else if (not lowIso and not lowMt) return 2; // fakes application region
  else if (    lowIso and not lowMt) return 3; // signal region
  return -1;  // should be impossible to get here, but just in case

}

Vec_f shiftPt_testPtScaleSystWmass(const float& pt, const float& eta, const unsigned int& nEtaBins=48, const double & binSize = 0.1, const double& etaMin = -2.4, const double& ptMin = -1.0, const double& ptMax = -1.0) {

  // return an array of size 2*nEtaBins, which contains pt in all bins except for the one corresponding to the value of eta
  // this embeds both positive and negative pt shifts in the same array
  // one can optionally specify ptMin and ptMax (setting both of them to a positive value) to make sure that the returned value is still in an allowed range
  Vec_f corPt(2*nEtaBins, pt); // initialize to nominal pT
  double shift = 0.0001;
  if      (fabs(eta) > 2.0) shift = 0.0003;
  else if (fabs(eta) > 1.2) shift = 0.0002;

  unsigned int binID = 0;
  for (unsigned int i = 1; i <= nEtaBins; ++i) {
    if (eta < (etaMin + binSize * i)) {
      binID = i - 1;
      if (ptMin > 0 and ptMax > 0) {
	corPt[binID]            = std::min(ptMax, std::max(ptMin, pt * (1 + shift) ) );
	corPt[nEtaBins + binID] = std::min(ptMax, std::max(ptMin, pt * (1 - shift) ) );
      } else {
	corPt[binID]            = pt * (1 + shift);
	corPt[nEtaBins + binID] = pt * (1 - shift);
      }
      break;
    }
  }
  return corPt;

}

}

#endif
