#ifndef WREMNANTS_FUNCTIONS_H
#define WREMNANTS_FUNCTIONS_H


// #include <stdio.h>
// #include <stdlib.h>
#include <iostream>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TRandom3.h"

#include "defines.h"

namespace wrem {

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

//// UTILITY FUNCTIONS NOT IN TFORMULA ALREADY

// ROOT::RDF::RResultPtr<double> getRDFcolumnSum(const ROOT::RDataFrame& rdf, const std::string& column) {
//   return rdf::Sum<double>(column);
// }

template <typename T>
T cropLargeValue(const T& wgt, const double max) {
  // max is already a positive number, no need to check or take absolute value here
  return static_cast<T>(std::copysign(1.0, wgt) * std::min<T>(std::abs(wgt), max));
}

template <typename T>
T removeLargeValue(const T& wgt, const T& max, const T& valueToZero = 0.0) {
    // max is already a positive number, no need to check or take absolute value here
    // std::cout << typeid(T).name() << std::endl;
    // std::cout << wgt << " " << max << " " << valueToZero << std::endl;
    return (std::abs(wgt) < max) ? wgt : valueToZero;
}


double genWeightLargeClipped(const double& wgt, const double& max) {

  // max is already a positive number, no need to check or take absolute value here
  return static_cast<double>(std::copysign(1.0, wgt) * std::min<double>(std::abs(wgt), max));

}

double genWeightLargeRemoved(const double& wgt, const double& max) {

  // max is already a positive number, no need to check or take absolute value here
  return (std::abs(wgt) < max) ? wgt : 0.0;

}

template <typename T>
Vec_i indices(const ROOT::VecOps::RVec<T>& vec, const int start = 0) {
    Vec_i res(vec.size(), 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

// Can return a std::array
Vec_i indices(const size_t size, const int start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

// Thanks stack overflow https://stackoverflow.com/questions/42749032/concatenating-a-sequence-of-stdarrays
template<typename T, int N, int M>
auto concatRVecsToArray(const ROOT::VecOps::RVec<T>& vec1, ROOT::VecOps::RVec<T>& vec2)
{
    std::array<T, N+M> result;
    std::copy (vec1.cbegin(), vec1.cend(), result.begin());
    std::copy (vec2.cbegin(), vec2.cend(), result.begin() + N);
    return result;
}

template<typename T>
auto concatRVecs(const ROOT::VecOps::RVec<T>& vec1, ROOT::VecOps::RVec<T>& vec2)
{
    ROOT::VecOps::RVec<T> result(vec1.size()+vec2.size());
    std::copy (vec1.cbegin(), vec1.cend(), result.begin());
    std::copy (vec2.cbegin(), vec2.cend(), result.begin()+vec1.size());
    return result;
}

TRandom3 *rand_smear = new TRandom3(0);

Vec_f shiftVar(const Vec_f& var, float shift = 0.0) {

  Vec_f res(var.size(),0);
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = var[i] * (1 + shift);
  }
  return res;

}

Vec_f smearAndShiftVar(const Vec_f& var, float shift = 0.0, float smear = 0.0) {

  Vec_f res(var.size(),0);
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = var[i] * (1 + shift + smear * rand_smear->Gaus(0.,1.));
  }
  return res;

}

double deltaPhi(float phi1, float phi2) {
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2.0*M_PI;
    while (result <= -1.0*M_PI) result += 2.0*M_PI;
    return result;
}

float if3(bool cond, float iftrue, float iffalse) {
    return cond ? iftrue : iffalse;
}

double deltaR2(float eta1, float phi1, float eta2, float phi2) {
    double deta = eta1-eta2;
    double dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}

double deltaR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

double minDeltaR2many(const Vec_f& eta, const Vec_f& phi, float etaProbe, float phiProbe) {

  double retmin = std::numeric_limits<double>::max();
  if (eta.size() == 0) return retmin;
  double dr2 = 0.0;
  for (unsigned int ij = 0; ij < eta.size(); ++ij) {
    dr2 = deltaR2(eta[ij], phi[ij], etaProbe, phiProbe);
    if (dr2 < retmin) retmin = dr2;
  }
  return retmin;
}

// returns minimum DR^2 between two vector collections
double minDeltaR2manyMore(const Vec_f& eta, const Vec_f& phi, const Vec_f& etaProbe, const Vec_f& phiProbe) {

  double retmin = std::numeric_limits<double>::max();
  if (eta.size() == 0 or etaProbe.size() == 0 ) return retmin;
  double dr2 = 0.0;
  for (unsigned int ij = 0; ij < eta.size(); ++ij) {
    for (unsigned int ik = 0; ik < etaProbe.size(); ++ik) {
      dr2 = deltaR2(eta[ij], phi[ij], etaProbe[ik], phiProbe[ik]);
      if (dr2 < retmin) retmin = dr2;
    }
  }
  return retmin;
}

// returns vector of bools filtering input collection based on DR^2 match with respect to a test collection
Vec_b matchDeltaR2(const Vec_f& eta, const Vec_f& phi, const Vec_f& etaTest, const Vec_f& phiTest, double dr2match = 0.01) {

  Vec_b res(eta.size(), false); // initialize to 0 (no match)
  if (eta.size() == 0 or etaTest.size() == 0 ) return res;
  double dr2 = 0.0;
  for (unsigned int ij = 0; ij < eta.size(); ++ij) {
    double retmin = std::numeric_limits<double>::max();
    for (unsigned int ik = 0; ik < etaTest.size(); ++ik) {
      dr2 = deltaR2(eta[ij], phi[ij], etaTest[ik], phiTest[ik]);
      if (dr2 < retmin) retmin = dr2;
    }
    if (retmin < dr2match) res[ij] = true;
  }
  return res;
}


float Hypot(float x, float y) {
  return hypot(x,y);
}

float pt_2(float pt1, float phi1, float pt2, float phi2) {
    phi2 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2), pt2*std::sin(phi2));
}

float rapidity_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Rapidity();
}

float mt_2(float pt1, float phi1, float pt2, float phi2) {
    return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

float mass_2_ene(float ene1, float eta1, float phi1, float m1, float ene2, float eta2, float phi2, float m2) {
  //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector unitp41(1.0,eta1,phi1,m1);
    PtEtaPhiMVector unitp42(1.0,eta2,phi2,m2);
    double theta1 = unitp41.Theta();
    double theta2 = unitp42.Theta();
    double pt1 = ene1*fabs(sin(theta1));
    double pt2 = ene2*fabs(sin(theta2));
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float mass_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float invariantmass(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).mass();
}

// Vec_f invariantmasses(float pt1, float eta1, float phi1, float m1, const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
float invMLepPhotons(const Vec_f& pt1, const Vec_f& eta1, const Vec_f& phi1, float m1, const Vec_f& pt, const Vec_f& eta, const Vec_f& phi) {
  // PtEtaPhiMVector p1(pt1[0], eta1[0], phi1[0], m1);
  PtEtaPhiMVector psum(0,0,0,0);
  for (unsigned int i = 0; i < pt1.size(); ++i) {
    PtEtaPhiMVector p1(pt1[i], eta1[i], phi1[i], m1);
    psum += p1;
  }
  for (unsigned int i = 0; i < pt.size(); ++i) {
    PtEtaPhiMVector p2(pt[i], eta[i], phi[i], 0.);
    psum += p2;
  }
  return psum.mass();
}

float rapidity(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).Rapidity();
}

float phi_pair(const Vec_f& pt, const Vec_f& phi) {
    TVector2 p1(pt[0]*std::cos(phi[0]), pt[0]*std::sin(phi[0]));
    TVector2 p2(pt[1]*std::cos(phi[1]), pt[1]*std::sin(phi[1]));
    return TVector2::Phi_mpi_pi((p1+p2).Phi()); // in the ntuples phi is defined between -pi and pi, but TVector2::Phi returns in 0,2pi
}


float transversemomentum(const Vec_f& pt, const Vec_f& phi) {
  float phidiff = phi[1] - phi[0];
  return hypot(pt[0] + pt[1] * std::cos(phidiff), pt[1]*std::sin(phidiff));
}

Vec_b chargedParticleByEventParity(const ULong64_t& event, const Vec_b& plus, const Vec_b& minus) {

  if (isOddEvent(event)) return plus;
  else                   return minus;

}

Vec_b goodMuonTriggerCandidate(const Vec_i& TrigObj_id, const Vec_f& TrigObj_pt, const Vec_f& TrigObj_l1pt, const Vec_f& TrigObj_l2pt, const Vec_i& TrigObj_filterBits) {

   Vec_b res(TrigObj_id.size(),false); // initialize to 0
   for (unsigned int i = 0; i < res.size(); ++i) {
       if (TrigObj_id[i]  != 13 ) continue;
       if (TrigObj_pt[i]   < 24.) continue;
       if (TrigObj_l1pt[i] < 22.) continue;
       if (! (( TrigObj_filterBits[i] & 8) || (TrigObj_l2pt[i] > 10. && (TrigObj_filterBits[i] & 2) )) ) continue;
       res[i] = true;
   }
   // res will be goodTrigObjs in RDF
   // e.g. RDF::Define("goodTrigObjs","goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
   return res;
}

Vec_b hasTriggerMatch(const Vec_f& eta, const Vec_f& phi, const Vec_f& TrigObj_eta, const Vec_f& TrigObj_phi) {

   Vec_b res(eta.size(),false); // initialize to 0
   for (unsigned int i = 0; i < res.size(); ++i) {
      for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
	  // use deltaR*deltaR < 0.3*0.3, to be faster
          if (deltaR2(eta[i], phi[i], TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.09) {
              res[i] = true;
              break; // exit loop on trigger objects, and go to next muon
          }
      }
   }
   // res will be triggerMatchedMuons in RDF, like
   // RDF::Define("triggerMatchedMuons","hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
   return res;

}

bool hasTriggerMatch(const float& eta, const float& phi, const Vec_f& TrigObj_eta, const Vec_f& TrigObj_phi) {

  for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
    if (deltaR2(eta, phi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.09) return true;
  }
  return false;

}

Vec_b cleanJetsFromMuons(const Vec_f& Jet_eta, const Vec_f& Jet_phi, const Vec_f& Muon_eta, const Vec_f& Muon_phi) {

   Vec_b res(Jet_eta.size(), true); // initialize to true and set to false whenever the jet overlaps with a muon
   for (unsigned int ij = 0; ij < res.size(); ++ij) {
     for (unsigned int im = 0; im < Muon_eta.size(); ++im) {
       if (deltaR2(Jet_eta[ij], Jet_phi[ij], Muon_eta[im], Muon_phi[im]) < 0.16) { // cone DR = 0.4
	 res[ij] = false;
	 break;
       }
     }
   }

   return res;
}

Vec_b cleanJetsFromLeptons(const Vec_f& Jet_eta, const Vec_f& Jet_phi, const Vec_f& Muon_eta, const Vec_f& Muon_phi, const Vec_f& Electron_eta, const Vec_f& Electron_phi) {

   Vec_b res(Jet_eta.size(), true); // initialize to true and set to false whenever the jet overlaps with a muon

   for (unsigned int ij = 0; ij < res.size(); ++ij) {

     for (unsigned int im = 0; im < Muon_eta.size(); ++im) {
       if (deltaR2(Jet_eta[ij], Jet_phi[ij], Muon_eta[im], Muon_phi[im]) < 0.16) { // cone DR = 0.4
	 res[ij] = false;
	 break;
       }
     }

     if (res[ij]) {
       for (unsigned int ie = 0; ie < Electron_eta.size(); ++ie) {
	 if (deltaR2(Jet_eta[ij], Jet_phi[ij], Electron_eta[ie], Electron_phi[ie]) < 0.16) { // cone DR = 0.4
	   res[ij] = false;
	   break;
	 }
       }
     }

   }

   return res;
}

template <typename T>
ROOT::VecOps::RVec<T> absoluteValue(const ROOT::VecOps::RVec<T> & val) {

  ROOT::VecOps::RVec<T> res(val.size(), 0.0); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = std::abs(val[i]);
  }
  return res;

}

Vec_f absoluteValue(const Vec_f& val) {

  Vec_f res(val.size(),0.0); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = std::abs(val[i]);
  }
  return res;

}

Vec_i absoluteValue(const Vec_i& val) {

  Vec_i res(val.size(),0.0); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = std::abs(val[i]);
  }
  return res;

}

template <typename T>
ROOT::VecOps::RVec<T> truncateAndCropVector(const ROOT::VecOps::RVec<T>& vec, double maxVal, unsigned int numberOfCells = -1, unsigned int minIndex = 0) {

    if (numberOfCells < 0 or numberOfCells > vec.size())
        numberOfCells = vec.size();
    ROOT::VecOps::RVec<T> res(std::begin(vec)+minIndex, std::begin(vec)+numberOfCells);
    for (unsigned int i = 0; i < res.size(); ++i) {
        res[i] = cropLargeValue(res[i], maxVal);
    }
    return res;

}



TRandom3 *randy = NULL;

float mass_2_smeared(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, int isData) {
    float finalpt1;
    float finalpt2;
    if (isData){
        finalpt1 = pt1;
        finalpt2 = pt2;
    }
    else {
        if (!randy) randy = new TRandom3(42);
        finalpt1 = pt1*(1.+0.34/pt1/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt1;
        finalpt2 = pt2*(1.+0.34/pt2/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt2;
    }

//    std::cout << "initial pT1: " << pt1 << " corrected pT1: " << finalpt1 << std::endl;
//    std::cout << "initial pT2: " << pt2 << " corrected pT2: " << finalpt2 << std::endl;

    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(finalpt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(finalpt2,eta2,phi2,m2);

    PtEtaPhiMVector p43(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p44(pt2,eta2,phi2,m2);

    float finalmll = (p41+p42).M();
    float initialm = (p43+p44).M();

    //std::cout << "initial mll " << initialm << " final mll " << finalmll << std::endl;


    return (p41+p42).M();
}

float eta_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Eta();
}

float pt_3(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3) {
    phi2 -= phi1;
    phi3 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3), pt2*std::sin(phi2) + pt3*std::sin(phi3));
}

float mass_3(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    return (p41+p42+p43).M();
}

float pt_4(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3, float pt4, float phi4) {
    phi2 -= phi1;
    phi3 -= phi1;
    phi4 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3) + pt4 * std::cos(phi4), pt2*std::sin(phi2) + pt3*std::sin(phi3) + pt4*std::sin(phi4));
}

float mass_4(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float pt4, float eta4, float phi4, float m4) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    PtEtaPhiMVector p44(pt4,eta4,phi4,m4);
    return (p41+p42+p43+p44).M();
}

float mt_llv(float ptl1, float phil1, float ptl2, float phil2, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}

float mt_lllv(float ptl1, float phil1, float ptl2, float phil2, float ptl3, float phil3, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptl3*std::cos(phil3) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptl3*std::sin(phil3) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptl3+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}


float mtw_wz3l(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float mZ1, float met, float metphi)
{
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)) < 0.01) return mt_2(pt3,phi3,met,metphi);
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt2,phi2,met,metphi);
    if (abs(mZ1 - mass_2(pt2,eta2,phi2,m2,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt1,phi1,met,metphi);
    return 0;
}

float mt_lu_cart(float lep_pt, float lep_phi, float u_x, float u_y)
{
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float u = hypot(u_x,u_y);
    float uDotLep = u_x*lep_px + u_y*lep_py;
    return sqrt(2*lep_pt*sqrt(u*u+lep_pt*lep_pt+2*uDotLep) + 2*uDotLep + 2*lep_pt*lep_pt);
}

float u1_2(float met_pt, float met_phi, float ref_pt, float ref_phi)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_px + uy*ref_py)/ref_pt;
}
float u2_2(float met_pt, float met_phi, float ref_pt, float ref_phi)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_py - uy*ref_px)/ref_pt;
}


double parallelProjection(float phi, float ptRef, float phiRef) {
    // return u = a*ref/|a| where * is the scalar product
    // the magnitude of a is not necessary
    return ptRef * std::cos(phiRef - phi);
}

double parallelProjectionLepBoson(const Vec_f& pt, const Vec_f& phi, bool useSecondElement = true) {
    // special case of function above where reference is W/Z and a lepton is used
    // the boson is lep+lep2, so it is more convenient to pass the two leptons rather than the precooked boson

    int id1 = 0;
    int id2 = 1;
    if (useSecondElement) {
        id1 = 1;
        id2 = 0;
    }
    TVector2 lep_unity(std::cos(phi[id1]), std::sin(phi[id1]));
    TVector2 boson(pt[id2]*std::cos(phi[id2]), pt[id2]*std::sin(phi[id2])); // start defining boson as the pther lepton
    boson += (pt[id1]*lep_unity); // before summing get lepton with actual magnitude
    return (lep_unity*boson);

}

double parallelProjectionLepBoson(const float pt1, const float phi1, const float pt2, const float phi2) {
    // as above, but passing explicitly the components, so no ambiguity to select the lepton to project on the Z

    TVector2 lep_unity(std::cos(phi1), std::sin(phi2));
    TVector2 boson(pt2*std::cos(phi2), pt2*std::sin(phi2)); // start defining boson as the pther lepton
    boson += (pt1*lep_unity); // before summing get lepton with actual magnitude
    return (lep_unity*boson);

}


float met_cal(float met_pt, float met_phi, float lep_pt, float lep_phi, float u_coeff, float u_syst)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float ux = met_px + lep_px, uy = met_py + lep_py;
    float metcal_px = - u_coeff*ux*(1+u_syst) - lep_px, metcal_py = - u_coeff*uy*(1+u_syst) - lep_py;
    return hypot(metcal_px,metcal_py);
}


//==================================================


bool valueInsideRange(float value, float low, float high) {

  if (value > low and value < high) return true;
  else                              return false;

}


//==================================================

Vec_b valueInsideRange(const Vec_f& value, float low, float high) {

  Vec_b ret(value.size(), false);
  for (unsigned int i = 0; i < ret.size(); ++i) {
    if (value[i] >= low and value[i] < high) ret[i] = true;
  }
  return ret;

}

//==================================================

bool valueOutsideRange(float value, float low, float high) {

  if (value < low or value > high) return true;
  else                             return false;

}

//==================================================

Vec_b valueOutsideRange(const Vec_f& value, float low, float high) {

  Vec_b ret(value.size(), false);
  for (unsigned int i = 0; i < ret.size(); ++i) {
    if (value[i] < low or value[i] > high) ret[i] = true;
  }
  return ret;

}

//==================================================

float varLepPlusFromPair(float var1, int pdgid1, float var2, int pdgid2) {

  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0) return -9999.0;

  if (pdgid1 > 0) return var2;
  else            return var1;

}

//==================================================
float varLepMinusFromPair(float var1, int pdgid1, float var2, int pdgid2) {

  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0) return -9999.0;

  if (pdgid1 > 0) return var1;
  else            return var2;

}

//==================================================
float varChargedLepFromPair(int requestedCharge, float var1, int pdgid1, float var2, int pdgid2) {

  // requestedCharge must be > 0 for positive charge and < 0 otherwise

  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0)   return -9999.0;
  if (requestedCharge > 0) return (pdgid1 < 0) ? var1 : var2;
  else                     return (pdgid1 > 0) ? var1 : var2;

}

//==================================================
TRandom3 *randy_v2 = NULL;
double randomVarFromPair(float var1, float var2) {

  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return 0 if not
  if (!randy_v2) randy_v2 = new TRandom3(0);
  if (randy_v2->Rndm() > 0.5) return var1;
  else                        return var2;

}

double ptDiffCharge(float pt1, int charge1, float pt2, int charge2) {

  if (charge1 < 0) return pt1-pt2;
  else             return pt2-pt1;

}



void functions() {}

}

#endif
