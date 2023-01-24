#include <ROOT/RVec.hxx>
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include <eigen3/Eigen/Dense>
#include <memory>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <typeinfo>
#include <algorithm>

namespace wrem {

using ROOT::VecOps::RVec;

template<typename T>
ROOT::VecOps::RVec<ROOT::VecOps::RVec<T>> splitNestedRVec(const ROOT::VecOps::RVec<T> &vec, const ROOT::VecOps::RVec<int> &counts) {
  using ROOT::VecOps::RVec;

  RVec<RVec<T>> res;
  res.reserve(counts.size());

  int total = 0;
  for (unsigned int i = 0; i < counts.size(); ++i) {
    const int count = counts[i];
    // avoid copies by adopting memory
    res.emplace_back(const_cast<T*>(vec.begin()) + total, count);
    total += count;
  }

  return res;
}

template <std::ptrdiff_t NJac = 3, std::ptrdiff_t NReplicas = 0>
class CVHCorrectorSingle {
public:

  using V = ROOT::Math::PtEtaPhiM4D<double>;

  CVHCorrectorSingle(const std::string &filename) {
    TFile *fcor = TFile::Open(filename.c_str());

    TTree *idxmaptree = static_cast<TTree*>(fcor->Get("idxmaptree"));

    unsigned int idx;
    idxmaptree->SetBranchAddress("idx", &idx);

    idxmap_->reserve(idxmaptree->GetEntries());

    for (long long ientry=0; ientry < idxmaptree->GetEntries(); ++ientry) {
      idxmaptree->GetEntry(ientry);
      idxmap_->push_back(idx);
    }

    TTree *parmtree = static_cast<TTree*>(fcor->Get("parmtree"));

    float x;
    parmtree->SetBranchAddress("x", &x);

    std::array<float, NReplicas> xrep;
    if (NReplicas>0) {
      parmtree->SetBranchAddress("xreplicas", xrep.data());
    }

    x_->reserve(parmtree->GetEntries());
    xreplicas_->resize(parmtree->GetEntries(), NReplicas);

    for (long long ientry=0; ientry < parmtree->GetEntries(); ++ientry) {
      parmtree->GetEntry(ientry);
      x_->push_back(x);

      if (NReplicas > 0) {
        for (std::ptrdiff_t irep = 0; irep < NReplicas; ++irep) {
          (*xreplicas_)(ientry, irep) = xrep[irep];
        }
      }
    }
  }


  std::pair<V, int> operator() (float pt, float eta, float phi, int charge, const RVec<int> &idxs, const RVec<float> &jac) const {

    if (pt < 0.) {
      return std::make_pair<V, int>(V(), -99);
    }

    const Eigen::Matrix<double, 3, 1> curvmom = CurvMom(pt, eta, phi, charge);

    return CorMomCharge(curvmom, idxs, jac, *x_);
  }

protected:

  Eigen::Matrix<double, 3, 1> CurvMom(float pt, float eta, float phi, int charge) const {
    const double theta = 2.*std::atan(std::exp(-double(eta)));
    const double lam = M_PI_2 - theta;
    const double p = double(pt)/std::sin(theta);
    const double qop = double(charge)/p;

    return Eigen::Matrix<double, 3, 1>(qop, lam, phi);
  }

  template<typename U>
  std::pair<V, int> CorMomCharge(const Eigen::Matrix<double, 3, 1> &curvmom, const RVec<int> &idxs, const RVec<float> &jac, const U &parms) const {

    const auto nparms = idxs.size();

    const Eigen::Map<const Eigen::Matrix<float, NJac, Eigen::Dynamic, Eigen::RowMajor>> jacMap(jac.data(), NJac, nparms);

    Eigen::VectorXd xtrk(nparms);
    for (unsigned int i = 0; i < nparms; ++i) {
      xtrk[i] = parms[(*idxmap_)[idxs[i]]];
    }

    const Eigen::Matrix<double, 3, 1> curvmomcor = curvmom + jacMap.template topRows<3>().template cast<double>()*xtrk;

    if (curvmomcor.array().isNaN().any() || curvmomcor.array().isInf().any()) {
      return std::make_pair<V, int>(V(), -99);
    }

    const double qopcor = curvmomcor[0];
    const double lamcor = std::clamp(curvmomcor[1], -M_PI_2, M_PI_2);
    const double phicor = curvmomcor[2];

    const double pcor = 1./std::abs(qopcor);
    const int qcor = std::copysign(1., qopcor);

    const double ptcor = pcor*std::cos(lamcor);

    const double thetacor = std::clamp(M_PI_2 - lamcor, 0., M_PI);
    const double etacor = -std::log(std::tan(0.5*thetacor));

    return std::make_pair<V, int>(V(ptcor, etacor, phicor, wrem::muon_mass), int(qcor));

  }

  std::shared_ptr<std::vector<unsigned int>> idxmap_ = std::make_shared<std::vector<unsigned int>>();
  std::shared_ptr<std::vector<double>> x_ = std::make_shared<std::vector<double>>();
  std::shared_ptr<Eigen::MatrixXd> xreplicas_ = std::make_shared<Eigen::MatrixXd>();
};


class CVHCorrectorValidation : public CVHCorrectorSingle<5, 0> {
public:
  using base_t = CVHCorrectorSingle<5, 0>;

  using V = typename base_t::V;

  using base_t::base_t;

  std::pair<V, int> operator() (const RVec<float> &refParms, const RVec<unsigned int> &idxs, const RVec<float> &jac) const {

    const Eigen::Matrix<double, 3, 1> curvmom(refParms[0], refParms[1], refParms[2]);

    return CorMomCharge(curvmom, idxs, jac, *base_t::x_);
  }
};

template <std::ptrdiff_t NReplicas>
class CVHCorrectorReplicas : public CVHCorrectorSingle<5, NReplicas> {
public:
  using base_t = CVHCorrectorSingle<5, NReplicas>;

  using V = typename base_t::V;

  using base_t::base_t;

  std::array<std::pair<V, int>, NReplicas> operator() (const RVec<float> &refParms, const RVec<unsigned int> &idxs, const RVec<float> &jac) const {

    std::array<std::pair<V, int>, NReplicas> res;

    const Eigen::Matrix<double, 3, 1> curvmom(refParms[0], refParms[1], refParms[2]);

    for (std::ptrdiff_t irep = 0; irep < NReplicas; ++irep) {
      res[irep] = CorMomCharge(curvmom, idxs, jac, base_t::xreplicas_->col(irep));
    }

    return res;

  }
};


class CVHCorrector : public CVHCorrectorSingle<> {
public:


  CVHCorrector(const std::string &filename) : CVHCorrectorSingle(filename) {}
  CVHCorrector(const CVHCorrectorSingle &corsingle) : CVHCorrectorSingle(corsingle) {}

  RVec<std::pair<V, int>> operator () (const RVec<float> &ptv, const RVec<float> &etav, const RVec<float> &phiv, const RVec<int> &chargev, const RVec<RVec<int>> &idxsv, const RVec<RVec<float>> &jacv) {
    RVec<std::pair<V, int>> res;
    res.reserve(ptv.size());

    for (unsigned int i = 0; i < ptv.size(); ++i) {
      res.emplace_back(CVHCorrectorSingle::operator()(ptv[i], etav[i], phiv[i], chargev[i], idxsv[i], jacv[i]));
    }

    return res;
  }
};

template <typename HIST>
class calibration_uncertainty_helper {

public:

  using tensor_t = typename HIST::storage_type::value_type::tensor_t;
  static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
  static constexpr auto nvars = sizes[0];
  static constexpr auto nparms = sizes[1];

  using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nvars, 2>>;

  calibration_uncertainty_helper(HIST &&hist, double minGenPt, double maxWeight) : 
    hist_(std::make_shared<const HIST>(std::move(hist))),
    minGenPt_(minGenPt),
    maxWeight_(maxWeight) {}

  out_tensor_t operator() (
    double qop,
    double pt,
    float eta,
    float phi,
    int charge,
    RVec<float> &momcov,
    double genQop,
    float genPt,
    float genEta,
    float genPhi,
    int genCharge,
    double nominal_weight) const {

    //TODO move this into a helper

    out_tensor_t res;
    res.setConstant(nominal_weight);
    res *= scale_res_weight(qop, pt, eta, phi, charge, momcov, genQop, genPt, genEta, genPhi, genCharge);
    return res;
  }

private:

  out_tensor_t scale_res_weight(
    double qop, double pt, double eta, double phi, int charge, const RVec<float> &cov,
    double genQop, double genPt, double genEta, double genPhi, int genCharge,
    bool abQop = false, bool fullParam = false
  ) const {

    const double theta = 2.*std::atan(std::exp(-double(eta)));
    const double lam = M_PI_2 - theta;
    const double p = double(pt)/std::sin(theta);
    const double Qop = double(charge)/p;
    const Eigen::Vector3d parms(
      (abQop? qop : Qop),
      (fullParam? lam : 0),
      (fullParam? phi : 0)
    );

    const double gentheta = 2.*std::atan(std::exp(-double(genEta)));
    const double genlam = M_PI_2 - gentheta;
    const double genp = double(genPt)/std::sin(gentheta);
    const double genqop = double(genCharge)/genp;
    const Eigen::Vector3d genparms(
      (abQop? genQop : genqop),
      (fullParam? genlam : 0),
      (fullParam? genPhi : 0)
    );

    const Eigen::Vector3d deltaparms = parms - genparms;

    const Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>> covMap(cov.data(), 3, 3);

    Eigen::Matrix<double, 3, 3> covd = covMap.cast<double>();

    if (fullParam) {
      // fill in lower triangular part of the matrix, which is stored as zeros to save space
      covd.triangularView<Eigen::Lower>() = covd.triangularView<Eigen::Upper>().transpose();
    } else {
      covd.row(0) << covd(0,0), 0, 0;
      covd.row(1) << 0, 1, 0;
      covd.row(2) << 0, 0, 1;

    }

    const Eigen::Matrix<double, 3, 3> covinv = covd.inverse();
    const double covdet = covd.determinant();

    const double lnp = -0.5*deltaparms.transpose()*covinv*deltaparms;

    // no need to initialize since all elements are explicit filled
    out_tensor_t res;

    const auto &varparms = hist_->at(hist_->template axis<0>().index(genEta)).data();

    for (std::ptrdiff_t ivar = 0; ivar < nvars; ++ivar) {
      Eigen::Vector3d parmvar = Eigen::Vector3d::Zero();
      Eigen::Matrix<double, 3, 3> covvar = Eigen::Matrix<double, 3, 3>::Zero();

      const double A = varparms(ivar, 0);
      const double e = varparms(ivar, 1);
      const double M = varparms(ivar, 2);
      const double sig = varparms(ivar, 3);

      // scale
      parmvar[0] += (abQop? A*genQop : A*genqop);
      parmvar[0] += (abQop? -e*genQop/genPt : -e*genqop/genPt);
      parmvar[0] += (abQop? genCharge*M*genPt*genQop : genCharge*M*genPt*genqop);

      // resolution
      covvar += sig*covd;

      for (std::ptrdiff_t idownup = 0; idownup < 2; ++idownup) {

        const double dir = idownup == 0 ? -1. : 1.;

        const Eigen::Vector3d deltaparmsalt = deltaparms + dir*parmvar;
        const Eigen::Matrix<double, 3, 3> covdalt = covd + dir*covvar;

        const Eigen::Matrix<double, 3, 3> covinvalt = covdalt.inverse();
        const double covdetalt = covdalt.determinant();

        const double lnpalt = -0.5*deltaparmsalt.transpose()*covinvalt*deltaparmsalt;

        const double weight = std::sqrt(covdet/covdetalt)*std::exp(lnpalt - lnp);

        if (false) {
          if (std::isnan(weight) || std::isinf(weight) || std::fabs(weight) == 0.) {
            std::cout << "invalid weight: " << weight << std::endl;
            std::cout << "pt " << pt << std::endl;
            std::cout << "genPt " << genPt << std::endl;
            std::cout << "qop " << qop << std::endl;
            std::cout << "lam " << lam << std::endl;
            std::cout << "genlam " << genlam << std::endl;
            std::cout << "phi " << phi << std::endl;
            std::cout << "genlam " << genPhi << std::endl;
            std::cout << "covdet " << covdet << std::endl;
            std::cout << "covdetalt " << covdet << std::endl;
            std::cout << "lnp " << lnp << std::endl;
            std::cout << "lnpalt " << lnpalt << std::endl;
            std::cout << "covd\n" << covd << std::endl;
            std::cout << "covdalt\n" << covdalt << std::endl;
            std::cout << "covinv\n" << covinv << std::endl;
            std::cout << "covinvalt\n" << covinvalt << std::endl;
          }
        }

        // protect against outliers
        // if (weight > 0.9998 && weight < 1.0002) {cout << "smearing weight is " << weight << "covd is " << covd << std::endl;}
        res(ivar, idownup) = std::min(weight, maxWeight_);
      }
    }
    
    return res;
  }
  std::shared_ptr<const HIST> hist_;
  double minGenPt_;
  double maxWeight_;

};

    float smearGenPt(RVec<float> cov, int charge, float pt, float theta) {
        float sigma2_qop = cov[0];
        double sigma2_pt = pow((pt * pt / (charge * sin(theta))), 2) * sigma2_qop;
        return gRandom -> Gaus(pt, sqrt(sigma2_pt));
    }
    
    double smearGenQop(RVec<float> cov, double qop) {
        double sigma2_qop = cov[0];
        return gRandom -> Gaus(qop, sqrt(sigma2_qop));
    }

    int getGoodGenMuons0IdxInReco(RVec<bool> goodMuons, RVec<bool> goodMuonsByGenTruth) {
        int first_goodGenMuon_idx_in_goodMuons = 0;
        for (int i = 0; i < goodMuonsByGenTruth.size(); i++) {
            if (goodMuonsByGenTruth[i]) {break;}
            else {first_goodGenMuon_idx_in_goodMuons += 1;}
        }
        int goodGenMuons0_idx_in_recos = 0;
        int num_goodMuons_found = 0;
        for (int i = 0; i < goodMuons.size(); i++) {
            if (goodMuons[i]) {
                num_goodMuons_found += 1;
                if (num_goodMuons_found == first_goodGenMuon_idx_in_goodMuons + 1) {
                    break;
                }
            }
            goodGenMuons0_idx_in_recos += 1;
        }
        return goodGenMuons0_idx_in_recos;
    }

    RVec<float> getCovMatForGoodMuons0(
        RVec<float> covmat,
        RVec<int> covmat_counts,
        RVec<bool> goodMuons,
        RVec<bool> goodMuonsByGenTruth
    ) {
        int goodGenMuons0_idx_in_recos = wrem::getGoodGenMuons0IdxInReco(
            goodMuons, goodMuonsByGenTruth
        );
        int covmat_start_idx = 0;
        for (int i = 0; i < goodGenMuons0_idx_in_recos; i++) {
            covmat_start_idx += covmat_counts[i];
        }
        ROOT::VecOps::RVec<int> idxRange(9);
        for (int i = 0; i < 9; i++) {
            idxRange[i] = i + covmat_start_idx;
        }
        return Take(covmat, idxRange);
    }
  
template <typename T>
class BiasCorrectionsHelper {

public:

    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;
    static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
    static constexpr auto nUnc = sizes[sizes.size() - 1]; // 1 for cnetral value
    using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nUnc, 2>>;

    BiasCorrectionsHelper(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    // helper for bin lookup which implements the compile-time loop over axes
    template<typename... Xs, std::size_t... Idxs>
    const tensor_t &get_tensor_impl(std::index_sequence<Idxs...>, const Xs&... xs) {
        return correctionHist_->at(correctionHist_->template axis<Idxs>().index(xs)...).data();
    }

    // variadic templated bin lookup
    template<typename... Xs>
    const tensor_t &get_tensor(const Xs&... xs) {
        return get_tensor_impl(std::index_sequence_for<Xs...>{}, xs...);
    }
    
    // for central value of pt
    float operator() (float eta, float pt, int charge) {
        const double bias = get_tensor(eta, pt);
        return (1.0 + bias) * pt;
    }

private:
    std::shared_ptr<const T> correctionHist_;
};

}
