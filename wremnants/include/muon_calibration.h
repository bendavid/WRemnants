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
#include "defines.h"

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

  using V = ROOT::Math::PtEtaPhiMVector;

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
    Eigen::Vector3d parms(
      (abQop? qop : Qop),
      (fullParam? lam : 0),
      (fullParam? phi : 0)
    );

    const double gentheta = 2.*std::atan(std::exp(-double(genEta)));
    const double genlam = M_PI_2 - gentheta;
    const double genp = double(genPt)/std::sin(gentheta);
    const double genqop = double(genCharge)/genp;
    Eigen::Vector3d genparms(
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

    RVec<bool> filterRecoMuonsByGenTruth(
        RVec<bool> muonSelection,
        RVec<int> Muon_genPartIdx,
        RVec<int> GenPart_pdgId,
        RVec<int> GenPart_statusFlags,
        RVec<int> GenPart_status,
        bool requirePrompt = true
    ) {
        ROOT::VecOps::RVec<bool> res(muonSelection.size());
        for (int i = 0; i < muonSelection.size(); i++) {
            res[i] = (
                muonSelection[i] &&
                (!(requirePrompt) || GenPart_statusFlags[Muon_genPartIdx[i]] & 0x01) &&
                (abs(GenPart_pdgId[Muon_genPartIdx[i]]) == MUON_PDGID) &&
                (GenPart_status[Muon_genPartIdx[i]] == 1)
            );
        }
        return res;
    }

    RVec<bool> filterRecoMuonsByCovMat(
        RVec<bool> muonSelection,
        RVec<float> covMat,
        RVec<int> covMatCounts
    ) {
        ROOT::VecOps::RVec<bool> res(muonSelection.size());
        int covMat_start_idx = 0;
        for (int i = 0; i < muonSelection.size(); i++) {
            res[i] = (
                muonSelection[i] &&
                //covMatCounts[i] &&
                covMat[covMat_start_idx] > 0
            );
            covMat_start_idx += covMatCounts[i];
        }
        return res;
    }

    template<typename T> RVec<RVec<T>> getCovMatForSelectedRecoMuons(
        RVec<float> covmat,
        RVec<int> covmat_counts,
        RVec<bool> selection,
        int nMuons = -1 // limit the number of muons to be N; -1 to take all selected muons
    ) {
        if (covmat_counts.size() != selection.size()) {
            throw std::invalid_argument("selection masks should match the size of all RECO muons");
        }
        
        using ROOT::VecOps::RVec;
        RVec<RVec<T>> res;
        res.reserve(nMuons > 0 ? nMuons : 1);

        int covmat_start_idx = 0;
        for (int i = 0; i < covmat_counts.size(); i++) {
            if (selection[i]) {
                ROOT::VecOps::RVec<int> idxRange(covmat_counts[i]);
                for (int i = 0; i < idxRange.size(); i++) {
                    idxRange[i] = i + covmat_start_idx;
                }
                res.emplace_back(Take(covmat, idxRange));
            }
            covmat_start_idx += covmat_counts[i];
        }
        return res;
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
  
    using ROOT::VecOps::RVec;

    double calculateTheta(float eta) {
        return 2.*std::atan(std::exp(-double(eta)));
    }

    double calculateLam(float eta) {
        double theta = calculateTheta(eta);
        return M_PI_2 - theta;
    }

    double calculateQop(float pt, float eta, int charge) {
        double theta = calculateTheta(eta);
        return (charge * std::sin(theta) / double(pt));
    }

    double calculatePt(double qop, float eta, int charge) {
        double theta = calculateTheta(eta);
        return (charge * std::sin(theta) / qop);
    }

    double calculateQopUnc(float eta, int charge, double KUnc) {
        double theta = calculateTheta(eta);
        return (charge * std::sin(theta) * KUnc);
    }

    double calculateQopUnc(float pt, float eta, int charge, double ptUnc) {
        double theta = calculateTheta(eta);
        return ((-1. * charge * std::sin(theta)) / pow(pt, 2)) * ptUnc;
    }

// jpsi correction central value for one muon
template <typename T>
class JpsiCorrectionsHelper {

public:
    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;

    JpsiCorrectionsHelper(T&& corrections) :
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
    float operator() (float cvhPt, float cvhEta, int charge) {
        const auto &params = get_tensor(cvhEta);
        const double A = params(0);
        const double e = params(1);
        const double M = params(2);
        double k = 1.0 / cvhPt;
        double magnetic = 1.0 + A;
        double material = -1.0 * e * k;
        double alignment = charge * M;
        double kCrctd = (magnetic + material) * k + alignment;
        return (1.0 / kCrctd);
    }

private:
    std::shared_ptr<const T> correctionHist_;
};

// jpsi corrections central value for multiple muons
template <typename T>
class JpsiCorrectionsRVecHelper : public JpsiCorrectionsHelper<T> {

using base_t = JpsiCorrectionsHelper<T>;

public:
    //inherit constructor
    using base_t::base_t;

    RVec<float> operator() (const RVec<float>& pts, const RVec<float> etas, RVec<int> charges) {
        RVec<float> corrected_pt(pts.size(), 0.);
        assert(etas.size() == pts.size() && etas.size() == charges.size());
        for (size_t i = 0; i < pts.size(); i++) {
            corrected_pt[i] = JpsiCorrectionsHelper<T>::operator()(pts[i], etas[i], charges[i]);
        }

        return corrected_pt;
    }
};

// unlike the JpsiCorrectionHelper which only deals with one muon,
// and uses inheritance to work on a vector of muons,
// this helper is already vectorized
template <typename T>
class JpsiCorrectionsUncHelper {

public:
    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;
    static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
    static constexpr auto nUnc = sizes[sizes.size() - 1]; // 1 for central value
    using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nUnc, 2>>;

    JpsiCorrectionsUncHelper(T&& corrections) :
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

    // for smearing weights derived from propagating uncs on A, e, M to uncs on qop

    out_tensor_t operator() (
        const RVec<double> &genQops, 
        const RVec<float> &genPhis,
        const RVec<float> &genEtas,
        const RVec<double> &recoQops,
        const RVec<float> &recoPhis,
        const RVec<float> &recoEtas,
        const RVec<int> &recoCharges,
        const RVec<double> &recoPts,
        const RVec<RVec<float>> &covs, // for sigma on the Gaussian
        double nominal_weight = 1.0,
        bool fullParam = false
    ) {
        out_tensor_t res;
        res.setConstant(nominal_weight);

        const std::size_t nmuons = recoPts.size();

        for (std::size_t i = 0; i < nmuons; ++i) {
            res *= smearingWeight_oneMuon(
                genQops[i], genPhis[i], genEtas[i],
                recoQops[i], recoPhis[i], recoEtas[i], recoCharges[i], recoPts[i],
                covs[i], fullParam
            );
        }

        return res;
    }

private:
    out_tensor_t smearingWeight_oneMuon(
        double genQop,  float genPhi,  float genEta,
        double recoQop, float recoPhi, float recoEta, int recoCharge, float recoPt,  
        const RVec<float> &cov, bool fullParam = false
    ) {
        Eigen::Vector3d parms(
            recoQop,
            (fullParam? calculateLam(recoEta) : 0),
            (fullParam? recoPhi: 0)
        ); // (qop, lam, phi)

        Eigen::Vector3d genparms(
            genQop,
            (fullParam? calculateLam(genEta) : 0),
            (fullParam? genPhi: 0)
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

        const auto &params = get_tensor(recoEta);

        // no need to initialize since all elements will be explicitly filled
        out_tensor_t res;
    
        for (std::ptrdiff_t ivar = 0; ivar < nUnc; ++ivar) {
            const double AUnc = params(0, ivar);
            const double eUnc = params(1, ivar);
            const double MUnc = params(2, ivar);
            double recoK = 1.0 /recoPt;
            double recoKUnc = (AUnc - eUnc * recoK) * recoK + recoCharge * MUnc;
            Eigen::Vector3d parmvar = Eigen::Vector3d::Zero();   
            parmvar[0] = recoCharge * std::sin(calculateTheta(recoEta)) * recoKUnc;

            for (std::ptrdiff_t idownup = 0; idownup < 2; ++idownup) {
                const double dir = idownup == 0 ? -1. : 1.;
                const Eigen::Vector3d deltaparmsalt = deltaparms + dir*parmvar;
                const Eigen::Matrix<double, 3, 3> covdalt = covd; //+ dir*covvar;
    
                const Eigen::Matrix<double, 3, 3> covinvalt = covdalt.inverse();
                const double covdetalt = covdalt.determinant();
    
                const double lnpalt = -0.5*deltaparmsalt.transpose()*covinvalt*deltaparmsalt;
    
                const double weight = std::sqrt(covdet/covdetalt)*std::exp(lnpalt - lnp);
    
                res(ivar, idownup) = weight;
            }
        }
        return res;
    }

    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class ZNonClosureChargeDepHelper {

public:
    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;
    using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<2>>;

    ZNonClosureChargeDepHelper(T&& corrections) :
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

    out_tensor_t operator() (
        double genQop, float genEta,
        double recoQop, float recoEta, int recoCharge, double recoPt,  
        const RVec<float> &cov, double nominal_weight = 1.0,
        std::vector<std::string> vars =  {"A", "M"}
    ) {
        const auto &params = get_tensor(recoEta);
        double recoK = 1.0 /recoPt;
        double recoKUnc = 0.0;

        if (std::find(vars.begin(), vars.end(), "A") != vars.end()) {
            const double AUnc = params(0);
            recoKUnc += AUnc * recoK;
        }
        if (std::find(vars.begin(), vars.end(), "e") != vars.end()) {
            const double eUnc = params(1);
            recoKUnc += -1.0 * eUnc * recoK * recoK;
        }
        if (std::find(vars.begin(), vars.end(), "M") != vars.end()) {
            const double MUnc = params(2);
            recoKUnc += recoCharge * MUnc;
        }

        double recoQopUnc = recoCharge * std::sin(calculateTheta(recoEta)) * recoKUnc;

        const Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>> covMap(cov.data(), 3, 3);
        Eigen::Matrix<double, 3, 3> covd = covMap.cast<double>();
        const double sigma2Qop = covd(0,0);
        const double sigma2QopAlt = sigma2Qop;

        const double lnp = -0.5 * pow((recoQop - genQop), 2) / sigma2Qop;

        out_tensor_t res;
        res.setConstant(nominal_weight);

        for (std::ptrdiff_t idownup = 0; idownup < 2; ++idownup) {
            const double dir = idownup == 0 ? -1. : 1.;
            const double lnpAlt = -0.5 * pow((recoQop + dir * recoQopUnc - genQop), 2) / sigma2QopAlt;
            const double weight = std::sqrt(sigma2Qop / sigma2QopAlt) * std::exp(lnpAlt - lnp);
            res(idownup) *= weight; 
        }
        return res;
    }

private:
    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class ZNonClosureChargeIndHelper {

public:
    using hist_t = T;
    using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<2>>;

    ZNonClosureChargeIndHelper(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    // smaering weights for Z non-closure number on pt 
    out_tensor_t operator() (
        double genQop, float genPt, float genEta, int genCharge,
        double recoQop, double recoPt, float recoEta, int recoCharge,
        const RVec<float> &cov, double nominal_weight = 1.0
    ) {
        out_tensor_t res;
        res.setConstant(nominal_weight);
        /*
        if (recoPt <= 25.0 || recoPt >= 100.0 || abs(recoEta) >= 2.4 ) {
            return res;
        }
        */
        unsigned int iEta = correctionHist_->template axis<0>().index(recoEta);
        unsigned int iPt = correctionHist_->template axis<1>().index(recoPt);
        const double nonClosure = correctionHist_->at(iEta, iPt).value();
        const double recoKUnc = (nonClosure - 1) * (1 / recoPt);
        const double recoQopUnc = calculateQopUnc(recoEta, recoCharge, recoKUnc);

        const Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>> covMap(cov.data(), 3, 3);
        Eigen::Matrix<double, 3, 3> covd = covMap.cast<double>();
        const double sigma2Qop = covd(0,0);
        const double sigma2QopAlt = sigma2Qop;

        const double lnp = -0.5 * pow((recoQop - genQop), 2) / sigma2Qop;

        for (std::ptrdiff_t idownup = 0; idownup < 2; ++idownup) {
            const double dir = idownup == 0 ? -1. : 1.;
            const double lnpAlt = -0.5 * pow((recoQop + dir * recoQopUnc - genQop), 2) / sigma2QopAlt;
            const double weight = std::sqrt(sigma2Qop / sigma2QopAlt) * std::exp(lnpAlt - lnp);
            res(idownup) *= weight; 
        }
        return res;
    }

private:
    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class CorrectionHelperBase {

public:
    CorrectionHelperBase(unsigned int nSlots, T&& corrections) :
        myRndGens(nSlots),
        correctionHist_(std::make_shared<const T>(std::move(corrections)))
        {
            init_random_generators();
        }


    CorrectionHelperBase(unsigned int nSlots, T&& corrections, T&& uncertainties) : 
        myRndGens(nSlots),
        correctionHist_(std::make_shared<const T>(std::move(corrections))),
        uncertaintyHist_(std::make_shared<const T>(std::move(uncertainties))
        ){
            init_random_generators();
        }

    void init_random_generators(){
        int seed = 1; // not 0 because seed 0 has a special meaning
        for (auto &&gen : myRndGens)
        {
            gen.SetSeed(seed++);
        }
    }

    // helper for bin lookup which implements the compile-time loop over axes
    template<typename... Xs, std::size_t... Idxs>
    const float get_value_impl(std::index_sequence<Idxs...>, const Xs&... xs) {
      return correctionHist_->at(correctionHist_->template axis<Idxs>().index(xs)...);
    }

    // variadic templated bin lookup
    template<typename... Xs>
    const float get_value(const Xs&... xs) {
        return get_value_impl(std::index_sequence_for<Xs...>{}, xs...);
    }

    // helper for bin lookup which implements the compile-time loop over axes
    template<typename... Xs, std::size_t... Idxs>
    const float get_error_impl(std::index_sequence<Idxs...>, const Xs&... xs) {
      return uncertaintyHist_->at(uncertaintyHist_->template axis<Idxs>().index(xs)...);
    }

    // variadic templated bin lookup
    template<typename... Xs>
    const float get_error(const Xs&... xs) {
        return get_error_impl(std::index_sequence_for<Xs...>{}, xs...);
    }

    virtual float get_correction(unsigned int slot, float pt, float eta) {return 0;}

    float get_random(unsigned int slot, float mean, float std){
        return myRndGens[slot].Gaus(mean, std);
    }


    RVec<float> operator() (unsigned int slot, const RVec<float>& pts, const RVec<float>& etas) {
        RVec<float> corrected_pt(pts.size(), 0.);
        assert(etas.size() == pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            corrected_pt[i] = get_correction(slot, pts[i], etas[i]);
        }
        return corrected_pt;
    }


private:
    std::vector<TRandom3> myRndGens; 

    std::shared_ptr<const T> correctionHist_;
    std::shared_ptr<const T> uncertaintyHist_;
};


template <typename T>
class BiasCalibrationHelper : public CorrectionHelperBase<T> {

using base_t = CorrectionHelperBase<T>;

public:
    //inherit constructor
    using base_t::base_t;

    float get_correction(unsigned int slot, float pt, float eta) override {
        const double bias = base_t::get_value(eta, pt);
        const double error = base_t::get_error(eta, pt);

        if(error>0.)
            return pt*(1.0 + base_t::get_random(slot, bias, error));
        else 
            return pt*(1.0+bias);
    }
};

class SmearingHelper{

public:
    SmearingHelper(unsigned int nSlots, TH2D&& smearings) :
        myRndGens(nSlots),
        hsmear(std::make_shared<const TH2D>(std::move(smearings)))
        {
            init_random_generators();
        }

    void init_random_generators(){
        int seed = 1; // not 0 because seed 0 has a special meaning
        for (auto &&gen : myRndGens)
        {
            gen.SetSeed(seed++);
        }
    }

    float get_random(unsigned int slot, float mean, float std){
        return myRndGens[slot].Gaus(mean, std);
    }


    RVec<float> operator() (unsigned int slot, const RVec<float>& pts, const RVec<float>& etas) {
        RVec<float> corrected_pt(pts.size(), 0.);
        assert(etas.size() == pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            corrected_pt[i] = get_correction(slot, pts[i], etas[i]);
        }
        return corrected_pt;
    }


    float get_correction(unsigned int slot, float pt, float eta) {

        const double xlow = hsmear->GetXaxis()->GetBinLowEdge(1);
        const double ylow = hsmear->GetYaxis()->GetBinLowEdge(1);

        const double xhigh = hsmear->GetXaxis()->GetBinUpEdge(hsmear->GetNbinsX());
        const double yhigh = hsmear->GetYaxis()->GetBinUpEdge(hsmear->GetNbinsY());

        double pt_cap = pt;
        double eta_cap = eta;

        // If eta is outside the range, set the interpolated value to the interpolated value of the closest x bin edge
        if (eta <= xlow){
          eta_cap = xlow + 0.00001;
        } 
        else if (eta >= xhigh)
        {
          eta_cap = xhigh - 0.00001;
        }

        // If pt is outside the range, set the interpolated value to the interpolated value of the closest y bin edge
        if (pt <= ylow){
          pt_cap = ylow + 0.00001;
        }
        else if (pt >= yhigh)
        {
          pt_cap = yhigh - 0.00001;
        }

        const double sigma = hsmear->Interpolate(eta_cap, pt_cap); // this is sigma_p/p

        if(sigma>0.)
            return 1. / (1./pt + get_random(slot, 0., sigma/pt));
        else 
            return pt;
    }

private:
    std::vector<TRandom3> myRndGens; 
    std::shared_ptr<const TH2D> hsmear;
};

}
