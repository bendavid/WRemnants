#include <ROOT/RVec.hxx>
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TFile.h"
#include "TTree.h"
#include <eigen3/Eigen/Dense>
#include <memory>

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

    x_->reserve(parmtree->GetEntries());

    for (long long ientry=0; ientry < parmtree->GetEntries(); ++ientry) {
      parmtree->GetEntry(ientry);
      x_->push_back(x);
    }
  }


  std::pair<V, int> operator() (float pt, float eta, float phi, int charge, const RVec<int> &idxs, const RVec<float> &jac) {

    if (pt < 0.) {
      return std::make_pair<V, int>(V(), -99);
    }

    const double theta = 2.*std::atan(std::exp(-double(eta)));
    const double lam = M_PI_2 - theta;
    const double p = double(pt)/std::sin(theta);
    const double qop = double(charge)/p;

    const Eigen::Matrix<double, 3, 1> curvmom(qop, lam, phi);

    const auto nparms = idxs.size();

    const Eigen::Map<const Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::RowMajor>> jacMap(jac.data(), 3, nparms);

    Eigen::VectorXd xtrk(nparms);
    for (unsigned int i = 0; i < nparms; ++i) {
      xtrk[i] = (*x_)[(*idxmap_)[idxs[i]]];
    }

    const Eigen::Matrix<double, 3, 1> curvmomcor = curvmom + jacMap.cast<double>()*xtrk;

    const double qopcor = curvmomcor[0];
    const double lamcor = curvmomcor[1];
    const double phicor = curvmomcor[2];

    const double pcor = 1./std::abs(qopcor);
    const int qcor = std::copysign(1., qopcor);

    const double ptcor = pcor*std::cos(lamcor);

    const double thetacor = M_PI_2 - lamcor;
    const double etacor = -std::log(std::tan(0.5*thetacor));

    return std::make_pair<V, int>(V(ptcor, etacor, phicor, wrem::muon_mass), int(qcor));

  }

protected:

  std::shared_ptr<std::vector<unsigned int>> idxmap_ = std::make_shared<std::vector<unsigned int>>();
  std::shared_ptr<std::vector<double>> x_ = std::make_shared<std::vector<double>>();
};


class CVHCorrector : public CVHCorrectorSingle {
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

  calibration_uncertainty_helper(HIST &&hist) : hist_(std::make_shared<const HIST>(std::move(hist))) {}

  out_tensor_t operator() (const RVec<float> &pts,
                      const RVec<float> &etas,
                      const RVec<float> &phis,
                      const RVec<int> &charges,
                      const RVec<int> &matchedGenIdxs,
                      const RVec<RVec<float>> &momcovs,
                      const RVec<int> &selection,
                      const RVec<float> &genPts,
                      const RVec<float> &genEtas,
                      const RVec<float> &genPhis,
                      const RVec<int> &genPdgIds,
                      const RVec<int> &genStatusFlags,
                      double nominal_weight = 1.0) const {

    //TODO move this into a helper

    const std::size_t nmuons = pts.size();

    out_tensor_t res;
    res.setConstant(nominal_weight);

    for (std::size_t i = 0; i < nmuons; ++i) {
      if (!selection[i]) {
        continue;
      }

      if (pts[i] < 0.) {
        continue;
      }

      const int genidx = matchedGenIdxs[i];

      if (genidx < 0) {
        continue;
      }

      // matched genparts should be status muons, but we need to explicitly check if they are prompt
      if (!(genStatusFlags[genidx] & 0x01)) {
        continue;
      }

      const int gencharge = genPdgIds[genidx] > 0 ? -1 : 1;

      res *= scale_res_weight(pts[i], etas[i], phis[i], charges[i], genPts[genidx], genEtas[genidx], genPhis[genidx], gencharge, momcovs[i]);
    }

    return res;

  }

private:

  out_tensor_t scale_res_weight(double pt, double eta, double phi, int charge, double genPt, double genEta, double genPhi, int genCharge, const RVec<float> &cov) const {

    const double theta = 2.*std::atan(std::exp(-double(eta)));
    const double lam = M_PI_2 - theta;
    const double p = double(pt)/std::sin(theta);
    const double qop = double(charge)/p;

    const Eigen::Vector3d parms(qop, lam, phi);

    const double gentheta = 2.*std::atan(std::exp(-double(genEta)));
    const double genlam = M_PI_2 - gentheta;
    const double genp = double(genPt)/std::sin(gentheta);
    const double genqop = double(genCharge)/genp;
    const double genqopt = double(genCharge)/genPt;

    const Eigen::Vector3d genparms(genqop, genlam, genPhi);

    const Eigen::Vector3d deltaparms = parms - genparms;

    const Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>> covMap(cov.data(), 3, 3);

    Eigen::Matrix<double, 3, 3> covd = covMap.cast<double>();
    // fill in lower triangular part of the matrix, which is stored as zeros to save space
    covd.triangularView<Eigen::Lower>() = covd.triangularView<Eigen::Upper>().transpose();

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
      parmvar[0] += A*genqop;
      parmvar[0] += -e*genqop/genPt;
      parmvar[0] += genCharge*M*genPt*genqop;

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

        res(ivar, idownup) = weight;

      }
    }

    return res;
  }

  std::shared_ptr<const HIST> hist_;

};


}
