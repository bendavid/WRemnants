#include <boost/histogram.hpp>
#include <stdlib.h>

namespace wrem {
template <typename T>
class JpsiCorrectionsHelper {

public:

    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;
    static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
    static constexpr auto nUnc = sizes[sizes.size() - 1]; // 1 for cnetral value
    using out_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nUnc, 2>>;

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
    double operator() (float cvh_eta, float cvh_pt, int charge) {
        const auto &params = get_tensor(cvh_eta);
        const double A = params(0, 0);
        const double e = params(1, 0);
        const double M = params(2, 0);
        double k = 1.0 / cvh_pt;
        double magnetic = 1.0 + A;
        double material = -1.0 * e * k;
        double alignment = charge * M;
        double k_crctd = (magnetic + material) * k + alignment;
        return (1.0 / k_crctd);
    }

    // for uncertainties on pt
    out_tensor_t operator() (
        float cvh_eta, float cvh_pt, int charge, float jpsi_crctd_pt
    ) {
        const auto &params = get_tensor(cvh_eta);
        double k = 1.0 / cvh_pt;
        out_tensor_t res;
        for (int i = 0; i < nUnc; i++) {
            const double A_unc = params(0, i);
            const double e_unc = params(1, i);
            const double M_unc = params(2, i);
            double k_crctd = 1.0 / jpsi_crctd_pt;
            double k_shift = (A_unc - e_unc * k) * k + charge * M_unc;
            double k_crctd_up = k_crctd + k_shift, k_crctd_down = k_crctd - k_shift;
            res(i, 0) = (1.0 / k_crctd_up) - cvh_pt;
            res(i, 1) = (1.0 / k_crctd_down) - cvh_pt;
        }
        return res;
    }

    // for smearing weights derived from qop
    out_tensor_t operator() (
        float cvh_eta, float cvh_pt, int cvh_charge, float jpsi_crctd_pt,
        double qop, double pt, double eta, double phi, int charge
    ) {
        const auto &params = get_tensor(cvh_eta);
        double k = 1.0 / cvh_pt;
        out_tensor_t res;
        for (int i = 0; i < nUnc; i++) {
            const double A_unc = params(0, i);
            const double e_unc = params(1, i);
            const double M_unc = params(2, i);
            double k_crctd = 1.0 / jpsi_crctd_pt;
            double k_shift = (A_unc - e_unc * k) * k + cvh_charge * M_unc;
            double k_crctd_up = k_crctd + k_shift, k_crctd_down = k_crctd - k_shift;
            res(i, 0) = (1.0 / k_crctd_up) - cvh_pt;
            res(i, 1) = (1.0 / k_crctd_down) - cvh_pt;
        }
        return res;
    }

private:

  out_tensor_t smearing_weight(
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
    std::shared_ptr<const T> correctionHist_;

};

}
