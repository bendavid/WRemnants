#include <boost/histogram.hpp>
#include <stdlib.h>

namespace wrem {

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

    double calculateDeltaQop(float pt, float deltaPt, float eta, int charge) {
        double deltaK = -deltaPt / (pt * (pt + deltaPt));
        return charge * std::sin(calculateTheta(eta)) * deltaK;
    }

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
    /*
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
            double k_unc = (A_unc - e_unc * k) * k + charge * M_unc;
            double k_crctd_up = k_crctd + k_unc, k_crctd_down = k_crctd - k_unc;
            res(i, 0) = cvh_pt - (1.0 / k_crctd_down);
            res(i, 1) = (1.0 / k_crctd_up) - cvh_pt;
        }
        return res;
    }

    // for smearing weights derived from qop
    out_tensor_t operator() (
        double genQop, double genPhi, int genCharge, double genEta, double genPt, //for GEN params
        double qop, double phi, int cvhCharge,
        float cvhEta, float cvhPt, float jpsiCrctdPt, // for RECO params
        const RVec<float> &cov, // for sigma on the Gaussian
        out_tensor_t &deltaPts, // for the variations
        bool abQop = false, bool fullParam = false
    ) {
        const double cvhLam = calculateLam(cvhEta);
        const double cvhQopAbPt = calculateQop(cvhPt, cvhEta, cvhCharge);
        const Eigen::Vector3d parms(
            (abQop? qop : cvhQopAbPt),
            (fullParam? cvhLam : 0),
            (fullParam? phi : 0)
        );
    
        const double genLam = calculateLam(genEta);
        const double genQopAbPt = calculateQop(genPt, genEta, genCharge);
        const Eigen::Vector3d genparms(
            (abQop? genQop : genQopAbPt),
            (fullParam? genLam : 0),
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
    
        for (std::ptrdiff_t ivar = 0; ivar < nUnc; ++ivar) {

            double deltaQopDown = calculteDeltaQop(cvhPt, deltaPts(ivar, 0), cvhEta, cvhCharge);
            double deltaQopUp = calculteDeltaQop(cvhPt, deltaPts(ivar, 1), cvhEta, cvhCharge);
    
            for (std::ptrdiff_t idownup = 0; idownup < 2; ++idownup) {
                Eigen::Vector3d parmvar = Eigen::Vector3d::Zero();   
                const double dir = idownup == 0 ? -1. : 1.;
                parmvar[0] =  idownup == 0 ? deltaQopDown : deltaQopUp;
                const Eigen::Vector3d deltaparmsalt = deltaparms + dir*parmvar;
                const Eigen::Matrix<double, 3, 3> covdalt = covd; //+ dir*covvar;
    
                const Eigen::Matrix<double, 3, 3> covinvalt = covdalt.inverse();
                const double covdetalt = covdalt.determinant();
    
                const double lnpalt = -0.5*deltaparmsalt.transpose()*covinvalt*deltaparmsalt;
    
                const double weight = std::sqrt(covdet/covdetalt)*std::exp(lnpalt - lnp);
    
            // protect against outliers
            // if (weight > 0.9998 && weight < 1.0002) {cout << "smearing weight is " << weight << "covd is " << covd << std::endl;}
                res(ivar, idownup) = weight;
            }
        }
        return res;
    }
    */
private:
    std::shared_ptr<const T> correctionHist_;
};

}
