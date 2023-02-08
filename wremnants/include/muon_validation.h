#include <boost/histogram.hpp>
#include <stdlib.h>
#include <defines.h>

namespace wrem {

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

// jpsi correction central value for one muon
template <typename T>
class JpsiCorrectionsHelper {

public:
    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;
    static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
    static constexpr auto nUnc = sizes[sizes.size() - 1]; // 1 for central value
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
    float operator() (float cvhPt, float cvhEta, int charge) {
        const auto &params = get_tensor(cvhEta);
        const double A = params(0, 0);
        const double e = params(1, 0);
        const double M = params(2, 0);
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
        double genQop, float genPhi, float genEta,
        double recoQop, float recoPhi, int recoCharge, float recoEta, float recoPt,
        const RVec<float> &cov, // for sigma on the Gaussian
        double nominal_weight, bool fullParam = false
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
    
            // protect against outliers
            // if (weight > 0.9998 && weight < 1.0002) {cout << "smearing weight is " << weight << "covd is " << covd << std::endl;}
                res(ivar, idownup) = nominal_weight * weight;
            }
        }
        return res;
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

}
