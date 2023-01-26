#include <boost/histogram.hpp>
#include <stdlib.h>

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
    float operator() (float cvh_pt, float cvh_eta, int charge) {
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

    private:
    std::shared_ptr<const T> correctionHist_;
};

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
