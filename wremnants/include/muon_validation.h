#include <boost/histogram.hpp>
#include <stdlib.h>

namespace wrem {
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

    double operator() (float cvh_eta, float cvh_pt, int charge) {
        const auto &params = get_tensor(cvh_eta);
        const double A = params(0);
        const double e = params(1);
        const double M = params(2);
        double k = 1.0 / cvh_pt;
        double magnetic = 1.0 + A;
        double material = -1.0 * e * k;
        double alignment = charge * M;
        double k_crctd = (magnetic + material) * k + alignment;
        system("test");
        return (1.0 / k_crctd);
    }
private:
    std::shared_ptr<const T> correctionHist_;
};

}
