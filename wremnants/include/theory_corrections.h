#ifndef WREMNANTS_THEORY_CORRECTIONS_H
#define WREMNANTS_THEORY_CORRECTIONS_H

#include <boost/histogram.hpp>
#include "utils.h"
#include <cmath>
#include "traits.h"
#include "csVariables.h"
#include "theoryTools.h"

namespace wrem {

const size_t NHELICITY = 9;
typedef Eigen::TensorFixedSize<double, Eigen::Sizes<NHELICITY>> helicity_tensor;

helicity_tensor csAngularFactors(const CSVars& csvars) {
    const double sinThetaCS = csvars.sintheta;
    const double cosThetaCS = csvars.costheta;
    const double sinPhiCS = csvars.sinphi;
    const double cosPhiCS = csvars.cosphi;

    const double sin2ThetaCS = 2.*sinThetaCS*cosThetaCS;
    const double sin2PhiCS = 2.*sinPhiCS*cosPhiCS;
    const double cos2ThetaCS = 1. - 2.*sinThetaCS*sinThetaCS;
    const double cos2PhiCS= 1. - 2.*sinPhiCS*sinPhiCS;
    helicity_tensor angular;
    angular(0) = 1.+cosThetaCS*cosThetaCS;
    angular(1) = 0.5*(1. - 3.*cosThetaCS*cosThetaCS);
    angular(2) = sin2ThetaCS*cosPhiCS;
    angular(3) = 0.5*sinThetaCS*sinThetaCS*cos2PhiCS;
    angular(4) = sinThetaCS*cosPhiCS;
    angular(5) = cosThetaCS;
    angular(6) = sinThetaCS*sinThetaCS*sin2PhiCS;
    angular(7) = sin2ThetaCS*sinPhiCS;
    angular(8) = sinThetaCS*sinPhiCS;
    return angular;
}

template <typename T>
class TensorCorrectionsHelper {

public:

    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;

    TensorCorrectionsHelper(T&& corrections) :
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

    //tensor_t operator() (double x, double y, double z, double nominal_weight = 1.0) {
    //    return nominal_weight*get_tensor(x, y, z);
    //}

    tensor_t operator() (double x1, double x2, double x3, int x4, double nominal_weight) {
        return nominal_weight*get_tensor(x1, x2, x3, x4);
    }
private:
    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class QCDScaleByHelicityCorrectionsHelper : public TensorCorrectionsHelper<T> {

using base_t = TensorCorrectionsHelper<T>;

using tensor_t = typename T::storage_type::value_type::tensor_t;
static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
static constexpr auto nhelicity = sizes[0];
static constexpr auto nmur = sizes[1];
static constexpr auto nmuf = sizes[2];


using small_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nhelicity, 1, 1>>;


public:

    //inherit constructor
    using base_t::base_t;

    // TODO: Surely this can be done better
    small_tensor_t csAngularFactorsTensor(const CSVars &csvars) {
        small_tensor_t angular;
        // n.b. setValues is error-prone for multi-dim tensors
        auto angular1D = csAngularFactors(csvars);
        for (size_t i = 0; i < NHELICITY; i++) {
            angular(i, 0, 0) = angular1D[i];
        }
        return angular;
    }

    tensor_t operator() (double mV, double yV, double ptV, int qV, const CSVars &csvars, const scale_tensor_t &scale_tensor, double nominal_weight) {
        // pure angular terms without angular coeffs multiplied through
        const auto angular = csAngularFactorsTensor(csvars);

        static_assert(sizes.size() == 3);
        static_assert(nhelicity == NHELICITY);

        constexpr std::array<Eigen::Index, 3> broadcastscales = { 1, nmur, nmuf };
        // now multiplied through by angular coefficients (1.0 for 1+cos^2theta term)
        const tensor_t angular_with_coeffs = angular.broadcast(broadcastscales)*base_t::get_tensor(mV, yV, ptV, qV);

        // denominator for each scale combination
        constexpr std::array<Eigen::Index, 1> helicitydims = { 0 };
        constexpr std::array<Eigen::Index, 3> broadcasthelicities = { nhelicity, 1, 1 };
        constexpr std::array<Eigen::Index, 3> reshapeden = { 1, nmur, nmuf };
        auto denominator = angular_with_coeffs.sum(helicitydims).reshape(reshapeden).broadcast(broadcasthelicities);

        constexpr std::array<Eigen::Index, 3> reshapescale = { 1, nmur, nmuf };
        auto scale = scale_tensor.reshape(reshapescale).broadcast(broadcasthelicities);

        return nominal_weight*scale*angular_with_coeffs/denominator;
    }

};

template <typename T>
class CentralCorrByHelicityHelper : public TensorCorrectionsHelper<T> {
using base_t = TensorCorrectionsHelper<T>;

using tensor_t = typename T::storage_type::value_type::tensor_t;
static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
static constexpr auto nhelicity = sizes[0];

public:
    using base_t::base_t;

    tensor_t operator() (double mV, double yV, double ptV, int qV, const CSVars &csvars, double nominal_weight) {
        const auto angular = csAngularFactors(csvars);
        const auto coeffs = base_t::get_tensor(mV, yV, ptV, qV);
        auto uncorr = coeffs.chip(0, 0)*angular;
        auto corr = coeffs.chip(0, 1)*angular;

        return corr.sum()/uncorr.sum();
    }
};

}

#endif
