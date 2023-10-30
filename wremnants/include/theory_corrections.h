#ifndef WREMNANTS_THEORY_CORRECTIONS_H
#define WREMNANTS_THEORY_CORRECTIONS_H

#include <boost/histogram.hpp>
#include "utils.h"
#include <cmath>
#include "traits.h"
#include "csVariables.h"
#include "theoryTools.h"

namespace wrem {

template <typename T>
class TensorCorrectionsHelper {

public:

    using hist_t = T;
    using tensor_t = typename T::storage_type::value_type::tensor_t;

    TensorCorrectionsHelper(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    template<typename... Xs>
    const tensor_t &get_tensor(const Xs&... xs) {
        return narf::get_value(*correctionHist_, xs...).data();
    }
    tensor_t operator() (double x1, double x2, double x3, int x4, double nominal_weight) {
        return nominal_weight*get_tensor(x1, x2, x3, x4);
    }


private:
    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class TensorCorrectionsHelper2D : public TensorCorrectionsHelper<T> {

using base_t = TensorCorrectionsHelper<T>;
using tensor_t = typename T::storage_type::value_type::tensor_t;

public:
    //inherit constructor
    using base_t::base_t;

    tensor_t operator() (double x1, int charge, double nominal_weight) {
        return nominal_weight*base_t::get_tensor(x1, charge);
    }
};

template <typename T>
class TensorCorrectionsHelper3D : public TensorCorrectionsHelper<T> {

using base_t = TensorCorrectionsHelper<T>;
using tensor_t = typename T::storage_type::value_type::tensor_t;

public:

    //inherit constructor
    using base_t::base_t;

    tensor_t operator() (double x1, double x2, int charge, double nominal_weight) {
        return nominal_weight*base_t::get_tensor(x1, x2, charge);
    }
};

template <typename T>
class TensorCorrectionsHelperWeighted4D : public TensorCorrectionsHelper<T> {

using base_t = TensorCorrectionsHelper<T>;
using tensor_t = typename T::storage_type::value_type::tensor_t;

public:

    //inherit constructor
    using base_t::base_t;

    tensor_t operator() (double x1, double x2, double x3, int charge, tensor_t nominal_weights) {
        return nominal_weights*base_t::get_tensor(x1, x2, x3, charge);
    }
};

template <typename T>
class QCDScaleByHelicityCorrectionsHelper : public TensorCorrectionsHelper<T> {

using base_t = TensorCorrectionsHelper<T>;

using tensor_t = typename T::storage_type::value_type::tensor_t;
static constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
// static constexpr auto nhelicity = sizes[0];
static constexpr auto nhelicity = 10;
static constexpr auto nmur = sizes[1];
static constexpr auto nmuf = sizes[2];


using small_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<nhelicity, 1, 1>>;


public:

    //inherit constructor
    using base_t::base_t;

    tensor_t operator() (double mV, double yV, double ptV, int qV, const CSVars &csvars, const scale_tensor_t &scale_tensor, double nominal_weight) {
        // denominator for each scale combination
        constexpr std::array<Eigen::Index, 1> helicitydims = { 0 };
        constexpr std::array<Eigen::Index, 3> broadcasthelicities = { nhelicity-1, 1, 1 };
        constexpr std::array<Eigen::Index, 3> reshapeden = { 1, nmur, nmuf };

        // pure angular terms without angular coeffs multiplied through
        const auto &angular = csAngularFactors(csvars);
        const auto &coeffs = base_t::get_tensor(mV, yV, ptV, qV);

        static_assert(sizes.size() == 3);
        static_assert(nhelicity == NHELICITY);

        constexpr std::array<Eigen::Index, 3> broadcastscales = { 1, nmur, nmuf };
        // now multiplied through by angular coefficients (1.0 for 1+cos^2theta term)
        tensor_t angular_with_coeffs;
        Eigen::TensorFixedSize<double, Eigen::Sizes<1, nmur, nmuf>> denominator;

        for(unsigned int imur = 0; imur < 3;imur++) {
            for(unsigned int imuf = 0; imuf < 3;imuf++) {
                angular_with_coeffs(0,imur,imuf) = coeffs(0,imur,imuf)*angular(9); // 1.0 times 1+cos^2theta for all scales
                denominator(0,imur,imuf) = coeffs(0,imur,imuf)*angular(9);
                for(unsigned int ihel = 1; ihel < NHELICITY-1;ihel++) {
                    angular_with_coeffs(ihel,imur,imuf) = coeffs(ihel,imur,imuf)*angular(ihel);
                    denominator(0,imur,imuf)+=coeffs(ihel,imur,imuf)*angular(ihel);
                }
            }
        }

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
static constexpr auto ncorrs = sizes[1];
static constexpr auto nvars = sizes[2];
const size_t NHELICITY = 10;

// TODO: Can presumably get the double type from the template param
typedef Eigen::TensorFixedSize<double, Eigen::Sizes<nvars>> var_tensor_t;

public:
    using base_t::base_t;

    var_tensor_t operator() (double mV, double yV, double ptV, int qV, const CSVars &csvars, double nominal_weight) {
        static_assert(sizes.size() == 3);
        static_assert(nhelicity == NHELICITY);
        static_assert(ncorrs == 2);

        const auto &angular = csAngularFactors(csvars);
        const auto &coeffs = base_t::get_tensor(mV, yV, ptV, qV);

        constexpr std::array<Eigen::Index, 3> reshapedims = {nhelicity, 1, 1};
        constexpr std::array<Eigen::Index, 1> reduceddims = {0};

        Eigen::TensorFixedSize<double, Eigen::Sizes<nhelicity-1, 1, 1>> coeffs_with_angular = coeffs*angular;

        coeffs_with_angular(0) = coeffs(0) * angular(9);
        for(unsigned int i = 1; i < NHELICITY-1;i++) {
            coeffs_with_angular(i) = coeffs(i) * angular(i);
        }

        auto uncorr_hel = coeffs_with_angular.chip(0, nhelicity-1).reshape(reshapedims).broadcast(sizes).chip(0, 1);
        auto corr_hel = coeffs_with_angular.chip(0, nhelicity-1).reshape(reshapedims).broadcast(sizes).chip(1, 1);

        var_tensor_t corr_weight_vars = corr_hel.sum(reduceddims)/uncorr_hel.sum(reduceddims)*nominal_weight;

        return corr_weight_vars;
    }
};

}

#endif
