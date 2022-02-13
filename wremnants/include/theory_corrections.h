#ifndef WREMNANTS_THEORY_CORRECTIONS_H
#define WREMNANTS_THEORY_CORRECTIONS_H

#include <boost/histogram.hpp>
#include "utils.h"
#include <cmath>
#include "traits.h"

namespace wrem {

template <typename T>
class TensorCorrectionsHelper3D {

typedef typename T::storage_type::value_type::tensor_t tensor_t;

public:
    TensorCorrectionsHelper3D(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    tensor_t operator() (double x, double y, double z) {
        auto xbin = correctionHist_->template axis<0>().index(x);
        auto ybin = correctionHist_->template axis<1>().index(y);
        auto zbin = correctionHist_->template axis<2>().index(z);
        return correctionHist_->at(xbin, ybin, zbin).data();
    }
private:
    std::shared_ptr<const T> correctionHist_;
};

template <typename T>
class QCDScaleByHelicityCorrectionsHelper : public TensorCorrectionsHelper3D<T> {

typedef typename T::storage_type::value_type::tensor_t tensor_t;

public:
    QCDScaleByHelicityCorrectionsHelper(T&& corrections) : TensorCorrectionsHelper3D<T>(std::move(corrections)) {}

    tensor_t csAngularFactors(const double sinThetaCS, const double cosThetaCS, 
                const double sinPhiCS, const double cosPhiCS) {
        constexpr auto sizes = narf::tensor_traits<tensor_t>::sizes;
        const double NORM = 3./16.*M_PI;
        const double sin2ThetaCS = 2.*sinThetaCS*cosThetaCS;
        const double cos2ThetaCS = 1-2*sinThetaCS*sinThetaCS;
        const double cos2PhiCS= 1-2*sinPhiCS*sinPhiCS;
        Eigen::TensorFixedSize<double, Eigen::Sizes<sizes[sizes.size()-1]>> angular;
        angular.setValues({(1.+cosThetaCS*cosThetaCS),
            0.5*(1.-3.*cosThetaCS*cosThetaCS),
            sin2ThetaCS*cosPhiCS, 
            0.5*sinThetaCS*sinThetaCS*cos2PhiCS,
            sinThetaCS*cosPhiCS,
            cosThetaCS,
            sinThetaCS*sinThetaCS*sin2ThetaCS,
            sin2ThetaCS*sinPhiCS,
            sinThetaCS*sinPhiCS
        });
        auto angscale = angular*NORM;
        // Eigen syntax for broadcasting is a bit ridiculous, first need to reshape
        // with new axes of size 1, then give the number of times to duplicate each axis
        std::array<Eigen::Index, sizes.size()> shape;
        shape.fill(1.);
        shape[sizes.size()-1] = sizes[sizes.size()-1];
        auto angresize = angular.reshape(shape);
        std::array<Eigen::Index, sizes.size()> broadcast = sizes;
        broadcast[sizes.size()-1] = 1;
        return angresize.broadcast(broadcast);
    }

    tensor_t operator() (double ptV, double yV, double qV, std::array<double, 4>& csSineCosThetaPhi) {
        auto angular = csAngularFactors(csSineCosThetaPhi[0], csSineCosThetaPhi[1], csSineCosThetaPhi[2], csSineCosThetaPhi[3]);
        return(angular*TensorCorrectionsHelper3D<T>::operator()(ptV, yV, qV));
    }

};

}

#endif
