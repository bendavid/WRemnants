#ifndef WREMNANTS_HISTO_SCALING_H
#define WREMNANTS_HISTO_SCALING_H

#include "utils.h"

namespace wrem {

    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> constantScaling(double nominal_weight, double scale) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
            
        // Down weight, then up weight
        outWeights(0) = nominal_weight / scale;
        outWeights(1) = nominal_weight * scale;
        return outWeights;
        
    }

    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> twoPointScaling(double nominal_weight, double scaleDown, double scaleUp, double scaleCentral = 1.0) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
            
        // Down weight, then up weight, scaleCentral might generally be 1.0 in case it was not already included in nominal_weight for all histograms
        outWeights(0) = nominal_weight * scaleDown / scaleCentral;
        outWeights(1) = nominal_weight * scaleUp   / scaleCentral;
        return outWeights;
        
    }

}

#endif
