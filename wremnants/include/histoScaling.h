#ifndef WREMNANTS_HISTO_SCALING_H
#define WREMNANTS_HISTO_SCALING_H

#include "utils.h"

namespace wrem {

    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> dummyScaling(double nominal_weight, double scale) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
        outWeights.setConstant(nominal_weight);
            
        // Down weight, then up weight
        outWeights(0) = (1./scale) * nominal_weight;
        outWeights(1) = scale      * nominal_weight;
        return outWeights;
        
    }

}

#endif
