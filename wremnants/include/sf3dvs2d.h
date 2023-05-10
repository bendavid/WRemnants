#ifndef WREMNANTS_SF3DVS2D_H
#define WREMNANTS_SF3DVS2D_H

#include "utils.h"

namespace wrem {

    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> SF3DVS2D(double up, double down) {

        Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
            
        // Down weight, then up weight
        outWeights(0) = down;
        outWeights(1) = up;
        return outWeights;
        
    }

}

#endif
