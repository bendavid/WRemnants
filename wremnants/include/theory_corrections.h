#ifndef WREMNANTS_THEORY_CORRECTIONS_H
#define WREMNANTS_THEORY_CORRECTIONS_H

#include <boost/histogram.hpp>
#include "utils.h"

namespace wrem {

template <typename T>
class ScetlibCorrectionsHelper {

typedef typename T::storage_type::value_type::tensor_t tensor_t;

public:
    ScetlibCorrectionsHelper(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    tensor_t operator() (double mass, double y, double pt) {
        auto massbin = correctionHist_->template axis<0>().index(mass);
        auto ybin = correctionHist_->template axis<1>().index(y);
        auto ptbin = correctionHist_->template axis<2>().index(pt);
        return correctionHist_->at(massbin, ybin, ptbin).data();
    }
private:
    std::shared_ptr<const T> correctionHist_;
};

}

#endif
