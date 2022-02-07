#ifndef WREMNANTS_THEORY_CORRECTIONS_H
#define WREMNANTS_THEORY_CORRECTIONS_H

#include <boost/histogram.hpp>
#include "utils.h"

namespace wrem {

template <typename T>
class ScetlibCorrectionsHelper {
public:
    ScetlibCorrectionsHelper(T&& corrections) :
        correctionHist_(std::make_shared<const T>(std::move(corrections))) {}

    const auto& operator() (float mass, float y, float pt) {
        size_t massbin = correctionHist_->template axis<0>().index(mass);
        size_t ybin = correctionHist_->template axis<1>().index(y);
        size_t ptbin = correctionHist_->template axis<2>().index(pt);
        return correctionHist_->at(massbin, ybin, ptbin);
    }
private:
    std::shared_ptr<const T> correctionHist_;
};

}

#endif
