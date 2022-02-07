#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <ROOT/RVec.hxx>
#include "utils.h"

namespace wrem {

Eigen::TensorFixedSize<int, Eigen::Sizes<2>> prefsrLeptons(const wrem::EigenRVecView<int>& status,
        const wrem::EigenRVecView<int>& statusFlags, const wrem::EigenRVecView<int>& pdgId, const wrem::EigenRVecView<int>& motherIdx) {
    auto leptons = pdgId.abs() >= 11 && pdgId.abs() <= 16;
    auto status746 = status == 746;
    auto status23 = status == 23;
    auto motherV = pdgId(motherIdx) == 23 || pdgId(motherIdx).abs() == 24;
    auto fromHardProcess = statusFlags.unaryExpr([](int x) { return x & (1 << 8); }).cast<bool>();

    // TODO: Is there a way to relax the fromHardProcess condition?
    auto others = leptons && (motherV || status23) && fromHardProcess;
    
    // If there are status = 746 leptons, they came from photos and are pre-FSR
    // (but still need to check the mother in case photos was applied to other particles in the
    // event, and in case of radiation from taus and their decay products)
    auto photos = leptons && status746 && motherV;
    auto all = photos || others;

    Eigen::Array<Eigen::Index, 2, 1> selected;
    if (photos.count() == 2) {
        selected = wrem::make_nonzero(photos);
    }
    else if (all.count() == 2) {
        selected = wrem::make_nonzero(all);
    }
    else {
        throw std::range_error("Expected to find 2 pre-FSR leptons, but found " + std::to_string(all.count()) + ", " + std::to_string(photos.count()));
    }
    auto ids = pdgId(selected);
    bool partIdx = ids[0] > 0;
    Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
    out(0) = selected(partIdx);
    out(1) = selected(!partIdx);
    return out;
}

} 


