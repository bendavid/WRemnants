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
    auto photos = leptons && status746;
    auto all = photos || others;

    Eigen::Array<Eigen::Index, 2, 1> selected;
    if (photos.count() == 2) {
        selected = wrem::make_nonzero(photos);
    }
    else {
        selected = wrem::make_nonzero(all);
    }
    auto ids = pdgId(selected);
    bool partIdx = ids[0] > 0;
    Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
    out(0) = selected(partIdx);
    out(1) = selected(!partIdx);
    return out;
}

} 


