#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <ROOT/RVec.hxx>
#include "utils.h"

namespace wrem {

Eigen::TensorFixedSize<int, Eigen::Sizes<2>> prefsrLeptons(const ROOT::VecOps::RVec<int>& status,
        const ROOT::VecOps::RVec<int>& statusFlags, const ROOT::VecOps::RVec<int>& pdgId, const ROOT::VecOps::RVec<int>& motherIdx) {

  const std::size_t ngenparts = status.size();


  std::array<std::size_t, 2> photos_idxs;
  std::array<std::size_t, 2> all_idxs;
  std::size_t nphotos = 0;
  std::size_t nall = 0;
  for (std::size_t i = 0; i < ngenparts; ++i) {
    const int &istatus = status[i];
    const int &istatusFlags = statusFlags[i];
    const int &ipdgId = pdgId[i];
    const int &imotherIdx = motherIdx[i];

    const int absPdgId = std::abs(ipdgId);

    const bool is_lepton = absPdgId >= 11 && absPdgId <= 16;
    const bool is_status746 = istatus == 746;
    const bool is_status23 = istatus == 23;
    const int &motherPdgId = pdgId[imotherIdx];
    const bool is_motherV = motherPdgId == 23 || std::abs(motherPdgId) == 24;

    const bool is_photos = is_lepton && is_status746 && is_motherV;


    const bool is_fromHardProcess = istatusFlags & ( 1 << 8 );

    // TODO: Is there a way to relax the fromHardProcess condition?
    const bool is_other = is_lepton && (is_motherV || is_status23) && is_fromHardProcess;

    // If there are status = 746 leptons, they came from photos and are pre-FSR
    // (but still need to check the mother in case photos was applied to other particles in the
    // event, and in case of radiation from taus and their decay products)
    const bool is_all = is_photos || is_other;

    if (is_photos) {
      if (nphotos < 2) {
        photos_idxs[nphotos] = i;
      }
      ++nphotos;
    }

    if (is_all) {
      if (nall < 2) {
        all_idxs[nall] = i;
      }
      ++nall;
    }
  }

  std::array<std::size_t, 2> selected_idxs;
  if (nphotos == 2) {
    selected_idxs = photos_idxs;
  }
  else if (nall == 2) {
    selected_idxs = all_idxs;
  }
  else {
    throw std::range_error("Expected to find 2 pre-FSR leptons, but found " + std::to_string(nall) + ", " + std::to_string(nphotos));
  }

  std::array<int, 2> selected_pdgids = { pdgId[selected_idxs[0]], pdgId[selected_idxs[1]] };
  const bool partIdx = selected_pdgids[0] > 0;
  Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
  out(0) = selected_idxs[partIdx];
  out(1) = selected_idxs[!partIdx];
  return out;

}

Eigen::TensorFixedSize<int, Eigen::Sizes<2>> prefsrLeptonsTensor(const ROOT::VecOps::RVec<int>& status_in,
        const ROOT::VecOps::RVec<int>& statusFlags_in, const ROOT::VecOps::RVec<int>& pdgId_in, const ROOT::VecOps::RVec<int>& motherIdx_in) {

    auto const status = tensor_view(status_in);
    auto const statusFlags = tensor_view(statusFlags_in);
    auto const pdgId = tensor_view(pdgId_in);
    auto const motherIdx = tensor_view(motherIdx_in);

    auto const leptons = pdgId.abs() >= 11 && pdgId.abs() <= 16;
    auto const status746 = status == 746;
    auto const status23 = status == 23;
    auto const motherPdgId = fancy_index(pdgId, motherIdx);
    auto const motherV = motherPdgId == 23 || motherPdgId.abs() == 24;
    auto const fromHardProcess = statusFlags.unaryExpr([](int x) { return x & (1 << 8); }).cast<bool>();

    // TODO: Is there a way to relax the fromHardProcess condition?
    auto const others = leptons && (motherV || status23) && fromHardProcess;

    // If there are status = 746 leptons, they came from photos and are pre-FSR
    // (but still need to check the mother in case photos was applied to other particles in the
    // event, and in case of radiation from taus and their decay products)
    auto const photos = leptons && status746 && motherV;
    auto const all = photos || others;

    const bool is_photos = tensor_count_eval(photos) == 2;

    Eigen::TensorFixedSize<Eigen::Index, Eigen::Sizes<2>> selected;
    if (tensor_count_eval(photos) == 2) {
      selected = make_nonzero_tensor(photos);
    }
    else if (tensor_count_eval(all) == 2) {
      selected = make_nonzero_tensor(all);
    }
    else {
      throw std::range_error("Expected to find 2 pre-FSR");
    }

    const Eigen::TensorFixedSize<int, Eigen::Sizes<2>> ids = fancy_index(pdgId, selected);

    const bool partIdx = ids[0] > 0;
    Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
    out(0) = selected(partIdx);
    out(1) = selected(!partIdx);
    return out;

//     auto const is_photos = tensor_count(photos) == std::size_t(2);
//
//     auto const prefsr = scalar_select(is_photos, photos, all);
//
//     if (tensor_count_eval(prefsr) != std::size_t(2)) {
//       throw std::range_error("Expected to find 2 pre-FSR leptons, but found " + std::to_string(tensor_count_eval(prefsr)));
//     }
//
//     const Eigen::TensorFixedSize<Eigen::Index, Eigen::Sizes<2>> selected = make_nonzero_tensor(prefsr);
//
//     const Eigen::TensorFixedSize<int, Eigen::Sizes<2>> ids = fancy_index(pdgId, selected);
//
//     const bool partIdx = ids[0] > 0;
//     Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
//     out(0) = selected(partIdx);
//     out(1) = selected(!partIdx);
//     return out;
}

Eigen::TensorFixedSize<int, Eigen::Sizes<2>> prefsrLeptonsArr(const wrem::EigenRVecView<int>& status,
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


