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

} 


