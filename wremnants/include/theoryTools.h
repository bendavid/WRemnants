#ifndef WREMNANTS_THEORYTOOLS_H
#define WREMNANTS_THEORYTOOLS_H

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

using scale_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>>;

scale_tensor_t makeScaleTensor(const Vec_f &scale_weights) {
  // from nanoaod doc lines (do I trust them?)
//[0] is mur=0.5 muf=0.5; [1] is mur=0.5 muf=1; [2] is mur=0.5 muf=2; [3] is mur=1 muf=0.5 ; [4] is mur=1 muf=1; [5] is mur=1 muf=2; [6] is mur=2 muf=0.5; [7] is mur=2 muf=1 ; [8] is mur=2 muf=2)*

  //ordering of the tensor axes are mur, muf with elements ordered 0.5, 1.0, 2.0

  scale_tensor_t res;
  res(0, 0) = scale_weights[0]; //mur=0.5 muf=0.5;
  res(0, 1) = scale_weights[1]; //mur=0.5 muf=1.0;
  res(0, 2) = scale_weights[2]; //mur=0.5 muf=2.0;
  res(1, 0) = scale_weights[3]; //mur=1.0 muf=0.5;
  res(1, 1) = scale_weights[4]; //mur=1.0 muf=1.0;
  res(1, 2) = scale_weights[5]; //mur=1.0 muf=2.0;
  res(2, 0) = scale_weights[6]; //mur=2.0 muf=0.5;
  res(2, 1) = scale_weights[7]; //mur=2.0 muf=1.0;
  res(2, 2) = scale_weights[8]; //mur=2.0 muf=2.0;

  return res;


}

using helicity_scale_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<9, 3, 3>>;

helicity_scale_tensor_t makeHelicityMomentScaleTensor(const CSVars &csvars, const scale_tensor_t &scale_tensor, double original_weight = 1.0) {

  constexpr Eigen::Index nhelicity = 9;
  constexpr Eigen::Index nmur = 3;
  constexpr Eigen::Index nmuf = 3;

  const double sinThetaCS = csvars.sintheta;
  const double cosThetaCS = csvars.costheta;
  const double sinPhiCS = csvars.sinphi;
  const double cosPhiCS = csvars.cosphi;

  const double sin2ThetaCS = 2.*sinThetaCS*cosThetaCS;
  const double sin2PhiCS = 2.*sinPhiCS*cosPhiCS;
  const double cos2ThetaCS = 1. - 2.*sinThetaCS*sinThetaCS;
  const double cos2PhiCS= 1. - 2.*sinPhiCS*sinPhiCS;

  // computing moments e.g. as used in arxiv:1708.00008 eq. 2.13
  Eigen::TensorFixedSize<double, Eigen::Sizes<nhelicity, 1, 1>> moments;
  moments(0, 0, 0) = 1.;
  moments(1, 0, 0) = cosThetaCS*cosThetaCS;
  moments(2, 0, 0) = sin2ThetaCS*cosPhiCS;
  moments(3, 0, 0) = sinThetaCS*sinThetaCS*cos2PhiCS;
  moments(4, 0, 0) = sinThetaCS*cosPhiCS;
  moments(5, 0, 0) = cosThetaCS;
  moments(6, 0, 0) = sinThetaCS*sinThetaCS*sin2PhiCS;
  moments(7, 0, 0) = sin2ThetaCS*sinPhiCS;
  moments(8, 0, 0) = sinThetaCS*sinPhiCS;

  constexpr std::array<Eigen::Index, 3> broadcastscales = { 1, nmur, nmuf };
  constexpr std::array<Eigen::Index, 3> broadcasthelicities = { nhelicity, 1, 1 };
  constexpr std::array<Eigen::Index, 3> reshapescale = { 1, nmur, nmuf };

  return original_weight*scale_tensor.reshape(reshapescale).broadcast(broadcasthelicities)*moments.broadcast(broadcastscales);


}

} 

#endif

