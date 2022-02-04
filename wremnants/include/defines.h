#ifndef WREMNANTS_DEFINES_H
#define WREMNANTS_DEFINES_H

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

namespace wrem {

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

typedef enum {isoTrigPlus=0, isoTrigMinus, isoNotrig, noisoTrigPlus, noisoTrigMinus, noisoNotrig, antiisoTrigPlus, antiisoTrigMinus, antiisoNotrig, trackingReco, isoTrigPlusStepOfTnp, isoNotrigStepOfTnP, tracking, altTracking, reco, altReco, global, isoOnly, antiisoOnly, isoNotrigOnly, antiisoNotrigOnly} ScaleFactorType; // to use product
typedef enum {BToH=0, BToF, B, C, D, E, F, F_preVFP, F_postVFP, G, H, GToH} DataEra;  // keep sorted with postVFP after preVFP (BToH is an exception but we don't use it anyway, leave it as 0)
typedef enum {MC=0, Data} DataType;

bool isOddEvent(ULong64_t evt) {

  return (evt%2) ? 1 : 0;

}

bool isEvenEvent(ULong64_t evt) {

  return (evt%2) ? 0 : 1;

}

}

#endif
