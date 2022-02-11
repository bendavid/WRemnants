import time
time0 = time.time()

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=["WplusmunuPostVFP"])
parser.add_argument("--useBoost", type=bool, help="use boost histograms", default=False)
parser.add_argument("--useTensor", type=bool, help="use tensors", default=False)
parser.add_argument("--nHists", type=int, help="number of hist copies", default=1)
parser.add_argument("--useSingles", type=bool, help="use tensors", default=False)


args = parser.parse_args()

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
if not args.nThreads:
    ROOT.ROOT.EnableImplicitMT()
elif args.nThreads != 1:
    ROOT.ROOT.EnableImplicitMT(args.nThreads)

#ROOT.gErrorIgnoreLevel = ROOT.kPrint
#ROOT.ROOT.Experimental.RLogManager.Get().SetVerbosity(ROOT.ROOT.Experimental.ELogLevel.kInfo)

import pickle
import gzip

import narf
import wremnants
import hist
import lz4.frame
import logging
import numpy as np
import numba

filt = lambda x,filts=args.filterProcs: any([f in x.name for f in filts])
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles, filt=filt if args.filterProcs else None)

for dataset in datasets:
  print(dataset.name)


wprocs = ["WplusmunuPostVFP", "WminusmunuPostVFP", "WminustaunuPostVFP", "WplustaunuPostVFP"]
zprocs = ["ZmumuPostVFP", "ZtautauPostVFP"]

# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]
ptV_axis = hist.axis.Variable([0.0, 2.9, 4.7, 6.7, 9.0, 11.8, 15.3, 20.1, 27.2, 40.2, 13000.0], name="genPtV")

axis_pdfsyst_idx = hist.axis.Integer(0, 103, name = "pdfsyst_idx")

#assert(0)

ROOT.gInterpreter.Declare("""

namespace wtest {

using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_f = ROOT::VecOps::RVec<float>;

unsigned int count_fiducial(const Vec_i &pdgId, const Vec_i &status, const Vec_i &statusFlags, const Vec_f &pt, const Vec_f &eta) {

  unsigned int n_fiducial = 0;
  for (unsigned int i = 0; i < pdgId.size(); ++i) {

    if ( (abs(pdgId[i]) == 11 || abs(pdgId[i]) == 13) && status[i] == 1 && abs(eta[i])<2.4 && pt[i] > 25.) {
      ++n_fiducial;
    }
  }

  return n_fiducial;
}

}

template<typename V>
auto tensor_view(const V &vec, std::size_t start = 0) {
  return Eigen::TensorMap<const Eigen::Tensor<typename V::value_type, 1>>(vec.data() + start, vec.size() - start);
}

template <typename ArgType, typename = std::enable_if_t<std::is_same_v<typename ArgType::Scalar, bool>>>
auto tensor_count(const ArgType &arg) {
  return arg.template cast<std::size_t>().sum();
}

template <typename ArgType, typename = std::enable_if_t<std::is_same_v<typename ArgType::Scalar, bool>>>
std::size_t tensor_count_eval(const ArgType &arg) {
  return Eigen::TensorFixedSize<std::size_t, Eigen::Sizes<>>(tensor_count(arg))();
}



""")

@ROOT.Numba.Declare(["RVec<int>", "RVec<int>", "RVec<int>", "RVec<float>", "RVec<float>" ], "int")
def count_leptons(pdgId, status, statusFlags, pt, eta):
    abspdg11 = np.abs(pdgId) == 11
    abspdg13 = np.abs(pdgId) == 13
    return np.count_nonzero((abspdg11 | abspdg13) & status == 1 & (np.abs(eta)<2.4) & (pt > 25.))

@ROOT.Numba.Declare(["RVec<int>", "RVec<int>", "RVec<int>", "RVec<float>", "RVec<float>" ], "int")
def count_leptons_loop(pdgId, status, statusFlags, pt, eta):
    count = 0
    for i in range(len(pdgId)):
        abspdg11 = abs(pdgId[i]) == 11
        abspdg13 = abs(pdgId[i]) == 13
        if (abspdg11 or abspdg13) and status[i] == 1 and (abs(eta[i])<2.4) and (pt[i] > 25.):
            count += 1
    return count


c_sig = numba.types.intc(numba.types.CPointer(numba.types.intc),
                   numba.types.CPointer(numba.types.intc),
                   numba.types.CPointer(numba.types.intc),
                   numba.types.CPointer(numba.types.float32),
                   numba.types.CPointer(numba.types.float32),
                   numba.types.intc)

@numba.cfunc(c_sig)
def count_leptons_cfunc(pdgId_raw, status_raw, statusFlags_raw, pt_raw, eta_raw, nparts):
    pdgId = numba.carray(pdgId_raw, nparts)
    status = numba.carray(status_raw, nparts)
    statusFlags = numba.carray(statusFlags_raw, nparts)
    pt = numba.carray(pt_raw, nparts)
    eta = numba.carray(eta_raw, nparts)

    abspdg11 = np.abs(pdgId) == 11
    abspdg13 = np.abs(pdgId) == 13
    return np.count_nonzero((abspdg11 | abspdg13) & status == 1 & (np.abs(eta)<2.4) & (pt > 25.))

print(count_leptons_cfunc.address)

ROOT.gInterpreter.Declare(f"""

const std::ptrdiff_t count_leptons_cfunc_addr = {count_leptons_cfunc.address};

using func_ptr_t = int(*)(int*, int*, int*, float*, float*, int);

const func_ptr_t count_leptons_cfunc_ptr = reinterpret_cast<func_ptr_t>(count_leptons_cfunc_addr);

""")

ROOT.gInterpreter.Declare("""


using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_f = ROOT::VecOps::RVec<float>;

int count_leptons_raw(const Vec_i &pdgId, const Vec_i &status, const Vec_i &statusFlags, const Vec_f &pt, const Vec_f &eta) {

   return count_leptons_cfunc_ptr(const_cast<int*>(pdgId.data()), const_cast<int*>(status.data()), const_cast<int*>(statusFlags.data()), const_cast<float*>(pt.data()), const_cast<float*>(eta.data()), pdgId.size());

}

""")

#assert(0)


def build_graph(df, dataset):
    results = []

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")





    for i in range(1):
        #df = df.Define(f"fiducial_leptons_{i}", "(abs(GenPart_pdgId) == 11 || abs(GenPart_pdgId) == 13) && GenPart_status == 1 && (GenPart_statusFlags & 0x1) && abs(GenPart_eta)<2.4 && GenPart_pt > 25.")
        #df = df.Define(f"n_fiducial_leptons_{i}", f"Sum(fiducial_leptons_{i})")


        #df = df.Define(f"n_fiducial_leptons_{i}", "Sum((abs(GenPart_pdgId) == 11 || abs(GenPart_pdgId) == 13) && GenPart_status == 1 && abs(GenPart_eta)<2.4 && GenPart_pt > 25.)")
        #sum_fiducial = df.Sum(f"n_fiducial_leptons_{i}")



        #df = df.Define(f"n_fiducial_leptons_{i}", "int(tensor_count_eval((tensor_view(GenPart_pdgId).abs() == 11 || tensor_view(GenPart_pdgId).abs() == 13) && tensor_view(GenPart_status) == 1 && tensor_view(GenPart_statusFlags).unaryExpr([](int x) {return x & 0x1;}) && tensor_view(GenPart_eta).abs() < 2.f4 && tensor_view(GenPart_pt) > 25.f))")


        #ROOT.gInterpreter.Declare("""
        #template<typename V>
        #auto array_view(const V &vec) {
          #return Eigen::Map<const Eigen::Array<typename V::value_type, Eigen::Dynamic, 1>>(vec.data(), vec.size());
        #}
        #""")

        #df = df.Define(f"n_fiducial_leptons_{i}", "((array_view(GenPart_pdgId).abs() == 11 || array_view(GenPart_pdgId).abs() == 13) && array_view(GenPart_status) == 1 && array_view(GenPart_statusFlags).unaryExpr([](int x) {return x & 0x1;}) && array_view(GenPart_eta).abs() < float(2.4) && array_view(GenPart_pt) > float(25.)).count()")
        #sum_fiducial = df.Sum(f"n_fiducial_leptons_{i}")



        #df = df.Define(f"n_fiducial_leptons_{i}", "Numba::count_leptons(GenPart_pdgId, GenPart_status,GenPart_statusFlags, GenPart_pt, GenPart_eta)")
        df = df.Define(f"n_fiducial_leptons_{i}", "Numba::count_leptons_loop(GenPart_pdgId, GenPart_status,GenPart_statusFlags, GenPart_pt, GenPart_eta)")
        #df = df.Define(f"n_fiducial_leptons_{i}", "count_leptons_raw(GenPart_pdgId, GenPart_status,GenPart_statusFlags, GenPart_pt, GenPart_eta)")


        sum_fiducial = df.Sum(f"n_fiducial_leptons_{i}")
        results.append(sum_fiducial)


        #df = df.Define(f"n_fiducial_leptons_{i}", "((array_view(GenPart_pdgId).abs() == 11 || array_view(GenPart_pdgId).abs() == 13) && array_view(GenPart_status) == 1 &&  array_view(GenPart_eta).abs() < 2.4 && array_view(GenPart_pt) > 25.).count()")

        #df = df.Define(f"fiducial_leptons_{i}", "tensor_view(GenPart_pt) > 25.f")
        #df = df.Define(f"n_fiducial_leptons_{i}", f"tensor_count_eval(fiducial_leptons_{i})")

        #df = df.Define(f"n_fiducial_leptons_{i}", "wtest::count_fiducial(GenPart_pdgId, GenPart_status,GenPart_statusFlags, GenPart_pt, GenPart_eta)")

        #sum_fiducial = df.Sum(f"n_fiducial_leptons_{i}")
        #results.append(sum_fiducial)


    return results, weightsum

time_to_graph_setup = time.time()

resultdict = narf.build_and_run(datasets, build_graph)

print("init time:", time_to_graph_setup-time0)

