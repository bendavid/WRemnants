#pragma once

#include <ROOT/RDataFrame.hxx>
#include <eigen3/Eigen/Dense>

#include "matrix_utils.hpp"

namespace wrem {

class GradHelper : public ROOT::Detail::RDF::RActionImpl<GradHelper> {

public:
  using grad_t = float;

  using Result_t = std::vector<double>;

  GradHelper(unsigned int nparms)
      : nparms_(nparms), grad_(std::make_shared<Result_t>()) {}
  //   GradHelper(unsigned int nparms, std::shared_ptr<Result_t> grad) :
  //   nparms_(nparms), grad_(grad) {} GradHelper(GradHelper && other) =
  //   default; GradHelper(GradHelper && other) : nparms_(other.nparms_),
  //   grad_(other.grad_) {} GradHelper(const GradHelper &other) :
  //   nparms_(other.nparms_), grad_(other.grad_) {}

  std::shared_ptr<Result_t> GetResultPtr() const { return grad_; }

  void Exec(unsigned int slot, ROOT::VecOps::RVec<grad_t> const &vec,
            ROOT::VecOps::RVec<unsigned int> const &idxs) {
    Exec(slot, vec, idxs, 1.);
  }

  void Exec(unsigned int slot, ROOT::VecOps::RVec<grad_t> const &vec,
            ROOT::VecOps::RVec<unsigned int> const &idxs, double w) {
    std::vector<double> &grad = gradtmp_[slot];

    for (unsigned int i = 0; i < vec.size(); ++i) {
      const unsigned int &idx = idxs[i];
      grad[idx] += w * vec[i];
    }
  }
  void InitTask(TTreeReader *, unsigned int) {}

  void Initialize() {
    //     const unsigned int nslots = ROOT::IsImplicitMTEnabled() ?
    //     ROOT::GetImplicitMTPoolSize() : 1;
    const unsigned int nslots =
        ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
    gradtmp_.clear();
    gradtmp_.resize(nslots, std::vector<double>(nparms_, 0.));

    if (grad_->empty()) {
      grad_->clear();
      grad_->resize(nparms_, 0.);
    }
  }

  void Finalize() {
    for (auto const &grad : gradtmp_) {
      for (unsigned int i = 0; i < grad_->size(); ++i) {
        (*grad_)[i] += grad[i];
      }
    }
  }

  std::string GetActionName() { return "GradHelper"; }

private:
  unsigned int nparms_;
  std::vector<std::vector<double>> gradtmp_;
  std::shared_ptr<Result_t> grad_;
};

class HessHelper : public ROOT::Detail::RDF::RActionImpl<HessHelper> {

public:
  using grad_t = float;

  using Result_t = SymMatrixAtomic;
  //   using Result_t = std::vector<std::atomic<double> >;
  //   using Result_t = std::vector<double>;
  //   using Data_t = std::vector<std::atomic<double> >;

  HessHelper(unsigned int nparms) : grad_(std::make_shared<Result_t>(nparms)) {}

  //   HessHelper(unsigned int nparms) : nparms_(nparms),
  //   grad_(std::make_shared<Result_t>()) {} HessHelper(unsigned int nparms,
  //   std::shared_ptr<Result_t> grad) : nparms_(nparms), grad_(grad) {}
  //   HessHelper(HessHelper && other) = default;
  //   HessHelper(HessHelper && other) : nparms_(other.nparms_),
  //   grad_(other.grad_) {} HessHelper(const HessHelper &other) :
  //   nparms_(other.nparms_), grad_(other.grad_) {}

  std::shared_ptr<Result_t> GetResultPtr() const { return grad_; }
  //   std::shared_ptr<Result_t> GetResultPtr() const {
  //     return
  //     std::shared_ptr<Result_t>(reinterpret_cast<Result_t*>(gradatom_.get()));
  //   }

  void Exec(unsigned int slot, ROOT::VecOps::RVec<grad_t> const &vec,
            ROOT::VecOps::RVec<unsigned int> const &idxs) {
    Exec(slot, vec, idxs, 1.);
  }

  void Exec(unsigned int slot, ROOT::VecOps::RVec<grad_t> const &vec,
            ROOT::VecOps::RVec<unsigned int> const &idxs, double w) {

    unsigned int k = 0;
    for (unsigned int i = 0; i < idxs.size(); ++i) {
      const unsigned int iidx = idxs[i];
      for (unsigned int j = i; j < idxs.size(); ++j) {
        const unsigned int jidx = idxs[j];
        //         const double val = vec[k];

        const double val =
            (iidx == jidx && i != j) ? 2. * w * vec[k] : w * vec[k];

        grad_->fetch_add(iidx, jidx, val);

        ++k;
      }
    }
  }
  void InitTask(TTreeReader *, unsigned int) {}

  void Initialize() { timestamp_ = std::chrono::steady_clock::now(); }

  void Finalize() {
    //     std::vector<double>& tmp =
    //     *reinterpret_cast<std::vector<double>*>(&gradatom_);
    //     grad_->swap(tmp);

    //     double *data = reinterpret_cast<double*>(gradatom_.get());
    //     double *data = reinterpret_cast<double*>(gradatom_);
    //     std::vector<double> tmp(std::move(data),
    //     std::move(data+nparms_*nparms_)); grad_->swap(tmp); grad_.reset(new
    //     std::vector<double>(std::move(data),
    //     std::move(data+nparms_*nparms_)));

    //     std::cout << "val0 = " << (*grad_)[0] << std::endl;
    //     std::cout << "val1 = " << (*grad_)[1] << std::endl;

    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();

    auto timediff =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - timestamp_);

    std::cout << "Elapsed time = " << timediff.count() << std::endl;
  }

  std::string GetActionName() { return "HessHelper"; }

private:
  //   unsigned long long nparms_;
  std::shared_ptr<Result_t> grad_;
  //   std::unique_ptr<std::atomic<double>[]> gradatom_;
  //   std::atomic<double>* gradatom_;
  //   std::vector<std::vector<unsigned long long> > tmpidxs_;
  //   std::vector<unsigned long long> offsets_;
  std::chrono::steady_clock::time_point timestamp_;
};

} // namespace wrem