#ifndef WREMNANTS_UTILS_H
#define WREMNANTS_UTILS_H


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>


namespace wrem {


template<std::ptrdiff_t N, typename V>
auto vec_to_tensor(const V &vec, std::size_t start = 0) {
  Eigen::TensorFixedSize<typename V::value_type, Eigen::Sizes<N>> res;
  std::copy(vec.begin() + start, vec.begin() + start + N, res.data());
  return res;
}

template<typename T, std::ptrdiff_t N, typename V>
auto vec_to_tensor_t(const V &vec, std::size_t start = 0) {
  Eigen::TensorFixedSize<T, Eigen::Sizes<N>> res;
  std::copy(vec.begin() + start, vec.begin() + start + N, res.data());
  return res;
}

template<typename T>
T clip_tensor(const T &tensor, const typename T::Scalar &thres) {
  return tensor.cwiseMax(-thres).cwiseMin(thres);
}

// like std::make_shared but detach the TH1-derived object from the current directory
template< class T, class... Args >
std::shared_ptr<T> make_shared_TH1( Args&&... args ) {
  using hist_t = std::decay_t<T>;
  hist_t *hist = new hist_t(std::forward<Args>(args)...);
  hist->SetDirectory(nullptr);
  return std::shared_ptr<T>(hist);
}

template <typename T>
class EigenRVecView : public Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>> {

private:
  using base_t = Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>;

public:
  EigenRVecView(const ROOT::VecOps::RVec<T> &vec) : base_t(vec.data(), vec.size()) {}

};

template<class ArgType>
struct nonzero_helper {
  using ArrayType =  Eigen::Array<typename Eigen::Index,
                 Eigen::Dynamic,
                 1,
                 Eigen::ColMajor,
                 ArgType::MaxSizeAtCompileTime,
                 1>;
};

template<class ArgType>
class nonzero_functor {
  const ArgType &m_vec;
public:
  using ArrayType = typename nonzero_helper<ArgType>::ArrayType;

  nonzero_functor(const ArgType& arg) : m_vec(arg) {}

  typename Eigen::Index operator() (Eigen::Index row) const {
    const bool cached = lastrow_ == (row - 1);
    lastrow_ = row;
    if (cached) {
      for (Eigen::Index i = lastidx_ + 1; i < m_vec.rows(); ++i) {
        if (m_vec[i] != 0) {
          lastidx_ = i;
          return i;
        }
      }
    }
    else {
      for (Eigen::Index i = 0, count = 0; i < m_vec.rows(); ++i) {
        if (m_vec[i] != 0) {
          if (count++ == row) {
            lastidx_ = i;
            return i;
          }
        }
      }
    }
    return -1;
  }

private:
  mutable Eigen::Index lastrow_ = -1;
  mutable Eigen::Index lastidx_ = -1;
};

template <class ArgType>
Eigen:: CwiseNullaryOp<nonzero_functor<ArgType>, typename nonzero_functor<ArgType>::ArrayType>
make_nonzero(const Eigen::ArrayBase<ArgType>& arg)
{
  using ArrayType = typename nonzero_helper<ArgType>::ArrayType;
  static_assert(ArrayType::ColsAtCompileTime == 1);
  std::size_t size;
  if constexpr (std::is_same_v<typename ArrayType::Scalar, bool>) {
    size = arg.count();
  }
  else {
    size = (arg != 0).count();
  }
  return ArrayType::NullaryExpr(size, 1, nonzero_functor<ArgType>(arg.derived()));
}

}

#endif
