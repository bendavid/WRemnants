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

template<typename V>
auto array_view(const V &vec, std::size_t start = 0) {
  return Eigen::Map<const Eigen::Array<typename V::value_type, Eigen::Dynamic, 1>>(vec.data() + start, vec.size() - start);
}

template<typename V>
auto array_view(V &vec, std::size_t start = 0) {
  return Eigen::Map<Eigen::Array<typename V::value_type, Eigen::Dynamic, 1>>(vec.data() + start, vec.size() - start);
}

template<typename V>
auto tensor_view(const V &vec, std::size_t start = 0) {
  return Eigen::TensorMap<const Eigen::Tensor<typename V::value_type, 1>>(vec.data() + start, vec.size() - start);
}

template<typename V>
auto tensor_view(V &vec, std::size_t start = 0) {
  return Eigen::TensorMap<Eigen::Tensor<typename V::value_type, 1>>(vec.data() + start, vec.size() - start);
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

template <typename ArgType, typename = std::enable_if_t<std::is_same_v<typename ArgType::Scalar, bool>>>
auto tensor_count(const ArgType &arg) {
  return arg.template cast<std::size_t>().sum();
}

template <typename ArgType, typename = std::enable_if_t<std::is_same_v<typename ArgType::Scalar, bool>>>
std::size_t tensor_count_eval(const ArgType &arg) {
  return Eigen::TensorFixedSize<std::size_t, Eigen::Sizes<>>(tensor_count(arg))();
}


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

  nonzero_functor(const ArgType& arg) : m_vec(arg) {}

  typename Eigen::Index operator() (Eigen::Index row) const {
    if (row == lastrow_) {
      return lastidx_;
    }
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

template<class ArgType>
class nonzero_tensor_functor {
public:

  nonzero_tensor_functor(const ArgType& arg) : arg_(arg) {}

  typename Eigen::Index operator() (Eigen::Index row) const {
    if (row == lastrow_) {
      return lastidx_;
    }
    const bool cached = lastrow_ == (row - 1);
    lastrow_ = row;
    if (cached) {
      for (Eigen::Index i = lastidx_ + 1; i < arg_.size(); ++i) {
        if (arg_(i)) {
          lastidx_ = i;
          return i;
        }
      }
    }
    else {
      for (Eigen::Index i = 0, count = 0; i < arg_.size(); ++i) {
        if (arg_(i)) {
          if (count++ == row) {
            lastidx_ = i;
            return i;
          }
        }
      }
    }
    lastidx_ = -1;
    return -1;
  }

private:
  const Eigen::TensorRef<Eigen::Tensor<typename ArgType::Scalar, 1>> arg_;
  mutable Eigen::Index lastrow_ = -1;
  mutable Eigen::Index lastidx_ = -1;
};

template <class ArgType>
Eigen:: CwiseNullaryOp<nonzero_functor<ArgType>, typename nonzero_helper<ArgType>::ArrayType>
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

template <typename ArgType>
auto make_nonzero_tensor(const ArgType &arg) {
  if constexpr (std::is_same_v<typename ArgType::Scalar, bool>) {
    auto asints = arg.template cast<Eigen::Index>();
    auto size = asints.sum();
    Eigen::TensorFixedSize<Eigen::Index, Eigen::Sizes<>> sizetensor = size;
    const Eigen::array<Eigen::Index, 1> offsets = {0};
    const Eigen::array<Eigen::Index, 1> extents = {sizetensor()};
    auto slice = asints.slice(offsets, extents);
    return slice.nullaryExpr(nonzero_tensor_functor(arg));
  }
  else {
    auto notzero = arg != static_cast<typename ArgType::Scalar>(0);
    auto asints = notzero.template cast<Eigen::Index>();
    auto size = asints.sum();
    Eigen::TensorFixedSize<Eigen::Index, Eigen::Sizes<>> sizetensor = size;
    const Eigen::array<Eigen::Index, 1> offsets = {0};
    const Eigen::array<Eigen::Index, 1> extents = {sizetensor()};
    auto slice = asints.slice(offsets, extents);
    return slice.nullaryExpr(nonzero_tensor_functor(notzero));
  }
}

template<class ArgType, class IndexType>
class fancy_index_tensor_functor {
public:
  fancy_index_tensor_functor(const ArgType &arg, const IndexType &idxs) : arg_(arg), idxs_(idxs) {}

  typename Eigen::Index operator() (Eigen::Index row) const {
    return arg_(idxs_(row));
  }


private:
  const Eigen::TensorRef<Eigen::Tensor<typename ArgType::Scalar, 1>> arg_;
  const Eigen::TensorRef<Eigen::Tensor<typename IndexType::Scalar, 1>> idxs_;
};

template<class ArgType, class IndexType>
auto fancy_index(const ArgType &arg, const IndexType &idxs) {
  return idxs.template cast<typename ArgType::Scalar>().nullaryExpr(fancy_index_tensor_functor(arg, idxs));
}

template<class ArgType, class MaskType>
auto bool_index(const ArgType &arg, const MaskType &mask) {
  return fancy_index(arg, make_nonzero_tensor(mask));
//   return make_nonzero_tensor(mask);
}

template<typename ArgTypeIf, typename ArgTypeThen, typename ArgTypeElse>
auto scalar_select(const ArgTypeIf &cond, const ArgTypeThen &arg0, const ArgTypeElse &arg1) {
  Eigen::TensorRef<Eigen::Tensor<typename ArgTypeThen::Scalar, 1>> arg0ref(arg0);
  Eigen::array<Eigen::Index, 1> shape{1};
  Eigen::array<Eigen::Index, 1> broadcast{arg0ref.size()};
  return cond.reshape(shape).broadcast(broadcast).select(arg0, arg1);
}

}


#endif
