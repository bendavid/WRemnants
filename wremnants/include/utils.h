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

}

#endif
