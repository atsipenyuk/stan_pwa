#ifndef PWA_STAN__SRC__MODEL_HPP
#define PWA_STAN__SRC__MODEL_HPP

#include "model_def.hpp"

namespace mc = stan_pwa::complex;

namespace stan_pwa {

  template <typename T>
  C_t<T> Model::amplitude(unsigned int i, const Var_t<T>& y) {
    return this->amplitudes_[i].value(y);
  };


  template <typename T>
  CV_t<T> Model::amplitude_vector(const Var_t<T>& y) {
    CV_t<T> res(2, (Eigen::Matrix<T,Eigen::Dynamic,1> (this->num_res_)));
    for (unsigned int i = 0; i < this->num_res_; i++) {
      res[i] = this->amplitude(i,y);
    }
    return res;
  };


  template <typename T>
  C_t<T> Model::amplitude_sym(unsigned int i, const Var_t<T>& y) {
    return this->amplitudes_[i].value_sym(y);
  };


  template <typename T>
  CV_t<T> Model::amplitude_vector_sym(const Var_t<T>& y) {
    CV_t<T> res(2, (Eigen::Matrix<T,Eigen::Dynamic,1> (this->num_res_)));
    for (int i = 0; i < this->num_res_; i++) {
      res[i] = this_amplitude_sym(i,y);
    }
    return res;
  };


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar
  //Model::f_genfit(const CV_t<T0>& amp, const CV_t<T1>& theta) {
  Model::f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& A_r,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta) {

      return mc::scalar::abs2(
               mc::vector::sum(
                 mc::vector::mult(A_r, theta)));
  };


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type
  Model::norm(const CV_t<T0>& theta, const CV_t<T0>& I) {
      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < NUM_RES; i++) {
	for (int j = 0; j < NUM_RES; j++) {
	  // Complex multiplication
	  tmp[0] = I[0](i,j) * theta[0](j) - I[1](i,j) * theta[1](j);
	  tmp[1] = I[0](i,j) * theta[1](j) + I[1](i,j) * theta[0](j);
          // Keep only the real part of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
	  res = res + theta[0](i) * tmp[0] + theta[1](i) * tmp[1];
        }
      }

      return res;
  };

} // end of pwa_stan

#endif









