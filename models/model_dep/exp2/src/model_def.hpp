#ifndef SRC__MODEL_DEF_HPP
#define SRC__MODEL_DEF_HPP

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/typedefs.h>
#include <stan_pwa/src/exp_structures/res.hpp>

namespace mc = stan_pwa::complex;

namespace stan_pwa {

class Model {
public:
  Model(std::vector<Base_3d_res> resonances, unsigned int num_var) :
    resonances_(resonances),
    num_res_(resonances.size()),
    num_var_(num_var) {};

    ///> returns the n-th resonance value (complex number)
    template <typename T>
    //C_t<typename boost::math::tools::promote_arg<T>::type> 
    std::vector<typename boost::math::tools::promote_arg<T>::type> amplitude(unsigned int, const Var_t<T>&);

    ///> returns all resonance values
    template <typename T>
    std::vector< Eigen::Matrix<typename boost::math::tools::promote_arg<T>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Var_t<T>&);

  template <typename T0, typename T1>
  static
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const CV_t<T0>& x, const CV_t<T1>& y);

  ///> Calculates the normalization integral for the fitting
  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type
  norm(const CV_t<T0>&, const CV_t<T0>&);

  std::vector<Base_3d_res> resonances_;

  int get_num_res() {return this->num_res_;};
  int get_num_var() {return this->num_var_;};

private:
  int num_res_;
  int num_var_;
};


  // External cpp files here cause problems in the STAN linker
  // They may prabably be resolved adjusting CFLAGS in Stan makefile
  template <typename T>
  std::vector<typename boost::math::tools::promote_arg<T>::type>
  //  C_t<typename boost::math::tools::promote_arg<T>::type> 
  Model::amplitude(unsigned int i, const Var_t<T>& y) {
    return this->resonances_[i].value(y(0), y(1));
  };


  template <typename T>
  std::vector< Eigen::Matrix<typename boost::math::tools::promote_arg<T>::type, Eigen::Dynamic, 1> >
  //CV_t<typename boost::math::tools::promote_arg<T>::type> 
  Model::amplitude_vector(const Var_t<T>& y) {
    CV_t<typename boost::math::tools::promote_arg<T>::type> 
      res(2, (Eigen::Matrix<T,Eigen::Dynamic,1> (this->num_res_)));

    for (unsigned int i = 0; i < this->num_res_; i++) {
      res[i] = this->resonances_[i].value(y(0), y(1));
    }
    return res;
  };


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  Model::f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& x,
		  const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& y)
  {
    return mc::scalar::abs2(mc::vector::sum(mc::vector::mult(x,y)));
  }


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type
  Model::norm(const CV_t<T0>& theta, const CV_t<T0>& I) {
      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < this->num_res_; i++) {
	for (int j = 0; j < this->num_res_; j++) {
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

}
#endif

