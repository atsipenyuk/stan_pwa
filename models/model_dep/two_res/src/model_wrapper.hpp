#ifndef STAN_PWA__SRC__MODEL_WRAPPER_HPP
#define STAN_PWA__SRC__MODEL_WRAPPER_HPP

#include "model_def.hpp"
#include "model_inst.hpp"
// Wrap user's input in model_inst.hpp to Stan-usable form

namespace stan {
  namespace math {

    template <typename T>
    stan_pwa::Wrapper<T> ModelWrapper(stan_pwa::sym_flag,
				      stan_pwa::amp_ptr<T>,
				      stan_pwa::amp_vector_ptr<T>);

    template <typename T>
    inline
    std::vector<typename boost::math::tools::promote_arg<T>::type>
    amplitude(const unsigned int &res_id, 
	      const Eigen::Matrix<T, Eigen::Dynamic,1>& y) 
    {
      return (stan_pwa::MyModel.*stan_pwa::amp_ptr<T>)(res_id, y);
    }


    template <typename T0>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) 
    {
      return (stan_pwa::MyModel.*stan_pwa::amp_vector_ptr<T0>)(y);
    }


    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& A_r,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta) {
    
      /*    f_genfit(const CV_t<T0>& A_r, const CV_t<T1>& theta) {*/
      return stan_pwa::MyModel.f_genfit(A_r, theta);
      }


    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    norm(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& theta,
         const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> >& I) {
      return stan_pwa::MyModel.norm(theta, I);
    }


    inline int num_resonances() {
      return stan_pwa::MyModel.get_num_res();
    }


    inline int num_variables() {
      return stan_pwa::MyModel.get_num_var();
    }

  }
}
#endif
