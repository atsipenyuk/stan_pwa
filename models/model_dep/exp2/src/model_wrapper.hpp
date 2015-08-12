#ifndef SRC__STRUCTURES__MODEL_MODEL_WRAPPER_HPP
#define SRC__STRUCTURES__MODEL_MODEL_WRAPPER_HPP

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/model_inst.hpp>
#include <stan_pwa/src/typedefs.h>

namespace stan {
  namespace math {

  template <typename T>
  inline
  std::vector<typename boost::math::tools::promote_arg<T>::type>
  amplitude(const unsigned int &res_id, const Var_t<T>& y) {
//    return (stan_pwa::My_model.*stan_pwa::amp_ptr<T>)(res_id, y);
return stan_pwa::My_model.amplitude(res_id, y);
  }


    template <typename T0>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) 
    {
//      return (stan_pwa::My_model.*stan_pwa::amp_vector_ptr<T0>)(y);
      return stan_pwa::My_model.amplitude_vector<T0>(y);
    }


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const CV_t<T0>& x, const CV_t<T1>& y) {
    return stan_pwa::My_model.f_genfit(x,y);
  };


    
    inline int num_resonances() {
      return stan_pwa::My_model.get_num_res();
    }


    inline int num_variables() {
      return stan_pwa::My_model.get_num_var();
    }

  }
}






#endif


