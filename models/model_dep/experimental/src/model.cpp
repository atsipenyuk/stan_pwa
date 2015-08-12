/*#ifndef SRC__MODEL_CPP
#define SRC__MODEL_CPP

#include <stan_pwa/src/model_def.hpp>

namespace stan_pwa {

  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& x,
	   const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& y) {
          return mc::scalar::abs2(mc::vector::sum(mc::vector::mult(x,y)));
  };

}
#endif*/
