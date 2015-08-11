#ifndef SRC__STRUCTURES__MODEL_MODEL_WRAPPER_HPP
#define SRC__STRUCTURES__MODEL_MODEL_WRAPPER_HPP


//#include "model_inst.hpp"

namespace stan {
  namespace math {

  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& x,
	   const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& y) {
    return x+y;
  }

  }
}






#endif
