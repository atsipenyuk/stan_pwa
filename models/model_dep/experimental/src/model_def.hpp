#ifndef SRC__STRUCTURES__SQR_HPP
#define SRC__STRUCTURES__SQR_HPP

#include "resonances/base_3d_res.hpp"
#include <src/typedefs.h>

template <typename T>
class Model {
public:
  Model(std::vector<Base_3d_res<T> > resonances) :
    resonances_(resonaces) {};

  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& x,
	   const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& y) {
    return this->resonances_[0].value(x,y)
  }

  std::vector<Base_3d_res<T> > resonances_;
}


#endif
