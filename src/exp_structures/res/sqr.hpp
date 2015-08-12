#ifndef SRC__STRUCTURES__SQR_HPP
#define SRC__STRUCTURES__SQR_HPP

#include "base_3d_res.hpp"

namespace stan_pwa {

  class Sqr : public Base_3d_res {
  public:
    Sqr(particle P, particle a, particle b, particle c) :
      Base_3d_res(P, a, b, c) {};
    ~Sqr() {};

    template <typename T>    
    C_t<promoted_t<T> > value(const T& x, const T& y) {
      C_t<T> res = {x*y, 0};
      return res;
    }

    template <typename T>
    C_t<promoted_t<T> > value_sym(const T& x, const T& y) {
      C_t<T> res = {x*y + y*x, 0};
      return res;
    }
  };

}
#endif
