#ifndef SRC__STRUCTURES__BASE_3D_RES_HPP
#define SRC__STRUCTURES__BASE_3D_RES_HPP

#include <stan_pwa/src/exp_structures/particles_def.hpp>
#include <stan_pwa/src/typedefs.h>

namespace stan_pwa {

  class Base_3d_res {
  public:
    Base_3d_res(particle P, particle a, particle b, particle c) :
      P_(P), a_(a), b_(b), c_(c) {};

    ~Base_3d_res() {};

    template <typename T>
    C_t<promoted_t<T> > value(const T& x, const T& y);

    template <typename T>
    C_t<promoted_t<T> > value_sym(const T& x, const T& y);

    const particle P_; ///> Parent particle
    const particle a_; ///> Final state particles
    const particle b_; 
    const particle c_;
  };

}

#endif
