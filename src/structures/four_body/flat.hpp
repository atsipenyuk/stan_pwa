#ifndef STAN_PWA__SRC__STRUCTURES__FOUR_BODY__FLAT_HPP
#define STAN_PWA__SRC__STRUCTURES__FOUR_BODY__FLAT_HPP

/*#include <cmath> // sqrt

#include <stan_pwa/lib/include/fct.hpp>
#include <stan_pwa/lib/include/complex.hpp>
#include "base.hpp" // base class

namespace resonances {

// deprecated???

// Non-resonant (flat) resonance, 5-dimensional 4 body decay.
struct flat_4 : resonances::resonance_base_4
{
  // Constructor
  flat_4(particle _P,
      particle _a, particle _b, particle _c, particle _d) :
        resonance_base_4(_P, _a, _b, _c, _d) {};

  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  std::vector<typename
  boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type>
  value(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
      const T3 &m2_34, const T4& m2_13)
  {

    typedef typename
        boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type T_res;


    std::vector<T_res> res(2, 0.0);

    return res;
  }
};

}*/

#endif
