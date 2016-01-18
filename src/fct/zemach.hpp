#ifndef STAN_PWA__SRC__FCT__ZEMACH_HPP
#define STAN_PWA__SRC__FCT__ZEMACH_HPP

#include <vector>
#include <cmath>
#include <complex>

#include <stan_pwa/src/fct.hpp>
namespace mfct = stan_pwa::fct;

namespace stan_pwa {
namespace fct {

  // Zemach - as described in 'The Physics of the B Factories', 
  // Chapter 13, p. 151.
  template <typename T0, typename T1, typename T2>
  typename boost::math::tools::promote_args<T0,T1,T2>::type
  zemach(int J, const T0& m2_ab, const T1& m2_bc, const T2& m2_R,
         const particle &a, const particle &b, const particle &c) {

    if (J == 0) return 1;
    if (J == 1) {
      return m2_R + a.m2 + b.m2 + c.m2 - m2_ab - 2 * m2_bc - 
          (m2_R - c.m2) * (a.m2 - b.m2) / m2_ab;
    }

    if (J == 2) {
      return mfct::zemach(1, m2_ab, m2_bc, m2_R, a, b, c) -
	(m2_ab - 2.*m2_R - 2.*c.m2 + (m2_R - c.m2) * (m2_R - c.m2) / m2_ab) *
	(m2_ab - 2.*a.m2 - 2.*b.m2 + (a.m2 - b.m2) * (a.m2 - b.m2) / m2_ab) /
	3.0;
    }

    return 0;
  }


  // Zemach - as described in 'Covariant spin tensors
  // in meson spectroscopy', Filippini, Fontana, Rotondi,
  // Phys. Rev. D, Vol. 51, Nr. 5, 2247-2261, published 1995.
  // For different channels J -> j + l the functions are
  // called zemach(J,j,l,...), respectively.
  // Note: to see, how the invariant mass variables m2_ij
  // are converted to z2, cos2_theta, see, e.g., 
  // src/struct/four_body/struct_four_particle_decay_channel.hpp
  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type
  zemach(const int J, const int j, const int l, 
	 const T0& z2, const T1& cos2_theta) {

    typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
    T_res res(0);

    switch (J) {
    case 0:
      // 0 -> 0 + 0
      if (j == 0 && l == 0) {
        res = 1.;
        break;
      }
      
      // 0 -> 1 + 1
      if (j == 1 && l == 1) {
        res = (1.0 + z2) * cos2_theta;
        break;
      }


    case 1:
      // 1 -> 1 + 0
      if (j == 1 && l == 0) {
        res = 1.0 + z2 * cos2_theta;
        break;
      }

      // 1 -> 1 + 1
      if (j == 1 && l == 1) {
        res = 1.0 - cos2_theta;
        break;
      }

      // 1 -> 1 + 2
      if (j == 1 && l == 2) {
        res = 1.0 + (3. + 4. * z2) * cos2_theta;
        break;
      }

      // 1 -> 2 + 1
      if (j == 2 && l == 1) {
        res = (1. + z2) * (1. + 3.*cos2_theta + 9.*z2 * pow(cos2_theta - 1./3., 2) );
        break;
      }


    case 2:
      // 2 -> 1 + 1
      if (j == 1 && l == 1) {
        res = 3. + (1. + 4.*z2) * cos2_theta;
        break;
      }

      // 2 -> 1 + 2
      if (j == 1 && l == 2) {
        res = 1. - cos2_theta;
        break;
      }

      // 2 -> 2 + 1
      if (j == 2 && l == 1) {
        res = 1. + z2/9. + (z2/3. - 1.) * cos2_theta - z2 * pow(cos2_theta - 1./3., 2);
        break;
      }

      // 2 -> 2 + 0
      if (j == 2 && l == 0) {
        res = 1. + z2/3. + z2 * cos2_theta + z2*z2 * pow(cos2_theta - 1./3., 2);
        break;
      }

    default:
      std::cerr << "fct::zemach - combination of J, j and l not handled. Return 0." << std::endl;
      res = 0.;
    }    

   
    return res;
  }
}
}
#endif
