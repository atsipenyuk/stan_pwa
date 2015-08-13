#ifndef STAN_PWA__SRC__FCT__BLATT_WEISSKOPF_HPP
#define STAN_PWA__SRC__FCT__BLATT_WEISSKOPF_HPP

#include <vector>
#include <math.h>

#include <stan_pwa/src/fct/breakup_momentum.hpp>

namespace stan_pwa {
namespace fct {

  /**
   * Return floating-point Blatt-Weisskopf form factor.
   *
   * Implemented as in: arxiv:1406.6311v2, p. 151, eq. (13.2.8).
   *
   * @param J_R resonance spin
   * @param r2_p parent Particle squared radius
   * @param m2_ab Dalitz plot variable (squared mass)
   * @param m_a 1st daughter Particle mass of m2_ab
   * @param m_b 2nd daughter Particle mass of m2_ab
   * @return Blatt-Weisskopf form factor
   */
  template <typename T0, typename T1, typename T2>
  inline
  typename boost::math::tools::promote_args<T0,T1,T2>::type
  blatt_weisskopf(int J_R, double r2_P, 
                  const T0 &m2_ab, const T1& m_a, const T2 &m_b) {

    if (J_R == 0 or J_R > 2) return 1;

    typedef typename boost::math::tools::promote_args<T0,T1,T2>::type T;

    T p2 = fct::breakup_momentum::p2(m2_ab, m_a, m_b) * r2_P;
    if (J_R == 1) {
      return sqrt(1.0 / (1.0 + p2));
    }
    if (J_R == 2) {
      return sqrt(1.0 / (9.0 + 3.0 * p2 + p2 * p2));
    }

    return 0;
  }

}
}
#endif
