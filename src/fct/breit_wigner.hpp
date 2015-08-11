#ifndef STAN_PWA__SRC__FCT__BREIT_WIGNER_HPP
#define STAN_PWA__SRC__FCT__BREIT_WIGNER_HPP

#include <vector>
#include <cmath>

#include <stan_pwa/src/fct/breakup_momentum.hpp>
#include <stan_pwa/src/fct/blatt_weisskopf.hpp>
#include <stan_pwa/src/complex.hpp>
namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;

namespace stan_pwa {
namespace fct {
  namespace breit_wigner {

    /**
     * Return complex Breit-Wigner form factor.
     *
     * @param M_R resonance mass
     * @param m2_ab Dalitz plot variable (squared mass)
     * @param width_m2_ab resonance width
     * @return Breit-Wigner dynamical form factor
     */
    template <typename T0, typename T1, typename T2>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2>::type>
    value(const T0& M_R, const T1& m2_ab, const T2& width_m2_ab) {

      typedef typename boost::math::tools::promote_args<T0,T1,T2>::type T_res;

      std::vector<T_res> res(2);
      res = mc::scalar::complex(M_R * M_R - m2_ab, - M_R * width_m2_ab);
      return mc::scalar::inverse(res);

    }


    /**
     * Return Relativistic Breit Wigner resonance width.
     *
     * Implemented as in: arxiv:1406.3611v2, p.150, eq. (13.2.4).
     *
     * @param M_R resonance mass
     * @param W_R resonance width
     * @param J_R resonance spin
     * @param r_R resonance radius
     * @param m2_ab Dalitz plot variable
     * @param m_a 1st daughter mass
     * @param m_b 2nd daughter mass
     */
    template <typename T0, typename T1, typename T2>
    typename boost::math::tools::promote_args<T0,T1,T2>::type
    relativistic_width(double M_R, double W_R, double J_R, double r_R,
		       const T0& m2_ab, const T1& m_a, const T2& m_b) {

      typedef typename boost::math::tools::promote_args<T0,T1,T2>::type T_res;
      T_res res;
      res = W_R * M_R / sqrt(m2_ab) *
        pow(mfct::breakup_momentum::r2(m2_ab, M_R*M_R, m_a, m_b), J_R + 0.5) *
        pow(mfct::blatt_weisskopf(J_R, r_R*r_R, m2_ab, m_a, m_b), 2) /
        pow(mfct::blatt_weisskopf(J_R, r_R*r_R, M_R*M_R, m_a, m_b), 2);

      return res;
    }


  }
}
}
#endif
