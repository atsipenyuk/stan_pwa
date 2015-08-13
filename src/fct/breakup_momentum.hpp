#ifndef STAN_PWA__SRC__FCT__BREAKUP_MOMENTUM_HPP
#define STAN_PWA__SRC__FCT__BREAKUP_MOMENTUM_HPP

#include <vector>
#include <cmath>

#include <stan_pwa/src/complex.hpp>
namespace mc = stan_pwa::complex;

namespace stan_pwa {
namespace fct {

  namespace breakup_momentum {

    /**
     * Return squared breakup momentum (R -> ab).
     *
     * @param m2_R decaying Particle squared mass
     * @param m_a 1st daughter mass
     * @param m_b 2nd daughter mass
     * @return breakup momentum
     */
    template <typename T0, typename T1, typename T2>
    inline
    typename boost::math::tools::promote_args<T0,T1,T2>::type
    p2(const T0& m2_R, const T1& m_a, const T2& m_b) {

      if (m_a == m_b) {
        return m2_R / 4.0 - m_a * m_a;
      }

      return (m2_R - (m_a + m_b) * (m_a + m_b)) * 
             (m2_R - (m_a - m_b) * (m_a - m_b)) / m2_R / 4.0;
    }


    /**
     * Return complex breakup momentum.
     *
     * @param m2_R decaying Particle squared mass
     * @param m_a 1st daughter mass
     * @param m_b 2nd daughter mass
     * @return breakup momentum
     */
    template <typename T0, typename T1, typename T2>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2>::type>
    complex_p(const T0& m2_R, const T1& m_a, const T2& m_b) {

      typedef typename boost::math::tools::promote_args<T0,T1,T2>::type T_res;

      T_res p2 = stan_pwa::fct::breakup_momentum::p2(m2_R, m_a, m_b);
      std::vector<T_res> res(2);      

      if (p2 >= 0)
        res = mc::scalar::complex(sqrt(p2), 0.);
      else
        res = mc::scalar::complex(0., sqrt(-p2));

      return res;
    }


    /**
     * Ratio of squared breakup momenta (N -> ab) / (D -> ab).
     *
     * @param m2_N 1st decaying Particle squared mass
     * @param m2_D 2nd decaying Particle squared mass
     * @param m_a 1st daughter mass
     * @param m_b 2nd daughter mass
     * @return breakup momentum ratio
     */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    typename boost::math::tools::promote_args<T0,T1,T2,T3>::type
    r2(const T0& m2_N, const T1& m2_D, const T2& m_a, const T3& m_b) {

      typedef typename boost::math::tools::promote_args<T0,T2,T3>::type T_r1;
      typedef typename boost::math::tools::promote_args<T1,T2,T3>::type T_r2;

      T_r1 p2_N = stan_pwa::fct::breakup_momentum::p2(m2_N, m_a, m_b);
      T_r2 p2_D = stan_pwa::fct::breakup_momentum::p2(m2_D, m_a, m_b);

      return p2_N / p2_D;
    }



  }
}
}
#endif
