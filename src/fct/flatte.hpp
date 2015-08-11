#ifndef STAN_PWA__SRC__FCT__FLATTE_HPP
#define STAN_PWA__SRC__FCT__FLATTE_HPP

#include <vector>
#include <cmath>

#include <stan_pwa/src/structures/particles.hpp> // particles::pi, k
#include <stan_pwa/src/fct/breakup_momentum.hpp> // breakup_momentum::complex_p
#include <stan_pwa/src/complex.hpp>
namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;


namespace stan_pwa {
namespace fct {
  namespace flatte {

    /**
     * Return complex Flatte form factor.
     *
     * TODO: Eliminate the dependence on specific particles (pi, k)
     * Implemented as in: arxiv:1406.6311v2, p. 150, eq. (13.2.5)
     *
     * @param M_R resonance mass
     * @param m2 Dalitz plot variable (squared mass)
     * @param gpp phase-space factor*coupling const**2 of the channel -> pi+pi
     * @param gkk -//-                                 of the channel -> K+K
     * @return Flatte dynamical form factor
     */
    template <typename T0, typename T1, typename T2, typename T3>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3>::type>
    value(const T0& M_R, const T1& m2_ab, const T2& gpp, const T3& gkk) {

      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3>::type T_res;
      std::vector<T_res> pp(2);
      std::vector<T_res> kk(2);
      std::vector<T_res> temp(2);

      pp = mc::scalar::mult(mc::scalar::complex(gpp * gpp, 0.0),
             mfct::breakup_momentum::complex_p(m2_ab, particles::pi.m, 
					      particles::pi.m));

      kk = mc::scalar::mult(mc::scalar::complex(gkk * gkk, 0.0),
             mfct::breakup_momentum::complex_p(m2_ab, particles::k.m,
					      particles::k.m));
     
      temp = mc::scalar::mult(mc::scalar::complex(0.0, 1.0), 
				   mc::scalar::add(pp,kk));

      return mc::scalar::inverse(
        mc::scalar::subtract(
	mc::scalar::complex(M_R*M_R - m2_ab, 0.0), 
	mc::scalar::mult(mc::scalar::complex(2. / sqrt(m2_ab),0.0),
			      temp)));
    }
  }
}
}
#endif
