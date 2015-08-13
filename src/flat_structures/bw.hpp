#ifndef STAN_PWA__SRC__FLAT_STRUCTURES__THREE_BODY__BW_HPP
#define STAN_PWA__SRC__FLAT_STRUCTURES__THREE_BODY__BW_HPP


#include <cmath> // sqrt

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/src/flat_structures/particles_def.hpp>
#include <stan_pwa/src/typedefs.h>

namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;

namespace stan_pwa {

  // Breit-Wigner resonance, 3-body decay
  template <typename T, ExternalParticles4 const &Ext, 
	    Particle_w const &R>
  C_t<T>
  breit_wigner(const T& m2_ab, const T& m2_bc) {

    if (mfct::valid(m2_ab, m2_bc, Ext.P, Ext.a, Ext.b, Ext.c) == true) {
      
      T m_ab = sqrt(m2_ab);
      // Form factor P -> Rc
      T F_P = mfct::blatt_weisskopf(R.J, Ext.P.r2, Ext.P.m2, m_ab, Ext.c.m) /
	mfct::blatt_weisskopf(R.J, Ext.P.r2, Ext.P.m2, R.m, Ext.c.m);

      // Form factor R -> ab
      T F_R = mfct::blatt_weisskopf(R.J, R.r2, m2_ab, Ext.a.m, Ext.b.m)/
	mfct::blatt_weisskopf(R.J, R.r2, R.m2, Ext.a.m, Ext.b.m);

      T width = mfct::breit_wigner::relativistic_width(R.m, R.W, R.J, R.r,
						       m2_ab,Ext.a.m,Ext.b.m);
      
      std::vector<T> T_R = mfct::breit_wigner::value(R.m,m2_ab,width);
      // If the parent particle does not have spin 0, some adjustments
      // must be performed in this Zemach function (use angular orbital
      // momentum between P and R instead of R.J)
      T Z = mfct::zemach(R.J, m2_ab, m2_bc, Ext.P.m, Ext.a, Ext.b, Ext.c);
      
      C_t<T> res(2);
      res = mc::scalar::mult(F_P * F_R * Z, T_R);
      return res;
    }
    else {
      std::vector<T> res(2, 0.0);
      return res;
    }
  }


  // Evaluates the resonance at the given point in the Dalitz plot
  // for the decay P -> ABC (symmetrized, i.e. A==C)
  template <typename T, ExternalParticles4 const &Ext, 
	    Particle const &R, double const &W>
  C_t<T>
  breit_wigner_sym(const T& m2_ab, const T& m2_bc) {
    return mc::scalar::add(breit_wigner<T,Ext,R,W>(m2_ab, m2_bc),
			   breit_wigner<T,Ext,R,W>(m2_bc, m2_ab));
  };
  


} // end of stan_pwa
#endif
