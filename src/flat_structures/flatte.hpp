#ifndef STAN_PWA__SRC__FLAT_FLAT_STRUCTURES__THREE_BODY__FLATTE_HPP
#define STAN_PWA__SRC__FLAT_FLAT_STRUCTURES__THREE_BODY__FLATTE_HPP

#include <cmath> // sqrt

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/src/flat_flat_structures/three_body/base.hpp>
namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;
namespace mresonances = stan_pwa::resonances;

namespace stan_pwa {
namespace resonances {

  // Flatte Resonance for the 3-body decay
  struct flatte : public mresonances::resonance_base_3
  {
    // Flatte has same properties as a particle + 2 widths
    const particle R;
    const double G_pp;
    const double G_kk;

    flatte(particle _P, particle _a, particle _b, particle _c,
	   particle _R, double _G_pp, double _G_kk) :
      resonance_base_3(_P, _a, _b, _c), R(_R), G_pp(_G_pp), G_kk(_G_kk) {};

    // Returns the amplitude of the decay P->abc via Flatte resonance.
    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc) 
    {
      if (mfct::valid(m2_ab, m2_bc, 
		      this->P, this->a, this->b, this->c) == true) {

	T m_ab = sqrt(m2_ab);

        // Form factor P -> Rc
        T F_P = mfct::blatt_weisskopf(this->R.J, this->P.r2, 
				      this->P.m2, m_ab, this->c.m) /
	  mfct::blatt_weisskopf(this->R.J, this->P.r2,
				this->P.m2, this->R.m, this->c.m);

        // Form factor R -> ab
        T F_R = mfct::blatt_weisskopf(this->R.J, this->R.r2, 
				      m2_ab, this->a.m, this->b.m)/
	  mfct::blatt_weisskopf(this->R.J, this->R.r2,
				this->R.m2, this->a.m, this->b.m);

	std::vector<T> T_R = mfct::flatte::value(this->R.m, m2_ab,
						 this->G_pp, this->G_kk);
        T Z = mfct::zemach(this->R.J, m2_ab, m2_bc, 
			   this->P.m, this->a, this->b, this->c);

	std::vector<T> res(2);
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
    template <typename T>
    inline
    std::vector<T>
    value_sym(const T& m2_ab, const T& m2_bc) {
      return mc::scalar::add(this->value(m2_ab, m2_bc),
			     this->value(m2_bc, m2_ab));
    }

  };
}
}
#endif
