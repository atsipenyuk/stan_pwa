#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODY__BW_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODY__BW_HPP


#include <cmath> // sqrt

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/src/structures/three_body/base.hpp>
namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;
namespace mresonances = stan_pwa::resonances;


namespace stan_pwa {
namespace resonances {

  // Breit-Wigner resonance, 3-body decay
  struct breit_wigner : public mresonances::resonance_base_3
  {
    // A BW resonance has the same properties as a particle, and a width
    const particle R; // "Resonance = particle + width"
    const double W; // Width of the resonance

    breit_wigner(particle _P, particle _a, particle _b, particle _c, 
		 particle _R, double _W) :
      resonance_base_3(_P, _a, _b, _c), R(_R), W(_W) {};


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABC (not symmetrized)
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

        T width = mfct::breit_wigner::relativistic_width(this->R.m, this->W,
							 this->R.J, this->R.r,
							 m2_ab, this->a.m, 
							 this->b.m);

	std::vector<T> T_R = mfct::breit_wigner::value(this->R.m,m2_ab,width);
	// If the parent particle does not have spin 0, some adjustments
	// must be performed in this Zemach function (use angular orbital
	// momentum between P and R instead of R.J)
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
