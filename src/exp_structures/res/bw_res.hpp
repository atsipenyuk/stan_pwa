#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODY__BW_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODY__BW_HPP


#include <cmath> // sqrt

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/src/exp_structures/res/base_3d_res.hpp>

namespace mc = stan_pwa::complex;
namespace mfct = stan_pwa::fct;


namespace stan_pwa {

  // Breit-Wigner resonance, 3-body decay
  class Bw_res : public Base_3d_res {
  public:
    Bw_res(particle P, particle a, particle b, 
	   particle c, particle R, double W) :
      Base_3d_res(P, a, b, c), R_(R), W_(W) {};


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABC (not symmetrized)
    template <typename T>
    C_t<promoted_t<T> >
    //    std::vector<T>
    value(const T& m2_ab, const T& m2_bc) {

      // Check if we are in the correct region of the Dalitz plot
      if (mfct::valid(m2_ab, m2_bc, 
		      this->P_, this->a_, this->b_, this->c_) == true) {

	T m_ab = sqrt(m2_ab);

        // Form factor P -> Rc
        T F_P = mfct::blatt_weisskopf(this->R_.J, this->P_.r2, 
				      this->P_.m2, m_ab, this->c_.m) /
	  mfct::blatt_weisskopf(this->R_.J, this->P_.r2, 
				this->P_.m2, this->R_.m, this->c_.m);

        // Form factor R -> ab
        T F_R = mfct::blatt_weisskopf(this->R_.J, this->R_.r2, 
				      m2_ab, this->a_.m, this->b_.m)/
                mfct::blatt_weisskopf(this->R_.J, this->R_.r2, 
				      this->R_.m2, this->a_.m, this->b_.m);

        T width = mfct::breit_wigner::relativistic_width(this->R_.m, this->W_,
							 this->R_.J, 
							 this->R_.r,
							 m2_ab, this->a_.m, 
							 this->b_.m);

	std::vector<T> T_R = mfct::breit_wigner::value(this->R_.m,
						       m2_ab, width);
	// If the parent particle does not have spin 0, some adjustments
	// must be performed in this Zemach function (use angular orbital
	// momentum between P and R instead of R.J)
        T Z = mfct::zemach(this->R_.J, m2_ab, m2_bc, 
			   this->P_.m, this->a_, this->b_, this->c_);

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
    C_t<promoted_t<T> >
    //    std::vector<T>
    value_sym(const T& m2_ab, const T& m2_bc) {
      return mc::scalar::add(this->value(m2_ab, m2_bc),
			     this->value(m2_bc, m2_ab));
    }

  private:
    // A BW resonance has the same properties as a particle, and a width
    const particle R_; // "Resonance = particle + width"
    const double W_; // Width of the resonance

  };

}

#endif
