#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODY__FLAT_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODY__FLAT_HPP


#include <cmath> // sqrt

#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/structures/three_body/base.hpp>
namespace mfct = stan_pwa::fct;
namespace mresonances = stan_pwa::resonances;


namespace stan_pwa {
namespace resonances {

  // Non-resonant (flat) resonance, 3-body decay
  struct flat_3 : public mresonances::resonance_base_3
  {
    // Constructor
    flat_3(particle _P, particle _a, particle _b, particle _c) : 
      resonance_base_3(_P, _a, _b, _c) {};

  
    // Returns 1 if we are within Dalitz plot bounds, 0 else.
    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc) {
      std::vector<T> res(2, 0.0);
      if (mfct::valid(m2_ab, m2_bc, 
		     this->P, this->a, this->b, this->c) == true) {
	res[0] = 1.0;
      }
      return res;
    }


    // Returns 1 if we are within Dalitz plot bounds, 0 else.
    // For flat background symmetrized and non-symmetrized functions are
    // the same.
    template <typename T>
    std::vector<T>
    value_sym(const T& m2_ab, const T& m2_bc) {
      std::vector<T> res(2, 0.0);
      if (mfct::valid(m2_ab, m2_bc, 
		     this->P, this->a, this->b, this->c) == true) {
	res[0] = 1.0;
      }
      return res;
    }


  };

}
}
#endif
