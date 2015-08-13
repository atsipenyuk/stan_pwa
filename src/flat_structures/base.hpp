#ifndef STAN_PWA__SRC__FLAT_FLAT_STRUCTURES__THREE_BODY__BASE_HPP
#define STAN_PWA__SRC__FLAT_FLAT_STRUCTURES__THREE_BODY__BASE_HPP

#include <boost/any.hpp>

#include <stan_pwa/src/flat_flat_structures/particles_def.hpp>

namespace stan_pwa {
namespace resonances {

  ///> Base struct for model-dependent 3-body-decay resonances
  ///> ALL 3-body resonances MUST be derived from this struct.
  class resonance_base_3
  {
  public:
    resonance_base_3(stan_pwa::particle _P, stan_pwa::particle _a, 
		     stan_pwa::particle _b, stan_pwa::particle _c) :
      P(_P), a(_a), b(_b), c(_c) {};

    template <typename T>    
    std::vector<T> value(const T&, const T&);
    template <typename T>
    std::vector<T> value_sym(const T&, const T&);

    const stan_pwa::particle P; ///> Parent particle
    const stan_pwa::particle a; ///> Final state particles
    const stan_pwa::particle b; 
    const stan_pwa::particle c;
  };

}
}
#endif
