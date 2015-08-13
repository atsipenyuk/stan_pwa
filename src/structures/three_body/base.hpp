#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODYaoeu__BASE_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODYaoeu__BASE_HPP

#include <boost/any.hpp>

#include <stan_pwa/src/flat_structures/particles_def.hpp>

namespace stan_pwa {
namespace resonances {

  ///> Base struct for model-dependent 3-body-decay resonances
  ///> ALL 3-body resonances MUST be derived from this struct.
  class resonance_base_3
  {
  public:
    resonance_base_3(stan_pwa::Particle _P, stan_pwa::Particle _a, 
		     stan_pwa::Particle _b, stan_pwa::Particle _c) :
      P(_P), a(_a), b(_b), c(_c) {};

    template <typename T>    
    std::vector<T> value(const T&, const T&);
    template <typename T>
    std::vector<T> value_sym(const T&, const T&);

    const stan_pwa::Particle P; ///> Parent Particle
    const stan_pwa::Particle a; ///> Final state Particles
    const stan_pwa::Particle b; 
    const stan_pwa::Particle c;
  };

}
}
#endif
