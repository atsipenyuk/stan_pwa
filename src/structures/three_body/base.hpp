#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODY__BASE_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODY__BASE_HPP

#include <stan_pwa/src/structures/struct_particles.hpp>

namespace stan_pwa {
namespace resonances {

  // Base struct for model-dependent 3-body-decay resonances
  struct resonance_base_3
  {

    const stan_pwa::particle P; // Parent particle
    const stan_pwa::particle a; // Final state particles
    const stan_pwa::particle b; 
    const stan_pwa::particle c;

    resonance_base_3(stan_pwa::particle _P, stan_pwa::particle _a, 
		     stan_pwa::particle _b, stan_pwa::particle _c) :
      P(_P), a(_a), b(_b), c(_c) {};
  };

}
}
#endif
