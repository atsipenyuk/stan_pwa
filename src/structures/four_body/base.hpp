#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__BASE_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__BASE_HPP


#include <meson_deca/lib/c_lib/structures/struct_particles.hpp> // Particles

namespace stan_pwa {
namespace resonances {

  struct resonance_base_4 // Base struct for the 4 particle decay
  {

    // Particle masses and radii are stored in the structure 'particle'
    const particle P; // Parent particle
    const particle a; // Final state particles a,b,c,d
    const particle b;
    const particle c;
    const particle d;

    resonance_base_4(particle _P, particle _a, particle _b,
		     particle _c, particle _d) :
      P(_P), a(_a), b(_b), c(_c), d(_d) {};
  };

}
}

#endif
