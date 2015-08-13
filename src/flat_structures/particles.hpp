#ifndef STAN_PWA__SRC__FLAT_STRUCTURES__PARTICLES_HPP
#define STAN_PWA__SRC__FLAT_STRUCTURES__PARTICLES_HPP

#include <stan_pwa/src/flat_structures/particles_def.hpp>

namespace stan_pwa {
namespace particles {

  // "Particles"
  Particle pi(0.13957018, 5., 0); // Charged pion
  Particle k(0.493677, 5., 0);  // Charged kaon
  Particle d(1.86484, 5., 0); // Charged D meson
  Particle D0(1.86961, 5., 0); // neutral D meson

  // "Resonances"
  /* The description of resonances as particles needs certain 
   * clarification. Namely, if we describe resonances as particles,
   * why do we need the radius, and what about the width of the
   * resonance? Well, we need the radius, because in the four
   * body decay we may use resonance as a parent particle.
   * We do not define the widths here, because width is more
   * model-dependent (in a certain manner); therefore, it seems
   * more convenient to define width of a resonance while
   * instantiating a given decay channel (see 'resonances.hpp').
   */
  // Spin 0
  Particle toy0_1000(1.000, 5., 0); // Toy resonance
  Particle toy0_1200(1.200, 5., 0); // Toy resonance

  Particle f0_500(0.475, 5., 0); // a.k.a sigma
  Particle f0_980(0.990, 5., 0);
  Particle f0_1370(1.350, 5., 0);
  Particle f0_1500(1.505, 5., 0);

  // Spin 1
  Particle a1(1.230, 5., 1);
  Particle rho_770(0.77526, 5., 1);
  Particle omega_782(0.78265, 5., 1);

  // Spin 2
  Particle f2_1270(1.2751, 5., 2);

}
}
#endif
