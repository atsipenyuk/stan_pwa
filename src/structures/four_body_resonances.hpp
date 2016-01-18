#ifndef STAN_PWA__SRC__STRUCTURES__FOUR_BODY_RESONANCES_HPP
#define STAN_PWA__SRC__STRUCTURES__FOUR_BODY_RESONANCES_HPP

#include <stan_pwa/src/structures/struct_resonances.hpp>
namespace mresonances = stan_pwa::resonances;

namespace stan_pwa {
namespace resonances {

  // 4-body resonances

  // Consecutive decays
  // To be adjusted: width of a1
  mresonances::P_R1d_R2cd_abcd D_a_rho_S_wave(particles::D0,
					     particles::pi, particles::pi, 
					     particles::pi, particles::pi,
					     1, 0, 1,
					     particles::a1, 
					     particles::rho_770, 
					     0.1, 0.1491);

  mresonances::P_R1d_R2cd_abcd D_a_rho_D_wave(particles::D0,
					     particles::pi, particles::pi, 
					     particles::pi, particles::pi,
					     1, 2, 1,
					     particles::a1, 
					     particles::rho_770, 
					     0.1, 0.1491); 

  // To be adjusted: width of a1 and sigma
  mresonances::P_R1d_R2cd_abcd D_a_sigma(particles::D0,
               particles::pi, particles::pi,
               particles::pi, particles::pi,
               1, 1, 0, // ????
               particles::a1,
               particles::f0_500, // = sigma
               0.1, 0.550);

  // R1 R2 decays
  mresonances::P_R1R2_abcd D_rho_rho(particles::D0,
               particles::pi, particles::pi,
               particles::pi, particles::pi,
               1, 1, 1,
               particles::rho_770, particles::rho_770,
               0.1491, 0.1491);

  mresonances::P_R1R2_abcd D_omega_omega(particles::D0,
               particles::pi, particles::pi,
               particles::pi, particles::pi,
               1, 1, 1,
               particles::omega_782, particles::omega_782,
               0.00849, 0.00849);

  mresonances::P_R1R2_abcd D_rho_omega(particles::D0,
               particles::pi, particles::pi,
               particles::pi, particles::pi,
               1, 1, 1,
               particles::rho_770, particles::omega_782,
               0.1491, 0.00849);

}
} 
#endif
