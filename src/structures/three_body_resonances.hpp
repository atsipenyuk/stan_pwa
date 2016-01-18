#ifndef STAN_PWA__SRC__STRUCTURES__THREE_BODY_RESONANCES_HPP
#define STAN_PWA__SRC__STRUCTURES__THREE_BODY_RESONANCES_HPP

#include <stan_pwa/src/structures/struct_resonances.hpp>
namespace mresonances = stan_pwa::resonances;

namespace stan_pwa {
namespace resonances {

  // 3-body decay, Spin 0
  mresonances::flat_3 flat_D3pi(particles::d, particles::pi, particles::pi, 
				particles::pi);

  mresonances::breit_wigner toy0_1000(particles::d, particles::pi,
				      particles::pi, particles::pi,
				      particles::toy0_1000, 0.1);

  mresonances::breit_wigner toy0_1200(particles::d, particles::pi,
				      particles::pi, particles::pi,
				      particles::toy0_1200, 0.1);

  mresonances::flatte toy0_flatte(particles::d, particles::pi,
				  particles::pi, particles::pi,
				  particles::toy0_1000, 0.329, 2*0.329);

  mresonances::flatte f0_980(particles::d, particles::pi,
			     particles::pi, particles::pi,
			     particles::f0_980, 0.329, 2*0.329);

  mresonances::breit_wigner f0_500(particles::d, particles::pi,
				   particles::pi, particles::pi,
				   particles::f0_500, 0.800);

  mresonances::breit_wigner f0_1370(particles::d, particles::pi,
				    particles::pi, particles::pi,
				    particles::f0_1370, 0.350);

  mresonances::breit_wigner f0_1500(particles::d, particles::pi,
				    particles::pi, particles::pi,
				    particles::f0_1370, 0.109);


  // 3-body decay, Spin 1
  mresonances::breit_wigner rho_770(particles::d, particles::pi,
				    particles::pi, particles::pi,
				    particles::rho_770, 0.1491);

  mresonances::breit_wigner_only rho_770_bw_only(particles::d, particles::pi,
						 particles::pi, particles::pi,
						 particles::rho_770, 0.1491);


  // 3-body decay, Spin 2
  mresonances::breit_wigner f2_1270(particles::d, particles::pi,
				    particles::pi, particles::pi,
				    particles::f2_1270, 0.1852);

}
} 
#endif
