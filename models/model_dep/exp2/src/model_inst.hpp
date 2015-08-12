#ifndef SRC__STRUCTURES__MODEL_INST_HPP
#define SRC__STRUCTURES__MODEL_INST_HPP


#include <stan_pwa/src/exp_structures/particles.hpp>
#include <stan_pwa/src/exp_structures/res.hpp>
#include <stan_pwa/src/model_def.hpp>

namespace stan_pwa {

  particle P = particles::d;
  particle a = particles::pi;
  particle b = particles::pi;
  particle c = particles::pi;

  Sqr my_res_1(P, a, b, c);//, particles::rho_770, 0.1491);
  Sqr my_res_2(P, a, b, c);//, particles::f0_1370, 0.350);
  
  std::vector<Base_3d_res> my_res_list = {my_res_1, my_res_2};
  bool sym_flag = 1;
  
  Model My_model(my_res_list, 2);

}
#endif
