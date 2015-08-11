#ifndef PWA_STAN__SRC__MODEL_INST_HPP
#define PWA_STAN__SRC__MODEL_INST_HPP

#include <boost/any.hpp>

#include <stan_pwa/src/structures.hpp>
#include "model_def.hpp"

namespace stan_pwa {

  /* EDIT THE FOLLOWING SECTION -- YOU NEED TO EDIT TWO THINGS**********/
  /**
   * DO THIS (1): Declare your model-dependent resonances. (The vector
   * resonance_list is instantiated with them.)
   * (Need suggestions? 
   *  Look at stan_pwa/src/structures/three_body_resonances.hpp)
   * (Want to know which particles are defined? 
   *  Look at stan_pwa/src/structures/particles.hpp)
   */
  // rho_770
  resonances::breit_wigner rho_770 = 
    resonances::breit_wigner(particles::d, particles::pi,
			     particles::pi, particles::pi,
			     particles::rho_770, 0.1491);
  // f0_1370
  resonances::breit_wigner f0_1370 = 
    resonances::breit_wigner(particles::d, particles::pi,
			     particles::pi, particles::pi,
			     particles::f0_1370, 0.350);
    
  std::vector<resonances::resonance_base_3> resonance_list = 
  {rho_770, f0_1370};
  
  ///> DO THIS (2): Should your model be symmetrized? If yes, set sym_flag
  ///> to 1. Else, set to 0.
  bool sym_flag = 1;
  
  /* END of edited section **********************************************/
  
  
  // Declare the model
  // First argument tells how many variables we have 
  // (two for 3-body decay, five for 4-body-decay).
  Model MyModel = Model(2,sym_flag,resonance_list);

  // Define pointers to vector amplitudes depending on the symmetry of the
  // model
  template <typename T>
  using APTR_t = C_t<T> (Model::*) (const unsigned int, const Var_t<T>&);

  template <typename T>
  using AVPTR_t = CV_t<T> (Model::*) (const Var_t<T>&);

  // By default, point to the symmetric model
  template <typename T>
  APTR_t<T> amp_ptr=&Model::amplitude_sym;
  template <typename T>
  AVPTR_t<T> amp_vector_ptr=&Model::amplitude_vector_sym;

  template <typename T>
  class Wrapper {
  public:
    Wrapper(bool a_sym_flag,
	    APTR_t<T> a_amp_ptr,
	    AVPTR_t<T> a_amp_vector_ptr)
    {
      if (a_sym_flag == 0) {
	  a_amp_ptr = &Model::amplitude;
	  a_amp_vector_ptr = &Model::amplitude_vector;
      } // end if
    };    
  };


}
#endif 









