#ifndef PWA_STAN__SRC__RESONANCE_WRAPPER_HPP
#define PWA_STAN__SRC__RESONANCE_WRAPPER_HPP


#include <stan_pwa/src/structures.hpp>

// Workaround to the fact that we can not create pointers
// to member functions (such as resonencos).

namespace stan_pwa {

  // Pointers to model-dependent Breit-Wigner resonances
  typedef std::vector<boost::any> (resonances::resonance_base_3::*val_ptr_t)
  (const boost::any& y0, const boost::any& y1);

  val_ptr_t val_ptr = &resonances::resonance_base_3::value;
  val_ptr_t val_sym_ptr = &resonances::resonance_base_3::value_sym;


  /*  std::vector<boost::any> 
  val(const resonances::resonance_base_3& MyRes, 
      const boost::any& y0, const boost::any& y1)pp
  {
    return MyRes.value(y0,y1);
  }

  std::vector<boost::any> 
  val_sym(const resonances::resonance_base_3& MyRes, 
	  const boost::any& y0, const boost::any& y1)
  {
    return MyRes.value_sym(y0,y1);
    }*/
}

#endif
