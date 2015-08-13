#ifndef SRC__MODEL_DEF_HPP
#define SRC__MODEL_DEF_HPP

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/typedefs.h>
#include <stan_pwa/src/exp_structures/res.hpp>

namespace mc = stan_pwa::complex;

namespace stan_pwa {
  const int num_res_ =2;
  const int num_var_ = 2;

  template <typename T>
  struct Model {
    Model() :
      amplitudes(0) {};

    void add_amplitude(Amp_ptr_2<T> a)
    {
      this->amplitudes.push_back(a);
    }

    std::vector<Amp_ptr_2<T> > amplitudes; 
  };

  //Model My_Model();
}
#endif

