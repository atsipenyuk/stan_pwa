#ifndef STAN_PWA__SRC__PHASE_SPACE_GEN_HPP
#define STAN_PWA__SRC__PHASE_SPACE_GEN_HPP

#include <stan_pwa/src/structures/three_body_resonances.hpp>

namespace mresonances = stan_pwa::resonances;

namespace stan {
  namespace math {

    /**
     * scalar phase_space_gen(vector)
     *
     * Takes an event vector; if this event lies in the phase space,
     * then return 1; otherwise, return 0. Basically, it is the 
     * probability function that allows us to sample the phase space region.
     *
     */
    template <typename T>
    T
    phase_space_gen(const Eigen::Matrix<T, Eigen::Dynamic,1>& y) {
      T res;
      res = stan_pwa::resonances::flat_D3pi.value(y(0,0), y(1,0))[0];
      return res;
    };

  }
}
#endif
