#ifndef STAN_PWA__SRC__STRUCTURES__FOUR_BODY__P_Rcd_abcd_HPP
#define STAN_PWA__SRC__STRUCTURES__FOUR_BODY__P_Rcd_abcd_HPP

#include <cmath> // sqrt
#include <math.h> // isnan

#include <stan_pwa/src/fct.hpp> // Breit-Wigner, Blatt-Weisskopf, etc.
#include <stan_pwa/src/complex.hpp> // Complex numbers

#include "base.hpp" // base class

#include <assert.h>

namespace mfct = stan_pwa::fct;
namespace mcomplex = stan_pwa::complex;
namespace mresonances = stan_pwa::resonances;


namespace stan_pwa {
namespace resonances {

  // Amplitude function for the 4 particle decay
  //   P -> R c d -> a b c d
  //
  // Model-dependent description is glued from
  // two Breit-Wigner 3-body decays
  //   P   -> R c d
  //   R   -> a b
  //
  // For detailed description, see
  //   https://www.overleaf.com/2630893ckqcxt#/6945226/
  // or feel free to contact me at arseniy.tsipenyuk@gmail.com.
  struct P_Rcd_abcd : resonances::resonance_base_4
  {

    const int l_1; // Orbital angular momentum between R_1 and d
    const int l_2; // Orbital angular momentum between R_2 and c
    const int l_3; // Orbital angular momentum between a and b
    const particle R; // resonance
    const double W_R; // Width of the resonance
 

    // Default constructor
    P_R1d_R2cd_abcd(particle _P, particle _a, particle _b, 
		    particle _c, particle _d,
		    int _l_1, int _l_2, int _l_3,  
		    particle _R_1, particle _R_2,
		    double _W_R_1, double _W_R_2) : 
      resonance_base_4(_P,_a, _b, _c, _d), 
      l_1(_l_1), l_2(_l_2), l_3(_l_3),
      R_1(_R_1), R_2(_R_2), W_R_1(_W_R_1), W_R_2(_W_R_2) {};


  private:

    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABCD (not symmetrized)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    // m2_12 is the invariant square mass of particles a and b.
    // Analogously, m2_34 is i.sq.m. of c and d, m2_23 - of b and c, etc.
    value_dynamic(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13) {

      typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;


    }



    template <typename T0>
    std::vector<T0>
    value_angular(theta_z_values v)
    {
      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;



      return 0;
    }



    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13, theta_z_values v) {
      return mcomplex::scalar::mult(value_dynamic(m2_12,m2_14,m2_23,m2_34,m2_13),
          value_angular(v));
    }


  public:

    // Evaluates the resonance at the given point
    // for the decay P -> ABCD (symmetrized, i.e. A==C, B==D)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_sym(const T0& m2_12, const T1& m2_14, const T2& m2_23,
	      const T3& m2_34, const T4& m2_13, std::vector<theta_z_values> theta_z_values_sym) {

      return mcomplex::scalar::add(this->value(m2_12,m2_14,m2_23,m2_34,m2_13, theta_z_values_sym[0]),
             mcomplex::scalar::add(this->value(m2_23,m2_34,m2_12,m2_14,m2_13, theta_z_values_sym[1]), // 1 <-> 3
             mcomplex::scalar::add(this->value(m2_14,m2_12,m2_34,m2_23,m2_13, theta_z_values_sym[2]), // 2 <-> 4
                                   this->value(m2_34,m2_23,m2_14,m2_12,m2_13, theta_z_values_sym[3]))));  // 1 <-> 3, 2 <-> 4
    }


  };

}
}

#endif

