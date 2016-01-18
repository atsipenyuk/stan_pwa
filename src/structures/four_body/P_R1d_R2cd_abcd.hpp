#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__P_R1D_R2cd_abcd_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__P_R1D_R2cd_abcd_HPP

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
  //   P -> R_1 d -> R_2 c d -> a b c d
  //
  // Model-dependent description is glued from
  // two Breit-Wigner 3-body decays
  //   P   -> R_2 c d
  //   R_1 -> a   b c
  //
  // For detailed description, see
  //   https://www.overleaf.com/2630893ckqcxt#/6945226/
  // or feel free to contact me at arseniy.tsipenyuk@gmail.com.
  struct P_R1d_R2cd_abcd : mresonances::resonance_base_4
  {

    const int l_1; // Orbital angular momentum between R_1 and d
    const int l_2; // Orbital angular momentum between R_2 and c
    const int l_3; // Orbital angular momentum between a and b
    const particle R_1; // First decay resonance (e.g. a_1)
    const particle R_2; // 2nd order decay resonance (e.g. rho_0)
    const double W_R_1, W_R_2; // Width of the 1st, 2nd resonance
 

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

      const T_123 m2_123 = m2_12 + m2_13 + m2_23 - a.m2 - b.m2 - c.m2;

      // Form factor P -> R_1 d
      const T_123 F_P = mfct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2,
          sqrt(m2_123), this->d.m) /
          mfct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2,
              this->R_1.m, this->d.m);

      // Form factor R_1 -> R_2 c
      const T_123 F_R_1 = mfct::blatt_weisskopf(this->l_2, this->R_1.r2, m2_123,
          this->R_2.m, this->c.m) /
          // POSSIBLY m_12 instead of R_2.m above and below
          mfct::blatt_weisskopf(this->l_2, this->R_1.r2, this->R_1.m2,
              this->R_2.m, this->c.m);

      // Form factor R_2 -> a b
      const T0 F_R_2 = mfct::blatt_weisskopf(this->l_3, this->R_2.r2, m2_12,
          this->a.m, this->b.m) /
          mfct::blatt_weisskopf(this->l_3, this->R_2.r2, this->R_2.m2,
              this->a.m, this->b.m);

      // Dynamical (Breit-Wigner) form factor of the first resonance
      const T_123 width_R_1 = mfct::breit_wigner::relativistic_width(this->R_1.m, W_R_1,
                    this->l_2, this->R_1.r,
                    m2_123, this->R_2.m2,
                    c.m2);
      const std::vector<T_123> T_R_1 = mfct::breit_wigner::value(this->R_1.m,
                m2_123,width_R_1);

      // Dynamical (Breit-Wigner) form factor of the 2nd resonance
      const T_123 width_R_2 = mfct::breit_wigner::relativistic_width(this->R_2.m, W_R_2,
                    this->l_3, this->R_2.r,
                    m2_12, a.m2, b.m2);
      const std::vector<T_123> T_R_2 = mfct::breit_wigner::value(this->R_2.m,
                m2_12, width_R_2);


      // Combine the factors to the decay amplitude
      return mcomplex::scalar::mult(F_P * F_R_1 * F_R_2,
				mcomplex::scalar::mult(T_R_1,T_R_2));
    }



    template <typename T0>
    std::vector<T0>
    value_angular(mfct::theta_z_values<T0> v)
    {
      std::vector<T0> A(2, 0.0);

      if (!(v.cos2_theta_1 >= 0. && v.cos2_theta_1 <= 1.) ||
          !(v.cos2_theta_2 >= 0. && v.cos2_theta_2 <= 1.)) {
        std::cout << "cos2_theta out of range; Return 0 \n";
        return A;
      }

      //assert(v.cos2_theta_1 >= 0. && v.cos2_theta_1 <= 1.);
      assert(v.z2_1 >= 0.);
      assert(v.z2_2 >= 0.);

      T0 Z_1 = mfct::zemach(this->P.J,   this->R_1.J, l_1, v.z2_1, v.cos2_theta_1);
      T0 Z_2 = mfct::zemach(this->R_1.J, this->R_2.J, l_2, v.z2_2, v.cos2_theta_2);

      return mcomplex::scalar::complex(Z_1 * Z_2, 0.);
    }



    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13,
        mfct::theta_z_values<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type> v) {
      return mcomplex::scalar::mult(value_dynamic(m2_12,m2_14,m2_23,m2_34,m2_13),
          this->value_angular(v));
    }


  public:

    // Evaluates the resonance at the given point
    // for the decay P -> ABCD (symmetrized, i.e. A==C, B==D)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_sym(const T0& m2_12, const T1& m2_14, const T2& m2_23,
	      const T3& m2_34, const T4& m2_13,
	      std::vector<mfct::theta_z_values<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type> > theta_z_values_sym) {

      return mcomplex::scalar::add(this->value(m2_12,m2_14,m2_23,m2_34,m2_13, theta_z_values_sym[0]),
             mcomplex::scalar::add(this->value(m2_23,m2_34,m2_12,m2_14,m2_13, theta_z_values_sym[1]), // 1 <-> 3
             mcomplex::scalar::add(this->value(m2_14,m2_12,m2_34,m2_23,m2_13, theta_z_values_sym[2]), // 2 <-> 4
                                   this->value(m2_34,m2_23,m2_14,m2_12,m2_13, theta_z_values_sym[3]))));  // 1 <-> 3, 2 <-> 4
    }


  };

}
}

#endif

