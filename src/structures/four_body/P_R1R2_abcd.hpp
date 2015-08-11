#ifndef STAN_PWA__SRC__STRUCTURES__FOUR_BODY__P_R1R2_abcd_HPP
#define STAN_PWA__SRC__STRUCTURES__FOUR_BODY__P_R1R2_abcd_HPP

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
  //   P -> R_1 R_2 -> a b c d
  //
  //   P   -> R_1 R_2
  //   R_1 -> a b
  //   R_2 -> c d
  //
  // For detailed description, see
  //   https://www.overleaf.com/2630893ckqcxt#/6945226/
  // or feel free to contact me at arseniy.tsipenyuk@gmail.com.
  struct P_R1R2_abcd : mresonances::resonance_base_4
  {
    const int l_1; // Orbital angular momentum between R_1 and R_2
    const int l_2; // Orbital angular momentum between a and b
    const int l_3; // Orbital angular momentum between c and d
    const particle R_1; //
    const particle R_2; //
    const double W_R_1, W_R_2; // Widths

    // Default constructor
    P_R1R2_abcd(particle _P, particle _a, particle _b,
		    particle _c, particle _d,
		    int _l_1, int _l_2, int _l_3,
		    particle _R_1, particle _R_2,
		    double _W_R_1, double _W_R_2) :
      resonance_base_4(_P,_a, _b, _c, _d),
      l_1(_l_1), l_2(_l_2), l_3(_l_3),
      R_1(_R_1), R_2(_R_2),
      W_R_1(_W_R_1), W_R_2(_W_R_2) {};


  private:

    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABCD (not symmetrized)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_dynamic(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13) {

      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

      const T_res F_P = mfct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2,
                sqrt(m2_12), sqrt(m2_34));

      const T_res relativistic_width_R_1 =
          mfct::breit_wigner::relativistic_width(R_1.m, W_R_1, l_2, R_1.r, m2_12, a.m, b.m);
      const std::vector<T_res> BW_R_1 = mfct::breit_wigner::value(R_1.m, m2_12, relativistic_width_R_1);

      const T_res F_R_1 = mfct::blatt_weisskopf(this->l_2, R_1.r2, R_1.m2, a.m, b.m);
      const T_res relativistic_width_R_2 =
          mfct::breit_wigner::relativistic_width(R_2.m, W_R_2, l_3, R_2.r, m2_34, c.m, d.m);
      const std::vector<T_res> BW_R_2 = mfct::breit_wigner::value(R_2.m, m2_34, relativistic_width_R_2);
      const T_res F_R_2 = mfct::blatt_weisskopf(this->l_3, R_2.r2, R_2.m2, c.m, d.m);


      // Combine the factors to the decay amplitude
      return mcomplex::scalar::mult(F_P * F_R_1 * F_R_2,
          mcomplex::scalar::mult(BW_R_1, BW_R_2) );
    }


    template <typename T0>
    std::vector<T0>
    value_angular_parallel(mfct::helicity_angles<T0> a)
    {
      return mcomplex::scalar::complex(1./sqrt(2.) * cos(a.chi) * sin(a.theta_1) * sin(a.theta_2), 0);
    }

    template <typename T0>
    std::vector<T0>
    value_angular_perpendicular(mfct::helicity_angles<T0> a)
    {
      return mcomplex::scalar::complex(0, 1./sqrt(2.) * sin(a.chi) * sin(a.theta_1) * sin(a.theta_2));
    }

    template <typename T0>
    std::vector<T0>
    value_angular_longitudinal(mfct::helicity_angles<T0> a)
    {
      return mcomplex::scalar::complex(cos(a.theta_1) * cos(a.theta_2), 0);
    }


    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_parallel(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13, mfct::helicity_angles<T0> v) {
      return mcomplex::scalar::mult(value_dynamic(m2_12,m2_14,m2_23,m2_34,m2_13),
          this->value_angular_parallel(v));
    }

    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_perpendicular(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13, mfct::helicity_angles<T0> v) {
      return mcomplex::scalar::mult(value_dynamic(m2_12,m2_14,m2_23,m2_34,m2_13),
          this->value_angular_perpendicular(v));
    }

    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_longitudinal(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13, mfct::helicity_angles<T0> v) {
      return mcomplex::scalar::mult(value_dynamic(m2_12,m2_14,m2_23,m2_34,m2_13),
          this->value_angular_longitudinal(v));
    }


  public:

    // Evaluates the resonance at the given point
    // for the decay P -> ABCD (symmetrized, i.e. A==C, B==D)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_sym_parallel(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13,
        std::vector<mfct::helicity_angles<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type> > helicity_angles_sym) {

      return mcomplex::scalar::add(this->value_parallel(m2_12,m2_14,m2_23,m2_34,m2_13, helicity_angles_sym[0]),
             mcomplex::scalar::add(this->value_parallel(m2_23,m2_34,m2_12,m2_14,m2_13, helicity_angles_sym[1]), // 1 <-> 3
             mcomplex::scalar::add(this->value_parallel(m2_14,m2_12,m2_34,m2_23,m2_13, helicity_angles_sym[2]), // 2 <-> 4
                                   this->value_parallel(m2_34,m2_23,m2_14,m2_12,m2_13, helicity_angles_sym[3]))));  // 1 <-> 3, 2 <-> 4
    }


    // Evaluates the resonance at the given point
    // for the decay P -> ABCD (symmetrized, i.e. A==C, B==D)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_sym_perpendicular(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13,
        std::vector<mfct::helicity_angles<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type> > helicity_angles_sym) {

      return mcomplex::scalar::add(this->value_perpendicular(m2_12,m2_14,m2_23,m2_34,m2_13, helicity_angles_sym[0]),
             mcomplex::scalar::add(this->value_perpendicular(m2_23,m2_34,m2_12,m2_14,m2_13, helicity_angles_sym[1]), // 1 <-> 3
             mcomplex::scalar::add(this->value_perpendicular(m2_14,m2_12,m2_34,m2_23,m2_13, helicity_angles_sym[2]), // 2 <-> 4
                                   this->value_perpendicular(m2_34,m2_23,m2_14,m2_12,m2_13, helicity_angles_sym[3]))));  // 1 <-> 3, 2 <-> 4
    }


    // Evaluates the resonance at the given point
    // for the decay P -> ABCD (symmetrized, i.e. A==C, B==D)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    value_sym_longitudinal(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13,
        std::vector<mfct::helicity_angles<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type> > helicity_angles_sym) {

      return mcomplex::scalar::add(this->value_longitudinal(m2_12,m2_14,m2_23,m2_34,m2_13, helicity_angles_sym[0]),
             mcomplex::scalar::add(this->value_longitudinal(m2_23,m2_34,m2_12,m2_14,m2_13, helicity_angles_sym[1]), // 1 <-> 3
             mcomplex::scalar::add(this->value_longitudinal(m2_14,m2_12,m2_34,m2_23,m2_13, helicity_angles_sym[2]), // 2 <-> 4
                                   this->value_longitudinal(m2_34,m2_23,m2_14,m2_12,m2_13, helicity_angles_sym[3]))));  // 1 <-> 3, 2 <-> 4
    }

  };

}
}

#endif
