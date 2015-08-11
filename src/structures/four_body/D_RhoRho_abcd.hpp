#ifndef STAN_PWA__SRC__STRUCTURES__FOUR_BODY__D_RhoRho_abcd_HPP
#define STAN_PWA__SRC__STRUCTURES__FOUR_BODY__D_RhoRho_abcd_HPP

#include <cmath> // sqrt
#include <math.h> // isnan

#include <stan_pwa/src/fct.hpp> // Breit-Wigner, Blatt-Weisskopf, etc.
#include <stan_pwa/src/complex.hpp> // Complex numbers

#include <stan_pwa/src/structures/four_body/base.hpp> // base class

#include <assert.h>

namespace resonances {

  // Amplitude function for the 4 particle decay
  //   P -> rho rho -> a b c d
  //
  // Model-dependent description is glued from
  // two Breit-Wigner 3-body decays
  //   P   -> R_1 R_2
  //   R_1 -> a b
  //   R_2 -> c d
  //
  // For detailed description, see
  //   https://www.overleaf.com/2630893ckqcxt#/6945226/
  // or feel free to contact me at arseniy.tsipenyuk@gmail.com.
  struct P_RhoRho_abcd : resonances::resonance_base_4
  {

    const double W_rho, W_omega; // Widths
    double phi; // phase between rho and omega

    // Default constructor
    P_RhoRho_abcd(particle _P, particle _a, particle _b,
		    particle _c, particle _d) :
      resonance_base_4(_P,_a, _b, _c, _d)
      // TODO: constructor
      {
        phi = (particles::omega.m * W_omega - particles::rho_770.m * W_rho) /
            (particles::omega.m2 - particles::rho_770.m2);

      };


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABCD (not symmetrized)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    // m2_12 is the invariant square mass of particles a and b.
    // Analogously, m2_34 is i.sq.m. of c and d, m2_23 - of b and c, etc.
    value(const T0& m2_12, const T1& m2_14, const T2& m2_23,
        const T3& m2_34, const T4& m2_13) {

      typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

      std::vector<T_res> A(2, 0.0);


      // rho rho lineshape from Crystal Barrel
      
      // modified Breit Wigner for rho from 12
      T_res p = sqrt( fct::breakup_momentum::p2(P.m2, sqrt(m2_12), sqrt(m2_34)) ); // rho rho breakup momentum
      T_res q_12 = sqrt( fct::breakup_momentum::p2(m2_12, a.m, b.m) ); // pi+_1 pi-_2 breakup momentum
      T_res q_rho = sqrt( fct::breakup_momentum::p2(particles::rho_770.m2, a.m, b.m) ); // rho -> pi+ pi- breakup momentum
      T_res q_omega = sqrt( fct::breakup_momentum::p2(particles::omega.m2, a.m, b.m) ); // omega -> pi+ pi- breakup momentum

      T_res_relativistic_width_rho_12 = fct::relativistic_width(particles::rho_770.m, W_rho, 1, particles::rho_770.r,
          m2_12, a.m, b.m);
      T_res Dp_rho = p    * particles::rho_770.r * sqrt( 2. / (pow(p    * particles::rho_770.r, 2) + 1.) );
      T_res Dq_rho = q_12 * particles::rho_770.r * sqrt( 2. / (pow(q_12 * particles::rho_770.r, 2) + 1.) );
      T_res rho_rho = 1. / particles::rho_770.m * sqrt( particles::rho_770.m2 - 4.*a.m*b.m);
      std::vector<T_res> BW_rho_12 = Dp_rho * particles::rho_770.m * W_rho / rho_rho * Dp_rho / Dq_rho *
          fct::breit_wigner::value(particles::rho_770.m, m2_12, T_res_relativistic_width_rho_12);

      T_res_relativistic_width_omega_12 = fct::relativistic_width(particles::omega.m, W_omega, 1, particles::omega.r,
          m2_12, a.m, b.m);
      T_res Dp_omega = p    * particles::omega.r * sqrt( 2. / (pow(p    * particles::omega.r, 2) + 1.) );
      T_res Dq_omega = q_12 * particles::omega.r * sqrt( 2. / (pow(q_12 * particles::omega.r, 2) + 1.) );
      T_res rho_omega = 1. / particles::omega.m * sqrt( particles::omega.m2 - 4.*a.m*b.m);
      std::vector<T_res> BW_omega_12 = Dp * particles::omega.m * W_omega / rho_omega * Dp / Dq_omega *
          fct::breit_wigner::value(particles::omega.m, m2_12, T_res_relativistic_width_omega_12);



      complex::scalar::complex corrFac = complex::scalar::mult(complex::scalar::complex(cos(phi), sin(phi)),
          complex::scalar::complex(delta * (particles::rho_770.m + particles::omega.m)), 0.);
      complex::scalar::mult( corrFac, complex::scalar::inverse(
          complex::scalar::complex( particles::omega.m2 - particles::rho_770.m2,
          -1.*particles::omega.m*W_omega + particles::rho_770.m * W_rho) ) );








      T_res q_34 = sqrt( fct::breakup_momentum::p2(m2_34, c.m, d.m) ); // pi+_3 pi-_4 breakup momentum
      T_res_relativistic_width_rho_34 = fct::relativistic_width(particles::rho_770.m, W_rho, 1, particles::rho_770.r,
          m2_34, c.m, d.m);



    }

  };

}

#endif
