#ifndef STAN_PWA__SRC__FCT__P_V1V2_ANGLES_HPP
#define STAN_PWA__SRC__FCT__P_V1V2_ANGLES_HPP

#include <cmath> // sqrt
#include <math.h> // isnan

#include "breakup_momentum.hpp"


namespace mfct = stan_pwa::fct;

namespace stan_pwa {
namespace fct {

template <typename T0>
struct helicity_angles {
  T0 theta_1;
  T0 theta_2;
  T0 chi;
};

/**
 * Calculate helicity angles theta_1 and theta_2, and acoplanarity angle chi
 * from m2_12, m2_14, m2_23, m2_34, m2_13.
 *
 * See LINK paper p. 30
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4>
helicity_angles<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type>
P_V1V2_angles(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
    const T3 &m2_34, const T4& m2_13,
    const Particle& P, const Particle& a, const Particle& b, const Particle& c, const Particle& d)
{
  typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

  helicity_angles<T_res> angles;

  /**
   * This function has been successfully tested against ROOTPWA's nBodyPhaseSpaceGen.
   */

  // helicity angle theta_1
  const T_res m2_24  = P.m2 + 2.0 * (a.m2 + b.m2 + c.m2 + d.m2) - (m2_12 + m2_14 + m2_23 + m2_34 + m2_13);
  const T_res m2_234 = m2_23 + m2_24 + m2_34 - (b.m2 + c.m2 + d.m2);

  const T_res p2_2 = mfct::breakup_momentum::p2(m2_12, a.m, b.m); // in R1
  const T_res E_2 = sqrt(b.m2 + p2_2);

  const T_res p2_34 = mfct::breakup_momentum::p2(P.m2, sqrt(m2_12), sqrt(m2_34)); // in P
  const T_res v2_R_1 = p2_34 / (m2_12 + p2_34);
  const T_res v_R_1 = sqrt(v2_R_1);
  const T_res gamma_R_1 = 1. / sqrt(1. - v2_R_1);

  const T_res p_34_in_R_1 = gamma_R_1 * (sqrt(p2_34) + v_R_1 * sqrt(m2_34 + p2_34));
  const T_res p2_34_in_R_1 = p_34_in_R_1 * p_34_in_R_1;
  const T_res E_34_in_R_1 = sqrt(m2_34 + p2_34_in_R_1);
  const T_res p_2_dot_p_34_in_R_1 = E_2 * E_34_in_R_1 - 0.5 * (m2_234 - b.m2 - m2_34);

  const T_res cos_theta_1 = p_2_dot_p_34_in_R_1 / sqrt(p2_2 * p2_34_in_R_1);
  angles.theta_1 = acos(cos_theta_1);


  // helicity angle theta_2 // same as theta_1, but particles 1 <-> 3 and 2 <-> 4 swapped
  const T_res m2_124 = m2_14 + m2_24 + m2_12 - (a.m2 + b.m2 + d.m2);

  const T_res p2_4 = mfct::breakup_momentum::p2(m2_34, c.m, d.m); // in R2
  const T_res E_4 = sqrt(d.m2 + p2_4);

  const T_res p2_12 = p2_34; // in P
  const T_res v2_R_2 = p2_12 / (m2_34 + p2_12);
  const T_res v_R_2 = sqrt(v2_R_2);
  const T_res gamma_R_2 = 1. / sqrt(1. - v2_R_2);

  const T_res p_12_in_R_2 = gamma_R_2 * (sqrt(p2_12) + v_R_2 * sqrt(m2_12 + p2_12));
  const T_res p2_12_in_R_2 = p_12_in_R_2 * p_12_in_R_2;
  const T_res E_12_in_R_2 = sqrt(m2_12 + p2_12_in_R_2);
  const T_res p_4_dot_p_12_in_R_2 = E_4 * E_12_in_R_2 - 0.5 * (m2_124 - d.m2 - m2_12);

  const T_res cos_theta_2 = p_4_dot_p_12_in_R_2 / sqrt(p2_4 * p2_12_in_R_2);
  angles.theta_2 = acos(cos_theta_2);


  // acoplanarity angle chi
  const T_res p2_1 = p2_2; // in R1
  const T_res E_1 = sqrt(a.m2 + p2_1); // in R1

  const T_res p2_1_in_P = pow(gamma_R_1 * (-cos_theta_1 * sqrt(p2_1) - E_1 * v_R_1), 2)
      + (1. - cos_theta_1*cos_theta_1) * p2_1;

  const T_res p2_2_in_P = pow(gamma_R_1 * ( cos_theta_1 * sqrt(p2_2) - E_2 * v_R_1), 2)
          + (1. - cos_theta_1*cos_theta_1) * p2_2;


  const T_res p2_3 = p2_4;
  const T_res E_3 = sqrt(c.m2 + p2_3);

  const T_res p2_3_in_P = pow(gamma_R_2 * (-cos_theta_2 * sqrt(p2_3) - E_3 * v_R_2), 2)
      + (1. - cos_theta_2*cos_theta_2) * p2_3;

  const T_res p2_4_in_P = pow(gamma_R_2 * ( cos_theta_2 * sqrt(p2_4) - E_4 * v_R_2), 2)
          + (1. - cos_theta_2*cos_theta_2) * p2_4;

  const T_res E_1_in_P = sqrt(a.m2 + p2_1_in_P);
  const T_res E_2_in_P = sqrt(b.m2 + p2_2_in_P);
  const T_res E_3_in_P = sqrt(c.m2 + p2_3_in_P);
  const T_res E_4_in_P = sqrt(d.m2 + p2_4_in_P);

  const T_res p_1_dot_p_2_in_P = E_1_in_P * E_2_in_P - 0.5 * (m2_12 - a.m2 - b.m2);
  const T_res p_1_dot_p_3_in_P = E_1_in_P * E_3_in_P - 0.5 * (m2_13 - a.m2 - c.m2);
  const T_res p_1_dot_p_4_in_P = E_1_in_P * E_4_in_P - 0.5 * (m2_14 - a.m2 - d.m2);
  const T_res p_2_dot_p_3_in_P = E_2_in_P * E_3_in_P - 0.5 * (m2_23 - b.m2 - c.m2);
  const T_res p_2_dot_p_4_in_P = E_2_in_P * E_4_in_P - 0.5 * (m2_24 - b.m2 - d.m2);
  const T_res p_3_dot_p_4_in_P = E_3_in_P * E_4_in_P - 0.5 * (m2_34 - c.m2 - d.m2);

  const T_res p_1_cross_p_2_dot_p_3_cross_p_4 = p_1_dot_p_3_in_P * p_2_dot_p_4_in_P
      - p_2_dot_p_3_in_P * p_1_dot_p_4_in_P;

  const T_res p_1_cross_p_2_abs2 = p2_1_in_P * p2_2_in_P - p_1_dot_p_2_in_P * p_1_dot_p_2_in_P;
  const T_res p_3_cross_p_4_abs2 = p2_3_in_P * p2_4_in_P - p_3_dot_p_4_in_P * p_3_dot_p_4_in_P;

  const T_res cos_chi = p_1_cross_p_2_dot_p_3_cross_p_4 / sqrt(p_1_cross_p_2_abs2 * p_3_cross_p_4_abs2);
  angles.chi = acos(cos_chi);


  return angles;
}


}
}

#endif
