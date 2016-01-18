#ifndef MESON_DECA__LIB__C_LIB__FCT__P_R1D_R2CD_HPP
#define MESON_DECA__LIB__C_LIB__FCT__P_R1D_R2CD_HPP

#include <cmath> // sqrt
#include <math.h> // isnan

#include "breakup_momentum.hpp"


namespace mfct = stan_pwa::fct;

namespace stan_pwa {
namespace fct {

template <typename T0>
struct theta_z_values {
  T0 cos2_theta_1;
  T0 z2_1;
  T0 cos2_theta_2;
  T0 z2_2;
};

/**
 * Calculate theta1 (angle between c and d in R1), theta2 (angle between b and c in R2)
 * and z1 (p_d / sqrt(s)) and z2 (p_c / sqrt(s))
 * from m2_12, m2_14, m2_23, m2_34, m2_13.
 *
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4>
theta_z_values<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type>
P_R1d_R2cd_theta_z(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
    const T3 &m2_34, const T4& m2_13,
    const particle& P, const particle& a, const particle& b, const particle& c, const particle& d)
{
  typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
  typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

  theta_z_values<T_res> v;

  // Invariant square mass of particles a,b,c together
  const T_123 m2_123 = m2_12 + m2_13 + m2_23 - a.m2 - b.m2 - c.m2;

  const T_123 m_123 = sqrt(m2_123);

  // Zemach tensors
  // Calculate transformed variables z2, cos2_theta for
  // the decay D-> R_1 d -> R_2 c d
  // in the rest frame of R_1
  const T_123 p2_c = mfct::breakup_momentum::p2(m2_123, sqrt(m2_12), c.m); // in R1
  assert(p2_c >= 0.);

  const T_123 E_c = sqrt(c.m2 + p2_c); // in R1

  const T_123 p2_d_rest_frame_of_P = mfct::breakup_momentum::p2(P.m2, // in P
      sqrt(m2_123), d.m);
  assert(p2_d_rest_frame_of_P >= 0.);
  const T_123 p_d_rest_frame_of_P = sqrt(p2_d_rest_frame_of_P); // in P

  const T_123 E_d_rest_frame_of_P = sqrt(d.m2 + p2_d_rest_frame_of_P); // in P
  // R_1 and d have the same abs. momenta |p2| in rest frame of P
  const T_123 E_R_1_rest_frame_of_P = sqrt(m2_123 + p2_d_rest_frame_of_P); // in P

  // Compute p2_d in the rest frame of R_1 by
  // performing a Lorentz boost in the direction -v_R_1
  const T_123 v_R_1 = p_d_rest_frame_of_P / E_R_1_rest_frame_of_P; // in P // p2_d_rest_frame_of_P == p2_R1_rest_frame_of_P
  const T_123 gamma = 1.0 / sqrt(1.0 - v_R_1 * v_R_1); // in P

  const T_123 p_d = gamma * ( p_d_rest_frame_of_P + E_d_rest_frame_of_P * v_R_1 ); // in R1 // in direction of -v_R_1
  const T_123 p2_d = p_d*p_d; // in R1
  const T_123 E_d = sqrt(d.m2 + p2_d); // in R1

  // m2_34 = (E_c + E_d)^2 - (p2_c + 2*p_c*p_d + p2_d)
  const T_res p_c_dot_p_d = (-0.5) * (m2_34 - c.m2 - d.m2 - 2.0 * E_c * E_d); // in R1
  //const T_res p_c_dot_p_d = 0.5 * (pow(E_c + E_d, 2) - m2_34 - p2_c - p2_d); // gives same result

  v.cos2_theta_1 = p_c_dot_p_d * p_c_dot_p_d / p2_c / p2_d; // in R1

  const T_123 s = m2_123 + d.m2 + 2.0 * m_123 * E_d; // sqrt(s) = E_abcd in R1
  v.z2_1 = p2_d / s;




  // Calculate transformed variables z2, cos2_theta for
  // the decay R_1 -> R_2 c -> a b c
  // in the rest frame of R_2

  // p2_c above is calculated in the rest frame of R_1.
  // We want it in the rest frame of R_2, so we perform a
  // Lorentz boost again.
  const T_123 p_c = sqrt(p2_c); // in R1
  const T_123 E_R_2 = sqrt(m2_12 + p2_c); // in R1

  const T_123 v_R_2 = p_c / E_R_2; // in R1
  const T_123 gamma_2 = 1.0 / sqrt(1.0 - v_R_2 * v_R_2); // in R1
  const T_123 p_c_rest_frame_of_R_2 = gamma_2 * (p_c + E_c * v_R_2); // in R2
  assert(p_c_rest_frame_of_R_2 >= 0.);
  const T_123 p2_c_rest_frame_of_R_2 = p_c_rest_frame_of_R_2 * p_c_rest_frame_of_R_2; // in R2

  const T_123 E_c_rest_frame_of_R_2 = sqrt(c.m2 + p2_c_rest_frame_of_R_2); // in R2

  const T_123 p2_b = mfct::breakup_momentum::p2(m2_12, a.m, b.m); // in R2
  assert(p2_b >= 0.);
  const T_123 E_b = sqrt(b.m2 + p2_b); // in R2

  const T_res p_b_dot_p_c = (-0.5) * (m2_23 - b.m2 - c.m2 - 2.0 * E_b * E_c_rest_frame_of_R_2); // in R2
  v.cos2_theta_2 = p_b_dot_p_c * p_b_dot_p_c / p2_b / p2_c_rest_frame_of_R_2;

  const T0 m_12 = sqrt(m2_12);

  const T_123 s_2 = m2_12 + c.m2 + 2.0 * m_12 * E_c_rest_frame_of_R_2;
  v.z2_2 = p2_c / s_2;

  assert(v.z2_2 >= 0.);


  return v;
}


}
}

#endif
