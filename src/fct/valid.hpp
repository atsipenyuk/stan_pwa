#ifndef STAN_PWA__SRC__FCT__VALID_HPP
#define STAN_PWA__SRC__FCT__VALID_HPP

#include <cmath> // sqrt

#include <stan_pwa/src/structures/particles_def.hpp> 
// class particle

/*
 * Check whether we are in the energetically allowed region of the decay.
 * 
 * FUNCTIONS
 * valid(m2_ab, m2_bc, p, a, b, c) - Check m2_ab, m2_bc for p->abc decay
 * valid(m2_ab, m2_bc, m2_p, m2_a, m2_b, m2_c) - As above, other input type
 * valid_5d() // WRONG. To be corrected.
 * 
 */

namespace stan_pwa {
namespace fct {

  /**
   * bool valid(m2_ab, m2_bc, p, a, b, c)
   *
   * Determines whether we are in an energetically allowed
   * phase space region of the decay p-> a + b + c.
   *
   * m2_ab, m2_bc are the invariant square masses of a and b, b and c 
   * (scalars). p, a, b, c are decay particles, usually specified in
   * c_lib/structures/particles.hpp .
   *
   * Since we do not need the momenta of p, a, b, c, we can use their
   * masses as arguments - the overloaded function may be found below.
   * 
   */
  template <typename T>
  inline
  bool valid(const T &m2_ab, const T &m2_bc, 
             const particle &p, const particle &a, 
             const particle &b, const particle &c)
  {
    if ( (m2_ab < (a.m2 + b.m2 + 2. * sqrt(a.m2 * b.m2))) ||
         (m2_ab > (p.m2 + c.m2 - 2. * sqrt(p.m2 * c.m2)))   ) {
      return false;
    }

   typedef typename boost::math::tools::promote_args<T, double>::type T_res;

    T_res E_b = (m2_ab - a.m2 + b.m2) / 2. / sqrt(m2_ab);
    T_res E_c = (p.m2 - m2_ab - c.m2) / 2. / sqrt(m2_ab);
    T_res P_b = sqrt(E_b * E_b - b.m2);
    T_res P_c = sqrt(E_c * E_c - c.m2);

    if ((fabs(m2_bc - b.m2 - c.m2 - 2. * E_b * E_c)) <= (2. * P_b * P_c)){
      return true;
    }

    return false;
  }


  /*
   * Overloaded input arguments to allow particle masses.
   */
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline
  bool valid(const T0 &m2_ab, const T0 &m2_bc, 
             const T1 &m2_p, const T2 &m2_a, const T3 &m2_b, const T4 &m2_c)
  {
    if ( (m2_ab < (m2_a + m2_b + 2. * sqrt(m2_a * m2_b))) ||
         (m2_ab > (m2_p + m2_c - 2. * sqrt(m2_p * m2_c)))   ) {
      return false;
    }

    typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4,double>::type T_res;

    T_res E_b = (m2_ab - m2_a + m2_b) / 2. / sqrt(m2_ab);
    T_res E_c = (m2_p - m2_ab - m2_c) / 2. / sqrt(m2_ab);
    T_res P_b = sqrt(E_b * E_b - m2_b);
    T_res P_c = sqrt(E_c * E_c - m2_c);

    /* verbose debugging
    std::cout << "Valid range: from " 
	      << (E_b + E_c) * (E_b + E_c) - (P_b + P_c)
	      << " to "
	      << (E_b + E_c) * (E_b + E_c) - (P_b - P_c)
	      << "\n";
    */

    if ((fabs(m2_bc - m2_b - m2_c - 2. * E_b * E_c)) <= (2. * P_b * P_c)){
      return true;
    }

    return false;
  }


  /**
   * bool valid_5d(m2_12, ..., P, a, b, c)
   *
   * Determines whether we are in an energetically allowed
   * phase space region of the decay P -> a b c d.
   */
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  bool valid_5d(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
    const T3 &m2_34, const T4& m2_13,
    const particle &Parent,
    const particle &a, const particle &b,
    const particle &c, const particle &d)
  {
    typedef typename boost::math::tools::promote_args<T1,T3,T4>::type T_134;
    typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

    // 5D hypercube lower boundaries
    if ( m2_12 < pow(a.m + b.m, 2) ||
         m2_14 < pow(a.m + d.m, 2) ||
         m2_23 < pow(b.m + c.m, 2) ||
         m2_34 < pow(c.m + d.m, 2) ||
         m2_13 < pow(a.m + c.m, 2) ) {
      return false;
    }

    // 5D hypercube upper boundaries
    if ( m2_12 > pow(Parent.m - c.m - d.m, 2) ||
         m2_14 > pow(Parent.m - b.m - c.m, 2) ||
         m2_23 > pow(Parent.m - a.m - d.m, 2) ||
         m2_34 > pow(Parent.m - a.m - b.m, 2) ||
         m2_13 > pow(Parent.m - b.m - d.m, 2) ) {
      return false;
    }

    const T_res m2_24 = (Parent.m2 + 2.*(a.m2 + b.m2 + c.m2 + d.m2)) - m2_12 - m2_14 - m2_23 - m2_34 - m2_13;

    if ( sqrt(m2_12) + sqrt(m2_34) > Parent.m ||
         sqrt(m2_14) + sqrt(m2_23) > Parent.m ||
         sqrt(m2_13) + sqrt(m2_24) > Parent.m ) {
      return false;
    }
    



    const T0 M = m2_12;
    const T0 M2 = M*M;

    const T3 N = m2_34;
    const T3 N2 = N*N;

    const T_res P = m2_12 + m2_14 + m2_24 - a.m2 - b.m2 - d.m2; // m2_124
    const T_res P2 = P*P;

    const T_134 Q = m2_13 + m2_14 + m2_34 - a.m2 - c.m2 - d.m2; // m2_134;
    const T_134 Q2 = Q*Q;

    const T1 R = m2_14;
    const T1 R2 = R*R;

    const double m = c.m2;
    const double m2 = m*m;

    const double n = b.m2;
    const double n2 = n*n;
    
    const double p = a.m2;
    const double p2 = p*p;

    const double q = d.m2;
    const double q2 = q*q;

    const double r = Parent.m2; // E^2
    const double r2 = r*r;

    
    T_res B = (M2*Q2 + N2*P2 + M2*R2 + N2*R2 + P2*Q2) - 2.*(M2*Q*R + N2*P*R + M*N*R2 + M*P*Q2 + N*P2*Q)
        + 2.*(M*N*P*Q + M*N*P*R + M*N*Q*R + M*P*Q*R + N*P*Q*R) - 2.*(M2*Q*m + N2*P*n + M2*R*m + N2*R*n + M*Q2*q
            + N*P2*p + M*R2*r + N*R2*r + P2*Q*p + P*Q2*q) - 2.*(M*N*P*m + M*N*Q*n + M*P*R*p + N*Q*R*q + P*Q*R*r)
        + 2.*(M*N*R*(m + n - 2.*r) + M*P*Q*(m + p - 2.*q) + N*Q*P*(n + q - 2.*p) + Q*R*M*(q + r - 2.*m) + P*R*N*(p + r - 2.*n))
        + (M2*m2 + N2*n2 + P2*p2 + Q2*q2 + R2*r2) + 2.*(M*N*m*n + M*P*m*p + N*Q*n*q + P*R*p*r + Q*R*q*r)
        + 2.*(M*Q*(m*q + m*n + q*n + m*p + q*r - p*r) + N*P*(n*p + n*m + p*m + p*r + n*q - q*r) + M*R*(m*r + m*p + r*p + m*n + r*q - n*q)
            + N*R*(n*r + n*q + r*q + n*m + r*p - m*p) + P*Q*(p*q + p*r + q*r + p*m + q*n - m*n) ) - 2.*(M*m*(m*p + m*n + q*r - p*r - n*q + 2.*n*p)
            + N*n*(n*m + n*q + p*r - p*m - q*r + 2.*m*q) + P*p*(p*m + p*r + n*q - m*n - q*r + 2.*m*r) + Q*q*(q*n + q*r + m*p - m*n - p*r + 2.*n*r)
            + R*r*(r*p + r*q + m*n - m*p - n*q + 2.*p*q) ) + (m2*n2 + m2*p2 + n2*q2 + p2*r2 + q2*r2) - 2.*(m2*n*p + m*n2*q + m*p2*r + n*q2*r + p*q*r2)
            + 2.*(m*n*p*q + m*n*p*r + m*n*q*r + m*p*q*r + n*p*q*r);

  
    if (B < 0.)
      return true;

    return false;

  }

}
}
#endif
