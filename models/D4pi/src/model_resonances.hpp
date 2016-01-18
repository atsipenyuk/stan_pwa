#ifndef STAN_PWA__SRC__MODEL_RESONANCES_HPP
#define STAN_PWA__SRC__MODEL_RESONANCES_HPP

#include <vector>

namespace mfct = stan_pwa::fct;
namespace mcomplex = stan_pwa::complex;
namespace mres = stan_pwa::resonances;


// Predefined constants for our model

// These variables should be adjusted manually
const int NUM_RES=6; // Number of PWA resonances
const int NUM_VAR=5; // Number of independent masses(e.g., 2 for 3-body-decay)
//const int NUM_BCKGR=2; // Number of background amplitudes

const stan_pwa::particle P = stan_pwa::particles::D0;
const stan_pwa::particle a = stan_pwa::particles::pi;
const stan_pwa::particle b = stan_pwa::particles::pi;
const stan_pwa::particle c = stan_pwa::particles::pi;
const stan_pwa::particle d = stan_pwa::particles::pi;


template <typename T>
class struct_model_resonances {

private:
  // Invariant mass squares
  T m2_12_, m2_14_, m2_23_, m2_34_, m2_13_;

  // Declare here all the derived members such as helicity angles
  // and Breit-Wigner function values that should be used several times
  //[...DECLARE DERIVED MEMBERS... theta, f, etc]
  bool valid_;
  std::vector< mfct::helicity_angles<T> > helicity_angles_sym; // helicity angles theta1, theta2, acoplanarity angle chi
  std::vector< mfct::theta_z_values<T> > theta_z_values_sym;


public:
  // Constructor - does nothing
  struct_model_resonances() :
    m2_12_(0), m2_14_(0), m2_23_(0), m2_34_(0), m2_13_(0),
    valid_(false)
  {
    helicity_angles_sym.resize(4);
    theta_z_values_sym.resize(4);
  };


private:
  // Update variables and derived variables (such as helicity angles)
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  void update(int i, const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
      const T3 &m2_34, const T4& m2_13)
  {
    // Write here how all the derived members are calculated

    helicity_angles_sym[i] = mfct::P_V1V2_angles(m2_12, m2_14, m2_23, m2_34, m2_13,
        P, a, b, c, d);

    theta_z_values_sym[i] = mfct::P_R1d_R2cd_theta_z(m2_12, m2_14, m2_23, m2_34, m2_13,
        P, a, b, c, d);

  };


public:
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  void updateSymmetrized(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
      const T3 &m2_34, const T4& m2_13)
  {
    m2_12_ = m2_12;
    m2_14_ = m2_14;
    m2_23_ = m2_23;
    m2_34_ = m2_34;
    m2_13_ = m2_13;

    valid_ = mfct::valid_5d(m2_12_, m2_14_, m2_23_, m2_34_, m2_13_,
        P, a, b, c, d);

    if (valid_) {
      update(0, m2_12_, m2_14_, m2_23_, m2_34_, m2_13_);
      update(1, m2_23_, m2_34_, m2_12_, m2_14_, m2_13_); // 1 <-> 3
      update(2, m2_14_, m2_12_, m2_34_, m2_23_, m2_13_); // 2 <-> 4
      update(3, m2_34_, m2_23_, m2_14_, m2_12_, m2_13_); // 1 <-> 3, 2 <-> 4
    }

  }



  /**
   * complex_scalar A_c(vector)
   *
   * Takes the data vector y as an argument, returns the corresponding PWA 
   * amplitude (complex number). The number res_id tells, which resonance
   * to use.
   *
   * @tparam T0__ Scalar type of the data vector
   */
  //inline
  std::vector<T>
  A_c(const int &res_id) {
    
    if (!valid_)
      return std::vector<T>(2, 0.0);

    switch (res_id) {
    // This resonance list must be adjusted manually

    //  D->a1 + something
    case  0: return mres::D_a_rho_S_wave.value_sym(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, theta_z_values_sym);
    case  1: return mres::D_a_rho_D_wave.value_sym(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, theta_z_values_sym);

    case  2: return mres::D_a_sigma.value_sym(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, theta_z_values_sym);

    // D->rho+rho
    case  3: return mres::D_rho_rho.value_sym_parallel(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case  4: return mres::D_rho_rho.value_sym_perpendicular(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case  5: return mres::D_rho_rho.value_sym_longitudinal(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
      
    /*case  6: return mres::D_omega_omega.value_sym_parallel(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case  7: return mres::D_omega_omega.value_sym_perpendicular(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case  8: return mres::D_omega_omega.value_sym_longitudinal(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);

    case  9: return mres::D_rho_omega.value_sym_parallel(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case 10: return mres::D_rho_omega.value_sym_perpendicular(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
    case 11: return mres::D_rho_omega.value_sym_longitudinal(m2_12_,m2_14_,m2_23_,m2_34_,m2_13_, helicity_angles_sym);
*/
    // 3 body: D-> R + pi + pi
    /*case 12: // f_0(980)
    case 13: // f_2(1270)
    case 14: // sigma*/
      
      
    default: {
      std::cout << "Fatal error: Unknown resonance occured.";
      return std::vector<T>(2, 0.0);
    }
    }
  }


  void testfunc() {std::cout<<"test\n";};
};

//template <typename T>
//struct_model_resonances<T> MakeStruct_model_resonances(T t) { return struct_model_resonances(t); }





#endif
