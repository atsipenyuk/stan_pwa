#ifndef MESON_DECA__LIB__C_LIB__MODEL_HPP
#define MESON_DECA__LIB__C_LIB__MODEL_HPP

#include <vector>

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>

#include <meson_deca/lib/c_lib/complex.hpp>
#include <meson_deca/lib/c_lib/real.hpp>
#include <meson_deca/lib/c_lib/structures/resonances.hpp>


// These variables should be adjusted manually
const int NUM_RES=7; // Number of PWA resonances
const int NUM_VAR=2; // Number of independent masses (e.g., 2 for 3-body-decay)
const int NUM_BCKGR=2; // Number of background amplitudes

namespace stan {
  namespace math {

    /**
     * complex_scalar A_c(vector)
     *
     * Takes the data vector y as an argument, returns the corresponding PWA 
     * amplitude (complex number). The number res_id tells, which resonance
     * to use.
     *
     * @tparam T0__ Scalar type of the data vector
     */
    template <typename T0__>
    //inline
    std::vector<typename boost::math::tools::promote_arg<T0__>::type>
    A_c(const int &res_id, const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        switch (res_id) {
	// This resonance list must be adjusted manually
        case 1: return resonances::flat_D3pi.value(y(0,0), y(1,0));

        case 2: return resonances::f0_980.value_sym(y(0,0), y(1,0));

        case 3: return resonances::f0_600.value_sym(y(0,0), y(1,0));

        case 4: return resonances::f0_1370.value_sym(y(0,0), y(1,0));

        case 5: return resonances::f0_1500.value_sym(y(0,0), y(1,0));

        case 6: return resonances::rho_770.value_sym(y(0,0), y(1,0));

        case 7: return resonances::f2_1270.value_sym(y(0,0), y(1,0));

        default: {
            std::cout << "Fatal error: Unknown resonance occured.";
            return complex::scalar::one(y(0,0)); // Dummy return
        }
        }
    }


    /**
     * complex_scalar A_c_background(vector)
     *
     * Takes the data vector y as an argument, returns the corresponding 
     * background amplitude (complex number). 
     * The number res_id tells, which resonance
     * to use.
     *
     * @tparam T0__ Scalar type of the data vector
     */
    template <typename T0__>
    //inline
    std::vector<typename boost::math::tools::promote_arg<T0__>::type>
    A_c_background(const int &res_id, const Eigen::Matrix<T0__, 
		   Eigen::Dynamic,1>& y) {

        switch (res_id) {
	// This resonance list must be adjusted manually
        case 1: return resonances::flat_D3pi.value(y(0,0), y(1,0));

        case 2: return resonances::rho_770.value_sym(y(0,0), y(1,0));

        default: {
            std::cout << "Fatal error: Unknown resonance occured.";
            return complex::scalar::one(y(0,0)); // Dummy return
        }
        }
    }


    /**
     * complex_vector A_cv(vector)
     *
     * Takes the data vector y as an argument, returns 
     * complex vector [A(1,y) ... A(NUM_RES, y)] of PWA amplitudes.
     */
    template <typename T0__>
    //inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1> >
    A_cv(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0__>::type T2;

        // Somewhat convoluted initialization of the return
        std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T2,Eigen::Dynamic,1> (NUM_RES)));
        for (int i = 0; i < NUM_RES; i++) {
            std::vector<T0__> tmp;
            tmp = A_c(i+1, y);
            res[0](i) = tmp[0];
            res[1](i) = tmp[1];
        }
        return res;
    }


    /**
     * vector A_v_background_abs2(vector)
     *
     * Takes the data vector y as an argument, returns 
     * real vector [|A(1,y)|^2 ... |A(NUM_BCKGR, y)|^2] 
     + of PWA background amplitudes squared.
     */
    template <typename T0__>
    //inline
    Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic,1>
    A_v_background_abs2(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0__>::type T2;

	Eigen::Matrix<T2, Eigen::Dynamic, 1> res(NUM_BCKGR);
        for (int i = 0; i < NUM_BCKGR; i++) {
	  res(i) = complex::scalar::abs2(A_c_background(i+1, y));
        }
        return res;
    }


    /**
     *
     * double f_model(vector A_y[2], vector theta[2],
     *                vector A_y_background_abs2, vector theta_background_abs2)
     *
     * Takes four comlex vectors, returns 
     *  |A_y * theta|^2 + A_y_background_abs2 * theta_background_abs2)
     *
     */
    template <typename T1, typename T2, typename T3, typename T4>
    typename boost::math::tools::promote_args<T1,T2,T3,T4>::type
    f_model(const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& A_r,
	    const std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> >& theta,
	    const Eigen::Matrix<T3, Eigen::Dynamic, 1>& A_v_background_abs2_,
	    const Eigen::Matrix<T4, Eigen::Dynamic, 1>& theta_background_abs2_) {

      //typename boost::math::tools::promote_args<T0,T1>::type res = 0;

      return 
	// Sum coherent amplitudes
	complex::scalar::abs2(
          complex::vector::sum(
	    complex::vector::mult(A_r, theta)))

	// Add background
	+ real::vector::mult(A_v_background_abs2_,theta_background_abs2_);
    }


    /**
     *
     * double Norm(vector theta[2], matrix I[2],
     *             vector theta_background_abs2_, vector I_background_abs2_)
     *
     * Takes complex vector theta and complex matrix I, real vectors
     * theta_background_abs2_, I_background_abs2_,
     * returns 
     *   conj(theta)' * I * theta + theta_background_abs2_ * I_background_abs2_.
     */
    template <typename T0, typename T1, typename T2, typename T3>
    typename boost::math::tools::promote_args<T0,T1,T2,T3>::type
    Norm(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& theta,
         const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> >& I,
	 const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta_background_abs2_,
	 const Eigen::Matrix<T3, Eigen::Dynamic, 1>& I_background_abs2_) {

      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < NUM_RES; i++) {
      for (int j = 0; j < NUM_RES; j++) {
        // Complex multiplication
        tmp[0] = I[0](i,j) * theta[0](j) - I[1](i,j) * theta[1](j);
        tmp[1] = I[0](i,j) * theta[1](j) + I[1](i,j) * theta[0](j);
          // Keep only the real palt of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
        res = res + theta[0](i) * tmp[0] + theta[1](i) * tmp[1];
        }
      }

      res = res + real::vector::mult(theta_background_abs2_,I_background_abs2_);

      return res;
    }


    /**
     *
     * [double, double] = Norm(vector theta[2], matrix I[2],
     *                         vector theta_background_abs2_, 
     *                         vector I_background_abs2_)
     *
     * Previous function is overloaded here to return not the
     *   total integral, but the coherent and incoherent integrals separately.
     */
    template <typename T0, typename T1, typename T2, typename T3>
    Eigen::Matrix<typename boost::math::tools::promote_args<T0,T1,T2,T3>::type, Eigen::Dynamic, 1>
    Norm(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& theta,
         const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> >& I,
	 const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta_background_abs2_,
	 const Eigen::Matrix<T3, Eigen::Dynamic, 1>& I_background_abs2_) {

      typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
      Eigen::Matrix<T_res, Eigen::Dynamic, 1> res(2,0);
      // I * theta holder, real and imaginary part
      T_res tmp[2];

      for (int i = 0; i < NUM_RES; i++) {
      for (int j = 0; j < NUM_RES; j++) {
        // Complex multiplication
        tmp[0] = I[0](i,j) * theta[0](j) - I[1](i,j) * theta[1](j);
        tmp[1] = I[0](i,j) * theta[1](j) + I[1](i,j) * theta[0](j);
          // Keep only the real palt of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
        res(0) = res(0) + theta[0](i) * tmp[0] + theta[1](i) * tmp[1];
        }
      }

      res(1) = res(1) + real::vector::mult(theta_background_abs2_,I_background_abs2_);

      return res;
    }


    /**
     * int num_resonances()
     *
     * Returns the number of resonances
     */
    inline int num_resonances() {
        return NUM_RES;
    }

    /**
     * int num_variables()
     *
     * Returns the number of variables
     */
    inline int num_variables() {
        return NUM_VAR;
    }



    /**
     * int num_background()
     *
     * Returns the number of background amplitudes.
     */
    inline int num_background() {
        return NUM_BCKGR;
    }
  }
}

#endif
