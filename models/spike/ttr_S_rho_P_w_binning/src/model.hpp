#ifndef STAN_PWA__SRC__MODEL_HPP
#define STAN_PWA__SRC__MODEL_HPP

#include <vector>

#include <stan/math.hpp>
#include <boost/math/tools/promotion.hpp>

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/real.hpp>
#include <stan_pwa/src/structures/three_body_resonances.hpp>

/* MODEL_HPP
 * 
 * This model uses two different model functions:
 *  - f_gen is used for data generation (model-dependent);
 *  - f_fit is used for amplitude fitting (S wave is binned).
 *
 * The former is defined in this file. The latter is defined in
 * model_fit.hpp.
 */

const int NUM_RES=4; // Number of PWA resonances
const int NUM_VAR=2; // Number of independent masses (e.g.,2 for 3-body-decay)
const int NUM_BCKGR=1; // Number of background amplitudes

#include <stan_pwa/models/ttr_S_rho_P_w_binning/src/phase_space_gen.hpp>

namespace mcomplex = stan_pwa::complex;
namespace mreal = stan_pwa::real;
namespace mresonances = stan_pwa::resonances;

namespace stan {
  namespace math {

    /**
     * complex_scalar A_c(vector) - the PWA amplitudes are defined here.
     *
     * Takes the data vector y as an argument, returns the corresponding PWA 
     * amplitude (complex number). The number res_id tells, which resonance
     * to use.
     *
     * CAVEAT. All the S-wave resonances must be defined BEFORE all other
     * resonances.
     * 
     * @tparam T0__ Scalar type of the data vector
     */
    template <typename T0__>
    std::vector<typename boost::math::tools::promote_arg<T0__>::type>
    A_c(const int &res_id, const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        switch (res_id) {
        // Adjust this list according to your model
        case 1: return mresonances::f0_980.value_sym(y(0,0), y(1,0));

        case 2: return mresonances::f0_1370.value_sym(y(0,0), y(1,0));

        case 3: return mresonances::rho_770.value_sym(y(0,0), y(1,0));

        case 4: return mresonances::f2_1270.value_sym(y(0,0), y(1,0));

         // Dummy return. May be reworked to raise an error.
        default: {
	  std::cout << "Fatal error: Unknown resonance occured." << res_id;
	  return mcomplex::scalar::one(y(0,0)); // Dummy return.
          }
        }
    }


    /**
     * complex_scalar A_c_background(vector) - the PWA background amplitudes
     *   are defined here.
     *
     * Takes the data vector y as an argument, returns the corresponding 
     * background amplitude (complex number). The number res_id tells, which 
     * resonance to use.
     *
     * @tparam T0__ Scalar type of the data vector
     */
    template <typename T0__>
    std::vector<typename boost::math::tools::promote_arg<T0__>::type>
    A_c_background(const int &res_id, 
		   const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        switch (res_id) {
	// This resonance list must be adjusted manually
        //case 1: return mresonances::flat_D3pi.value(y(0,0), y(1,0));

        case 1: return mresonances::rho_770_bw_only.value_sym(y(0,0), y(1,0));

        default: {
            std::cout << "Fatal error: Unknown resonance occured.";
            return mcomplex::scalar::one(y(0,0)); // Dummy return
          }
        }
    }


    /**
     * complex_vector amplitude_vector(vector)
     *
     * Takes the data vector y as an argument, returns 
     * complex vector [A(1,y) ... A(NUM_RES, y)] of PWA amplitudes.
     */
    template <typename T0__>
    //inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {
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
     * vector background_vector(vector)
     *
     * Takes the data vector y as an argument, returns 
     * real vector [|A(1,y)|^2 ... |A(NUM_BCKGR, y)|^2] 
     + of PWA background amplitudes squared.
     */
    template <typename T0__>
    //inline
    Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic,1>
    background_vector(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0__>::type T2;

	Eigen::Matrix<T2, Eigen::Dynamic, 1> res(NUM_BCKGR);
        for (int i = 0; i < NUM_BCKGR; i++) {
	  res(i) = mcomplex::scalar::abs2(A_c_background(i+1, y));
        }
        return res;
    }


    /**
     *
     * double f_gen(vector A_y[2], vector theta[2],
     *                vector A_y_background_abs2, vector theta_background_abs2)
     *
     * Takes four comlex vectors, returns 
     *  |A_y * theta|^2 + A_y_background_abs2 * theta_background_abs2)
     *
     */
    template <typename T1, typename T2, typename T3, typename T4>
    typename boost::math::tools::promote_args<T1,T2,T3,T4>::type
    f_gen(const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& A_r,
	  const std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> >& theta,
	  const Eigen::Matrix<T3, Eigen::Dynamic, 1>& background_vector_,
	  const Eigen::Matrix<T4, Eigen::Dynamic, 1>& theta_background_abs2_) {

      //typename boost::math::tools::promote_args<T0,T1>::type res = 0;

      return 
	// Sum coherent amplitudes
	mcomplex::scalar::abs2(
          mcomplex::vector::sum(
	    mcomplex::vector::mult(A_r, theta)))

	// Add background
	// Use STAN vector dot multiplication.
	+ mreal::vector::mult(background_vector_,theta_background_abs2_);
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

#include <stan_pwa/models/ttr_S_rho_P_w_binning/src/model_fit.hpp>

#endif













