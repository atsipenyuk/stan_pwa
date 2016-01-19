#ifndef STAN_PWA__SRC__MODEL_HPP
#define STAN_PWA__SRC__MODEL_HPP

#include <vector>
#include <stan/math.hpp>
#include <boost/math/tools/promotion.hpp>

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/structures.hpp>

namespace mc = stan_pwa::complex;
namespace mresonances = stan_pwa::resonances;

// These variables should be adjusted manually
const int NUM_RES=2; // Number of PWA resonances
const int NUM_VAR=2; // Number of independent masses (e.g., 2 for 3-body-decay)

namespace stan {
  namespace math {

    /**
     * complex_scalar A_c(int res_id, vector y)
     *
     * Takes the data vector y as an argument, returns the corresponding PWA 
     * amplitude (complex number). The number res_id tells, which resonance
     * to use.
     *
     */
    template <typename T>
    inline
    std::vector<typename boost::math::tools::promote_arg<T>::type>
    A_c(const int &res_id, const Eigen::Matrix<T, Eigen::Dynamic,1>& y) {

        switch (res_id) {
	// This resonance list must be adjusted manually
        case 1: return mresonances::toy0_1000.value_sym(y(0,0), y(1,0));
	case 2: return mresonances::toy0_1200.value_sym(y(0,0), y(1,0));

        default: {
            std::cout << "Fatal error: Unknown resonance occured.";
            return mc::scalar::one(y(0,0)); // Dummy return
        }
        }
    }


    /**
     * complex_vector amplitude_vector(vector)
     *
     * Takes the data vector y as an argument, returns 
     * complex vector [A(1,y) ... A(NUM_RES, y)] of PWA amplitudes.
     */
    template <typename T0>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0>::type T2;

        // Somewhat convoluted initialization of the return
        std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T2,Eigen::Dynamic,1> (NUM_RES)));
        for (int i = 0; i < NUM_RES; i++) {
            std::vector<T0> tmp;
            tmp = A_c(i+1, y);
            res[0](i) = tmp[0];
            res[1](i) = tmp[1];
        }
        return res;
    }



    /**
     *
     * double f_genfit(vector A_y[2], vector theta[2]
     *
     * Takes two comlex vectors, returns |A_y * theta|^2
     *
     */
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& A_r,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta) {

      //typename boost::math::tools::promote_args<T0,T1>::type res = 0;

      return mc::scalar::abs2(
               mc::vector::sum(
                 mc::vector::mult(A_r, theta)));
    }


    /**
     *
     * double norm(vector theta[2], matrix I[2])
     *
     * Takes complex vector theta and complex matrix I,
     * returns conj(theta)' * I * theta.
     */
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    norm(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& theta,
         const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> >& I) {

      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < NUM_RES; i++) {
      for (int j = 0; j < NUM_RES; j++) {
        // Complex multiplication
        tmp[0] = I[0](i,j) * theta[0](j) - I[1](i,j) * theta[1](j);
        tmp[1] = I[0](i,j) * theta[1](j) + I[1](i,j) * theta[0](j);
          // Keep only the real part of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
        res = res + theta[0](i) * tmp[0] + theta[1](i) * tmp[1];
        }
      }

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
  }
}

#endif
