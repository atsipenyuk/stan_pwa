#ifndef STAN_PWA__SRC__MODEL_HPP
#define STAN_PWA__SRC__MODEL_HPP

#include <vector>

#include <stan/math.hpp>
#include <boost/math/tools/promotion.hpp>

#include <stan_pwa/src/complex.hpp>
#include <stan_pwa/src/structures.hpp>

#include <stan_pwa/src/fct.hpp>
#include <stan_pwa/models/D4pi/src/model_resonances.hpp>


namespace mfct = stan_pwa::fct;
namespace mcomplex = stan_pwa::complex;
namespace mresonances = stan_pwa::resonances;

extern const int NUM_RES;
extern const int NUM_VAR;


namespace stan {
  namespace math {

    /**
     * complex_vector amplitude_vector(vector)
     *
     * Takes the data vector y as an argument, returns 
     * complex vector [A(1,y) ... A(NUM_RES, y)] of PWA amplitudes.
     */
    template <typename T0__>
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1> >
    amplitude_vector(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0__>::type T2;

        assert(y.rows() == NUM_VAR && y.cols() == 1);

        static struct_model_resonances<T0__>* model_resonances = new struct_model_resonances<T0__>();
        model_resonances->updateSymmetrized(y(0),y(1),y(2),y(3),y(4));

        // Somewhat convoluted initialization of the return
        std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T2,Eigen::Dynamic,1> (NUM_RES)));
        for (int i = 0; i < NUM_RES; i++) {
            std::vector<T0__> tmp;
            tmp = model_resonances->A_c(i);
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

      return mcomplex::scalar::abs2(
               mcomplex::vector::sum(
                 mcomplex::vector::mult(A_r, theta)));
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
          // Keep only the real palt of the product; imaginary part
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
