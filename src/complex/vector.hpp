#ifndef STAN_PWA__SRC__COMPLEX__VECTOR_HPP
#define STAN_PWA__SRC__COMPLEX__VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

/*
 *  Introduce complex vector operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    See stan_pwa/src/complex/complex.hpp
 *
 *  FUNCTIONS
 *    complex_vector mult(complex_vector, complex_vector)
 *    complex_vector mult(vector, complex_vector)
 *    complex_vector sum(complex_vector, complex_vector)
 */

namespace stan_pwa {
namespace complex {
  namespace vector {

    /**
     * complex_vector mult(complex_vector, complex_vector)
     *
     * Multiplication of two complex vectors (represented by a 2d array
     * of vectors - 1st element of the array is the real part of the complex 
     * vector, 2nd element is the imag. part of the complex vector).
     *
     * @tparam T Scalar vector type
     */
    template <typename T0, typename T1>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0,T1>::type, Eigen::Dynamic, 1> >
    mult(const std::vector<Eigen::Matrix<T0,Eigen::Dynamic,1> > &v1, 
            const std::vector<Eigen::Matrix<T1,Eigen::Dynamic,1> > &v2) {

        // check size of v1 and v2
        int v_len = v1[0].rows();
        if (v_len != v2[0].rows()) {
            std::cout << "Arugment size mismatch in complex::vector::mult.";
        }

	typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
        std::vector<Eigen::Matrix<T_res, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T_res, Eigen::Dynamic, 1> (v_len)));
        for (int i = 0; i < v_len; i++) {
            res[0](i) = v1[0](i) * v2[0](i) - v1[1](i) * v2[1](i);
            res[1](i) = v1[0](i) * v2[1](i) + v1[1](i) * v2[0](i);
        }
        return res;
    }



    /**
     * complex_vector mult(complex_vector, complex_vector)
     *
     * Multiplication of real vector with a complex vector 
     * (represented by a 2d array of vectors - 1st element 
     * of the array is the real part of the complex vector, 
     * 2nd element is the imag. part of the complex vector).
     *
     * @tparam T Scalar vector type
     */
    template <typename T0, typename T1>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0,T1>::type, Eigen::Dynamic, 1> >
    mult(const Eigen::Matrix<T0,Eigen::Dynamic,1> &v1, 
         const std::vector<Eigen::Matrix<T1,Eigen::Dynamic,1> > &v2) {

        // check size of v1 and v2
        int v_len = v1.rows();
        if (v_len != v2[0].rows()) {
            std::cout << "Arugment size mismatch in complex::vector::mult.";
        }

	typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
        std::vector<Eigen::Matrix<T_res, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T_res, Eigen::Dynamic, 1> (v_len)));
        for (int i = 0; i < v_len; i++) {
            res[0](i) = v1(i) * v2[0](i);
            res[1](i) = v1(i) * v2[1](i);
        }
        return res;
    }


    /**
     * complex_number sum(complex_vector)
     *
     * Sum of all the elements of a complex vector.
     *
     * @tparam T Scalar matrix type
     */
    template <typename T>
    inline
    std::vector<T> sum(const std::vector<Eigen::Matrix<T,Eigen::Dynamic,1> > &v) {

        std::vector<T> res(2);
        res[0] = 0.0;
        res[1] = 0.0;

        for (int i = 0; i < v[0].rows(); i++) {
            res[0] += v[0](i);
            res[1] += v[1](i);
        }

        return res;
    }
  }
}
}
#endif
