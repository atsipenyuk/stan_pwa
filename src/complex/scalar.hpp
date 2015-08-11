#ifndef STAN_PWA__SRC__COMPLEX__SCALAR_HPP
#define STAN_PWA__SRC__COMPLEX__SCALAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>


/*
 *  Introduce complex number operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    See stan_pwa/src/complex/complex.hpp
 *
 *  FUNCTIONS
 *    complex_scalar abs2(complex_scalar)
 *    complex_scalar add(complex_scalar, complex_scalar)
 *    complex_scalar complex(scalar, scalar)
 *    complex_scalar inverse(complex_scalar)
 *    complex_scalar one(scalar)
 *    complex_scalar one_i(scalar)
 *    complex_scalar mult(complex_scalar, complex_scalar)
 *    complex_scalar mult(scalar, complex_scalar)
 *    complex_scalar subtract(complex_scalar, complex_scalar)
 */


namespace stan_pwa {
namespace complex {
  namespace scalar {

    /**
     * scalar abs2(complex_scalar)
     *
     * Square magnitude of a complex number.
     * Takes the vector (a, b), returns a**2 + b**2.
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline T
    abs2(const std::vector<T> &v) {
        return v[0] * v[0] + v[1] * v[1];
    }


    /**
     * complex_scalar add(complex_scalar, complex_scalar)
     *
     * Sum of two complex numbers.
     *
     * @tparam T Scalar type
     */
    template <typename T0, typename T1>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1>::type> 
    add(const std::vector<T0> &v1, const std::vector<T1> &v2) {
        typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
        std::vector<T_res> res(2);
        res[0] = v1[0] + v2[0];
        res[1] = v1[1] + v2[1];
        return res;
    }


    /**
     * complex_scalar complex(scalar, scalar)
     *
     * Complex number constructor.
     * Takes two scalars a,b returns the vector (a, b).
     *
     * @tparam T0, T1 Scalar type
     */
    template <typename T0, typename T1>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1>::type>
    complex(const T0& re, const T1& im) {
      typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
      std::vector<T_res> res(2);
      res[0] = re;
      res[1] = im;
      return res;
    }


    /**
     * complex_scalar inverse(complex_scalar)
     *
     * Inverse of a complex number.
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> inverse(const std::vector<T>& y) {

        std::vector<T> res(2);
        T norm = y[0] * y[0] + y[1] * y[1];

        res[0] = y[0] / norm;
        res[1] = -y[1] / norm;
        return res;
    }


    /**
     * complex_scalar one(scalar)
     *
     * Returns the complex number 1.
     *
     * @param y Dummy scalar variable
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> one(const T& y) {
        std::vector<T> res(2);
        res[0] = 1.0;
        res[1] = 0.0;
        return res;
    }


    /**
     * complex_scalar one(scalar)
     *
     * Returns the imaginary unit 1*i.
     *
     * @param y Dummy scalar variable
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> one_i(const T& y) {
        std::vector<T> res(2);
        res[0] = 0.0;
        res[1] = 1.0;
        return res;
    }


    /**
     * complex_scalar mult(complex_scalar, complex_scalar) [Overl. op.]
     *
     * Multiplication of two complex numbers.
     *
     * @tparam T0,T1 Scalar types
     */
    template <typename T0, typename T1>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1>::type> 
    mult(const std::vector<T0> &v1, const std::vector<T1> &v2) {
      typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
      std::vector<T_res> res(2);
      res[0] = v1[0] * v2[0] - v1[1] * v2[1];
      res[1] = v1[1] * v2[0] + v1[0] * v2[1];
      return res;
    }


    /**
     * complex_scalar mult(scalar, complex_scalar) [Overloaded operator]
     *
     * Multiplication of a real and a complex number.
     *
     * @tparam T0,T1 Scalar types
     */
    template <typename T0, typename T1>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1>::type> 
    mult(const T0 &v1, const std::vector<T1> &v2) {
      typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
      std::vector<T_res> res(2);
      res[0] = v1 * v2[0];
      res[1] = v1 * v2[1];
      return res;
    }


    /**
     * complex_scalar subtract(complex_scalar, complex_scalar)
     *
     * Subtraction of two complex numbers.
     *
     * @tparam T Scalar type
     */
    template <typename T0, typename T1>
    inline
    std::vector<typename boost::math::tools::promote_args<T0,T1>::type> 
    subtract(const std::vector<T0> &v1, const std::vector<T1> &v2) {
        typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
        std::vector<T_res> res(2);
        res[0] = v1[0] - v2[0];
        res[1] = v1[1] - v2[1];
        return res;
    }

  }
}
}

#endif
