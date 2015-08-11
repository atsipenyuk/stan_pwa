#ifndef STAN_PWA__SRC__REAL__VECTOR_HPP
#define STAN_PWA__SRC__REAL__VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

/*
 *  Introduce real-valued vector operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    See stan_pwa/src/real.hpp
 *
 *  FUNCTIONS
 *    scalar mult(vector, vector)
 */

namespace stan_pwa {
namespace real {
  namespace vector {

    /**
     * scalar mult(vector, vector)
     *
     * Dot products of two vectors.
     *
     * @tparam T Scalar vector type
     */
    template <typename T0, typename T1>
    inline
    typename boost::math::tools::promote_args<T0,T1>::type
    mult(const Eigen::Matrix<T0,Eigen::Dynamic,1> &v1, 
            const Eigen::Matrix<T1,Eigen::Dynamic,1> &v2) {

        // check size of v1 and v2
        int v_len = v1.rows();
        if (v_len != v2.rows()) {
            std::cout << "Arugment size mismatch in real::vector::mult.";
        }

	typedef typename boost::math::tools::promote_args<T0,T1>::type T_res;
        T_res res = 0;
        for (int i = 0; i < v_len; i++) {
	  res =  res + v1(i) * v2(i);
        }
        return res;
    }

  }
}
}
#endif
