#ifndef STAN_PWA__SRC__COMPLEX__MATRIX_HPP
#define STAN_PWA__SRC__COMPLEX__MATRIX_HPP

#include <stan/math.hpp>
#include <vector>

/*
 *  Introduce complex matrix operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    See meson_deca/lib/c_lib/complex/complex.hpp
 *
 *  FUNCTIONS
 *    complex_matrix ct(complex_matrix)
 */


namespace stan_pwa {
namespace complex {
  namespace matrix {

    /**
     * complex_matrix ct(complex_matrix)
     *
     * Returns the complex conjugate of the transposed matrix.
     *
     */
    template <typename T>
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T>::type, Eigen::Dynamic, Eigen::Dynamic> >
    ct(std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T>::type, Eigen::Dynamic, Eigen::Dynamic> > &m) {

      // Initialize return - I know, looks pretty horrible
      std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T>::type, Eigen::Dynamic, Eigen::Dynamic> > res(2,  (Eigen::Matrix<T, Eigen::Dynamic, 1> (m[0].cols(), m[0].rows())));

      res[0] = m[0].transpose();      // Transpose the real part
      res[1] = -m[1].transpose();     // Transpose and conj the im part

      return res;
    };

  }
}
}
#endif
