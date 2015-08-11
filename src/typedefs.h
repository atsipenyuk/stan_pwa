#ifndef PWA_STAN__SRC__TYPEDEFS_H
#define PWA_STAN__SRC__TYPEDEFS_H

// These types are quite useful, so I define them in both namespaces

namespace stan_pwa {

  ///> Complex variable type ('templated typename')
  template <typename T>
    using C_t = typename std::vector<typename boost::math::tools::promote_arg<T>::type>;

  ///> Complex vector type
  template <typename T>
    using CV_t = typename std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T>::type, Eigen::Dynamic, 1> >;
  
  ///> Variable type (Real-valued vector)
  template <typename T>
    using Var_t = typename Eigen::Matrix<T, Eigen::Dynamic,1>;
  
}



namespace stan {
  namespace math {

  ///> Complex variable type ('templated typename')
  template <typename T>
    using C_t = typename std::vector<typename boost::math::tools::promote_arg<T>::type>;

  ///> Complex vector type
  template <typename T>
    using CV_t = typename std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T>::type, Eigen::Dynamic, 1> >;
  
  ///> Variable type (Real-valued vector)
  template <typename T>
    using Var_t = typename Eigen::Matrix<T, Eigen::Dynamic,1>;
  
}
}

#endif
