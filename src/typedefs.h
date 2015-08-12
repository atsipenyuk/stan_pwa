#ifndef PWA_STAN__SRC__TYPEDEFS_H
#define PWA_STAN__SRC__TYPEDEFS_H

// These types are quite useful, so I define them in both namespaces

namespace stan_pwa {

  ///> Complex variable type ('templated typename')
  template <typename T>
    using C_t = std::vector<T>;

  ///> Complex vector type
  template <typename T>
    using CV_t = std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >;

  ///> Variable type (Real-valued vector)
  template <typename T>
    using Var_t = Eigen::Matrix<T, Eigen::Dynamic,1>;

  ///> Variable type (Real-valued vector)
  template <typename T>
    using promoted_t = typename boost::math::tools::promote_args<T>::type;
  
}



namespace stan {
  namespace math {

  ///> Complex variable type ('templated typename')
  template <typename T>
    using C_t = std::vector<T>;

  ///> Complex vector type
  template <typename T>
    using CV_t = std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >;

  ///> Variable type (Real-valued vector)
  template <typename T>
    using Var_t = typename Eigen::Matrix<T, Eigen::Dynamic,1>;


  ///> Variable type (Real-valued vector)
  template <typename T>
    using promoted_t = typename boost::math::tools::promote_args<T>::type;
  
}
}

#endif
