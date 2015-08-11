#ifndef PWA_STAN__SRC__MODEL_DEF_HPP
#define PWA_STAN__SRC__MODEL_DEF_HPP

#include <vector>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/any.hpp>

#include <stan_pwa/src/structures.hpp>
#include <stan_pwa/src/typedefs.h>

namespace stan_pwa {

  /**
   * This class wraps together multiple resonances to single PWA model.
   */
  class Model {
  public:
    ///> Set number of variables: 2 for 3-body-decay, 5 for 4-body-decay.
    Model(unsigned int num_var, bool sym_flag, 
	  std::vector<resonances::resonance_base_3> amplitudes) : 
      num_var_(num_var), 
      num_res_(amplitudes.size()),
      sym_flag_(sym_flag),
      amplitudes_(amplitudes)
    {};
    ~Model() {};

    ///> returns the n-th resonance value (complex number)
    template <typename T>
    C_t<T> amplitude(unsigned int, const Var_t<T>&);

    ///> returns all resonance values
    template <typename T>
    CV_t<T> amplitude_vector(const Var_t<T>&);

    ///> Same as above, but symmetrized
    template <typename T>
    C_t<T> amplitude_sym(unsigned int, const Var_t<T>&);

    template <typename T>
    CV_t<T> amplitude_vector_sym(const Var_t<T>&);

    ///> Combines vectors and amplitudes to (un-normalized) likelihood fct
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar
    //f_genfit(const CV_t<T0>&, const CV_t<T1>&);
    f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >&,
	     const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >&);

    ///> Calculates the normalization integral for the fitting
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    norm(const CV_t<T0>&, const CV_t<T0>&);

    // get_num_res
    int get_num_res() {return num_res_;}

    // get_num_var
    int get_num_var() {return num_var_;};

    bool get_sym_flag() {return sym_flag_;}

  private:

    unsigned int num_res_; ///> Number of resonances
    unsigned int num_var_; ///> Number of variables

    ///>  whether the model must be symmetrized or not
    bool sym_flag_;

    ///> Vector containing PWA amplitude functions
    std::vector<resonances::resonance_base_3> amplitudes_;
  };


} // end of stan_pwa
#endif






