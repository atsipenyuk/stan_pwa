#ifndef SRC__STRUCTURES__MODEL_INST_HPP
#define SRC__STRUCTURES__MODEL_INST_HPP


#include <stan_pwa/src/flat_structures/particles.hpp>
#include <stan_pwa/src/model_def.hpp>
#include <stan_pwa/src/flat_structures/bw.hpp>
#include <stan_pwa/src/structures/three_body/bw.hpp>

namespace stan_pwa {

  ExternalParticles4 Ext(particles::d,
			 particles::pi,
			 particles::pi,
			 particles::pi);

  Particle_w R(particles::f0_1370, 0.350);

  stan_pwa::resonances::breit_wigner breit_wigner1(particles::d,
					  particles::pi,
					  particles::pi,
					  particles::pi,
					  particles::f0_1370, 0.2);

  // External cpp files here cause problems in the STAN linker
  // They may prabably be resolved adjusting CFLAGS in Stan makefile
  //template<typename T>
  Model<double> My_model();
  
  //template<typename T>
  //My_model.amplitudes.push_back(*breit_wigner<double,Ext,R>)

  template <typename T> 
  C_t<T> (*my_ptr) (const T&, const T&) = &breit_wigner<T,Ext,R>;
  
  template <typename T>
  //  std::vector<typename boost::math::tools::promote_arg<T>::type>
  //C_t<T>
  C_t<typename boost::math::tools::promote_arg<T>::type> 
  amplitude(unsigned int i, const Var_t<T>& y) {
    return breit_wigner<T,Ext,R>(y(0), y(1));
  };

 
  template <typename T>
  std::vector< Eigen::Matrix<typename boost::math::tools::promote_arg<T>::type, Eigen::Dynamic, 1> >
  //CV_t<typename boost::math::tools::promote_arg<T>::type> 
  amplitude_vector(const Var_t<T>& y) {
    CV_t<typename boost::math::tools::promote_arg<T>::type> 
      res(2, (Eigen::Matrix<T,Eigen::Dynamic,1> (num_res_)));
        for (int i = 0; i < num_res_; i++) {
            std::vector<T> tmp;
            tmp =  my_ptr<T>(y(0), y(1));
            //tmp =  breit_wigner1.value(y(0),y(1));
            res[0](i) = tmp[0];
            res[1](i) = tmp[1];
        }
    return res;
  };


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type ///> return scalar   
  f_genfit(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& x,
		  const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& y)
  {
    return mc::scalar::abs2(mc::vector::sum(mc::vector::mult(x,y)));
  }


  template <typename T0, typename T1>
  typename boost::math::tools::promote_args<T0,T1>::type
  norm(const CV_t<T0>& theta, const CV_t<T0>& I) {
      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < num_res_; i++) {
	for (int j = 0; j < num_res_; j++) {
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
  };

  
  //Sqr my_res_1(P, a, b, c);//, particles::rho_770, 0.1491);
  /*  Sqr my_res_2(P, a, b, c);//, particles::f0_1370, 0.350);
  
  std::vector<Base_3d_res> my_res_list = {my_res_1, my_res_2};
  bool sym_flag = 1;
  
  Model My_model(my_res_list, 2);*/

}
#endif
