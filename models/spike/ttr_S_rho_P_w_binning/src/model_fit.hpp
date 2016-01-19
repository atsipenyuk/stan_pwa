#ifndef STAN_PWA__SRC__MODEL_FIT_HPP
#define STAN_PWA__SRC__MODEL_FIT_HPP

#include <stan_pwa/src/structures/particles.hpp>

const int NUM_BINS_Y1=100; // Number of bins for the S wave for the 1st var
const int NUM_BINS_Y2=NUM_BINS_Y1; // We use symmetric model
const int NUM_NON_S_RES=2; // Number of non-S-wave resonances

namespace mcomplex = stan_pwa::complex;
namespace mreal = stan_pwa::real;

namespace stan {
  namespace math {


    // Describe the binning boundaries
    double bin_y1_min = 0.0;
    double bin_y1_max = 3.0;
    double bin_y2_min = 0.0;
    double bin_y2_max = 3.0;


    /**
     * float bin_y1_length()
     *
     * Returns the length of the bin in the phase space.
     */
    double bin_y1_length() {
      return (bin_y1_max - bin_y1_min) / double(NUM_BINS_Y1);
    }


    /**
     * float bin_volume()
     *
     * Returns the length of the bin in the phase space.
     */
    double bin_y2_length() {
      return (bin_y2_max - bin_y2_min) / double(NUM_BINS_Y2);
    }


    /**
     * integer bin_y1(double)
     *
     * Binning in the first Dalitz plot coordinate. Given a square mass, return the
     * number of the corresponding bin
     */
    template <typename T0>
    typename boost::math::tools::promote_args<int>::type
    bin_y1(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {
      // These values are decay-specific and should be defined in a more
      // appropriate way, e.g., in resonances.hpp.
      typename boost::math::tools::promote_args<int>::type res;

      if (y(0) >= bin_y1_max) {
	std::cout << "Dalitz plot variable out of bounds!";
	return false;
      }
      if (y(0) < bin_y1_min) {
	std::cout << "Dalitz plot variable out of bounds!";
	return false;
      }

      // This function is defined for symmetric decays. With other words, it 
      // suffices to consider binning only in the first coordinate.
      // This is a uniform binning.
      // res is from 1 .. NUM_BINS.
      res = 1 + int(floor((y(0) - bin_y1_min) / (bin_y1_max - bin_y1_min) * double(NUM_BINS_Y1)));
      return res;
    }


    /**
     * integer bin_y2(double)
     *
     * Binning in the first Dalitz plot coordinate. Given a square mass, return the
     * number of the corresponding bin
     */
    template <typename T0>
    typename boost::math::tools::promote_args<int>::type
    bin_y2(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

      typename boost::math::tools::promote_args<int>::type res;

      if (y(1) >= bin_y2_max) {
	std::cout << "Dalitz plot variable out of bounds!";
	return false;
      }
      if (y(1) < bin_y2_min) {
	std::cout << "Dalitz plot variable out of bounds!";
	return false;
      }

      res = 1 + int(floor((y(1) - bin_y2_min) / (bin_y2_max - bin_y2_min) * double(NUM_BINS_Y2)));
      return res;
    }


    /**
     * std::vector<integer> bin(std::vector<double>
     *
     * Given a coordinate in the phase space, return the corresponding bin numbers.
     */
    template <typename T0>
    std::vector<typename boost::math::tools::promote_args<int>::type>
    bin(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

      std::vector<typename boost::math::tools::promote_args<int>::type> res(2);
      res[0] = bin_y1(y);
      res[1] = bin_y2(y);

      return res;
    }


    /**
     * complex_vector A_FFZ_y1(vector). Decay-dependent.
     *
     * Returns the PWA S-wave amplitude of an event.
     */
    template <typename T0>
    typename boost::math::tools::promote_arg<T0>::type
    A_FFZ_y1(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

      typedef typename boost::math::tools::promote_arg<T0>::type T_res;

      //particle P = particles::DO;
      //particle c = particles::pi;

      //T_res m_ab = sqrt(y(0,0));

      // Form factor P->Rc, S-wave resonances have spin 0
      T_res F_P = 1.0;

      // Form factor R->ab
      T_res F_R = 1.0;

      // Zemach tensor (S-wave): 1 -> 0 + 
      T_res Z = 1.0;

      return F_P * F_R * Z;
    }


    /**
     * complex_vector A_FFZ_y2(vector). Decay-dependent.
     *
     * Returns the PWA S-wave amplitude of an event.
     */
    template <typename T0>
    typename boost::math::tools::promote_arg<T0>::type
    A_FFZ_y2(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

      typedef typename boost::math::tools::promote_arg<T0>::type T_res;

      // Form factor P->Rc, S-wave resonances have spin 0
      T_res F_P = 1.0;

      // Form factor R->ab
      T_res F_R = 1.0;

      // Zemach tensor (S-wave): 1 -> 0 + 
      T_res Z = 1.0;

      return F_P * F_R * Z;
    }


    /**
     * complex_vector amplitude_vector_non_S(vector)
     *
     * Returns the PWA waves of amplitude_vector that are fitted in a model-dependent way.
     * Relies on a correctly defined amplitude_vector function (S-wave resonances MUST come first).
     */
    template <typename T0>
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0>::type, Eigen::Dynamic, 1> >
    amplitude_vector_non_S(const Eigen::Matrix<T0, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0>::type T2;

        // Somewhat convoluted initialization of the return
        std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T2,Eigen::Dynamic,1> (NUM_NON_S_RES)));

        for (int i = 0; i < NUM_NON_S_RES; i++) {

            std::vector<T0> tmp;
            tmp = stan::math::A_c((NUM_RES - NUM_NON_S_RES) + i + 1, y);
            res[0](i) = tmp[0];
            res[1](i) = tmp[1];
        }
        return res;
    }


    /**
     * double f_fit( 
     *              scalar bin_FFZ_y1, scalar theta_bin_y1[2],
     *              scalar bin_FFZ_y2, scalar theta_bin_y2[2],
     *              vector A_y[2], vector theta[2],
     *              vector A_y_background_abs2, 
     *              vector theta_background_abs2)
     *
     *
     * Returns the model function for the given event 'i'.
     *
     * Parameters:
     *   bin_FFZ - a real number; bin_FFZ[i] contains the non-dynamical 
     *     description of the resonance at the event 'i'. With other
     *      words,
     *        bin_FFZ = Blatt-Weisskopf(Parent particle, at 'i') x
     *                  Blatt-Weisskopf(Resonance, at 'i') x
     *                  Zemach-Tensor(decay, at 'i').
     *   theta_bin - complex parameter, fitting the amplitude of the bin.
     *
     *   Since for two Dalitz variables we have two variables, we also bin
     *   in two axes. If the decay is symmetric, for example, D->3pi, then
     *   theta_bin_y1 and theta_bin_y2 represent bin contents of the same
     *   binning vector.
     *
     *   A_y - complex vector, contans model-dependent resonance amplitudes
     *     (P,D waves) evaluated at the event 'i'.
     *   theta - complex vector, contains fitting parameters for
     *     model-dependent parameter amplitudes (P,D waves).
     *
     *   A_y_background_abs2 - real vector, contains model-dependent squared
     *     absolute values of the background resonance amplitudes.
     *   theta_background_abs2 - real vector, fitting the background
     *     intensities.
     */
    template <typename T2, typename T3, typename T4,
	      typename T5, typename T6, typename T7>
    typename boost::math::tools::promote_args<T2,T3,T4,T5,T6,T7>::type
    f_fit(
      const T2& bin_TTZ_y1, const std::vector<T3>& theta_bin_y1,
      const T2& bin_TTZ_y2, const std::vector<T3>& theta_bin_y2,
      const std::vector<Eigen::Matrix<T4, Eigen::Dynamic, 1> >& A_r,
      const std::vector<Eigen::Matrix<T5, Eigen::Dynamic, 1> >& theta,
      const Eigen::Matrix<T6, Eigen::Dynamic, 1>& background_vector_,
      const Eigen::Matrix<T7, Eigen::Dynamic, 1>& theta_background_abs2_){

      return
        // Sum coherent amplitudes
        mcomplex::scalar::abs2(
          mcomplex::scalar::add(mcomplex::scalar::mult(bin_TTZ_y1, theta_bin_y1),
          mcomplex::scalar::add(mcomplex::scalar::mult(bin_TTZ_y2, theta_bin_y2),
				mcomplex::vector::sum(mcomplex::vector::mult(A_r, theta)))))
        // Incoherently add background amplitudes
        + mreal::vector::mult(background_vector_,theta_background_abs2_);
    }


    /**
     *
     * double norm_bin(complex_vector theta_bin_y1,
     *                complex_vector theta_bin_y2,
     *                complex_vector theta_res,
     *                real_vector I_y1y1,
     *                real_vector I_y2y2,
     *                complex_matrix I_y1y2,
     *                complex_matrix I_y1res,
     *                complex_matrix I_y2res,
     *                complex_matrix I_resres,
     *                real_vector theta_background_abs2_,
     *                real_vector I_background_abs2_)
     *
     * Takes complex vector theta and complex matrix I, real vectors
     * theta_background_abs2_, I_background_abs2_,
     * returns 
     *   conj(theta)' * I * theta + theta_background_abs2_ * 
     *                              I_background_abs2_,
     * where
     *          ____________________________________________________________
     *          |                      |                 |                 |
     *          |   diag (I_y1y1)      |     I_y1y2      |    I_y1res      |
     *          |______________________|_________________|_________________|
     *          |                      |                 |                 |
     *      I = | conj(transp(I_y1y2)) |  diag (I_y2y2)  |    I_y2res      |
     *          |______________________|_________________|_________________|
     *          |                      |                 |                 |
     *          |        (I_y1res)'*   |    (I_y2res)'*  |    I_resres     |
     *          |______________________|_________________|_________________|
     * 
     * and
     * 
     *              theta = [ theta_bin_y1 | theta_bin_y2 | theta_res ].
     */
    template <typename T1, typename T2, typename T3,
	      typename T7, typename T8>
    typename boost::math::tools::promote_args<T1,T2,T3,T7,T8>::type
    norm_bin(
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta_bin_y1,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta_bin_y2,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta_res,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& I_y1y1,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& I_y2y2,
      const std::vector<Eigen::Matrix<T3, Eigen::Dynamic, Eigen::Dynamic> >& I_y1y2,
      const std::vector<Eigen::Matrix<T3, Eigen::Dynamic, Eigen::Dynamic> >& I_y1res,
      const std::vector<Eigen::Matrix<T3, Eigen::Dynamic, Eigen::Dynamic> >& I_y2res,
      const std::vector<Eigen::Matrix<T3, Eigen::Dynamic, Eigen::Dynamic> >& I_resres,
      const Eigen::Matrix<T7, Eigen::Dynamic, 1>& theta_background_abs2_,
      const Eigen::Matrix<T8, Eigen::Dynamic, 1>& I_background_abs2_) {

      typedef typename boost::math::tools::promote_args<T1,T2,T3,T7,T8>::type T_res;

      T_res res = 0;

      T_res temp;
      std::vector<T_res> current_theta_bin(2);

      // I * theta holder, real and imaginary part
      T_res tmp[2];

      // The sum has 9 matrix entries + background: 
      //   1) theta_y1  I_y1y1    theta_y1
      //   2) theta_y2  I_y2y2    theta_y2
      //   3) theta_res I_resres  theta_res

      //   4) theta_y1  I_y1y2    theta_y2 
      //   5) theta_y1  I_y1res   theta_res
      //   6) theta_y2  I_y2res   theta_res

      //   7) theta_y2 (I_y1y2)'* theta_y1
      //   8) theta_res(I_y1res)'*theta_y1
      //   9) theta_res(I_y2res)'*theta_y2

      //   and then background. Note that cases
      //   4 and 7 are similar and can be treated simultaneously,
      //   as well as cases 5 and 8, 6 and 9.

      // Add theta_y1 * diag(I_y1y1) * theta_y1 (add real part;
      // imag. part should be zero
      for (int i = 0; i < NUM_BINS_Y1; i++) {
        current_theta_bin[0] = theta_bin_y1[0](i);
	current_theta_bin[1] = theta_bin_y1[1](i);
	temp =  mcomplex::scalar::mult(I_y1y1(i),
	        mcomplex::scalar::mult(current_theta_bin, 
		                       current_theta_bin))[0];
      res = res + temp;
      }

      // Add 2) analogously
      for (int i = 0; i < NUM_BINS_Y2; i++) {
        current_theta_bin[0] = theta_bin_y2[0](i);
	current_theta_bin[1] = theta_bin_y2[1](i);
	temp =  mcomplex::scalar::mult(I_y2y2(i),
	        mcomplex::scalar::mult(current_theta_bin, 
		                       current_theta_bin))[0];
      res = res + temp;
      }
 
      // Add 3) theta_res * I_resres * theta_res
      for (int i = 0; i < NUM_NON_S_RES; i++) {
      for (int j = 0; j < NUM_NON_S_RES; j++) {
        // Complex multiplication
        tmp[0] = I_resres[0](i,j) * theta_res[0](j) - 
	         I_resres[1](i,j) * theta_res[1](j);
        tmp[1] = I_resres[0](i,j) * theta_res[1](j) + 
	         I_resres[1](i,j) * theta_res[0](j);
          // Keep only the real palt of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
        res = res + theta_res[0](i) * tmp[0] + theta_res[1](i) * tmp[1];
        }
      }


     // Add 4) theta_y1 * I_y1y2 * theta_y2 and 7 (its complex conjugate).
      for (int i = 0; i < NUM_BINS_Y1; i++) {
      for (int j = 0; j < NUM_BINS_Y2; j++) {
        // Complex multiplication
        tmp[0] = I_y1y2[0](i,j) * theta_bin_y2[0](j) - 
	         I_y1y2[1](i,j) * theta_bin_y2[1](j);
        tmp[1] = I_y1y2[0](i,j) * theta_bin_y2[1](j) + 
	         I_y1y2[1](i,j) * theta_bin_y2[0](j);
	// Keep only the real palt of the product; imaginary part
	// should be 0 (+- float calculation errors).
	// Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1].
	// The factors 2.0 account for the complex conjugate part.
        res = res + 2.0 * theta_bin_y1[0](i) * tmp[0] 
	          + 2.0 * theta_bin_y1[1](i) * tmp[1];
        }
      }


     // Add 5) theta_y1 * I_y1res * theta_res and 8 (its complex conjugate).
      for (int i = 0; i < NUM_BINS_Y1; i++) {
      for (int j = 0; j < NUM_NON_S_RES; j++) {
        // Complex multiplication
        tmp[0] = I_y1res[0](i,j) * theta_res[0](j) - 
	         I_y1res[1](i,j) * theta_res[1](j);
        tmp[1] = I_y1res[0](i,j) * theta_res[1](j) + 
	         I_y1res[1](i,j) * theta_res[0](j);
	// Keep only the real palt of the product; imaginary part
	// should be 0 (+- float calculation errors).
	// Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1].
	// The factors 2.0 account for the complex conjugate part.
        res = res + 2.0 * theta_bin_y1[0](i) * tmp[0] 
	          + 2.0 * theta_bin_y1[1](i) * tmp[1];
        }
      }


     // Add 6) theta_y2 * I_y2res * theta_res and 9 (its complex conjugate).
      for (int i = 0; i < NUM_BINS_Y2; i++) {
      for (int j = 0; j < NUM_NON_S_RES; j++) {
        // Complex multiplication
        tmp[0] = I_y2res[0](i,j) * theta_res[0](j) - 
	         I_y2res[1](i,j) * theta_res[1](j);
        tmp[1] = I_y2res[0](i,j) * theta_res[1](j) + 
	         I_y2res[1](i,j) * theta_res[0](j);
	// Keep only the real palt of the product; imaginary part
	// should be 0 (+- float calculation errors).
	// Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1].
	// The factors 2.0 account for the complex conjugate part.
        res = res + 2.0 * theta_bin_y2[0](i) * tmp[0] 
	          + 2.0 * theta_bin_y2[1](i) * tmp[1];
        }
      }


      // Add background integral component
      res = res + mreal::vector::mult(theta_background_abs2_,I_background_abs2_);

      return res;
    }


    /**
     * int num_bins_y1()
     *
     * Returns the number of bins in the first coordinate.
     */
    inline int num_bins_y1() {
      return NUM_BINS_Y1;
    }


    /**
     * int num_bins_y2()
     *
     * Returns the number of bins in the second coordinate.
     */
    inline int num_bins_y2() {
      return NUM_BINS_Y2;
    }


    /**
     * int num_non_S_res()
     *
     * Returns the number of model-dependent resonances (resonances whose
     * amplitude will be fitted not via bins).
     */
    inline int num_non_S_res() {
      return NUM_NON_S_RES;
    }

  }
}


#endif // MESON_DECA__LIB__C_LIB__MODEL_FIT_HPP
