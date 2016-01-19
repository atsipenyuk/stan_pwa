#include<iostream>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <stan_pwa/src/model.hpp>


// Python wrapper for functions specified in model.hpp

namespace stan {
  namespace math {

    /**
     * vector _amplitude_vector_py_wrapper(vector)
     *
     * Argument wrapper for the function amplitude_vector to call amplitude_vector from python.
     *
     * We need to convert python list to Eigen::Matrix, pass it as an 
     * argument to amplitude_vector, and then convert the result back to some 
     * pythonian type.
     *
     *
     * Boost::Python can convert std::vector to Python; more explicitely
     * - to some python class; in our case, this class is specified further 
     * below and it is called StdVrDouble.
     */
     inline
     std::vector<double>
     _amplitude_vector_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        // Call amplitude_vector
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > res;
        res = amplitude_vector(y);
        if (res.size() != 2)
	  std::cout << "Something's wrong here - dimension mismatch.";

        // Return a flatten version of res - an std::vector
        std::vector<double> std_res(2*NUM_RES);
        for (int i=0; i < NUM_RES ; i++) {
	  std_res[2*i] = res[0](i,0);
          std_res[2*i + 1] = res[1](i,0);
        }
        return std_res;
     }

    /**
     * vector _amplitude_vector_non_S_py_wrapper(vector)
     *
     * Argument wrapper for the function amplitude_vector to call amplitude_vector from python.
     */
     inline
     std::vector<double>
     _amplitude_vector_non_S_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        // Call amplitude_vector
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > res;
        res = amplitude_vector_non_S(y);
        if (res.size() != 2)
	  std::cout << "Something's wrong here - dimension mismatch.";

        // Return a flatten version of res - an std::vector
        std::vector<double> std_res(2*NUM_NON_S_RES);
        for (int i=0; i < NUM_NON_S_RES ; i++) {
	  std_res[2*i] = res[0](i,0);
          std_res[2*i + 1] = res[1](i,0);
        }
        return std_res;
     }



    /**
     * vector _A_FFZ_y1_py_wrapper(vector)
     *
     * Argument wrapper for the function A_FFZ_y1.
     */
     inline
     double
     _A_FFZ_y1_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        return  A_FFZ_y1(y);
     }


    /**
     * vector _A_FFZ_y2_py_wrapper(vector)
     *
     * Argument wrapper for the function A_FFZ_y1.
     */
     inline
     double
     _A_FFZ_y2_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        return A_FFZ_y2(y); 
     }


    /**
     * vector _A_v_backgr_py_wrapper(vector)
     *
     * Argument wrapper for the function background_vector to call it from 
     * python as background_vector.
     *
     * Proceed as in the function above with mutatis mutandis siplifications.
     */
     inline
     std::vector<double>
     _A_v_backgr_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        // Call A_v and convert Eigen to std::vector
	std::vector<double> res(num_background());
	Eigen::Matrix<double, Eigen::Dynamic, 1> res_eigen = background_vector(y);
	for (int i=0; i<num_background(); i++) {
	  res[i] = res_eigen(i);
	}

        return res;
     }


    /**
     * vector _bin_py_wrapper(vector)
     *
     * Argument wrapper for the function bin().
     */
     inline
     std::vector<double>
     _bin_py_wrapper(int y_len, boost::python::list mapping) {

        // Convert python list to Eigen::Matrix
        Eigen::Matrix<double, Eigen::Dynamic, 1> y(y_len);
        for (int i=0; i<y_len; i++) {
            y[i] = boost::python::extract<double>(mapping[i]);
        }

        return stan::math::bin(y);
     }

  }
}


BOOST_PYTHON_MODULE(model)
{
    using namespace boost::python;
    class_<std::vector<double> >("StdVrDouble")
        .def(vector_indexing_suite<std::vector<double> >() );

    def("A_FFZ_y1", stan::math::_A_FFZ_y1_py_wrapper, args("x","y"));
    def("A_FFZ_y2", stan::math::_A_FFZ_y2_py_wrapper, args("x","y"));
    def("amplitude_vector", stan::math::_amplitude_vector_py_wrapper, args("x","y"));
    def("amplitude_vector_non_S", stan::math::_amplitude_vector_non_S_py_wrapper, args("x","y"));
    def("background_vector", stan::math::_A_v_backgr_py_wrapper, args("x","y"));
    def("bin", stan::math::_bin_py_wrapper);
    def("bin_y1_length", stan::math::bin_y1_length);
    def("bin_y2_length", stan::math::bin_y2_length);
    def("num_background", stan::math::num_background);
    def("num_bins_y1", stan::math::num_bins_y1);
    def("num_bins_y2", stan::math::num_bins_y2);
    def("num_non_S_res", stan::math::num_non_S_res);
    def("num_resonances", stan::math::num_resonances);
    def("num_variables", stan::math::num_variables);
}







