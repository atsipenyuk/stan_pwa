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
     * Boost::Python can convert std::vector's to Python; more explicitely
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
  }
}


BOOST_PYTHON_MODULE(model)
{
    using namespace boost::python;
    class_<std::vector<double> >("StdVrDouble")
        .def(vector_indexing_suite<std::vector<double> >() );

    def("amplitude_vector", stan::math::_amplitude_vector_py_wrapper, args("x","y"));
    def("num_resonances", stan::math::num_resonances);
    def("num_variables", stan::math::num_variables);
}


