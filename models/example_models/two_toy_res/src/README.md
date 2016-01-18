src
===

The file 'model.hpp' contains functions describing PWA. For example, for
data generation:

 - f_genfit is the likelihood function, f(y,theta)=|theta*amplitude_vector(y)|^2;
 - amplitude_vector is the complex vector-valued function returning the partial wave amplitudes;
 - theta is the parameter vector describing complex numbers in front of the amplitudes.

The variable y=(m2_ab, m2_bc, ...) denotes the vector of variables that 
parametrize our decay.

For data fitting, it is necessary to define the normalization function  
`norm(y,theta) = \int f_genfit(y, theta) dy`  
By substituting the definition of 'f_genfit' in the integral, one can rewrite 'norm' as
`norm(theta,I) = theta'* I theta`,
where `I[i,j] = \int amplitude_vector[i](y)* amplitude_vector[j](y) dy`. 

To calculate these integrals 'I', we use the script bin/data_analysis__root_to_data_R.py. This is a python script, in which we use the functions
f_genfit, amplitude_vector, etc; hence, we need to convert these functions from C++ to
Python. To do so, we define the necessary wrappers in 'py_wrapper.cpp'.
This latter file may be compiled to a python module using bin/py_wrapper_setup.py, or simply by calling './../../../wrap_python.py' from the two_toy_res
directory. 
