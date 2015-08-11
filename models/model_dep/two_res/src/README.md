src
===

The file 'model.hpp' contains functions describing PWA. For example, for
data generation:

 - f_model is the likelihood function, f(y,theta)=|theta*A_cv(y)|^2;
 - A_cv is the complex vector-valued function returning the partial wave amplitudes;
 - theta is the parameter vector describing complex numbers in front of the amplitudes.

The variable y=(m2_ab, m2_bc, ...) denotes the vector of variables that 
parametrize our decay.

For data fitting, it is necessary to define the normalization function  
`Norm(y,theta) = \int f_model(y, theta) dy`  
By substituting the definition of 'f_model' in the integral, one can rewrite 'Norm' as
`Norm(theta,I) = theta'* I theta`,
where `I[i,j] = \int A_cv[i](y)* A_cv[j](y) dy`. 

To calculate these integrals 'I', we use the script bin/data_analysis__root_to_data_R.py. This is a python script, in which we use the functions
f_model, A_cv, etc; hence, we need to convert these functions from C++ to
Python. To do so, we define the necessary wrappers in 'py_wrapper.cpp'.
This latter file may be compiled to a python module using bin/py_wrapper_setup.py, or simply by calling './../../../wrap_python.py' from the two_toy_res
directory. 
