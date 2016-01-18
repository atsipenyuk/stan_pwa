#ifndef STAN_PWA__SRC__COMPLEX_HPP
#define STAN_PWA__SRC__COMPLEX_HPP

#include <stan_pwa/src/complex/scalar.hpp>
#include <stan_pwa/src/complex/vector.hpp>
#include <stan_pwa/src/complex/matrix.hpp>

/*
 *  Introduce complex number operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    Complex numbers are described as std::vector<double>
 *    (STAN array of reals; also referred to as 'complex_scalar').
 * 
 *    Complex vectors are described as 
 *      std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >.
 *    (STAN array of vectors; also referred to as 'complex_vector').
 *
 *    Complex matrices are described as 
 *      std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >.
 *    (STAN array of matrices; also referred to as 'complex_matrix').
 *
 *    There are no new type definitions for complex-valued objects described
 *    above. None of the operators are overloaded. The reason is that STAN
 *    would not be able to differentiate between a complex vector and an array
 *    of vectors. Therefore, BEWARE the following: under such operators as +, 
 *    *, etc., all complex-valued objects behave as std::vectors. For example,
 *    if you want to multiply two complex numbers, do not use '*' operator, 
 *    but rather the function complex::scalar::mult. But using the operator 
 *    '+' is okay, since (in Cartesian coordinates) the sum of two complex 
 *    numbers is isomorphic to the sum of two vectors.
 *
 *    Since no new types/classes are really introduced, all distinction between
 *    complex objects is implemented via namespaces. 
 *
 *  FUNCTIONS
 *    Are currently listed in particular files - scalar.hpp, vector.hpp, matrix.hpp.
 */

#endif
