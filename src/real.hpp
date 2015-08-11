#ifndef STAN_PWA__SRC__REAL_HPP
#define STAN_PWA__SRC__REAL_HPP

#include <stan_pwa/src/real/vector.hpp>

/*
 *  Introduce real-valued number operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    Just as for complex numbers, we need to redefine some functions
 *    like a dot product in a general, STAN-friendly way.
 * 
 *    You might think: wait, this is ridiculous! STAN already provides
 *    nice functions for real-valued vectors and matrices. And this is
 *    true: for example, in CmdSTAN 2.6.2, look at 
 *    stan/src/stan/math/prim/mat/fun/multiply.hpp.
 *    And you shall see that they are defined for vectors and matrices
 *    of scalar type 'double'. Unfortunately, I can not make their
 *    code work for, say, dot multiplication of Eigen::Matrix<double, ...>
 *    with Eigen::Matrix<stan::agrad::var, ...>. The trouble I am having
 *    probably lies in the instantiation order or conversion order...
 *    Well, anyway, I am working around that without really understanding
 *    why the STAN templates do not work, so that is probably not the 
 *    smartest idea. My bad!
 *
 *  FUNCTIONS
 *    Are currently listed in particular files - vector.hpp.
 */

#endif
