#ifndef SRC__STRUCTURES__SQR_HPP
#define SRC__STRUCTURES__SQR_HPP

#include "base_3d_res.hpp"

template <typename T>
class Sqr : public Base_3d_res<T> {
public:
  Sqr() {};
  Sqr() {};

  T value(const T& x, const T& y) {
    return x*y;
  }

  T value_sym(const T& x, const T& y) {
    return x*y + y*x;
  }
}
#endif
