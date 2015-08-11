#ifndef SRC__STRUCTURES__BASE_3D_RES_HPP
#define SRC__STRUCTURES__BASE_3D_RES_HPP


template <typename T>
class Base_3d_res {
public:
  base_3d_res() {};
  ~base_3d_res() {};

  virtual T value(const T& x, const T& y);
  virtual T value_sym(const T& x, const T& y);
}
#endif
