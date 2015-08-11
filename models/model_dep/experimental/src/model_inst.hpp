#ifndef SRC__STRUCTURES__MODEL_INST_HPP
#define SRC__STRUCTURES__MODEL_INST_HPP


#include "resonances/square_res.hpp"
#include "model_def.hpp"

template <typename T>
Sqr<T> my_res_1;
template <typename T>
Sqr<T> my_res_2;

std::vector<Base_3d_res<T> > my_res_list = {my_res_1, my_res_2};

Model My_model(my_res_list);


#endif
