#!/bin/bash

# relink_model.sh
#   Links
#
#    *    stan_pwa/src/model_def.hpp
#    *    stan_pwa/src/model.cpp
#    *    stan_pwa/src/model_inst.hpp
#
#   to
#
#    *    <PWD>/src/model_def.hpp
#    *    <PWD>/src/model.cpp
#    *    <PWD>/src/model_inst.hpp

###### FUNCTIONS
function cd_stan_pwa
{
  while [[ $PWD != '/' && ${PWD##*/} != 'stan_pwa' ]]; do cd ..; done
}

###### MAIN
# Define model directory and meson_deca directory
MODEL_DIR=$(pwd)
cd_stan_pwa
MDECA_DIR=$(pwd)
cd $MODEL_DIR

rm $MDECA_DIR/src/model_def.hpp
rm $MDECA_DIR/src/model.cpp
rm $MDECA_DIR/src/model_inst.hpp
rm $MDECA_DIR/src/model_wrapper.hpp


ln -s $MODEL_DIR/src/model_def.hpp $MDECA_DIR/src/model_def.hpp
#ln -s $MODEL_DIR/src/model.cpp $MDECA_DIR/src/model.cpp
ln -s $MODEL_DIR/src/model_inst.hpp $MDECA_DIR/src/model_inst.hpp
ln -s $MODEL_DIR/src/model_wrapper.hpp $MDECA_DIR/src/model_wrapper.hpp
