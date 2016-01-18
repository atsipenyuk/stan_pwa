#!/bin/bash

# relink_model.sh
#   Links
#
#    *    stan_pwa/src/model.hpp
#
#   to
#
#    *    <PWD>/src/model.hpp


###### FUNCTIONS
function cd_stan_pwa
{
  while [[ $PWD != '/' && ${PWD##*/} != 'stan_pwa' ]]; do cd ..; done
}

###### MAIN
# Define model directory and meson_deca directory
MODEL_DIR=$(pwd)
cd_stan_pwa
STAN_PWA_DIR=$(pwd)
cd $MODEL_DIR

rm -f $STAN_PWA_DIR/src/model.hpp
ln -s $MODEL_DIR/src/model.hpp $STAN_PWA_DIR/src/model.hpp
