#!/bin/bash

# wrap_python.sh
#   This script rebuilds src/py_wrapper.cpp to a python module and copies
#   that module into the 'build' folder


###### FUNCTIONS
function cd_stan_pwa
{
  while [[ $PWD != '/' && ${PWD##*/} != 'stan_pwa' ]]; do cd ..; done
}

###### MAIN
# Define model directory and meson_deca directory
MODEL_FOLDER=$(pwd)
cd_stan_pwa
MDECA_FOLDER=$(pwd)
cd "$MODEL_FOLDER"

rm -rf build/lib*
rm -rf build/temp*
rm -f build/model.so

python "${MDECA_FOLDER}/bin/py_wrapper_setup.py" "build"
cp build/lib*/model.so build
