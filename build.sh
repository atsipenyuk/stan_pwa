#!/bin/bash

# build.sh
#   Builds
#
#    *    build/STAN_data_generator
#    *    build/STAN_amplitude_fitting
#
#   from the corresponding *.stan files stored in
#    *    stan/STAN_data_generator.stan
#    *    stan/STAN_amplitude_fitting.stan

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

mkdir -p build
mkdir -p output

python $MDECA_DIR/scripts/make_stan.py stan/STAN_data_generator.stan build/STAN_data_generator
#python $MDECA_DIR/scripts/make_stan.py stan/STAN_amplitude_fitting.stan build/STAN_amplitude_fitting 

if [ -f stan/STAN_phase_space_gen.stan ];
then
    python $MDECA_DIR/scripts/make_stan.py stan/STAN_phase_space_gen.stan build/STAN_phase_space_gen
fi



