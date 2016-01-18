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
# 
# An additional file STAN_phase_space_gen.stan is build,
# if it is defined. This file is supposed to generate
# (uniformely distributed or incorporating detector 
# efficiency) events in the phase space. (These events
# are then used to evaluate Monte Carlo Integration.)

###### FUNCTIONS
function cd_stan_pwa
{
  while [[ $PWD != '/' && ${PWD##*/} != 'stan_pwa' ]]; do cd ..; done
}

###### MAIN
# Define model directory and stan_pwa directory
MODEL_DIR=$(pwd)
cd_stan_pwa
STAN_PWA_DIR=$(pwd)
cd $MODEL_DIR

mkdir -p build
mkdir -p output

python $STAN_PWA_DIR/bin/make_stan.py stan/STAN_data_generator.stan build/STAN_data_generator
python $STAN_PWA_DIR/bin/make_stan.py stan/STAN_amplitude_fitting.stan build/STAN_amplitude_fitting 

if [ -f stan/STAN_phase_space_gen.stan ];
then
    python $STAN_PWA_DIR/bin/make_stan.py stan/STAN_phase_space_gen.stan build/STAN_phase_space_gen
fi



