#!/bin/bash

# generate.sh NUM_SAMPLES
#
# USAGE
#
# ./../../../generate.sh NUM_SAMPLES
#
# where NUM_SAMPLES is an integer, for example, 10000.
#
# Generate the data using build/STAN_data_generator executable, with
# NUM_SAMPLES events. The samples are converted to a '.root' file.
# The results are saved in output/generated_data.csv and 
# output/generated_data.root.
#
# CAVEAT: run from the model folder.

###### FUNCTIONS
function cd_stan_pwa
{
  while [[ $PWD != '/' && ${PWD##*/} != 'stan_pwa' ]]; do cd ..; done
}

###### MAIN
# Define locations of necessary files
MODEL_DIR=$PWD
cd_stan_pwa
STAN_PWA_DIR=$PWD
cd $MODEL_DIR

# If no arguments are supplied, generate 1000 samples
if [ $# -eq 0 ]
  then
    NUM_SAMPLES=1000
  else
    NUM_SAMPLES=$1
fi

# Generate data
./build/STAN_data_generator sample num_samples=$NUM_SAMPLES data file=stan/STAN_data_generator.data.R output file=output/generated_data.csv

# Create root file containing generated events
${STAN_PWA_DIR}/bin/csv_to_root.py output/generated_data.csv output/generated_data.root

# Run data analysis: evaluate Monte Carlo integrals, evaluate amplitude at the events
${MODEL_DIR}/bin/data_analysis__root_to_dataR.py

