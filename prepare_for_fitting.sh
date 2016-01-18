#!/bin/bash

# prepare_for_fitting.sh
#
# USAGE
#
# ./../../../prepare_for_fitting.sh
#
# Launches bin/data_analysis__root_to_dataR.py --
# the script that takes event from .root tree and 
# prepares them for fitting (evaluates amplitudes at 
# given data events) and calculates Monte Carlo 
# integrals necessary for normalization.
#
# CAVEAT: run from the model folder.

./bin/data_analysis__root_to_dataR.py output/generated_data.root output/STAN_amplitude_fitting.data.R

