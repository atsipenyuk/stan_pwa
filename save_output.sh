#!/bin/bash

# NAME
#   save_output.sh - copies output.root and other files to specified location
#
# SYNOPSIS
#   save_output.sh FOLDER_NAME
#
# DESCRIPTION
#   Copies
#
#    *    output/output.root
#    *    output/generated_data.root
#    *    output/STAN_amplitude_fitting.data.R
#    *    stan/STAN_data_generator.data.R
#
#   to the folder FOLDER_NAME.
#
# CAVEAT
#  Uses relative paths. MUST be run from the model directory.

FOLDER_NAME=$1
mkdir -p $FOLDER_NAME
cp output/output.root $FOLDER_NAME
cp output/generated_data.root $FOLDER_NAME
cp output/STAN_amplitude_fitting.data.R $FOLDER_NAME
cp stan/STAN_data_generator.data.R $FOLDER_NAME
cp background_to_signal_ratio.py $FOLDER_NAME
cp normalization_integral.py $FOLDER_NAME
