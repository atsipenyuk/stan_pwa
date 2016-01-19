#!/bin/bash

# merge_output_chains.sh
#   Merges files
#
#    *    output/output?.csv
#
#   to a single file
#    *    output/output.csv
#   and translates the latter one to output/output.root

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

# Merge
grep lp__ output/output1.csv > output/output.csv
sed '/^[#l]/d'  output/output?.csv >> output/output.csv

# Translate to root
$STAN_PWA_DIR/bin/csv_to_root.py $MODEL_DIR/output/output.csv $MODEL_DIR/output/output.root






