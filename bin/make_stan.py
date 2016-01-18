#!/usr/bin/env python
# make_stan.py
#
# NAME
#    make_stan.py - build an executable from STAN model. 
#
# SYNOPSIS
#    ./make_stan.py MODEL_STAN MODEL_EXEC
#
# DESCRIPTION
#    This is a small utility script. It changes the directory upwards until
#    the directory 'stan_pwa' is encountered, and goes one directory up (landing
#    in the cmdstan folder). Then, it calls the cmdstan 'make' to build the
#    executable MODEL_STAN.


import argparse
import os
import subprocess

import utils

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to build *.stan to an executable.')

parser.add_argument('f_in',
                     nargs='?')

parser.add_argument('f_out',
                    nargs='?')

args = parser.parse_args()

model_path = os.getcwd()
cmdstan_path = os.path.split(utils.get_path('stan_pwa'))[0]

# Relative path from cmdstan to model
i = model_path.find(cmdstan_path) + len(cmdstan_path) + 1
rel_model_path = model_path[i:]

# Relative path from cmdstan to f_in and f_out
rel_f_in_path_stan = os.path.join(rel_model_path, args.f_in)
rel_f_out_path = os.path.join(rel_model_path, args.f_out)

# Check that the passed f_in file has '.stan' extension
if rel_f_in_path_stan[-5:] != '.stan':
    raise InputError('First input argument must be a \'.stan\' file. Got ' + \
                     args.f_in + 'instead.\n')

# Get rid of the .stan ending in rel_f_in_path_stan
rel_f_in_path = rel_f_in_path_stan[:-5]

os.chdir(cmdstan_path)
subprocess.Popen("make " + rel_f_in_path, shell=True).wait()

subprocess.Popen("mv" + " " + rel_f_in_path + \
                 " " + rel_f_out_path, shell=True).wait()

print('Done. \'Generated ' + os.path.join(cmdstan_path, rel_f_out_path) + '\'.\n')
