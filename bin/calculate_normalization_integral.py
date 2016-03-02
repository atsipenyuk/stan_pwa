#!/usr/bin/env python
# calculate_normalization_integral.py

# NAME
#    calculate_normalization_integral.py - make a file containing I[i,j]
#
# SYNOPSIS
#    ./calculate_normalization_integral.py [OPTIONS...] l_bd h_bd l_bd h_bd
#
# EXAMPLE for 3-body decay usage, model with background:
#    ./calculate_normalization_integral.py -b 0 3 0 3
#
# DESCRIPTION
#    calculate_normalization_integral.py calculates the integrals
#    necessary to perform the correct normalization of the
#    log-likelihood function during STAN_amplitude_fitting
#    sampling.
#
#    This means the following: the function model.A_cv (python-wrapped
#    analogon of the STAN function A_cv) contains R =
#    model.num_resonances() amplitudes. This script calculates
#
#        I[i,j] = \int conj(A_r(i)) A_r(j) \dy,
#
#    where y is the model variable - for examle, the invariant mass of
#    the decay. The integral is calculated using Monte-Carlo
#    integration and is stored in 'normalization_integral.py'. The
#    integration bounds may be specified in [OPTIONS] - just pass the
#    arguments in the form y1_min, ... , yR_min, y1max, ... , yR_max.
#
# CAVEAT. This script MUST be called from the folder containing the
#    module model.so corresponding to the described model.

import argparse
import numpy as np
import os
import sys


import utils # integral for Monte Carlo integration
             # convert_to_vector - translates A_cv results to usable form
             # vector_to_string  - to save results in a .py file

# Import PWA data from current model
sys.path.insert(1, os.path.join(os.getcwdu(),'build'))
import model  # Function to be integrated

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to calculate the integral leading to the Normalization intergral.')

# Parse arguments - two boundaries for each resonance!
# y1min, y1max, ... yRmin, yRmax
N = model.num_variables()
R = model.num_resonances()

parser.add_argument('bounds',
                    nargs=2*N,
                    type=float,
                    help="integration bounds.")
parser.add_argument('--background','-b',
                    action='store_true',
                    help="Flag: check for background amplitudes.")

parser.add_argument('--num_int_samples','-n',
                    default=1000000,
                    type=int,
                    help="Number of Monte Carlo integration samples.")

args = parser.parse_args()

I = np.zeros([R,R], dtype=complex)

# Complex vector containing all PWA amplitudes
def func(y):
    return utils.convert_to_complex_vector(model.amplitude_vector(N, y))


# Repack the integration bounds
bounds = [[args.bounds[2*n], args.bounds[2*n+1]] for n in range(N)]

# Calculate the integral matrix
I = utils.integral_of_tensor_product(func, bounds, N=args.num_int_samples)

# If necessary, calculate background amplitude normalization
if args.background == 1:
    B = model.num_background()
    print "Calculating background integral...\n"
    def func_backgr(y):
        return utils.convert_to_vector(model.background_vector(B, y))
    I_background = np.zeros(B, dtype=float)
    I_background = utils.integral(func_backgr, bounds, N=args.num_int_samples)


f_py = open('normalization_integral.py', 'w')
f_py.write('I_ = ' + utils.array_to_string(I[0]) + '\n')
f_py.write('I_err_ = ' + utils.array_to_string(I[1]) + '\n') # error estimate
if args.background == 1:
    f_py.write('I_background_ = ' + utils.vector_to_string(I_background[0]) + '\n')
    f_py.write('I_background_err_ = ' + utils.vector_to_string(I_background[1]) + '\n')
f_py.close()
