#!/usr/bin/env python
# polar_data_R_to_re_im_data_R.py
#
# DESCRIPTION
#    polar_data_R_to_re_im_data_R.py - read a polar.data.R file. Treat all 
#    variables  with '*_polar' name as complex vectors written as 
#    (magnitude1, phase1, magnitude2, phase2, ...); compute the corresponding
#    vector (real1, imag1, real2, imag2, ...) and save it to the 
#    .data.R file. Phase must be written in DEGREES.
#
# SYNOPSIS
#     polar_data_R_to_re_im_data_R.py POLAR_DATA_R
#

import argparse
import numpy as np
import os
from itertools import chain

import utils

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert STAN output.csv to a .ROOT file.')

parser.add_argument('--f_in',
                    default = 'stan/STAN_data_generator.polar.data.R',
                    type=argparse.FileType('r'))

args = parser.parse_args()


#
print('polar_data_R_to_re_im_data_R.py: Reading {0}...'.format(args.f_in.name))

# Open CSV file
#print('Importing data...')
f_in = args.f_in

# Open ROOT file
f_out = open(args.f_in.name[:-len('.polar.data.R')] + '.data.R', 'w')


for line in f_in:
    # Parse the variables that end with '_polar'
    if '_polar' in line:
        lhs, rhs = line.split('<-')
        lhs = lhs.split('_polar')[0]

        if 'structure' in rhs:
            rhs = rhs.split('structure(')[-1] # throw away 'structure'
            rhs = rhs[:-1] # Throw away the ')' corresponding to 'structure'
            arr, dim = rhs.split(', .Dim =') # split array and dimension
            dim_string = '.Dim = ' + dim

            # throw away the 'c(' and ')'
            arr = arr.split('c(')[-1]
            arr = arr.rsplit(')')[0]
            dim = dim.split('c(')[-1]
            dim = dim.rsplit(')')[0]

            # get the numbers
            arr = np.asarray(arr.split(','), dtype=float)
            dim = np.asarray(dim.split(','), dtype=int)

            # arr = arr.reshape(dim, order='F')
            # Convert to polar vector
            mag = arr[::2]
            phase = arr[1::2]
            real = mag * np.cos(phase / 180 * np.pi)
            imag = mag * np.sin(phase / 180 * np.pi)
            # Zip mag and phase together
            res = list(chain(*zip(real,imag)))
            print(res)

            # Write results to file
            f_out.write(lhs + ' <- ' + 'structure(')
            f_out.write('c(')
            for val in res[:-1]:
                f_out.write(str(val) + ', ')
            f_out.write(str(res[-1]) + ')')
            
            # Write the dimension to file
            f_out.write(', ' + dim_string)
            # Close structure
            f_out.write(')\n')

        else:
            raise ValueError;

    else:
        f_out.write(line)

print("polar_data_R_to_re_im_data_R.py: Done. Data saved in {0}.".format(f_out.name))

f_out.close()
f_in.close()


