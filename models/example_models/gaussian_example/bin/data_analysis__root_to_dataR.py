#!/usr/bin/env python
#
# data_analysis__root_to_dataR.py
#
# DESCRIPTION
#    The events for our model are generated as a real-valued vector
#    y[] (for example, a vector containing invariant mass
#    squares). However, when we do the sampling over our events, it is
#    convenient to work with an array A_r(y) containing
#    model-dependent PWA amplitudes of our events. With other words,
#    to each event y_i = (m2_ab(i), m2_bc(i), ...), we would like to
#    compute complex numbers A_1(y_i), A_2(y_i), ... A_R(y_i) and
#    store them in a file that STAN can read - with other words, as
#    an array A_r_data dumped in the file
#    STAN_amplitude_fitting.data.R (or in the f_out argument).
#
#    To sum up, this script: opens a .ROOT file containing a Tree
#    with branches y.1 ... y.R (R meaning the number of variables);
#    gets the model-dependent amplitude function A_r(y) from model.so;
#    applies it to all events in the tree and dumps the result to 
#    *.data.R file.
#
# USAGE
#    data_analysis__root_to_dataR.py f_in f_out
#
#    Takes the ROOT file f_in ('output/generated_data.root' by default),
#    saves results in f_out ('output/STAN_amplitude_fitting.data.R' by 
#    default).

import argparse
import numpy as np
import ROOT
import os

import utils # convert_to_matrix - Translates A_cv results to usable form
             # stan_rdump - saves data to R dump file

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert a .ROOT file to STAN data.R file.')


parser.add_argument('f_in',
                     default='output/generated_data.root',
                     nargs='?',
                     type=argparse.FileType('r'))

parser.add_argument('f_out',
                    default='output/STAN_amplitude_fitting.data.R',
                    nargs='?',
                    type=argparse.FileType('w')
)


args = parser.parse_args()
MODEL_FOLDER = os.getcwdu()

# display what's happening!
print("data_analysis__root_to_dataR.py: Reading {0}.".format(args.f_in.name))

# Load the branches y from *.root into python cache
f_in = ROOT.TFile(MODEL_FOLDER + '/' + args.f_in.name)
t = f_in.Get("t")

# Declare the variables that will contain tree entries
y = np.asarray(0, dtype=float)

# Aliase the nave of the variable to the tree branch
branch_name = "y"
t.SetBranchAddress(branch_name, y)

# D_ is the number of events
D_ = t.GetEntries()
y_data_ = np.zeros(D_, dtype = float)


### FILL DATA FROM TREE ###
for d in range(D_):
    t.GetEntry(d)
    y_data_[d] = y
data = dict(N = D_, y = y_data_)

utils.stan_rdump(data, MODEL_FOLDER + '/' + args.f_out.name)
print("data_analysis__root_to_dataR.py: Done. Data dumped in {0}.".format(args.f_out.name))
