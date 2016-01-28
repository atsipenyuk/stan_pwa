#!/usr/bin/env python
# csv_to_root_w_renormalized_background.py

# NAME
#    csv_to_root_w_renormalized_background.py - make a ROOT tree from the 'output.csv' STAN file
#
# SYNOPSIS
#    ./csv_to_root_w_renormalized_background.py OUTPUT_CSV OUTPUT_RENORMALIZED_BACKGROUND_ROOT[OPTIONS...]
#
# DESCRIPTION
#    csv_to_root.py parses the specified file (by default, 
#    'generated_data.csv'), ignoring all lines starting with '#', and writes 
#    the entries of the file in a tree. Also, we introduce a 3-dim. symplectic vector b, 
#    containing the renormalized background entries. Specifically,
#      b.2 = b0 * theta_background_flat_abs2,
#      b.3 = b0 * theta_background_rho_770_abs2,
#      b.1 = 1 - b.2 - b.3,
#    and b0 = (1 + theta_background_flat_abs2 + theta_background_rho_770_abs2)^{-1}.




import argparse
import numpy as np
import os
import ROOT


# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert STAN output.csv to a .ROOT file.')

parser.add_argument('--tree_name', 
                    default='t',
                    help="Name of the tree in the .ROOT file.")

parser.add_argument('--tree_title', 
                    default='t',
                    help="Title of the tree in the .ROOT file.")

parser.add_argument('f_in', 
                     default='generated_data.csv',
                     nargs='?',
                     type=argparse.FileType('r'))

parser.add_argument('f_out', 
                    default='generated_data.root',
                    nargs='?',
                    type=argparse.FileType('w')
)

parser.add_argument('-p', '--print_dalitz_plot',
                    default=0,
                    type=float)

args = parser.parse_args()


#
print('csv_to_root.py: Reading {0}...'.format(args.f_in.name))

# Open CSV file
#print('Importing data...')
f_in = args.f_in

# Open ROOT file
f_out = ROOT.TFile(args.f_out.name, "recreate")
t = ROOT.TTree(args.tree_name, args.tree_title)


for line in f_in:
    # Ignore lines containing STAN commands and STAN info...
    if line[0] in ['#', '\n', ' ']:
        pass

    # Parse the line containing names of the parameters...
    elif line[0] == 'l':
        param_names = line.split(",")
        param_names = [x.replace("\n", "") for x in param_names if x != '']
        # Add the extra branches for b.1, b.2, b.3
        param_names = param_names + ["b.1", "b.2", "b.3"]
        # Parameter values will be stored in param[0], param[1], ... etc.
        param = [np.zeros(1, dtype=float) for i in range(len(param_names))]
        param = np.asarray(param)
        for i in range(len(param_names)):
            t.Branch(param_names[i], param[i], param_names[i] + '/D')

    # Parse the lines containing parameter values...
    else:
        a = line.split(",")
        a = [float(x) for x in a if x != '']
        # Add renormalized background
        # The background parameters are the two last parameters
        b0 = 1.0 / (1.0 + a[-1] + a[-2])
        b_1 = b0
        b_2 = b0 * a[-2]
        b_3 = b0 * a[-1]
        a = a + [b_1, b_2, b_3]
        param[:,0] = [x for x in a]
        t.Fill()

# Make a *.pdf drawing of the parameters 'm2_ab', 'm2_bc'
num_bins = str(200)
if args.print_dalitz_plot == 1:
   c = ROOT.TCanvas("y.1:y.2", "Dalitz plot")
   t.Draw("y.1:y.2>>hh(" + num_bins + ", 0, 3.1, " + num_bins + ", 0, 3.1)", "", "COLZ")
   c.Print(args.f_out.name[:-5] + '.pdf' )


f_out.Write()
f_out.Close()

f_in.close()


print("csv_to_root.py: Done. Data saved in {0}.".format(args.f_out.name))
