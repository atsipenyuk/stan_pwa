#!/usr/bin/env python
# fit.py
#
# NAME
#    fit.sh - fit parameters to the Stan model.
#
# SYNOPSIS
#     fit.sh [OPTION]... [STAN_EXEC_FILE]... [OUTPUT_FILE]..
#
# DESCRIPTION
#     Fit the model using build/STAN_amplitude_fitting and convert results
#     to a .root file.
#
#     -I, --input_file
#         exec file to be called (by default, build/STAN_amplitude_fitting
#
#     -O, --output_file
#         .csv output file filled with samples (defaults to output/output.csv)
#
#     -c, --chains=INT
#         how many sampling chains should be ran at parallel: for example,
#         if you have 3 free cores, you may run 3 parallel chains in back-
#         ground (defaults to 0, meaning that there is one chain and it
#         does not run in background)
#
#     -w, --warmup_only
#         starts one chain with random initial values; generates 1000
#         samples, but keeps only the last sample that is saved to
#         output/STAN_amplitude_fitting.init.data.R
#
#     -i, --init_data=FILENAME
#         starts the chain with initial values specified in FILENAME;
#         by default, starts with random initial values
#
#     -d, --data_file=FILENAME
#         sampling data file, defaults to 
#         output/STAN_amplitude_fitting.data.R
#
#     -n, --num_samples=INT

import argparse
import os
import subprocess

import utils

# Parse the arguments
parser = argparse.ArgumentParser(description='Fit parameters using sampler build/STAN_amplutide_fitting')

parser.add_argument('--input_file','-I',
                    default="build/STAN_amplitude_fitting")

parser.add_argument('--output_file','-O',
                    default="output/output.csv")

parser.add_argument('--chains','-c',
                    default=0,
                    type=int)

parser.add_argument('--num_samples','-n',
                    default=1000,
                    type=int)

parser.add_argument('--warmup_only','-w',
                    action='store_true')

parser.add_argument('--init_data','-i',
                    default='none_given')

parser.add_argument('--data_file','-d',
                    default="output/STAN_amplitude_fitting.data.R")

args = parser.parse_args()

model_path = os.getcwd()
stan_pwa_path = utils.get_path('stan_pwa')

# Warmup only?
if args.warmup_only:
    # Sample
    subprocess.Popen("./" + args.input_file + \
                          " sample num_warmup=" + str(args.num_samples) + \
                          " num_samples=1 " + \
                          " data file=" + args.data_file + \
                          " output file=" + args.output_file + " " + \
                          " diagnostic_file=" + args.output_file[:-4] + \
                          "_diagnostics.csv", shell=True).wait()
    # Convert to *.init.data.R
    subprocess.Popen(os.path.join(stan_pwa_path,'bin/csv_to_init_data_R.py ')+ \
                     args.output_file + ' ' + \
                     args.data_file[:-7] + '.init.data.R' + '\n', shell=True).wait()

# OK, real sampling
else:
    # Initial values - defined or random
    if args.init_data=='none_given':
        init = ''
        num_warmup = 1000
    else:
        init = ' init=' + args.init_data + ' '
        num_warmup = 1000

    # Single chain in foreground
    if args.chains==0:
        # Sample
        subprocess.Popen("./" + args.input_file + \
                             " sample num_warmup=" + str(num_warmup) + \
                             " num_samples=" + str(args.num_samples) + \
                             " data file=" + args.data_file + \
                             " output file=" + args.output_file + \
                             init + "\n", shell=True).wait()
        # Translate to *.root
        subprocess.Popen(os.path.join(stan_pwa_path,'bin/csv_to_root.py ') + \
                             args.output_file + ' ' + \
                             args.output_file[:-4] + '.root \n',shell=True).wait()

    # Multiple chains in background
    # No more than 9 chains though (Stan is fine with it, but
    # correct the merge outputs part, if you need more than 9 chains).
    else:
        # Sample
        subprocess.Popen('for i in {1..' + str(args.chains) + '}\n' + \
                             'do \n' + \
                             "./" + args.input_file + \
                             " sample num_warmup="+ str(num_warmup) + \
                             " num_samples=" + str(args.num_samples) + \
                             ' id=$i ' \
                             " data file=" + args.data_file + \
                             " output file=" + args.output_file[:-4] + \
                             "$i.csv" + \
                             init + " &\n" + \
                             'done \n', shell=True).wait()
        ## In plain bash, the command above translates to following:
        # 
        # for i in {1..1}
        # do
        # ./build/STAN_amplitude_fitting sample num_warmup=1000 num_samples=1000 id=$i data file=output/STAN_amplitude_fitting.data.R output file=output/output$i.csv init=output.STAN_amplitude_fitting.init.data.R &
        # done
        #
        ##

        # Merge chain outputs (DEPRECATED) - is called too early
        #subprocess.Popen("grep lp__ " + args.output_file[:-4] + "1.csv" + \
        #                 " > " + args.output_file, shell=True).wait()
        #subprocess.Popen("sed '/^[#l]/d' " + args.output_file[:-4] + \
        #                 "?.csv >> " + \
        #                 args.output_file,shell=True).wait()

        # Translate to *.root
        #subprocess.Popen(os.path.join(stan_pwa_path,'bin/csv_to_root.py ') + \
        #                     args.output_file + ' ' + \
        #                     args.output_file[:-4] + '.root \n',shell=True).wait()

