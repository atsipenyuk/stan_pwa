#!/usr/bin/env python
#
# py_wrapper_setup.py
#
# NAME
#     py_wrapper_setup.py - build python module from a c++ file using Boost library
#
# SYNOPSIS
#     py_wrapper_setup.py cmd [OPTIONS]
#
# DESCRIPTION
#     Uses boost library to build a python module from c++ file.
#     Standard parameters that are accepted by distutils.core setup 
#     scripts apply. The names of the source/output files are specified
#     in the code below

import argparse
import os
from distutils.core import setup
from distutils.extension import Extension

import utils

# Change these strings if you want to change the name of the input/output file
input_file='src/py_wrapper.cpp'
output_file='model'

# Change these lines (e.g. to "g++") if you want to use another c++ compiler
os.environ["CC"] = "clang++"
os.environ["CXX"] = "clang++"

# set cmdstan directory (as the first directory containing '/stan_pwa/')
model_path = os.getcwd()
cmdstan_path = os.path.split(utils.get_path('stan_pwa'))[0]

setup(name="Model_Dep_Functions",
      ext_modules=[
          Extension(output_file, [input_file],
          libraries = ["boost_python"],
          include_dirs=[cmdstan_path + "/stan_2.9.0/lib/stan_math_2.9.0",
                        cmdstan_path + "/stan_2.9.0/lib/stan_math_2.9.0/lib/eigen_3.2.4",
                        #                    cmdstan_path + "/stan_2.9.0/lib/stan_math_2.9.0/boost_1.58.0",
                        #                    cmdstan_path + "/stan_2.9.0/src",
                        cmdstan_path],
          undef_macros=['NDEBUG']
          )
      ])
