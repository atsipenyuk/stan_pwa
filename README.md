# stan_pwa
CmdStan module for meson decay PWA

The project <b>stan_pwa</b> is an extension module to CmdStan (https://github.com/stan-dev/cmdstan). 
Its purpose is to simplify modelling partial wave analysis 
of heavy meson physics. We implement such functions as
Breit-Wigner/Flatte dynamical form factors, 
Blatt-Weisskopf functions and Zemach tensors
and bundle them together to model-dependent PWA amplitudes 
that describe 3- and 4-body meson decays.  
  
We introduce only a handful of functions to Stan language;
the core of the module is build from C++ functions that are not
directly callable from Stan (however, they are templated and 
can be exposed to Stan parser, should the need present itself).

The module also contains some small python scripts and modules
that are aimed at debugging the C++ code and plotting/analyzing
the results of the Stan fitting.
  
This code may be interesting to you if you want to:  

 * See an example how to expose your functions to Stan (look at the `install` and `reload_libraries` functions of the `makefile`);
 * Fit a function that looks like "f(y,theta) = |A(y) * theta|^2", where "A(y)" and "theta" are complex vectors (look at `models/example_models/two_toy_res/src/model.hpp`);
 * Use complex numbers in Stan (look at `src/complex`);
 * Use some model-dependent PWA functions (templated C++: `src/fct`).

### Licensing

The core Stan C++ code and CmdStan are licensed under the new BSD. For more information, see
https://github.com/stan-dev/cmdstan/blob/develop/LICENSE. This project is licensed under 
the MIT license.

### Download and install 
Download CmdStan-2.9.0. Unzip it, download and install stan_pwa files:
```bash
 ..$ unzip cmdstan-2.9.0.zip
 ..$ cd cmdstan-2.9.0
 ..$ git clone https://github.com/atsipenyuk/stan_pwa.git
 ..$ cd stan_pwa
 ..$ make -s install
```
You also need to add the following line to your .bashrc file:
```bash
export PYTHONPATH=/your_path_to_stan_pwa/stan_pwa/lib/py:$PYTHONPATH
```  

### Dependencies 
The module requires following dependencies:

 * libboost-python-dev;  
 * PyROOT (you should be able to call import ROOT during a python session).
 * numpy

The PWA model is written in C++; however, it is often more convenient to python for tasks not directly related to sampling. Such tasks may include data analysis or calculation of the LF (likelihood function) normalization integrals. This is why wo use the boost library: it allows conversion of the c++ code to python modules. To get libboost-python-dev, run  
```bash
 ..$ sudo apt-get install libboost-python-dev
```
or the corresponding command for your system.  

The PyROOT package is used to wrap Stan output files to Root trees and vice versa. Check out the [installation guide](https://root.cern.ch/drupal/content/pyroot). Note that on E18 Linux, it suffices to check that the following line
```bash
export PYTHONPATH=$PYTHONPATH:/nfs/mnemosyne/sys/slc6/sw/root/x86-64/5.34.21/root/lib
```
is present in your .bashrc file.

### Example
The following commands demonstrate how to generate data and sample one complex parameter modeling D->3pi decay via two fictitious Breit-Wigner resonances.
```bash
..$ cd cmdstan-2.9.0/stan pwa/models/two_toy_res  
..$ ./../../../relink model.sh # Use the correct model.hpp file  
..$ ./../../../build.sh # Build executable files  
..$ ./../../../generate.sh 100000 # Generate 100 000 events  
..$ root output/generated data.root # You may check the generated data    
root [1] t->Draw("y.2:y.1>>hh(100,0,3,100,0,3)","","COLZ",20000, 0);  
root [2] .q  
..$ # To evaluate amplitudes, python library of our model must be made  
..$ ./../../../wrap python.sh # Creates build/model.so  
..$ # Calculate normalization integrals; bounds of Dalitz plot are 0 and 3 in both axes  
..$ ./../../../calculate_normalization_integrals.py 0 3 0 3  
..$ # Evaluate amplitudes and Monte Carlo integrals  
..$ ./../../../prepare_for_fitting.sh  
..$ # Fitting 1000 warmup + 1000 samples using this model  
..$ # with 100 000 data pts requires ca. 45 min. on a home computer  
..$ ./../../../fit.py -c 2 # Run two sampling chains, 1000 samples  
..$ ./../../../merge_output_chains.sh    
..$ root output/output.root  
root [1] t->Draw("y.2:y.1","","COLZ", 20000, 0);  
root [2] .q  
```
