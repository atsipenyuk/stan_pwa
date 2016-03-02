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

Reminder: the core Stan C++ code and CmdStan are licensed under the new BSD.

### Download and install 
Download CmdStan-2.6.2. Unzip it, download and install stan_pwa files:
```bash
 ..$ unzip cmdstan-2.6.2.zip
 ..$ cd cmdstan-2.6.2
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

