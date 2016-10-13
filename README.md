

REDUNDENT. The AC_tools example code held here have been moved to AC_tools' folder of examples/Scripts (https://github.com/tsherwen/AC_tools/Scripts). All (and more) of "MChem_tools" functionality is present in AC_tools.

Please use AC_tools instead of MChem_tools.

This code is no longer maintained. 

MChem_tools
===========

Example analysis of chemical transport model (CTM) output for
 working with GEOS-Chem binary (ctm.bpch) and csv (planeflight.dat.YYYYMMDD). 


- MChem_bpch_plotter.py :
Plots up GEOS-Chem ctm.bpch as 2D maps using PyGChem.
This programme can alsoo calculate tropospheric burdens.  

- MChem_planeflight_plotter.py :
Plots up GEOS-Chem planeflight output as timeseries plots.

- bpch2netCDF.py : 
Programme to convert multiple bpch files to a single netCDF.  This is 
portable and can be called as post run script or directly at command line.

- Prod_Loss_4_spec.py :
Tools to allow for automated analysis of prod/loss species in by reading
globchem.dat and smv2.log files. This allows for debugging/checking of 
tagged simulations

- MChem_tools.py :
Module containing functions for programmes above. Note: this module 
now mostly just imports from the AC_tools submodule and will be 
removed in the near future. 


Requires:
PyGChem: https://github.com/benbovy/PYGChem

Scientific python modules e.g. Annoconda ( https://store.continuum.io/cshop/anaconda/  )


Notes for use:

As this repository uses a submodule (AC_tools), please recursively clone URL: 

git clone --recursive https://github.com/tsherwen/MChem_tools.git

And to update:

git submodule update --recursive

To be able to access these modules from any place, add the following to your $HOME/.bashrc file:
Where MChem_tools is the location you downloaded the MChem_tools folder.

To add MChem_tools to python path, add the line below to bashrc
export PYTHONPATH=${PYTHONPATH}:$HOME/MChem_tools

This module and the submodule (AC_tools ) it uses are setup to be portable, 
so when using functions contained within these modules a direct import into 
a new script would be recommended.

e.g. 

from AC_tools.funcs4time import *

or

from MChem_tools.MChem_tools import *


Monthly Run
===========

Script to split GEOS-Chem run into monthly jobs is in directory 'monthly_run'. 
See sepereate README in dricetory for more details.
