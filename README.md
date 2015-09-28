MChem_tools
===========

Tools for MJE MChems

- MChem_bpch_plotter.py

GEOS-Chem ctm.bpch plotter using PyGChem 

- MChem_planeflight_plotter.py

GEOS-Chem planeflight output plotter 

- bpch2netCDF.py

Tool for converting multiple bpch files to a single netCDF 

- Prod_Loss_4_spec.py 

GEOS-Chem smvgear analyser of prod/loss for reaction tag or "spec" 

- MChem_tools.py

Module containing functions for programmes above 



Requires:

PyGChem: https://github.com/benbovy/PYGChem

Scientific python modules e.g. Annoconda ( https://store.continuum.io/cshop/anaconda/  )


Notes for use:

To be able to access these modules from any place, add the following to your $HOME/.bashrc file:
Where MChem_tools is the location you downloaded the MChem_tools folder.

To add MChem_tools to python path, add the line below to bashrc
export PYTHONPATH=${PYTHONPATH}:$HOME/MChem_tools

Monthly Run
===========

Script to split GEOS-Chem run into monthly jobs is in directory 'monthly_run'. See sepereate README in dricetory for more details.