MChem_tools
===========

Tools for MJE MChems

- GEOS-Chem ctm.bpch plotter using PyGChem (MChem_plotter.py)

- GEOS-Chem planeflight output plotter ( MChem_planeflight_plotter.py )

- Tool for converting multiple bpch files to a single netCDF (bpch2netCDF.py)

- Module containing functions for programmes above ( MChem_tools.py )

<<<<<<< HEAD
=======
- GEOS-Chem smvgear analyser of prod/loss for reaction tag or "spec" ( Prod_Loss_4_spec.py ) 

Requires:

PyGChem: https://github.com/benbovy/PYGChem

>>>>>>> 7d8b30732a6eb73eee495b0ce085eef1af406aa1

To be able to access these modules from any place, add the following to your $HOME/.bashrc file:
Where MChem_tools is the location you downloaded the MChem_tools folder.

# Add MChem_tools
export PYTHONPATH=${PYTHONPATH}:$HOME/MChem_tools

