""" Note: this module now essentially just imports from AC_tools.

    remaining Prod/Loss functions will also be updated to use 
    AC_tools in future, making this module redundant. """
# ------------  MChem Tools - tms -------------------------------------
# --------------  
# ---- Section 0 ----- Modules required

# I/O
import os
import glob
import csv
import sys
from netCDF4 import Dataset
# import PyGChem ( + allow for back compatibility )
import pygchem
if pygchem.__version__ == '0.2.0':
    import pygchem.diagnostics as gdiag
else:
    # standard code or dev branch?
    try:
        from pygchem import datasets
    except:
        import pygchem.datafields as datasets


# Plotting/Analysis
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime as datetime 
import numpy as np

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Core Programmes
# 1.01 - Open ctm.bpch  ( pygchem version 0.2.0)
# 1.02 - Get np array (4D) of ctm.bpch ( lon, lat, alt, time)
# 1.03 - Get Air mass 
# 1.04 - Plot up GEOS-Chem slice
# 1.06 - Process time/date to CV days equivilent - mje
# 1.07 - incremental increase datetime by given months - credit: Dave Webb
# 1.08 - processes provided files to extract data/names
# 1.09 - date specific (Year,Month,Day) planeflight output reader 
# 1.10 - Get var data for model run ( pygchem version > 0.3.0)

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Generic Tools
# 2.01 - find nearest values in array
# 2.02 - Get GEOS-Chem longitude index for a given latitude
# 2.03 - Get GEOS-Chem latitude index for a given latitude
# 2.04 - Get GC datetime
# 2.05 - Get GEOS-Chem altitude index for a given latitude

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 3 ----- Prod-Loss processing
# 3.01 - extract reactions tracked by prod loss diag for a given p/l family
# 3.02- Get tags for reactions
# 3.03 - Get prod loss reactions for a given family.
# 3.04 - Extract reactions to form a dictionary of active reactions

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 4 ----- Species and Reference GEOS-Chem data
# 4.01 - Class for holding extra Geos-Chem (GC) species data.
# 4.02 - GEOS-Chem species in latex form
# 4.03 - What GEOS-Chem (GC) Species am i? takes TRA_## & returns GC ID - setup for tms iodine tracer # > 53 (v9-01-03/v9-2 - iodine branch)
# 4.04 - Retrieve lat, lon and alt for a given resolution and region.
# 4.05 - Get CTM (GEOS-Chem) array dimension for a given resolution
# 4.06 - Convert gamap category/species name to Iris/bpch name
# 4.07 - Returns tracers unit and scale (if requested)
# 4.99 - Reference data (from GChem - credit: Gerrit Kuhlmann)

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 5 ----- Core Programmes
# 4.01 - Planeflight dat file generation... 
# 4.02 -  


# --------------                                                                                 
# 1.01 - open ctm.bpch using PyGChem <= REDUNDENT
# --------------                                                                                
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import open_ctm_bpch

# --------------                                                                                 
# 1.02 - get np array (4D) of ctm.bpch ( lon,lat , alt,time) <= REDUNDENT                         
# --------------                                                                                
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_gc_data_np

# --------------
# 1.03 - get air mass (4D) numpy array <= REDUNDENT 
# -------------       
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_air_mass_np

# -----                                                                                                            
# 1.04 - plot geos slice                                                                                           
# -----                                                                                                            
# <= update to use map_plot from AC_Tools?
def map_plot(scalar, species=None, unit=None, res='4x5', **Kwargs):

    # Setup slices                                                                                                 
    # Grid/Mesh values for Lat, lon, & alt                                                                         
    lon, lat, alt = get_latlonalt4res( res=res )

    # Setup mesh grids                                                                                             
    x, y = np.meshgrid(lon,lat)
    print len(x), len(y)

    # Add Labels for axis + title  if species provided
    plt.ylabel('Latitude', fontsize = 20)
    plt.xlabel('Longitude',fontsize = 20)
    if not isinstance( species, type(None) ):
        plt.title( '{} / {} '.format( latex_spec_name( species ), unit ) )

    # Setup map ("m") using Basemap                                                                                
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-182.5,\
                    urcrnrlon=177.5,\
                    resolution='c')
    m.drawcoastlines()

    parallels = np.arange(-90,91,15)
    meridians = np.arange(-180,151,30)

    plt.xticks(meridians) # draw meridian lines
    plt.yticks(parallels) # draw parrelel lines
#    m.drawparallels(parallels) # add to map
                                                                                                                   
    # Create meshgrid to plot onto                                                                             
    x, y = np.meshgrid(*m(lon, lat))
    print len(x), len(y)

    plt.xlim(-180,175)
    plt.ylim(-89,89)

    poly = m.pcolor(lon, lat, scalar, cmap = plt.cm.Blues)

    # Add labels/annotations                                                                                       
    cb = plt.colorbar(poly, ax = m.ax,shrink=0.4)

    return plt , cb 

# --------------
# 1.06 - Process time/date to CV days equivilent - mje
# -------------
""" translate year to "since2006" function """
# Now import function directly from AC_tools 
from AC_tools.funcs4time import year_to_since_2006

# -------------
# 1.07 -  incremental increase datetime by given months - credit: Dave Webb
# -------------
# Now import function directly from AC_tools 
from AC_tools.funcs4time import add_months
    
# --------------
# 1.08 - processes provided files to extract data/names
# -------------    
# Now import function directly from AC_tools 
from AC_tools.funcs4pf import process_files_to_read

# --------------
# 1.09 - date specific (Year,Month,Day) planeflight output reader 
# -------------
# 
# Now import function directly from AC_tools 
from AC_tools.funcs4pf import readfile

# ----
# 1.10 - Get var data for model run
# ---
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_GC_output

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Generic Tools
# --------   
# 2.01 - Find nearest
# --------
# Now import function directly from AC_tools 
from AC_tools.funcs4generic import find_nearest

# ----
# 2.02 - Get Longitube in GC grid box number
# ---- 
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_gc_lon

# ----
# 2.03 - Get Latitude in GC grid box number
# ---- 
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_gc_lat

# ----
# 2.04 - Get gc datetime
# -----
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_gc_datetime

# -------------- 
# 2.05 - get_gc_alt  ( KM => box num )   
# -------------
# Now import function directly from AC_tools 
from AC_tools.funcs4GEOSC import get_gc_alt


# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 3 ----- Prod-Loss processing

# -------------
# 3.01 - extract reactions tracked by prod loss diag for a given p/l family
# ------------- 
# <= update to use AC_tools package
def rxns_in_pl( wd, spec='LOX' ):
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if all( [ i in row for i in 'Family','coefficient' , 'rxns', spec] ):
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'Family' 
    rxns = [ i for i in rxns if (  'Family' not in i ) ]
    n = [int(rxn[1]) for rxn in rxns ]
    rxns = [rxn[2:] for rxn in rxns ]

    rdict = dict( zip(n, rxns) )
    return rdict

# -------------
# 3.02 - Get tags for reactions
# ------------- 
# <= update to use AC_tools package
def get_p_l_tags( rxns ):

    # (PD??, RD??, LO3_??, PO3_??, LR??)
    for rxn in rxns:
#        print rxn
#        print [ i for i in rxn if  any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        tags =  [i for i in rxn if any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        try:
            tagsl.append( tags)
        except:
            tagsl = [tags]

    return tagsl

# -------------
# 3.03 - Get prod loss reactions for a given family.
# ------------- 
# <= update to use AC_tools package
def prod_loss_4_spec( wd, spec ):
    # ---  Get Dict of all reactions, Keys = #s
    rdict = rxn_dict_from_smvlog(wd)

    # ---  Get reaction # tracked by p/l diag for spec and coefficient.
    rxns = rxns_in_pl(wd, spec)
    nums =  rxns.keys() 
    Coe = [ rxn[-1] for rxn in rxns.values() ]

    # --- get all details from full reaction dictionary
    rxns =  [ rdict[i] for i in nums ]

    # --- get tags for tracked reactions, state where reactions are un tracked
    tags = get_p_l_tags( rxns )

    return nums, rxns, tags, Coe

# -------------
# 3.04 - extract reactions to form a dictionary of active reactions
# ------------- 
# <= update to use AC_tools package
def rxn_dict_from_smvlog( wd, spec='LOX' ):
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if 'NMBR' in row:
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'NMBR'
    rxns = [ i for i in rxns if (  'NMBR' not in i ) ]
    n = [int(rxn[0]) for rxn in rxns ]
    rxns = [rxn[1:] for rxn in rxns ]
    rdict = dict( zip(n, rxns) )
    return rdict
    


# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 4 ----- Species and Reference GEOS-Chem data

# ------------
# 4.01 - Class for holding extra Geos-Chem (GC) species data. Input the string into the class init function and it can give several forms of data out.
# e.g. 
#  O3 = species('O3')
#  print O3.Latex
#  print O3.RMM
#  Output: 'O_3'   
# -------------

class species:

   def __repr__(self):
        return "This is a class to hold chemical species information"
   def __str__(self):
        try: 
            print "name = " + self.name
        except:
            pass
        try: 
            print "GeosChem group = " + self.group
        except:
            pass
        try:
             print "formula = " + self.formula
        except:
            pass
        try:
             print "InchI = " + self.InChI
        except:
            pass
        try:
             print "smiles = " + self.smiles
        except:
            pass
        try:
             print "RMM = " + str(self.RMM)
        except:
            pass
        try:
             print "Latex = " + self.Latex
        except:
            pass
        return "Class for holding species data"

   def __init__(self, name):
      self.name = name
      self.help = ("""This is a class to get information on species from a local CSV folder
   It might contain the following information:
   self.RMM       = The Mean Mass of the species.
   self.latex     = The latex name of the species.
   self.smiles    = The smiles string of the species.
   self.InChI     = The InChI string of the species.
   """)
      species_filename = os.path.dirname(__file__) + "/Species.csv"

      try:
         species_file     = open(species_filename, 'rb')
      except IOError:
         print "Error: Species.csv does not appear to exist."
      species_csv      = csv.reader(species_file)

      if (name in ['OH', 'HO2']):
         self.group = 'CHEM-L=$'
      else:
         self.group = 'IJ-AVG-$'

      species_in_csv=False
      for row in species_csv:
         if (str(self.name) == row[0].strip()):
            species_in_csv=True
            self.formula   = row[1]
            self.InChI     = row[2]
            self.smiles    = row[3]
            self.RMM       = float(row[4])
            if row[5].isspace() or row[5] == "": # Check if the Latex is defined (not whitespace or empty)
               self.Latex = name
            else: 
              self.Latex     = row[5]
      if not species_in_csv:               
         print "Warning:"
         print "In MChem_tools.py class Species: Species " + name + " not found in species CSV file"
         return

# --------------
# 4.02 -  GEOS-Chem species in latex form
# --------------
# Now import function directly from AC_tools 
from AC_tools.funcs_vars import latex_spec_name 

# --------------
# 4.03 - What GEOS-Chem (GC) Species am i? takes TRA_## & returns GC ID - setup for tms iodine tracer # > 53 (v9-01-03)
# -------------
# Now import function directly from AC_tools 
from AC_tools.funcs_vars import what_species_am_i 

# ----
#  4.04 - returns lon/lat/alt for res
# ----
# Now import function directly from AC_tools 
from AC_tools.funcs4core import get_latlonalt4res

# --------   
# 4.05 - Get CTM (GEOS-Chem) array dimension for a given resolution
# --------
# Now import function directly from AC_tools 
from AC_tools.funcs4core import get_dims4res

# --------   
# 4.06 - Convert gamap category/species name to Iris/bpch name
# --------
# Now import function directly from AC_tools 
from AC_tools.funcs_vars import diagnosticname_gamap2iris

# --------------
#  4.07 - Returns tracer's unit (and scale if requested)
# --------------
# Now import function directly from AC_tools 
from AC_tools.funcs_vars import tra_unit

# --------------
# 4.99 - Reference data, (inc. grid data) from GChem 
# ------------- 
""" Updated to dictionary from gchemgrid (credit: Gerrit Kuhlmann ) with 
     addition grid adds """       
# Now import function directly from AC_tools 
from AC_tools.funcs4core import gchemgrid
                                                                  
                                                                                 

