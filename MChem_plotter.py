# modules
from MChem_tools_II import *
import numpy as np

# setup
species       = 'O3'#'CO2'
RMM_species   = 16.*3.
wd            = '/work/home/ts551/data/all_model_simulations/standard_code_runs/standard_run_code_13_01_17_1_year/run'
#wd           = '/home/tt561/MChem/SMV_run/run'

# get data  ( lon,lat , alt,time) 
ctm_f         = open_ctm_bpch(wd) 
data          = get_gc_data_np(ctm_f, species)
air_mass      = get_air_mass_np(ctm_f)  

# calc the total & mean mass of a speices - select data you want to calc & print
air_moles     = (air_mass*1E3) / (.78*28.+.22*32.)  # mass (convert to g) / RMM air
species_moles = air_moles * data            # moles * v/v
species_mass  = species_moles * RMM_species
print species_mass.shape
print np.sum( species_mass )  / 1E12 # Convert to Tg
print np.sum( np.mean(species_mass, axis=3) ) / 1E12 # Convert to Tg

#  select data you want to plot -  mean() of time                                                                                                                                            
data          = data[:,:,:,:].mean(axis=3)

# select surface
print data.shape
data          = data[:,:,0]
print data.shape
data          = np.transpose(data)
print data.shape

# plot surface
plt, cb       = plot_geos_alt_slice(data)

plt.show()
