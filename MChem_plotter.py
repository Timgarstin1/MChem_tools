# modules
from MChem_tools_II import *
import numpy as np

# setup
species = 'O3'#'CO2'
RMM_species = 12 + 16 + 16
wd= '/home/ts551/labbook/all_model_simulations/iodine_runs/iodine_module_WIP_14_01_01_1_year_I_Br_all_I_org_II_I2_r_p_l_jones_MUTD_end_at_I2O2/run'
wd = '/home/tt561/MChem/SMV_run/run'

# get data  ( lon,lat , alt,time) 
ctm_f = open_ctm_bpch(wd) 
data = get_gc_data_np(ctm_f, species)
air_mass = get_air_mass_np(ctm_f)  

# calc the total mass of a speices - select data you want to calc & print
air_moles = (air_mass*1E3) / (.8*28+.2*32)  # mass (convert to g) / RMM air
species_moles = air_moles * data            # moles * v/v
species_mass = species_moles * RMM_species
print species_mass.shape
print np.sum(species_mass)

#  select data you want to plot
# mean () of time                                                                                                                                            
data = data[:,:,:,:].mean(axis=3)

# select surface
print data.shape
data = data[:,:,0]
print data.shape
data = np.transpose(data)
print data.shape

# plot surface
plt, cb = plot_geos_alt_slice(data)

plt.show()





