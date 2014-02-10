## --------------------------------- Planeflight Plotter - tms -------------------------------------
# --------------  
import matplotlib.pyplot as plt 
import matplotlib.font_manager as font_manager
from matplotlib import cm
from MChem_tools_II import *

#---------------------------------------- SET PLOTTING HERE ----------------------------------------
# Plot out data?
#print_diags=True

# Where are the model files? And What Species do you want to plot?
wd='/home/jts507/MChem/2yearrun/Run/plane_flight_logs/plane.log.200*'

# Which species to plot?  (must be in list form) 
species_to_plot= 'ClNO2'#'Cl'

# What years, months, days to plot?  (must be in list form) 
years_to_use, months_to_use, days_to_use  = ['2009'], ['02'], [ "{0:0>2}".format(i)  for i in range(10, 20, 1)  ] # set day_to_use by adjusting range
print years_to_use, months_to_use, days_to_use 

# Which locations? (where 10 = 1000 Hpa, 08 = 800 Hpa etc ... ) (must be in list form) 
locations=['TX1','LA1'] # must be as a list of strings
print locations

# ---------------------------------------- START PLOTTING HERE -------------------------------------------------------------
# -------------
species_to_plot=[what_species_am_i(species_to_plot)]

fig = plt.figure(figsize=(15,6), dpi=80, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 8})

for i,site in enumerate(locations):
    model, names = readfile( sorted(glob.glob(wd)), site, years_to_use, months_to_use, days_to_use)
    k=names.index(species_to_plot[0])
    plt.plot(  year_to_since_2006(model), model[:,k]*1e9*1E3,color=cm.jet(1.*i/len(locations)),label='{0}'.format(locations[i])  )    

plt.xlabel('Time/ CV days')
plt.ylabel('Conc./ p.p.t.v.')
plt.legend(loc='upper right',prop=font_manager.FontProperties(size=10))
plt.grid(b=None, which='major', axis='both', alpha=0.3)
plt.show()
