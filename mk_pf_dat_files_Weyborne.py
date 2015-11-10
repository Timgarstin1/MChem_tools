# --- Import packages
from AC_tools.funcs4GEOSC import *
from AC_tools.funcs4pf import * 
import numpy as np
from time import gmtime, strftime
import time
import glob

# --- Settings
#wd = '/work/home/ts551/data/'#IO_obs'
GAWfile = 'Weyborne.dat'
start_year, end_year = 2004,2006
debug = False
# Version - Nic, this is the version you are using for 0.25*0.3125
ver = '1.7'  #ver = '1.6.1'

# --- Read in site ('GAWfile') Detail
numbers, lats, lons, pres, locs = readin_gaw_sites( GAWfile )
# make sure numbers are floats
lats, lons, pres = [ np.float64(i) for i in lats, lons, pres ]
print lats[0:4]

# --- Set Variables
slist = pf_var('slist_REAs_all_OH_extras', ver=ver )
MUTDwd =  MUTD_runs()[0]
rdict = rxn_dict_from_smvlog(MUTDwd)

nvar=len(slist)
yr = range(start_year, end_year )
m  = range(01,13)
da  = range(01,32,1)
h  = range(00,24,1)
mi = range(00,60,1)
minute = 0

# --- loop dates and make Planeflight log files for given points
# 
end_str = '999999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
loc_str = '{:>6}   {:<4}{:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
for year in yr:
    for b in m: 
        for c in da:

            if debug:
                print year, b, c
            a=open('Planeflight.dat.'+str(year)+'{:0>2}{:0>2}'.format( b, c), 'w')  

            print >>a, 'Planeflight.dat -- input file for ND40 diagnostic GEOS_FP'
            print >>a, 'Tomas Sherwen'
            print >>a, strftime("%B %d %Y", gmtime())
            print >>a, '-----------------------------------------------'
            print >>a, '{:<4}'.format(nvar)        ,'! Number of variables to be output'
            print >>a, '-----------------------------------------------'

            for i in range(0,len(slist)):
                print >>a, slist[i]
            
            print >>a, '-------------------------------------------------'
            print >>a, 'Now give the times and locations of the flight'
            print >>a, '-------------------------------------------------'
            print >>a, 'Point   Type DD-MM-YYYY HH:MM     LAT     LON   PRESS'
                    
            counter=0
            for d in h:
                
                for i in range(len(lats)):
                    print >>a, loc_str.format( counter,locs[i], c, b,  year, d, minute, lats[i], lons[i], pres[i] )
                    counter+=1
                   
            print >>a, end_str
            a.close()

