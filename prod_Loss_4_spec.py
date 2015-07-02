from MChem_tools import *
import re
import sys

# ---- Get inputs from command line/defaults
try:    # chcck if a directory was given ad command line
    wd    = sys.argv[1]
except: # Otherwise use path below
    wd    = '<insert GEOS-Chem run direcotory path here>'

# --- Settings
#spec = 'LOX' 
spec = 'POX' 
#spec =  'PIOx'#'LIOx'#'L2OI'#'L_Iy' #'P_Iy'
nums, rxns, tags, Coe = prod_loss_4_spec( wd,  spec )
rdict =  rxn_dict_from_smvlog( wd )
diagnose = False
#diagnose = True

# ---  compare with tracked families in tms diags
for n, rxn in enumerate( rxns ):
    print n, nums[n], rxn[4:][:5], tags[n], Coe[n]
    # -- ouput rxns without tags
    if diagnose:
        if len( tags[n] ) < 1 :
            print 'tracking needed for this reaction: {}'.format( nums[n]  )
            try:
                extra.append( nums[n]  )
            except:
                extra = [ nums[n] ]
                
# ---  print detail on untracked reactions to screen 
if diagnose:            
    print [ ( len(i), i ) for i in [extra]]
    for rxn in extra:
        print rxn , rdict[rxn]
