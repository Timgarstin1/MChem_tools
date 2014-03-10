# ----------
import numpy as np
import glob
import datetime as datetime 
import csv

#  processes provided files to extract data/names
def process_files_to_read(files, location, big, names):
    print 'process_files_to_read called'
    print files
    reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace = True)
    for row in reader:
#        print location , (row[0] == 'POINT'), (row[1] == location) ,  len(row) , len(big), len(names)#, big.shape#, (row[2:])[0], (row[-1]) , l
        if row[1] == location: 
            new=row[2:]
            try:    
                big.append(new)
            except:
                big=[new]
        if row[0] == 'POINT':
            names = row[2:]
    return big, names


# --------------
# 1.01 - date specific (Year,Month,Day) planeflight output reader - tms
# -------------

def readfile(filename, location,  years_to_use, months_to_use, days_to_use, plot_all_data=False,debug=True, **kwargs):
    print 'readfile called'
    big, names = [],[]
             # sort for choosen years/months
    for files in filename:
             # loop through list on choosen years
        if (not plot_all_data):
            lll = 0
            for year in range(len(years_to_use)):
                if (("{0}".format(years_to_use[year])) in (files)) :
                    # is it not the last year specificied?
                    if (debug):
                        print 'years_to_use[year]', years_to_use[year], 'years_to_use[-1]', years_to_use[-1]
                    if (not (years_to_use[year] == years_to_use[-1])):
                        # just read all years upto point uptyo final year
                        big, names=process_files_to_read(files, location,big, names)
                        print 'i got to line 91'
                    # If last year selected, then only plot the given months & days
                    if (years_to_use[year] == years_to_use[-1]):
                        # Plot months exceot last one
                        for month in range(len(months_to_use)):                                                                                                  
                            if (debug):
                                print 'months_to_use[month]', months_to_use[month], 'months_to_use[-1]', months_to_use[-1], 'months_to_use', months_to_use, 'type(months_to_use)', type(months_to_use)
                            if (("{0}{1}".format(years_to_use[year],months_to_use[month])) in files) :  
                                if (not (months_to_use[month] == months_to_use[-1])):
                                    big, names=process_files_to_read(files, location,big, names)
                                    print 'i got to line 100',  'month',month,'in',len(months_to_use),  'year', year, 'in' , len(years_to_use)
                                if (months_to_use[month] == months_to_use[-1]):
                                        # For last month, plot days upto last day
                                    for day in range(len(days_to_use)):                                                                                          
                                        if (("{0}{1}{2}".format(years_to_use[year],months_to_use[month],days_to_use[day])) in files) : 
                                            if (debug):
                                                print 'days_to_use[day]', days_to_use[day], 'days_to_use[-1]', days_to_use[-1]
                                            big, names=process_files_to_read(files, location,big, names)
                                            if (debug):
                                                print 'i got to line 108'
                                                print 'readfile read big of size: ', len(big)
                                                
        if (plot_all_data):
            big, names=process_files_to_read(files, location,big, names)
            print 'reading all data'

    big=np.float64(big)             
    print 'readfile read big of size: ', len(big)
    return big, names

# --------------
# 1.02 - Process time/date to CV days equivilent - mje
# -------------
# translate year to "since2006" function
def year_to_since_2006(model):
            year=(model[:,0]//10000)
            month=((model[:,0]-year*10000)//100)
            day=(model[:,0]-year*10000-month*100)
            hour=model[:,1]//100
            min=(model[:,1]-hour*100)
            doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                        np.int(hour[i]),np.int(min[i]),0)- \
                      datetime.datetime(2006,1,1,0,0,0) \
                      for i in range(len(year))]
            since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
            return since2006

# --------------
# 1.03 - What GEOS-Chem (GC) Specie am i? takes TRA_## & returns GC ID or other wayround 
# -------------
#      1st - gc_to_tra = True (GC to TRA_##) or (TRA_## to GC) - HASHED OUT
#        2nd - input specie
def what_species_am_i(input_species) :
#(gc_to_tra)
#            if gc_to_tra :
# tracer library
            tracer_library={'O3':'O3','CO':'CO','NO':'NO','TRA_68': 'I2O', 'TRA_79': 'BrCl', 'TRA_71': 'I2O4', 'TRA_70': 'I2O3', 'TRA_59': 'IONO', 'TRA_58': 'HI', 'TRA_75': 'Cl', 'TRA_74': 'Cl2', 'TRA_77': 'ClO', 'TRA_76': 'HOCl', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_50': 'BrNO3', 'TRA_57': 'OIO', 'TRA_56': 'IO', 'TRA_55': 'HOI', 'TRA_54': 'I2', 'TRA_69': 'INO', 'TRA_62': 'CH3I', 'TRA_63': 'CH2I2', 'TRA_60': 'IONO2', 'TRA_61': 'I2O2', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_64': 'IBr', 'TRA_65': 'ICl', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_73': 'AERI', 'TRA_72': 'I2O5', 'TRA_78': 'OClO', 'TRA_66': 'I', 'TRA_67': 'HIO3','CH3Br': 'REA_53', 'HOBr': 'REA_47', 'CHBr3': 'REA_51', 'Br2': 'REA_44', 'BrO': 'REA_46', 'Br': 'REA_45', 'CH2Br2': 'REA_52', 'BrNO2': 'REA_49', 'BrNO3': 'REA_50', 'HBr': 'REA_48','NO': 'NO', 'O3': 'O3', 'IO': 'TRA_56', 'CH3Br': 'TRA_53', 'HI': 'TRA_58', 'Br': 'TRA_45', 'IONO': 'TRA_59', 'Cl': 'TRA_75', 'BrO': 'TRA_46', 'HIO3': 'TRA_67', 'OClO': 'TRA_78', 'CH3I': 'TRA_62', 'CHBr3': 'TRA_51', 'ClO': 'TRA_77', 'I2O2': 'TRA_61', 'HOI': 'TRA_55', 'BrNO2': 'TRA_49', 'BrNO3': 'TRA_50', 'I': 'TRA_66', 'I2O': 'TRA_68', 'OIO': 'TRA_57', 'Cl2': 'TRA_74', 'BrCl': 'TRA_79', 'CH2Br2': 'TRA_52', 'ICl': 'TRA_65', 'CH2I2': 'TRA_63', 'IBr': 'TRA_64', 'I2O5': 'TRA_72', 'CO': 'CO', 'HBr': 'TRA_48', 'HOCl': 'TRA_76', 'HOBr': 'TRA_47', 'Br2': 'TRA_44', 'I2': 'TRA_54', 'I2O4': 'TRA_71', 'AERI': 'TRA_73', 'IONO2': 'TRA_60', 'I2O3': 'TRA_70', 'INO': 'TRA_69','REA_327': 'HO + CH3IT => H2O + I (CH2I)', 'REA_378': 'HI => .5I2', 'REA_325': 'IO + CH3O2 =(M)> I + HO2 + HCHO', 'REA_324': 'OIO + OH => HIO3', 'REA_373': 'I2O3 + I2O4 => 4AERI', 'REA_322': 'IO + IO (O2) =>I2O2', 'REA_363': 'I2O2 =(M)> IO + IO', 'REA_320': 'IO + IO (O2) => I + OIO', 'REA_309': 'I + O3 => IO + O2', 'REA_362': 'OIO + OIO => I2O4', 'REA_360': 'IO + OIO =(M)> I2O3', 'REA_367': 'I2O2 + O3 => I2O3 +O2', 'REA_369': 'I2O4 + O3 => I2O5 + O2', 'REA_365': 'I2O2 =(M)> OIO + I', 'REA_328': 'I + NO <=(M)> INO', 'REA_368': 'I2O3 + O3 => I2O4 +O2', 'REA_444': 'IONO =(hv)> I + NO2', 'REA_370': 'I2O4 =(M)> 2OIO', 'REA_446': 'I2O2 => IO + IO', 'REA_440': 'I2 =(hv)> 2I', 'REA_451': 'INO => I + NO', 'REA_447': 'CH3IT =>  <CH3 +I>', 'REA_441': 'HOI =(hv)> I + OH', 'REA_448': 'CH2I2 => <CH2 + I +I >', 'REA_332': 'I + NO2 <=(M)> IONO', 'REA_379': 'IONO2 => 0.5I2', 'REA_449': 'IBr =(hv)> I + Br', 'REA_445': 'IONO2 =(hv)> I + NO3', 'REA_375': 'I2O4 + I2O4 => 4AERI', 'REA_352': 'IO + BrO => Br + I + O2', 'REA_353': 'IO + BrO => Br +OIO', 'REA_350': 'IO + ClO => I + Cl + O2', 'REA_351': 'IO + ClO => I + OClO', 'REA_356': 'ICl + OH => HOCl  +I', 'REA_357': 'IBr + Br => I + Br2', 'REA_354': 'ICl + Cl => Cl2 + I', 'REA_355': 'ICL + Br => BrCl + I', 'REA_358': 'IBr + OH => HOI + Br', 'REA_359': 'IBr + OH => HOBr + I', 'REA_372': 'I2O3 + OIO => 3AERI', 'REA_334': 'IONO =(delta)> I + NO2', 'REA_335': 'IONO + IONO => I2 + 2NO2', 'REA_336': 'I2 + NO3 => I + IONO2', 'REA_337': 'IO + NO2 => IONO2', 'REA_330': 'INO =(delta)> NO + I', 'REA_331': 'INO + INO => I2 + 2NO', 'REA_318': 'IO + HO2 => HOI + O2', 'REA_319': 'IO + NO => I + NO2', 'REA_316': 'HI + OH => I + H2O', 'REA_317': 'HOI + OH => IO + H2O', 'REA_314': 'I + I2O => IO + I2', 'REA_315': 'I2 + OH => HOI + I', 'REA_374': 'I2O4 + OIO => 3AERI', 'REA_313': 'I + IO => I2O', 'REA_310': 'I + HO2 => HI + O2', 'REA_311': 'I + I =(M)> I2', 'REA_376': 'OIO + NO => NO2 + IO', 'REA_377': 'IO => .5I2', 'REA_442': 'IO =(hv)> I + O3', 'REA_450': 'ICl =(hv)> I + Cl', 'REA_339': 'IONO2 + I = I2 + NO3', 'REA_345': 'I2 + Br => I + IBr', 'REA_344': 'I2 + Cl => I + ICl', 'REA_347': 'IO + Cl => I + ClO', 'REA_346': 'I2 +BrO => IO + IBr', 'REA_340': 'IONO2 =(M)> IO + NO2', 'REA_343': 'I + BrO => IO + Br', 'REA_342': 'I + Br2 => Br + IBr', 'REA_381': '"HIO3 => AERI ( ""aerosol"")"', 'REA_380': '"OIO => AERI (""aerosol"")"', 'REA_383': 'HOI =>0.5I2', 'REA_382': 'I2O5 => 2AERI', 'REA_349': 'IO + ClO => ICl + O2', 'REA_348': 'IO + Br => I + BrO', 'REA_443': 'OIO =(hv)> I + O2(M)','PD59': 'IONO2 => 0.5I2', 'PD58': 'HI => .5I2', 'PD57': 'IO => .5I2', 'PD56': 'OIO + NO => NO2 + IO', 'PD55': 'I2O4 + I2O4 => 4AERI', 'PD54': 'I2O4 + OIO => 3AERI', 'PD53': 'I2O3 + I2O4 => 4AERI', 'PD52': 'I2O3 + OIO => 3AERI', 'PD51': 'I2O4 =(M)> 2OIO', 'PD50': 'I2O4 + O3 => I2O5 + O2', 'PD70': 'I2O2 => IO + IO', 'PD63': 'HOI =>0.5I2', 'PD69': 'IONO2 =(hv)> I + NO3', 'PD67': 'OIO =(hv)> I + O2(M)', 'PD62': 'I2O5 => 2AERI', 'PD71': 'CH3IT =>  <CH3 +I>', 'PD49': 'I2O3 + O3 => I2O4 +O2', 'PD39': 'ICL + Br => BrCl + I', 'PD38': 'ICl + Cl => Cl2 + I', 'PD60': '"OIO => AERI (""aerosol"")"', 'PD01': 'I + O3 => IO + O2', 'PD02': 'I + HO2 => HI + O2', 'PD03': 'I + I =(M)> I2', 'PD04': 'I + IO => I2O', 'PD05': 'I + I2O => IO + I2', 'PD06': 'I2 + OH => HOI + I', 'PD07': 'HI + OH => I + H2O', 'PD08': 'HOI + OH => IO + H2O', 'PD09': 'IO + HO2 => HOI + O2', 'PD24': 'IONO2 + I = I2 + NO3', 'PD25': 'IONO2 =(M)> IO + NO2', 'PD22': 'I2 + NO3 => I + IONO2', 'PD23': 'IO + NO2 => IONO2', 'PD20': 'IONO =(delta)> I + NO2', 'PD21': 'IONO + IONO => I2 + 2NO2', 'PD68': 'IONO =(hv)> I + NO2', 'PD47': 'I2O2 =(M)> OIO + I', 'PD48': 'I2O2 + O3 => I2O3 +O2', 'PD28': 'I2 + Cl => I + ICl', 'PD75': 'INO => I + NO', 'PD44': 'IO + OIO =(M)> I2O3', 'PD45': 'OIO + OIO => I2O4', 'PD46': 'I2O2 =(M)> IO + IO', 'PD29': 'I2 + Br => I + IBr', 'PD40': 'ICl + OH => HOCl  +I', 'PD41': 'IBr + Br => I + Br2', 'PD42': 'IBr + OH => HOI + Br', 'PD43': 'IBr + OH => HOBr + I', 'PD26': 'I + Br2 => Br + IBr', 'PD64': 'I2 =(hv)> 2I', 'PD27': 'I + BrO => IO + Br', 'PD65': 'HOI =(hv)> I + OH', 'PD74': 'ICl =(hv)> I + Cl', 'PD33': 'IO + ClO => ICl + O2', 'PD61': '"HIO3 => AERI ( ""aerosol"")"', 'PD72': 'CH2I2 => <CH2 + I +I >', 'PD32': 'IO + Br => I + BrO', 'PD66': 'IO =(hv)> I + O3', 'PD73': 'IBr =(hv)> I + Br', 'PD13': 'OIO + OH => HIO3', 'PD12': 'IO + IO (O2) =>I2O2', 'PD11': 'IO + IO (O2) => I + OIO', 'PD10': 'IO + NO => I + NO2', 'PD17': 'INO =(delta)> NO + I', 'PD16': 'I + NO <=(M)> INO', 'PD15': 'HO + CH3IT => H2O + I (CH2I)', 'PD14': 'IO + CH3O2 =(M)> I + HO2 + HCHO', 'PD31': 'IO + Cl => I + ClO', 'PD30': 'I2 +BrO => IO + IBr', 'PD19': 'I + NO2 <=(M)> IONO', 'PD18': 'INO + INO => I2 + 2NO', 'PD35': 'IO + ClO => I + OClO', 'PD34': 'IO + ClO => I + Cl + O2', 'PD37': 'IO + BrO => Br +OIO', 'PD36': 'IO + BrO => Br + I + O2','IONO2_hv': 'REA_445', 'ICl_hv': 'REA_450', 'CH2I2_hv': 'REA_448', 'IBr_hv': 'REA_449', 'IONO_hv': 'REA_444', 'OIO_hv': 'REA_443', 'IO_hv': 'REA_442', 'I2O2_hv': 'REA_446', 'INO_hv': 'REA_451', 'CH3IT_hv': 'REA_447', 'HOI_hv': 'REA_441', 'I2_hv': 'REA_440', 'NO2_hv': 'REA_385', 'BrNO2_hv': 'REA_438', 'NO3_hv_II': 'REA_394', 'BrNO3_hv_II': 'REA_437', 'O3_hv': 'REA_384', 'NO3_hv': 'REA_393', 'CH3Br_hv': 'REA_439', 'HOBr_hv': 'REA_435', 'HONO_hv': 'REA_391', 'Br2_hv': 'REA_433', 'BrNO3_hv': 'REA_436', 'BrO_hv': 'REA_434','REA_247': 'O3 emission', 'REA_391': 'HONO =(hv)>', 'REA_150': 'O3 + PRPE =>', 'REA_249': 'CH3Br emission', 'REA_152': 'O3 + PMN =>', 'REA_437': 'BrNO3 =(hv)>BrO +NO2', 'REA_436': 'BrNO3 =(hv)> Br + NO3', 'REA_435': 'HOBr =(hv)>', 'REA_434': 'BrO =(hv)>', 'REA_433': 'Br2 =(hv)>', 'REA_439': 'CH3Br =(hv)>', 'REA_438': 'BrNO2 =(hv)>Br + NO2', 'REA_394': 'NO3 =(hv)> ONO + O2', 'REA_393': 'NO3 =(hv)> NO2 + O3', 'REA_1': 'O3 + NO =>', 'REA_2': 'O3 + OH =>', 'REA_3': 'O3 +HO2 =>', 'REA_4': 'O3 + NO2 =>', 'REA_251': 'Br2 emission', 'REA_197': 'O3 + IALD =>', 'REA_169': 'O3 + MACR =>', 'REA_168': 'O3 + MVK =>', 'REA_254': 'I2 emission', 'REA_253': 'CH2I2 emission', 'REA_252': 'CH3IT emission', 'REA_277': 'Br + O3 =>', 'REA_250': 'CH2Br2 emission', 'REA_279': 'Br + OH2 =>', 'REA_278': 'Br + OH =>', 'REA_167': 'O3 + ISOP =>', 'REA_385': 'NO2 =(hv)>', 'REA_384': 'O3 =(hv)> OH + OH', 'REA_248': 'HNO3 emission', 'TRA_80':'CH2ICl', 'TRA_81': 'CH2IBr', 'TRA_83': 'C3H5I','TRA_82': 'C3H7I'}


            output_species=tracer_library[input_species]
            return output_species
