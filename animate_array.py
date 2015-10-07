from matplotlib import animation
import matplotlib.pyplot as plt
import time
import numpy as np
import gc
import datetime

debug=True

def main( spec='NO' , pcent=False, fixcb=None, limit_by_dates=False, \
        extend='neither', debug=False ):
    """ Extract data array from given location and make Animation"""

    # Get data in array ( time, lat, lon ) and dates (datetime.datetime )
    arr, dates = get_data_dates( spec=spec, limit_by_dates=limit_by_dates,
        debug=debug)
    if debug:
        print [ (i[:5], i.shape ) for i in arr, dates]
    
    # Get titles and other varibles for run (e.g. res )
    lat, lon, units, fname, res, title, scale = get_run_info( spec=spec )
    arr = arr*scale

    # Setup figure and axis
    fig, ax = setup_figure_and_axis( )

    # manual limit colormap ( colorbar )
#    fixcb = np.array([ 0, 100])
#    extend ='both'
    
    # Setup first frame for re-use ( inc. basemap ) and get key variables 
    cmap, specplt, clevs, cnorm, minc, maxc, m = setup_plot2animate( arr,\
            fig=fig, ax=ax, lat=lat, lon=lon, units=units, res=res, \
            fixcb=fixcb, debug=debug )

    # setup figure ascetics
    setup_figure_ascetics( dates, cmap=cmap, cnorm=cnorm, maxc=maxc, \
        minc=minc, units=units, fig=fig, title=title, extend=extend, \
        arr=arr, fixcb=fixcb, debug=debug)
    
    # animate the array and save as 
    animate_array( arr, dates, specplt, clevs=clevs, cnorm=cnorm, \
        cmap=cmap, debug=debug, fig=fig, m=m, lon=lon, lat=lat, \
        spec=spec, fname=fname )

def get_data_dates( spec='O3', dates_variable='time', \
            fill_invalid_with_mean=True, limit_by_dates=False, 
            sdate=datetime.datetime(2005, 01, 01), 
            edate=datetime.datetime(2005, 01, 07), debug=False ):
    """ Extracts dates and data from a given location """

    import numpy as np
    from pandas import DataFrame
    from netCDF4 import Dataset
    import datetime
    from AC_tools.funcs_vars import get_dir

    # Set file to use
    wd = get_dir('npwd')
#    wd = '/work/home/ts551/temp/'
    f= 'pf_iGEOSChem_1.7_v10_G5_EU_run.0.25x0.3125.2012.1week.sucess_3D.nc'

#    f= 'pf_iGEOSChem_1.7_v10_G5_EU_run_3D.nc'
    print wd + f
    # Extract data
    with Dataset( wd+f , 'r' ) as rootgrp:

        if debug:
            print [i for i in rootgrp.variables ]

        # return data as an array
#        arr = np.ma.array( rootgrp.variables[ spec ]   )
        arr = np.array( rootgrp.variables[ spec ]   )
        print rootgrp.variables[ spec ] 
        print np.array( rootgrp.variables[ spec ]  )

        # get dates
        dates = np.ma.array( rootgrp.variables[ dates_variable ]   )

    if debug:
        print [ ( type( i ), i.shape ) for i in arr, dates ]
        print [ (i.min(), i.max(), i.mean()) for i in arr, dates ]

    # mask with invalid values ( fill values in array) 
    arr = np.ma.masked_invalid( arr  )
    if debug:
        print [ (i.min(), i.max(), i.mean()) for i in [arr] ]
        
    # Make sure dates are as datetime and in a numpy array
    dates = [datetime.datetime.fromtimestamp(i) for i in dates ]    
    dates = np.array(dates )

    if debug:
        print [ i.shape for i in arr, dates ]
    
    print edate, sdate, dates[0]
    print [ type(i) for i in edate, sdate, dates[0] ]
    
    print dates > sdate 
    print dates < edate 
    
#    print dates[ np.where( (dates >= sdate)  ) ] 
    
    # Limit to given dates ( e.g. 1st month 2005)
    if limit_by_dates:
        dates = dates[ np.where( dates >= sdate ) ] 
        print [i.shape for i in arr, dates ]
        dates = dates[ np.where( dates < edate ) ] 
        print [i.shape for i in arr, dates ]
        # Kludge, remove 1st dimension added by method.  <= improve this. 
        arr = arr[ np.where( dates >= sdate ), ... ][0,...] 
        print [i.shape for i in arr, dates ]
        arr = arr[ np.where( dates < edate ), ... ][0,...] 
        print [i.shape for i in arr, dates ]
        
        
    # The corner point is a NaN for 0.25 output.  <= improve this. 
    # Set this ( not visible on the plot window ) to mean to allow for save
    if fill_invalid_with_mean:
        np.ma.set_fill_value(arr, arr.mean())
        arr = arr.filled()

    return arr, dates

def get_run_info( spec='O3', res='0.25x0.3125', region='EU', fname='', \
            scale=1, pcent=False ):
    """ get descriptive variables for run period ( e.g. res ) """
    from AC_tools.funcs4core import get_latlonalt4res
    from AC_tools.funcs_vars import tra_unit, latex_spec_name
    
    # Set variables (e.g. res) or automation of variable setting here
#    res = get_run_descriptors()
#    res='0.5x0.666'
    
    # Get lat and lon of GC grid for given resolution
    lon, lat, NIU = get_latlonalt4res( res=res )

    # Set units based on species name + get scaling 
    if pcent:
        units ='%'        
    else:
        units, scale = tra_unit( spec, scale=True)

    # Set filename from run detail
    fname += '{}_{}_{}_{}_{}.mp4'.format( region, res, 
        spec, units,  time.strftime("%y_%m_%d_%H_%M")   )

    # setup plot title
    title = 'Surface {} / {}'.format( latex_spec_name( spec ), units )
    
    return lat, lon, units, fname, res, title, scale


def setup_figure_and_axis( ):
    """ Initialise figure and axis """

    # Setup plot
    fig, ax = plt.subplots(figsize=[16,9])
    ax.set_aspect('equal')
    ax.autoscale_view(False)
    
    return fig, ax

def setup_figure_ascetics(  dates, f_size=10, title=None, cmap=None, \
            minc=None, units=None, cnorm=None, fig=None, maxc=None,\
            format='%.0f', extend='neither', arr=None, fixcb=None, debug=False):
    """ Add colorbar, logos and titles to figure """

    # if on Univeristy of York/NCAS servers, add logos
    if platform.platform() == 'Linux-3.0.101-0.47.52-default-x86_64-with-SuSE-11-x86_64':
        from funcs4plotting_special import add_logos_NCAS_york_bottom, mk_cb

        # add title and logos for NCAS/NERC
        fig = add_logos_NCAS_york_bottom( fig)
        fig.suptitle( title, fontsize=f_size*2, x=.55 , y=.95  )

    # Set min and max of colorbar if not given
    if isinstance( fixcb, type(None) ):
        fixcb = ( minc, maxc )

    # Add colorbar  
    mk_cb(fig, units=units, left=0.925, bottom=0.2, width=0.015, height=0.6, \
        cmap=cmap, vmin=fixcb[0], vmax=fixcb[1], f_size=f_size, format=format, \
        norm=cnorm, extend=extend,  ) 

    # Add text showing model time on figure 
    plt.figtext(x=0.5,y=0.875, s=dates[0], fontsize=f_size)

def setup_plot2animate( arr, fig=None, ax=None, lat=None, lon=None, \
        units=None, contour=False, res=None, f_size=20, min_change=0.5,  \
        period=1, adjust_window=3, alpha=0.85, npoints=50, \
        everyother=1, interval=1, resolution='f', drawcountries=True, \
        cnorm = None, clevs=None, fixcb=None, debug=False ):

    from funcs4plotting_special import get_basemap, get_colormap

    # --- settings
    plt.ioff() # turn off interactive plotting
    global specplt

    # Get colormap for fixed length or min/max or array
    if isinstance( fixcb, type(None) ):
        fixcb = arr.copy() 
    cmap = get_colormap( fixcb )
    
    print 3, 'detail on output: ', [ [ np.ma.min(i), np.ma.max(i), \
        np.ma.mean(i), type(i),i.shape ] for i in [arr] ]

    # Setup basemap
    m = get_basemap( lat=lat, lon=lon, resolution=resolution, res=res,\
        everyother=everyother, interval=interval, f_size=f_size, \
        drawcountries=drawcountries )

    # adjust window size
    m.ax = ax
    plt.xlim( lon[0+adjust_window], lon[-1-adjust_window] )
    plt.ylim(lat[0+adjust_window], lat[-1-adjust_window])

    # --- Setup grid and Pcolor/contour
    x, y = np.meshgrid(lon,lat)

    # Make sure array is a masked numpy array
    arr  = np.ma.array( arr )

    # Setup colour or pcolor plot
    if contour:
        minc=arr.min()
        maxc=arr.max()

        clevs = np.linspace(minc, maxc+min_change, npoints )
        cnorm = mpl.colors.Normalize(vmin=minc, vmax=maxc)
        specplt = m.contourf(x,y,arr[0,:,:],clevs,cmap=cmap, 
            latlon=True, norm=cnorm,antialiased=True)

        fname = 'contour_'+fname

    else:
        minc=arr.min()
        maxc=arr.max()
        specplt = m.pcolor( lon, lat, arr[0,:,:], cmap=cmap)
        #, antialiased=True )

    if debug:
        print [len(i) for i in lat, lon, x, y ]
        print arr.shape, type( arr) 
        print arr.shape, type( arr) 

    return cmap, specplt, clevs, cnorm, minc, maxc, m 

def animate_array( arr, dates, specplt, spec='O3', min_change=0.5, \
        period=1, adjust_window=3, alpha=0.85, npoints=50, wd='./', \
        clevs=None, cnorm=None, cmap=None, contour=False, \
        fig=None, m=None, lon=None, lat=None, fname=None, debug=False ):
    """ Animates array with specplt as first frame"""

    # clean memory
    gc.collect()

    # Function to loop frame and save to animation
    def plotgc(i):
        global specplt
        print i, 'for {} @ model time of {} @  real time of {}'.format( spec, \
             dates[i], time.strftime('%Y %m %d %H:%M:%S') )
        if contour:
            # remove previous plot data, but retain basemap....
            for c in specplt.collections:
                c.remove()
            # Fill basemap 
            specplt = m.contourf(x,y,arr[i,:,:], clevs, cmap=cmap, 
                latlon=True, norm=cnorm, alpha=alpha, antialiased=True)
        else:
            # remove previous plot data, but retain basemap....
            specplt.remove()
            # Fill basemap/ plot
            specplt = m.pcolor( lon, lat, arr[i,:,:], cmap=cmap)

        # clean memory
        gc.collect()

        # update date
        fig.texts[-1].set_text( dates[i] )
        return specplt

    # Animate
    ani=animation.FuncAnimation(fig, plotgc, \
            frames=np.arange(len(dates)-2),blit=True)

    # Save
    ani.save( wd + fname, 'ffmpeg', fps=24, extra_args=['-vcodec', 
            'libx264','-pix_fmt', 'yuv420p'])
    print 'Video saved & Closed as/at: ', fname

if __name__ == "__main__":
    main( debug=debug )
