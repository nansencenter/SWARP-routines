import os,sys
import mod_reading as MR
import MIZchar as MC
from netCDF4 import Dataset as ncopen
import numpy as np
from datetime import datetime
import pyproj
from skimage import measure as msr

TEST_PLOT   = 1
STEP        = 4 # don't take all points (eg STEP=4 gives every 4th)

# stereographic projection
map   = pyproj.Proj(proj='stere',lon_0=-45,lat_0=90,lat_ts=70,ellps='WGS84')

# file info
fcdir = '/work/users/timill/RealTime_Models/results/TP4a0.12/ice_only'
ncfil = fcdir+'/20170126/final_output/SWARPiceonly_forecast_start20170126T000000Z.nc'
nci   = MR.nc_getinfo(ncfil)
print('\nReading\n'+ncfil)

# global attributes
glob_atts   = nci.get_global_atts()
del(glob_atts['history'])
#print(glob_atts)
#print('\n')

# lon-lat
LON,LAT  = nci.get_lonlat()
X,Y      = map(LON,LAT)

# time info
time_index  = 0
tval        = nci.get_time()[time_index]
time_atts   = nci.get_var_atts('time')
#print(tval,time_atts)
#print('\n')

# variable info
vbl      = 'icec'
V_atts   = nci.get_var_atts(vbl)
fillval  = None
if '_FillValue' in V_atts:
   fillval  = V_atts['_FillValue']
   del(V_atts['_FillValue'])


# output file
dto      = nci.datetimes[time_index]
outdir   = 'out/contours/STEP'+str(STEP)
if not os.path.exists(outdir):
   os.mkdir(outdir)
ncfil2   = outdir+'/'+nci.basename+'_'+vbl+dto.strftime('_contours_%Y%m%dT%H%M%SZ.nc')
print('\nWriting\n'+ncfil2+'\n')


if TEST_PLOT:
   if 1:
      # zoom
      DC       = 'k'#date color
      #reg      = 'gre'# greenland sea
      reg      = 'FR1'# greenland sea
      figname  = ncfil2.replace('.nc','_'+reg+'.png')
   else:
      DC       = 'r'
      reg      = 'Arctic'
      figname  = ncfil2.replace('.nc','.png')
   print('Saving figure: \n'+figname+'\n')
   var_opts    = MR.make_plot_options(vbl,ice_mask=True)
   pobj,bmap   = nci.plot_var(var_opts,show=False,HYCOMreg=reg,\
                  date_label=2,date_color=DC,\
                  clim=[0,1],add_cbar=True,clabel='sea ice concentration')
   fig   = pobj.fig
   ax    = pobj.ax


# ---------------------------------------------------------------------------
# get the contours

print('Getting contours...')
V  = nci.get_var(vbl,time_index=time_index)
if 0:
   # get vbl range
   Vmin     = V.min()
   Vmax     = V.max()
   dV       = (Vmax-Vmin)/20.
   Vlevels  = np.arange(Vmin,Vmax+dV,dV)
   print(Vmin,Vmax)
else:
   # manual vbl range
   Vmin     = .15
   Vmax     = 1.
   dV       = .1
   Vlevels  = np.arange(Vmin,Vmax+dV,dV)

Nlevels  = len(Vlevels)
Ntot     = 0
Vconts   = []
for Vlev in Vlevels:
   # create a binary corresponding to the level
   B                 = np.zeros(V.shape)
   good              = np.logical_not(V.values.mask)
   data              = V.values.data[good]
   Bgood             = np.zeros(len(data))
   Bgood[data>=Vlev] = 1.
   B[good]           = Bgood
   #print(len(Bgood[Bgood>.5]))
   #print(data)

   # find the contours
   vcont = sorted(msr.find_contours(B,.5), key=len)
   Vcont = []
   for cont in vcont:
      ij_in_rows  = np.array(cont)
      x,y         = MC.ij2xy(ij_in_rows,X,Y).transpose()
      if len(x)<10:
         continue
      elif len(x)<20:
         lon,lat  = map(x,y,inverse=True)
      else:
         # reduce points according to STEP
         lon,lat  = map(x[::STEP],y[::STEP],inverse=True)

      length   = len(lon)
      Vcont.append({'lon':lon,'lat':lat,'contour_value':Vlev,'length':length})
      Ntot += length

      if TEST_PLOT:
         bmap.plot(lon,lat,latlon=True,ax=ax,linewidth=2)

   #*Vcont is list of dictionaries corresponding to a single level of the variable
   #*Vconts is list of all the "Vcont"s for all the levels of the variable
   Vconts.append(Vcont)
# ---------------------------------------------------------------------------


if TEST_PLOT:
   print('Saving figure...')
   fig.savefig(figname)
   MR.plt.show(fig)
   ax.cla()
   MR.plt.close(fig)
   
cno_vec  = np.zeros(Ntot,dtype='int')  # contour numbers
c_vec    = np.zeros(Ntot)              # values of the contour
lon_vec  = np.zeros(Ntot)              # lon position
lat_vec  = np.zeros(Ntot)              # lat position

position       = 0
contour_number = 0
for Vcont in Vconts:
   for D in Vcont:
      N     = D['length']
      C     = D['contour_value']
      lonc  = D['lon']
      latc  = D['lat']
      cno_vec[position:position+N]  = contour_number
      c_vec  [position:position+N]  = C
      lon_vec[position:position+N]  = lonc
      lat_vec[position:position+N]  = latc

      #increment position
      position         += N
      contour_number   += 1


print('Writing output netcdf file...')
nc = ncopen(ncfil2,'w', format='NETCDF4') #'w' stands for write

"""
The above line creates a netCDF file called "sample.nc" in the current folder. nc is a netCDF Dataset object that provides methods for storing data to the
file.
"""

if 0:
   """
   nc also doubles as the root group. A netCDF group is basically a directory or folder within the netCDF dataset.
   This allows you to organize data as you would in a unix file system. Let's create a group for the heck of it:
   """
   tempgrp = nc.createGroup('Temp_data')
else:
   # write directly to the root group
   tempgrp = nc

"""
Specifying dimensions

The next step is to specify the dimensions of the data. If you plan to save a multidimensional array of data, each dimension of that array needs to be given a name and a length:
"""

tempgrp.createDimension('point_index', Ntot)
tempgrp.createDimension('time', 1)

"""
Building variables

This step essentially pre-allocates NetCDF variables for data storage. NetCDF variables are very similar to numpy arrays. To construct them, you use the createVariable method:
"""
#idx   = tempgrp.createVariable('point_index'    ,'i4' ,'point_index',zlib=True)
Cno   = tempgrp.createVariable('contour_number' ,'i4' ,'point_index',zlib=True)
Lon   = tempgrp.createVariable('longitude'      ,'f4' ,'point_index',\
            zlib=True,least_significant_digit=8)
Lat   = tempgrp.createVariable('latitude'       ,'f4' ,'point_index',\
            zlib=True,least_significant_digit=8)
Vnc   = tempgrp.createVariable(vbl              ,'f4' ,'point_index',\
            fill_value=fillval,zlib=True,least_significant_digit=4)
Time  = tempgrp.createVariable('time'           ,'i4' ,'time',zlib=True)

"""
Passing data into variables

Here, you simply pass your data into the variables you just created:
"""
Cno[:]   = cno_vec #The "[:]" at the end of the variable instance is necessary
Lon[:]   = lon_vec
Lat[:]   = lat_vec
Vnc[:]   = c_vec
Time[0]  = tval

"""
Adding attributes

NetCDF attributes can be used to provide additional information about the dataset (i.e. metadata). You can add attributes to variables, groups and the dataset itself. This is optional but considered good practice:
"""

#Add global attributes
for att in glob_atts:
   setattr(tempgrp,att,glob_atts[att])


#Add time attributes
for att in time_atts:
   setattr(Time,att,time_atts[att])

# lon,lat attributes
lon_atts = nci.get_var_atts('lon')
for att in lon_atts:
   setattr(Lon,att,lon_atts[att])
lat_atts = nci.get_var_atts('lat')
for att in lat_atts:
   setattr(Lat,att,lat_atts[att])

# main variable
for att in V_atts:
   setattr(Vnc,att,V_atts[att])

Cno.long_name  = 'contour_number'
Cno.units      = ''

tempgrp.comment   = 'Positions correspond to isolines of '\
      +vbl+' ('+V_atts['standard_name']+') for '+dto.strftime('%Y-%m-%d %H:%M:%S (')\
                     +'record number '+str(time_index+1)+') in ' +nci.basename+'.nc'
nc.close()
