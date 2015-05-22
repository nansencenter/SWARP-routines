# Purpose of the script is:
# 1)Reproject daily binary file (TODO even the .nc maybe) into OSI-SAF grid
# 2)Calculate the Area of disagreement between the two

from netCDF4 import Dataset
import sys,os
import glob
import numpy as np
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure as ms

pdir  = '../py_funs'
if pdir not in sys.path:
   sys.path.append(pdir)

import mod_reading   as Mrdg
import fns_plotting  as Mplt

############################################################################
# BM MAP
# setup stereographic basemap
# - for plotting and also for limiting search area
# lat_ts is latitude of true scale.
# (lon_0,lat_0) is central point -> (x,y)=(0,0)
rad   = 18.          # approx radius of image (degrees)
xmax  = rad*111.e3   # half width of image [m]
ymax  = rad*111.e3   # half height of image [m]
cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
#
lat_ts   = 74. # deg N
lon_0    = 50. # deg E
lat_0    = 74. # deg N
#
bm = Basemap(width=2*xmax,height=2*ymax,\
             resolution=cres,projection='stere',\
						              lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)

############################################################################
def finish_map(bm):
	# finish_map(bm)
	# *bm is a basemap
	bm.drawcoastlines()
	bm.fillcontinents(color='gray')
	
	# draw parallels and meridians.
	bm.drawparallels(np.arange(60.,91.,10.),\
		labels=[True,False,True,True]) # labels = [left,right,top,bottom]
	bm.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
		labels=[True,False,False,True])
	bm.drawmapboundary() # fill_color='aqua')
	
	return
############################################################################

##################################################################
# READ IN ARRAYS FROM TP4DAILY BINARIES (.a)
# lon/lat from grid binary
tdir  = 'data'
gfil  = tdir+'/regional.grid.a'
plon  = Mrdg.get_array_from_HYCOM_binary(gfil,1)
nx,ny = plon.shape
plat  = Mrdg.get_array_from_HYCOM_binary(gfil,2,dims=(nx,ny))

if 0:
	plt.imshow(plon)
	plt.colorbar()
	plt.show()
	sys.exit('L22')
 
afil	= ''.join(glob.glob(tdir+'/TP4*.a'))
vlst  = {'ficem':1,'hicem':2,'dfloe':3,'swh':4,'mwp':5}
vbl   = 'ficem'
recno = vlst[vbl]
V     = Mrdg.get_array_from_HYCOM_binary(afil,recno,dims=(nx,ny))
##################################################################
date		= '20150520'
# READ TP4 DAILY
ncfil		= ''.join( glob.glob('./data/TP4DAILY_start*_dump'+date+'.nc'))
slon		= 'longitude'
slat		= 'latitude'
sconc		= 'fice'
lon			= Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
lat			= Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
conc		= Mrdg.nc_get_var(ncfil,sconc,time_index=0)
tlon		=
Z				= conc[:,:].data
mask		= conc[:,:].mask
Z[mask]	= np.NaN
##################################################################
# READ IN OSI-SAF FILE
ncfil2	= './data/ice_conc_nh_polstere-100_multi_'+date+'1200.nc'
clon		= 'lon'
clat		= 'lat'
cconc		= 'ice_conc'
lon2		= Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
lat2		= Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
conc		= Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
olon		= lon2[:,:]
olat		= lat2[:,:]
Z2			= conc[:,:].data
mask2			= conc[:,:].mask
Z2[mask2]	= np.NaN
##################################################################

# Need a rank 1 array for x and y of original file



# PRINT OUT NAME OF VARIABLES - TO BE USED ONLY UNTIL WE WORK IN IPYTHON
print('For model: plon,plat,V - For OSI: olon,olat,Z2')
