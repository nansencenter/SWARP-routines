# Python script for reprojecting model product into OSI-SAF grid

from netCDF4 import Dataset
import sys,os
import glob
import numpy as np
import subprocess
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure as ms
from scipy.interpolate import griddata as grd

sys.path.append('../py_funs')
import mod_reading as Mrdg

############################################################################
def finish_map(m):
	# *m is a basemap
	m.drawcoastlines()
	m.fillcontinents(color='gray')
	# draw parallels and meridians.
	m.drawparallels(np.arange(60.,91.,10.),\
		labels=[True,False,True,True]) # labels = [left,right,top,bottom]
	m.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
		labels=[True,False,False,True])
	m.drawmapboundary() # fill_color='aqua')
	return
############################################################################

# DATE
date = raw_input('Insert date (YYYYMMDD):	')
fcday	= raw_input('Forecast day 0,1,2 (#):	')
typo	= 'waves'
subprocess.call(['./data/fetch_daily.sh',date,fcday,typo])

# DEFINING THE BASEMAP
m = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)

# READ TP4 DAILY
ncfil = ''.join( glob.glob('./data/TP4DAILY_start*_dump'+date+'.nc'))
print('TP4DAILY ice_only file = ' +ncfil+'\n')
slon		= 'longitude'
slat		= 'latitude'
sconc		= 'fice'
lon			= Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
lat			= Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
conc		= Mrdg.nc_get_var(ncfil,sconc,time_index=0)
X,Y      = m(lon[:,:],lat[:,:],inverse=False)
Z        = conc[:,:].data
mask     = conc[:,:].mask
Z[mask]  = np.NaN

# READ IN OSI-SAF FILE
ncfil2 = './data/ice_conc_nh_polstere-100_multi_'+date+'1200.nc'
print('OSISAF file = '+ncfil2+'\n')
clon     = 'lon'
clat     = 'lat'
cconc    = 'ice_conc'
lon2     = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
lat2     = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
conc		 = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
xc			 = Mrdg.nc_get_var(ncfil2,'xc')
yc			 = Mrdg.nc_get_var(ncfil2,'yc')
X2,Y2    = m(lon2[:,:],lat2[:,:],inverse=False)
Z2       = conc[:,:].data
mask2    = conc[:,:].mask
Z2[mask2] = np.NaN
ZO				= Z2/100

# REPROJECTION
X3 = X.reshape(X.size)
Y3 = Y.reshape(Y.size)
Z3 = Z.reshape(Z.size)
C = [X3,Y3]
C = np.array(C)
C = C.T

ZN = grd(C,Z3,(X2,Y2),method='nearest')
ZL = grd(C,Z3,(X2,Y2),method='linear')
ZC = grd(C,Z3,(X2,Y2),method='cubic')

# BINARY
BO					= np.copy(ZO)
BO[BO<.15]	= 0
BO[BO>=.15]	= 1
thenans			= np.isnan(BO)
BO[thenans]	= 0

BN					= np.copy(ZN)
BN[BN<.15]	= 0
BN[BN>=.15]	= 1
thenans			= np.isnan(BN)
BN[thenans]	= 0

BL					= np.copy(ZL)
BL[BL<.15]	= 0
BL[BL>=.15]	= 1
thenans			= np.isnan(BL)
BL[thenans]	= 0

# PLOT RESULTS
#plt.figure(0) 
#plt.imshow(ZO)
#plt.title('OSI')
#
#plt.figure(1)
#plt.imshow(ZN)
#plt.title('NRS')
#
#plt.figure(2)
#plt.imshow(ZL)
#plt.title('LIN')
#
#plt.figure(3)
#plt.imshow(ZC)
#plt.title('CUB')
#
#plt.show()


