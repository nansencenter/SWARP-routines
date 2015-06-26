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
############################################################################
def binary_mod(data,thresh):
   ndata                = np.copy(data)
   ndata[ndata<thresh]	= 0
   ndata[ndata>=thresh]	= 1
   #thenans		= np.isnan(ndata)
   #ndata[thenans]      = 0
   return(ndata)
############################################################################
############################################################################
def binary_diff(data1,data2,thresh):
   mdata                   = np.copy(data1)
   mdata[mdata<thresh]	= 0
   mdata[mdata>=thresh]	= 1
   odata			= np.copy(data2)
   odata[odata<thresh]	= 0
   odata[odata>=thresh]	= 1	
   ddata			= odata - mdata
   #thenans		= np.isnan(ddata)
   #ddata[thenans]		= 0
   return(ddata)
###########################################################################
###########################################################################
def figure_save(X,Y,Z,name,m):
   f = plt.figure()
   m.pcolor(X,Y,Z,cmap='YlGnBu',vmin=-1,vmax=1)
   finish_map(m)
   fname = './outputs/aod/'+name+'.png'
   plt.colorbar()
   plt.title(name)
   plt.savefig(fname,format='png',dpi=1000)
   plt.close()
   f.clf()
###########################################################################
###########################################################################
def get_stats(data,name):
   stat0	= (data == 0).sum()
   stat1	= (data == 1).sum()
   stat2	= (data == -1).sum()
   stat	= data.size
   pone	= (stat1 / float(stat)) * 100
   pmone	= (stat2 / float(stat)) * 100
   print "Number of elements:	",stat
   print "Hit:			",stat0
   print "Underprediction(+1):	",stat1,pone
   print "Overprediction(-1):	",stat2,pmone
   print " " 
   f = open('./outputs/aod/'+name+'.txt','w')
   f.write("Number of elements:	%s \n" % stat)
   f.write("Hit:			%s \n" % stat0)
   f.write("Underprediction(+1):	%s - %s \n" % (stat1,pone))
   f.write("Overprediction(-1):	%s - %s \n" % (stat2,pmone))
   f.write(" ")
   f.close
   return()
###########################################################################
###########################################################################
def binary_cont(X,Y,D,O,M):
   ND = np.copy(D)
   for m,el in enumerate(X2):
      for n,el in enumerate(X2[m]):
         if D[m][n] == 1:
            around = ((m,n+1),(m,n-1),(m+1,n),(m-1,n),(m+1,n+1),(m-1,n+1),(m+1,n-1),(m-1,n-1))
            for hor,ver in around:
               if O[hor][ver] == M[hor][ver] == 1:
                  ND[hor][ver] = 2
               elif O[hor][ver] == M[hor][ver] == 0:
                  ND[hor][ver] = -2
               elif np.isnan(O[hor][ver]) or np.isnan(M[hor][ver]):
                  ND[hor][ver] = 3
         elif D[m][n] == -1:
            around = ((m,n+1),(m,n-1),(m+1,n),(m-1,n),(m+1,n+1),(m-1,n+1),(m+1,n-1),(m-1,n-1))
            for hor,ver in around:
               if O[hor][ver] == M[hor][ver] == 1:
                  ND[hor][ver] = -2
               elif O[hor][ver] == M[hor][ver] == 0:
                  ND[hor][ver] = 2
               elif np.isnan(O[hor][ver]) or np.isnan(M[hor][ver]):
                  ND[hor][ver] = 3
   return(ND)
###########################################################################
##########################################################################

# DEFINING THE BASEMAP
m = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)

# CHOOSE
FIGURE = raw_input('1 for figure, [Enter] for not:	')

# READ TP4 DAILY
ncfil = ''.join( glob.glob('./data/TP4DAILY*.nc'))
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

# DATE
dadate		= ncfil[-11:-3]

# READ IN OSI-SAF FILE
ncfil2 = ''.join( glob.glob('./data/ice_conc_nh_polstere-100_multi_*.nc'))
print('OSISAF file = '+ncfil2+'\n')
clon     = 'lon'
clat     = 'lat'
cconc    = 'ice_conc'
lon2     = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
lat2     = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
conc		 = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
xc			 = Mrdg.nc_get_var(ncfil2,'xc')
yc			 = Mrdg.nc_get_var(ncfil2,'yc')
X2,Y2		 = m(lon2[:,:],lat2[:,:],inverse=False)
XO			 = np.copy(X2)
YO			 = np.copy(Y2)
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
BO	= binary_mod(ZO,.15)
BN	= binary_mod(ZN,.15)
BL	= binary_mod(ZL,.15)

# DIFFERENCES
DN = binary_diff(ZN,ZO,.15)
DL = binary_diff(ZL,ZO,.15)

print "DN stats"
get_stats(DN,'DNstats'+dadate)
print "DL stats"
get_stats(DL,'DLstats'+dadate)

NDN = binary_cont(XO,YO,DN,BO,BN)

if FIGURE:
	figure_save(X2,Y2,DN,'DN'+dadate,hqm)
	figure_save(X2,Y2,DL,'DL'+dadate,hqm)

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


