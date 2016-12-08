import os,sys        # system stuff
import numpy as np # numerical stuff
from shapely.wkt import dumps, loads
from shapely.geometry import Polygon, Point 
from mpl_toolkits.basemap import Basemap
import datetime as dtm

#  takes .[a,b] files, calcs waves-in-ice region
import mod_reading as mr   # for reading files > $SWARP_ROUTINES repo: add $SWARP_ROUTINES/py_funs   to $PYTHONPATH
import WIM2d_f2py as Mwim  # WIMmodel + attenuation calc > $WIM2d repo:           add $WIM2D_PATH/fortran/bin   to $PYTHONPATH (infterfaces .F code)
import geometry_planar as GP

 
# reference time
tref  = dtm.datetime(1970,1,1,0,0,0)
# TODO find closest topaz file 00.00 + dt=3h/6h
 

################################################################
# Use projection from Dr Fab's polygons
# projection - EPSG:3413
lonc     = -45.      # central meridian (x=0)
latc     = 90.       # central lat
lat_ts   = 70.       # standard parallel (true scale lat)
ellps    = 'WGS84'   # shape of earth
width    = 1000.e3   # width of domain in m
height   = 1000.e3   # height of domain in m

# basemap corresponding to the projection
bmap     = Basemap(projection='stere',ellps=ellps,\
                     lon_0=lonc,lat_0=latc,lat_ts=lat_ts,\
                     width=width,height=height)
################################################################

# parameters in WIM mode
#  TODO put in as a range of values ?
young     = 5.49e9
visc_rp   = 13

# TODO get T from Dr. Fabs wave periods
T = 20
################################################################
def atten_nondim(h,T,young,visc_rp):
   # attenuation per floe
   # - for scalars (ie single grid cell)

   # h        = inputs(1)
   # om       = inputs(2)
   # young    = inputs(3)
   # visc_rp  = inputs(4)

   om       = 2*np.pi/T
   inputs   = np.array([h,om,young,visc_rp])

   outs  = Mwim.atten_youngs(inputs) # calls interface to fortran code
   # outputs  = (/damping,kice,kwtr,int_adm,
   #              ac,modT,argR,argT/)

   k_visc   = outs[0]
   kice     = outs[1]
   kwtr     = outs[2]
   alp_scat = outs[4]

   return k_visc,kice,kwtr,alp_scat
################################################################


################################################################
def atten_dim(k_visc,alp_scat,c,Dmax):
# TODO change Dmax to Dmin, should be able to calculate from Dmax
# from floe scaling in topaz code ...
   alp   = 2*c*k_visc+alp_scat*c/Dmax
  
   return alp
################################################################


# waves-ice file
#afil  = os.getenv('HOME')+'/test_lucia/TP4archv_wav.2015_351_120000.a'
afil = '/home/nersc/timill/test_lucia/TP4archv_wav.2015_351_120000.a'
hi    = mr.HYCOM_binary_info(afil)

# TODO get time
# hitime = hi.get_var('time').values
fice  = hi.get_var('fice').values
hice  = hi.get_var('hice').values
Dmax  = hi.get_var('dfloe').values
lon,lat  = hi.get_lonlat()

nx,ny = hi.Nx,hi.Ny
# remap to Dr Fabs projection
hix,hiy = bmap(lon,lat)

# mask for polygon and for larger square
pmask = np.zeros_like(fice)
smask = np.zeros_like(fice)

# loop over grid:
# - get c,h,Dmax, convert to local atten (alp_dim) 
# - also need T from Dr Fab's netcdf file (TODO)(also in dropbox)

# TODO  get corners from Dr Fab's polygon (dropbox)
# TODO  get a) max/min Lon/Lat get hi poins inside, or
#           b) hi points inside Dr Fab's polygon
#        
# POLYGON filename ex:
# ODL_wavesinice_poly20151012T100000Z.txt
# POLYGON file format ex:
# Polygon | Lon | Lat | Type
# 0    -135.505860677717465    66.696926716725855    2
# 0    -135.999887551280068    66.789449259182035    2
polyfile = 'ODL_wavesinice_poly20151012T100000Z.txt'
polydir  = '/home/nersc/bergh/python_script/WIM/ODL/ODL_polys/'
file = open(polydir+polyfile, "r")
head = file.readline()
data = np.loadtxt(file,skiprows=0)
file.close()
pnr  = data[:,0]
plon = data[:,1]
plat = data[:,2]
ptype= data[:,3]
plonmax=max(plon)
plonmin=min(plon)
platmax=max(plat)
platmin=min(plat)
# get polygon lon,lat into proj map x,y 
px,py = bmap(plon,plat)
pxpy=[(px[i],py[i]) for i in range(len(px))]
# pp=Polygon(pxpy)

# TODO limit i,j positions within polygon
# mask loop before/after when Dmax(i,j) = 0, => "RuntimeWarning: divide by zero" 
inside=GP.maskgrid_outside_polygon(hix,hiy,pxpy)
mask=np.logical_and(inside,Dmax>0.) # only want atten inside ice
alp_dim=np.zeros_like(Dmax)

for i in range(nx):
   for j in range(ny):
      if mask[i,j]:
         c  = fice[i,j]
         h  = hice[i,j]
         dmax  = Dmax[i,j]
         k_visc,kice,kwtr,alp_scat  = atten_nondim(h,T,young,visc_rp)
         alp_dim[i,j]               = atten_dim(k_visc,alp_scat,c,dmax)





