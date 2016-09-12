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

# GET DR. FABS POLYGON AREA AND DATA
# polygon file and datetime
polyfile = 'ODL_wavesinice_poly20151012T100000Z.txt'
polydir  = '/home/nersc/bergh/python_script/WIM/ODL/ODL_polys/'
# get thsi from file or file-header ?
year     = 2015
month    = 10
date     = 12
hh       = 10 #find closest arond 10
pdate = dtm.datetime(year,month,date,hh,00,00)
print 'ODL Date Time: ',pdate.__str__()

################################################################
# Use projection from Dr Fab's polygons
# projection - EPSG:3413
lonc     = -45.      # central meridian (x=0)
latc     = 90.       # central lat
lat_ts   = 70.       # standard parallel (true scale lat)
ellps    = 'WGS84'   # shape of earth
width    = 1000.e3   # width of domain in m
height   = 1000.e3   # height of domain in m

#) basemap corresponding to the projection
bmap     = Basemap(projection='stere',ellps=ellps,\
                     lon_0=lonc,lat_0=latc,lat_ts=lat_ts,\
                     width=width,height=height)
################################################################

# TODO Do we need this ?
################################################################
# wave data from ERA_Interim
eraidir='/work/shared/nersc/ERA-I_0.125deg/'
eraifile='SWH_'+str(year)+'.nc'
ncfil=eraidir+eraifile
nci=mr.nc_getinfo(ncfil)
################################################################

nci_nearestdate = nearestDate(nci.datetimes,pdate)
rec=nci.datetimes.index(nci_nearestdate)
print 'Nearest Date Time: ',nci.datetimes[rec].__str__(),' record: ',rec


# TODO FIND CLOSEST TIME AND GET FILE ?
# TODO LOOP OVER SOME PERIOD AROUND EVENT
# GET TOPAZ DATA 
adir = '/home/nersc/timill/test_lucia/'
afil = 'TP4archv_wav.2015_351_120000.a'
hi    = mr.HYCOM_binary_info(adir+afil)

hidt  = hi.datetime
fice  = hi.get_var('fice').values
hice  = hi.get_var('hice').values
Dmax  = hi.get_var('dfloe').values
swh   = hi.get_var('swh').values
mwp   = hi.get_var('mwp').values
mwd   = hi.get_var('mwd').values
lon,lat  = hi.get_lonlat()
#nx,ny = hi.Nx,hi.Ny
# remap to bmap projection
# hix,hiy = bmap(lon,lat)

# mask for polygon and for larger square
pmask = np.zeros_like(fice)
smask = np.zeros_like(fice)# waves-ice file

# Read in polygon TODO ? read in dr fabs x,y directly
file  = open(polydir+polyfile, "r")
head  = file.readline()
data  = np.loadtxt(file,skiprows=0)
file.close()
pnr   = data[:,0]
plon  = data[:,1]
plat  = data[:,2]
ptype = data[:,3]
# get polygon lon,lat into proj map x,y 
px,py     = bmap(plon,plat)
poly_pxpy = Polygon(zip(px,py))
# Get centroid and bounding box of polygon
pxc,pyc = poly_pxpy.centroid.xy
xll,yll,xur,yur = poly_pxpy.bounds
print(px[0],py[0])
#print(xll,yll)
#print(xur,yur)
plonc,platc=bmap(pxy[0],pyc[0],inverse=True)

################################################################
# NEW projection centered around polygon area
dist=100  # hor and ver size of area in km
# projection - EPSG:3413
lonc     = plonc      # central meridian (x=0)
latc     = platc       # central lat
lat_ts   = 70.       # standard parallel (true scale lat)
ellps    = 'WGS84'   # shape of earth
width    = dist*1.e3   # width of domain in m
height   = dist*1.e3   # height of domain in m

#) basemap corresponding to the projection
plotmap     = Basemap(projection='stere',ellps=ellps,\
                     lon_0=lonc,lat_0=latc,lat_ts=lat_ts,\
                     width=width,height=height)
################################################################
hix,hiy = plotmap(lon,lat)
px,py   = plotmap(plon,plat)




################################################################
# FUNCTION FIND NEAREST TIME RECORD TO PIVOT TIME 
################################################################
# dates: list of dates in datetime format
# pivot: ONE date in datetime format
def nearestDate(dates, pivot):
  return min(dates, key=lambda x: abs(x - pivot))
################################################################


