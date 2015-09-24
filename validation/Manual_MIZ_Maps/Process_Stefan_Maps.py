import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import rtree.index	as Rindex
from datetime import datetime
import shapely.geometry as SHgeom

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import fns_plotting as Fplt
import geometry_planar as GP
import geometry_sphere as GS
import fns_Stefan_Maps as FSM
import Laplace_eqn_solution as Leqs
import MIZchar as mizc

# Choose method (METH)
METH  = 2
# METH  = 4
"""
0     : direct Laplacian with specified boundary flags
1     : direct Laplacian with boundary flags determined from PCA
2,3   : direct Laplacian to simplified polygon (lower fractal dimension index),
         with boundary flags determined from PCA for new shape
 * 2 > get convex hull
 * 3 > dilation to get less complicated shape (in between original and convex hull)
4     : Use PCA to determine direction to get the width in,
         then just take straight lines across (in stereographic projection space)
"""

# save/don't save results to shapefile (1/0)
SV_SHP   = 1

if 1:
   fmon        = '201402'
   fix_invalid = True
   check_ends  = True
   use_thresh  = False
   if 0:
      vsn   = 'SWARP_MIZ_maps_def'
   elif 0:
      vsn   = 'SWARP_MIZ_maps_def_v2'
   else:
      vsn   = 'SWARP_MIZ_maps_def_new'
else:
   fmon        = '201309'
   vsn         = 'SWARP_MIZ_maps_def'
   fix_invalid = True
   check_ends  = False
   use_thresh  = True

# rootdir  = '/Volumes/sim/tim/Projects/SWARP/Stefan'
rootdir  = '/Users/timill/Documents/SWARP/Stefan'
indir    = rootdir+'/'+fmon+'/'+vsn+'/polygons/txt'

# location of outputs
outdir   = rootdir+'/'+fmon+'/'+vsn+'/polygons/out'
if not os.path.exists(outdir):
   os.mkdir(outdir)
figdir   = outdir+'/figs'
if not os.path.exists(figdir):
   os.mkdir(figdir)
shpdir   = outdir+'/shps'
if not os.path.exists(shpdir):
   os.mkdir(shpdir)

#######################################################################
# make standard stereographic basemap
# - for plotting and also for limiting search area
# lat_ts is latitude of true scale.
# (lon_0,lat_0) is central point -> (x,y)=(0,0)
rad   = 10.          # approx radius of image (degrees)
xmax  = rad*111.e3   # half width of image [m]
ymax  = rad*111.e3   # half height of image [m]
cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
#
if fmon=='201402':
   lon_0    = 36.    # deg E
elif fmon=='201309':
   lon_0    = 5.    # deg E
lat_0    = 80.    # deg N
lat_ts   = lat_0  # deg N
bmap     = Basemap(width=2*xmax,height=2*ymax,\
              resolution=cres,projection='stere',\
              lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)

if 1:
   day0  = 1
   #day0  = 25
   day1  = 28
   show  = False
else:
   # day   = 1
   # day   = 5
   # day   = 8
   # day   = 17
   # day   = 18
   day   = 26
   # day   = 28
   day0  = day
   day1  = day

   show  = True
   # show  = False

for iday in range(day0,day1+1):
   cday  = '%2.2d' %(iday)
   cdate = fmon+cday
   fname = indir+'/'+cdate+'_polys.txt'
   if os.path.exists(fname):
      print('\n'+60*'--'+'\nReading '+fname+'\n')
      fig      = plt.figure()
      ax1      = fig.add_subplot(1,1,1)
      Psolns   = mizc.single_file(fname,bmap,pobj=[fig,ax1],\
                                    cdate=cdate,METH=METH)

      ############################################################
      # finish off fig and show/save
      Fplt.finish_map(bmap,ax=ax1)
      if show:
         plt.show(fig)
      else:
         figname  = figdir+'/'+cdate+"_polys.png"
         print('>'+figname)
         fig.savefig(figname)
      ax1.cla()
      plt.close(fig)
      ############################################################

      if SV_SHP:
         #save shapefile
         # sname = 'out/test.shp'
         sname = shpdir+'/'+cdate+'_polys.shp'
         mizc.save_shapefile(Psolns,filename=sname)
