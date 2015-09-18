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
import fns_Stefan_Maps as FSM

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

rootdir  = '/Volumes/sim/tim/Projects/SWARP/Stefan'
indir    = rootdir+'/'+fmon+'/'+vsn+'/polygons/txt'

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
bmap  = Basemap(width=2*xmax,height=2*ymax,\
                resolution=cres,projection='stere',\
                lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)
#######################################################################

if 0:
   day0  = 1
   day1  = 28
   show  = False
else:
   # day   = 1
   # day   = 5
   # day   = 8
   day   = 17
   # day   = 18
   # day   = 25
   # day   = 28
   day0  = day
   day1  = day
   show  = True

for iday in range(day0,day1+1):
   cday  = '%2.2d' %(iday)
   cdate = fmon+cday
   fname = indir+'/'+cdate+'_polys.txt'
   if os.path.exists(fname):
      print('\n'+60*'--'+'\nReading '+fname+'\n')

      ############################################################
      # get polys as "poly_info" objects
      Pols  = FSM.read_txt_file_polys(fname)
      Polys = []
      for llc,fvals in Pols:
         Poly  = FSM.poly_info(cdate,llc,bmap,fvals)
         Polys.append(Poly)
      ############################################################

      """
      TODO: run Laplacian on Polys
      - save as .shp files with info
      - plot on big basemap + indicate width on poly
      """

      if 1:
         # TEST PLOT - with polygons
         #
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,1,1)
         lstil = ['-','--']

         for Poly in Polys:
            x,y   = np.array(Poly.xy_coords).transpose()
            shp   = SHgeom.Polygon(Poly.xy_coords)
            if shp.is_valid:
               bmap.plot(x,y,'k',linewidth=2,ax=ax1)
            else:
               bmap.plot(x,y,'c',linewidth=2,ax=ax1)

         Fplt.finish_map(bmap,ax=ax1)

         if 1:#show:
            plt.show(fig)
         else:
            figname  = figdir+'/'+cdate+"_polys.png"
            print('>'+figname)
            fig.savefig(figname)
         ax1.cla()
         plt.close(fig)
