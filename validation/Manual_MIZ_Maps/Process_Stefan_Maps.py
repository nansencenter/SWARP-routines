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
METH     = 1
"""
0     : direct Laplacian with specified boundary flags
1     : direct Laplacian with boundary flags determined from PCA
2,3   : direct Laplacian to simplified polygon (lower fractal dimension index),
         with boundary flags determined from PCA for new shape
 * 2 > get convex hull
 * 3 > dilation to get less complicated shape (in between original and convex hull)
4     : Use PCA to determine direction to get the width in
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


#######################################################################
def Simplify(lons,lats,bmap,res=10000.,method='ConvexHull'):

   x,y   = bmap(lons,lats)
   xy    = np.array([x,y]).transpose()
   xyc   = [tuple(xyi) for xyi in xy]
   shp   = SHgeom.Polygon(xyc)

   if method=='ConvexHull':
      # get convex hull
      shp2  = shp.convex_hull
      x2,y2 = shp2.exterior.coords.xy
   else:
      # get convex hull
      shp2  = mizc.covering_polygon(shp)
      x2,y2 = shp2.exterior.coords.xy

   # increase resolution (m) (this increases the number of points):
   x3,y3       = mizc.fill_poly(x2,y2,res=res)
   lons2,lats2 = bmap(x3,y3,inverse=True)
   
   # apply Laplacian soln to simplified polygon
   Psoln = Leqs.get_MIZ_widths(lons2,lats2,basemap=bmap)

   ####################################################################
   # restrict contour lines to within original poly
   MIZlines = []
   for llc in Psoln.area_info.lonlat_contours:
      lonv,latv   = np.array(llc).transpose()
      xx,yy       = bmap(lonv,latv)
      xyv         = np.array([xx,yy]).transpose()
      xyv         = [tuple(xyi) for xyi in xyv]
      #
      LS    = SHgeom.LineString(xyv)
      if LS.intersects(shp):
         LSi   = LS.intersection(shp)
         MIZlines.append(mizc.MIZcont(LSi,bmap))
   ####################################################################

   return Psoln,MIZlines
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
         Poly  = FSM.poly_info(llc,bmap,cdate=cdate,func_vals=fvals)
         Polys.append(Poly)
      ############################################################

      ############################################################
      fig      = plt.figure()
      ax1      = fig.add_subplot(1,1,1)
      Psolns   = []
      for Poly in Polys:
         lons,lats   = np.array(Poly.ll_coords).transpose()
         if METH<2:
            if METH==0:
               # direct Laplacian soln
               # - use fvals
               print('\nUsing Laplacian on original polygon, with boundary flags\n')
               Psoln = Leqs.get_MIZ_widths(lons,lats,fvals=Poly.func_vals,basemap=bmap)
            else:
               print('\nUsing Laplacian on original polygon, with PCA\n')
               # - use PCA
               Psoln = Leqs.get_MIZ_widths(lons,lats,basemap=bmap)

            cbar  = (Psolns==[])
            Psoln.plot_soln(pobj=[fig,ax1],bmap=bmap,cbar=cbar)
            Psolns.append(Psoln)

         elif METH<4:
            # apply Laplacian method to simpler covering polygon
            if METH==2:
               method   = 'ConvexHull'
            else:
               method   = 'Buffer'

            print('\nUsing Laplacian on simplified polygon ('+method\
                  +'), with PCA\n')
            Psoln,MIZlines = Simplify(lons,lats,bmap,method=method)
            cbar           = (Psolns==[])
            Psoln.plot_soln(pobj=[fig,ax1],bmap=bmap,cbar=cbar)
            Psolns.append(Psoln)
            #
            x,y   = np.array(Poly.xy_coords).transpose()
            bmap.plot(x,y,'k',linewidth=2,ax=ax1)
            #
            for MIZc in MIZlines:
               MIZc.plot_lines(bmap,ax=ax1,color='c')

         elif METH==4:
            print('\nUsing PCA without Laplacian solution\n')
            x,y   = np.array(Poly.xy_coords).transpose()
            bmap.plot(x,y,'k',linewidth=2,ax=ax1)
            #
            PCA      = mizc.pca_mapper(Poly.xy_coords)
            MIZinfo  = PCA.get_MIZ_lines(bmap)
            MIZinfo.plot_soln(bmap,ax=ax1,color='c')
            Psolns.append(MIZinfo)
      ############################################################

      Fplt.finish_map(bmap,ax=ax1)
      if 1:#show:
         plt.show(fig)
      else:
         figname  = figdir+'/'+cdate+"_polys.png"
         print('>'+figname)
         fig.savefig(figname)
      ax1.cla()
      plt.close(fig)

      if SV_SHP:
         #save shapefile
         mizc.save_shapefile(Psolns,filename='out/test.shp')
