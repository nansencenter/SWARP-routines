import os,sys
import numpy as np
import shapely.geometry as shgeom
from matplotlib import pyplot as plt

sys.path.append('../py_funs')

import geometry_planar as GP
import shapefile_utils as SFU

SHOW        = 0
PLOTunionB  = 1

if 0:
   rect0 = [(0,0.1),(2,0),(3,1),(0,1),(0,0.1)]
else:
   rect0 = [(0,0),(3.2,0),(3,3),(0,3),(0,2),(2,2),(2,1),(0,1),(0,0)]

pol   = shgeom.Polygon(rect0)

# set line to split polygon
L     = 33.4142135624
if 0:
   # split into 2 poly's
   x0       = 1
   y0       = 1.5
   linedir  = np.pi/3.
   lin      = [(x0,y0),(x0+L*np.cos(linedir),y0+L*np.sin(linedir))]
elif 1:
   # split into 3 poly's
   x0       = -.5
   y0       = -.5
   linedir  = np.pi/3.
   lin      = [(x0,y0),(x0+L*np.cos(linedir),y0+L*np.sin(linedir))]
elif 1:
   # start on boundary
   x0       = 1
   y0       = 1
   linedir  = np.pi/2.
   lin      = [(x0,y0),(x0+L*np.cos(linedir),y0+L*np.sin(linedir))]
else:
   x0       = 3
   y0       = 3
   linedir  = np.pi/3.
   lin      = [(x0,y0),(x0+L*np.cos(linedir),y0+L*np.sin(linedir))]

ll = shgeom.LineString(lin)
ii = ll.intersection(pol)

if 1:
   #########################################################################
   if ii.geom_type=='LineString':
      Ni       = 1 # 2 intersection points (split into 2 polygons)
      MPlist   = GP.line_splits_polygon_in2(pol,ii)

   else:
      MPlist      = []
      Ni          = len(ii.geoms)   # 2.Ni intersection points (split into Ni+1 polygons)
      pol1,pol2   = GP.line_splits_polygon_in2(pol,ii.geoms[0])
         # split pol into 2 poly's

      ################################################################
      # loop over remaining geometries
      for i in range(1,Ni):
         pol3,pol4   = GP.line_splits_polygon_in2(pol,ii.geoms[i])
         # split pol into another 2 poly's
         if not pol3.intersects(pol1):
            # pol3,pol1 are disjoint polygons
            MPlist.extend([pol3,pol1])

            # => pol2,pol4 are big ones
            # => reduce size of big polygon "pol"
            pol   = pol2.intersection(pol4)

         elif not pol4.intersects(pol1):
            # pol4,pol1 are disjoint polygons
            MPlist.extend([pol4,pol1])

            # => pol2,pol3 are big ones
            # => reduce size of big polygon "pol"
            pol   = pol2.intersection(pol3)

         elif not pol4.intersects(pol2):
            # pol4,pol2 are disjoint polygons
            MPlist.extend([pol4,pol2])

            # => pol1,pol3 are big ones
            # => reduce size of big polygon "pol"
            pol   = pol1.intersection(pol3)

         else:
            # pol3,pol2 are disjoint polygons
            MPlist.extend([pol3,pol2])

            # => pol1,pol4 are big ones
            # => reduce size of big polygon "pol"
            pol   = pol1.intersection(pol4)
      ################################################################

      # finished looping: pol is last polygon
      MPlist.insert(0,pol)
   #########################################################################


   if 1:
      for n,poly in enumerate(MPlist):
         print(poly,)
         print(poly.is_valid)
         SFU.plot_poly(poly)

      if 0:
         xl,yl = ll.coords.xy
         plt.plot(xl,yl,'k')

      if Ni==1:
         xi,yi = ii.coords.xy
         plt.plot(xi,yi,'.r')
      else:
         for i in range(Ni):
            xi,yi = ii.geoms[i].coords.xy
            plt.plot(xi,yi,'.')

      plt.xlim([-1, 4])
      plt.ylim([-1, 4])
      if SHOW:
         plt.show()
      else:
         figname  = 'out/test2.png'
         print('saving to '+figname+'\n')
         plt.savefig(figname)
         plt.close()
