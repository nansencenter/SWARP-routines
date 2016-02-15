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
   if ii.geom_type=='LineString':
      Ni = 1 # 2 intersection points (split into 2 polygons)

      bb       = pol.boundary.union(ii)
      MPlist   = []
      MPlist2  = []
      P0       = pol.boundary.coords[0]

      for n,gg in enumerate(bb.geoms):

         print(gg)
         print(gg!=ii)
         if gg!=ii:
            if P0 not in gg.coords:
               # normal polygon
               print('normal polygon')
               MPlist.append(shgeom.Polygon(gg))
            else:
               print('contains '+str(P0))
               # boundary can be artificially split by 1st point also
               # this list should be merged into one polygon later
               MPlist2.extend(list(gg.coords))

         if PLOTunionB:
            xb,yb = gg.coords.xy
            plt.plot(xb,yb)

      if PLOTunionB:
         plt.show()
         sys.exit()

   else:
      Ni = len(ii.geoms)   # 2.Ni intersection points (split into Ni+1 polygons)
      bb = pol.boundary
      for i in range(len(ii.geoms)):
         bb = bb.union(ii.geoms[i])

      MPlist   = []
      MPlist2  = []
      P0       = pol.boundary.coords[0]

      for n,gg in enumerate(bb.geoms):

         print(gg)
         if gg not in ii.geoms:
            if P0 not in gg.coords:
               # normal polygon
               print('normal polygon')
               MPlist.append(shgeom.Polygon(gg))
            else:
               print('contains '+str(P0))
               # boundary can be artificially split by 1st point also
               # this list should be merged into one polygon later
               MPlist2.extend(list(gg.coords))

         if PLOTunionB:
            xb,yb = gg.coords.xy
            plt.plot(xb,yb)

      MPlist.append(shgeom.Polygon(MPlist2))

      if PLOTunionB:
         plt.xlim([-1, 4])
         plt.ylim([-1, 4])
         plt.show()
         sys.exit()

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


if 0:
   if 0:
      # use intersection calc' outside of "GP.line_intersection_polygon"
      isec  = ii
   else:
      isec  = GP.line_intersection_polygon(rect0,lin[0],linedir)

   if isec.geom_type=='LineString':
      # 1 intersection
      print(isec.coords)
      Ni = 1
   elif isec.geom_type=='Point':
      Ni = -1
   else:
      Ni = len(isec.geoms)
      if Ni>0:
         # just use 1 pair (1st/last):
         c0    = isec.geoms[0].coords[0]
         c1    = isec.geoms[-1].coords[-1]
         isec  = shgeom.LineString([c0,c1])

   xp,yp = GP.coords2xy(rect0)
   plt.plot(xp,yp,'b')

   xl,yl = ll.coords.xy
   plt.plot(xl,yl,'k')

   xi,yi = isec.coords.xy
   plt.plot(xi,yi,'.g')

   if Ni>0:
      if 1:
         print('shortest')
         clst  = GP.line_splits_polygon(pol,isec,shortest=True)
      else:
         print('longest')
         clst  = GP.line_splits_polygon(pol,isec,shortest=False)
      print(clst)

      xs,ys = GP.coords2xy(clst)
      plt.plot(xs,ys,'r',linewidth=2)
   elif Ni==-1:
      print('Point intersection')
   else:
      print('No intersection')

   plt.xlim([-1, 4])
   plt.ylim([-1, 4])

   if SHOW:
      plt.show()
   else:
      figname  = 'out/test.png'
      print('saving to '+figname+'\n')
      plt.savefig(figname)
      plt.close()

   r,th  = GP.xy2polar_coords(xs,ys)
   plt.plot(r,th/np.pi)
   plt.xlabel('r')
   plt.ylabel('$\\theta/\pi$')
   plt.show()
