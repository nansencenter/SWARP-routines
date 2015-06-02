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

pol      = shgeom.Polygon(rect0)
pol_orig = pol

# set line to split polygon
L     = 33.4142135624
if 0:
   # split into 2 poly's
   x0       = 1
   y0       = 1.5
   linedir  = np.pi/3.
elif 1:
   # split into 3 poly's
   L        = 28.6728315753
   x0       = -.5
   y0       = -.5
   linedir  = np.pi/3.
elif 1:
   # split into 3 poly's
   L        = 28.6728315753
   x0       = 1.
   y0       = -.5
   linedir  = np.pi/2.
elif 1:
   # start on boundary
   x0       = 1
   y0       = 1
   linedir  = np.pi/2.
else:
   x0       = 3
   y0       = 3
   linedir  = np.pi/3.

lin      = [(x0,y0),(x0+L*np.cos(linedir),y0+L*np.sin(linedir))]
print(lin)
ll = shgeom.LineString(lin)
print(ll)
ii = ll.intersection(pol)

class arctan2_polygon:

   ##################################################################################################
   def __init__(self,pol,coords,linedir):
      import geometry_planar as GP
      import shapefile_utils as SFU

      self.polygon      = pol
      self.branch_point = coords
      self.branch_dir   = linedir
      #
      MPlist,self.intersection_types,\
         self.branch_line = GP.line_splits_polygon(pol,coords,linedir)

      self.connected_polygon  = MPlist[-1]
      self.disjoint_polygons  = MPlist[:-1]

      # get corrections for each disjoint polygon
      self._get_arctan2_constants()

      return
   ##################################################################################################

   ##################################################################################################
   def _get_arctan2_constants(self):

      # disjoint polygons
      MPlist      = self.disjoint_polygons
      Ncl         = len(MPlist)
      isec_types  = self.intersection_types

      # connected polygons
      pol_conn          = self.connected_polygon
      const_atan2       = []
      const_atan2_bdy   = []

      #####################################################################################################
      # loop over disjoint polygons
      for i in range(Ncl):
         pol_disj = MPlist[i]
         int_bdy  = isec_types.crossing_lines[i]
         Li       = int_bdy.length

         # to see jump across int_bdy, get 2 points
         # - one from just inside pol_conn
         # - one from just inside pol_disj
         pic   = int_bdy.centroid
         dsk_c = pic.buffer(1.e-4*Li/2.).intersection(pol_conn)
         dsk_d = pic.buffer(1.e-4*Li/2.).intersection(pol_disj)
         xi,yi = GP.coords2xy(pic.coords)
         xc,yc = GP.coords2xy(dsk_c.centroid.coords)
         xd,yd = GP.coords2xy(dsk_d.centroid.coords)

         Ac = GP.arctan2_branch(yc-y0,xc-x0,branch_dir=linedir)
         Ad = GP.arctan2_branch(yd-y0,xd-x0,branch_dir=linedir)
         Ai = GP.arctan2_branch(yi-y0,xi-x0,branch_dir=linedir)

         ##################################################################################################
         # correction inside the disjoint polygon
         if Ad<Ac-np.pi:
            const_atan2.append(1)   # need to add 2*pi to Ad to make it continuous across the boundary
         elif Ad>Ac+np.pi:
            const_atan2.append(-1)  # need to add -2*pi to Ad to make it continuous across the boundary

         # correction on the boundary itself
         if Ai<Ac-np.pi:
            const_atan2_bdy.append(1)   # need to add 2*pi to Ai to make it continuous across the boundary
         elif Ai>Ac+np.pi:
            const_atan2_bdy.append(-1)  # need to add -2*pi to Ai to make it continuous across the boundary
         else:
            const_atan2_bdy.append(0)  # arctan2_branch is continous approaching the boundary
                                       # (the jump is just after the boundary)
         ##################################################################################################

         self.corrections_internal  = const_atan2
         self.corrections_boundary  = const_atan2_bdy
      # finish loop over disjoint polygons
      #####################################################################################################

      return
      ##################################################################################################

   ##################################################################################################
   def plot_polygons(self):
      import geometry_planar as GP
      import shapefile_utils as SFU


      poly  = self.polygon
      xlim  = [poly.bounds[0]-1.,poly.bounds[2]+1.]
      ylim  = [poly.bounds[1]-1.,poly.bounds[3]+1.]
      plt.xlim(xlim)
      plt.ylim(ylim)
      SFU.plot_poly(poly,color='m')

      for n,poly in enumerate(self.disjoint_polygons):
         print(poly)
         print(poly.is_valid)
         SFU.plot_poly(poly)

      if 1:
         xl,yl = self.branch_line.coords.xy
         plt.plot(xl,yl,'--r')

      for cl in self.intersection_types.crossing_lines:
         xi,yi = cl.coords.xy
         plt.plot(xi,yi,'.')


      if 1:
         plt.show()
      else:
         figname  = 'out/test2.png'
         print('saving to '+figname+'\n')
         plt.savefig(figname)
         plt.close()

      return
   ##################################################################################################

   # ##################################################################################################
   # def eval_solution(self,x,y):
   #    return
   # ##################################################################################################

   ##################################################################################################
   def plot_solution(self,pobj=None,boundary=False,add_corrections=False):

      import numpy as np

      if pobj is None:
         from matplotlib import pyplot as pobj

      poly  = self.polygon
      SFU.plot_poly(poly,color='k')

      if not boundary:
         # plot everywhere
         SFU.plot_poly(poly,color='k',linewidth=2)
         bbox  = poly.bounds

         # get a grid to plot F on:
         x0    = bbox[0]-1
         y0    = bbox[1]-1
         x1    = bbox[2]+1
         y1    = bbox[3]+1
         eps   = np.max([x1-x0,y1-y0])/200.
         x     = np.arange(x0,x1+eps,eps)
         y     = np.arange(y0,y1+eps,eps)

         # make pcolor/contour plot
         X,Y   = np.meshgrid(x,y)
         xc,yc = self.branch_point

         if not add_corrections:
            F  = GP.arctan2_branch(Y-yc,X-xc,self.branch_dir)
         else:
            F  = self.eval_solution(X,Y)

         pobj.pcolor(X,Y,F/np.pi)
         pobj.colorbar()

      if 1:
         plt.show()

      return
   ##################################################################################################


a2p   = arctan2_polygon(pol,(x0,y0),linedir)
# a2p.plot_polygons()
a2p.plot_solution()

if 0:
   MPlist,isec_types = GP.line_splits_polygon(pol,(x0,y0),linedir)

   if 0:
      # plot the polygons
      for n,poly in enumerate(MPlist):
         print(poly,)
         print(poly.is_valid)
         SFU.plot_poly(poly)

      if 0:
         xl,yl = ll.coords.xy
         plt.plot(xl,yl,'k')

      for cl in isec_types.crossing_lines:
         xi,yi = cl.coords.xy
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
   else:
      Ncl               = len(MPlist)-1
      pol_conn          = MPlist[-1]
      const_atan2       = []
      const_atan2_bdy   = []

      for i in range(Ncl):
         pol_disj = MPlist[i]
         int_bdy  = isec_types.crossing_lines[i]
         Li       = int_bdy.length

         # to see jump across int_bdy, get 2 points
         # - one from just inside pol_conn
         # - one from just inside pol_disj
         pic   = int_bdy.centroid
         dsk_c = pic.buffer(1.e-4*Li/2.).intersection(pol_conn)
         dsk_d = pic.buffer(1.e-4*Li/2.).intersection(pol_disj)
         xi,yi = GP.coords2xy(pic.coords)
         xc,yc = GP.coords2xy(dsk_c.centroid.coords)
         xd,yd = GP.coords2xy(dsk_d.centroid.coords)

         Ac = GP.arctan2_branch(yc-y0,xc-x0,linedir=linedir)
         Ad = GP.arctan2_branch(yd-y0,xd-x0,linedir=linedir)
         Ai = GP.arctan2_branch(yi-y0,xi-x0,linedir=linedir)

         # print('**************************************************')
         # print('getting arctan2 constants')
         # print(pol_disj)
         # print(pol_conn)
         # print(pic)
         # print((xc,yc))
         # print((xd,yd))
         # print(' ')
         # print(Ad,Ai,Ac)
         # print('**************************************************\n')

         ##################################################################################################
         # correction inside the disjoint polygon
         if Ad<Ac-np.pi:
            const_atan2.append(1)   # need to add 2*pi to Ad to make it continuous across the boundary
         elif Ad>Ac+np.pi:
            const_atan2.append(-1)  # need to add -2*pi to Ad to make it continuous across the boundary

         # correction on the boundary itself
         if Ai<Ac-np.pi:
            const_atan2_bdy.append(1)   # need to add 2*pi to Ai to make it continuous across the boundary
         elif Ai>Ac+np.pi:
            const_atan2_bdy.append(-1)  # need to add -2*pi to Ai to make it continuous across the boundary
         else:
            const_atan2_bdy.append(0)  # arctan2_branch is continous approaching the boundary
                                       # (the jump is just after the boundary)
         ##################################################################################################

      print(const_atan2)
      print(const_atan2_bdy)
