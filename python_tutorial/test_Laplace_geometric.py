import os,sys
import shapely.geometry as shgeom
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt

sys.path.append('../py_funs')
import Laplace_eqn_solution as Leqs

circ  = shgeom.Point((0,0)).buffer(5.e3)
xy0   = np.array(circ.boundary.coords).transpose()
x0,y0 = xy0

if 1:
   # rotate by 45 degrees, and stretch x axis
   s2    = np.sqrt(2)
   Rmat  = np.array([[1/s2,-1/s2],[1/s2,1/s2]])
   Dmat  = np.diag([1.,2.])
   xy2   = Rmat.dot(Dmat.dot(xy0))
   x0,y0 = xy2
   x0    = x0+15.e3

bmap        = Basemap(projection='stere',lat_0=90.,lon_0=-45.,lat_ts=90.,height=100.e3,width=100.e3)
lons,lats   = bmap(x0,y0,inverse=True)

####################################################################################
if 1:
   AI,fun_sol,stream = Leqs.get_MIZ_widths(lons,lats)
else:
   import geometry_sphere as GS # from py_funs
   import f_vals_smoother as smt	

   R    = 6378273.# default is hyc2proj radius (m)
   AL   = GS.arc_length(lons,lats,R=R,radians=False,closed=True)
   P    = AL[-1]
   #
   latc = np.mean(lats)
   lonc = np.mean(lons)   # TODO 1st should make lon continuous around polygon
                        # - to make sure lonc is inside polygon
   basemap = Basemap(lon_0=lonc,lat_0=latc,lat_ts=latc,\
                   projection='stere',rsphere=[R,R],\
                   width=P,height=P)
   # x,y         = basemap(lons,lats)
   x,y         = bmap(lons,lats)
   xy_coords2  = np.array([x,y]).transpose()
   #
   PCA    = smt.pca_mapper(xy_coords2)
   fvals2 = PCA.set_func_vals()

   if 0:
      # test forward map
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      ax.plot(PCA.X,PCA.Y)
      # x2,y2 = PCA.map(x0,y0,inverse=False)
      x2,y2 = PCA.map(PCA.x,PCA.y,inverse=False)
      ax.plot(x2,y2,'--r')
      ax.set_aspect('equal')
      fig.show()

   if 0:
      # test backwards map
      fig2  = plt.figure()
      ax2   = fig2.add_subplot(111)
      ax2.plot(PCA.x,PCA.y)
      # ax2.plot(x0,y0,'--m')
      x1,y1 = PCA.map(PCA.X,PCA.Y,inverse=True)
      ax2.plot(x1,y1,'--r')
      ax2.set_aspect('equal')
      fig2.show()

   if 0:
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      ax.plot(fvals2)
      fig.show()

   # solve Laplace's eqn (2x)
   fun_sol  = Leqs.dirichlet_fund_soln(xy_coords2,fvals2)
   stream   = Leqs.dirichlet_stream_func(potential=fun_sol)
   AI       = stream.get_contour_lengths()
####################################################################################

if 1:
   fig   = plt.figure()
   ax    = fig.add_subplot(111)
   fun_sol.plot_solution(plot_boundary=False,pobj=[fig,ax],show=False)
   #
   for cc in AI.xy_contours:
      xx,yy = np.array(cc).transpose()
      ax.plot(xx/1.e3,yy/1.e3,'k',linewidth=2)

   ax.title.set_text('median contour length (km): '+str(AI.length_median/1.e3))
   ax.set_aspect('equal')
   fig.show()
