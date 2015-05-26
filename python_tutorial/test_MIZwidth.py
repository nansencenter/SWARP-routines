import os,sys
import numpy            as np
import shapely.geometry as shgeom
import shapely.ops      as shops
from matplotlib import pyplot as plt

if '../py_funs' not in sys.path:
   sys.path.append('../py_funs')

import shapefile_utils as SFU
import geometry_planar as GP

##################################################
class multipole:
   # make a test function that satisfies Laplace's equation
   # - of form r^n*cos(n\theta) or r^n*sin(n\theta)
   def __init__(self,radius,center,even=True):
      self.radius = radius
      self.center = center
      self.even   = even
      return

   def eval_solution(self,x,y,orders,coeffs,polar=False):

      Nx    = x.size
      shp   = x.shape
      out   = np.zeros(Nx)

      if not polar:
         # reshape and make them relative to centre
         x  = x.reshape(Nx)-self.center[0]
         y  = y.reshape(Nx)-self.center[1]
         r  = np.sqrt(x*x+y*y)
         th = np.arctan2(y,x)
      else:
         # just reshape
         r  = x.reshape(Nx)
         th = y.reshape(Nx)
      
      a  = self.radius # scale r by this:
      # a  = 1. # scale r by this:
      if self.even:
         for m,order in enumerate(orders):
            out   = out+coeffs[m]*pow(r/a,order)*np.cos(order*th)
      else:
         for m,order in enumerate(orders):
            out   = out+coeffs[m]*pow(r/a,order)*np.sin(order*th)

      return out.reshape(shp)
##################################################

##################################################
def get_grid(bbox,nx=50,ny=50):
   x0    = bbox[0]
   y0    = bbox[1]
   x1    = bbox[2]
   y1    = bbox[3]
   #
   dx = (x1-x0)/float(nx)
   dy = (y1-y0)/float(ny)
   x  = np.arange(x0,x1+dx,dx)
   y  = np.arange(y0,y1+dy,dy)

   return np.meshgrid(x,y)
##################################################

p1 = shgeom.Point((0,0)).buffer(3)
p2 = shgeom.Point((4,0)).buffer(3)
#
poly  = shops.unary_union([p1,p2])
x,y   = poly.exterior.coords.xy
x     = np.array(x)
y     = np.array(y)
if 0:
   SFU.plot_poly(poly)
   plt.show()

# make a test function that satisfies Laplace's equation
mpol_centre = (2.5,0)
mpol     = multipole(6,mpol_centre)
if 0:
   orders   = [2,3,5]
   coeffs   = [.2,1.3,.5]
else:
   orders   = [3]
   coeffs   = [1.]

# solve Dirichlet problem
func_vals   = mpol.eval_solution(x,y,orders,coeffs)
coords      = GP.xy2coords(x,y)
fun_sol     = GP.dirichlet_fund_soln(coords,func_vals)

fdir  = 'out/funsol'
Vmin  = -.2
Vmax  = .2
dc    = (Vmax-Vmin)/20.
clev  = np.arange(Vmin,Vmax+dc,dc)
if not os.path.exists(fdir):
   os.mkdir(fdir)

if 1:
   fig   = plt.figure()

   # plot test solution:
   circ  = shgeom.Point(mpol_centre).buffer(6)
   SFU.plot_poly(circ,color='k',linewidth=2)
   X,Y   = get_grid(circ.bounds)
   G     = mpol.eval_solution(X,Y,orders,coeffs)
   plt.pcolor(X,Y,G,vmin=Vmin,vmax=Vmax)
   plt.contour(X,Y,G,clev,colors='k')

   # plot boundary of region where Laplace's eqn will be solved
   SFU.plot_poly(poly,color='b',linewidth=2)

   # plot locations of sing's
   xy2   = np.array(fun_sol.singularities)
   plt.plot(xy2[:,0],xy2[:,1],'.k')

   plt.colorbar()
   plt.savefig(fdir+'/test_mpol.png')
   plt.close()
   fig.clf()

if 1:
   fvals_ap = fun_sol.eval_solution(x,y)
   fig   = plt.figure()
   plt.plot(func_vals)
   plt.plot(fvals_ap,'--r')
   plt.savefig(fdir+'/test_bdy.png')
   plt.close()
   fig.clf()

if 1:
   # plot "fun_sol" solution:
   fig   = plt.figure()
   circ  = shgeom.Point(mpol_centre).buffer(6)
   SFU.plot_poly(circ,color='k',linewidth=2)
   X,Y   = get_grid(circ.bounds)
   F     = fun_sol.eval_solution(X,Y)
   plt.pcolor(X,Y,F,vmin=Vmin,vmax=Vmax)
   plt.contour(X,Y,G,clev,colors='k')

   # plot boundary of region where Laplace's eqn will be solved
   SFU.plot_poly(poly,color='b',linewidth=2)

   # plot locations of sing's
   xy2   = np.array(fun_sol.singularities)
   plt.plot(xy2[:,0],xy2[:,1],'.k')

   plt.colorbar()
   plt.savefig(fdir+'/test_full.png')
   plt.close()
   fig.clf()
