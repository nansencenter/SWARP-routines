import os,sys
import numpy                  as np
import shapely.geometry       as shgeom
import shapely.ops            as shops
from matplotlib import pyplot as plt
from skimage import measure   as msr

if '../py_funs' not in sys.path:
   sys.path.append('../py_funs')

import shapefile_utils as SFU
import geometry_planar as GP
import Laplace_eqn_solution as Lap


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
if 0:
   orders   = [2,3,5]
   coeffs   = [.2,1.3,.5]
   even     = True
else:
   orders   = [1]
   coeffs   = [1.]
   even     = False
mpol        = Lap.multipole(6,mpol_centre,even=even)

# solve Dirichlet problem
func_vals   = mpol.eval_solution(x,y,orders,coeffs)
coords      = GP.xy2coords(x,y)
fun_sol     = Lap.dirichlet_fund_soln(coords,func_vals)

fdir  = 'out/funsol'
Vmin  = func_vals.min()
Vmax  = func_vals.max()

print('\nValue range:')
print(Vmin,Vmax)
print('')

dc    = (Vmax-Vmin)/20.
clev  = np.arange(Vmin,Vmax+dc,dc)
if not os.path.exists(fdir):
   os.mkdir(fdir)

if 0:
   fig   = plt.figure()

   # plot test solution:
   circ  = shgeom.Point(mpol_centre).buffer(6)
   SFU.plot_poly(circ,color='k',linewidth=2)
   X,Y   = get_grid(circ.bounds)
   G     = mpol.eval_solution(X,Y,orders,coeffs)
   plt.pcolor(X,Y,G,vmin=Vmin,vmax=Vmax)
   plt.colorbar()

   plt.contour(X,Y,G,clev,colors='k')

   # plot boundary of region where Laplace's eqn will be solved
   SFU.plot_poly(poly,color='b',linewidth=2)

   # plot locations of sing's
   xy2   = np.array(fun_sol.singularities)
   plt.plot(xy2[:,0],xy2[:,1],'.k')

   figname  = fdir+'/test_mpol.png'
   print('\nSaving to figure '+figname+'\n')
   plt.savefig(figname)
   plt.close()
   fig.clf()

if 0:
   # plot boundary values
   fig   = plt.figure()
   fun_sol.plot_solution(plot_boundary=True,pobj=plt,show=False)
   #   x2,y2       = np.array(fun_sol.coords).transpose() # coords can be reversed
   #   fvals_ap1   = fun_sol.eval_solution(x2,y2)
   #   fvals_ap    = fun_sol.func_vals_approx
   #   fig         = plt.figure()
   #   ss          = fun_sol.get_arc_length()[:-1]
   #   plt.plot(ss,fun_sol.func_vals,'.k')
   #   plt.plot(ss,fvals_ap,'-b')
   #   plt.plot(ss,fvals_ap1,'--r')
   figname  = fdir+'/test_bdy.png'
   print('\nSaving to figure '+figname+'\n')
   plt.savefig(figname)
   plt.close()
   fig.clf()

if 0:
   # plot "fun_sol" solution:
   fig   = plt.figure()
   circ  = shgeom.Point(mpol_centre).buffer(6)
   SFU.plot_poly(circ,color='k',linewidth=2)
   X,Y   = get_grid(circ.bounds)
   F     = fun_sol.eval_solution(X,Y)
   plt.pcolor(X,Y,F,vmin=Vmin,vmax=Vmax)
   plt.colorbar()

   # contours
   plt.contour(X,Y,F,clev,colors='k')

   # plot boundary of region where Laplace's eqn will be solved
   SFU.plot_poly(poly,color='b',linewidth=2)

   # plot locations of sing's
   xy2   = np.array(fun_sol.singularities)
   plt.plot(xy2[:,0],xy2[:,1],'.k')

   figname  = fdir+'/test_full.png'
   print('\nSaving to figure '+figname+'\n')
   plt.savefig(figname)
   plt.close()
   fig.clf()

if 0:
   # test tangent derivative
   F_s   = fun_sol.tangent_deriv
   ss    = np.cumsum(fun_sol.spacings)

   fig   = plt.figure()
   plt.plot(ss,F_s)
   # 
   F_s2  = (fun_sol.func_vals[1:]-fun_sol.func_vals[:-1])/fun_sol.spacings[:-1]
   plt.plot(ss[:-1],F_s2,color='r',linestyle='--')

   figname  = fdir+'/testTanDeriv.png'
   print('\nSaving to figure '+figname+'\n')
   plt.savefig(figname)
   plt.close()
   fig.clf()

if 0:
   # test gradient:
   if not even:
      dy = .1
      yt = np.arange(-2,2+dy,dy)
      xt = 4+0*yt
      f1 = fun_sol.eval_solution(xt,yt)
      plt.plot(yt,f1)
      #
      fx,fy = fun_sol.eval_gradient(xt,yt)
      f2 = f1[0]+0*f1
      for n in range(1,len(yt)):
         f2[n] = f2[n-1]+dy*fy[n]
      plt.plot(yt,f2,'--r')

      f3 = mpol.eval_solution(xt,yt,orders,coeffs)
      plt.plot(yt,f3,'.g')

      figname  = fdir+'/testFy.png'
      print('\nSaving to figure '+figname+'\n')
      plt.savefig(figname)
      plt.close()
      fig.clf()
   else:
      dx = .1
      xt = np.arange(0.5,4+dx,dx)
      yt = 0*xt
      f1 = fun_sol.eval_solution(xt,yt)
      plt.plot(xt,f1)
      #
      fx,fy = fun_sol.eval_gradient(xt,yt)
      # integrate
      f2 = np.cumsum(fx*dx)
      f2 = f2-f2[0]+f1[0]
      plt.plot(xt,f2,'--r')

      f3 = mpol.eval_solution(xt,yt,orders,coeffs)
      plt.plot(xt,f3,'.g')
      #
      figname  = fdir+'/testFx.png'
      print('\nSaving to figure '+figname+'\n')
      plt.savefig(figname)
      plt.close()
      fig.clf()

if 0:
   #test isoline extraction
   if 0:
      merged_levels = fun_sol.get_isolines()

      fig   = plt.figure()
      print('\nplotting isolines...\n')
      poly  = fun_sol.shapely_polygon
      bbox  = poly.bounds
      eps   = fun_sol.resolution/2.
      #
      x0    = bbox[0]
      y0    = bbox[1]
      x1    = bbox[2]
      y1    = bbox[3]
      xp = np.arange(x0-eps/2.,x1+1.5*eps,eps)
      yp = np.arange(y0-eps/2.,y1+1.5*eps,eps)
      SFU.plot_poly(poly,pobj=plt,color='k',linewidth=2)

      for cont in merged_levels:
         xx,yy = np.array(cont).transpose()
         plt.plot(xx,yy)

      plt.show()
      fig.clf()
      # sys.exit()

   elif 0:
      # evaluate solution on the grid
      print('extracting isolines...\n')
      poly  = fun_sol.shapely_polygon
      bbox  = poly.bounds
      eps   = fun_sol.resolution/2.
      #
      x0    = bbox[0]
      y0    = bbox[1]
      x1    = bbox[2]
      y1    = bbox[3]
      #
      xv    = np.arange(x0,x1+eps,eps)
      yv    = np.arange(y0,y1+eps,eps)
      Nx    = len(xv)
      Ny    = len(yv)
      X,Y   = np.meshgrid(xv,yv)
      F     = fun_sol.eval_solution(X,Y)

      # get contours
      nlevels        = fun_sol.number_of_points/3
      vmin           = fun_sol.func_vals.min()
      vmax           = fun_sol.func_vals.max()
      dv             = (vmax-vmin)/float(nlevels)
      vlev           = np.arange(vmin+dv/2.,vmax+dv/2.,dv)
      print(str(nlevels)+'contours, for isolines between '+\
            str(vmin)+' and '+str(vmax))
      #
      contours  = []
      for m,V in enumerate(vlev):
         print('level '+str(m)+' ('+str(V)+'), out of '+str(nlevels))
         #
         B              = np.zeros(F.shape)
         B[F>=V]        = 1.
         B[np.isnan(F)] = np.nan
         conts0         = msr.find_contours(B,.5)  # list of [ivec,jvec] arrays
         conts          = []                       # list of (xvec,yvec) arrays

         ##################################################
         #convert conts0->conts, or (i,j)->(x,y)
         for m2,cont in enumerate(conts0):
            ivec  = cont[:,0]
            jvec  = cont[:,1]
            #
            xvec  = xv[0]+jvec*eps# j=0->x[0],j=1->x[1]??
            yvec  = yv[0]+ivec*eps
            keep  = np.logical_and([not bol for bol in np.isnan(xvec)],\
                                    [not bol for bol in np.isnan(yvec)])

            xk    = xvec[keep]
            yk    = yvec[keep]
            Nok   = len(xk)

            if Nok>1:
               # convert to list
               xy = list(np.array([xk,yk]).transpose())
               xy = [(xx,yy) for xx,yy in xy]

               # end points of contour
               c0 = xy[0]
               cl = xy[-1]

               # find nearest points on boundary to end points of contour
               i0 = list(fun_sol.coord_index.nearest(c0,1))[0]
               il = list(fun_sol.coord_index.nearest(cl,1))[0]
               #
               xy.insert(0,fun_sol.coords[i0])
               xy.append(fun_sol.coords[il])
               conts.append(xy)

               if 0:#m==1:
                  xp = np.arange(x0-eps/2.,x1+1.5*eps,eps)
                  yp = np.arange(y0-eps/2.,y1+1.5*eps,eps)
                  plt.pcolor(xp,yp,B,vmin=-.5,vmax=1.5)
                  plt.plot(xk,yk)
                  plt.show()
                  sys.exit()

         contours.append(conts)
         ##################################################

      ##################################################
      fig   = plt.figure()
      print('\nplotting isolines...\n')
      xp = np.arange(x0-eps/2.,x1+1.5*eps,eps)
      yp = np.arange(y0-eps/2.,y1+1.5*eps,eps)
      SFU.plot_poly(poly,pobj=plt,color='k',linewidth=2)
      plt.pcolor(xp,yp,F,vmin=vmin,vmax=vmax)
      plt.colorbar()

      for conts in contours:
         for cont in conts:
            xx,yy = np.array(cont).transpose()
            plt.plot(xx,yy)

      plt.show()
      fig.clf()
      ##################################################

if 1:
   # test stream function
   stream   = Lap.dirichlet_fund_soln(1*fun_sol.coords,1*fun_sol.stream_func_bdy)
   AI       = stream.get_contour_lengths(pobj=plt)
