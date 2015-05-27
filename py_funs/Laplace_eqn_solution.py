#########################################################
class dirichlet_fund_soln:
   # solve Dirichlet problem to get coefficients of expansion
   # F(z)=\sum_n.a_n.log|z-z_n|
   # (F satisfies Laplace's eqn)
   # this functions solves it and makes the mapping F
   
   #######################################################
   # INITIALISATION
   #######################################################

   #######################################################
   def __init__(self,coords,func_vals,solve_exactly=True):
      # initialise object

      import numpy as np
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html
      import geometry_planar  as GP     # also in py_funs

      x,y         = GP.coords2xy(coords)
      area        = GP.area_polygon_euclidean(x,y)
      self.area   = abs(area)

      # want to go round curve anti-clockwise
      if area<0:
         print("Curve traversed in clockwise direction - reversing arrays' order")
         coords.reverse
         func_vals   = list(func_vals)
         func_vals.reverse

      # don't want to repeat the last point
      if coords[0]!=coords[-1]:
         self.coords    = coords
         self.func_vals = np.array(func_vals)
      else:
         self.coords    = coords[:-1]
         self.func_vals = np.array(func_vals[:-1])

      N0             = len(self.coords)
      self.length    = N0

      # gets spacings,directions between points, and perimeter
      self.perimeter,self.resolution,\
            self.spacings,self.tangent_dirn  = GP.curve_info(self.coords)
      #
      self.shapely_polygon  = shgeom.Polygon(coords)

      # get points around boundary of polygon, then
      # expand points by an amount related to the spacings of the coords
      self.buffer_resolution  = 16 # default buffer resolution
                                   # (number of segments used to approximate a quarter circle around a point.)

      self.solve_exactly   = solve_exactly
         # if True: number of singularities and number of boundary points should be the same

      get_sings = True
      while get_sings:
         # May need a couple of repetitions if too few singularities are found
         get_sings   = self._get_singularities()

      # solve Laplace's eqn (get a_n)
      self.solve_laplace_eqn()

      # # evaluate error on boundary:
      print('\nCalculating error on the boundary...\n')
      self._eval_solution_boundary()

      # evaluate normal derivative -> stream function on boundary:
      print('\nCalculating normal derivative -> stream function on the boundary...\n')
      self._eval_derivs_boundary()

      return
   #######################################################

   #######################################################
   # FUNCTIONS FOR INITIALISATION OF THE OBJECT BELOW:
   #######################################################

   #######################################################
   def _get_singularities(self):

      print('Getting singularities...\n')

      bufres   = self.buffer_resolution
      eps      = self.resolution/2.
      poly     = self.shapely_polygon.buffer(eps,resolution=bufres)

      # put singularities (z_n) on the boundary of expanded polygon
      self.singularities            = poly.exterior.coords
      N1                            = len(self.singularities)
      self.number_of_singularities  = N1

      # check the number of singularities:
      do_check = self._check_singularities()
      # do_check = False # accept automatically

      return do_check
   #######################################################

   # #######################################################
   # move this outside class to make it more general
   # def calc_perimeter(self):
   #    import numpy as np
   #    import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html

   #    p0       = shgeom.Point(self.coords[0])
   #    Nc       = len(self.coords)
   #    P        = 0
   #    spacings = []
   #    for n in range(1,Nc):
   #       p1 = shgeom.Point(self.coords[n])
   #       ds = p0.distance(p1)
   #       P  = P+ds
   #       spacings.append(ds)
   #       p0 = p1

   #    self.spacings     = np.array(spacings)
   #    self.resolution   = np.mean(self.spacings)
   #    self.perimeter    = P
   #    return
   # #######################################################

   #######################################################
   def solve_laplace_eqn(self):
      import numpy as np
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html

      print("\nSolving Laplace's equation...\n")
      Nc = len(self.coords)
      Ns = len(self.singularities)
      M  = np.zeros((Nc,Ns))
      for m,c0 in enumerate(self.coords): # rows of matrix (eqn's)
         p0 = shgeom.Point(c0)
         for n,c1 in enumerate(self.singularities): # cols of matrix (var's)
            p1       = shgeom.Point(c1)
            R        = p0.distance(p1)
            M[m,n]   = np.log(R)

      # solve in a least squares sense:
      # M*a=self.func_vals
      self.sing_coeffs  = np.linalg.lstsq(M,self.func_vals)[0]
      return
   #######################################################

   #######################################################
   # EXTERNAL FUNCTIONS THE CLASS PROVIDES
   #######################################################

   #######################################################
   def _eval_derivs_boundary(self):
      import numpy as np
      import geometry_planar as GP

      # get x,y for boundary
      x,y      = GP.coords2xy(self.coords)
      F_x,F_y  = self.eval_gradient(x,y)

      # get tangent derivative
      ctd   = np.cos(self.tangent_dirn)
      std   = np.sin(self.tangent_dirn)
      F_s   = ctd*F_x+std*F_y

      # get normal derivative
      F_n   = -std*F_x + ctd*F_y

      # get stream function on boundary:
      # G_s = -F_n -> integrate
      G  = np.zeros(len(x))
      Nc = len(x)
      for n in range(1,Nc):
         Gs    = -F_n[n]
         ds    = self.spacings[n]
         G[n]  = G[n-1]+Gs*ds

      self.tangent_deriv   = F_s
      self.normal_deriv    = F_n
      self.stream_func     = G
      return
   #######################################################

   #######################################################
   def _eval_solution_boundary(self):
      import geometry_planar as GP
      
      # get x,y for boundary
      x,y   = GP.coords2xy(self.coords)

      # evaluate func on boundary
      self.func_vals_approx   = self.eval_solution(x,y)

      # calculate error
      self.boundary_error  = GP.vector_norm(self.func_vals-self.func_vals_approx)/GP.vector_norm(self.func_vals)

      return
   #######################################################

   #######################################################
   def eval_solution(self,x,y):

      import numpy as np
      from shapely.prepared import prep   # want "contains" function
      import shapely.geometry as shgeom

      Nx    = x.size
      F     = np.zeros(Nx)
      shp   = x.shape
      x     = x.reshape(Nx)
      y     = y.reshape(Nx)

      poly2 = prep(self.shapely_polygon)
      for m in range(Nx):
         # loop over points to evaluate F at
         p0 = shgeom.Point((x[m],y[m]))

         # check if p0 is inside the domain
         if poly2.intersects(p0):
            # add contribution from each singularity
            for n,c1 in enumerate(self.singularities): # singularities
               p1    = shgeom.Point(c1)
               R     = p0.distance(p1)
               F[m]  = F[m]+self.sing_coeffs[n]*np.log(R)
         else:
            F[m]  = np.nan

      return F.reshape(shp)
   #######################################################

   #######################################################
   def eval_gradient(self,x,y):

      import numpy as np
      from shapely.prepared import prep   # want "contains" function
      import shapely.geometry as shgeom

      Nx    = x.size
      shp   = x.shape
      F_x   = np.zeros(Nx)
      F_y   = np.zeros(Nx)
      x     = x.reshape(Nx)
      y     = y.reshape(Nx)

      poly2 = prep(self.shapely_polygon) # need to test if (x,y) are inside region
      for m in range(Nx):
         # loop over points to evaluate F_x,F_y at
         p0 = shgeom.Point((x[m],y[m]))

         # check if p0 is inside the domain
         if poly2.intersects(p0):
            # add contribution from each singularity
            for n,c1 in enumerate(self.singularities): # singularities
               dx       = x[m]-c1[0]
               dy       = y[m]-c1[1]
               R2       = dx*dx+dy*dy # R^2
               F_x[m]   = F_x[m]+self.sing_coeffs[n]*dx/R2
               F_y[m]   = F_y[m]+self.sing_coeffs[n]*dy/R2
         else:
            F_x[m]  = np.nan
            F_y[m]  = np.nan

      return F_x.reshape(shp),F_y.reshape(shp)
   #######################################################

   #######################################################
   def plot_solution(self,pobj=None,plot_boundary=True,show=True):
      import numpy            as np
      import shapefile_utils  as SFU
      import geometry_planar  as GP

      if pobj is None:
         # set a plot object if none exists
         from matplotlib import pyplot as pobj
      
      poly  = self.shapely_polygon
      if plot_boundary:
         # just plot values of F at boundary
         pobj.plot(self.spacings,self.func_vals,'.k',markersize=1.5)
         # 
         x,y   = poly.exterior.boundary.coords.xy
         F     = self.eval_function(x,y)
         pobj.plot(ss,self.func_vals,'b')
      else:
         # plot F everywhere
         SFU.plot_poly(poly,color='k',linewidth=2)
         bbox  = poly.bounds
         eps   = self.resolution/2.

         # get a grid to plot F on:
         x0 = bbox[0]
         y0 = bbox[1]
         x1 = bbox[2]
         y1 = bbox[3]
         x  = np.arange(x0,x1+eps,eps)
         y  = np.arange(y0,y1+eps,eps)

         # make pcolor/contour plot
         nlevels  = 10
         vmin     = self.func_vals.min()
         vmax     = self.func_vals.max()
         dv       = (vmax-vmin)/(nlevels)
         vlev     = np.arange(vmin,vmax+dv,dv)
         #
         X,Y      = np.meshgrid(x,y)
         F        = self.eval_solution(X,Y)
         pobj.pcolor(X,Y,F,vmin=vmin,vmax=vmax)
         pobj.colorbar()
         pobj.contour(X,Y,F,vlev,colors='k')
         #
         x,y   = GP.coords2xy(self.singularities)
         pobj.plot(x,y,'.k',markersize=1.5)

      if show:
         pobj.show()

      return
   #######################################################

   #######################################################
   def _check_singularities(self):
      # don't want too many singularities
      import numpy as np
      import geometry_planar  as GP

      N0 = self.length
      N1 = self.number_of_singularities

      # set limits for N1
      Nthresh  = 100
      frac     = .2

      if self.solve_exactly:
         # number of singularities and number of boundary points should be the same
         Ntarget  = N0
      elif N0<=Nthresh:
         # if N0<=Nthresh, try to get N1~N0
         Ntarget  = N0
      else:
         # try to get N1~frac*N0, if frac*N0<Nthresh
         Ntarget  = int(np.max(np.round(frac*N0),Nthresh))

      print('\nChecking singularities...\n')
      print('Number of boundary points       : '+str(N0))
      print('Number of singularities         : '+str(N1))
      print('Desired number of singularities : '+str(Ntarget))

      check_again = False
      # NB this applies to the case where N1==Ntarget
      # NB also if N1==Ntarget, we can reduce N1 to Ntarget exactly in 1 go

      if N1>Ntarget:

         ##########################################################################
         # reduce the number of sing's by increasing spacing between points
         coords                        = self.singularities
         perimeter,resolution,\
            spacings,tangent_dirn  = GP.curve_info(coords)
         #
         s_target    = N1/float(Ntarget)*resolution # increase mean spacing between points
         new_coords  = [coords[0]]
         #
         ss = 0
         s0 = 0
         s1 = s0+s_target
         for n,c0 in enumerate(coords[1:]):
            ss = ss+spacings[n]
            if ss>=s1:
               new_coords.append(c0)
               s0 = s1
               s1 = s0+s_target

         # no need to call _get_singularities again in this case
         N2 = len(new_coords)
         if not self.solve_exactly:

            #######################################################################
            # update list of singularities:
            self.singularities            = new_coords
            self.number_of_singularities  = N2
            #######################################################################

         elif N2>N0:

            #######################################################################
            # if self solve_exactly:
            # delete a few points randomly from new_coords to get right number
            Ndel  = N2-N0
            while Ndel>0:
               idel        = np.random.randint(N2)
               new_coords  = list(new_coords).remove(new_coords[idel])
               N2          = len(new_coords)
               Ndel        = N2-N0

            # update list of singularities
            self.singularities            = np.array(new_coords)
            self.number_of_singularities  = N2
            #######################################################################

         ##########################################################################
         else:

            #######################################################################
            # if self solve_exactly:
            # get some more points from self.singularities
            # *this adds more at start preferentially
            # but shouldn't matter since it won't be too many
            Nadd  = N0-N2
            iadd  = 0

            for coord in self.singularities:
               if coord in new_coords:
                  # in there so now need to insert before next element of new_coords
                  iadd  = iadd+1
               else:
                  # adds coord before new_coords[iadd]
                  new_coords  = list(new_coords).insert(iadd,coord)
                  N2          = len(new_coords)
                  Nadd        = N0-N2

                  # new element so still need to increase point at which to place 
                  iadd  = iadd+1

               if Nadd==0:
                  break

            # update list of singularities
            self.singularities            = np.array(new_coords)
            self.number_of_singularities  = N2
            #######################################################################

         ##########################################################################
         print('New number of singularities     : '+str(Ntarget)+'\n')
         ##########################################################################

      elif N1<Ntarget:

         ##########################################################################
         # set up for another iteration
         fac                     = int(np.ceil(1.2*Nmax/float(N1)))
         self.buffer_resolution  = fac*self.buffer_resolution

         # need to call _get_singularities again,
         # to reset the singularities and check them
         check_again = True
         print('Trying again to get singularities (too few)...\n')
         ##########################################################################

      return check_again
   #######################################################

##########################################################

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

      import numpy as np

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
