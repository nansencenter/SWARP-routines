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
   def __init__(self,coords,func_vals,singularities=None):
      # initialise object

      import numpy as np
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html
      import geometry_planar  as GP     # also in py_funs

      #######################################################
      if (coords[0][0]!=coords[-1][0]) and (coords[0][1]!=coords[-1][1]):
         # can't repeat points during Laplace's eqn solution
         self.func_vals = np.array(func_vals)

         # BUT need last and first coordinate the same for shapely polygon
         coords.append(coords[0])
         self.coords = coords[:-1]
         print(len(self.coords))
         print(len(self.func_vals))
      else:
         self.coords    = coords[:-1]
         self.func_vals = np.array(func_vals[:-1])
      #######################################################

      # make shapely polygon
      # (coords now has end-point repeated,
      #  self.coords doesn't)
      self.shapely_polygon  = shgeom.Polygon(coords)
      
      #######################################################
      # want to go round curve anti-clockwise
      x,y         = GP.coords2xy(self.coords)
      area        = GP.area_polygon_euclidean(x,y)
      self.area   = abs(area)

      if area<0:
         print("Curve traversed in clockwise direction - reversing arrays' order")
         self.func_vals = list(self.func_vals)
         self.func_vals.reverse()
         self.func_vals = np.array(self.func_vals)
         #
         self.coords    = list(self.coords)
         self.coords.reverse()
         self.coords    = self.coords
      #######################################################

      self.number_of_points   = len(self.coords)

      # gets spacings,directions between points, and perimeter
      self.perimeter,self.resolution,\
            self.spacings,self.tangent_dirn  = GP.curve_info(self.coords)
      #

      # get points around boundary of polygon, then
      # expand points by an amount related to the spacings of the coords
      self.buffer_resolution  = 16 # default buffer resolution
                                   # (number of segments used to approximate a quarter circle around a point.)

      # self.solve_exactly   = solve_exactly
      # # if True: number of singularities and number of boundary points should be the same

      if singularities is None:
         get_sings = True
      else:
         get_sings                     = False
         self.singularities            = singularities 
         self.number_of_singularities  = len(self.singularities)
         print('Number of boundary points : '+str(self.number_of_points))
         print('Number of singularities   : '+str(self.number_of_singularities)+'\n')

      while get_sings:
         # May need a couple of repetitions if too few singularities are found
         get_sings   = self._get_singularities()

      # solve Laplace's eqn (get a_n)
      self._solve_laplace_eqn()

      # # evaluate error on boundary:
      print('\nCalculating error on the boundary...')
      self._eval_solution_boundary()
      print(str(self.boundary_error)+'\n')

      if 0:
         # evaluate normal derivative -> stream function on boundary:
         print('\nCalculating normal derivative -> stream function on the boundary...\n')
         self._eval_derivs_boundary()

      elif 1:
         # use analytical definition of stream function
         # - for log, this is arctan2 (+const)
         self._get_stream_func_bdy()

      return
   #######################################################

   #######################################################
   # FUNCTIONS FOR INITIALISATION OF THE OBJECT BELOW:
   #######################################################

   #######################################################
   def _get_stream_func_bdy(self):

      # calculate stream function on boundary
      # - continuous on boundary (also should be periodic)
      
      import numpy as np
      import geometry_planar as GP

      # coordinates of polygon
      x,y   = np.array(self.coords).transpose()
      Nx    = len(x)
      nvec  = np.arange(Nx)
      sfun  = 0*x

      # for each singularity (just outside polygon)
      perimeter,resolution,\
            spacings,tangent_dirn  = GP.curve_info(self.singularities)
      for i,sing in enumerate(self.singularities):
         an          = self.sing_coeffs[i]
         branch_dir  = np.pi/2.+tangent_dirn[i]
         atan2       = GP.arctan2_branch(y,x=x,branch_point=sing,branch_dir=branch_dir)
         atan2       = GP.make_arctan2_cts(atan2)
         sfun        = sfun+an*atan2

      self.stream_func_bdy = sfun
      return
   #######################################################

   #######################################################
   def _get_singularities(self):

      print('Getting singularities...\n')

      bufres   = self.buffer_resolution
      eps      = self.resolution/2.
      poly     = self.shapely_polygon.buffer(eps,resolution=bufres)

      # put singularities (z_n) on the boundary of expanded polygon
      self.singularities   = list(poly.exterior.coords)

      if self.singularities[0]==self.singularities[-1]:
         self.singularities   = self.singularities[:-1]

      self.number_of_singularities  = len(self.singularities)

      if 1:
         # check the number of singularities:
         do_check = self._check_singularities()
      else:
         # accept automatically
         do_check = False
         print('Warning: not checking singularities')
         print('Number of boundary points : '+str(self.number_of_points))
         print('Number of singularities   : '+str(self.number_of_singularities)+'\n')

      return do_check
   #######################################################

   #######################################################
   def _check_singularities(self):
      # don't want too many singularities
      import numpy as np
      import geometry_planar  as GP

      N0 = self.number_of_points
      N1 = self.number_of_singularities

      # set limits for N1
      Nthresh  = 100
      frac     = 1.2
      # frac     = .2

      # if self.solve_exactly:
      #    # number of singularities and number of boundary points should be the same
      #    # - this doesn't work too well - need more sing's
      #    Ntarget  = N0
      if N0<=Nthresh:
         # if N0<=Nthresh, try to get N1~N0
         Ntarget  = N0+30
      else:
         # try to get N1~frac*N0, if frac*N0<Nthresh
         Ntarget  = int(np.max(np.round(frac*N0),Nthresh))

      print('\nChecking singularities...\n')
      print('Number of boundary points       : '+str(N0))
      print('Number of singularities         : '+str(N1))
      print('Desired number of singularities : '+str(Ntarget))

      check_again = False
      # NB this applies to the case where N1==Ntarget
      # NB also if N1>=Ntarget, we can reduce N1 to Ntarget exactly in 1 go

      if N1>Ntarget:

         ##########################################################################
         # reduce the number of sing's by increasing spacing between points
         coords   = self.singularities
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

         #######################################################################
         # update list of singularities:
         N2                            = len(new_coords)
         self.singularities            = new_coords
         self.number_of_singularities  = N2
         #######################################################################

         if N2>Ntarget:

            #######################################################################
            # delete a few points randomly from new_coords to get right number
            Ndel  = N2-Ntarget
            while Ndel>0:
               idel  = np.random.randint(N2)
               new_coords.remove(new_coords[idel])
               N2    = len(new_coords)
               Ndel  = N2-Ntarget

            # update list of singularities
            self.singularities            = new_coords
            self.number_of_singularities  = N2
            #######################################################################

         ##########################################################################
         elif N2<Ntarget:

            #######################################################################
            # get some more points from self.singularities
            # *this adds more at start preferentially
            # but shouldn't matter since it won't be too many
            Nadd  = Ntarget-N2
            iadd  = 0

            for coord in self.singularities:
               if coord in new_coords:
                  # in there so now need to insert before next element of new_coords
                  iadd  = iadd+1
               else:
                  # adds coord before new_coords[iadd]
                  new_coords  = list(new_coords)
                  new_coords.insert(iadd,coord)
                  N2          = len(new_coords)
                  Nadd        = Ntarget-N2

                  # new element so still need to increase point at which to place 
                  iadd  = iadd+1

               if Nadd==0:
                  break

            # update list of singularities
            self.singularities            = new_coords
            self.number_of_singularities  = N2
            #######################################################################

         ##########################################################################
         print('New number of singularities     : '+str(Ntarget)+'\n')
         ##########################################################################

      elif N1<Ntarget:

         ##########################################################################
         # set up for another iteration
         # TODO this may not be the best way
         # - perhaps add more points on line segments of polygon
         fac                     = int(np.ceil(1.2*Nmax/float(N1)))
         self.buffer_resolution  = fac*self.buffer_resolution

         # need to call _get_singularities again,
         # to reset the singularities and check them
         check_again = True
         print('Trying again to get singularities (too few)...\n')
         ##########################################################################

      return check_again
   #######################################################

   #######################################################
   def _solve_laplace_eqn(self):
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
   # EXTERNAL FUNCTIONS THE CLASS PROVIDES
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
   def get_arc_length(self):

      import numpy   as np

      ss = list(np.cumsum(self.spacings))
      ss.insert(0,0)
      ss = np.array(ss)

      return ss
   #######################################################

   #######################################################
   def plot_solution(self,pobj=None,plot_boundary=True,show=True):
      import numpy            as np
      import shapefile_utils  as SFU
      import geometry_planar  as GP

      if pobj is None:
         # set a plot object if none exists
         from matplotlib import pyplot as pobj
      
      if plot_boundary:
         # just plot values of F at boundary
         ss = self.get_arc_length()[:-1]
         pobj.plot(ss,self.func_vals_approx,'b')
         pobj.plot(ss,self.func_vals,'.k')

         x2,y2 = np.array(self.coords).transpose() # coords can be reversed
         f2    = self.eval_solution(x2,y2)
         pobj.plot(ss,f2,'--r')
         pobj.xlabel('arc length from '+str(self.coords[0]))
         pobj.ylabel('values of target function')

      else:
         # plot F everywhere
         poly  = self.shapely_polygon
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
   def plot_stream_func_bdy(self,pobj=None,show=True,i=None):
      import numpy            as np
      import shapefile_utils  as SFU
      import geometry_planar  as GP

      if pobj is None:
         # set a plot object if none exists
         from matplotlib import pyplot as pobj
      
      ss = list(np.cumsum(self.spacings[:-1]))
      ss.insert(0,0)
      ss = np.array(ss)

      if i is None:
         # just plot values of stresm function at boundary
         pobj.plot(ss,self.stream_func_bdy,'b')
         out   = None

      else:

         # coordinates of polygon
         x,y   = np.array(self.coords).transpose()
         Nx    = len(x)
         nvec  = np.arange(Nx)

         # for each singularity (just outside polygon)
         sing        = self.singularities[i]
         branch_dir  = np.pi/2.+self.tangent_dirn[i]
         atan2       = GP.arctan2_branch(y,x=x,branch_point=sing,branch_dir=branch_dir)
         pobj.plot(ss,atan2/np.pi,'r')
         out   = [atan2]
         atan2 = GP.make_arctan2_cts(atan2)
         pobj.plot(ss,atan2/np.pi,'b')
         out.append(atan2)

      if show:
         pobj.show()

      return out
   #######################################################
##########################################################

#######################################################
def dirichlet_stream_func(coords=None,potential=None,func_vals=None):
   # CALL: stream = dirichlet_stream_func(potential=potential)
   # inputs:
   # > potential: dirichlet_fund_soln object
   #     - solution to Dirichlet problem with boundary values potential.func_vals
   #
   # OR:
   #
   # CALL: stream = dirichlet_stream_func(coords=coords,func_vals=None):
   # inputs:
   # > coords: list of points (polygon boundary)
   # > func_vals: value of function at each point in coords
   # use these to calculate "potential" object
   #     - solution to Dirichlet problem with boundary values func_vals
   # 
   # outputs:
   # > stream: dirichlet_fund_soln object
   #     -  solution to Dirichlet problem with boundary values potential.stream_func_bdy

   if potential is None:
      if (func_vals is None) or (coords is None):
         raise ValueError('Please pass in either potential or func_vals and coords')
      else:
         potential   = dirichlet_fund_soln(coords,func_vals)

   stream   = dirichlet_fund_soln(1*potential.coords,\
      1*potential.stream_func_bdy,singularities=1*potential.singularities)

   return stream
#######################################################

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
