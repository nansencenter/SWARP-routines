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
      
      import numpy   as np
      
      #######################################################
      # check enough points, anticlockwise rotation, evenly spaced,...
      # - also define area, perimeter...
      # - define rtree index
      self._check_coords(coords,func_vals,do_check=(singularities is None))
      #######################################################
      
      #######################################################
      # GET SINGULARITIES
      self.buffer_resolution	= 8 # default buffer resolution
         # (number of segments used to approximate a quarter circle around a point.)
      
      if singularities is None:
         self._get_singularities()
      else:
         self.singularities            = singularities 
         self.number_of_singularities  = len(self.singularities)
         print('Number of boundary points : '+str(self.number_of_points))
         print('Number of singularities   : '+str(self.number_of_singularities)+'\n')
      #######################################################

      #######################################################
      # SOLVE LAPLACE'S EQN (get a_n)
      self._solve_laplace_eqn()

      # # evaluate error on boundary:
      print('\nCalculating error on the boundary...')
      self._eval_solution_boundary()
      print(str(self.boundary_error)+'\n')
      #######################################################

      #######################################################
      if 0:
         # evaluate normal derivative -> stream function on boundary:
         print('\nCalculating normal derivative -> stream function on the boundary...\n')
         self._eval_derivs_boundary()

      elif 1:
         # use analytical definition of stream function
         # - for log, this is arctan2 (+const)
         self._get_stream_func_bdy()
      #######################################################

      return
   #######################################################
   
   #######################################################
   # FUNCTIONS FOR INITIALISATION OF THE OBJECT BELOW:
   #######################################################

   #############################################################
   def _check_coords(self,coords,fvals,do_check=True):

      import rtree.index	as Rindex
      import numpy            as np
      import MIZchar          as mizc
      import geometry_planar  as GP
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html

      #############################################################
      def do_fill(pcoords,res=None,f=None):

         #############################################################
         if res is None: 
            print('Too few points - adding more')
         else:
            print('Checking boundary points are close enough together')

         xc,yc    = np.array(pcoords).transpose()
         xc2,yc2  = mizc.fill_poly(xc,yc,res=res)
         pc       = np.array([xc2,yc2]).transpose()   # arrays of arrays: coords = rows
         pc       = [tuple(xy) for xy in pc]          # list of tuples
         #############################################################

         #############################################################
         if f is not None:
            fv = []
            i0 = -1
            for i,cc in enumerate(pc):
               if cc in pcoords:
                  # NB this includes 1st and last points
                  i0 = i0+1
                  fv.append(f[i0])
               else:
                  x,y   = cc
                  x0,y0 = pcoords[i0]
                  x1,y1 = pcoords[i0+1]
                  wt    = np.sqrt(pow(x-x0,2)+pow(y-y0,2))/np.sqrt(pow(x1-x0,2)+pow(y1-y0,2))
                  #
                  f0 = f[i0]
                  f1 = f[i0+1]
                  fv.append(f0+wt*(f1-f0))
         #############################################################

         return pc,fv
      #############################################################

      pcoords  = list(coords)
      fvals    = list(fvals)
      pcoords  = [tuple(cc) for cc in pcoords]

      #######################################################
      # check for periodicity
      if (pcoords[0]!=pcoords[-1]):
         # not periodic
         # - BUT need last and first coordinate the same for shapely polygon
         print('not periodic')
         pcoords.append(pcoords[0])
         fvals.append(fvals[0])
      #######################################################

      if do_check:

         #######################################################
         # want to go round curve anti-clockwise
         x,y         = np.array(pcoords).transpose()
         area        = GP.area_polygon_euclidean(x,y)
         self.area   = abs(area)
         if area<0:
            print("Curve traversed in clockwise direction - reversing arrays' order")
            pcoords.reverse()
            fvals.reverse()
         #######################################################

         #############################################################
         #check that there are enough points 
         Nmin  = 40
         while len(pcoords)<Nmin:
            #double no of points
            pcoords,fvals  = do_fill(pcoords,f=fvals)
         #############################################################

         #############################################################
         #check the points are evenly spaced
         spc            = GP.curve_info(pcoords)[2]
         res            = np.mean(spc)
         pcoords,fvals  = do_fill(pcoords,res=res,f=fvals)
         #######################################################
         

      #######################################################
      # make shapely polygon
      # (coords now has end-point repeated,
      #	self.coords doesn't)
      self.shapely_polygon = shgeom.Polygon(pcoords)
      self.coords          = pcoords[:-1]
      self.func_vals       = np.array(fvals[:-1])
      #######################################################


      #######################################################
      # gets spacings,directions between points, and perimeter
      self.number_of_points   = len(self.coords)
      self.perimeter,self.resolution,\
      	 self.spacings,self.tangent_dirn  = GP.curve_info(self.coords)
      #######################################################

      #######################################################
      # make rtree index
      idx   = Rindex.Index()
      for i,(xp,yp) in enumerate(self.coords):
         idx.insert(i,(xp,yp,xp,yp)) # a point is a rectangle of zero side-length
      self.coord_index	= idx
      #######################################################

      return
   #######################################################
   
   #######################################################
   def _get_stream_func_bdy(self):

      # calculate stream function on boundary
      # - continuous on boundary (also should be periodic)

      import numpy            as np
      import geometry_planar  as GP
      
      # coordinates of polygon
      x,y   = np.array(self.coords).transpose()
      Nx    = len(x)
      nvec  = np.arange(Nx)
      sfun  = 0*x
      
      # for each singularity (just outside polygon)
      perimeter,resolution,\
         spacings,tangent_dirn	= GP.curve_info(self.singularities)

      for i,sing in enumerate(self.singularities):
         an	         = self.sing_coeffs[i]
         branch_dir  = np.pi/2.+tangent_dirn[i]
         atan2	   = GP.arctan2_branch(y,x=x,branch_point=sing,branch_dir=branch_dir)
         atan2	   = GP.make_arctan2_cts(atan2)
         sfun        = sfun+an*atan2
      
      self.stream_func_bdy = sfun
      return
   #######################################################

   #######################################################
   def _get_singularities(self):
      import MIZchar as mizc
      import numpy as np

      print('Getting singularities...\n')

      bufres   = self.buffer_resolution
      eps      = self.resolution/2.
      poly     = self.shapely_polygon.buffer(eps,resolution=bufres)

      # make sure sing's are evenly spaced:
      xs,ys = np.array(poly.exterior.coords).transpose()
      xs,ys = mizc.fill_poly(xs,ys,1.5*eps)
      xys   = np.array([xs,ys]).transpose()
      #
      self.singularities            = [tuple(xyi) for xyi in xys]
      self.number_of_singularities  = len(xys)

      # put singularities (z_n) on the boundary of expanded polygon
      self.singularities   = list(poly.exterior.coords)

      if self.singularities[0]==self.singularities[-1]:
         self.singularities = self.singularities[:-1]

      self.number_of_singularities  = len(self.singularities)
      
      do_check = True
      while do_check:
         # check the number of singularities:
         do_check = self._check_singularities()
      # else:
      #    # accept automatically
      #    do_check = False
      #    print('Warning: not checking singularities')
      #    print('Number of boundary points : '+str(self.number_of_points))
      #    print('Number of singularities	: '+str(self.number_of_singularities)+'\n')
         
      return
   #######################################################
   
   #######################################################
   def _check_singularities(self):
      # don't want too many singularities
      import numpy as np
      import geometry_planar	as GP
      
      N0 = self.number_of_points
      N1 = self.number_of_singularities
      
      # set limits for N1
      Nthresh  = 100
      frac     = 1.3
      # frac   = .2
      
      # if self.solve_exactly:
      #		# number of singularities and number of boundary points should be the same
      #		# - this doesn't work too well - need more sing's
      #		Ntarget	= N0
      if N0 <= Nthresh:
         # if N0<=Nthresh, try to get N1~N0
         Ntarget  = N0+30
      else:
         # try to get N1~frac*N0, if frac*N0<Nthresh
         Ntarget	= int(np.max([np.round(frac*N0),Nthresh]))
      
      print('\nChecking singularities...\n')
      print('Number of boundary points	     : '+str(N0))
      print('Number of singularities         : '+str(N1))
      print('Desired number of singularities : '+str(Ntarget))
      
      check_again = False
      # NB this applies to the case where N1==Ntarget
      # NB also if N1>=Ntarget, we can reduce N1 to Ntarget exactly in 1 go
      
      if N1>Ntarget:
      
         ##########################################################################
         # reduce the number of sing's by increasing spacing between points
         coords	= self.singularities
         perimeter,resolution,\
               spacings,tangent_dirn   = GP.curve_info(coords)
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
               idel	= np.random.randint(N2)
               new_coords.remove(new_coords[idel])
               N2       = len(new_coords)
               Ndel     = N2-Ntarget
            
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
            Nadd	= Ntarget-N2
            iadd	= 0
            
            for coord in self.singularities:
               if coord in new_coords:
               	  # in there so now need to insert before next element of new_coords
               	  iadd	= iadd+1
               else:
               	  # adds coord before new_coords[iadd]
               	  new_coords  = list(new_coords)
               	  new_coords.insert(iadd,coord)
               	  N2    = len(new_coords)
               	  Nadd  = Ntarget-N2

               	  # new element so still need to increase point at which to place 
               	  iadd	= iadd+1

               if Nadd==0:
                  break

            # update list of singularities
            self.singularities            = new_coords
            self.number_of_singularities  = N2
            #######################################################################

         ##########################################################################
         print('New number of singularities		: '+str(Ntarget)+'\n')
         ##########################################################################

      elif N1<Ntarget:

         ##########################################################################
         # set up for another iteration
         if 0:
            # increase buffer_resolution
            # - not so good since this only adds more points round corners
            fac                     = int(np.ceil(1.2*Ntarget/float(N1)))
            self.buffer_resolution  = fac*self.buffer_resolution
         else:
            import MIZchar as mizc
            xs,ys = np.array(self.singularities).transpose()

            # double the number of sings
            # TODO this may not be the best way
            # - perhaps specify resolution
            xs,ys = mizc.fill_poly(xs,ys)
            xys   = np.array([xs,ys]).transpose()
            #
            self.singularities            = [tuple(xyi) for xyi in xys]
            self.number_of_singularities  = len(xys)
         ##########################################################################

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
      M	= np.zeros((Nc,Ns))
      for m,c0 in enumerate(self.coords): # rows of matrix (eqn's)
         p0 = shgeom.Point(c0)
         for n,c1 in enumerate(self.singularities): # cols of matrix (var's)
            p1       = shgeom.Point(c1)
            R        = p0.distance(p1)
            M[m,n]   = np.log(R)
      
      # solve in a least squares sense:
      # M*a=self.func_vals
      self.sing_coeffs	= np.linalg.lstsq(M,self.func_vals)[0]
      return
   #######################################################
   
   #######################################################
   def _eval_solution_boundary(self):
      import geometry_planar as GP
      
      # get x,y for boundary
      x,y   = GP.coords2xy(self.coords)
      
      # evaluate func on boundary
      self.func_vals_approx   = self.eval_solution_fast(x,y)
      
      # calculate error
      self.boundary_error  = GP.vector_norm(self.func_vals-self.func_vals_approx)\
      			      /GP.vector_norm(self.func_vals)

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
      G	 = np.zeros(len(x))
      Nc = len(x)
      for n in range(1,Nc):
         Gs    = -F_n[n]
         ds    = self.spacings[n]
         G[n]  = G[n-1]+Gs*ds
      
      self.tangent_deriv	= F_s
      self.normal_deriv		= F_n
      self.stream_func		= G
      return
   #######################################################
   
   #######################################################
   # EXTERNAL FUNCTIONS THE CLASS PROVIDES
   #######################################################
   
   #######################################################
   def _maskgrid_outside_polygon(self,x,y):
      # use shapely
      
      import numpy as np
      import shapely.geometry as shgeom
      
      Nx    = x.size
      shp   = x.shape
      x	    = x.reshape(Nx)
      y	    = y.reshape(Nx)
      mask  = np.zeros(Nx)
      
      poly2 = self.shapely_polygon
      for m in range(Nx):
         # loop over points to evaluate F at
         p0 = shgeom.Point((x[m],y[m]))
         
         # check if p0 is inside the domain
         mask[m]  = (poly2.intersects(p0))
      
      return np.array(mask.reshape(shp),dtype=bool)
   #######################################################
   
   #######################################################
   def eval_solution(self,x,y):
      # check if inside the polygon before evaluating solution

      import numpy as np
      import geometry_planar as GP

      Nx		= x.size
      F		= np.zeros(Nx)+np.nan
      shp	= x.shape
      x		= x.reshape(Nx)
      y		= y.reshape(Nx)

      if 0:
         # use shapely to get points in polygon
         mask  = self._maskgrid_outside_polygon(x,y)
      else:
         # use matplotlib to get points in polygon
         mask  = GP.maskgrid_outside_polygon(x,y,1*self.coords)

      if 0:
         nvec  = np.arange(Nx)
         nc    = nvec[mask]
         print(mask)
         print(nc)
         print(Nx,len(nc))

      F[mask]  = self.eval_solution_fast(x[mask],y[mask])
      return F.reshape(shp)
   #######################################################
   
   #######################################################
   def eval_solution_fast(self,x,y):
      # eval solution without checking if inside the polygon
      
      import numpy as np
      
      Nx    = x.size
      F	    = np.zeros(Nx)
      shp   = x.shape
      x	    = x.reshape(Nx)
      y	    = y.reshape(Nx)
      
      for m in range(Nx):
         # loop over points to evaluate F at
         # - add contribution from the singularities in 1 step
         x0,y0 = x[m],y[m]
         xs,ys = np.array(self.singularities).transpose()
         Rsq   = pow(xs-x0,2)+pow(ys-y0,2)
         F[m]  = F[m]+(.5*np.log(Rsq)).dot(self.sing_coeffs)

      return F.reshape(shp)
   #######################################################

   #######################################################
   def eval_gradient(self,x,y):

      import numpy as np
      from shapely.prepared import prep	# want "contains" function
      import shapely.geometry as shgeom

      Nx    = x.size
      shp   = x.shape
      F_x   = np.zeros(Nx)
      F_y   = np.zeros(Nx)
      x	    = x.reshape(Nx)
      y	    = y.reshape(Nx)

      poly2 = prep(self.shapely_polygon) # need to test if (x,y) are inside region
      for m in range(Nx):
         # loop over points to evaluate F_x,F_y at
         p0 = shgeom.Point((x[m],y[m]))

         # check if p0 is inside the domain
         if poly2.intersects(p0):
            # add contribution from each singularity
            for n,c1 in enumerate(self.singularities): # singularities
               dx	= x[m]-c1[0]
               dy	= y[m]-c1[1]
               R2	= dx*dx+dy*dy # R^2
               F_x[m]	= F_x[m]+self.sing_coeffs[n]*dx/R2
               F_y[m]	= F_y[m]+self.sing_coeffs[n]*dy/R2
         else:
            F_x[m]	= np.nan
            F_y[m]	= np.nan

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
   def plot_solution(self,pobj=None,plot_boundary=True,show=True,\
         cbar=True,**kwargs):

      import numpy            as np
      import shapefile_utils  as SFU
      import geometry_planar  as GP
      from matplotlib import cm
      from matplotlib import pyplot as plt
      
      if pobj is None:
         # set a plot object if none exists
         fig = plt.figure()
         ax  = fig.add_subplot(1,1,1)
      else:
         fig,ax   = pobj
      
      if plot_boundary:
         # just plot values of F at boundary
         ss = self.get_arc_length()/1.e3 # km
         ax.plot(ss,self.func_vals_approx,'b',**kwargs)
         ax.plot(ss,self.func_vals,'.k',**kwargs)

         x2,y2  = np.array(self.coords).transpose() # coords can be reversed
         f2     = self.eval_solution_fast(x2,y2)
         ax.plot(ss,f2,'--r',**kwargs)
         xc  = np.round(10.*self.coords[0][0]/1.e3)/10. # km (1dp)
         yc  = np.round(10.*self.coords[0][1]/1.e3)/10. # km (1dp)
         ax.set_xlabel('arc length (km) from '+str((xc,yc))+' (km)')
         ax.set_ylabel('values of target function')

      else:
         # plot F everywhere
         bbox	= self.shapely_polygon.bounds
         eps	= self.resolution/2.

         # get a grid to plot F on:
         x0  = bbox[0]
         y0  = bbox[1]
         x1  = bbox[2]
         y1  = bbox[3]
         x   = np.arange(x0,x1+eps,eps)
         y   = np.arange(y0,y1+eps,eps)
         xp  = np.arange(x0-.5*eps,x1+1.5*eps,eps)/1.e3 #km
         yp  = np.arange(y0-.5*eps,y1+1.5*eps,eps)/1.e3

         # make pcolor/contour plot
         nlevels  = 10
         vmin	  = self.func_vals.min()
         vmax	  = self.func_vals.max()
         dv	  = (vmax-vmin)/(nlevels)
         vlev	  = np.arange(vmin,vmax+dv,dv)
         #
         X,Y	= np.meshgrid(x,y)
         F	= self.eval_solution(X,Y)
         #
         cmap	= cm.jet
         cmap.set_bad(color='w')
         Fm  = np.ma.array(F,mask=np.isnan(F))

         ##################################################################
         if 'basemap' not in str(type(ax)):
            # plotter is not a basemap (eg a normal axis object)
            PC  = ax.pcolor(xp,yp,Fm,vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)
            if cbar:
               fig.colorbar(PC)
            ax.contour(X/1.e3,Y/1.e3,F,vlev,colors='k',**kwargs)

            # plot polygon boundary
            x,y   = GP.coords2xy(self.coords)
            ax.plot(x/1.e3,y/1.e3,'k',linewidth=2,**kwargs)

            # plot singularities
            x,y   = GP.coords2xy(self.singularities)
            ax.plot(x/1.e3,y/1.e3,'.k',markersize=5)

            # axes labels
            ax.set_xlabel('x, km')
            ax.set_ylabel('y, km')
         else:
            # plotter is a basemap

            # NB need units in m for basemap
            PC = ax.pcolor(xp*1e3,yp*1e3,Fm,vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)
            if cbar:
               fig.colorbar(PC)

            # plot outline
            x,y   = GP.coords2xy(self.coords)
            ax.plot(x,y,'k',linewidth=2,**kwargs)

            # plot singularities
            x,y   = GP.coords2xy(self.singularities)
            ax.plot(x,y,'.k',markersize=5,**kwargs)
         ##################################################################

      if show:
         plt.show(fig)

      return
   #######################################################
   
   #######################################################
   def plot_stream_func_bdy(self,pobj=None,show=True,i_test=None):
      import numpy            as np
      import shapefile_utils  as SFU
      import geometry_planar  as GP
      
      if pobj is None:
         # set a plot object if none exists
         from matplotlib import pyplot as pobj
      
      ss = self.get_arc_length()/1.e3 # km

      if i_test is None:
         # just plot values of stresm function at boundary
         pobj.plot(ss,self.stream_func_bdy,'b')
         out   = None
         xc    = np.round(10.*self.coords[0][0]/1.e3)/10. # km (1dp)
         yc    = np.round(10.*self.coords[0][1]/1.e3)/10. # km (1dp)
         pobj.xlabel('arc length (km) from '+str((xc,yc))+' (km)')
         pobj.ylabel('values of stream function')

      else:

      	# coordinates of polygon
      	x,y    = np.array(self.coords).transpose()
      	Nx     = len(x)
      	nvec   = np.arange(Nx)

      	# for each singularity (just outside polygon)
      	sing         = self.singularities[i_test]
      	branch_dir   = np.pi/2.+self.tangent_dirn[i_test]
      	atan2	     = GP.arctan2_branch(y,x=x,branch_point=sing,branch_dir=branch_dir)
      	pobj.plot(ss,atan2/np.pi,'r')
      	out    = [atan2]
      	atan2  = GP.make_arctan2_cts(atan2)
      	pobj.plot(ss,atan2/np.pi,'b')
      	out.append(atan2)
      
      if show:
         pobj.show()

      return out
   #######################################################

   #######################################################
   def get_isolines(self,test_function=None,pobj=None,show=True,func_vals_orig=None):
      from skimage import measure as msr
      import numpy as np
      from matplotlib import cm
      import time

      poly	= self.shapely_polygon
      bbox	= poly.bounds
      eps	= self.resolution/3.

      # get a grid to plot F on:
      x0 = bbox[0]
      y0 = bbox[1]
      x1 = bbox[2]
      y1 = bbox[3]
      #
      xv = np.arange(x0,x1+eps,eps)
      yv = np.arange(y0,y1+eps,eps)
      Nx = len(xv)
      Ny = len(yv)
      
      t0 = time.clock()
      # evaluate solution on the grid
      print('evaluating solution on the grid...')
      X,Y   = np.meshgrid(xv,yv)
      F	    = self.eval_solution(X,Y)
      t1    = time.clock()
      print('time taken (mins): '+str((t1-t0)/60.)+'\n')
      
      # get contours
      print('extracting isolines...\n')
      nlevels  = self.number_of_points/4
      vmin     = self.func_vals.min()
      vmax     = self.func_vals.max()
      dv       = (vmax-vmin)/float(nlevels)
      vlev     = np.arange(vmin,vmax+dv,dv)
      print(str(nlevels)+' contours, for isolines between '+\
      	 str(vmin)+' and '+str(vmax))
      #
      contours          = []
      Fmin              = np.nanmin(F)
      F2                = F.copy() 
      F2[np.isnan(F)]   = Fmin-1
      for V in vlev:
         B              = np.zeros(F.shape)
         B[F2>=V]       = 1.
         B[np.isnan(F)] = np.nan
         conts0	        = msr.find_contours(B,.5)  # list of [ivec,jvec] arrays
         conts	        = []                       # list of (xvec,yvec) 
         ##################################################
         #convert conts0->conts, or (i,j)->(x,y)
         for m2,cont in enumerate(conts0):
            ivec  = cont[:,0] # NB these indices correspond to center of pixels
            jvec  = cont[:,1]
            #
            xvec  = xv[0]+jvec*eps
            yvec  = yv[0]+ivec*eps
            keep  = np.logical_and([not bol for bol in np.isnan(xvec)],\
                                    [not bol for bol in np.isnan(yvec)])
            xk    = xvec[keep]
            yk    = yvec[keep]
            Nok   = len(xk)

            ##################################################
            # if contour is bigger than 1 pixel,
            # add it to list, and also include the nearest points
            # on the polygon boundary at the end
            if Nok>1:
               # convert to list
               xy = list(np.array([xk,yk]).transpose())
               xy = [(xx,yy) for xx,yy in xy]
               
               # end points of contour
               c0 = xy[0]
               cl = xy[-1]
               
               # find nearest points on boundary to end points of contour
               i0 = list(self.coord_index.nearest(c0,1))[0]
               il = list(self.coord_index.nearest(cl,1))[0]
               
               if test_function is not None:
               	  # check if contours cross from "0" to "1"
               	  keep	= test_function(i0,il)
               else:
                  keep	= True

               if keep:
                  xy.insert(0,self.coords[i0])
                  xy.append(self.coords[il])
                  conts.append(xy)
            ##################################################

         contours.append(conts)
         ##################################################

      t2 = time.clock()
      print('time taken (s): '+str(t2-t1)+'\n')

      ##################################################
      # merge contours from different vlevels
      merged_levels	= []
      for conts in contours:
         merged_levels.extend(conts)
      ##################################################

      ##################################################
      # if plot object passed in, make a plot
      if pobj is not None:
         from matplotlib import cm
         cmap	= cm.jet
         # cmap2 = cm.coolwarm
         cmap2 = cm.PRGn
         fig,ax = pobj

         print('plotting isolines...\n')
         xb,yb = np.array(poly.boundary.coords).transpose()
         ax.plot(xb,yb,color='k',linewidth=2)

         if func_vals_orig is not None:
            c	= func_vals_orig
         else:
            c	= self.stream_func_bdy
         
         if (len(c)==len(xb)+1):
            c	= c[:len(xb)]
         elif (len(c)==len(xb)-1):
            c	= np.concatenate([c,np.array([c[0]])])

         ax.scatter(xb, yb, marker='o', s=150, linewidths=0.5,\
                     c=c, cmap=cmap2)
         #
         xp = np.arange(x0-eps/2.,x1+1.5*eps,eps)
         yp = np.arange(y0-eps/2.,y1+1.5*eps,eps)
         #
         cmap.set_bad(color='w')
         Fm = np.ma.array(F,mask=np.isnan(F))
         PC = ax.pcolor(xp,yp,Fm,vmin=vmin,vmax=vmax,cmap=cmap)
         fig.colorbar(PC)
         
         for cont in merged_levels:
            xx,yy = np.array(cont).transpose()
            ax.plot(xx,yy)
         
         if show:
            fig.show()
      ##################################################

      return merged_levels
   #######################################################

   #######################################################
   def get_contour_lengths(self,bmap=None,pobj=None,show=True,test_function=None,\
                              func_vals_orig=None):

      import numpy as np
      ###########################################################################################
      class area_info:

         ########################################################################################
         def __init__(self,xy_bdy_coords,xy_conts,\
                        area,perimeter,lengths,\
                        spherical_geometry=False,\
                        ll_bdy_coords=None,ll_contours=None,\
                        func_vals=None,stream_func=None):
         
            import MIZchar as mizc
            # NB use "1*" to remove pointers to the arrays outside the function
            # - like copy, but works for lists also
            self.xy_bdy_coords   = 1*xy_bdy_coords # (x,y) coordinates of boundary (ie in projected space)
            self.xy_contours     = 1*xy_conts      # (x,y) coordinates of each contour
            self.lengths         = 1*lengths	   # lengths of each contour
            self.area	         = area		   # area of polygon
            self.perimeter       = perimeter	   # perimeter of polygon
            self.FDI             = mizc.frac_dim_index([self.area,self.perimeter])
            
            # lon-lat info if present
            self.spherical_geometry = spherical_geometry
            if spherical_geometry:
               if ll_contours is None:
                  raise ValueError('"ll_contours" not given')
               if ll_bdy_coords is None:
                  raise ValueError('"ll_bdy_coords" not given')
               self.lonlat_contours = 1*ll_contours   # (lon,lat) coordinates of each contour
               self.ll_bdy_coords   = 1*ll_bdy_coords # (lon,lat) coordinates of boundary

            if func_vals is not None:
               self.func_vals	 = 1*func_vals	   # value of function used by Laplace's equation
            if stream_func is not None:
               self.stream_func	 = 1*stream_func   # value of stream function produced by Laplace's equation

            # some summarising info about "lengths"
            lens                       = np.array(lengths)
            self.length_mean           = np.mean(lens)
            self.length_median         = np.median(lens)
            self.length_percentile05   = np.percentile(lens,5)
            self.length_percentile95   = np.percentile(lens,95)

            return
      ###########################################################################################

      if bmap is not None:
         # boundary/area information
         # - use routines for sphere
         import geometry_sphere as GS

         print('getting contour lengths on sphere...\n')

         x,y            = np.array(self.coords).transpose()
         lons,lats      = bmap(x,y,inverse=True)
         lst            = list(np.array([lons,lats]).transpose())
         ll_bdy_coords  = [(lo,la) for lo,la in lst]
         area           = GS.area_polygon_ellipsoid(lons,lats,radians=False)
         arclen         = GS.arc_length(lons,lats,radians=False,closed=True)
         perimeter      = arclen[-1]

         # get isolines of function (plot if pobj is not None)
         contours = self.get_isolines(pobj=pobj,show=show,test_function=test_function,\
                                       func_vals_orig=func_vals_orig)

         ll_conts = []
         lengths  = []

         for cont in contours:
            # list of coords (tuples)
            x,y         = np.array(cont).transpose()
            lons,lats	= bmap(x,y,inverse=True)
            arclen      = GS.arc_length(lons,lats,radians=False,closed=False)
            #
            lengths.append(arclen[-1])	#perimeter
            lst	  = list(np.array([lons,lats]).transpose())
            tups  = [(lo,la) for lo,la in lst]
            ll_conts.append(tups)

         # output object with all the info
         AI = area_info(self.coords,contours,\
                        area,perimeter,lengths,\
                        func_vals=self.func_vals,stream_func=self.stream_func_bdy,\
                        ll_bdy_coords=ll_bdy_coords,ll_contours=ll_conts,\
                        spherical_geometry=True)
      else:
         # boundary/area information
         # - just Euclidean routines
         import geometry_planar as GP

         print('getting contour lengths in the plane...\n')

         x,y         = np.array(self.coords).transpose()
         area        = self.shapely_polygon.area
         perimeter   = self.shapely_polygon.length

         # get isolines of function (plot if pobj is not None)
         contours = self.get_isolines(pobj=pobj,show=show,\
                                       test_function=test_function,\
                                       func_vals_orig=func_vals_orig)
         lengths  = []

         for cont in contours:
            # list of coords (tuples)
            P  = GP.curve_info(cont,closed=False)[0]	# perimeter
            lengths.append(P)

         # output object with all the info
         # area_info(xy_bdy_coords,xy_conts,\
         #			area,perimeter,lengths,\
         #			ll_bdy_coords=None,ll_conts=None,\
         #			func_vals=None,stream_func=None,):
         AI = area_info(self.coords,contours,\
                        area,perimeter,lengths,\
                        func_vals=self.func_vals,stream_func=self.stream_func_bdy,\
                        spherical_geometry=False)

      return AI
   #######################################################

##########################################################

#######################################################
def dirichlet_stream_func(coords=None,potential=None,func_vals=None):
   # CALL: stream = dirichlet_stream_func(potential=potential)
   # inputs:
   # > potential: dirichlet_fund_soln object
   #		- solution to Dirichlet problem with boundary values potential.func_vals
   #
   # OR:
   #
   # CALL: stream = dirichlet_stream_func(coords=coords,func_vals=None):
   # inputs:
   # > coords: list of points (polygon boundary)
   # > func_vals: value of function at each point in coords
   # use these to calculate "potential" object
   #		- solution to Dirichlet problem with boundary values func_vals
   # 
   # outputs:
   # > stream: dirichlet_fund_soln object
   #		-	solution to Dirichlet problem with boundary values potential.stream_func_bdy
   
   if potential is None:
      if (func_vals is None) or (coords is None):
         raise ValueError('Please pass in either potential or func_vals and coords')
      else:
         potential   = dirichlet_fund_soln(1*coords,1*func_vals)
   
   stream   = dirichlet_fund_soln(1*potential.coords,\
   	1*potential.stream_func_bdy,singularities=1*potential.singularities)
   
   return stream
#######################################################

#######################################################################
class MIZ_soln:

   #######################################################
   def __init__(self,AI,potential,stream):
      self.area_info = AI
      self.potential = potential
      self.stream    = stream

      # record for shapefile
      self.record = {}
      self.record.update({'Area'                      : AI.area})
      self.record.update({'Perimeter'                 : AI.perimeter})
      self.record.update({'Fractal_dimension_index'   : AI.FDI})
      self.record.update({'Width_mean'                : AI.length_mean})
      self.record.update({'Width_median'              : AI.length_median})
      self.record.update({'Width_percentile05'        : AI.length_percentile05})
      self.record.update({'Width_percentile95'        : AI.length_percentile95})
      return
   #######################################################

   #######################################################
   def parts(self):
      return [self.area_info.ll_bdy_coords]
   #######################################################

   #######################################################
   def plot_soln(self,option='Potential',pobj=None,show=False,bmap=None,\
         cbar=False,title=False):

      import numpy as np

      if pobj is None:
         from matplotlib import pyplot as plt
         fig   = plt.figure()
         ax    = fig.add_subplot(111)
      else:
         fig,ax   = pobj

      #################################################
      if option=='Potential':
         # plot potential
         plot_fun = self.potential.plot_solution
      elif option=='Stream':
         # plot stream function
         plot_fun = self.stream.plot_solution
      #################################################

      #################################################
      if bmap is not None:
         plot_fun(plot_boundary=False,pobj=[fig,bmap],show=False,
               cbar=cbar,ax=ax)

         # plot streamlines (m)
         for cc in self.area_info.xy_contours:
            xx,yy = np.array(cc).transpose()
            bmap.plot(xx,yy,'m',linewidth=2,ax=ax)

      else:
         plot_fun(plot_boundary=False,pobj=[fig,ax],show=False,cbar=cbar)

         # plot streamlines (km)
         for cc in self.area_info.xy_contours:
            xx,yy = np.array(cc).transpose()
            ax.plot(xx/1.e3,yy/1.e3,'m',linewidth=2)

         ax.set_aspect('equal')
      #################################################

      if title:
         ss = str(self.area_info.length_median/1.e3)
         ax.title.set_text('median contour length (km): '+ss)

      if show:
         fig.show()
      return pobj
   ####################################################################
#######################################################################

##################################################
class multipole:
   # make a test function that satisfies Laplace's equation
   # - of form r^n*cos(n\theta) or r^n*sin(n\theta)
   def __init__(self,radius,center,even=True):
      self.radius = radius
      self.center = center
      self.even	  = even
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
      
      a	= self.radius # scale r by this:
      # a   = 1. # scale r by this:
      if self.even:
         for m,order in enumerate(orders):
            out	= out+coeffs[m]*pow(r/a,order)*np.cos(order*th)
      else:
         for m,order in enumerate(orders):
            out	= out+coeffs[m]*pow(r/a,order)*np.sin(order*th)

      return out.reshape(shp)
##################################################

##################################################
def get_MIZ_widths(lons,lats,fvals=None,name=None,fig_outdir=None,basemap=None,xy_coords2=None):

   # Apparently the modules have to been called again 
   import time
   import sys,os
   import numpy as np 
   from matplotlib import pyplot as plt
   
   sys.path.append('../py_funs')
   import f_vals_smoother as smt	
   
   if xy_coords2 is None:

      ######################################################################################
      # make basemap if needed
      if basemap is None:
         from mpl_toolkits.basemap import Basemap
         import geometry_sphere as GS # from py_funs
         #
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
         ###################################################################################

      x,y         = basemap(lons,lats)
      xy_coords2  = np.array([x,y]).transpose()
      ######################################################################################

   if fvals is not None:
      fvals2 = smt.smoother(fvals)
      CSopt  = 0
   else:
      import MIZchar as mizc
      # use principal components to set
      # the function values on the boundary
      PCA    = mizc.pca_mapper(xy_coords2)
      fvals2 = PCA.set_func_vals()
      CSopt  = 1

   #print(fvals2)
   #print(xy_coords2)
   t0 = time.clock()
   print('\n**********************************************************************')
   print('Calculating potential...\n')
   fun_sol  = dirichlet_fund_soln(xy_coords2,fvals2)#,bmap=basemap)
   t1       = time.clock()
   print('\nTime to get potential (s): '+str(t1-t0))
   print('**********************************************************************\n')

   print('\n**********************************************************************')
   print('Calculating stream function...\n')
   stream   = dirichlet_stream_func(potential=fun_sol)
   t2       = time.clock()
   print('\nTime to get stream function (s): '+str(t2-t1))
   print('**********************************************************************\n')

   ###########################################################################################
   #add a test function to eliminate contours that don't cross from "0" to "1":
   class contour_selection:

      def __init__(self,func_vals):
      	 self.func_vals = 1*func_vals
         return

      def selector_binary(self,i0,il):
   	 # assume func_vals are a binary function
   	 # keep contours that end on the opposite value
   	 # OR the "unknown" part (between 0,1)
   	 f0 = self.func_vals[i0]
   	 fl = self.func_vals[il]
   	 if f0==1.:
   	    keep  = (fl<1.)
   	 elif f0==0.:
   	    keep  = (fl>0.)
   	 else:
   	    #keep = False
   	    keep  = np.logical_or(fl==0.,fl==1.)
   	 return keep

      def selector_binary_v2(self,i0,il):
         # assume func_vals are a binary function
         # remove contours that end on the same value (or group of values)
         f0 = self.func_vals[i0]
         fl = self.func_vals[il]
         if f0==1.:
            keep  = (fl<1.)
         elif f0==0.:
            keep  = (fl>0.)
         else:
            keep  = ((fl==0.) or (fl==1.))
         return keep

      def selector_opp_sign(self,i0,il):
         # remove contours that end on values with opposite signs
            # - good for sinusoidal type functions
         f0     = self.func_vals[i0]
         fl     = self.func_vals[il]
         keep   = (f0*fl<0) # opposite signs
         return keep
   ###########################################################################################

   CS = contour_selection(fun_sol.func_vals)
   if CSopt==1:
      selector_function   = CS.selector_opp_sign
      fstr	          = '_vOpp'
   else:
      if 1:
         selector_function = CS.selector_binary
         fstr	           = '_v1'
      elif 0:
         selector_function = CS.selector_binary_v2
         fstr  	           = '_v2'
      else:
         selector_function = None
         fstr	           = ''
   
   if fig_outdir is None:
      # do not make a figure
      pobj = None
   else:
      # make a figure
      fig     = plt.figure()
      ax      = fig.add_subplot(111)
      pobj    = [fig,ax]
      outdir  = './outputs/aod/'+str(fig_outdir)
   
      if not os.path.exists(outdir):
         os.mkdir(outdir)
   
      outdir = outdir+'/'+str(name)+'_laplacian'
      if not os.path.exists(outdir):
         os.mkdir(outdir)
   
   print('\n**********************************************************************')
   print('Getting streamlines...\n')
   if 0:
      # euclidean space
      AI = stream.get_contour_lengths(pobj=pobj,show=False,\
         test_function=selector_function,\
         func_vals_orig=1*fun_sol.func_vals)
      if pobj is not None:
         figname  = outdir+'/test_Laplacian_planar'+fstr+'.png'
   else:
      # spherical stuff
      AI = stream.get_contour_lengths(pobj=pobj,bmap=basemap,show=False,\
                                       test_function=selector_function,\
                                       func_vals_orig=1*fun_sol.func_vals)
      if pobj is not None:
      	figname	= outdir+'/test_Laplacian_spherical'+fstr+'.png'
   	
   ttl	= 'Median length (km) '+str(np.round(10.*AI.length_median/1.e3)/10.)
   if pobj is not None:
      fig,ax  = pobj
      ax.title.set_text(ttl)
      ax.set_aspect('equal')

      if 0:
         # show for testing
         fig.show()
      else:
         # make figure
         #print('Saving plot to figure '+figname)
         fig.savefig(figname)
         fig.clf()
   
   t3 = time.clock()
   print('\nTime to get streamlines (mins): '+str((t3-t2)/60.))
   print(ttl)
   print('**********************************************************************\n')
   

   Psoln = MIZ_soln(AI,fun_sol,stream)
   # Psoln = [AI,fun_sol,stream]

   return Psoln
