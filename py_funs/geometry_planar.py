#########################################################
def area_polygon_euclidean(x,y):
   # area of a polygon in Euclidean space
   A  = x[:-1]*y[1:]-x[1:]*y[:-1]

   return .5*np.sum(A)
#########################################################

#########################################################
def vector_norm(v):
   import numpy as np
   return np.sqrt(v.dot(v))
#########################################################

#########################################################
def unit_vector(v):
   norm  = vector_norm(v)
   return v/norm
#########################################################

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
   def __init__(self,coords,func_vals):
      # initialise object

      import numpy as np
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html
      from shapely.prepared import prep

      self.coords    = coords
      self.func_vals = func_vals
      self.calc_perimeter() # gets spacings between points and perimeter

      # get points around boundary of polygon, then
      # expand points by an amount related to the spacings of the coords
      self.shapely_polygon  = prep(shgeom.Polygon(coords))
      eps   = self.resolution/2.
      poly  = self.shapely_polygon.buffer(eps)

      # put singularities (z_n) on the boundary of expanded polygon
      self.singularities   = poly.exterior.boundary.coords

      # solve Laplace's eqn (get a_n)
      self.solve_laplace_eqn()

      return
   #######################################################

   #######################################################
   # FUNCTIONS FOR INITIALISATION OF THE PLANE BELOW:
   #######################################################

   #######################################################
   def calc_perimeter(self):
      import numpy as np
      import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html

      p0       = shgeom.Point(coords[0])
      Nc       = len(coords)
      P        = 0
      spacings = []
      for n in range(1,Nc):
         p1 = shgeom.Point(coords[n])
         ds = p0.distance(p1)
         P  = P+ds
         spacings.append(ds)

      self.spacings     = np.array(spacings)
      self.resolution   = np.mean(self.spacings)
      self.perimeter    = P
      return
   #######################################################

   #######################################################
   def solve_laplace_eqn(self):
      import numpy as np

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
   def eval_solution(self,x,y):

      import numpy as np
      Nx    = x.size
      F     = np.zeros(Nx)+np.nan
      shp   = x.shape
      x     = x.reshape(Nx)
      y     = y.reshape(Nx)

      for m in range(Nx): # points to evaluate F at
         p0 = shgeom.Point((x[m],y[m]))
         if self.shapely_polygon.contains(p0):
            for n,c1 in enumerate(self.singularities): # singularities
               p1    = shgeom.Point(c1)
               R     = p0.distance(p1)
               F[m]  = F[m]+self.sing_coeffs[n]*np.log(R)

      return F.reshape(shp)
   #######################################################

   #######################################################
   def plot_solution(self,pobj=None,plot_boundary=True):
      import numpy as np
      import shapefile_utils as SFU

      if pobj is None:
         # set a plot object if none exists
         from matplotlib import pyplot as pobj
      
      poly  = self.shapely_polygon
      if plot_boundary:
         # just plot values of F at boundary
         ss = np.concatenate([[0],self.spacings])
         pobj.plot(ss,self.func_vals,'.k',markersize=1.5)
         # 
         x,y   = poly.exterior.boundary.coords.xy
         F     = self.eval_function(x,y)
         pobj.plot(ss,self.func_vals,'b')
      else:
         # plot F everywhere
         SFU.plot_polygon(poly,color='k',linewidth=2)
         bbox  = poly.bbox
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
         X,Y      = np.meshgrid(x,y)
         F        = self.eval_solution(X,Y)
         pobj.pcolor(X,Y,F)
         pobj.contour(X,Y,F,nlevels)

      return
   #######################################################

##########################################################

##########################################################
class plane:


   #######################################################
   # INITIALISATION
   #######################################################

   #######################################################
   def __init__(self,norm_vec,const=0):
      # initialise object

      import numpy as np

      # define unit normal vector and constant
      norm  = np.sqrt(norm_vec.dot(norm_vec))
      if norm==0.:
         raise ValueError('normal vector must be non-zero')

      self.normal_vector   = norm_vec/norm
      self.constant        = const/norm

      # calculate basis vectors for plane
      self.basis_vectors   = self._plane_basis()
      return
   #######################################################

   #######################################################
   # FUNCTIONS FOR INITIALISATION OF THE PLANE BELOW:
   #######################################################

   #######################################################
   # def _plane_basis(self,norm_vec,const):
   def _plane_basis(self,bv=None):
      # return a basis for the plane

      import numpy as np
      norm_vec = self.normal_vector
      const    = self.constant

      u  = [np.zeros(3),np.zeros(3)]
      for i2 in range(2,-1,-1):

         ###########################################
         # look for 2 vectors in the plane
         n2 = norm_vec[i2]
         if n2!=0.:

            for i in range(2):
               j        = np.mod(i2-2+i,3)
               u[i][j]  = 1 # the two values of x[j] are free parameters
               cc       = norm_vec[j]/n2
               u[i][i2] = const/n2-cc

            break # have found 2 vectors in the plane
         ###########################################

      # normalise 1st vec:
      v     = unit_vector(u[0]-const*norm_vec)
      B     = [v] # 1st member of basis

      # get unit vector orthog to v and norm_vec
      c     = B[0].dot(u[1])
      v     = unit_vector(u[1]-const*norm_vec-c*B[0])
      B.append(v)

      # vector in plane has expansion: v=const*norm_vec+a0*B[0]+a1*B[1]
      return B 
   #######################################################

   #######################################################
   # FUNCTIONS FOR THE PLANE BELOW:
   #######################################################

   #######################################################
   def basis_coordinates(self,v):
      # return coordinates of vector relative to plane's basis
      import numpy as np

      M  = np.array(self.basis_vectors)
      x  = M.dot(v)
      return x
   #######################################################

   #######################################################
   def basis_expansion(self,a):
      # return absolute coordinates of vector,
      # given coordinates in terms of plane's basis
      import numpy as np

      M  = np.array(self.basis_vectors).transpose()
      v  = M.dot(a) + self.constant*self.normal_vector
      return v
   #######################################################

   #######################################################
   def reset_basis(self,v,abs_coords=True,index=0):
      # reset_basis(self,v,coords):
      # change basis so that v (inside the plane) is a basis vector (BV)
      # *index=0, v=1st BV
      # *index=1, v=2nd BV
      # can give absolute coords of v or expansion in terms of current basis

      import numpy as np
      
      # check if v is in the plane
      nn = self.normal_vector
      cc = self.constant
      MB = np.array(self.basis_vectors).transpose() # vectors in columns

      if abs_coords:
         if nn.dot(v)!=cc:
            raise ValueError('v is not in the plane')
         B  = [unit_vector(v-cc*nn)]
         a0 = self.basis_coordinates(v)
         a0 = unit_vector(a0) # normalise
      else:
         a0 = unit_vector(v) # normalise
         B  = [MB.dot(a0)]

      if index==0:
         # rotate by pi/2
         a1 = np.array([-a0[1],a0[0]]) # (1,0)->(0,1),(0,1)->(-1,0)
         B.append(MB.dot(a1))
      else:
         # rotate by -pi/2
         a1 = np.array([a0[1],-a0[0]]) # (1,0)->(0,-1),(0,1)->(1,0)
         B  = [MB.dot(a1),B[0]]

      # # now normalise the vectors
      # for n in range(2):
      #    norm  = np.sqrt(B[n].dot(B[n]))
      #    print('hi')
      #    print(norm)
      #    B[n]  = B[n]/norm
      #    print(B)

      # update the basis
      self.basis_vectors   = B
      return
   #######################################################

   #######################################################
   def intersection_with_line(self,line_points,abs_coords=False):

      import numpy as np

      v,b      = line_points  # numpy vectors, or (3 x N) numpy arrays: line goes through these
      norm_vec = self.normal_vector
      const    = self.constant

      # line x(t)=a*t+b
      # t=0 -> b, t=1 -> v
      a  = v-b

      ####################################################
      if len(a.shape)>1:
         m  = np.array(norm_vec.dot(a))
         c  = np.array(norm_vec.dot(b))
         t  = (const-c[m!=0])/m[m!=0]
         #
         x        = np.zeros(b.shape)+np.nan
         x[m!=0]  = a[m!=0]*t+b[m!=0] # absolute coordinates of intersection point
      else:
         m  = norm_vec.dot(a)
         c  = norm_vec.dot(b)
         if m!=0:
            t  = (const-c)/m
            x  = a*t+b # absolute coordinates of intersection point
         else:
            x  = np.nan+np.zeros(3)
      ####################################################

      if not abs_coords:
         # return intersection point in coordinates relative to basis vectors
         M  = np.array(self.basis_vectors)
         x  = M.dot(x)

      return x
   #######################################################

   #######################################################
   # END OF CLASS "plane"
   #######################################################

##########################################################