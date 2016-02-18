#######################################################
def maskgrid_outside_polygon(x,y,coords):
   # use matplotlib.path to test if multiple points are contained
   # inside a polygon
   # INPUTS:
   # x,y: numpy arrays with x and y coordinates to test respectively
   # coords: list of coordinates of polygon's boundary (as tuples or lists)
   # [(x0,y0),(x1,y1),...] or [[x0,y0],[x1,y1],...]
   # RETURNS:
   # mask of same shape as x and y
   # > mask element is True/False if corresponding point is inside/outside the polygon

   import numpy as np
   from matplotlib import path
   
   Nx    = x.size
   shp   = x.shape
   x     = x.reshape(Nx)
   y     = y.reshape(Nx)

   bbPath   = path.Path(np.array(coords))

   # points to test [(xi,yi)]
   # coords2  = [(x[i],y[i]) for i in range(Nx*Ny)]
   coords2  = np.array([x,y]).transpose()
   mask     = np.array(bbPath.contains_points(coords2),dtype=bool)

   return mask.reshape(shp)
#######################################################

#########################################################
def area_polygon_euclidean(x,y):
   # area of a polygon in Euclidean space
   # - negative if traveling in a clockwise direction
   # http://www.wikihow.com/Calculate-the-Area-of-a-Polygon
   import numpy as np

   # repeat last point
   if x[0]!=x[-1]:
      x  = list(x)
      x.append(x[0])
      x  = np.array(x)
      #
      y  = list(y)
      y.append(y[0])
      y  = np.array(y)

   A  = x[:-1]*y[1:]-x[1:]*y[:-1]
   return .5*np.sum(A)
#########################################################

#########################################################
def arctan2_branch(y,x=None,branch_point=None,branch_dir=None):
   # returns an angle in interval: linedir+(-pi,pi]

   import numpy as np

   if x is None:
      # y is coord seq
      print(y)
      x,y      = np.array(y).transpose()
      reshape  = 0
   else:
      Nx       = x.size
      shp      = x.shape
      x        = x.reshape((Nx))
      y        = y.reshape((Nx))
      reshape  = 1

   if branch_dir is None:
      branch_dir  = np.pi

   if branch_point is not None:
      x  = x-branch_point[0]
      y  = y-branch_point[1]

   rot   = np.pi-branch_dir # anti-clockwise rotation (towards negative real axis)
   w     = np.log(np.exp(1j*rot)*(x+1j*y))-1j*rot

   if reshape:
      out   = w.imag.reshape(shp)
   else:
      out   = w.imag
   
   return out
#########################################################

#########################################################
def make_arctan2_cts(atan2):
   import numpy as np

   # on a polygon boundary (for example),
   # arctan2 should be continuous (for it to be analytic inside the polygon)
   
   Nx    = len(atan2)
   nvec  = np.arange(Nx)

   #######################################################################
   # check if branch cut intersects polygon (atan2 jumps by 2\pi):
   diff        = 0*atan2
   diff[:-1]   = atan2[1:]-atan2[:-1]
   diff[-1]    = atan2[0] -atan2[-1]
   #
   jj = nvec[diff>1.5*np.pi]
   for j in jj:
      jp1         = np.mod(j+1,Nx)
      atan2[jp1:] = atan2[jp1:]-2*np.pi

   # check if branch cut intersects polygon (atan2 drops by 2\pi):
   diff[:-1]   = atan2[1:]-atan2[:-1]
   diff[-1]    = atan2[0] -atan2[-1]
   #
   jj = nvec[diff<-1.5*np.pi]
   for j in jj:
      jp1         = np.mod(j+1,Nx)
      atan2[jp1:] = atan2[jp1:]+2*np.pi
   #######################################################################

   return atan2
#########################################################

#########################################################
def line_splits_polygon_boundary(poly,isec,shortest=True):

   import shapely.geometry as shgeom
   import numpy as np

   # line should split polygon into
   if 'shapely' in str(type(isec)):
      isec  = list(isec.coords)

   c0 = isec[0]
   c1 = isec[-1]
   p0 = shgeom.Point(c0)
   p1 = shgeom.Point(c1)

   i0    = None
   i1    = None
   lst   = list(poly.boundary.coords)
   N     = len(lst)

   if c0 in poly.boundary.coords: 
      i0    = lst.index(c0)
      i0m1  = lst.index(c0)
   if c1 in poly.boundary.coords: 
      i1    = lst.index(c1)
      i1p1  = i1

   
   if i0 is None:
      # c0 not already present:
      # find the member of list that is closest:
      d0 = p0.distance(poly)
      for i,cc in enumerate(lst):
         ip1   = np.mod(i+1,N)
         lin   = shgeom.LineString([cc,lst[ip1]])

         if p0.distance(lin)==d0:
            i0    = ip1
            i0m1  = i
            break

   if i1 is None:
      # c1 not already present:
      # find the member of list that is closest:
      d1 = p1.distance(poly)
      for i,cc in enumerate(lst):
         ip1   = np.mod(i+1,N)
         lin   = shgeom.LineString([cc,lst[ip1]])

         if p1.distance(lin)==d1:
            i1    = i
            i1p1  = ip1
            break

   if i0<i1:
      clst1 = lst[i0:i1+1]       # go anti-clockwise between points
   else:
      clst1 = lst[i0:]
      clst1.extend(lst[:i1+1])

   if i1p1<i0m1:
      clst2 = lst[i1p1:i0m1+1]
   else:
      clst2 = lst[i1p1:]         # get rest of points (exclude i1, include i1p1)
      clst2.extend(lst[:i0m1+1]) # get rest (include i0m1, not i0)

   # reverse direction so i=0 corresp's to c0
   clst2.reverse()

   L1 = shgeom.LineString(clst1).length
   L2 = shgeom.LineString(clst2).length

   if shortest:
      # choose shortest line
      clst  = clst1
      # if len(clst2)<len(clst):
      if L2<L1:
         clst  = clst2
   else:
      # choose longest line
      clst  = clst1
      # if len(clst2)>len(clst):
      if L2>L1:
         clst  = clst2

   # make sure the initial points are included
   if c0 not in clst:
      clst.insert(0,c0)
   if c1 not in clst:
      clst.append(c1)
   
   return clst
#########################################################

#########################################################
def line_splits_polygon(pol,line_pt1,linedir):

   isec_types,\
      branch_line = line_intersection_polygon(pol,line_pt1,linedir)

   icl         = isec_types.crossing_lines
   Ncl         = len(icl)
   MPlist      = []

   if Ncl>=1:
      pol1,pol2   = line_splits_polygon_in2(pol,icl[0])

      if Ncl==1:

         ################################################################
         # don't need to split further
         # - just sort in decreasing size
         if pol1.length>=pol2.length:
            MPlist   = [pol1,pol2]
         else:
            MPlist   = [pol2,pol1]
         ################################################################

      else:

         ################################################################
         # loop over remaining crossing lines
         for i in range(1,Ncl):
            pol3,pol4   = line_splits_polygon_in2(pol,icl[i])
            # split pol into another 2 poly's
            if not pol3.intersects(pol1):
               # pol3,pol1 are disjoint polygons
               # NB ordering corresponds
               MPlist.extend([pol1,pol3])

               # => pol2,pol4 are big ones
               # => reduce size of big polygon "pol"
               pol   = pol2.intersection(pol4)

            elif not pol4.intersects(pol1):
               # pol4,pol1 are disjoint polygons
               MPlist.extend([pol1,pol4])

               # => pol2,pol3 are big ones
               # => reduce size of big polygon "pol"
               pol   = pol2.intersection(pol3)

            elif not pol4.intersects(pol2):
               # pol4,pol2 are disjoint polygons
               MPlist.extend([pol2,pol4])

               # => pol1,pol3 are big ones
               # => reduce size of big polygon "pol"
               pol   = pol1.intersection(pol3)

            else:
               # pol3,pol2 are disjoint polygons
               MPlist.extend([pol2,pol3])

               # => pol1,pol4 are big ones
               # => reduce size of big polygon "pol"
               pol   = pol1.intersection(pol4)

            # list may have repeated members so reduce it each timestep
            MPlist   = list(set(MPlist))
         ################################################################

         ################################################################
         # loop through again to match ordering of MPlist to icl
         MPlist1  = []
         for lin in icl:
            for poly in MPlist:
               if lin.intersects(poly):
                  MPlist1.append(poly)
                  MPlist.remove(poly)
                  break
         ################################################################

         ################################################################
         # finished looping: pol is last polygon
         # (connecting polygon)
         MPlist1.append(pol)
         ################################################################

   return MPlist1,isec_types,branch_line
#########################################################

#########################################################
def line_splits_polygon_in2(pol,isec):

   import shapely.geometry as shgeom
   import numpy as np

   # line should split polygon into 2
   if 'shapely' not in str(type(isec)):
      isec  = shgeom.LineString(isec)

   bb       = pol.boundary.union(isec)
   MPlist   = []
   MPlist2  = []
   P0       = pol.boundary.coords[0]

   for n,gg in enumerate(bb.geoms):

      print(gg)
      print(gg!=isec)
      if gg!=isec:
         if P0 not in gg.coords:
            # normal polygon
            print('normal polygon')
            MPlist.append(shgeom.Polygon(gg))
         else:
            print('contains '+str(P0))
            # boundary can be artificially split by 1st point also
            # this list should be merged into one polygon later
            MPlist2.extend(list(gg.coords))

   # turn MPlist2 into shapely polygon and append to MPlist
   MPlist.append(shgeom.Polygon(MPlist2))

   return MPlist
#########################################################

#########################################################
def line_intersection_polygon(pol,line_pt1,line_dir):
   import shapely.geometry as shgeom
   import numpy as np

   #####################################################################
   class intersection_types:
      # sort out intersections into points, lines and tangent lines

      def __init__(self,isec,pol):
         self.tangent_points  = []
         self.tangent_lines   = []
         self.crossing_lines  = []

         if isec.geom_type=='Point':
            # intersection is a single point (line just touches curve)
            # print('tangent point')
            # print(ii)
            self.tangent_points.append(isec)

         elif isec.geom_type=='LineString':
            # intersection is a single line
            # print('line')

            if isec.intersection(pol.boundary)==isec:
               # intersection is part of the boundary
               # (doesn't make a new polygon)
               # print('tangent line')
               # print(ii)
               self.tangent_lines.append(isec)

            else:
               # intersection is not part of the boundary
               # (will make a new polygon)
               # print('crossing line')
               # print(ii)
               self.crossing_lines.append(isec)

         else:
            # collection
            # print('collection')
            # print(isec)

            for ii in isec.geoms:
               # print(ii.geom_type)

               if ii.geom_type=='Point':
                  # intersection is a single point (line just touches curve)
                  # print('tangent point')
                  # print(ii)
                  self.tangent_points.append(ii)

               elif ii.geom_type=='LineString':
                  # intersection is a single line
                  # print('line')

                  if ii.intersection(pol.boundary)==ii:
                     # intersection is part of the boundary
                     # (doesn't make a new polygon)
                     # print('tangent line')
                     # print(ii)
                     self.tangent_lines.append(ii)
                  else:
                     # intersection is not part of the boundary
                     # (will make a new polygon)
                     # print('crossing line')
                     # print(ii)
                     self.crossing_lines.append(ii)

         # more information
         self.number_of_crossing_lines = len(self.crossing_lines)
         self.number_of_tangent_lines  = len(self.tangent_lines)
         self.number_of_tangent_points = len(self.tangent_points)
         self.number_of_geometries     = self.number_of_tangent_points+\
                                         self.number_of_tangent_lines+\
                                         self.number_of_crossing_lines

         return
      # finish of __init__
      #####################################################################

   # end of class
   #####################################################################

   if type(line_pt1)!=type((0,0)):
      raise ValueError('line_intersection_polygon: input line_pt1 should be a tuple')
   else:
      x0 = line_pt1[0]
      y0 = line_pt1[1]
      p0 = shgeom.Point(line_pt1)

   # make pol a shapely polygon if it isn't already
   if 'shapely' not in str(type(pol)):
      pol   = shgeom.Polygon(pol)

   P     = pol.length        # perimeter
   D     = p0.distance(pol)  # shortest distance to rectangle

   # make a linestring
   L     = 2*(D+P)            # length of line - if it is in right direction it should intersect
   u     = np.array([np.cos(line_dir),np.sin(line_dir)])
   x1    = x0+u[0]*L
   y1    = y0+u[1]*L
   lin   = shgeom.LineString([line_pt1,(x1,y1)])

   # intersection
   isec        = lin.intersection(pol)
   isec_sorted = intersection_types(isec,pol)

   return isec_sorted,lin
#########################################################

#########################################################
def xy2polar_coords(x,y=None):
   import numpy as np

   use_tups = (y is None)
   if use_tups:
      # x is a tuple list
      xy = np.array(x)
      x  = xy[:,0]
      y  = xy[:,1]

   r  = np.sqrt(x*x+y*y)
   th = np.arctan2(y,x)

   if not use_tups:
      out   = (r,th)
   else:
      # return list of tuples
      # (in same way that input was passed in)
      out   = xy2coords(r,th)

   return out
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
def calc_perimeter(coords,closed=False):
   import numpy as np

   cc    = list(coords)
   x0,y0 = cc[0]
   if closed:
      if cc[-1] is not cc[0]:
         cc.append(cc[0])

   Nc    = len(cc)
   P     = 0 # perimeter

   for n in range(1,Nc):
       x1,y1   = cc[n]
       dx      = x1-x0
       dy      = y1-y0
       #
       ds = np.sqrt(dx*dx+dy*dy)
       P  = P+ds
       #
       x0,y0   = x1,y1

   return P
#########################################################

#########################################################
def arc_length(coords,closed=False):
   import numpy as np

   cc    = list(coords)
   x0,y0 = cc[0]
   if closed:
      if cc[-1] is not cc[0]:
         cc.append(cc[0])

   Nc    = len(cc)
   S     = np.zeros(Nc) # arc length

   for n in range(1,Nc):
       x1,y1   = cc[n]
       dx      = x1-x0
       dy      = y1-y0
       #
       ds   = np.sqrt(dx*dx+dy*dy)
       S[n] = S[n-1]+ds
       #
       x0,y0   = x1,y1

   return S
#########################################################

#########################################################
def curve_info(coords,closed=True):
   import numpy as np

   x0,y0 = coords[0]
   Nc    = len(coords)
   P     = 0   # perimeter

   th_vec   = []  # direction of tangent to curve
   spacings = []  # distance between points
   for n in range(1,Nc):
       x1,y1   = coords[n]
       dx      = x1-x0
       dy      = y1-y0
       #
       ds = np.sqrt(dx*dx+dy*dy)
       th = np.arctan2(dy,dx)
       P  = P+ds
       spacings.append(ds)
       th_vec.append(th)
       x0,y0   = x1,y1

   if closed:
      # last point on curve should be first
      x1,y1   = coords[0]
      dx      = x1-x0
      dy      = y1-y0
      #
      ds = np.sqrt(dx*dx+dy*dy)
      th = np.arctan2(dy,dx)
      if ds>0.:
         # last and first are different,
         # so add them
         P  = P+ds
         spacings.append(ds)
         th_vec.append(th)

   resolution   = np.mean(spacings)
   return P,resolution,np.array(spacings),np.array(th_vec)
#########################################################

#########################################################
def xy2coords(x,y):
   # return list of tuples: [(x[0],y[0]),(x[1],y[1]),...]
   import numpy as np
   v     = np.vstack([x,y]).transpose()
   tup   = [tuple(xy) for xy in v]
   return tup
#########################################################

#########################################################
def coords2xy(coords):
   # convert list of tuples: [(x[0],y[0]),(x[1],y[1]),...]
   # to x,y numpy arrays
   import numpy as np
   x,y   = np.array(coords).transpose()
   return x,y
#########################################################

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
