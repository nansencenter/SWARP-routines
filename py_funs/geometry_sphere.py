#######################################################
def rsphere_authalic(ellipsoid_mat):
   # translation of rsphere.m
   # 'authalic' case (equal surface area sphere)
   # returns radius of the sphere with same radius
   import numpy as np

   a = ellipsoid_mat[0]
   e = ellipsoid_mat[1]
   if e > 0.:
      f1 = a*a/2.
      f2 = (1 - e*e) / (2.*e)
      f3 = np.log((1.+e) / (1.-e))
      r  = np.sqrt(f1 * (1. + f2 * f3));
   else:
      r = a
   
   return r
#######################################################

#######################################################
def convertlat_geodetic2authalic(e,latin):
   # latout = convertlat_geodetic2authalic(e,latin)
   # convert latitude to equivalent latitude on authalic sphere
   # (equal area sphere)
   # translation of 'geodetic2authalic' option of convertlat.m
   import numpy as np

   if e>.5:
      print('warning (geodetic2authalic): e>.5')
      print('Auxiliary sphere approximation weakens with high eccentricity')
      print(e)

   fact1 = pow(e,2)/3. + 31.*pow(e,4)/180. + 59.*pow(e,6)/560.
   fact2 = 17.*pow(e,4)/360. + 61.*pow(e,6)/1260.
   fact3 = 383.*pow(e,6)/45360.

   #  Truncated series expansion.
   latout = latin - \
            fact1*np.sin(2.*latin) + \
            fact2*np.sin(4.*latin) - \
            fact3*np.sin(6.*latin)

   return latout
#######################################################

#######################################################
def convertlat_geodetic2conformal(e,phi,inverse=False,radians=True):
   # chi=  convertlat_geodetic2conformal(ecc,phi):
   # phi=latitude,e=eccentricity,
   # chi=conformal latitude of Snyder (1982), eqn (3-1,3-2), p18
   import numpy as np

   if not radians:
      # to radians
      phi   = phi*np.pi/180.

   if not inverse:
      fact1 = 1. - e*np.sin(phi)
      fact2 = 1. + e*np.sin(phi)
      chi   = 2.*np.arctan(np.tan(np.pi/4. + phi/2.)*pow(fact1/fact2,e/2.)) - np.pi/2.
   else:
      fact1 = pow(e,2)/2. + 5.*pow(e,4)/24. + pow(e,6)/12.
      fact2 = 7.*pow(e,4)/48. + 29.*pow(e,6)/240.
      fact3 = 7.*pow(e,6)/120.
      latin = phi # actually chi
      chi   = latin + \
            fact1*np.sin(2.*latin) + \
            fact2*np.sin(4.*latin) + \
            fact3*np.sin(6.*latin) # phi

   if not radians:
      # back to degrees
      chi   = chi*180./np.pi

   return chi
#######################################################

#######################################################
def greatcircledist(lat1, lon1, lat2, lon2, R=None,radians=True):

   import numpy as np

   # Calculate great circle distance between points on a sphere using the
   # Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
   # length and has the same units as the radius of the sphere, R.  (If R is
   # 1, then RNG is effectively arc length in radians.)

   R0 = 6378273.# default is hyc2proj radius (m)
   if R is None:
      R  = R0
   else:
      R  = float(R)

   if not radians:
      lat1  = lat1*np.pi/180.
      lat2  = lat2*np.pi/180.
      lon1  = lon1*np.pi/180.
      lon2  = lon2*np.pi/180.

   # Haversine formula
   a     = pow( np.sin((lat2-lat1)/2.),2) + np.cos(lat1)*np.cos(lat2)*pow(np.sin((lon2-lon1)/2.),2)
   rng   = R*2.*np.arctan2(np.sqrt(a),np.sqrt(1 - a))

   # more accurate formula for close points
   if type(rng)==type(np.zeros(1)):
      # arrays
      nvec  = np.arange(len(rng))
      for n in nvec[(rng/R)<(5.e3/R0)]:
         dphi  = lat2[n]-lat1[n]
         dlam  = lon2[n]-lon1[n]
         c     = pow(np.sin(dphi/2.),2)+\
                  np.cos(lat1[n])*np.cos(lat2[n])*pow(np.sin(dlam/2.),2)
         #
         rng[n]  = R*2*np.arcsin(np.sqrt(c))
   else:
      # numbers
      if (rng/R)<(5.e3/R0):
         dphi  = lat2-lat1
         dlam  = lon2-lon1
         rng   = R*2*np.arcsin(np.sqrt(\
            pow(np.sin(dphi/2),2)+\
            np.cos(lat1)*np.cos(lat2)*pow(np.sin(dlam/2),2)\
               ))
   
   return rng
#######################################################

#######################################################
def arc_length(lons,lats,R=None,radians=False,closed=False):
   # given vectors or lists lons,lats
   # returns a numpy array of the same size,
   # ranging from 0 to P,
   # where P is the perimeter

   import numpy as np
   Nl       = len(lons)
   arc_len  = np.zeros(Nl)

   if closed and (lons[0]!=lons[-1] or lats[0]!=lats[-1]):
      # repeat last point
      lons  = list(lons)
      lats  = list(lats)
      lons.append(lons[0])
      lats.append(lats[0])

   for n in range(1,Nl):
      lon0        = lons[n-1]
      lat0        = lats[n-1]
      lon1        = lons[n]
      lat1        = lats[n]
      arc_len[n]  = arc_len[n-1]+greatcircledist(lat1, lon1, lat0, lon0, R=R,radians=radians)
   
   return arc_len
#######################################################

#######################################################
def perimeter(lons,lats,R=None,radians=False,closed=False):
   # given vectors or lists lons,lats
   # returns the perimeter

   import numpy as np

   if closed and (lons[0]!=lons[-1] or lats[0]!=lats[-1]):
      # repeat last point
      lons  = list(lons)
      lats  = list(lats)
      lons.append(lons[0])
      lats.append(lats[0])

   P  = 0
   for n in range(1,len(lons)):
      lon0  = lons[n-1]
      lat0  = lats[n-1]
      lon1  = lons[n]
      lat1  = lats[n]
      P     = P+greatcircledist(lat1, lon1, lat0, lon0, R=R,radians=radians)
   
   return P
#######################################################

#######################################################
def area_polygon_ellipsoid(lon,lat,ellipsoid_mat=None,ellipsoid=None,radians=False):
   # area of polygon on ellipsoid
   # translation of areaint.m
   import numpy as np

   # default ellipsoid
   if ellipsoid_mat is None:
      if ellipsoid is None:
         # ellipsoid=[a,b,flattening], a>b, flattening=(b-a)/a
         # same as hyc2proj (spherical earth)
         ellipsoid   = [6378273,6378273,0]

      a              = ellipsoid[0] # semi-major axis
      b              = ellipsoid[1]
      ecc            = np.sqrt(1-pow(b/a,2)) # eccentricity
      ellipsoid_mat  = [a,ecc]      # ellipsoid parameters used by matlab

   if not radians:
      # convert lons,lats to radians
      lon   = lon*np.pi/180.
      lat   = lat*np.pi/180.

   #########################################################################
   # define some auxiliary functions:
   #########################################################################

   #######################################################
   def greatcircleaz(lat1,lon1,lat2,lon2):

      import numpy as np
      # Inputs LAT1, LON1, LAT2, LON2 are in units of radians.

      y  = np.cos(lat2)*np.sin(lon2-lon1) #cos(lat2) .* sin(lon2-lon1)
      x  =   np.cos(lat1)*np.sin(lat2)\
           - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
           # cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1)

      az = np.arctan2(y,x)

      # Azimuths are undefined at the poles, so we choose a convention: zero at
      # the north pole and pi at the south pole.
      az[lat1 <= -np.pi/2] = 0.
      az[lat2 >=  np.pi/2] = 0.
      az[lat2 <= -np.pi/2] = np.pi
      az[lat1 >=  np.pi/2] = np.pi

      return az
   #######################################################

   #######################################################
   def distance(lat1, lon1, lat2, lon2) :

      ellipsoid_mat  = [1.,0.]

      # Start with spherical approximation even if using an ellipsoid
      rng = greatcircledist(lat1, lon1, lat2, lon2, ellipsoid_mat[0])

      # Azimuth on a sphere
      az = greatcircleaz(lat1, lon1, lat2, lon2)
      return rng,az
   #######################################################

   #######################################################
   def singleint(lon,lat):
      import numpy as np
      # Compute the area of a single polygon on a unit sphere:

      # Ensure the path segment closes upon itself by
      # repeating beginning point at the end.
      if (lat[-1]!=lat[0]) or (lon[-1]!=lon[0]):
         lat   = np.concatenate([lat,[lat[0]]],0)
         lon   = np.concatenate([lon,[lon[0]]],0)

      # Set origin for integration.  Any point will do, so (0,0) used
      lat0 = 0*lat
      lon0 = 0*lon

      # Get colatitude (a measure of surface distance as an angle)
      # and azimuth of each point in segment from the arbitrary origin
      colat,az = distance(lat0,lon0,lat,lon)
      # NB we have already mapped to a sphere with convertlat
      # [colat, az] = distance('gc',lat0,lon0,lat,lon,'radians');
      # gc -> use_geodesic=1,ellipsoid=[1,0] (unit sphere)

      # Calculate step sizes
      daz   = az[1:]-az[:-1]

      # keep daz in  [-pi,pi] interval
      daz[daz<-np.pi]  = daz[daz<-np.pi]+np.pi
      daz[daz>np.pi]   = daz[daz>np.pi]-np.pi

      # Determine average surface distance for each step
      deltas   = (colat[1:]-colat[:-1])/2.
      colat    = colat[:-1]+deltas

      # Integral over azimuth is 1-cos(colatitudes)
      integrands  = (1-np.cos(colat))*daz

      # Integrate and save the answer as a fraction of the unit sphere.
      # Note that the sum of the integrands will include a factor of 4\pi.
      area = abs(np.sum(integrands))/(4.*np.pi) # Could be area of inside or outside

      # Define the inside to be the side with less area.
      return np.min([area,1-area])
   #######################################################

   ###################################################################
   # start to do calculations
   ###################################################################

   # radius of sphere with same surface area
   radius = rsphere_authalic(ellipsoid_mat)

   # convert latitude to equivalent latitude on authalic sphere
   # lat on sphere with same surface area
   # NB this is not the "conformal lat" in Snyder (1982)
   # use convertlat_geodetic2conformal for this
   ecc   = ellipsoid_mat[1]
   lat   = convertlat_geodetic2authalic(ecc, lat)
   
   A  = singleint(lon,lat)
   A  = A * 4.*np.pi*pow(radius,2)
   return A
###################################################################

###################################################################
def lonlat2vec(lon,lat):
   import numpy as np
   z  = a*np.sin(lat)
   r  = b*np.cos(lat)
   x  = rc*np.cos(lon)
   y  = rc*np.sin(lon)
   v  = np.array([x,y,z])
   return v
###################################################################

###################################################################
def polar_stereographic_simple(lon_vec,lat_vec,NH=True,radius=6371.e3,inverse=False):
   """
   x,y      = polar_stereographic_simple(lon,lat,NH=True,radius=6371.e3,inverse=False)
   lon,lat  = polar_stereographic_simple(x,y,NH=True,radius=6371.e3,inverse=True)
   *x,y,lon,lat are numpy arrays
   *simple stereographic projection - when not worrying about area distortion etc,
   only conformity of mapping
   *NH=True: northern hemisphere - project from south pole;
   else: southern hemisphere - project from north pole
   *radius (of earth)
   """

   import numpy as np

   if type(lon_vec)!=type(np.array([0])):
      raise ValueError('lon_vec should be a numpy array')
   if type(lat_vec)!=type(np.array([0])):
      raise ValueError('lat_vec should be a numpy array')

   ###########################################
   def projection_cartesian(lon,lat,radius):

      a  = radius
      b  = radius
      z  = a*np.sin(np.pi/180.*lat)
      r  = b*np.cos(np.pi/180.*lat) # tan(lat)=(z*b)/(a*r)
      x  = r*np.cos(np.pi/180.*lon)
      y  = r*np.sin(np.pi/180.*lon) # tan(lon)=y/x

      return x,y,z
   ###########################################


   ###########################################
   if not inverse:

      x,y,z = projection_cartesian(lon_vec,lat_vec,radius)
      r     = np.sqrt(x**2+y**2)

      if NH:
         # project from South Pole
         r2 = r*radius/(z+radius)
         x2 = (r2/r)*x
         y2 = (r2/r)*y
      else:
         # project from North Pole
         r2 = r*radius/(-z+radius)
         x2 = (r2/r)*x
         y2 = (r2/r)*y

      return x2,y2

   ###########################################
   else:

      x2,y2 = lon_vec,lat_vec
      r2    = np.sqrt(x2**2+y2**2)

      if NH:
         # project from South Pole
         A  = radius**2+r2**2
         B  = 2*radius*r2**2
         C  = -radius**2*(radius**2-r2**2)
         z  = (-B+np.sqrt(B**2-4*A*C))/(2*A) # >0
         r  = r2/radius*(z+radius)
         x  = r/r2*x2
         y  = r/r2*y2
      else:
         # project from North Pole
         A  = radius**2+r2**2
         B  = -2*radius*r2**2
         C  = -radius**2*(radius**2-r2**2)
         z  = (-B-np.sqrt(B**2-4*A*C))/(2*A) # <0
         r  = r2/radius*(-z+radius)
         x  = r/r2*x2
         y  = r/r2*y2

      lon   = (180./np.pi)*np.arctan2(y,x)
      lat   = (180./np.pi)*np.arctan2(z,r)

      return lon,lat
   ###########################################
