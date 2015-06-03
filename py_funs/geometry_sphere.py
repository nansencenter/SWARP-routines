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
def greatcircledist(lat1, lon1, lat2, lon2, R):

   import numpy as np

   # Calculate great circle distance between points on a sphere using the
   # Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
   # length and has the same units as the radius of the sphere, R.  (If R is
   # 1, then RNG is effectively arc length in radians.)

   a     = pow( np.sin((lat2-lat1)/2),2) + np.cos(lat1)*np.cos(lat2)*pow(np.sin((lon2-lon1)/2.),2)
   rng   = R*2.*np.arctan2(np.sqrt(a),np.sqrt(1 - a))
   return rng
#######################################################

#######################################################
def area_polygon_ellipsoid(lon,lat,ellipsoid_mat=None,ellipsoid=None,radians=False):
   # translation of areaint.m
   import numpy as np

   # default ellipsoid
   if ellipsoid_mat is None:
      if ellipsoid is None:
         # ellipsoid=[a,b,flattening], a>b, flattening=(b-a)/a
         # same as hyc2proj (spherical earth)
         ellipsoid   = [6378273,6378.273,0]

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

# ###################################################################
# def polar_stereographic_ellipsoid(lon_vec,lat_vec,proj_params,ellipsoid=None,inverse=False):
# 
#    import numpy as np
# 
#    latc,lonc,lat_ts  = proj_params # degrees
#    # latc   = central lat
#    # lonc   = central lat
#    #  - lon=lonc     -> to X=0 in projection
#    #  - (lonc,latc)  -> (0,0)  in projection
#    # lat_ts = latitude of true scale
#    #  - when not using north or south pole as origin,
#    #    the [SHAPE=?] through (lonc,lat_ts) has distances at true scale
# 
#    if ellipsoid is not None:
#       a,b   = ellipsoid 
#       # a      = dist centre to north pole
#       # b      = radius at equator (a>b)
#    else:
#       import geopy.distance as geodist
#       ellipsoid   = geodist.ELLIPSOIDS['WGS-84'] # has extra parameter "flattening=(b-a)/a"
#       a,b         = ellipsoid[:-1]
# 
#    def projection(vp,a,b,lon,lat):
# 
#       z  = a*np.sin(lat)
#       r  = b*np.cos(lat)
#       x  = r*np.cos(lon)
#       y  = r*np.sin(lon)
#       v  = np.zeros((3,len(z)))
# 
#       
# 
#       
#       
# 
#       return X,Y
# 
#    NP    = (latc==90.)
#    SP    = (latc==-90.)
#    TSP   = (lat_ts==latc)
#    if TSP:
#       # dy/dlat=ds/dlat at pole
# 
#    # if NP:
#    #    # north pole is centre
#    #    lat_ref  = 80.
#    #    ang_ref  = -np.pi/2.# (lonc,lat_ref) should be on negative y axis 
#    # elif SP:
#    #    # south pole is centre
#    #    lat_ref  = -80.
#    #    ang_ref  = np.pi/2.# (lonc,lat_ref) should be on positive y axis 
#    # elif (TSP):
#    # else:
#    #    lat_ref  = lat_ts
#    #    ang_ref  = np.sign(lat_ts-latc)*np.pi/2.
#    #       # (lonc,lat_ts) should be on positive y axis if lat_ts>latc
#    #       # (lonc,lat_ts) should be on negative y axis if lat_ts<latc
# 
# 
# 
# 
#    if not inverse:
#       # convert from lon,lat to X,Y
# 
#       # convert to radians if necessary
#       latc     = latc*180./np.pi
#       lonc     = lonc*180./np.pi
#       lat_ts   = lat_ts*180./np.pi
#       lon      = (lon-lonc)*180./np.pi # redefine relative to lonc, so (lon=lonc) -> x=0 automatically
#       lat      = lat*180./np.pi
# 
#       # pole is opposite point to centre of proj in 3d:
#       vp = - lonlat2vec(lonc,latc)
# 
#       # projection plane
#       # plane through origin, normal to vp (vp\cdot x=0):
#       plane_params   = [vp,0]
# 
#       ########################################################################
#       # Get scale factor
# 
#       if NP:
#          # Get rotation and scale factor from lonc,lat_ts
#          v     = lonlat2vec(0,lat_ts)
# 
#          # line through v,vp:
#          line_coeffs = [vp-v,v] # [m,c], x=m*t+c: t=0 -> v, t=1 -> vp
# 
#          # intersection of line, projection plane:
#          X,Y   = intersection_line_plane(plane_params,line_coeffs)
# 
#          # rotation matrix:
#          Rang  = np.sign(lat_ts-latc)*np.pi/2.-np.atan2(Y,X)
#          Rmat  = [[np.cos(Rang),-np.sin(Rang)],
#                   [np.sin(Rang), np.cos(Rang)] ]
# 
# 
#       if lat_ts!=latc:
#          # Get rotation and scale factor from lonc,lat_ts
#          v     = lonlat2vec(lonc,lat_ts)
# 
#          # line through v,vp:
#          line_coeffs = [vp-v,v] # [m,c], x=m*t+c: t=0 -> v, t=1 -> vp
# 
#          # intersection of line, projection plane:
#          X,Y   = intersection_line_plane(plane_params,line_coeffs)
# 
#          # rotation matrix:
#          Rang  = np.sign(lat_ts-latc)*np.pi/2.-np.atan2(Y,X)
#          Rmat  = [[np.cos(Rang),-np.sin(Rang)],
#                   [np.sin(Rang), np.cos(Rang)] ]
#       else:
# 
#       ########################################################################
# 
#       for n in range(len(lon_vec)):
#          lon   = lon_vec[n]
#          lat   = lat_vec[n]
#          v     = lonlat2vec(lon,lat)
# 
#          # line through v,vp:
#          line_coeffs = [vp-v,v] # [m,c], x=m*t+c: t=0 -> v, t=1 -> vp
# 
#          # intersection of line, projection plane
#          X,Y   = intersection_line_plane(plane_params,line_coeffs)
# 
#    else:
#       # convert from X,Y to lon,lat 
#       X,Y   = lon,lat
