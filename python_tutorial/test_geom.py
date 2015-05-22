from matplotlib import pyplot as plt
import os,sys
import shapely.geometry as shgeom
import numpy as np

sys.path.append('../py_funs')
import shapefile_utils as SFU

#######################################################
def area_euclidean(x,y):
   # area of a polygon in Euclidean space
   A  = .5*(x[:-1]*y[1:]-x[1:]*y[:-1])
   print(A)
   return np.sum(A)
#######################################################

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
   # convert latitude to equivalent latitude on authalic sphere
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
def rng = greatcircledist(lat1, lon1, lat2, lon2, r):

   import numpy as np

   # Calculate great circle distance between points on a sphere using the
   # Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
   # length and has the same units as the radius of the sphere, R.  (If R is
   # 1, then RNG is effectively arc length in radians.)

   a     = pow( np.sin((lat2-lat1)/2),2) + np.cos(lat1) .* np.cos(lat2) .* pow(np.sin((lon2-lon1)/2.),2)
   rng   = r * 2. * np.atan2(np.sqrt(a),np.sqrt(1 - a))
   return rng
#######################################################

#######################################################
def az = greatcircleaz(lat1,lon1,lat2,lon2):

   import numpy as np
   # Inputs LAT1, LON1, LAT2, LON2 are in units of radians.

   az = atan2(np.cos(lat2) * np.sin(lon2-lon1),\
              np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2-lon1))

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

   ellipsoid_mat  = [1. 0.]

   # Start with spherical approximation even if using an ellipsoid
   rng = greatcircledist(lat1, lon1, lat2, lon2, ellipsoid_mat[0])

   # Azimuth on a sphere
   az = greatcircleaz(lat1, lon1, lat2, lon2)
   return rng,az
#######################################################

#######################################################
def singleint(lat,lon):
   # Compute the area of a single polygon on a unit sphere:

   # Ensure the path segment closes upon itself by
   # repeating beginning point at the end.

   if lat[-1]!=lat[0]:
      lat   = np.concatenate([lat,[lat[0]]],0)
      lon   = np.concatenate([lon,[lon[0]]],0)

   # Set origin for integration.  Any point will do, so (0,0) used
   lat0 = 0
   lon0 = 0

   # Get colatitude (a measure of surface distance as an angle)
   # and azimuth of each point in segment from the arbitrary origin

   colat,az = distance('gc',lat0,lon0,lat,lon,'radians')

   # Calculate step sizes
   daz   = az[1:]-az[:-1]

   # keep daz in  [-pi,pi] interval
   daz[daz<-np.pi]  = daz[daz<-np.pi]+np.pi
   daz[daz>np.pi]   = daz[daz>np.pi]-np.pi

   # Determine average surface distance for each step
   deltas   = (colat[1:]-colat[:-1])/2.
   colat    = colat[:-1]+deltas

   # Integral over azimuth is 1-cos(colatitudes)
   integrands  = (1-np.cos(colat))*daz;

   # Integrate and save the answer as a fraction of the unit sphere.
   # Note that the sum of the integrands will include a factor of 4\pi.
   area = abs(np.sum(integrands))/(4.*np.pi); # Could be area of inside or outside

   # Define the inside to be the side with less area.
   return np.min([area,1-area])

#######################################################
def area_ellipsoid(lon,lat,ellipsoid=None,radians=False):
   # translation of areaint.m

   if ellipsoid is None:
      import geopy.distance as geodist
      ellipsoid   = geodist.ELLIPSOIDS['WGS-84']

   a              = ellipsoid[0] # semi-major axis
   b              = ellipsoid[1]
   ecc            = np.sqrt(1-pow(b/a,2)) # eccentricity
   ellipsoid_mat  = [a,ecc]      # ellipsoid parameters used by matlab

   if not radians:
      # convert to radians
      lon   = lon*np.pi/180.
      lat   = lat*np.pi/180.

   # radius with same surface area
   radius = rsphere('authalic',ellipsoid_mat)

   # convert latitude to equivalent latitude on authalic sphere
   lat = convertlat_geodetic2authalic(ellipsoid, lat)
   
   A  = singleint(lon,lat)
   A  = A * 4.*np.pi*pow(radius,2)
   return A
#######################################################

# euclidean area of polygon
# http://www.wikihow.com/Calculate-the-Area-of-a-Polygon
x  = np.array([-3,-1,6,3,-4,-3]) 
y  = np.array([-2,4,1,10,9,-2])
p  = shgeom.Polygon(SFU.xy2tuple_list(x,y))
A  = p.area

fig   = plt.figure()
SFU.plot_poly(p,pobj=plt)

# triangle 1
x1 = np.array([-3,-1,-4,-3])
y1 = np.array([-2,4,9,-2])
p1 = shgeom.Polygon(SFU.xy2tuple_list(x1,y1))
A1 = p1.area
print(A1,area_euclidean(x1,y1))
SFU.plot_poly(p1,pobj=plt)

# triangle 2
x2 = np.array([-1,3,-4,-1]) 
y2 = np.array([4,10,9,4])
p2 = shgeom.Polygon(SFU.xy2tuple_list(x2,y2))
A2 = p2.area
print(A2,area_euclidean(x2,y2))
SFU.plot_poly(p2,pobj=plt)

# triangle 3
x3 = np.array([6,3,-1,6]) 
y3 = np.array([1,10,4,1])
p3 = shgeom.Polygon(SFU.xy2tuple_list(x3,y3))
A3 = p3.area
print(A3,area_euclidean(x3,y3))
SFU.plot_poly(p3,pobj=plt)

print(A,A1+A2+A3,area_euclidean(x,y))
plt.show()
