import os,sys
import numpy as np

# test area on a sphere
if '../py_funs' not in sys.path:
   sys.path.append('../py_funs')

import geometry_sphere as GS

if 0:
   #test area on an ellipsoid calc (compare to areaint.m)
   a              = 6378273.  # semi-major axis (m)
   ecc            = .1        # eccentricity of ellipse
   ellipsoid_mat  = [a,ecc]   # ellipsoid formatted ala matlab

   if 1:
      lons  = np.array([40,50,50,40])
      lats  = np.array([60,60,70,70])

      if 0:
         print('lons')
         print(lons)
         print(lons*np.pi/180.)
         #
         print('lats')
         print(lats)
         print(lats*np.pi/180.)

   area  = GS.area_polygon_ellipsoid(lons,lats,ellipsoid_mat=ellipsoid_mat)
   print('Area: '+str(area)+' m^2')

if 1:
   # test great circle dist
   a  = 6378273.  # semi-major axis (m)
   if 1:
      lat1,lon1   = 0.,0.
      # lat2,lon2   = 0.,1.
      lat2,lon2   = 0.,1.e-6
      dist2       = a*lon2*np.pi/180.

   if 1:
      lat1  = np.array([lat1])
      lon1  = np.array([lon1])
      lat2  = np.array([lat2])
      lon2  = np.array([lon2])

   dist  = GS.greatcircledist(lat1, lon1, lat2, lon2, R=a,radians=False)
   print(dist,dist2)



if 0:
   # test conformal lat calc + its inverse
   # (compare to convertlat.m)
   a              = 6378273   # semi-major axis (m)
   ecc            = .2        # eccentricity of ellipse
   ellipsoid_mat  = [a,ecc]   # ellipsoid formatted ala matlab

   phi   = 70
   chi   = GS.convertlat_geodetic2conformal(ecc,phi,inverse=False,radians=False)
   phi2  = GS.convertlat_geodetic2conformal(ecc,chi,inverse=True,radians=False)
   print(phi,chi,phi2)
