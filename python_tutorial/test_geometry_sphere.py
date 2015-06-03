import os,sys
import numpy as np

# test area on a sphere
if '../py_funs' not in sys.path:
   sys.path.append('../py_funs')

import geometry_sphere as GS

if 0:
   #test area on an ellipsoid calc (areaint.m)
   a              = 6378273   # semi-major axis
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
   # test conformal lat calc + its inverse
   a              = 6378273   # semi-major axis
   ecc            = .2        # eccentricity of ellipse
   ellipsoid_mat  = [a,ecc]   # ellipsoid formatted ala matlab

   phi   = 70
   chi   = GS.convertlat_geodetic2conformal(ecc,phi,inverse=False,radians=False)
   phi2  = GS.convertlat_geodetic2conformal(ecc,chi,inverse=True,radians=False)
   print(phi,chi,phi2)
   
