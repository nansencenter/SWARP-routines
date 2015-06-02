from matplotlib import pyplot as plt
import os,sys
import shapely.geometry as shgeom
import numpy as np

sys.path.append('../py_funs')
import shapefile_utils as SFU

import geom2

if 0:
   norm_vec = np.array([0,0,1])
   const    = 1
   plane    = geom2.plane(norm_vec,const)
   print('plane.normal_vector',plane.normal_vector)
   print('plane.constant',plane.constant)
   print('plane.basis_vectors',plane.basis_vectors)
else:
   norm_vec = np.array([1,0,0])
   const    = 2
   plane    = geom2.plane(norm_vec,const)
   print('plane.normal_vector',plane.normal_vector)
   print('plane.constant',plane.constant)
   print('plane.basis_vectors',plane.basis_vectors)
   #
   p1 = np.array([1,3,3])
   p2 = np.array([4,3,3])
   vi = plane.intersection_with_line([p1,p2],abs_coords=True)
   print(vi)
   #
   plane.reset_basis(vi)
   ai = plane.basis_coordinates(vi)
   print(ai)
   print(plane.basis_expansion(ai))
   print('plane.basis_vectors',plane.basis_vectors)
   # ai = plane.intersection_with_line([p1,p2],abs_coords=False)


sys.exit('L11')

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
