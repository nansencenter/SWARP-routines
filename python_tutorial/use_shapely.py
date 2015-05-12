from matplotlib import pyplot as plt
import shapely.geometry as shgeom
import shapely.ops      as shops
# how to manipulate groups of polygons with shapely
# shops.unary_union looks v useful

# poly 1
x1 = [1.,1.,0.,0.,1.]
y1 = [0.,1.,1.,0.,0.]

if 1:
   # disjoint
   xoffset   = 2.
   yoffset   = 2.
elif 0:
   # intersecting
   xoffset   = .2
   yoffset   = .2
else:
   # touching
   xoffset   = 1. # disjoint
   yoffset   = 0. # disjoint


# poly 2
x2 = [x1[i]+xoffset for i in range(5)]
y2 = [y1[i]+yoffset for i in range(5)]

# coords are list of tuples:
c1 = []
c2 = []
for i in range(5):
   tup   = (x1[i],y1[i])
   c1.append(tup)
   tup   = (x2[i],y2[i])
   c2.append(tup)

poly1 = shgeom.Polygon(c1)
poly2 = shgeom.Polygon(c2)
Mp1   = shgeom.MultiPolygon([poly1,poly2])

# can use either of these
if 1:
   uu = shops.unary_union(Mp1)
else:
   uu = shops.unary_union([poly1,poly2])

if hasattr(uu,'geoms')==0:
   # uu is a single polygon
   print("uu is a single polygon")
   coords   = uu.boundary.coords
   x        = []
   y        = []
   for m in range(len(coords)):
      p0 = coords[m][0]
      p1 = coords[m][1]
      x.append(p0)
      y.append(p1)
      print(p0,p1)

   plt.plot(x,y)
   plt.show()
else:
   # uu is made of many polygons
   print("uu is made of "+str(len(uu.geoms))+" polygons")
   Ndone = 0
   for vv in uu.geoms:
      coords   = vv.boundary.coords
      x        = []
      y        = []
      print('\nPoly index: '+str(Ndone))
      for m in range(len(coords)):
         p0 = coords[m][0]
         p1 = coords[m][1]
         x.append(p0)
         y.append(p1)
         print(p0,p1)

      plt.plot(x,y)
      Ndone = Ndone+1

   plt.show()


if 0:
   HAVE_UU  = 0
   if poly1.boundary.intersects(poly2.boundary):
      print('intersect')
      uu = shops.unary_union([poly1,poly2])
      HAVE_UU  = 1

   if HAVE_UU==1:
      coords   = uu.boundary.coords
      x        = []
      y        = []
      for m in range(len(coords)):
         x.append(coords[m][0])
         y.append(coords[m][1])

      plt.plot(x,y)
      plt.show()
elif 0:
      uu = shops.unary_union([poly1,poly2])
      coords   = uu.boundary.coords
      x        = []
      y        = []
      for m in range(len(coords)):
         x.append(coords[m][0])
         y.append(coords[m][1])

      plt.plot(x,y)
      plt.show()
# a  = poly1.boundary.symmetric_difference(poly2.boundary)

# contains     - inside poly/multipoly or not?
# unary_union  - get boundary of group of poly's if they are neighbours/overlapping
