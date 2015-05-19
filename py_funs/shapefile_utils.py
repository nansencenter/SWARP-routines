import shapefile
import shapely.geometry as shgeom
import shapely.ops      as shops
import numpy as np

############################################################################
def DMI_form_dictionary():

   # dictionary to map from strings to integers
   form_vals      = {'X':-1,'-9':np.NaN}
   form_vals_max  = 10
   for n in range(form_vals_max+1):
      form_vals.update({str(n):n})
   
   return form_vals
############################################################################

############################################################################
def xy2list(x,y):
   
   lst   = []
   for n in range(len(x)):
      lst.append([x[n],y[n]])

   return lst
############################################################################

############################################################################
def get_poly(shp,get_holes=True):
   # convert points from a shapefile shape into a shapely polygon

   pts   = shp.points
   Npts  = len(pts)
   #
   parts    = shp.parts
   Nparts   = len(parts)
   parts.append(Npts)

   bdy   = pts[:parts[1]]

   #########################################################################
   if get_holes and (Nparts>1):
      # add rings (holes) (if present and wanted)
      rings = []
      for n in range(1,Nparts):
         rpts  = pts[parts[n]:parts[n+1]]
         rings.append(rpts)

      poly  = shgeom.Polygon(bdy,rings)
   else:
      poly  = shgeom.Polygon(bdy)
   #########################################################################

   return poly
############################################################################

############################################################################
def plot_poly(poly,pobj=None,plot_holes=True,**kwargs):

   if pobj is None:
      from matplotlib import pyplot as plt
      pobj  = plt

   # external polygon
   q     = shgeom.Polygon(poly.exterior)
   x,y   = q.boundary.coords.xy
   pobj.plot(x,y,**kwargs)

   if plot_holes:
      # interior polygons
      qi    = poly.interiors
      for q in qi:
         x,y   = shgeom.Polygon(q).boundary.coords.xy
         pobj.plot(x,y,linestyle='--',**kwargs)

   return
############################################################################

############################################################################
def extract_shapefile_info(sfile,get_holes=True):

   Nrec     = sfile.numRecords
   fields   = sfile.fields[1:]
   Nfields  = len(fields)

   out   = []

   for n in range(Nrec):
      shp   = sfile.shape(n)
      rec   = sfile.record(n)

      # get field values:
      lut   = {}
      for m in range(Nfields):
         fld   = fields[m][0]
         val   = rec[m]
         lut.update({fld:val})

      # get shapely polygon representation of points
      poly  = get_poly(shp,get_holes=get_holes)

      # add lookup table and shapely polygon to out
      out.append([poly,lut])

   return out
############################################################################

# ############################################################################
# def prepMultiPolygon4ShapeFile(MP):
#    # convert shapely MultiPolygon to form usable by shapefile
#    # also calculate field values
#    Np          = len(MP.geoms)
#    parts_list  = Np*[0]
#    rec_list    = Np*[0]
# 
#    for n in range(Np):
#       parts = []
# 
#       ############################################################################
#       poly  = MP.geoms[n]
#       q     = poly.exterior
#       x,y   = q.boundary.coords.xy
#       lst   = xy2list(x,y)
# 
#       # record list: [Area]
#       # TODO do stereographic projection
#       #      > calc area,perimeter and width
#       rec_list.append([q.area])
#       
#       # pyshp: if points are anti-clockwise, the points are treated as exterior ones
#       # TODO check if shapely does the same
#       parts = [lst]
# 
#       Nint  = len(poly.interiors)
#       for m in range(Nint):
#          x,y   = poly.interiors[m].boundary.coords.xy
#          lst   = xy2list(x,y)
# 
#          # pyshp: if points are clockwise, the points are treated as interior ones
#          # TODO check if shapely does the same
#          parts.append(lst)
#       ############################################################################
#       
#       parts_list.append(parts)
#    return parts_list,rec_list
# ############################################################################
# 
# ############################################################################
# def MultiPolygon2ShapeFile(MP,filename):
#    # save MultiPolygon to shapefile
#    # 1. convert shapely MultiPolygon to form usable by shapefile
#    # 2. save shapefile
# 
# 
#    w  = shapefile.Writer(shapefile.POLYGON)
# 
#    # define attributes
#    w.field('Area','N','40') # name,type ('C'=character, 'N'=number), size (?)
# 
#    # TODO calculate these and add fields for them
#    # w.field('Perimeter','N','40')
#    # w.field('Width','N','40')
# 
#    parts_list,rec_list  = prepMultiPolygon4ShapeFile(MP)
# 
#    Nrec  = len(parts_list)
#    for n in range(Nrec):
#       # add polygons
#       w.poly(parts=parts_list[n])
#       w.record(rec_list[n])
# 
#    # save file
#    w.save(filename)
# 
#    return
# ############################################################################

# w  = shapefile.Writer(shapefile.POLYGON)
# 
# # define attributes
# w.field('Area','N','40') # name,type ('C'=character, 'N'=number), size (?)
# w.field('Perimeter','N','40')
# w.field('Width','N','40')
# 
# # add a polygon
# coords_list = [[1,5],[5,5],[5,1],[3,3],[1,1],[1,5]]
# w.poly(parts=[coords_list])
# w.record(Area=40,Perimeter=80,Width=5)
# 
# # add another polygon
# w.poly(parts=[coords_list])
# w.record(Area=50,Perimeter=70,Width=6)
# 
# # save file
# w.save('test')
# 
# if 1:
#    sf    = shapefile.Reader('test.shp')
#    Nrec  = sf.numRecords
#    print('number of records: '+str(Nrec))
#    # for i in range(Nrec):
#    for i in [0]:
#       shp   = sf.shape(i)
