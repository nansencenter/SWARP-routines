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
def xy2tuple_list(x,y):
   # make list of tuples from x,y
   
   lst   = []
   for n in range(len(x)):
      lst.append((x[n],y[n]))

   return lst
############################################################################

############################################################################
def get_poly(shp,get_holes=True):
   # convert points from a shapefile shape into a shapely polygon

   import shapely.geometry as shgeom

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

   import shapely.geometry as shgeom

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
def plot_Mpoly(MP,pobj=None,colors=None,plot_holes=True,**kwargs):
   # plot MultiPolygon MP
   import numpy as np
   import shapely.geometry as shgeom

   if pobj is None:
      from matplotlib import pyplot as plt
      pobj  = plt
   
   col   = None
   Ng    = len(MP.geoms)

   for m,poly in enumerate(MP.geoms):
      if colors is not None:
         m_    = np.mod(m,len(colors))
         col   = colors[m_]

      L  = len(shgeom.Polygon(poly.exterior).boundary.coords)
      plot_poly(poly,pobj=pobj,color=col,plot_holes=plot_holes,**kwargs)
      if ((not poly.is_valid) and (L>50)) or (L>50):
         print("polgon "+str(m)+' ('+str(Ng)+'): '+str(len(poly.interiors))+' interiors')
         print('Valid: '+str(poly.is_valid))
         print('Length of boundary: '+str(L) +'\n')

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

############################################################################
def prepMultiPolygon4ShapeFile(MP):
   # convert shapely MultiPolygon to form usable by shapefile
   # also calculate field values
   import shapely.geometry as shgeom

   Np          = len(MP.geoms)
   parts_list  = Np*[0]
   rec_list    = Np*[0]

   for n in range(Np):
      parts = []

      ############################################################################
      # record list: [AREA_ID]
      rec_list[n] = [n]
      # TODO do stereographic projection
      #      > calc area,perimeter and width and append
      #      > eg rec_list[n].append(Area)
      ############################################################################

      ############################################################################
      poly  = MP.geoms[n]
      q     = shgeom.Polygon(poly.exterior)
      x,y   = q.boundary.coords.xy
      lst   = xy2list(x,y)

      # pyshp: if points are anti-clockwise, the points are treated as exterior ones
      # TODO check if shapely does the same
      parts = [lst]

      Nint  = len(poly.interiors)
      for m in range(Nint):
         q     = shgeom.Polygon(poly.interiors[m])
         x,y   = q.boundary.coords.xy
         lst   = xy2list(x,y)

         # pyshp: if points are clockwise, the points are treated as interior ones
         # TODO check if shapely does the same
         parts.append(lst)
      ############################################################################
      
      parts_list[n]  = parts
   return parts_list,rec_list
############################################################################

############################################################################
def MultiPolygon2ShapeFile(MP,filename):
   # save MultiPolygon to shapefile
   # 1. convert shapely MultiPolygon to form usable by shapefile
   # 2. save shapefile

   import shapefile
   w  = shapefile.Writer(shapefile.POLYGON)

   # define attributes
   w.field('AREA_ID','N','40')

   # TODO calculate these and add fields for them
   # w.field('Area','N','40') # name,type ('C'=character, 'N'=number), size (?)
   # w.field('Perimeter','N','40')
   # w.field('Width','N','40')

   parts_list,rec_list  = prepMultiPolygon4ShapeFile(MP)

   Nrec  = len(parts_list)
   for n in range(Nrec):
      # add polygons
      w.poly(parts=parts_list[n])
      rec   = rec_list[n]
      w.record(AREA_ID=rec[0])

   # save file
   w.save(filename)

   return
############################################################################

############################################################################
def SKcontours2MP(contours,x,y):

   import shapely.geometry as shgeom
   import numpy as np

   plist = []
   for cont in contours:
      ivec  = np.array(cont[:,0],dtype=int)
      jvec  = np.array(cont[:,1],dtype=int)
      #
      xvec  = x[ivec,jvec]
      yvec  = y[ivec,jvec]
      tups  = xy2tuple_list(xvec,yvec)
      #
      poly  = shgeom.Polygon(tups)
      if 1:#poly.is_valid:
         plist.append(poly)

   MP = shgeom.MultiPolygon(plist)
   return MP
############################################################################
