import numpy as np
import shapefile

class inverse_map:
   def __init__(self,map):
      self.map = map
      return

   def __call__(self,x,y,inverse=False):
      """
      (x,y)->(lon,lat)
      inverse=True: (lon,lat)->(x,y)
      """
      if not inverse:
         return self.map(x,y,inverse=True)
      else:
         return self.map(x,y,inverse=False)

class coord_transformer:
   def __init__(self,invmap,fwdmap):
      self.invmap = invmap
      self.fwdmap = fwdmap
      return

   def __call__(self,x,y):
      """
      (x,y)->(x',y')
      """
      lon,lat  = self.invmap(x,y,inverse=True)
      return self.fwdmap(lon,lat,inverse=False)

# ===============================================================
def default_pyproj(hemi='N',chart_source='DMI'):

   import pyproj
   if chart_source=='NIC' and hemi=='N':
      # WGS84 ellipsoid
      a        = 6378137.0
      f        = 1./298.257223563 # flattening
      b        = a*(1-f)
      lon0     = 180.
      lat0     = 90.
      lat_ts   = 60.
   else:
      ecc      = 0.081816153
      a        = 6378.273e3 #equatorial radius in m
      b        = a*np.sqrt(1-pow(ecc,2))
      lon0     = -45.
      if hemi=='N':
         lat0     = 90.
         lat_ts   = 60.
      else:
         lat0     = -90.
         lat_ts   = -70.

   return pyproj.Proj(proj='stere',a=a,b=b,\
                  lon_0=lon0,lat_0=lat0,lat_ts=lat_ts)
# ===============================================================


# ===============================================================
def map_coords(coords,inverse=False,mapping=None,return_bboxes=False,**kwargs):
   import numpy as np

   if mapping is None:
      mapping  = default_pyproj(**kwargs)

   if not inverse:
      lon,lat  = zip(*coords)
      x,y      = mapping(lon,lat,inverse=False)
      cout     = zip(x,y)
   else:
      x,y      = zip(*coords)
      lon,lat  = mapping(x,y,inverse=True)
      cout     = zip(lon,lat)

   if not return_bboxes:
      return cout
   else:
      bbox_ll  = [np.min(lon),np.min(lat),np.max(lon),np.max(lat)]
      bbox_xy  = [np.min(x),np.min(y),np.max(x),np.max(y)]
      return cout,bbox_ll,bbox_xy
# ===============================================================


# ===============================================================
def new_bbox(b1,b2):
   import numpy as np

   if b1 is None and b2 is None:
      return
   elif b1 is None:
      return b2
   elif b2 is None:
      return b1
   else:
      return [np.min([b1[0],b2[0]]),\
              np.min([b1[1],b2[1]]),\
              np.max([b1[2],b2[2]]),\
              np.max([b1[3],b2[3]])]
# ===============================================================


# ===============================================================
def is_in_bbox(bbox,x,y):
   import numpy as np

   ok = np.logical_and(x>=bbox[0],x<=bbox[2])
   ok = np.logical_and(y>=bbox[1],ok)
   return np.logical_and(ok,y<=bbox[3])
# ===============================================================


# ==========================================================================
def form_dictionary():

   # dictionary to map from strings to integers
   form_vals      = {'X':-1,'-9':np.NaN,'99':np.NaN}
   form_vals_max  = 10
   for n in range(form_vals_max+1):
      form_vals.update({str(n):n})
   
   return form_vals
############################################################################


############################################################################
def form_cols(cmapname=None,keys=None,return_cmap=False):
   """
   make a dictionary mapping form to rgba color
   - for plotting
   """
   import matplotlib.pyplot as plt
   import matplotlib.cm as mplcm
   import matplotlib.colors as colors

   if cmapname is None:
      # see http://matplotlib.org/examples/color/colormaps_reference.html
      # for a full list of colormaps
      cmapname = 'Set1'

   cm       = plt.get_cmap(cmapname)
   addnan   = False

   if keys is None:
      # default is for plotting of form for DMI/AARI ice charts
      keys     = range(-1,15)
      addnan   = True
   else:
      cols        = {}

   num_colors  = len(keys)
   cNorm       = colors.Normalize(vmin=0, vmax=num_colors)
   scalarMap   = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
   col_list    = [scalarMap.to_rgba(i+1) for i in range(len(keys))] # add an extra color for nan

   for i,key in enumerate(keys):
      cols.update({key:col_list[i]})

   if addnan:
      cols.update({np.nan:col_list[0]})

   if return_cmap:
      return cols,cm
   else:
      return cols
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

def poly_bounds(poly,mapping=None):
   bbox_ll  = poly.bounds
   if mapping is None:
      return bbox_ll

   lon,lat  = poly.exterior.coords.xy
   x,y      = mapping(lon,lat)
   xmin     = np.min(x)
   ymin     = np.min(y)
   xmax     = np.max(x)
   ymax     = np.max(y)
   bbox_xy  = xmin,ymin,xmax,ymax
   return bbox_ll,bbox_xy


############################################################################
def get_poly(shp,mapping=None,get_holes=True,latlon=True,return_bboxes=False):
   # convert points from a shapefile shape into a shapely polygon

   import shapely.geometry as shgeom

   pts   = shp.points
   Npts  = len(pts)
   #
   parts    = shp.parts
   Nparts   = len(parts)
   parts.append(Npts)


   Nrings   = 0
   bdy      = pts[:parts[1]] # list of tuples - closed
   if return_bboxes:
      if not latlon:
         # points in x,y
         bbox_ll,bbox_xy  = map_coords(bdy,inverse=True,mapping=mapping,return_bboxes=True)[1:]
      else:
         # points in lon,lat
         bbox_ll,bbox_xy  = map_coords(bdy,inverse=False,mapping=mapping,return_bboxes=True)[1:]


   #########################################################################
   if get_holes and (Nparts>1):
      # add rings (holes) (if present and wanted)
      rings = []
      for n in range(1,Nparts):
         rpts  = pts[parts[n]:parts[n+1]]
         rings.append(rpts)

      # make a shapely polygon
      # - buffer(0) usually makes it valid if it's not already

      poly  = shgeom.Polygon(bdy,rings)
   else:
      # make a shapely polygon
      # - buffer(0) usually makes it valid if it's not already
      poly  = shgeom.Polygon(bdy)
   #########################################################################

   if not poly.is_valid:
      poly  = poly.buffer(0)

   # can sometimes make a multipolygon from making poly valid
   polys = []
   if hasattr(poly,'geoms'):
      polys = []
      for poly in poly.geoms:
         Nrings   = len(poly.interiors)
         polys.append([poly,Nrings])
   else:
      Nrings   = len(poly.interiors)
      polys.append([poly,Nrings])

   
   # try to ensure poly is valid
   if not return_bboxes:
      return polys
   else:
      return polys,bbox_ll,bbox_xy

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
def extract_shapefile_info(sfile,**kwargs):
   """
   shapefile_utils.extract_shapefile_info(sfile,**kwargs)

   inputs:
   sfile is a shapefile.Reader object
   kwargs for shapefile_utils.get_poly()
   eg extract_shapefile_info(sfile,latlon=True,get_holes=True,mapping=None)

   returns:
   out,bbox_ll,bbox_xy
   out      = list - each member is a list [poly,lut,Nrings]
   bbox_ll  = [lon_min,lat_min,lon_Max,lat_max]
   bbox_xy  = [x_min,y_min,x_Max,y_max]
   """

   Nrec     = sfile.numRecords
   fields   = sfile.fields[1:] # remove 'DeletionFlag'
   Nfields  = len(fields)

   out      = []
   bbox_ll  = None
   bbox_xy  = None
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
      polys,bbox_ll0,bbox_xy0 = get_poly(shp,return_bboxes=True,**kwargs)

      # update the bounding boxes
      bbox_ll  = new_bbox(bbox_ll,bbox_ll0)
      bbox_xy  = new_bbox(bbox_xy,bbox_xy0)

      # add lookup table and shapely polygon to out
      for poly,Nrings in polys:
         if (poly.exterior is not None and poly.is_valid):
            out.append([poly,lut,Nrings])

   return out,bbox_ll,bbox_xy
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


############################################################################
class shapes:
   def __init__(self,sf_info):
      self.polygons  = []
      self.records   = []
      self.Nholes    = []
      for poly,lut,Nrings in sf_info:
         self.polygons.append(poly)
         self.records .append(lut)
         self.Nholes  .append(Nrings)
      return
############################################################################

class shapefile_info:
   def __init__(self,fname_full,hemi='N',mapping=None,latlon=True,chart_source="DMI"):

      # for info
      self.shapefile    = fname_full
      self.chart_source = chart_source
      self.form_dict    = None # not all ice charts give info about "form" (floe size)
      self.latlon       = latlon

      # set mapping
      if mapping is None:
         self.mapping   = default_pyproj(hemi=hemi,chart_source=chart_source)
      else:
         self.mapping   = mapping

      import shapefile
      print('Reading '+fname_full+'...')
      sf = shapefile.Reader(fname_full)

      # ============================================================
      # get shapes and records from shapefile
      # self.shapes.records[n]
      # self.shapes.polygons[n]
      # give the info for the n-th shape
      # *poly is a shapely polygon
      # *record is the metadata
      sf_info,self.bbox_ll,self.bbox_xy   =\
            extract_shapefile_info(sf,get_holes=True,mapping=self.mapping,latlon=self.latlon)
      self.shapes = shapes(sf_info)
      self.Npolys = len(sf_info)

      self.fields = sf.fields[1:]
      self.sort_fields()

      print('\nExample record:')
      print(self.shapes.records[0])
      print('\n')
      # ============================================================


      # ============================================================
      # set limits
      # - mostly for plotting
      self.set_limits()
      return
   # ===============================================================


   # ===============================================================
   def set_limits(self):

      self.lon_min,self.lat_min,self.lon_max,self.lat_max   = \
            self.bbox_ll
      self.x_min,self.y_min,self.x_max,self.y_max   = \
            self.bbox_xy

      return
   # ============================================================


   # ============================================================
   def sort_fields(self):

      self.fieldnames   = []
      self.fieldnames_C = [] # characters
      self.fieldnames_N = [] # integers
      self.fieldnames_F = [] # floats

      dicC  = {}
      for i,fld in enumerate(self.fields):
         self.fieldnames.append(fld[0])
         if fld[1]=='C':
            dicC.update({fld[0]:[]})
            self.fieldnames_C.append(fld[0])
         elif fld[1]=='N':
            self.fieldnames_N.append(fld[0])
         elif fld[1]=='F':
            self.fieldnames_F.append(fld[0])

      # determine the range of values each character type field can take
      for rec in self.shapes.records:
         for fld in self.fieldnames_C:
            dicC[fld].append(rec[fld])

      self.fieldranges_C   = {}
      for fld in self.fieldnames_C:
         self.fieldranges_C.update({fld:list(set(dicC[fld]))})

      return
   # ============================================================


   # ============================================================
   def test_plot(self,vname=None,latlon=True,figname=None,\
         cmapname=None,cres='i',poly_num=None,plot_holes=True):

      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      fig      = plt.figure(figsize=(12,12))
      ax       = fig.add_subplot(111)

      if not self.latlon:
         # can be a mess if polygons not originally in lon,lat
         latlon   = False

      # ====================================================
      # limits of figure (not automatic with patches)
      if latlon:
         ax.set_xlim([self.lon_min,self.lon_max])
         ax.set_ylim([self.lat_min,self.lat_max])
         ax.set_xlabel('Longitude, $^\circ$E',fontsize=18)
         ax.set_ylabel('Latitude, $^\circ$N',fontsize=18)
         bmap  = None
         if self.latlon:
            # stay in lon,lat
            pmap  = None
         else:
            # x,y -> lon,lat
            pmap  = inverse_map(self.mapping)
      else:
         from mpl_toolkits.basemap import Basemap
         dic   = Fplt.pyproj_srs_to_dict(self.mapping.srs)
         dic.update({'width' :self.x_max-self.x_min})
         dic.update({'height':self.y_max-self.y_min})
         dic.update({'resolution':cres})

         # coastal resolution:
         # ``c`` (crude), ``l`` (low), ``i`` (intermediate),
         # ``h`` (high), ``f`` (fine)
         bmap  = Basemap(**dic)

         # draw the coast
         Fplt.finish_map(bmap,ax=ax)

         if self.latlon:
            # want to go from lon,lat to basemap coords
            pmap  = bmap
         else:
            # want to go from native to basemap coords
            pmap  = coord_transformer(self.mapping,bmap)
      # ====================================================


      # ====================================================
      # get colormap
      import matplotlib.patches as mpatches
      if vname is None:
         # AARI/DMI ice charts: plot form
         clist = form_cols(cmapname=cmapname)
         if poly_num is None:
            # plot all polygons
            for n,poly in enumerate(self.shapes.polygons):
               rec   = self.shapes.records[n]
               if rec['POLY_TYPE']=='W':
                  Fplt.plot_patches(ax,[poly],'c',mapping=pmap,plot_holes=plot_holes)
               else:
                  key   = self.form_dict[rec['FA']]
                  col   = clist[key]
                  Fplt.plot_patches(ax,[poly],col,mapping=pmap,plot_holes=plot_holes)
         else:
            # plot one polygon
            poly  = self.shapes.polygons[poly_num]
            rec   = self.shapes.records[poly_num]
            if rec['POLY_TYPE']=='W':
               Fplt.plot_patches(ax,[poly],'c',mapping=pmap,plot_holes=plot_holes)
            else:
               key   = self.form_dict[rec['FA']]
               col   = clist[key]
               Fplt.plot_patches(ax,[poly],col,mapping=pmap,plot_holes=plot_holes)

         # for legend
         Patches  = [mpatches.Patch(facecolor='c', edgecolor='k', label='W')]

      elif vname in self.fieldranges_C:
         frng     = self.fieldranges_C[vname]
         numcols  = len(frng)
         clist    = form_cols(keys=frng)

         if poly_num is None:
            # plot all
            for n,poly in enumerate(self.shapes.polygons):
               # print(n,poly)
               # print(poly.length,poly.area)
               rec   = self.shapes.records[n]
               # print(rec)
               val   = rec[vname]
               col   = clist[val]
               Fplt.plot_patches(ax,[poly],col,mapping=pmap,plot_holes=plot_holes)
         else:
            # plot one polygon
            poly  = self.shapes.polygons[poly_num]
            rec   = self.shapes.records[poly_num]
            val   = rec[vname]
            col   = clist[val]
            Fplt.plot_patches(ax,[poly],col,mapping=pmap,plot_holes=plot_holes)

         # for legend
         Patches  = []
      # ====================================================
      

      # ====================================================
      # legend
      for key in sorted(clist):
         col   = clist[key]
         if type(key)==type('hi'):
            lab   = key
         elif np.isnan(key):
            lab   = '-'
         elif key==-1:
            lab   = 'X'
         else:
            lab   = str(key)
         Patches.append(mpatches.Patch(facecolor=col,edgecolor='k',label=lab))

      ax.legend(handles=Patches,loc='lower right')
      # ====================================================


      for li in ax.get_xticklabels():
         li.set_fontsize(18)
      for li in ax.get_yticklabels():
         li.set_fontsize(18)


      # ====================================================
      # show or save fig
      if figname is None:
         plt.show(fig)
      else:
         print('\nSaving '+figname+'...\n')
         fig.savefig(figname,bbox_inches='tight')
         ax.cla()
         plt.close(fig)
      # ====================================================

      return
   # ============================================================


   # ============================================================
   def test_plot_grid(self,vname=None,figname=None,cres='i',\
         cmapname='bwr',res=10.e3, use_bmap=False,dic=None,\
         use_legend=True,**kwargs):

      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      fig      = plt.figure(figsize=(12,12))
      ax       = fig.add_subplot(111)

      # ====================================================
      # limits of figure (not automatic with patches)
      Dx = self.x_max-self.x_min
      Dy = self.y_max-self.y_min

      if use_bmap:
         from mpl_toolkits.basemap import Basemap
         Dic   = Fplt.pyproj_srs_to_dict(self.mapping.srs)
         Dic.update({'width' :Dx})
         Dic.update({'height':Dy})

         # coastal resolution:
         # ``c`` (crude), ``l`` (low), ``i`` (intermediate),
         # ``h`` (high), ``f`` (fine)
         Dic.update({'resolution':cres})
         bmap  = Basemap(**Dic)

         # draw the coast
         Fplt.finish_map(bmap,ax=ax)
      # ====================================================

      nx    = int(Dx/res)
      ny    = int(Dy/res)
      dx    = Dx/nx
      dy    = Dy/ny
      xvec  = np.linspace(self.x_min,self.x_max,nx) # corners (for pcolor)
      yvec  = np.linspace(self.y_min,self.y_max,ny)

      # X,Y   = np.meshgrid(xvec[1:]-.5*dx,yvec[1:]-.5*dy) # centres
      Y,X   = np.meshgrid(yvec[1:]-.5*dy,xvec[1:]-.5*dx) # centres

      # ====================================================
      # get colormap
      if use_legend:
         # for legend
         import matplotlib.patches as mpatches
         Patches  = []

      if vname is None:
         # AARI/DMI ice charts: plot form
         clist,cmap  = form_cols(cmapname=cmapname,return_cmap=True)

         if dic is None:
            keys  = self.form_dict.keys()
            dic   = {}
            i     = 1
            for key in keys:
               i       += 1
               dic.update({key:float(i)})

         Z,dic    = self.label_grid(X,Y,'FA',dic=dic,latlon=False)
         maskw    = self.is_cat(X,Y,'POLY_TYPE','W',latlon=True) #water
         Z[maskw] = 1.
         dic.update({'W':1.})

         if use_legend:
            Patches.append(mpatches.Patch(facecolor='c',edgecolor='k', label='W'))

      elif vname in self.fieldranges_C:
         if dic is None:
            frng  = self.fieldranges_C[vname]
         else:
            frng  = dic.keys()

         clist,cmap  = form_cols(keys=frng,cmapname=cmapname,\
               return_cmap=True)
         print(clist)
         Z,dic    = self.label_grid(X,Y,vname,dic=dic,latlon=False)
      # ====================================================
      
      if use_bmap:
         Yc,Xc       = np.meshgrid(yvec,xvec) # corners
         Lonc,Latc   = self.mapping(Xc,Yc,inverse=True)
         im          = bmap.pcolor(Lonc,Latc,Z,cmap=cmap,latlon=True,ax=ax,**kwargs)
         legloc      = 'lower right'
      else:
         im       = ax.imshow(Z.transpose(),cmap=cmap,origin='lower',**kwargs)
         legloc   = 'upper left'


      if use_legend:
         # ====================================================
         # legend
         for key in sorted(clist):
            col   = clist[key]
            if type(key)==type('hi'):
               lab   = key
            elif np.isnan(key):
               lab   = '-'
            elif key==-1:
               lab   = 'X'
            else:
               lab   = str(key)
            Patches.append(mpatches.Patch(facecolor=col,label=lab,edgecolor='k'))

         ax.legend(handles=Patches,loc=legloc)
         # ====================================================

      else:
         cbar  = fig.colorbar(im)


      # for li in ax.get_xticklabels():
      #    li.set_fontsize(18)
      # for li in ax.get_yticklabels():
      #    li.set_fontsize(18)



      # ====================================================
      # show or save fig
      if figname is None:
         plt.show(fig)
      else:
         print('\nSaving '+figname+'...\n')
         fig.savefig(figname,bbox_inches='tight')
         ax.cla()
         plt.close(fig)
      # ====================================================

      return Z,dic
   # ============================================================


   # =============================================================================
   def is_cat(self,X,Y,vname,cat,latlon=True):
      from matplotlib import path
      
      Nx    = X.size
      shp   = X.shape

      if latlon:
         # project from lon/lat into x,y
         X,Y   = self.mapping(X,Y)

      X        = X.reshape(Nx)
      Y        = Y.reshape(Nx)
      Mask     = np.zeros((Nx),dtype='bool')
      coords2  = np.array([X,Y]).transpose()

      # loop over shapes
      for n,poly in enumerate(self.shapes.polygons):
         
         # check the record
         rec   = self.shapes.records[n]
         if rec[vname]!=cat:
            continue

         lon,lat  = poly.exterior.coords.xy
         x,y      = self.mapping(lon,lat)
         coords   = np.array([x,y]).transpose()
         bbPath   = path.Path(coords)

         # test if inside external polygon
         mask  = np.array(bbPath.contains_points(coords2),dtype='bool')

         for qi in poly.interiors:
            # check not inside holes
            lon,lat  = qi.coords.xy
            x,y      = self.mapping(lon,lat)
            coords   = np.array([x,y]).transpose()
            bbPath   = path.Path(coords)
            inhole   = np.array(bbPath.contains_points(coords2),dtype='bool')
            mask     = np.logical_and(mask,np.logical_not(inhole))
            
         Mask[mask]  = True

      return Mask.reshape(shp)
   # =============================================================================

   def label_grid(self,Xin,Yin,vname,dic=None,latlon=True):
      from matplotlib import path
      import numpy as np

      
      if dic is None:
         dic   = {'None':0.}
         for i,key in enumerate(self.fieldranges_C[vname]):
            dic.update({key:float(i+1)})
         print(dic)
      elif 'None' not in dic:
         zmin  = 1.e30
         for key in dic:
            zmin  = np.min(zmin,dic[key])
         dic.update({'None':zmin-1.})

      print(vname)
      print(dic)


      Nx    = Xin.size
      shp   = Xin.shape
      Z0    = np.zeros(shp) + dic['None']
      if latlon:
         # project from lon/lat into x,y
         X0,Y0 = self.mapping(Xin,Yin)
      else:
         # just point to inputs
         # - no problem, since not changing them
         X0 = Xin
         Y0 = Yin


      # loop over shapes
      points_found   = 0
      for n,poly in enumerate(self.shapes.polygons):

         # # to test 1 polygon
         # if n!=43:
         #    continue

         print(n,len(self.shapes.records))

         # check the record
         rec   = self.shapes.records[n]
         key   = rec[vname]
         if key not in dic:
            # not interesting just assign 'None'
            continue
         else:
            zval  = dic[key]
            print(key,zval)

         # test if inside external polygon
         if self.latlon:
            lon,lat  = poly.exterior.coords.xy
            x,y      = self.mapping(lon,lat)
         else:
            x,y   = poly.exterior.coords.xy
         coords   = np.array([x,y]).transpose()
         bbPath   = path.Path(coords)
         bbox     = [np.min(x),np.min(y),np.max(x),np.max(y)]


         # reduce the search with polygon bbox
         in_bbox  = is_in_bbox(bbox,X0,Y0)
         X        = X0[in_bbox]
         Y        = Y0[in_bbox]
         coords2  = np.array([X,Y]).transpose()


         # reduce search further for holes
         inside   = np.array(bbPath.contains_points(coords2),dtype='bool')
         X1       = X[inside]
         Y1       = Y[inside]
         inholes  = np.zeros(X1.shape,dtype='bool')
         coords3  = np.array([X1,Y1]).transpose()

         for qi in poly.interiors:
            # check not inside holes
            if self.latlon:
               lon,lat  = qi.coords.xy
               x,y      = self.mapping(lon,lat)
            else:
               x,y   = qi.coords.xy

            coords            = np.array([x,y]).transpose()
            bbPath            = path.Path(coords)
            inhole            = np.array(bbPath.contains_points(coords3),dtype='bool')
            inholes[inhole]   = True

         inside[inside] = np.logical_not(inholes)
         Z              = Z0[in_bbox]
         Z[inside]      = zval
         Z0[in_bbox]    = Z

         points_found += len(Z[inside])

      print('found '+str(points_found)+', out of '+str(Nx))

      Z0[in_bbox] = Z
      return Z0,dic

############################################################################
class MIZ_from_shapefile:
   def __init__(self,fname_full,MIZ_criteria="FA_only",chart_source="DMI"):
      """
      MIZshp   = MIZ_from_shapefile(fname_full,MIZ_criteria="FA_only")
      """

      # for info
      self.shapefile    = fname_full
      self.chart_source = chart_source
      self.MIZ_criteria = MIZ_criteria

      sf = shapefile.Reader(fname_full)
      print('Processing '+fname_full+'...')

      # ============================================================
      # get shapes and records from shapefile
      # as list of lists:
      # sf_info[n]   = [poly,record]
      # where poly is a shapely polygon,
      # and record is all the metadata
      sf_info,self.bbox_ll,self.bbox_xy   =\
            extract_shapefile_info(sf,get_holes=True)
      self.shapes = shapes(sf_info)
      self.Npolys = len(sf_info)

      print('\nExample record:')
      print(self.shapes.records[0])
      print('\n')
      # ============================================================


      # ================================================================================
      # MIZ definition stuff
      self.MIZ_criteria = MIZ_criteria
      if MIZ_criteria=="FA_only":
         form_cats_MIZ     = ['FA']
         match_all_crits   = True
         #  in MIZ if (FA<=MIZthresh)
      elif MIZ_criteria=="FB_only":
         form_cats_MIZ     = ['FB']
         match_all_crits   = True
         #  in MIZ if (FB<=MIZthresh)
      elif MIZ_criteria=="FC_only":
         form_cats_MIZ     = ['FB']
         match_all_crits   = True
         #  in MIZ if (FC<=MIZthresh)
      else:
         form_cats_MIZ  = ['FA','FB','FC']
         if MIZ_criteria=="all_Fcats_small":
            match_all_crits   = True
            #  in MIZ if (FA<=MIZthresh) and (FB<=MIZthresh) and (FC<=MIZthresh)
         elif MIZ_criteria=="any_Fcats_small":
            match_all_crits   = False
            #  in MIZ if (FA<=MIZthresh) or (FB<=MIZthresh) or (FC<=MIZthresh)
         else:
            raise ValueError("unknown option for 'MIZ_criteria'\n"+\
                              "- valid options: 'FA_only', 'FB_only', 'FC_only',\n"+\
                              " 'all_Fcats_small' or 'any_Fcats_small'")

      form_vals   = form_dictionary()
      MIZthresh   = 4
      # ================================================================================


      # ============================================================
      # separate water, MIZ and pack:
      # find MIZ from values of FA,FB,FC
      self.MIZ_forms    = [] # MIZ polygons
      self.WTR_forms    = [] # water polygons
      self.PACK_forms   = [] # pack ice polygons

      print('\nTesting if polygons are in the MIZ...')
      for n,lut in enumerate(self.shapes.records):

         print(lut)
         if lut['POLY_TYPE']=='W':
            # water
            self.WTR_forms.append(n)
            # print(n,lut)
         else:
            # ========================================================
            # Try to define MIZ:
            fdct  = {}
            bools = np.zeros((len(form_cats_MIZ)))
            for i,key in enumerate(form_cats_MIZ):
               sval  = lut[key]        # string 
               val   = form_vals[sval] # integer from form_vals dictionary 
               fdct.update({key:val})
               
               bools[i] = val<=MIZthresh

            if match_all_crits:
               # define MIZ if all of the forms are <= threshhold
               isMIZ = np.all(bools)
            else:
               # define MIZ if any of the forms are <= threshhold
               isMIZ = np.any(bools)
            # ========================================================

            if isMIZ:
               self.MIZ_forms.append([n,fdct])
               # print(n,lut)
            else:
               # not MIZ but ice
               self.PACK_forms.append(n)
      # ============================================================


      # merge water,pack and MIZ polygons
      self.merge_polygons()

      # classify MIZ boundary
      self.classify_MIZ_bdy()

      # determine which of the original shapes
      # correspond to which of the merged polygons
      # - can then get areas of each merged one
      self.get_areas()

      return

   def write_text_file(self,tfil1,MIZtype=None):
      """
      write_text_file(self,textfil,MIZtype=None)
      MIZtype='inner','outer' or None (both in/out)
      """

      if self.N_MIZ==0:
         print("No MIZ polygons - not writing text file")
         return

      # ============================================================
      # write text files
      # tfil2 = Fname+'_MIZpolys_info.txt'
      print('\nSaving polygons + boundary flags to '+tfil1)
      # print('Saving metadata to '+tfil2)
      print('\n')
      tf1   = open(tfil1,'w')
      # tf2   = open(tfil2,'w')
      blk   = 3*' '

      Count = 0
      for npoly,dct in enumerate(self.MIZclass):
         npts  = len(dct['lon'])

         if MIZtype is not None:
            if dct['type']!=MIZtype:
               continue

         # ==================================================
         # headers
         if Count==0:
            ss = 'Polygon'+blk+\
                 'lon'+blk+\
                 'lat'+blk+\
                 'flag\n'
            tf1.write(ss)

            # keys  = dct['record'].keys()
            # ss    = 'Polygon'+blk
            # for key in keys[:-1]:
            #    ss   +=key+blk
            # ss   += keys[-1]+'\n'
            # tf2.write(ss)
         # ==================================================

         # # ==================================================
         # # metadata for poly
         # ss = str(npoly)+blk
         # for key in keys[:-1]:
         #    ss   += str(dct['record'][key])+blk
         # ss   += str(dct['record'][keys[-1]])+'\n'
         # tf2.write(ss)
         # # ==================================================


         # ==================================================
         # all points + classification
         for npt in range(npts):
            ss0   = [Count,\
                     dct['lon'][npt],
                     dct['lat'][npt],
                     dct['flags'][npt]]

            ss = ''
            for val in ss0[:-1]:
               ss   += str(val)+blk
            ss   += str(ss0[-1])+"\n"
            tf1.write(ss)
         # ==================================================

         Count+= 1

      tf1.close()
      # tf2.close()
      # ============================================================
      return

   def merge_polygons(self):
      import shapely.ops      as shops
      # ============================================================
      # merge pack polygons (also dissolving holes)
      MPlist_PACK = []
      Npack       = len(self.PACK_forms)

      print('Merging '+str(Npack)+' pack polygons...')
      bad   = []
      for n in self.PACK_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_PACK.append(P)
         else:
            bad.append(n)

      for n in bad:
         self.PACK_forms.remove(n)

      self.MP_PACK  = shops.unary_union(MPlist_PACK)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_PACK,'geoms'):
         # single polygon
         self.MPlist_PACK  = [self.MP_PACK]
      else:
         # many polygons
         self.MPlist_PACK  = []
         for poly in self.MP_PACK.geoms:
            self.MPlist_PACK.append(poly)

      self.N_PACK = len(self.MPlist_PACK)
      print(str(self.N_PACK)+ ' pack polygons after merging\n')
      # ============================================================


      # ============================================================
      # merge wtr polygons (also dissolving holes)
      MPlist_WTR  = []
      Nwtr        = len(self.WTR_forms)

      print('Merging '+str(Nwtr)+' water polygons...')
      bad   = []
      for n in self.WTR_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_WTR.append(P)
         else:
            bad.append(n)

      for n in bad:
         self.WTR_forms.remove(n)

      self.MP_WTR = shops.unary_union(MPlist_WTR)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_WTR,'geoms'):
         # single polygon
         self.MPlist_WTR   = [self.MP_WTR]
      else:
         # many polygons
         self.MPlist_WTR  = []
         for poly in self.MP_WTR.geoms:
            self.MPlist_WTR.append(poly)

      self.N_WTR  = len(self.MPlist_WTR)
      print(str(self.N_WTR)+ ' water polygons after merging\n')
      # ============================================================


      # ============================================================
      # merge MIZ polygons (also dissolving holes)
      MPlist_MIZ  = []
      N_MIZ       = len(self.MIZ_forms)

      print('Merging '+str(N_MIZ)+' MIZ polygons...')
      bad   = []
      for n,fdct in self.MIZ_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_MIZ.append(P)
         else:
            bad.append([n,fdct])

      for nd in bad:
         self.MIZ_forms.remove(nd)

      self.MP_MIZ = shops.unary_union(MPlist_MIZ)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_MIZ,'geoms'):
         # 1 polygon
         self.MPlist_MIZ  = [self.MP_MIZ]
      else:
         self.MPlist_MIZ  = []
         for poly in self.MP_MIZ.geoms:
            self.MPlist_MIZ.append(poly)

      self.N_MIZ  = len(self.MPlist_MIZ)
      print(str(self.N_MIZ)+ ' MIZ polygons after merging\n')
      # ===============================================================
      return
   
   def classify_MIZ_bdy(self):

      # ==================================================================
      # TODO match to original shapes to extract correct areas
      # TODO are holes important?
      N_MIZ_merged   = len(self.MPlist_MIZ)
      print('Sorting '+str(N_MIZ_merged)+' MIZ polygons...\n')

      self.MIZclass  = []
      if N_MIZ_merged==0:
         return

      count          = 0
      for poly in self.MPlist_MIZ:
         count += 1
         print("Progress: "+str(count)+"/"+str(N_MIZ_merged))
         llc            = list(poly.exterior.coords)
         lonc,latc      = poly.exterior.coords.xy
         Nc             = len(llc)
         boundary_flags = Nc*[2] # 2 = unknown, 0 = ice edge, 1 = MIZ-pack edge

         # =================================================
         # ice edge
         buf   = 1.e-5
         if poly.exterior.intersects(self.MP_WTR.buffer(buf)):
            IE = poly.exterior.intersection(self.MP_WTR.buffer(buf))
               # coords in ice edge

            if not hasattr(IE,'geoms'):
               LLc   = list(IE.coords)
            else:
               LLc   = []
               for IE_ in IE.geoms:
                  LLc.extend(list(IE_.coords))

            # reduce list of points to those in the original poly
            # - intersection gets filled somehow
            llc   = list(poly.exterior.coords)
            for cc in LLc:
               if cc in llc:
                  idx   = llc.index(cc)
                  boundary_flags[idx]  = 0
         # =================================================


         # =================================================
         # MIZ-pack edge
         if poly.exterior.intersects(self.MP_PACK.buffer(buf)):
            IE = poly.exterior.intersection(self.MP_PACK.buffer(buf))
               # coords in MIZ-pack edge

            if not hasattr(IE,'geoms'):
               LLc   = list(IE.coords)
            else:
               LLc   = []
               for IE_ in IE.geoms:
                  LLc.extend(list(IE_.coords))

            # reduce list of points to those in the original poly
            # - intersection gets filled somehow
            llc   = list(poly.exterior.coords)
            for cc in LLc:
               if cc in llc:
                  idx   = llc.index(cc)
                  boundary_flags[idx]  = 1
         # =================================================


         # ============================================================
         # finalise outputs
         # - which list to append poly to?
         # (inner MIZ or outer/proper MIZ)
         for flag in range(3):
            print('#'+str(flag)+': '+str(boundary_flags.count(flag)))
         print('Total number of points: '+str(Nc)+"\n")

         Dct   = {'lon'    :np.array(lonc),
                  'lat'    :np.array(latc),
                  'flags'  :np.array(boundary_flags)}

         if boundary_flags.count(0)>0:
            # does have an ice-water boundary
            Dct.update({'type':'outer'})
         else:
            # doesn't have an ice-water boundary
            Dct.update({'type':'inner'})

         self.MIZclass.append(Dct)
         # ============================================================
      return
      # ============================================================

   # ============================================================
   def get_areas(self):
   
      Nmiz                    = len(self.MPlist_MIZ)
      self.MIZ_areas          = np.zeros((Nmiz))
      self.MIZ_forms_sorted   = Nmiz*[[]]
      for i,Poly in enumerate(self.MPlist_MIZ):
         for n,fdct in self.MIZ_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.MIZ_areas[i]   += float(lut['AREA'])
               self.MIZ_forms_sorted[i].append(n)


      Npack                   = len(self.MPlist_PACK)
      self.PACK_areas         = np.zeros((Npack))
      self.PACK_forms_sorted  = Npack*[[]]
      for i,Poly in enumerate(self.MPlist_PACK):
         for n in self.PACK_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.PACK_areas[i]   += float(lut['AREA'])
               self.PACK_forms_sorted[i].append(n)


      Nwater                 = len(self.MPlist_WTR)
      self.WTR_areas         = np.zeros((Nwater))
      self.WTR_forms_sorted  = Nwater*[[]]
      for i,Poly in enumerate(self.MPlist_WTR):
         for n in self.WTR_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.WTR_areas[i]   += float(lut['AREA'])
               self.WTR_forms_sorted[i].append(n)

      return
   # ============================================================


   # ============================================================
   def test_plot(self,figname=None):

      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      fig      = plt.figure(figsize=(12,12))
      ax       = fig.add_subplot(111)

      # ====================================================
      # limits of figure (not automatic with patches)
      X0 = []
      X1 = []
      Y0 = []
      Y1 = []
      if self.N_WTR>0:
         x0,y0,x1,y1 = self.MP_WTR.bounds
         X0.append(x0)
         X1.append(x1)
         Y0.append(y0)
         Y1.append(y1)

      if self.N_PACK>0:
         x0,y0,x1,y1 = self.MP_PACK.bounds
         X0.append(x0)
         X1.append(x1)
         Y0.append(y0)
         Y1.append(y1)

      if self.N_MIZ>0:
         x0,y0,x1,y1 = self.MP_MIZ.bounds
         X0.append(x0)
         X1.append(x1)
         Y0.append(y0)
         Y1.append(y1)

      ax.set_xlim([np.min(X0),np.max(X1)])
      ax.set_ylim([np.min(Y0),np.max(Y1)])
      ax.set_xlabel('Longitude, $^\circ$E',fontsize=18)
      ax.set_ylabel('Latitude, $^\circ$N',fontsize=18)
      for li in ax.get_xticklabels():
         li.set_fontsize(18)
      for li in ax.get_yticklabels():
         li.set_fontsize(18)
      # ====================================================


      # ====================================================
      Fplt.plot_patches(ax,self.MPlist_WTR,'c')
      Fplt.plot_patches(ax,self.MPlist_PACK,'r')

      # distinguish between "inner" and "outer" MIZ
      ML1   = []
      ML2   = []
      for i,poly in enumerate(self.MPlist_MIZ):
         if self.MIZclass[i]['type']=='outer':
            ML1.append(poly)
         else:
            ML2.append(poly)

      Fplt.plot_patches(ax,ML1,'y')
      Fplt.plot_patches(ax,ML2,'g')
      # ====================================================


      # ====================================================
      # plot ice edges (flag=0) and inner edges (flag=1)
      lstils   = ['.b','.r']
      for MIZi in self.MIZclass:
         Lon   = MIZi['lon']
         Lat   = MIZi['lat']
         ax.plot(Lon,Lat,':')
         for flag in range(2):
            lon   = Lon[MIZi['flags']==flag]
            lat   = Lat[MIZi['flags']==flag]
            ax.plot(lon,lat,lstils[flag])
      # ====================================================


      # ====================================================
      # plot all original MIZ polygons
      for n,fdct in self.MIZ_forms:
         poly     = self.shapes.polygons[n]
         lon,lat  = poly.exterior.coords.xy
         ax.plot(lon,lat,':k')
      # ====================================================


      # ====================================================
      if figname is not None:
         print('Saving '+figname+'...\n')
         fig.savefig(figname,bbox_inches='tight')
      else:
         plt.show(fig)
      # ====================================================

      ax.cla()
      plt.close(fig)
      return
############################################################################
