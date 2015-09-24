import numpy as np
import geometry_planar as GP
import geometry_sphere as GS
import shapely.geometry as shg
from matplotlib import pyplot as plt

######################################################################
def fill_poly(x,y,res=None):
   # if resolution of a poly is too low, increase it artificially
   xy                            = np.array([x,y]).transpose()
   xyc                           = [tuple(xyi) for xyi in xy]
   P,resolution,spacings,th_vec  = GP.curve_info(xyc)

   ############################################################
   if res is None:
      # add a point in between each current one
      x2,y2 = xyc[0]
      x2    = [x2]
      y2    = [y2]
      for i in range(1,len(xyc)):
         x0,y0 = x2[-1],y2[-1]
         x1,y1 = xyc[i]
         x2.extend([.5*(x0+x1),x1])
         y2.extend([.5*(y0+y1),y1])
   else:
      # add points to make spacings <=res
      x2 = []
      y2 = []
      for i,spc in enumerate(spacings):
         x0,y0 = xyc[i]
         x1,y1 = xyc[i+1]
         dist  = np.sqrt(pow(x1-x0,2)+pow(y1-y0,2))
         if dist>res:
            N  = np.ceil(dist/float(res))
            xx = list(np.linspace(x0,x1,num=N))[:-1]
            yy = list(np.linspace(y0,y1,num=N))[:-1]
         else:
            xx = [x0]
            yy = [y0]

         x2.extend(xx)
         y2.extend(yy)

      # include last point
      x2.append(x1)
      y2.append(y1)
   ############################################################

   return np.array(x2),np.array(y2)
#######################################################################

########################################################################################
class MIZ_info:

   ########################################################################################
   def __init__(self,xy_coords,mapper,MIZlines,MIZwidths_int,MIZwidths_tot,func_vals=None):
   
      #################################################################
      x,y                  = np.array(xy_coords).transpose()
      lons,lats            = mapper(x,y,inverse=True)
      self.ll_bdy_coords   = [(lons[i],lats[i]) for i in range(len(lons))]

      # NB use "1*" to remove pointers to the arrays outside the function
      # - like copy, but works for lists also
      self.int_widths   = 1*MIZwidths_int
      self.tot_widths   = 1*MIZwidths_tot
      self.area	      = GS.area_polygon_ellipsoid(lons,lats) # area of polygon
      self.perimeter    = GS.perimeter(lons,lats)              # perimeter of polygon
      self.FDI          = frac_dim_index([self.area,self.perimeter])
      
      # lon-lat info if present
      self.spherical_geometry = True
      self.MIZlines           = 1*MIZlines   # (lon,lat) coordinates of each contour

      if func_vals is not None:
         self.func_vals	 = 1*func_vals	   # value of function used by Laplace's equation

      # some summarising info about "lengths"
      # - intersection widths
      lens                          = np.array(self.int_widths)
      self.int_width_mean           = np.mean(lens)
      self.int_width_median         = np.median(lens)
      self.int_width_percentile05   = np.percentile(lens,5)
      self.int_width_percentile95   = np.percentile(lens,95)

      # - total widths
      lens                          = np.array(self.tot_widths)
      self.tot_width_mean           = np.mean(lens)
      self.tot_width_median         = np.median(lens)
      self.tot_width_percentile05   = np.percentile(lens,5)
      self.tot_width_percentile95   = np.percentile(lens,95)

      # record for shapefile
      self.record = {}
      self.record.update({'Area'                      : self.area})
      self.record.update({'Perimeter'                 : self.perimeter})
      self.record.update({'Fractal_dimension_index'   : self.FDI})
      self.record.update({'Width_mean'                : self.int_width_mean})
      self.record.update({'Width_median'              : self.int_width_median})
      self.record.update({'Width_percentile05'        : self.int_width_percentile05})
      self.record.update({'Width_percentile95'        : self.int_width_percentile95})
      return
      ##############################################################

   #################################################################
   def parts(self):
      # give the "parts" needed by shapefile
      return [1*self.ll_bdy_coords]
   #################################################################

   #################################################################
   def plot_soln(self,bmap,**kwargs):
      for MIZc in self.MIZlines:
         MIZc.plot_lines(bmap,**kwargs)
      return
   #################################################################

   # end class MIZ_info
####################################################################

####################################################################
class comp_mapper:
   # if inverse=True:
   #    composite map: PCA coords -> basemap coords -> lon,lat
   #    (map1=pca_mapper.mapper,map2=basemap)
   # else:
   #    inverse composite: map PCA coords <- basemap coords <- lon,lat
   #    (map1=pca_mapper.mapper,map2=basemap)

   #################################################################
   def __init__(self,map1,map2):
      self.map1   = map1
      self.map2   = map2
      return
   #################################################################

   #################################################################
   def feval(self,X,Y,inverse=True):
      if inverse:
         x,y         = self.map1(X,Y,inverse=True) # eg PCA X,Y -> basemap x,y
         lons,lats   = self.map2(x,y,inverse=True) # eg basemap x,y -> lon,lat
         out         = lons,lats
      else:
         x,y   = self.map2(X,Y,inverse=False) # eg basemap x,y <- lon,lat
         X,Y   = self.map1(x,y,inverse=False) # eg PCA X,Y <- basemap x,y
         out   = X,Y
      return out
   #################################################################
####################################################################

####################################################################
class MIZline:

   #################################################################
   def __init__(self,lons,lats):
      arclen         = GS.arc_length(lons,lats,radians=False,closed=False)
      llc            = np.array([lons,lats]).transpose()
      self.ll_coords = [tuple(llci) for llci in llc]
      self.length    = arclen[-1]
      return
   #################################################################

   #################################################################
   def plot_line(self,bmap,**kwargs):
      lons,lats   = np.array(self.ll_coords).transpose()
      bmap.plot(lons,lats,latlon=True,**kwargs)
      return
   #################################################################

####################################################################

####################################################################
class MIZcont:

   #################################################################
   def __init__(self,LSi,mapper):

      if not hasattr(LSi,'geoms'):

         # LSi is a single line
         self.Nlines = 1
         xv,yv       = LSi.coords.xy
         lons,lats   = mapper(xv,yv,inverse=True)
         #
         self.intersection_length   = GS.perimeter(lons,lats)
         #
         MIZl              = MIZline(lons,lats)
         self.lines        = [MIZl]
         self.total_length = self.intersection_length

      else:
         # LSi is multiple line
         self.Nlines                = len(LSi.geoms)
         self.lines                 = []
         self.intersection_length   = 0
         self.total_length          = 0
         for i,Lsi in enumerate(LSi.geoms):
            xv,yv     = Lsi.coords.xy
            lons,lats = mapper(xv,yv,inverse=True)
            MIZl      = MIZline(lons,lats)
            self.lines.append(MIZl)
            d0 = GS.perimeter(lons,lats)
            #
            self.intersection_length   = self.intersection_length+d0

            # add distance between end of previous line and start of current line
            # to total length
            if i>0:
               x0,y0       = list(LSi.geoms[i-1].coords)[-1] # end of previous line
               x1,y1       = list(LSi.geoms[i].coords)[0]    # start of current line
               lon0,lat0   = mapper(x0,y0,inverse=True)
               lon1,lat1   = mapper(x1,y1,inverse=True)
               lon_ends    = np.array([lon0,lon1])
               lat_ends    = np.array([lat0,lat1])
               d0          = d0+GS.perimeter(lon_ends,lat_ends)

            self.total_length = self.total_length+d0
      return
   #################################################################

   #################################################################
   def plot_lines(self,bmap,**kwargs):
      for MIZl in self.lines:
         MIZl.plot_line(bmap,**kwargs)
      return
   #################################################################

####################################################################

##################################################################
def frac_dim_index(poly):
   # fractal dimension index of a polygon
   # circle ~ .78
   # square ~ 1.
   # increases with complexity (<2 in 2d space)
   if type(poly)==type([]):
      # just a list with P,A
      P,A   = poly
   else:
      # shapely poly
      P  = poly.length
      A  = poly.area

   if A==1:
      A  = 3
      P  = np.sqrt(3)*P
   return 2*np.log(P/4.)/np.log(A)
##################################################################

##################################################################
def covering_polygon(poly):
   # try to simplify shape by dilation
   # - reduce fractal dimension to below a "complexity threshold"
   # - will then be easier to apply Laplacian

   ###################################################################
   test = 0
   def test_plot(poly,fdi,poly0=None,figname=None):
      x,y = np.array(poly.exterior.coords).transpose()
      plt.plot(x,y)

      ss = 'fdi = %f' %(fdi)
      if poly0 is not None:
         x,y = np.array(poly0.exterior.coords).transpose()
         plt.plot(x,y)
         x,y = np.array(poly0.convex_hull.exterior.coords).transpose()
         plt.plot(x,y)
         fdi0  = frac_dim_index(poly0)
         fdi1  = frac_dim_index(poly0.convex_hull)
         ss = ss+' (%f,%f)' %(fdi0,fdi1)
      else:
         fdi1  = frac_dim_index(poly.convex_hull)
         ss    = ss+' (%f)' %(fdi1)
      plt.title(ss)

      if figname is not None:
         plt.savefig(figname)
      plt.show()
   ###################################################################

   fdi      = frac_dim_index(poly)
   thresh   = .5*(fdi+frac_dim_index(poly.convex_hull))
   if test:
      print('fdi = %f, thresh = %f' %(fdi,thresh))
      poly0 = poly
      test_plot(poly,fdi)

   count = 0

   xyc = list(poly.exterior.coords)
   P,res,spacings,th_vec   = GP.curve_info(xyc,closed=True)

   while fdi>thresh and count<10:
      poly  = poly.buffer(res)
      fdi   = frac_dim_index(poly)
      count = count+1

      if test:
         print('fdi = %f, thresh = %f' %(fdi,thresh))
         test_plot(poly,fdi,poly0)

   if test:
      print('fdi = %f, thresh = %f' %(fdi,thresh))
      test_plot(poly,fdi,poly0=poly0,figname='test.png')

   return poly
##################################################################

################################################################################################
class pca_mapper:

   #############################################################################################
   def __init__(self,xy_coords):
      import numpy as np
      #
      x,y   = np.array(xy_coords).transpose()
      #
      self.x0 = np.mean(x)
      self.y0 = np.mean(y)
      self.x  = x
      self.y  = y
      #
      xy_rel   = np.array([x-self.x0,y-self.y0]) # 2xN matrix
      cov      = xy_rel.dot(xy_rel.transpose()) # covariance (2x2 matrix)

      # reorder so 1st eig is bigger (so X represents the major axis)
      evals,evecs    = np.linalg.eig(cov)
      tmp            = sorted([(lam,i) for i,lam in enumerate(evals)],reverse=True)
      self.evals,jj  = np.array(tmp).transpose()
      jj             = [int(j) for j in jj]
      self.evecs     = evecs[:,jj]

      # coords of poly in transformed coords
      self.X,self.Y  = self.mapper(x,y,inverse=False)

      return
   #############################################################################################
      
   #############################################################################################
   def mapper(self,x,y,inverse=False):

      import numpy as np

      if not inverse:
         # in:  basemap coordinates
         # out: coordinates relative to principal components
         xy_rel   = np.array([x-self.x0,y-self.y0]) # 2xN matrix
         X,Y      = self.evecs.transpose().dot(xy_rel)
      else:
         # in:  coordinates relative to principal components
         # out: basemap coordinates
         xy    = np.array([x,y])
         X,Y   = self.evecs.dot(xy)
         X     = X+self.x0
         Y     = Y+self.y0

      return X,Y
   #############################################################################################

   #############################################################################################
   def set_func_vals(self):

      import geometry_planar as GP
      import numpy as np

      Nc       = len(self.X)
      coords   = np.array([self.X,self.Y]).transpose()
      ss       = GP.arc_length(coords,closed=True)
      P        = ss[-1] # perimeter
      ss       = ss[:-1] # drop last val
      #
      nvec     = range(Nc)
      fvals    = 0*ss

      # longest directions
      # - these can be the 2 zeros of a sine function
      Xsort = sorted([(x_,i) for i,x_ in enumerate(self.X)])
      i0    = Xsort[0][1]
      i1    = Xsort[-1][1]

      # orientation doesn't matter
      # - swap indices if inconvenient
      if i1<i0:
         i0,i1 = i1,i0

      # 1st half of polygon
      s_top       = ss[i0:i1]
      ntop        = range(i0,i1)
      L_top       = ss[i1]-ss[i0]
      fvals[ntop] = np.sin((np.pi/L_top)*(s_top-s_top[0]))

      # 2nd half of polygon
      s_bot = list(ss[i1:])
      nbot  = range(i1,Nc)
      s_bot.extend(list(P+ss[:i0]))
      nbot.extend(range(i0))
      L_bot       = ss[i0]+P-ss[i1]
      fvals[nbot] = np.sin(np.pi+(np.pi/L_bot)*(s_bot-s_bot[0]))

      return fvals
   #############################################################################################

   #############################################################################################
   def get_MIZ_lines(self,bmap):
      #
      X0 = self.X.min()
      X1 = self.X.max()
      Y0 = self.Y.min()
      Y1 = self.Y.max()
      #
      xyc   = np.array([self.X,self.Y]).transpose()
      xyc   = [tuple(xy) for xy in xyc]
      shp   = shg.Polygon(xyc)
      #
      P,resolution,spacings,th_vec  = GP.curve_info(xyc)
      #
      ny = 2*int((Y1-Y0)/resolution)
      nx = 2*int((X1-X0)/resolution)
      #
      YY = np.linspace(Y0,Y1,ny)
      XX = np.linspace(X0,X1,nx)
      #
      MIZlines       = []
      MIZwidths_int  = []
      MIZwidths_tot  = []
      mapper         = comp_mapper(self.mapper,bmap).feval

      for xi in XX:
         lin   = shg.LineString([(xi,yi) for yi in YY])
         if lin.intersects(shp):
            LSi   = lin.intersection(shp)
            MIZc  = MIZcont(LSi,mapper)
            #
            MIZlines.append(MIZc)
            MIZwidths_int.append(MIZc.intersection_length)
            MIZwidths_tot.append(MIZc.total_length)

      #TODO also collect MIZwidths into area_info object
      MIZi  = MIZ_info(xyc,mapper,MIZlines,MIZwidths_int,MIZwidths_tot)
      return MIZi
################################################################################################

################################################################################################
def save_shapefile(MIZpolys,filename='test.shp'):

   import shapefile
   w  = shapefile.Writer(shapefile.POLYGON)

   ###############################################################################
   # define attributes
   fields   = MIZpolys[0].record.keys()
   for fld in fields:
      # create field in shapefile
      w.field(fld,'N','40') # name,type ('C'=character, 'N'=number), size (?)
   ###############################################################################

   ###############################################################################
   for MIZi in MIZpolys:
      # add parts:
      parts = MIZi.parts()
      w.poly(parts=parts)

      # add record (dictionary):
      rec   = MIZi.record
      w.record(**rec)
   ###############################################################################

   # save file
   print('\nSaving polygons to shapefile: '+filename)
   w.save(filename)
   return
################################################################################################
