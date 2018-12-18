import os,sys
import numpy as np
import geometry_planar as GP
import geometry_sphere as GS




# MORPH = 'Open','Close' or None
# - 'Close' seems best
# MORPH = 'Close'
MORPH = None

# ======================================================
def augment_contour_ij(cont):
   # msr.contour only returns u,v points
   # - here we add the q points

   # ==================================================
   def get_loc(i,j):
      # determine if u or v point
      if i-np.floor(i)>.25:
         # u points have half-values for i
         return 'u'
      elif j-np.floor(j)>.25:
         # v points have half-values for j
         return 'v'
   # ==================================================

   Cont  = []
   Locs  = []
   N     = len(cont)
   for k in range(N-1):
      #closed curves repeat 1st point - so this works for open and closed

      # get current point
      c0    = cont[k]
      i,j   = c0
      loc   = get_loc(i,j)
      Cont.append(c0)
      
      # get next point
      kp1   = np.mod(k+1,N)
      i1,j1 = cont[kp1]
      loc1  = get_loc(i1,j1)
      if loc=='u' and loc1=='u':
         # vertical average
         i2 = i
         j2 = .5*(j+j1)
      elif loc=='v' and loc1=='v':
         # horizontal average
         j2 = j
         i2 = .5*(i+i1)
      elif loc=='u' and loc1=='v':
         # horizontally aligned with the v point (1),
         # vertically aligned with the u point   (0)
         # (nearest corner)
         i2 = i
         j2 = j1
      elif loc=='v' and loc1=='u':
         # horizontally aligned with the v point (0),
         # vertically aligned with the u point   (1)
         # (nearest corner)
         i2 = i1
         j2 = j

      # next point will always be a "q" (corner) point
      Cont.append((i2,j2))
      Locs.extend([loc,'q'])

   # add last point
   # if closed, this re-closes
   # - else gets the last point
   # TODO - do something more complicated at the array boundary?
   Cont.append((i1,j1))
   Locs.append(loc1)

   return Cont,Locs
# ======================================================


#########################################################
def get_Arctic_regions():
   return ['gre','bar','beau','lab','les']
#########################################################


#########################################################
def read_txt_file_polys(fname):

   # read in text file:
   # each line is:
   # [polygon number]  lon  lat  [function value]
   # 1st line should be a header line

   fid   = open(fname,'r')
   lins  = fid.readlines()
   fid.close()
   
   Polys = []
   llc   = []
   fvals = []
   Pn0   = int(lins[1].split()[0])

   for lin in lins[1:]:
      ss   = lin.split()
      Pno  = int(ss[0])
      lon  = float(ss[1])
      lat  = float(ss[2])
      fval = int(ss[3])

      if Pno==Pn0:
         llc  .append((lon,lat))
         fvals.append(fval)
      else:
         Polys.append([llc,fvals])
         llc     = [(lon,lat)]
         fvals   = [fval]
         Pn0     = Pno

   Polys.append([llc,fvals])
   return Polys
#########################################################


#######################################################################
class poly_info:
   """
   CALL: MIZchar.poly_info(ll_coords,func_vals=None,cdate=None,**kwargs)
   **kwargs: for geometry_sphere.polar_stereographic_simple,NH=True,radius=6371.e3
   make object with helpful info about polygons
   """
   def __init__(self,ll_coords,func_vals=None,cdate=None,**kwargs):

      import rtree.index as Rindex

      if type(cdate)==type('string'):
         self.datetime  = datetime.strptime(cdate,'%Y%m%d')
      else:
         self.datetime  = cdate # already a datetime object (or is None)

      if 0:
         # test coords
         ll_coords = [(0,77.),(1.,77.),(1.,78.),(0.,78.),(0,77.)]

      lon,lat     = np.array(ll_coords).transpose()
      self.length = len(ll_coords)

      # if type(bmap)==type([]):
      #    # already given list of x,y
      #    xy_coords   = 1*bmap
      # else:
      #    # get x,y from basemap
      #    x,y         = bmap(lon,lat)
      #    xy_coords   = list(np.array([x,y]).transpose())

      self.map = GS.polar_stereographic_simple 
      if 'inverse' in kwargs:
         del kwargs['inverse']
      x,y         = self.map(lon,lat,inverse=False,**kwargs)
      xy_coords   = list(np.array([x,y]).transpose())

      if func_vals is not None:
         fvals = list(func_vals)
      else:
         fvals = None
      xyc2  = [tuple(xyp) for xyp in xy_coords] # list of tuples (otherwise logical ops are difficult)
      llc2  = [tuple(llp) for llp in ll_coords] # list of tuples (otherwise logical ops are difficult)

      critter  = (len(ll_coords)!=len(xy_coords))
      flen     = 0
      if func_vals is not None:
         flen     = len(func_vals)
         critter  = critter\
                     or (len(ll_coords)!=len(func_vals))\
                     or (len(xy_coords)!=len(func_vals))

      #####################################################
      if critter:
         if flen>0:
            print('lengths (ll,xy,f):')
            print(len(ll_coords))
            print(len(xy_coords))
            print(flen)
         else:
            print('lengths (ll,xy):')
            print(len(ll_coords))
            print(len(xy_coords))
         raise ValueError('Inconsistent lengths of inputs')
      #####################################################

      # close polygon if necessary
      if xyc2[0]!=xyc2[-1]:
         llc2.append(llc2[0])
         xyc2.append(xyc2[0])
         if func_vals is not None:
            fvals.append(fvals[0])

      x,y   = np.array(xyc2).transpose()
      area  = GP.area_polygon_euclidean(x,y)
      if area<0:
         # reverse order (want anti-clockwise ordering)
         xyc2.reverse()
         llc2.reverse()
         if func_vals is not None:
            fvals.reverse()

      self.area      = abs(area)
      self.func_vals = fvals
      self.ll_coords = llc2
      self.xy_coords = xyc2
      self.perimeter = GP.calc_perimeter(xyc2,closed=True)

      # make rtree index
      idx   = Rindex.Index()
      for i in range(self.length):
         xp,yp = self.xy_coords[i]
         idx.insert(i,(xp,yp,xp,yp)) # a point is a rectangle of zero side-length
      self.index  = idx

      return
   ################################################################


   ################################################################
   def plot(self,pobj=None,latlon=False,show=True,check_flags=False,**kwargs):

      import mod_reading as MR
      from matplotlib import pyplot as plt

      if pobj is None:
         pobj  = MR.plot_object()

      
      if latlon:
         lon,lat  = np.array(self.ll_coords).transpose()
         pobj.ax.plot(lon,lat,**kwargs)
         if check_flags:
            FV = np.array(self.func_vals)
            pobj.ax.plot(lon[FV==1],lat[FV==1],'^',**kwargs)
            pobj.ax.plot(lon[FV==0],lat[FV==0],'o',**kwargs)
      else:
         lon,lat  = np.array(self.ll_coords).transpose()
         x,y      = self.map(lon,lat)
         pobj.ax.plot(x,y,**kwargs)
         if check_flags:
            FV = np.array(self.func_vals)
            pobj.ax.plot(x[FV==1],y[FV==1],'^',**kwargs)
            pobj.ax.plot(x[FV==0],y[FV==0],'o',**kwargs)

      if show:
         plt.show(pobj.fig)

      return pobj
   ################################################################


   ################################################################
   def basemap_plot(self,basemap,show=True,pobj=None,**kwargs):

      # add 'ax' to kwargs if not there already
      if pobj is None:
         if 'ax' not in kwargs:
            pobj  = MR.plot_object()
            kwargs.update({'ax':pobj.ax})
         # else 'ax' in kwargs and don't need to make an axis
      else:
         if 'ax' not in kwargs:
            kwargs.update({'ax':pobj.ax})
         else:
            raise ValueError('Conflicting axes in pobj and kwargs')

      lon,lat  = np.array(self.ll_coords).transpose()
      basemap.plot(lon,lat,latlon=True,**kwargs)

      if show and pobj is not None:
         from matplotlib import pyplot as plt
         import fns_plotting as Fplt
         Fplt.finish_map(basemap)
         plt.show(pobj.fig)

      return pobj
   ################################################################


   ################################################################
   def get_solution(self,METH=5):
      lons,lats   = np.array(self.ll_coords).transpose()
      if METH<2:
         import Laplace_eqn_solution as Leqs

         if METH==0:
            # direct Laplacian soln
            # - use fvals
            print('\nUsing Laplacian on original polygon, with boundary flags\n')
            Psoln = Leqs.get_MIZ_widths(lons,lats,fvals=self.func_vals,basemap=self.map)
         else:
            print('\nUsing Laplacian on original polygon, with PCA\n')
            # - use PCA
            Psoln = Leqs.get_MIZ_widths(lons,lats,basemap=self.map)

      elif METH<4:
         # apply Laplacian method to simpler covering polygon
         import Laplace_eqn_solution as Leqs

         if METH==2:
            method   = 'ConvexHull'
         else:
            method   = 'Buffer'

         print('\nUsing Laplacian on simplified polygon ('+method\
               +'), with PCA\n')
         Psoln = SimplifyPolygon(lons,lats,self.map,method=method)

      elif METH==4:
         print('\nUsing PCA without Laplacian solution\n')
         
         PCA   = pca_mapper(self.xy_coords)
         Psoln = PCA.get_MIZ_lines(self.map)

      elif METH==5:
         #print('\nUsing PCA without Laplacian solution')
         #print('\n - oriented wrt ice edge\n')

         if self.func_vals is None:
            raise ValueError("Need 'func_vals' input to use 'METH' = "+METH)

         subset   = (np.array(self.func_vals)==0)
         if not np.any(subset):
            subset   = None
         # TODO add selector function as extra input?

         # print(self.func_vals)
         # print(subset)
         # print(self.ll_coords)
         # print(self.xy_coords)
         PCA   = pca_mapper(self.xy_coords,subset=subset)
         Psoln = PCA.get_MIZ_lines(self.map)

      return Psoln
   ################################################################
#######################################################################

#################################################################
def get_region_v2(llc=None):
   """
   region = get_region_v2(llc=None)
   if llc is given:
   - output region of polygon (string)
   else:
   - plot the Arctic indicating the different regions
   """

   #########################################
   lat1  = 68.
   lat2  = 79.5
   lat7  = 66.
   lat8  = 72.
   lat9  = 78.2
   #
   lon1  = 25.
   lon2  = 16.
   lon3  = 82.
   lon4  = 140.
   lon5  = -123.
   lon6  = -44.
   lon7  = -63.5
   lon8  = -80.
   ##############################################################


   if llc is not None:
      lonc,latc   = llc

      #########################################
      if latc<0:
         # Antarctic
         region   = 'Antarctic'

      elif lonc<lon5 or lonc>lon4:
         # Beaufort
         region = 'beau'

      elif lonc>lon3:
         # Laptev & East Siberian Sea
         region = 'les'

      elif lonc>lon1:
         # Barents
         region = 'bar'

      elif lonc>lon6:
         # Greenland/Barents/Norway (fjord)
         if lonc<lon2 or latc>lat2:
            # Greenland
            region = 'gre'
         elif latc>lat1:
            # Barents
            region = 'bar'
         else:
            # Baltic Sea
            region = 'balt'

      elif lonc<lon8 or latc>lat9:
         # Canada (Hudson Bay, Canadian archipelago)
         region = 'can'

      elif latc>lat8 or lonc>lon7:
         # Labrador Sea
         region = 'lab'

      elif latc>lat7+(lat8-lat7)/(lon8-lon7)*(lonc-lon7):
         # Labrador Sea
         region = 'lab'

      else:
         # Canada (Hudson Bay, Canadian archipelago)
         region = 'can'

      return region
   else:
      # plot Arctic indicating different regions
      import mod_reading as mr
      import fns_plotting as Fplt
      from matplotlib import pyplot as plt
      bmap  = Fplt.start_HYCOM_map('Arctic')
      pobj  = mr.plot_object()

      # plot OSISAF concentration for 1 day:
      # cdate = '20150801'
      cdate = '20150307'
      ncfil = '/work/shared/nersc/msc/OSI-SAF/2015_nh_polstere/'+\
               'ice_conc_nh_polstere-100_multi_'+cdate+'1200.nc'

      nci      = mr.nc_getinfo(ncfil)
      var_opt  = mr.make_plot_options('fice',conv_fac=.01,ice_mask=True)
      nci.plot_var(var_opt,bmap=bmap,pobj=pobj,clabel='Ice concentration (OSISAF)',show=False)
      
      # plot regions
      fig   = pobj.fig
      ax    = pobj.ax
      ltk   = 2.5
      fsz   = 18
      lat0  = 40


      ############################################################
      # - Greenland Sea
      lats  = np.linspace(lat0,lat1)
      lons  = lon1+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lats  = np.linspace(lat0,lat1)
      lons  = lon2+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lons  = np.linspace(lon1,lon2)
      lats  = lat1+0*lons
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lats  = np.linspace(lat1,lat2)
      lons  = lon2+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lons  = np.linspace(lon2,lon1)
      lats  = lat2+0*lons
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lats  = np.linspace(lat2,90)
      lons  = lon1+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      x,y   = bmap(-8,75)
      ax.text(x,y,'GRE',fontsize=fsz)
      #
      x,y   = bmap(16,62)
      ax.text(x,y,'BALT',color='r',fontsize=fsz)
      ############################################################


      ############################################################
      # - Barents Sea
      lats  = np.linspace(lat0,90)
      lons  = lon3+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      x,y   = bmap(43,78)
      ax.text(x,y,'BAR',fontsize=fsz)
      ############################################################


      ############################################################
      # - Laptev & East Siberian Sea
      lats  = np.linspace(lat0,90)
      lons  = lon4+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      x,y   = bmap(110,70)
      ax.text(x,y,'LES',color='r',fontsize=fsz)
      ############################################################


      ############################################################
      # - Beaufort Sea
      lats  = np.linspace(lat0,90)
      lons  = lon5+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')

      x,y   = bmap(-165,75)
      ax.text(x,y,'BEAU',fontsize=fsz)
      ############################################################


      ############################################################
      # - Labrador Sea
      lats  = np.linspace(lat0,90)
      lons  = lon6+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')

      lats  = np.linspace(lat0,lat7)
      lons  = lon7+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')

      lons  = [lon7,lon8]
      lats  = [lat7,lat8]
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lats  = np.linspace(lat8,lat9)
      lons  = lon8+0*lats
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      lons  = np.linspace(lon8,lon6)
      lats  = lat9+0*lons
      bmap.plot(lons,lats,latlon=True,linewidth=ltk,color='k')
      #
      x,y   = bmap(-61,59.5)
      ax.text(x,y,'LAB',fontsize=fsz)
      #
      x,y   = bmap(-100,65)
      ax.text(x,y,'CAN',color='r',fontsize=fsz)
      ############################################################

      figname  = 'MIZregions.png'
      print('\nSaving to '+figname)
      fig.savefig(figname)
      # fig.show()

      plt.close(fig)
      return
#################################################################

################################################################################
def mask_region(MIZbins,lon,lat,region=None,vertices=None):

   shp         = MIZbins.shape
   sz          = MIZbins.size
   var_name    = MIZbins.var_name
   sitn_type   = MIZbins.sitn_type

   TEST_PLOT   = 0
   if TEST_PLOT:
      from matplotlib import pyplot as plt
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      ax.imshow(MIZbins.icemap,origin='lower')
      fig.show()
   MIZmask     = MIZbins.MIZmask .reshape(sz)
   ICEmask     = MIZbins.ICEmask .reshape(sz)
   PACKmask    = MIZbins.PACKmask.reshape(sz)

   lons  = lon.reshape(sz)
   lats  = lat.reshape(sz)

   if vertices is not None:
      print('\nUsing vertices to restrict MIZ calculation')
      jMIZ     = np.where(MIZmask==1)[0]
      MIZlons  = lons[jMIZ]
      MIZlats  = lats[jMIZ]
      #
      MM = MIZmask[jMIZ]
      IM = MIZmask[jMIZ]
      PM = MIZmask[jMIZ]

      # vertices
      vlons,vlats = np.array(vertices).transpose()
      NH          = (np.mean(vlats)>0)
      vx,vy       = GS.polar_stereographic_simple(vlons,vlats,NH=NH,inverse=False)
      xyverts     = [(vx[i],vy[i]) for i in range(len(vx))]
      # make sure polygon is closed
      if xyverts[0] != xyverts[-1]:
         print('- closing polygon')
         xyverts.append(xyverts[0])

      # simple (fast) stereographic projection
      xMIZ,yMIZ   = GS.polar_stereographic_simple(MIZlons,MIZlats,NH=NH,inverse=False)
      # ax1.plot(xMIZ,yMIZ,'.')
      
      # Mask out points not inside the polygon
      inside   = GP.maskgrid_outside_polygon(xMIZ,yMIZ,xyverts) #true if inside poly
      outside  = np.logical_not(inside)
      # ax1.plot(xMIZ[inside],yMIZ[inside],'or')
      # ax1.plot(vx,vy,'-')

      # if not in region set MIZ=0,ICE,PACK=NaN
      # - effectively this is land, so these boundaries become "unknown"
      MM[outside]    = 0.
      IM[outside]    = np.nan
      PM[outside]    = np.nan
      MIZmask [jMIZ] = MM
      ICEmask [jMIZ] = IM
      PACKmask[jMIZ] = PM

   elif region is None:
      raise ValueError('Need to provide vertices or name of region')

   elif region == 'Arctic':
      for i in np.where(MIZmask==1)[0]:
         # where there is MIZ, check if it's in the region
         # - if not set MIZ=0,ICE,PACK=NaN
         # - effectively this is land, so these boundaries become "unknown"
         llc         = (lons[i],lats[i])
         Preg        = get_region_v2(llc)
         gd_regions  = get_Arctic_regions()
            # discard Baltic Sea,
            # Canadian Archipelago and Hudson Bay,
            # Antarctic
         if Preg not in gd_regions:
            MIZmask [i] = 0
            ICEmask [i] = np.nan
            PACKmask[i] = np.nan
   else:
      for i in np.where(MIZmask==1)[0]:
         # where there is MIZ, check if it's in the region
         # - if not set MIZ=0,ICE,PACK=NaN
         # - effectively this is land, so these boundaries become "unknown"
         llc   = (lons[i],lats[i])
         Preg  = get_region_v2(llc)
         if Preg!=region:
            MIZmask [i] = 0
            ICEmask [i] = np.nan
            PACKmask[i] = np.nan

   MIZmask  = MIZmask .reshape(shp)
   ICEmask  = ICEmask .reshape(shp)
   PACKmask = PACKmask.reshape(shp)

   if TEST_PLOT:
      fig2  = plt.figure()
      ax2   = fig2.add_subplot(1,1,1)
      ax2.imshow(MIZmask,origin='lower')
      plt.show(fig2)

   return MIZbinaries(MIZmask,ICEmask,PACKmask,\
         var_name=var_name,sitn_type=sitn_type,region=region)
################################################################################

################################################################################
def shrink_poly(cont,thresh):
   """
   cont=[(i0,j0),(i1,j1),...]
   """
   import rtree
   import shapely.geometry as shg

   L  = len(cont)
   N  = np.ceil(L/float(thresh)) # no of new polygons
   
   shp0  = shg.Polygon(cont)
   Rind  = rtree.index.Index()
   for id_num,ij in enumerate(cont):
      Rind.insert(id_num, (ij[0],ij[1],ij[0],ij[1]))

   pca   = pca_mapper(cont)
   X,Y   = pca.X,pca.Y
   Xmin  = X.min()
   Xmax  = X.max()
   Ymin  = Y.min()
   Ymax  = Y.max()
   XX    = np.linspace(Xmin,Xmax,N+1)[1:-1]
   YY    = np.linspace(Ymin,Ymax,2)

   xyc   = np.array([X,Y]).transpose()
   xyc   = [tuple(xy) for xy in xyc]
   shp   = shg.Polygon(xyc).buffer(0)

   iLines = []
   for xi in XX:
      lin   = shg.LineString([(xi,yi) for yi in YY])
      LSi   = lin.intersection(shp)

      ##################################################################
      # sort intersection lines
      if not hasattr(LSi,'geoms'):
         # LSi is a single line
         iLines.append(LSi)
      else:
         # TODO this is the potentially problematic step
         Lmax  = 0
         for jj,lin in enumerate(LSi.geoms):
            # LSi is multiple lines
            # - choose longest
            Lj = lin.length
            if Lj>Lmax:
               Lmax  = Lj
               jmax  = jj
         iLines.append(LSi.geoms[jmax])
      ##################################################################


   #################################################
   # go back to (i,j) and use the rtree index
   # to make sure we don't add new points
   iLines_ij   = []
   for iLine in iLines:
      x,y      = np.array(list(iLine.coords)).transpose()
      I,J      = pca.mapper(x,y,inverse=True)
      iLine_ij = []
      for n,i_ in enumerate(I):
         ij1   = (i_,J[n])
         id    = list(Rind.nearest(ij1))[0]
         iLine_ij.append(cont[id])

      iLines_ij.append(shg.LineString(iLine_ij))
   #################################################


   #################################################
   # now split original polygon with boundary.union()
   # - this outputs N smaller polygons (KEEP)
   #   and also the lines in iLines_ij (CHUCK)
   mls   = shg.MultiLineString(iLines_ij)
   shps  = shp0.boundary.union(mls)
   conts = []
   for geom in shps.geoms:
      print(geom)
      if geom not in iLines:
         conts.append(list(geom.boundary.coords))
   #################################################

   return conts
################################################################################


################################################################################
def closing(data):
   from skimage import morphology as morph
   # apply closing to avoid small polynyas and clean up a little
   # - erosion then dilation of MIZ
   kernel = np.ones((3,3),np.uint8)
   data_cl = morph.closing(data,kernel)
   return(data_cl)
################################################################################


################################################################################
def opening(data):
   from skimage import morphology as morph
   # apply opening to avoid small polynyas and clean up a little
   # - dilation then erosion of MIZ
   kernel   = np.ones((3,3),np.uint8)
   data_o   = morph.opening(data,kernel)
   return(data_o)
################################################################################


################################################################################
class MIZbinaries:
   def __init__(self,MIZmask,ICEmask,PACKmask,var_name='dmax',sitn_type='MIZ',region=None):

      self.var_name  = var_name
      self.sitn_type = sitn_type
      self.region    = region

      self.MIZmask   = 1*MIZmask
      self.ICEmask   = 1*ICEmask
      self.PACKmask  = 1*PACKmask
      self.shape     = MIZmask.shape
      self.size      = MIZmask.size
      self.get_icemap()

      if sitn_type=='MIZ':
         self.conversion   = {'Pack':'Pack'  ,'Ice':'Ice'}

      elif sitn_type=='under':
         self.conversion   = {'Pack':'Model' ,'Ice':'Obs'}

      elif sitn_type=='over': 
         self.conversion   = {'Pack':'Obs'   ,'Ice':'Model'}

      else:
         raise ValueError('Unknown sitn_type: '+sitn_type)

      return
   #############################################################


   #############################################################
   def get_icemap(self):
      # gets icemap = 0: water; 1: MIZ; 2: pack; NaN: land
      # - only used for plotting I think

      # NB need to get pack in a consistent way from MIZmask:
      # - closing operation changed MIZ
      #   so initialise all ice to PACK (2), then add MIZ map (1)

      self.icemap = 2*self.ICEmask+0*self.PACKmask # this has either nans or a mask (land/missing vals)
      if 0:
         from matplotlib import pyplot as plt
         print(self.ICEmask)
         fig   = plt.figure()
         ax1   = fig.add_subplot(2,2,1)
         Im    = ax1.imshow(self.ICEmask)
         fig.colorbar(Im)
         ax2   = fig.add_subplot(2,2,2)
         Im    = ax2.imshow(self.PACKmask)
         fig.colorbar(Im)
         ax3   = fig.add_subplot(2,2,3)
         Im    = ax3.imshow(self.MIZmask)
         fig.colorbar(Im)
         plt.show(fig)
      # print(self.MIZmask==1)
      self.icemap[self.MIZmask==1]  = 1                              # no nans, only 0/1

      return
   #############################################################

################################################################################


################################################################################
def compobs2mizmap(DIFF,MODmask,OBSmask,sitn_type,var_name='fice'):
   # convert from over/under estimation to MIZ-like situation,
   # so arrays can be used in the MIZ width algorithms

   if sitn_type=='under':
      # Obs   is like the total ice: OBSmask -> ICEmask
      # Model is like the pack  ice: MODmask -> PACKmask
      # MIZbinaries(self,MIZmask,ICEmask,PACKmask,sitn_type='MIZ',var_name=var_name)
      return MIZbinaries(DIFF,OBSmask,MODmask,sitn_type='under',var_name=var_name)
   
   elif sitn_type=='over': 
      # Model is like the total ice: MODmask -> ICEmask
      # Obs   is like the pack  ice: OBSmask -> PACKmask
      # MIZbinaries(self,MIZmask,ICEmask,PACKmask,sitn_type='MIZ',var_name=var_name)
      return MIZbinaries(DIFF,MODmask,OBSmask,sitn_type='over',var_name=var_name)

   else:
      raise ValueError('Unknown option for sitn_type: '+sitn_type)

   return
################################################################################


################################################################################
def convert2mizmap(MODmask,OBSmask):
   # convert regions to MIZ-type 

   DIFF1 = np.zeros(MODmask.shape) # model overestimates  ice cover
   DIFF2 = np.zeros(MODmask.shape) # model underestimates ice cover
   DIFF  = MODmask-OBSmask
      # Difference MODEL - OSISAF
      # 1   : model=ice, obs=water
      # 0   : both agree
      # -1  : model=water, obs=ice

   # overestimation by model (mod-obs=1)
   # - obs mask is subset of mod mask
   #   so obs->pack and mod->ice
   good                 = np.isfinite(DIFF)
   Data                 = 0*DIFF[good]
   Data[DIFF[good]>0.5] = 1
   DIFF1[good]          = 1*Data

   # underestimation by model (mod-obs=-1)
   # - mod mask is subset of obs mask
   #   so mod->pack and obs->ice
   Data                    = 0*DIFF[good]
   Data[DIFF[good]<-0.5]   = 1
   DIFF2[good]             = 1*Data

   test_plot   = 0
   if test_plot:
      from matplotlib import pyplot as plt
      # plot original binaries
      fig   = plt.figure()
      ax1   = fig.add_subplot(3,2,1)
      IM    = ax1.imshow(DIFF)
      fig.colorbar(IM)
      ax2   = fig.add_subplot(3,2,3)
      IM    = ax2.imshow(DIFF1)
      fig.colorbar(IM)
      ax3   = fig.add_subplot(3,2,5)
      IM    = ax3.imshow(DIFF2)
      fig.colorbar(IM)

   if MORPH=='Open':
      # Opening operation (see def opening)
      # - can retrieve original MIZmask from ICEmask and PACKmask
      # - TODO make icemap consistent with new MIZmask?
      DIFF1 = opening(DIFF1)
      DIFF2 = opening(DIFF2)
   elif MORPH=='Close':
      # Closing operation (see def closing)
      # - can retrieve original MIZmask from ICEmask and PACKmask
      # - TODO make icemap consistent with new MIZmask?
      DIFF1 = closing(DIFF1)
      DIFF2 = closing(DIFF2)

   if test_plot:
      # binaries after closing
      ax4   = fig.add_subplot(3,2,4)
      IM    = ax4.imshow(DIFF1)
      fig.colorbar(IM)
      ax6   = fig.add_subplot(3,2,6)
      IM    = ax6.imshow(DIFF2)
      fig.colorbar(IM)
      plt.show(fig)

   # do conversion to MIZ-like situation
   over  = compobs2mizmap(DIFF1,MODmask,OBSmask,'over')
   under = compobs2mizmap(DIFF2,MODmask,OBSmask,'under')

   return over,under
###############################################################################


################################################################################
# make lon cts to make mean lon calc correct
def cts_lon(lon):
   
   # make lon cts:
   lon2  = np.array(lon)
   for i in range(1,len(lon2)):
      dl = lon2[i]-lon2[i-1]
      if dl>180:
         lon2[i:] = lon2[i:]-360.
      elif dl<-180:
         lon2[i:] = lon2[i:]+360.

   lon_mean = np.mean(lon2)
   while lon_mean>180.:
      lon_mean = lon_mean-360.
   while lon_mean<-180.:
      lon_mean = lon_mean+360.

   return lon2,lon_mean
################################################################################


#############################################################
# output: xy_coords
def ij2xy(cont,X,Y):
   # cont =np.array([(i0,j0),(i1,j1),...]).transpose()
   # so rows are (i,j) pairs
   # changes indexes to x and y
   # (NOTE i and j are inverted -> i = [:,1], j = [:,0])
   # ??
   x     = []
   y     = []
   if type(cont)==type([]):
      # if input list of coordinates,
      # convert to numpy array
      # - coords become rows
      cont  = np.array(cont)

   nx    = cont.shape[0]
   ni,nj = X.shape
   for n in range(nx):
      #########################################
      # bilinear interp to get x,y curves from X,Y arrays
      # -i:
      i_ = cont[n,0]
      if i_<=0:
         i0    = 0
         i1    = 0
         ifac  = 0
      elif i_>=ni-1:
         i0    = ni-1
         i1    = ni-1
         ifac  = 0
      else:
         i0    = int(np.floor(i_))
         i1    = i0+1
         ifac  = i_-i0

      # bilinear interp to get x,y curves from X,Y arrays
      # -j:
      j_ = cont[n,1]
      if j_<=0:
         j0    = 0
         j1    = 0
         jfac  = 0
      elif j_>=nj-1:
         j0    = nj-1
         j1    = nj-1
         jfac  = 0
      else:
         j0    = int(np.floor(j_))
         j1    = j0+1
         jfac  = j_-j0

      x_ = X[i0,j0]\
            +ifac*(X[i1,j0]-X[i0,j0])\
            +jfac*(X[i0,j1]-X[i0,j0])
      y_ = Y[i0,j0]\
            +ifac*(Y[i1,j0]-Y[i0,j0])\
            +jfac*(Y[i0,j1]-Y[i0,j0])
      x.append(x_)
      y.append(y_)
      #########################################

   xy_list = zip(x,y)
   xy_list = np.array(xy_list)
   return(xy_list)
#############################################################


################################################################################
# Baffin Bay reader
def baffin_bay():
   SR    = os.getenv('SWARP_ROUTINES')
   gdir  = SR+'/giacomo_hexagon/MIZwidth/geo_info'
   fil   = (gdir+'/baffin_bay.txt')

   bbay     = open(fil)
   bblon    = [] 
   bblat    = [] 
   bblone   = [] 
   bblate   = [] 
   bblonw   = [] 
   bblatw   = [] 

   for n,en in enumerate(bbay):
      nen = en.split(';')
      lon = float(nen[1])
      lat = float(nen[2])
      bblon.append(lon)
      bblat.append(lat)

   for n,en in enumerate(bblon):
      if n < 1266:
         bblone.append(bblon[n])
         bblate.append(bblat[n])
      else:
         bblonw.append(bblon[n])
         bblatw.append(bblat[n])

   return(bblone,bblate,bblonw,bblatw)
###############################################################################


###############################################################################
def plot_regions_v2():
   get_region_v2() # no input to the function and it makes a plot
   return
###############################################################################


###############################################################################
def array2binaries(data,threshm,threshM,ice_mask,var_name='dmax'):
   # takes in an array and min/max threshholds
   # returns MIZmask,ICEmask, PACKmask

   ICEmask  = np.array(ice_mask,dtype='float')
   MIZmask  = np.zeros(ICEmask.shape,dtype='float')

   # TODO error if ice_mask overlaps nans
   nans  = np.isnan(data)
   stop  = np.logical_and(nans,ice_mask).any()
   if stop:
      raise ValueError('data infinite/nan inside ice_mask')

   # ==================================================
   # MAKE MIZ MASK
   Data  = 1.+0*data[ice_mask] # start 1 where ice

   if threshm is not None:
      # remove data that is below the threshold
      Data[data[ice_mask]<threshm]   = 0.

   if threshM is not None:
      # remove data that is above the threshold
      Data[data[ice_mask]>threshM]  = 0.

   MIZmask[ice_mask] = Data # 0: <threshM (eg water+MIZ); 1: >=threshM (eg pack)
   # ==================================================

   # PACK = ICE-MIZ
   PACKmask       = ICEmask - MIZmask  # 1: MIZ (ice-pack); 0: water,pack,land(nans)
   ICEmask [nans] = np.nan #restore the nans
   PACKmask[nans] = np.nan #restore the nans

   if MORPH=='Open':
      # Opening operation (see def opening)
      # - can retrieve original MIZmask from ICEmask and PACKmask
      # - TODO make icemap consistent with new MIZmask?
      MIZmask  = opening(MIZmask)
   elif MORPH=='Close':
      # Closing operation (see def closing)
      # - can retrieve original MIZmask from ICEmask and PACKmask
      # - TODO make icemap consistent with new MIZmask?
      MIZmask = closing(MIZmask)

   return MIZbinaries(MIZmask,ICEmask,PACKmask,var_name=var_name,sitn_type='MIZ') 
###############################################################################


###############################################################################
def arrays2binary_diff(Zmod,Zobs,threshm):
   # takes in an arrays and min/max threshholds
   # returns MIZmask,ICEmask, PACKmask
   
   # generating the binary file, get the difference, divide into overprediction and
   # underprediction, apply the closing on both and finally re-map the land
   ICEmasks = []

   # get the 2 ice masks as masked arrays
   for Z in [Zmod,Zobs]:
      good     = np.logical_not(Z.mask)
      ICEmask  = np.ones(Z.shape)
      #
      Data                    = 1+0*Z[good] #stop errors from NaNs
      Data[Z[good]<threshm]   = 0
      ICEmask[good]           = 1*Data # 0: <threshm (eg water)    ; 1: >=threshm (eg ice)
      ICEmask[Z.mask]         = np.nan # 0: <threshm (eg water)    ; 1: >=threshm (eg ice)

      ICEmasks.append(1*ICEmask)
      
   over,under  = convert2mizmap(ICEmasks[0],ICEmasks[1])
   return over,under
############################################################################


############################################################################
def get_MIZ_poly(Z,lon,lat,fice=None,var_name='dmax',region=None,vertices=None):
   """
   get_MIZ_poly(ZM,var_name='dmax')
   *ZM is a masked array
   *var_name='dmax','fice','hice'
   """

   if fice is None:
      if var_name=='fice':
         fice  = Z
      else:
         raise ValueError('please pass in "fice" masked array for determination of ice edge')

   fmin           = .15 # min conc
   good           = np.logical_not(fice.mask) # only water/ice
   wtr_mask       = np.copy(good)
   ice_mask       = np.copy(good)
   wtr_mask[good] = (fice.data[good]<fmin)
   ice_mask[good] = (fice.data[good]>=fmin)
   
   # mask out water
   ZM                = np.ma.array(Z.data,mask=np.logical_or(Z.mask,wtr_mask))
   ZM.data[ZM.mask]  = np.nan # make sure masked val's corresp to nan's in data

   if var_name == 'fice':
      # conc MIZ
      thresh_min  = None # already applied with wtr_mask
      thresh_max  = .8
   elif var_name == 'hice':
      # thin ice (compare to SMOS)
      thresh_min  = None # already applied with wtr_mask
      thresh_max  = .45
   elif var_name == 'dmax':
      # FSD MIZ
      thresh_min  = None # already applied with wtr_mask
      thresh_max   = 250.
   elif var_name == 'swh':
      # waves-in-ice MIZ
      thresh_min  = .05
      thresh_max  = None # not needed

   MIZbins  = array2binaries(ZM,thresh_min,thresh_max,ice_mask,var_name=var_name)
      # icemap is 0: water; 1: MIZ; 2: pack; NaN: land

   if vertices is not None:
      MIZbins  = mask_region(MIZbins,lon,lat,vertices=vertices)
   elif region is not None:
      MIZbins  = mask_region(MIZbins,lon,lat,region=region)
   
   return MIZ_poly(MIZbins,lon,lat,region=region)
      # object with contours of MIZ and some methods 
      # - eg plot the binaries, get or write poly stats
############################################################################


############################################################################
def get_AOD_polys(Zmod,Zobs,lon,lat,region=None,vertices=None):
   ZM = np.ma.array(Zmod,mask=1-np.isfinite(Zmod))

   var_name    = 'fice'
   thresh_min  = .15
   over,under  = arrays2binary_diff(ZM,Zobs,thresh_min)
      # over.icemap  is 0: both water; 1: model ice/obs water;    2: both ice; NaN: land
      # under.icemap is 0: both water; 1: obs ice  /model water;  2: both ice; NaN: land

   if (vertices is not None) and (region is not None):
      raise ValueError('Cannot pass in both vertices and region')
   elif vertices is not None:
      over  = mask_region(over,lon,lat,vertices=vertices)
      under = mask_region(under,lon,lat,vertices=vertices)
   elif region is not None:
      over  = mask_region(over,lon,lat,region=region)
      under = mask_region(under,lon,lat,region=region)

   return MIZ_poly(over,lon,lat,region=region),\
          MIZ_poly(under,lon,lat,region=region)
      # objects with contours of AODs and some methods 
      # - eg plot the binaries, get or write poly stats
############################################################################


############################################################################
# This Class will find MIZ contours
class MIZ_poly:

   #########################################################################
   def __init__(self,MIZbins,lon,lat,region=None):

      self.Nx,self.Ny   = MIZbins.shape
      self.MIZbinaries  = MIZbins

      self.lon    = 1*lon
      self.lat    = 1*lat
      self.region = region
      
      self.poly_maker()

      # resolved regional classification by regional mask
      # TODO perhaps should still split large polys?
      # # check polys aren't too big
      # # - so ice edge orientation isn't too variable
      # # - so regional classification makes sense
      # #   (eg one poly can stretch right across one region)
      # self.perim_thresh = 150 #pixels: TP4 this is about 1900km
      # self.shrink_polys()

      return
      # __init__
   #############################################################


   #############################################################
   def poly_maker(self):
      # finding contours from the MIZ map
      # (i,j) coords (not always integers though)
      from skimage import measure as msr
      MIZcont  = msr.find_contours(self.MIZbinaries.MIZmask,.5)
      MIZcont  = sorted(MIZcont, key=len)
      shape    = self.Nx,self.Ny
      #
      self.MIZcont            = []
      self.MIZholes           = []
      self.MIZcont_locations  = []
      self.MIZhole_locations  = []
      for cont in MIZcont[::-1]:
         ii,jj = np.array(cont).transpose()
         A     = GP.area_polygon_euclidean(ii,jj) # area of polygon in pixels (can be <0)

         # we use opposite convention to msr.contour
         # * anticlockwise: external boundary
         #   - reverse order so these have positive areas
         # * clockwise: internal boundary
         #   - reverse order so these have negative areas
         Cont,grid_locs  = augment_contour_ij(cont[::-1])
         if A<0:
            self.MIZcont.append(Cont)
            self.MIZcont_locations.append(grid_locs)
         else:
            self.MIZholes.append(Cont)
            self.MIZhole_locations.append(grid_locs)

      # determine which poly each hole is inside
      self.MIZhole_indices = []
      for jh,hole in enumerate(self.MIZholes):
         for ci,cont in enumerate(self.MIZcont):
            inside   = msr.points_in_poly(hole[0:1],cont)
            if inside[0]:
               self.MIZhole_indices.append(ci)
               break

      return
   #############################################################


   # #############################################################
   # #TODO get this working
   # def shrink_polys(self):
   #    """
   #    Check if the poly's are not too big
   #    - if so, split them into smaller ones
   #    with shrink_poly(cont)
   #    """

   #    conts = []
   #    for cont in self.MIZcont:
   #       # [(i0,j0),(i1,j1),...]
   #       if len(cont)<=self.perim_thresh:
   #          conts.append(cont)
   #       else:
   #          conts.extend(shrink_poly(cont,self.perim_thresh))

   #    self.MIZcont   = sorted(conts,key=len)
   #    return
   # #############################################################
   

   #############################################################
   # output: xy_coords
   def ij2xy(self,cont,X,Y):
      # changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
      return ij2xy(cont,X,Y)
   #############################################################


   #############################################################
   def show_maps(self):
      # simple plotter of the different masks
      from matplotlib import pyplot as plt
      ttl   = ['icemap','MIZmask','ICEmask','PACKmask']
      fig   = plt.figure()

      for i,att in enumerate(ttl):
         arr   = getattr(self.MIZbinaries,att)
         arr2  = np.transpose(arr)
         ny,nx = arr2.shape
         ax    = fig.add_subplot(2,2,i+1)
         im    = ax.imshow(arr2,origin='upper')
         ax.set_title(ttl[i])
         # ax.axis([0,nx,0,ny])
         fig.colorbar(im)
         for cont in self.MIZcont:
            ivec  = cont[:,0]
            jvec  = cont[:,1]
            ax.plot(ivec,jvec)

      fig.show()
      return
   #############################################################


   #############################################################
   def get_poly_stat(self,number=0,**kwargs):
      cont  = self.MIZcont[number]
      kwargs.update({'number':number})

      PS = poly_stat(cont,self.MIZbinaries.icemap,self.lon,self.lat,**kwargs)
      return PS
   #############################################################


   #############################################################
   def write_poly_stats(self,outdir='.',\
         filename_start='MIZpolys',do_sort=True):

      if not os.path.exists(outdir):
         os.mkdir(outdir)

      if self.region is None:
         do_sort  = False

      # header of text file:
      ss = 'polygon    lon    lat    flag\n'
      if do_sort:
         ##########################################################
         # write to regional text files:
         # Nreg        = {'gre':0,'bar':0,'ncb':0,'les':0,'lab':0}
         Nreg        = {'gre':0,'bar':0,'can':0,'beau':0,'lab':0,'les':0,'balt':0,'Antarctic':0,'Arctic':0}
         reg_list    = Nreg.keys()
         tfil_list   = {}
         tfiles      = {}

         # open text files and write headers
         for reg in reg_list:
            fname = outdir+'/'+filename_start+'_'+reg+'.txt'
            tfiles.update({reg:fname})
            #
            tfil_list.update({reg:open(fname,'w')})
            tfil_list[reg].write(ss)
      else:
         reg         = 'all'
         fname       = outdir+'/'+filename_start+'.txt'
         tfiles      = {reg:fname}
         Nreg        = {reg:0}
         tfil_list   = {reg:open(fname,'w')}
         reg_list    = [reg]

         # open text file and write header
         tfil_list[reg].write(ss)

      for num,cont in enumerate(self.MIZcont):
         PS = poly_stat(cont,self.MIZbinaries.icemap,\
                        self.lon,self.lat,number=num)

         if do_sort:
            if self.region=='Arctic':
               # already been sorted to not analyse non-Arctic
               reg   = 'Arctic'
            else:
               reg   = PS.region
         else:
            reg   = 'all'

         Nreg[reg]   = Nreg[reg]+1
         for i in range(PS.cont_len):
            ss = '%d    %10.10f    %10.10f    %d\n'\
               %(Nreg[reg],PS.ll_list[i,0],PS.ll_list[i,1],PS.f_vals[i])
            tfil_list[reg].write(ss)

      for reg in reg_list:
         tfil_list[reg].close()
         if Nreg[reg]==0:
            os.remove(tfiles[reg])
            del tfiles[reg]

      return tfiles
   #############################################################
   
######################################################################


############################################################################
# This Class will analyse contours and adds polygon's info to the daily regional file
# -SINGLE POLYGON (need to loop over cont in MIZcont from MIZ_poly)
class poly_stat:
   def __init__(self,cont,icemap,X,Y,number=None,basemap=None,region=None):
   # def __init__(self,cont,icemap,X,Y,number=None,latlon=True,region=None):
      # basemap is None : X,Y are lon,lat
      # else            : X,Y are projected val's
      # icemap = 0, water; 1,MIZ; 2, pack 

      self.number    = number #get the number
      self.latlon    = (basemap is None)
      self.ij_list   = cont #get the contour
      self.class_def() #calculate class from contour (big,medium etc)

      if 0:
         xy_list  = self.ij2xy(X,Y) #ij -> xy (these can be lon/lat also)
         if self.latlon:
            # already lon,lat coords
            self.ll_list   = 1*xy_list
            lonl,latl      = xy_list.transpose()
            print('HEY')
         else:
            xl             = xy_list[:,0]
            yl             = xy_list[:,1]
            lonl,latl      = basemap(xl,yl,inverse=True)
            self.ll_list   = np.array([lonl,latl]).transpose()
      else:
         if self.latlon:
            # lat,lon given
            # project to x,y coords (lon is discts)
            # - GS.polar_stereographic_simple is faster than a basemap
            NH    = (np.mean(Y)>0)
            shp   = X.shape
            X,Y   = GS.polar_stereographic_simple(X[:],Y[:],\
                        NH=NH,radius=6371.e3,inverse=False)# routine outputs vectors

            # find x,y on boundaries of cont's
            X        = np.reshape(X,shp)
            Y        = np.reshape(Y,shp)
            xy_list  = self.ij2xy(X,Y) #ij -> xy (these can be lon/lat also)
            xl       = xy_list[:,0]    #np array
            yl       = xy_list[:,1]    #np array

            # transform back to lon,lat
            lonl,latl      = GS.polar_stereographic_simple(xl,yl,\
                                 NH=NH,radius=6371.e3,inverse=True)# routine outputs vectors
            self.ll_list   = np.array([lonl,latl]).transpose()
         else:
            # already projected coords (with basemap)
            # - find coords of boundary, project back to lon,lat
            xy_list        = self.ij2xy(X,Y)
            xl             = xy_list[:,0]
            yl             = xy_list[:,1]
            lonl,latl      = basemap(xl,yl,inverse=True)
            self.ll_list   = np.array([lonl,latl]).transpose()

      self.area	        = GS.area_polygon_ellipsoid(lonl,latl)       # area of polygon
      self.perimeter    = GS.perimeter(lonl,latl)                    # perimeter of polygon
      self.FDI          = frac_dim_index([self.area,self.perimeter]) # fractal dimension

      self.icemap = icemap
      f_vals      = self.split_cont() #

      # Defining the region
      lat_list       = 1*self.ll_list[:,1]
      latc           = np.mean(lat_list)
      self.lat_mean  = latc

      lon_list       = 1*self.ll_list[:,0]
      lon_list,lonc  = cts_lon(lon_list)
      self.lon_mean  = lonc

      if region is None:
         # self.get_region()
         self.region = get_region_v2((lonc,latc))
      else:
         self.region = region

      return

   # class definition
   def class_def(self):
      self.cont_len  = len(self.ij_list)
      if self.cont_len > 200:
         self.polygon_class = 'H'
      elif self.cont_len > 100 and self.cont_len <= 200:
         self.polygon_class = 'B'
      elif self.cont_len > 30 and self.cont_len <= 100:
         self.polygon_class = 'M'
      elif self.cont_len <= 30:
         self.polygon_class = 'S'
      return()

      #self.stat_chart(save=True)

   # output: xy_coords
   def ij2xy(self,X,Y):
      # changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
      return ij2xy(self.ij_list,X,Y)
   
   # output: f_vals
   def split_cont(self):
      # this function finds the different contours
      # NOTE if a point has non integer i coordinate is going to be a vertical edge,
      # if has non integer j coordinate is going to be a horizontal edge hence the use
      # of different arounds depending on the indexes of the point
      # NOTE if the polygon is an OVERESTIMATION the contour finding is inverted

      # get pack & total ice masks - NB these have NaNs, since icemap has them
      PACKmask = np.zeros(self.icemap.shape)+np.nan
      ICEmask  = np.zeros(self.icemap.shape)+np.nan
      good     = np.isfinite(self.icemap)

      # pack
      Data                       = 0*self.icemap[good]
      Data[self.icemap[good]==2] = 1
      PACKmask[good]             = 1*Data

      # ice
      Data                       = 0*self.icemap[good]
      Data[self.icemap[good]>=1] = 1
      ICEmask[good]              = 1*Data

      vs         = self.ij_list # list of (i,j) pixel indices
      in_cont    = []
      out_cont   = []
      unk_cont   = []
      func_vals  = []
      func_mod   = 0 # value of func_vals at model ice edge
      func_osi   = 1 # value of func_vals at OSISAF ice edge or MIZ edge
      func_unk   = 2 # value of func_vals at other parts of contour
      
      #################################################
      # loop over the contours ("vs" is (i,j) tuple)
      for n,el in enumerate(vs):
         #getting all the neighbours
         around2 = ((int(el[0]+.5),int(el[1]   )),(int(el[0]-.5),int(el[1]   ))) # vertical boundaries - OK!
         around1 = ((int(el[0]   ),int(el[1]+.5)),(int(el[0]   ),int(el[1]-.5))) # horizontal boundaries
         check_cont = 0
         if el[0]/int(el[0]) == 1:
            for h,v in around1:
               if PACKmask[h][v] == ICEmask[h][v] == 1:
                  in_cont.append(el)
                  func_val=func_mod
                  check_cont = 1
               elif PACKmask[h][v] == ICEmask[h][v] == 0:
                  out_cont.append(el)
                  func_val=func_osi
                  check_cont = 1
            if check_cont == 0:
                unk_cont.append(el)
                func_val=func_unk
            func_vals.append(func_val)
         else:
            for h,v in around2:
               if PACKmask[h][v] == ICEmask[h][v] == 1:
                  in_cont.append(el)
                  func_val=func_mod
                  check_cont = 1
               elif PACKmask[h][v] == ICEmask[h][v] == 0:
                  out_cont.append(el)
                  func_val=func_osi
                  check_cont = 1
            if check_cont == 0:
               unk_cont.append(el)
               func_val=func_unk
            func_vals.append(func_val)
      #################################################
      
      in_cont           = np.array(in_cont)
      out_cont          = np.array(out_cont)
      unk_cont          = np.array(unk_cont)
      self.f_vals       = np.array(func_vals)
      self.in_ice_edge  = in_cont
      self.out_ice_edge = out_cont
      self.unknown_edge = unk_cont
      return


   #########################################################
   # def get_region(self):
   #    # output: lon,lat,region of polygon
   #    # NOTE cont_list might be XY or IJ (no LONLAT)

   #    lat_list = 1*self.ll_list[:,1]
   #    latc     = np.mean(lat_list)

   #    lon_list       = 1*self.ll_list[:,0]
   #    lon_list,lonc  = cts_lon(lon_list)

   #    # Defining the region
   #    self.lon_mean,self.lat_mean   = lonc,latc

   #    # read in Baffin Bay outline
   #    bblone,bblate,bblonw,bblatw   = baffin_bay()

   #    if lonc < 90 and lonc > 16:
   #       # barents
   #       self.region = 'bar'
   #    elif lonc <= 16 and lonc > -44:
   #       # greenland
   #       self.region = 'gre'
   #    elif lonc <= -44 and lonc > -90:
   #       n = 0
   #       N = len(bblatw[:-1])
   #       if latc <= bblatw[-1] and lonc >= bblonw[-1]:
   #          self.region = 'lab'
   #       else:
   #          while n < N:
   #             if latc >= bblatw[n+1] and latc <= bblatw[n]:
   #                if lonc >= bblonw[n]: 
   #                   self.region = 'lab'
   #                   break
   #             n += 1
   #          else:
   #             # north canada/beaufort
   #             self.region = 'ncb'
   #    elif lonc < 180 and lonc > 90:
   #       # laptev/east siberian
   #       self.region = 'les'
   #    else:
   #       # north canada/beaufort
   #       self.region = 'ncb'
   #    
   #    return
   # #########################################################


   #########################################################
   # Function stat chart
   def stat_chart(self,save=False,METH=5):
      import fns_plotting as Fplt
      from matplotlib import pyplot as plt
      import matplotlib.gridspec as gridspec
      import matplotlib.lines    as mlines
      # The Statistical Chart is an output that tries to compress as many informations
      # as possible in a single figure.
      # The figure is:
      # 1) Top left arctic map with studied polygon highlighted
      # 2) Top right a close up of the polygon with contours and min. distances
      # 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
      # 4) Bottom right is a legend for the small image
      
      nmbr        = self.number
      pname       = 'Polygon_'+str(nmbr)
      self.name   = pname
      pclass      = self.polygon_class
      region      = self.region
      print('')
      print('Statistic Chart for ', pname)
      print('Class and Region ', pclass, region)
      # DN       = np.flipud(np.transpose(self.icemap))
      DN       = np.transpose(self.icemap)
      Ny,Nx    = DN.shape
      ij       = self.ij_list
      clon     = '%1.2f' %self.lon_mean
      clat     = '%1.2f' %self.lat_mean
      clonlat  = '{0}/{1}'.format(clon,clat)
      #
      inside_contour    = self.in_ice_edge
      outside_contour   = self.out_ice_edge
      unknown_contour   = self.unknown_edge

      area  = '%d' %(int(self.area/1.e6     ))  # convert to km^2
      perim = '%d' %(int(self.perimeter/1.e3))  # convert to km
      
      ##############################################################################
      # run MIZ width analysis (selective PCA method)
      bmap        = Fplt.start_HYCOM_map(self.region)
      Poly        = poly_info(self.ll_list,bmap,func_vals=self.f_vals,cdate=None)
      lons,lats   = np.array(Poly.ll_coords).transpose()
      if METH<2:
         if METH==0:
            # direct Laplacian soln
            # - use fvals
            print('\nUsing Laplacian on original polygon, with boundary flags\n')
            Psoln = Leqs.get_MIZ_widths(lons,lats,fvals=Poly.func_vals,basemap=bmap)
         else:
            print('\nUsing Laplacian on original polygon, with PCA\n')
            # - use PCA
            Psoln = Leqs.get_MIZ_widths(lons,lats,basemap=bmap)

      elif METH<4:
         # apply Laplacian method to simpler covering polygon
         if METH==2:
            method   = 'ConvexHull'
         else:
            method   = 'Buffer'

         print('\nUsing Laplacian on simplified polygon ('+method\
               +'), with PCA\n')
         Psoln = SimplifyPolygon(lons,lats,bmap,method=method)

      elif METH==4:
         print('\nUsing PCA without Laplacian solution\n')
         
         PCA      = pca_mapper(Poly.xy_coords)
         MIZinfo  = PCA.get_MIZ_lines(bmap)

      elif METH==5:
         print('\nUsing PCA without Laplacian solution')
         print('\n - oriented wrt ice edge\n')

         subset   = (np.array(Poly.func_vals)==0)
         if not np.any(subset):
            subset   = None

         PCA      = pca_mapper(Poly.xy_coords,subset=subset)
         MIZinfo  = PCA.get_MIZ_lines(bmap)
      ##############################################################################
      
      # Setting up the plot (2x2) and subplots
      fig   = plt.figure(figsize=(15,10))
      gs    = gridspec.GridSpec(2,2,width_ratios=[2,1],height_ratios=[4,2])
      plt.suptitle(pname+', class '+pclass+', '+region,fontsize=18)
      main  = plt.subplot(gs[0,0])
      polyf = plt.subplot(gs[0,1])
      tab   = plt.subplot(gs[1,0])
      leg   = plt.subplot(gs[1,1])
      tab.set_xticks([])
      leg.set_xticks([])
      tab.set_yticks([])
      leg.set_yticks([])
      tab.set_frame_on(False)
      leg.set_frame_on(False)
      
      # Main image on the top left
      main.imshow(DN)
      x1 = np.max([0,  np.min(ij[:,0])-10])
      x2 = np.min([Nx, np.max(ij[:,0])+10])
      y1 = np.max([0,  np.min(ij[:,1])-10])
      y2 = np.min([Ny, np.max(ij[:,1])+10])
      #
      Ymin  = 0
      Ymax  = Ny
      Y2    = (y2-Ymin)/(Ymax-Ymin)
      Y1    = (y1-Ymin)/(Ymax-Ymin)

      main.axvspan(x1,x2,ymin=Y1,ymax=Y2,\
         color='red',alpha=0.3)
      main.axis([0,Nx,Ymin,Ymax])
      
      # Polygon image on the top right
      polyf.imshow(DN)
      polyf.axis([x1,x2,y1,y2])
      if len(inside_contour) != 0:
         polyf.plot(inside_contour[:,0],inside_contour[:,1],'ro',markersize=4)
      if len(outside_contour) != 0:
         polyf.plot(outside_contour[:,0],outside_contour[:,1],'yo',markersize=4)
      if len(unknown_contour) != 0:
         polyf.plot(unknown_contour[:,0],unknown_contour[:,1],'go',markersize=4)
      
      # Legend on the bottom right
      mc = mlines.Line2D([],[],color='red',marker='o')      # inside:  MIZ-pack
      oc = mlines.Line2D([],[],color='yellow',marker='o')   # outside: Ice edge
      uc = mlines.Line2D([],[],color='green',marker='o')    # unknown
      leg.legend([mc,oc,uc],('MIZ-Pack Edge','Ice Edge','Unknown Cont.'))
               
      # Statistics text on the bottom left
      txt =  '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
             '2) Area = '+str(area)+' km^2\n'+ \
             '3) Perimeter = '+str(perim)+' km\n'+ \
             '4) Avg. Width = '+'\n'+ \
             '5) Widths SD = '
      tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
      if save:
         outdir = str(out_dir)
         if not os.path.exists(outdir):
             os.mkdir(outdir)
         valid_class = outdir+'/'+dadate
         if not os.path.exists(valid_class):
             os.mkdir(valid_class)
         if not os.path.exists(valid_class+'/'+region):
             os.mkdir(valid_class+'/'+region)
         fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
         plt.close()
      else:
         plt.show(False)
      print('Statistic chart done for '+str(pname))
      print('')
      return
####################################################################


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
   def __init__(self,xy_coords,mapper,MIZlines,func_vals=None):
   
      #################################################################
      x,y                  = np.array(xy_coords).transpose()
      lons,lats            = mapper(x,y,inverse=True)
      self.ll_bdy_coords   = [(lons[i],lats[i]) for i in range(len(lons))]

      # NB use "1*" to remove pointers to the arrays outside the function
      # - like copy, but works for lists also
      self.area	     = GS.area_polygon_ellipsoid(lons,lats)       # area of polygon
      self.perimeter = GS.perimeter(lons,lats,)                   # perimeter of polygon
      self.FDI       = frac_dim_index([self.area,self.perimeter]) # fractal dimension

      self.int_widths   = []
      self.tot_widths   = []
      for MIZc in MIZlines:
         self.int_widths.append(MIZc.intersection_length)
         self.tot_widths.append(MIZc.total_length)
      
      # lon-lat info if present
      self.spherical_geometry = True
      self.MIZlines           = 1*MIZlines   # (lon,lat) coordinates of each contour

      if func_vals is not None:
         self.func_vals	 = 1*func_vals	   # value of function used by Laplace's equation
      else:
         self.func_vals = None

      # some summarising info about "lengths"
      # - intersection widths
      lens                          = np.array(self.int_widths)
      self.int_width_mean           = np.mean(lens)
      self.int_width_median         = np.median(lens)
      self.int_width_percentile05   = np.percentile(lens,5)
      self.int_width_percentile95   = np.percentile(lens,95)
      self.int_width_max            = np.max(lens)

      # - total widths
      lens                          = np.array(self.tot_widths)
      self.tot_width_mean           = np.mean(lens)
      self.tot_width_median         = np.median(lens)
      self.tot_width_percentile05   = np.percentile(lens,5)
      self.tot_width_percentile95   = np.percentile(lens,95)
      self.tot_width_max            = np.max(lens)

      # record for shapefile
      self.record = {}
      self.record.update({'Area'                      : self.area})
      self.record.update({'Perimeter'                 : self.perimeter})
      self.record.update({'Fractal_dimension_index'   : self.FDI})
      self.record.update({'Width_mean'                : self.int_width_mean})
      self.record.update({'Width_median'              : self.int_width_median})
      self.record.update({'Width_percentile05'        : self.int_width_percentile05})
      self.record.update({'Width_percentile95'        : self.int_width_percentile95})
      self.record.update({'Width_max'                 : self.int_width_max})

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

   #################################################################
   def plot_representative_lines(self,bmap,plot_type='Width_mean',**kwargs):
      # locate representative curves for plotting

      Wavg  = None
      if plot_type=='Width_mean':
         Wavg  = self.record[plot_type]
      elif plot_type=='Width95':
         for key in self.record.keys():
            if 'idth_percentile95' in key:
               Wavg  = self.record[key]
               break

         if Wavg is None:
            print(self.record)
            raise ValueError('95th percentile not found in self.record')

      if Wavg is not None:
         count = 0
         for i,MIZc in enumerate(self.MIZlines):
            diff  = abs(Wavg-MIZc.intersection_length)/Wavg
            if (diff<.05) and (MIZc.Nlines==1):
               MIZc.plot_lines(bmap,**kwargs)
               count = count+1

      # if count==0:
      #    # TODO find median line

      return
   #################################################################

   #################################################################
   def bbox(self,bmap):
      lon,lat  = np.array(self.ll_bdy_coords).transpose()
      x,y      = bmap(lon,lat)
      return [x.min(),x.max(),y.min(),y.max()]
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
         x,y         = np.array(x),np.array(y)
         lons,lats   = self.map2(x,y,inverse=True) # eg basemap x,y -> lon,lat
         out         = lons,lats
      else:
         X,Y   = np.array(X),np.array(Y)
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
         xv,yv       = np.array(xv),np.array(yv)
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
            xv,yv     = np.array(xv),np.array(yv)
            lons,lats = mapper(xv,yv,inverse=True)
            MIZl      = MIZline(lons,lats)
            self.lines.append(MIZl)

            d0                         = GS.perimeter(lons,lats)
            self.intersection_length   = self.intersection_length+d0

            # add distance between end of previous line and start of current line
            # to total length
            if i>0:
               x0,y0     = list(LSi.geoms[i-1].coords)[-1] # end of previous line
               x1,y1     = list(LSi.geoms[i].coords)[0]    # start of current line
               x0,y0     = np.array(x0),np.array(y0)
               x1,y1     = np.array(x1),np.array(y1)
               lon0,lat0 = mapper(x0,y0,inverse=True)
               lon1,lat1 = mapper(x1,y1,inverse=True)
               lon_ends  = np.array([lon0,lon1])
               lat_ends  = np.array([lat0,lat1])
               d0        = d0+GS.perimeter(lon_ends,lat_ends)

            self.total_length = self.total_length+d0
      return
   #################################################################

   #################################################################
   def plot_lines(self,bmap,show=False,**kwargs):
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
      from matplotlib import pyplot as plt
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
   def __init__(self,xy_coords,subset=None):
      import numpy as np
      #
      x,y   = np.array(xy_coords).transpose()
      # print(subset)

      if subset is None or not np.any(subset):
         # use all points for the PCA analysis
         self.x0  = np.mean(x)
         self.y0  = np.mean(y)
         xy_rel   = np.array([x-self.x0,y-self.y0]) # 2xN matrix
      else:
         x_       = x[subset]
         y_       = y[subset]
         self.x0  = np.mean(x_)
         self.y0  = np.mean(y_)
         xy_rel   = np.array([x_-self.x0,y_-self.y0]) # 2xN matrix

      self.x   = x
      self.y   = y
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
      
   #######################################################################################
   def mapper(self,x,y,inverse=False):
      """
      *mapper to change from (x,y) to (x',y')
      - rotated coordinates relative to principal components
      *use inverse=True to go from (x',y') to (x,y)
      """

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
   #######################################################################################


   #######################################################################################
   def set_func_vals(self):
      """
      Set func vals using of arc length
      - zeros go at longest ends
      """

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

   #######################################################################################
   def get_MIZ_lines(self,bmap):
      """
      gets MIZ lines for a polygon
      """
      import shapely.geometry as shg
      #
      X0 = self.X.min()
      X1 = self.X.max()
      Y0 = self.Y.min()
      Y1 = self.Y.max()
      #
      xyc   = np.array([self.X,self.Y]).transpose()
      xyc   = [tuple(xy) for xy in xyc]
      # print(xyc)
      shp   = shg.Polygon(xyc).buffer(0)
      #
      P,resolution,spacings,th_vec  = GP.curve_info(xyc)


      #
      ny = np.max([2,2*int((Y1-Y0)/resolution)])
      nx = np.max([2,2*int((X1-X0)/resolution)])
      #
      YY = np.linspace(Y0,Y1,ny)
      XX = np.linspace(X0,X1,nx)
      #
      MIZlines       = []
      MIZwidths_int  = []
      MIZwidths_tot  = []
      mapper         = comp_mapper(self.mapper,bmap).feval

      # print(nx,ny)
      # print(X0,X1,(X1-X0),resolution)
      # print(Y0,Y1,(Y1-Y0),resolution)
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
      MIZi  = MIZ_info(xyc,mapper,MIZlines,func_vals=None)
         # PCA is geometric only, although by using a subset (eg ice edge)
         # we can use boundary info
      return MIZi
################################################################################################


################################################################################################
class MIZ_info_Lap:

   #################################################################
   def __init__(self,Psoln,MIZlines,lonlat):
      self.Laplacian_soln  = Psoln
      #
      lons,lats            = lonlat
      self.ll_bdy_coords   = [(lons[i],lats[i]) for i in range(len(lons))]
      #
      self.area      = GS.area_polygon_ellipsoid(lons,lats)
      self.perimeter = GS.perimeter(lons,lats,closed=True)
      self.FDI       = frac_dim_index([self.perimeter,self.area])
      #
      self.MIZlines        = 1*MIZlines
      self.int_widths   = []
      self.tot_widths   = []
      for MIZc in MIZlines:
         self.int_widths.append(MIZc.intersection_length)
         self.tot_widths.append(MIZc.total_length)

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
   #################################################################

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

   #################################################################
   def plot_representative_lines(self,bmap,plot_type='Width_mean',**kwargs):
      # locate representative curves for plotting
      Wavg  = self.record[plot_type]
      count = 0
      for i,MIZc in enumerate(self.MIZlines):
         diff  = abs(Wavg-MIZc.intersection_length)/Wavg
         if (diff<.05) and (MIZc.Nlines==1):
            MIZc.plot_lines(bmap,**kwargs)
            count = count+1

      # if count==0:
      #    # TODO find median line

      return
   #################################################################

   #################################################################
   def bbox(self,bmap):
      lon,lat  = np.array(self.ll_bdy_coords).transpose()
      x,y      = bmap(lon,lat)
      return [x.min(),x.max(),y.min(),y.max()]
   #################################################################

   # end class MIZ_info_Lap
################################################################################################

#######################################################################
def SimplifyPolygon(lons,lats,bmap,res=10000.,method='ConvexHull'):
   import shapely.geometry       as shg
   import Laplace_eqn_solution   as Leqs
   import geometry_sphere        as GS

   x,y   = bmap(lons,lats)
   xy    = np.array([x,y]).transpose()
   xyc   = [tuple(xyi) for xyi in xy]
   shp   = shg.Polygon(xyc).buffer(0)

   area        = GS.area_polygon_ellipsoid(lons,lats)
   perimeter   = GS.perimeter(lons,lats,closed=True)
   FDI         = frac_dim_index([perimeter,area])

   if method=='ConvexHull':
      # get convex hull
      shp2  = shp.convex_hull
      x2,y2 = shp2.exterior.coords.xy
   else:
      # get convex hull
      shp2  = covering_polygon(shp)
      x2,y2 = shp2.exterior.coords.xy

   # increase resolution (m) (this increases the number of points):
   x3,y3       = fill_poly(x2,y2,res=res)
   lons2,lats2 = bmap(x3,y3,inverse=True)
   
   # apply Laplacian soln to simplified polygon
   Lsoln = Leqs.get_MIZ_widths(lons2,lats2,basemap=bmap)

   ####################################################################
   # restrict contour lines to within original poly
   MIZlines = []
   for llc in Lsoln.area_info.lonlat_contours:
      lonv,latv   = np.array(llc).transpose()
      xx,yy       = bmap(lonv,latv)
      xyv         = np.array([xx,yy]).transpose()
      xyv         = [tuple(xyi) for xyi in xyv]
      #
      LS = shg.LineString(xyv)
      if LS.intersects(shp):
         LSi   = LS.intersection(shp)
         MIZlines.append(MIZcont(LSi,bmap))
   ####################################################################
   
   Psoln = MIZ_info_Lap(Lsoln,MIZlines,[lons,lats])
   return Psoln
#######################################################################


################################################################################################
def save_shapefile(MIZpolys,filename='test.shp'):

   import shapefile
   w  = shapefile.Writer(shapefile.POLYGON)

   ###############################################################################
   # define attributes
   fields   = MIZpolys.MIZ_info_objects[0].record.keys()
   for fld in fields:
      # create field in shapefile
      w.field(fld,'N','40') # name,type ('C'=character, 'N'=number), size (?)
   ###############################################################################

   ###############################################################################
   for MIZi in MIZpolys.MIZ_info_objects:
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


################################################################################################
def get_summary(MIZpolys,wt_meth='A'):

   # =============================================================================
   # Initialise outputs:

   # short names -> long names
   long_names  = {}
   long_names.update({'int_width_mean':'Mean_intersecting_width'})
   long_names.update({'tot_width_mean':'Mean_total_width'})
   long_names.update({'int_width_max' :'Maximum_intersecting_width'})
   long_names.update({'tot_width_max' :'Maximum_total_width'})
   long_names.update({'tot_perim'     :'Total_perimeter'})
   long_names.update({'tot_area'      :'Total_area'})

   # short names -> values
   summ  = {}
   for key in long_names.keys():
      summ.update({key:0})
   # =============================================================================


   # =============================================================================
   # Do analysis
   tot_perim   = 0.
   tot_area    = 0.
   for mp in MIZpolys.MIZ_info_objects:
      P           = mp.perimeter
      A           = mp.area
      tot_perim  += P
      tot_area   += A

      if wt_meth=='P':
         # use polygon perimeter as weight
         # TODO ice edge perimeter?
         wt = P
      elif wt_meth=='A':
         # use polygon area as weight
         wt = A

      # mean widths
      summ['int_width_mean'] += wt*mp.int_width_mean #sum of means
      summ['tot_width_mean'] += wt*mp.tot_width_mean #sum of means

      # max widths
      if 0:
         print(mp.int_width_percentile95/1.e3)
         print(mp.int_width_max/1.e3)
         print(mp.tot_width_percentile95/1.e3)
         print(mp.tot_width_max/1.e3)
      summ['int_width_max']   = np.max([summ['int_width_max'],mp.int_width_max])
      summ['tot_width_max']   = np.max([summ['tot_width_max'],mp.tot_width_max])

   summ['tot_perim'] = tot_perim
   summ['tot_area']  = tot_area
   if len(MIZpolys.MIZ_info_objects)>0:
      if wt_meth=='P':
         wt = float(tot_perim)
      else:
         wt = float(tot_area)

      summ['int_width_mean'] /= wt #divide sum of means by total perim/area to get mean of means
      summ['tot_width_mean'] /= wt #divide sum of means by total perim/area to get mean of means

   return summ,long_names
# ===========================================================================================


# ===========================================================================================
def save_summary(MIZpolys,filename,date_time=None,wt_meth='A'):

   # get summary info
   R,S   = get_summary(MIZpolys,wt_meth=wt_meth)

   # write to file
   print('\nWriting summary to '+filename+'...\n')

   w  = open(filename,'w')
   if date_time is not None:
      if type(date_time) == type('hey'):
         w.write('Date : '+date_time+'\n')
      else:
         w.write('Date : '+date_time.strftime('%Y%m%dT%H%M%SZ')+'\n')

   for fld in ['tot_perim','tot_area','int_width_mean','tot_width_mean','int_width_max','tot_width_max']:
      w.write(S[fld]+' : %16.0f\n' %(R[fld]))

   w.close()

   return R,S
# ===========================================================================================


#########################################################
def single_file(filename,cdate=None):
   """
   CALL: single_file(filename,cdate=None)
   *Read in text file "filename":
   >each line is:
   [polygon number] lon  lat  [function value]
   >First line is header

   returns a MIZchar.poly_info_list object
   """

   ############################################################
   # get polys as "poly_info" objects
   Pols  = read_txt_file_polys(filename)

   poly_infos = []
   for llc,fvals in Pols:
      Poly  = poly_info(llc,func_vals=fvals,cdate=cdate)
      poly_infos.append(Poly)

   # turn the list into a poly_info_list object
   return poly_info_list(poly_infos,cdate=cdate)
#########################################################


#########################################################
class MIZ_info_list:

   ######################################################
   def __init__(self,MIZ_infos,METH=5,cdate=None):
      self.METH               = METH
      self.MIZ_info_objects   = []
      for mi in MIZ_infos:
         self.MIZ_info_objects.append(mi)
      return
   ######################################################


   ######################################################
   def plot_solutions(self,basemap,pobj=None,figname=None):
      """
      MIZchar.MIZ_info_list.plot_solutions(basemap,pobj=None,figname=None)
      """
      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      bmap  = basemap
      if pobj is None:
         pobj  = Fplt.plot_object()
      fig,ax1  = pobj.fig,pobj.ax

      cbar  = False # colorbar not made yet
      for Psoln in self.MIZ_info_objects:

         # plot poly boundary
         lons,lats   = np.array(Psoln.ll_bdy_coords).transpose()
         bmap.plot(lons,lats,'k',latlon=True,linewidth=2,ax=ax1)

         #########################################################
         if self.METH<2:
            # soln from Laplacian method on original polygon
            Psoln.plot_soln(pobj=[fig,ax1],bmap=bmap,cbar=(not cbar),show=False)
            cbar  = True
         #########################################################
         

         #########################################################
         elif self.METH<4:
            # soln from Laplacian method on simpler covering polygon
            # plot Laplacian solution
            Psoln.Laplacian_soln.plot_soln(pobj=[fig,ax1],bmap=bmap,cbar=(not cbar),show=False)
            cbar  = True
            #
            for MIZc in Psoln.MIZlines:
               MIZc.plot_lines(bmap,ax=ax1,color='c')
         #########################################################


         #########################################################
         elif self.METH<6:
            # PCA without Laplacian solution
            # - 4: orientated wrt major axis of "ellipse"
            # - 5: oriented wrt ice edge
            Psoln.plot_soln(bmap,ax=ax1,color='c')
         #########################################################


      #########################################################
      # finish map and close before returning
      Fplt.finish_map(bmap)

      if figname is None:
         # fig.show()
         plt.show(fig)
      else:
         fig.savefig(figname)
         ax1.cla()
         plt.close(fig)
      #########################################################

         
      return pobj
   ######################################################
   
   
   ######################################################
   def save_shapefile(self,filename):
      save_shapefile(self,filename)
      return
   ######################################################


   ######################################################
   def save_summary(self,filename,**kwargs):
      """
      MIZpolys.save_summary(filename,**kwargs)
      returns MIZchar.save_summary(MIZpolys,filename,**kwargs)
      """
      return save_summary(self,filename,**kwargs)
   ######################################################


   ######################################################
   def get_summary(self,**kwargs):
      """
      MIZpolys.get_summary(**kwargs)
      returns MIZchar.get_summary(MIZpolys,**kwargs)
      """
      return get_summary(self,**kwargs)
   ######################################################

#########################################################


#########################################################
class poly_info_list:
   def __init__(self,poly_infos,cdate=None):
      """
      CALL: poly_info_list(filename,cdate=None)
      0     : direct Laplacian with specified boundary flags
      1     : direct Laplacian with boundary flags determined from PCA
      2,3   : direct Laplacian to simplified polygon (lower fractal dimension index),
               with boundary flags determined from PCA for new shape
       * 2 > get convex hull
       * 3 > dilation to get less complicated shape (in between original and convex hull)
      4     : Use PCA to determine direction to get the width in,
               then just take straight lines across (in stereographic projection space)
      5     : Use PCA to determine direction to get the width in,
               oriented wrt the ice edge
               then just take straight lines across (in stereographic projection space)

      *Read in text file "filename":
      >each line is:
      [polygon number] lon  lat  [function value]
      >First line is header
      """

      import fns_Stefan_Maps as FSM


      ############################################################
      # get polys as "poly_info" objects
      self.poly_info_objects  = 1*poly_infos
      self.cdate              = cdate
      

      # make a local stereographic projection
      from pyproj import Proj

      # start with NP stere
      map   = Proj(proj='stere',ellps='WGS84',lon_0=-45,lat_0=90,lat_ts=70)
      X,Y   = [],[]
      for pio in self.poly_info_objects:
         lons,lats   = np.array(pio.ll_coords).transpose()
         x,y         = map(lons,lats,inverse=False)#,NH=True)
         X.extend(list(x))
         Y.extend(list(y))

      # make a new proj centered around the mean coords
      xav         = np.mean(X)
      yav         = np.mean(Y)
      lon_0,lat_0 = map(xav,yav,inverse=True)
      self.map    = Proj(proj='stere',ellps='WGS84',\
                         lon_0=lon_0,lat_0=lat_0,lat_ts=lat_0)

      return
   ############################################################


   ############################################################
   def get_solutions(self,METH=5):

      Psolns   = []
      cbar     = False
      for Poly in self.poly_info_objects:
         Psoln = Poly.get_solution(METH=METH)
         Psolns.append(Psoln)

      # make a MIZ_info_list object and return
      Psolns   = MIZ_info_list(Psolns,METH=METH)

      return Psolns
      #########################################################

   def reduce_area(self,vertices):

      # ===============================================================
      # get info about vertices
      vlons,vlats = np.array(vertices).transpose()
      vx,vy       = self.map(vlons,vlats,inverse=False)
      xyverts     = [(vx[i],vy[i]) for i in range(len(vx))]

      # make sure polygon is closed
      if xyverts[0] != xyverts[-1]:
         print('- closing polygon')
         xyverts.append(xyverts[0])
      # ===============================================================


      # # test plot
      # fig   = plt.figure()
      # ax    = fig.add_subplot(1,1,1)
      # vx,vy = np.array(xyverts).transpose()
      # ax.plot(vx,vy,'r',linewidth=2)


      # ===============================================================
      # test individual polygons to see if they are inside the vertices
      pil_new  = []
      for pi_old in self.poly_info_objects:
         # print(pi_old.ll_coords)
         lons,lats   = np.array(pi_old.ll_coords).transpose()
         x,y         = self.map(lons,lats,inverse=False)
         fvals       = np.array(pi_old.func_vals)
         inside      = GP.maskgrid_outside_polygon(x,y,xyverts) #true if inside poly

         # add to test plot
         # ax.plot(x,y)

         # if any points inside, check size and add to new list
         if np.any(inside):
            lons,lats   = self.map(x[inside],y[inside],inverse=True)
            fvals       = 1*fvals[inside]
            N           = len(fvals)

            # if set of points are not too small add the polygon
            if N>3:
               llc      = [(lons[i],lats[i]) for i in range(N)]
               pi_new   = poly_info(llc,func_vals=fvals)
               pil_new.append(pi_new)
         # else:
         #    print(lons.min(),lons.max(),lats.min(),lats.max())
         #    print('not inside')


      # ===============================================================

      # show test plot
      # plt.show(fig)

      # ===============================================================
      if len(pil_new)>0:# convert to poly_info_list object and return
         return poly_info_list(pil_new)
      else:
         print('Reduced set of polygons is empty\n')
         return None
      # ===============================================================


   def plot_all(self,**kwargs):

      init  = 0
      if "pobj" not in kwargs:
         init  = 1

      if "show" in kwargs:
         show           = kwargs["show"]
         kwargs["show"] = False

      for PI in self.poly_info_objects:
         if init:
            init  = 0
            pobj  = PI.plot(**kwargs)
            kwargs.update({"pobj":pobj})
         else:
            pobj  = PI.plot(**kwargs)

      if show:
         from matplotlib import pyplot as plt
         plt.show(pobj.fig)

      return pobj


#########################################################
