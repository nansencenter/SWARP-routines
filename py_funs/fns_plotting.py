# ===============================================================
def pyproj_srs_to_dict(srs):
   lst      = srs.split('+')[1:]

   dic   = {}
   for ss in lst:
      key,val  = ss.strip().split('=')
      try:
         dic.update({key:float(val)})
      except:
         dic.update({key:val})

   if 'proj' in dic:
      dic.update({'projection':dic['proj']})
      del(dic['proj'])

   if 'units' in dic:
      if dic['units']!='m':
         raise ValueError('need units=m for conversion to basemap')
      del(dic['units'])

   if 'a' in dic and 'b' in dic:
      dic.update({'rsphere':[dic['a'],dic['b']]})
      del(dic['a'])
      del(dic['b'])

   elif 'a' in dic:
      dic.update({'rsphere':dic['a']})
      del(dic['a'])

   return dic
# ===============================================================


##########################################################
class plot_object:
   """
   create object with:
   pobj  = plot_object(fig=None,ax=None,cbar=None,axpos=None)
   fig  is a pyplot.figure instance
   ax   is a subplot axis of fig
   cbar is a colorbar associated with fig and a plot on ax
   - used in the plotting routines of mod_reading
   TODO move to fns_plotting
   """

   def __init__(self,fig=None,ax=None,cbar=None,axpos=None):
      from matplotlib import pyplot as plt

      if fig is None:
         self.fig   = plt.figure()
      else:
         self.fig   = fig

      if ax is None:
         self.ax = self.fig.add_subplot(1,1,1)
      else:
         self.ax  = ax

      self.cbar   = cbar
      self.axpos  = axpos

      return

   def close(self):
      plt.close(self.fig)
      return

   def get(self):
      return self.fig,self.ax,self.cbar

   def renew(self,axpos=None):
      
      pobj  = plot_object(fig=self.fig,ax=self.ax,cbar=self.cbar,axpos=axpos)

      return pobj
##########################################################


############################################################################
def start_HYCOM_map(region,cres='i'):

   if region is None or region=='Arctic':
      lonc     = -45.
      latc     = 90.
      rad      = 33 # radius in deg

   elif region=='Antarctic':
      lonc     = 180.
      latc     = -90.
      rad      = 60 # radius in deg

   elif region=='TP4':
      lonc     = -45.
      latc     = 85.
      rad      = 33 # radius in deg

   elif region=='BS1' or region=='bar':
      # Barents Sea model
      lonc     = 48.
      latc     = 74.
      rad      = 12 # radius in deg

   elif region=='FR1':
      # Fram Strait model
      lonc     = 0.5
      latc     = 78.75
      rad      = 6.0 # radius in deg

   elif region=='gre':
      # Greenland Sea - bit further south than FR1
      lonc     = -2.
      latc     = 74.
      rad      = 13.5 # radius in deg

   elif region=='svalbard':
      # Greenland Sea - bit further south than FR1
      lonc     = 2.
      latc     = 78.
      rad      = 5.4 # radius in deg

   elif region=='les':
      # Laptev/East Siberian Seas
      lonc     = 100.
      latc     = 74.
      rad      = 10.0 # radius in deg

   elif region=='beau':
      # Beaufort Sea
      lonc     = -175.
      latc     = 77.
      rad      = 23.0 # radius in deg

   elif region=='ncb':
      # North Canada & Beaufort Seas
      lonc     = -120.
      latc     = 74.
      rad      = 40.0 # radius in deg

   elif region=='lab':
      # Labrador Sea
      lonc     = -60.
      latc     = 65.
      rad      = 20.0 # radius in deg

   elif region=='can':
      # Labrador Sea
      lonc     = -80.
      latc     = 65.
      rad      = 20.0 # radius in deg

   elif region=='balt':
      # Baltic sea
      lonc     = 22.
      latc     = 64.
      rad      = 2.0 # radius in deg

   else:
      raise ValueError('Unknown region: '+region)

   return start_map_simple(lonc,latc,rad,cres=cres)
############################################################################

############################################################################
def start_map(bbox,cres='i'):

   # set plot domain from bbox
   lon0  = bbox[0]
   lat0  = bbox[1]
   lon1  = bbox[2]
   lat1  = bbox[3]
   #
   fac      = 1.25
   rad      = .5*fac*(lat1-lat0) #degrees
   #
   lon_av   = .5*(lon0+lon1)
   lat_av   = .5*(lat0+lat1)
   if 1:
      print('radius (deg lat) ='+str(rad))
      print('width/height (km) ='+str(2*rad*111))
      print('lon_av (deg E) ='+str(lon_av))
      print('lat_av (deg N) ='+str(lat_av))

   return start_map_simple(lon_av,lat_av,rad,cres=cres)
############################################################################

############################################################################
def start_map_simple(lonc,latc,rad,lat_ts=None,cres='i'):
   """
   start_map_simple(lonc,latc,rad,lat_ts=None,cres='i')
   *latc,lonc are central lat,lon
   *radius is radius of area in degrees
   *lat_ts is true scale latitude (default is latc)
   *cres is coastal resolution to use
   """
   from mpl_toolkits.basemap import Basemap

   fac      = 2
   width    = 2*rad*111.e3
   height   = 2*rad*111.e3

   if lat_ts is None:
      lat_ts   = latc

   return Basemap(width=width,height=height,\
                resolution=cres,projection='stere',\
                lat_ts=lat_ts,lat_0=latc,lon_0=lonc)
############################################################################

############################################################################
def get_nice_meridians(lonmin,lonmax):
   # get nice values of meridians for plotting gridlines
   # (also works for parallels)

   import numpy as np
   lonrng   = lonmax-lonmin

   if lonrng>80:
      dl = 20.
   elif lonrng>35:
      dl = 10.
   elif lonrng>15:
      dl = 5.
   elif lonrng>8:
      dl = 2.
   elif lonrng>4:
      dl = 1.
   elif lonrng>2:
      dl = .5
   elif lonrng>1:
      dl = .25
   else:
      dl    = lonrng/4.
      ldl   = np.log10(dl)
      dl0   = -int(np.fix(ldl))
      nd    = 1+dl0 # no of dp to round to
      dl    = np.round(dl,nd)

   l0 = (np.round(lonmin/dl)-2)*dl
   return np.arange(l0,lonmax,dl)
############################################################################

############################################################################
def finish_map(bm,do_fill=True,**kwargs):
   # finish_map(bm)
   # *bm is a basemap

   # coast/land
   bm.drawcoastlines(**kwargs)
   if do_fill:
      bm.fillcontinents(color='gray',**kwargs)

   # draw parallels and meridians.
   merids   = get_nice_meridians(bm.lonmin,bm.lonmax)
   parals   = get_nice_meridians(bm.latmin,bm.latmax)

   bm.drawparallels(parals,\
         labels=[True,False,True,True],**kwargs) # labels = [left,right,top,bottom]
   bm.drawmeridians(merids,latmax=90.,\
         labels=[True,False,False,True],**kwargs)
   bm.drawmapboundary(**kwargs)

   return
############################################################################


############################################################################
def plot_scalar(lon,lat,V,figname=None,text=None,\
      mask_lower_than=None,\
      HYCOM_region='Arctic',bmap=None,clim=None,clabel=None):

   pobj  = plot_object()
   if bmap is None:
      bmap  = start_HYCOM_map(HYCOM_region)

   if clim is None: 
      vmin  = None
      vmax  = None
   else:
      vmin,vmax   = clim

   if mask_lower_than is not None:
      import numpy as np
      good        = np.logical_not(V.mask)
      mask1       = np.zeros(V.shape,dtype='bool')
      mask1[good] = (V[good]<mask_lower_than)   # water
      mask1       = np.logical_or(V.mask,mask1)    # water or NaN
      V_          = np.ma.array(V.data,mask=mask1) # new array so we don't change V
   else:
      V_ = V # just pointer to V

   PC       = bmap.pcolor(lon,lat,V_,latlon=True,ax=pobj.ax,vmin=vmin,vmax=vmax)
   cbar     = pobj.fig.colorbar(PC)
   if clabel is not None:
      cbar.set_label(clabel,rotation=270,labelpad=20,fontsize=16)

   if text is not None:
      if HYCOM_region=='TP4':
         xyann = (0.05,.925)
      else:
         xyann = (0.4,.925)
      pobj.ax.annotate(text,xy=xyann,xycoords='axes fraction',fontsize=18)

   finish_map(bmap)

   if figname is not None:
      pobj.fig.savefig(figname)

      print('Saving '+figname+'\n')
      pobj.close()

   return pobj,bmap
############################################################################


# ===============================================================
def get_lists_for_patches(ccl,cent):
   """
   when doing complicated patches (eg to remove holes),
   need these lists (called from plot_patches)
   tst1,lst2=get_lists_for_patches(ccl,cent)
   inputs:
   ccl=[(x0,y0),(x1,y1),...]
    - list of coordinates of polygon to be drawn
   cent=(x_centroid,y_centroid)
    - centroid - used to close the polygon
   """
   from matplotlib import path
   # start the poly
   lst2  = [path.Path.MOVETO]

   # specify the poly
   lst1  = 1*ccl
   N     = len(ccl)-1
   lst2.extend(N*[path.Path.LINETO])

   # close the poly
   lst1.append(cent)
   lst2.append(path.Path.CLOSEPOLY)
   return lst1,lst2
# ===============================================================


# ===============================================================
def plot_patches(ax,MPlist,col,mapping=None,plot_holes=True):
   from matplotlib import patches,cm,collections
   from matplotlib import path
   patch_list  = []

   for poly in MPlist:
      if mapping is None:
         # plot in lon,lat
         ccl   = list(poly.exterior.coords)
      else:
         lon,lat  = poly.exterior.coords.xy
         x,y      = mapping(lon,lat)
         ccl      = [(x[i],y[i]) for i in range(len(x))]

      if len(poly.interiors)==0 or not plot_holes:
         # no holes
         # - simple way to plot
         patch_list.append(patches.Polygon(ccl,True))
      else:
         # ========================================================
         # remove holes
         x0,y0 = poly.centroid.coords[0]
         if mapping is not None:
            # lon/lat to x,y
            x0,y0 = mapping(x0,y0)
         cent  = (x0,y0)
         clist,Plist = get_lists_for_patches(ccl,cent)

         for qi in poly.interiors:
            if mapping is None:
               # plot in lon,lat
               ccl   = list(qi.coords)
               cent2 = poly.centroid.coords[0]
            else:
               lon,lat  = qi.coords.xy
               x,y      = mapping(lon,lat)
               ccl      = [(x[i],y[i]) for i in range(len(x))]
               x0,y0    = qi.centroid.coords[0]
               x0,y0    = mapping(x0,y0)
               cent2    = (x0,y0)

            clist2,Plist2 = get_lists_for_patches(ccl,cent2)
            clist.extend(clist2)
            Plist.extend(Plist2)

         patch_list.append(patches.PathPatch(path.Path(clist,Plist)))
         # ========================================================


   pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=1.) # alpha is transparency
   pc.set_facecolor(col)
   ax.add_collection(pc)
   return
############################################################################
