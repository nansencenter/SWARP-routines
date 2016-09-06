from matplotlib import pyplot as plt

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

   def get(self):
      return self.fig,self.ax,self.cbar

   def renew(self,axpos=None):
      
      pobj  = plot_object(fig=self.fig,ax=self.ax,cbar=self.cbar,axpos=axpos)

      return pobj
##########################################################


############################################################################
def start_HYCOM_map(region,cres='i'):

   from mpl_toolkits.basemap import Basemap

   if region=='Arctic':
      lonc     = -45.
      latc     = 90.
      lat_ts   = latc
      rad      = 33 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='Antarctic':
      lonc     = 180.
      latc     = -90.
      lat_ts   = latc
      rad      = 60 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='TP4':
      lonc     = -45.
      latc     = 85.
      lat_ts   = latc
      rad      = 33 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='BS1' or region=='bar':
      # Barents Sea model
      lonc     = 48.
      latc     = 74.
      lat_ts   = latc
      rad      = 12 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='FR1':
      # Fram Strait model
      lonc     = 0.5
      latc     = 78.75
      lat_ts   = latc
      rad      = 6.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='gre':
      # Greenland Sea - bit further south than FR1
      lonc     = -2.
      latc     = 74.
      lat_ts   = latc
      rad      = 13.5 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='les':
      # Laptev/East Siberian Seas
      lonc     = 100.
      latc     = 74.
      lat_ts   = latc
      rad      = 10.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='beau':
      # Beaufort Sea
      lonc     = -175.
      latc     = 77.
      lat_ts   = latc
      rad      = 23.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='ncb':
      # North Canada & Beaufort Seas
      lonc     = -120.
      latc     = 74.
      lat_ts   = latc
      rad      = 40.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='lab':
      # Labrador Sea
      lonc     = -60.
      latc     = 65.
      lat_ts   = latc
      rad      = 20.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='can':
      # Labrador Sea
      lonc     = -80.
      latc     = 65.
      lat_ts   = latc
      rad      = 20.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   elif region=='balt':
      # Baltic sea
      lonc     = 22.
      latc     = 64.
      lat_ts   = latc
      rad      = 2.0 # radius in deg
      width    = 2*rad*111.e3
      height   = 2*rad*111.e3
      #
      bm = Basemap(width=width,height=height,\
                   resolution='i',projection='stere',\
                   lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

   else:
      raise ValueError('Unknown region: '+region)

   return bm
############################################################################

############################################################################
def start_map(bbox,cres='i'):

   from mpl_toolkits.basemap import Basemap

   # set plot domain from bbox
   lon0  = bbox[0]
   lat0  = bbox[1]
   lon1  = bbox[2]
   lat1  = bbox[3]
   #
   fac      = 2
   rad      = .5*fac*(lat1-lat0) #degrees
   width    = 2*rad*111.e3
   height   = 2*rad*111.e3
   #
   lon_av   = .5*(lon0+lon1)
   lat_av   = .5*(lat0+lat1)
   if 1:
      print('radius (deg lat) ='+str(rad))
      print('width (km) ='+str(width/1.e3))
      print('height (km) ='+str(height/1.e3))
      print('lon_av (deg E) ='+str(lon_av))
      print('lat_av (deg N) ='+str(lat_av))
   #
   bm = Basemap(width=width,height=height,\
                resolution=cres,projection='stere',\
                lat_ts=lat_av,lat_0=lat_av,lon_0=lon_av)
   return bm
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
def finish_map(bm,**kwargs):
   # finish_map(bm)
   # *bm is a basemap

   # coast/land
   bm.drawcoastlines(**kwargs)
   bm.fillcontinents(color='gray',**kwargs)

   # draw parallels and meridians.
   merids   = get_nice_meridians(bm.lonmin,bm.lonmax)
   parals   = get_nice_meridians(bm.latmin,bm.latmax)

   bm.drawparallels(parals,\
         labels=[True,False,True,True],**kwargs) # labels = [left,right,top,bottom]
   bm.drawmeridians(merids,latmax=90.,\
         labels=[True,False,False,True],**kwargs)
   bm.drawmapboundary(**kwargs) # fill_color='aqua')

   return
############################################################################


############################################################################
def plot_anomaly(lon,lat,anom,anom_fig,text=None,\
      HYCOM_region='Arctic',clim=None,clabel=None):

   pobj     = plot_object()
   bmap     = start_HYCOM_map(HYCOM_region)

   if clim is None: 
      vmin  = None
      vmax  = None
   else:
      vmin,vmax   = clim

   PC       = bmap.pcolor(lon,lat,anom,latlon=True,ax=pobj.ax,vmin=vmin,vmax=vmax)
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
   pobj.fig.savefig(anom_fig)

   print('Saving '+anom_fig+'\n')
   plt.close(pobj.fig)
   return
############################################################################


############################################################################
def plot_patches(ax,MPlist,col):
   from matplotlib import patches,cm,collections
   patch_list  = []

   for poly in MPlist:
      ccl   = list(poly.exterior.coords)
      patch_list.append(patches.Polygon(ccl,True))

   pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=.5)
   pc.set_facecolor(col)
   ax.add_collection(pc)
   return
############################################################################
