############################################################################
def start_HYCOM_map(region,cres='i'):

   from mpl_toolkits.basemap import Basemap

   if region=='TP4':
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

   elif region=='BS1':
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
   bm.drawcoastlines()
   bm.fillcontinents(color='gray')

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
