import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

sys.path.append('../py_funs')
import fns_plotting as Fplt

def read_txt_file(fname):
   # reads in text file from Stefan,
   # outputs list of polygons

   fid   = open(fname,'r')
   lin   = fid.readline()    # 1 line
   fid.close()
   #
   lins     = lin.strip('\n').split('Line')[1:] # line 1 = ice edge, line 2 = MIZ-pack line
   Lins     = []

   # loop over line types
   for ll in lins:
      ll2   = ll.split('Part ')[1:]
      parts = []

      ################################################################################
      # loop over parts
      for part in ll2:
         ll3         = part.split('n')[1:-1] # now a list with members ['lon1 lat1','lon2 lat2',...]

         # loop over coords
         ll_coords   = [] # this will be a list of tuples [(lon1,lat1),(lon2,lat2),...]
         for cc in ll3:
            cc2      = cc.split()
            lonlat   = (float(cc2[0]),float(cc2[1]))
            ll_coords.append(lonlat)

         parts.append(ll_coords)# list of list of tuples
      ################################################################################

      Lins.append(parts)

   return Lins

fdir  = '/Volumes/sim/tim/Projects/SWARP/Stefan/Feb2014/SWARP_MIZ_maps_def/txt'
fname = fdir+'/20140201_lonlat.txt'
Lins  = read_txt_file(fname)

if 1:
   # TEST PLOT
   # setup stereographic basemap
   # - for plotting and also for limiting search area
   # lat_ts is latitude of true scale.
   # (lon_0,lat_0) is central point -> (x,y)=(0,0)
   rad   = 5.          # approx radius of image (degrees)
   xmax  = rad*111.e3   # half width of image [m]
   ymax  = rad*111.e3   # half height of image [m]
   cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
   #
   lon_0    = 35.    # deg E
   lat_0    = 80.    # deg N
   lat_ts   = lat_0  # deg N
   #
   fig   = plt.figure()
   ax1   = fig.add_subplot(1,1,1)
   bmap  = Basemap(width=2*xmax,height=2*ymax,\
                   resolution=cres,projection='stere',\
                   lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)

   lstil = ['-','--']
   for n,Lin in enumerate(Lins):
      for part in Lin:
         lon,lat  = np.array(part).transpose()
         x,y      = bmap(lon,lat)
         bmap.plot(x,y,linestyle=lstil[n],ax=ax1)

   Fplt.finish_map(bmap,ax=ax1)
   plt.show(fig)
   ax1.cla()
   plt.close(fig)
