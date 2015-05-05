import os,sys
import numpy as np
import shapefile
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt

############################################################################
def finish_map(bm):
   # finish_map(bm)
   # *bm is a basemap
   bm.drawcoastlines()
   bm.fillcontinents(color='gray')

   # draw parallels and meridians.
   bm.drawparallels(np.arange(60.,91.,10.),\
         labels=[True,False,True,True]) # labels = [left,right,top,bottom]
   bm.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
         labels=[True,False,False,True])
   bm.drawmapboundary() # fill_color='aqua')

   return
############################################################################

indir    = 'test_inputs'
fnames   = os.listdir(indir)
snames   = []

outdir   = 'test_outputs'
if not os.path.exists(outdir):
   os.mkdir(outdir)

for fname in fnames:
   ext   = os.path.splitext(fname)
   if ext[1]=='.shp':
      snames.append(ext[0])

for fname in snames:
   sf    = shapefile.Reader(indir+"/"+fname+'.shp')
   fig   = plt.figure()
   #
   if 1:
      # manually set plot domain
      rad   = 30.          # approx radius of image (degrees)
      xmax  = rad*111.e3   # half width of image [m]
      ymax  = rad*111.e3   # half height of image [m]
      cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
      #
      lat_ts   = 60. # deg N
      lon_0    = 00. # deg E
      lat_0    = 60. # deg N
      #
      bm = Basemap(width=2*xmax,height=2*ymax,\
                   resolution=cres,projection='stere',\
                   lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)
   else:
      # set plot domain from bbox
      bbox  = sf.bbox # [lonmin,latmin,lonmax,latmax]
      lon0  = bbox[0]
      lat0  = bbox[1]
      lon1  = bbox[2]
      lat1  = bbox[3]
      #
      lon_av   = .5*(lon0+lon1)
      lat_av   = .5*(lat0+lat1)
      bm       = Basemap(llcrnrlon=lon0,llcrnrlat=lat0,\
                         urcrnrlon=lon1,urcrnrlat=lat1,\
                         resolution='i',projection='stere',\
                         lat_ts=lat_av,lat_0=lat_av,lon_0=lon_av)

   fields   = sf.fields
   Nfields  = len(fields)

   ####################################################
   # get indices of the three "forms"
   dct      = {'FA':-1,'FB':-1,'FC':-1}
   Ncats    = len(dct.keys())
   MIZvals  = ['X','1','2','3','4'] # values of forms which will be in MIZ (floe size < 500m)
   MIZcols  = ['c','g','m','b','r'] # colours corresponding to each form categatory 

   for n in range(1,Nfields):
      for m in range(Ncats):
         key   = dct.keys()[m]
         if fields[n][0]==key:
            dct[key] = n-1
   ####################################################

   ####################################################
   # find MIZ from values of FA,FB,FC
   recs     = sf.shapeRecords()
   Npolys   = len(recs)

   MIZforms = []
   for n in range(Npolys):

      rec   = recs[n].record
      lst   = Ncats*[-1]
      isMIZ = 0

      for m in range(Ncats):
         key      = dct.keys()[m]
         ind_m    = dct[key]
         val      = rec[ind_m]
         if val in MIZvals[1:]:
            lst[m]   = int(val)
            isMIZ    = 1
         elif val==MIZvals[0]:
            lst[m]   = 0
            isMIZ    = 1

      if isMIZ==1:
         MIZforms.append([n,np.array(lst)])
   ####################################################

   ####################################################
   # plot outlines of polygons (colour-coded)

   if 0:
      # only plot single polygon as a test
      # tval     = 0#n=0
      # tval     = 10#n=14
      tval     = 29#n=73
      to_plot  = [tval]
      n        = MIZforms[tval][0]
      lst      = MIZforms[tval][1]

      figname  = outdir+'/'+fname+'_MIZpoly'+str(n)+'.png'
      ttl      = 'Polygon '+str(n)+': FA='

      if lst[0]>0:
         ttl   = ttl+str(lst[0])+', FB='
      else:
         ttl   = ttl+"'X', FB="

      if lst[1]>0:
         ttl   = ttl+str(lst[1])+', FC='
      else:
         ttl   = ttl+"'X', FC="

      if lst[2]>0:
         ttl   = ttl+str(lst[2])
      else:
         ttl   = ttl+"'X'"

   else:
      # plot all polygons in MIZ
      to_plot  = range(len(MIZforms))
      figname  = outdir+'/'+fname+'MIZ.png'
      ttl      = 'All polygons in MIZ'

   for m in to_plot:
      n     = MIZforms[m][0]
      lst   = MIZforms[m][1]
      #
      val   = lst.max()    # take val as max value over all 3 thickness categories
      print(n)
      print(lst)
      print(val)
      col   = MIZcols[val] # colour corresponding to this value
      #
      shp   = sf.shape(n)
      pts   = shp.points
      #
      x  = []
      y  = []
      for p in range(len(pts)):
         x.append(pts[p][0])
         y.append(pts[p][1])

      bm.plot(np.array(x),np.array(y),col,latlon=True)
   ####################################################

   finish_map(bm)
   plt.title(ttl,y='1.06')
   plt.savefig(figname)
   plt.close()
   fig.clf()
