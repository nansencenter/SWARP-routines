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
   fname_full  = indir+"/"+fname+'.shp'
   print('\nOpening '+fname_full+'...')
   sf    = shapefile.Reader(fname_full)
   print('Processing '+fname_full+'...')
   #

   fields   = sf.fields
   Nfields  = len(fields)

   ####################################################
   # get indices of the three "forms"
   form_cats   = ['FA','FB','FC']
   Ncats       = len(form_cats)
   form_vals   = {'X':-1,'-9':np.NaN}

   for m in range(11):
      form_vals.update({str(m):m})

   MIZvals     = ['X','0','1','2','3','4'] # values of forms which will be in MIZ (floe size < 500m)
   MIZcols     = ['c','k','g','m','b','r'] # colours corresponding to each form categatory 
   MIZthresh   = 4 # in MIZ if (FA<=MIZthresh) or (FA<=MIZthresh) or (FA<=MIZthresh)

   # make a dictionary from sf.fields:
   sdct  = {}
   for n in range(1,Nfields):
      sdct.update({fields[n][0]:n-1})
   ####################################################

   ####################################################
   # find MIZ from values of FA,FB,FC
   recs     = sf.shapeRecords()
   Npolys   = len(recs)

   MIZforms = []
   print('\nTesting if polygons are in the MIZ...')
   for n in range(Npolys):

      rec   = recs[n].record
      lst   = [n]
      isMIZ = False


      # Try to define MIZ:
      fdct  = {}
      for key in form_cats:
         ind_m = sdct[key]
         sval  = rec[ind_m]     # string 
         val   = form_vals[sval]# convert from string to integer with form_vals dictionary
         fdct.update({key:val})
         
         # define MIZ if any of the forms are <= threshhold
         isMIZ = (isMIZ or val<=MIZthresh)

      if isMIZ:
         MIZforms.append([n,fdct])
   ####################################################

   Nmiz  = len(MIZforms)
   print('Number of good polygons: '+str(Nmiz))
   if Nmiz==0:
      print('No polygons in MIZ are found')
      break
   elif 0:
      # test:
      for m in range(Nmiz):
         n     = MIZforms[m][0]
         rec   = recs[n].record
         ind_n = sdct['POLY_TYPE']
         print('\nPOLY_TYPE: '+rec[ind_n])
         print('Forms:')
         print(MIZforms[m][1])
      #return

   ####################################################
   # plot outlines of polygons (colour-coded)
   for key in form_cats:
      print('\nPlotting MIZ according to '+key+'...')
      fig   = plt.figure()
      if 1:
         # manually set plot domain
         rad   = 30.          # approx radius of image (degrees)
         xmax  = rad*111.e3   # half width of image [m]
         ymax  = rad*111.e3   # half height of image [m]
         cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
         #
         lat_ts   = 60. # deg N
         lon_0    = -20.# deg E
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

      if 0:
         # only plot single polygon as a test
         # tval     = 0#n=0
         # tval     = 10#n=14
         tval     = 14#n=33
         # tval     = 28#n=72
         # tval     = 29#n=73
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
         figname  = outdir+'/'+fname+'_MIZ_'+key+'.png'
         ttl      = 'All polygons in MIZ - using '+key
   
      ####################################################
      for m in to_plot:
         n     = MIZforms[m][0]
         fdct  = MIZforms[m][1]
         val   = fdct[key]
         # print(n)
         # print(lst)
         # print(val)
         if val<=MIZthresh:
            col   = MIZcols[val+1] # colour corresponding to this value
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
      print('Saving to '+figname)
