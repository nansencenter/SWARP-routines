import os,sys
import numpy as np
import shapefile
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
from matplotlib import lines as mlines

swarp = os.getenv('SWARP_ROUTINES')
pyswarp  = swarp+'/py_funs'
if pyswarp not in sys.path:
   sys.path.append(pyswarp)

import ordered_set as ordset

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

# indir    = 'test_inputs'
indir    = 'test_inputs2'
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
   # for key in ['FC']:
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
         PLOT_COMBINED  = 0

         # tval     = 0#n=0
         # tval     = 10#n=14
         tval     = 14#n=33
         # tval     = 28#n=72
         # tval     = 29#n=73
         to_plot  = [tval]
         n        = MIZforms[tval][0]
         lst      = MIZforms[tval][1]

         cdate    = fname[:8]
         if not os.path.exists(outdir+cdate):
            os.mkdir(outdir+cdate)
         figname  = outdir+'/'+cdate+'/'+fname+'_MIZpoly'+str(n)+'.png'
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
         cdate    = fname[:8]
         if not os.path.exists(outdir+cdate):
            os.mkdir(outdir+cdate)

         PLOT_COMBINED  = 0
         if PLOT_COMBINED==0:
            figname  = outdir+'/'+cdate+'/'+fname+'_MIZ_'+key+'.png'
         else:
            figname  = outdir+'/'+cdate+'/'+fname+'_MIZ_'+key+'2.png'
         ttl      = 'All polygons in MIZ - using '+key
   
      ####################################################
      Ncols = len(MIZcols)
      x_all = []
      y_all = []

      Ndone  = 0 # need to 
      for m in to_plot:
         n     = MIZforms[m][0]
         fdct  = MIZforms[m][1]
         val   = fdct[key]
         # print(n)
         # print(lst)
         # print(val)

         #########################################################
         Mpolys   = []
         if val<=MIZthresh:
            col   = MIZcols[val+1] # colour corresponding to this value
            #
            shp   = sf.shape(n)
            pts   = shp.points
            #
            x  = []
            y  = []
            for p in range(len(pts)):
               p0 = pts[p][0]
               p1 = pts[p][1]
               x.append(p0)
               y.append(p1)

               if PLOT_COMBINED==1:
                  # test if point in (x_all,y_all) already
                  i0 = -1
                  if (p0 in x_all):
                     i0 = x_all.index(p0)
                  i1 = -1
                  if (p1 in y_all):
                     i1 = y_all.index(p1)

                  # if x in x_all and y in y_all already,
                  # check if they have the same index
                  # if not add to list
                  if not((i0==i1) and (i0>=0)):
                     x_all.append(p0)
                     y_all.append(p1)

         # #  TODO use shapely to convert (x,y) to geometric object 
         # # usable in operations (eg disjoint,
         # # unary union = merge overlapping polys into single poly)
         # import shapely.geometry as shgeom
         # import shapely.ops      as shops

         # poly  = shgeom.Polygon(pts)
         # if Ndone==0:
         #    Mp = shgeom.MultiPolygon([poly])
         # else:
         #    Mp = shops.unary_union([Mp,poly]) # this should automatically merge neighbouring polygons

         # Ndone = Ndone+1
         #########################################################

   
         #########################################################
         if PLOT_COMBINED==0:
            bm.plot(np.array(x),np.array(y),col,latlon=True)
         #########################################################

      #########################################################
      if PLOT_COMBINED==1:
         print("Doing plot now...")
         bm.plot(np.array(x_all),np.array(y_all),'.k',latlon=True)
      ####################################################

      ####################################################
      # Draw legend
      # MIZvals     = ['X','0','1','2','3','4'] # values of forms which will be in MIZ (floe size < 500m)
      # MIZcols     = ['c','k','g','m','b','r'] # colours corresponding to each form categatory 
      Ncols    = len(MIZcols)
      handles  = Ncols*[0]
      for m in range(Ncols):
         handles[m]  = mlines.Line2D([], [], color=MIZcols[m],label=MIZvals[m]) # also eg marker='*',markersize='15'
      plt.legend(handles=handles)
      ####################################################
   
      finish_map(bm)
      plt.title(ttl,y='1.06')
      plt.savefig(figname)
      plt.close()
      fig.clf()
      print('Saving to '+figname)
