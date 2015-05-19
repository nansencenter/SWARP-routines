import os,sys
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
from matplotlib import lines as mlines

# pyshp->shapefile info:
# https://pypi.python.org/pypi/pyshp
import shapefile

# shapely info:
# official docs: http://toblerity.org/shapely/shapely.geometry.html
# more demos: http://toblerity.org/shapely/manual.html
import shapely.geometry as shgeom
import shapely.ops      as shops

swarp = os.getenv('SWARP_ROUTINES')
pyswarp  = swarp+'/py_funs'
if pyswarp not in sys.path:
   sys.path.append(pyswarp)

import shapefile_utils as SFU

############################################################################
def start_map(bbox,cres='i'):
   # set plot domain from bbox
   lon0  = bbox[0]
   lat0  = bbox[1]
   lon1  = bbox[2]
   lat1  = bbox[3]
   #
   fac      = 1.5
   rad      = .5*fac*(lat1-lat0) #degrees
   width    = 2*rad*111.e3
   height   = 2*rad*111.e3
   #
   lon_av   = .5*(lon0+lon1)
   lat_av   = .5*(lat0+lat1)
   if 0:
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
def finish_map(bm):
   # finish_map(bm)
   # *bm is a basemap

   def get_nice_meridians(lonmin,lonmax):
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

   # coast/land
   bm.drawcoastlines()
   bm.fillcontinents(color='gray')

   # draw parallels and meridians.
   merids   = get_nice_meridians(bm.lonmin,bm.lonmax)
   parals   = get_nice_meridians(bm.latmin,bm.latmax)

   bm.drawparallels(parals,\
         labels=[True,False,True,True]) # labels = [left,right,top,bottom]
   bm.drawmeridians(merids,latmax=90.,\
         labels=[True,False,False,True])
   bm.drawmapboundary() # fill_color='aqua')

   return
############################################################################

############################################################################
def xy2coords(x,y):
   # convert x,y list to list of coordinates (tuples)

   pts   = []
   for n in range(len(x)):
      pts.append((x[n],y[n]))

   return pts
############################################################################

############################################################################
def coords2xy(pts):
   # convert list of coordinates (tuples) to x,y list

   x  = []
   y  = []
   for n in range(len(pts)):
      x.append(pts[n][0])
      y.append(pts[n][1])

   return np.array(x),np.array(y)
############################################################################

############################################################################
def plot_coords(bm,pts,**kwargs):

   x,y   = coords2xy(pts)
   bm.plot(x,y,**kwargs)

   return
############################################################################

############################################################################
def make_valid_poly(shp):

   pts      = shp.points
   Nparts   = len(shp.parts)
   if Nparts==1:
      poly  = shgeom.Polygon(pts)
   else:
      polys    = []

      # all parts but last
      for n in range(Nparts-1):
         i0 = shp.parts[n]
         i1 = shp.parts[n+1]
         pn = shgeom.Polygon(pts[i0:i1])
         polys.append(pn)

      # last part
      for n in [Nparts-1]:
         i0 = shp.parts[n]
         i1 = len(pts)
         pn = shgeom.Polygon(pts[i0:i1])
         polys.append(pn)

      if 1:
         # just take 1st poly
         poly  = polys[0]
      else:
         # include all parts
         # TODO include these?
         # TODO does unary_union remove holes
         poly  = shgeom.MultiPolygon(polys)

   return poly
############################################################################

indir    = 'test_inputs'
# indir    = 'test_inputs2'
fnames   = os.listdir(indir)
snames   = []

outdir   = 'test_outputs'
if not os.path.exists(outdir):
   os.mkdir(outdir)

for fname in fnames:
   ext   = os.path.splitext(fname)
   if ext[1]=='.shp':
      snames.append(ext[0])

# MIZ definition stuff
form_cats   = ['FA','FB','FC']
form_vals   = SFU.DMI_form_dictionary()
MIZthresh   = 4 # in MIZ if (FA<=MIZthresh) or (FB<=MIZthresh) or (FC<=MIZthresh)

# string values (for plotting later)
MIZvals  = ['X'] # values of forms which will be in MIZ (floe size < 500m)
for n in range(MIZthresh+1):
   MIZvals.append(str(n))
MIZcols     = ['c','k','g','m','b','r'] # colours corresponding to each form categatory 

for fname in snames:
   fname_full  = indir+"/"+fname+'.shp'
   print('\nOpening '+fname_full+'...')
   sf    = shapefile.Reader(fname_full)
   print('Processing '+fname_full+'...')
   #
   sf_info  = SFU.extract_shapefile_info(sf,get_holes=True)

   ####################################################
   # find MIZ from values of FA,FB,FC
   Npolys   = len(sf_info)
   MIZforms = []

   print('\nTesting if polygons are in the MIZ...')
   for n in range(Npolys):

      shp_info = sf_info[n]
      lut      = shp_info[1]
      lst      = [n]
      isMIZ    = False

      # Try to define MIZ:
      fdct  = {}
      for key in form_cats:
         sval  = lut[key]        # string 
         val   = form_vals[sval] # convert from string to integer with form_vals dictionary
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
   # for key in form_cats:
   for key in ['FC']:

      print('\nPlotting MIZ according to '+key+'...')
      fig   = plt.figure()

      PLOT_ALL = 1
      if PLOT_ALL==0:
         # only plot single polygon as a test
         PLOT_COMBINED  = 0

         # tval     = 0#n=0
         # tval     = 10#n=14
         # tval     = 14#n=33
         tval     = 28#n=72 #has holes (islands)
         # tval     = 29#n=73
         to_plot  = [tval]
         n        = MIZforms[tval][0]
         dct      = MIZforms[tval][1]

         cdate    = fname[:8]
         if not os.path.exists(outdir+'/'+cdate):
            os.mkdir(outdir+'/'+cdate)
         figname  = outdir+'/'+cdate+'/'+fname+'_MIZpoly'+str(n)+'.png'

         ttl   = 'Polygon '+str(n)+': '
         for key in dct.keys():
            ttl   = ttl+key+'='
            val   = dct[key]
            if val>=0:
               ttl   = ttl+str(val)+', '
            else:
               ttl   = ttl+"'X', "
         ttl   = ttl[:-2]

         # bounding box of the single shape
         bbox  = sf.shape(n).bbox
         bm    = start_map(bbox,cres='f')
         print('bbox=')
         print(bbox)

      else:
         # plot all polygons in MIZ
         to_plot  = range(len(MIZforms))
         cdate    = fname[:8]
         if not os.path.exists(outdir+'/'+cdate):
            os.mkdir(outdir+'/'+cdate)

         PLOT_COMBINED  = 1 # merge neighbouring MIZ polygons together
         if PLOT_COMBINED==0:
            figname  = outdir+'/'+cdate+'/'+fname+'_MIZ_'+key+'.png'
         else:
            figname  = outdir+'/'+cdate+'/'+fname+'_MIZ_'+key+'2.png'
         ttl      = 'All polygons in MIZ - using '+key

         # set plot domain from overall bbox
         bbox  = sf.bbox # [lonmin,latmin,lonmax,latmax]
         bm    = start_map(bbox)
      ####################################################
   
      ####################################################
      Ncols = len(MIZcols)
      x_all = []
      y_all = []

      MPlist   = []
      for m in to_plot:
         n     = MIZforms[m][0]
         fdct  = MIZforms[m][1]
         val   = fdct[key]

         #########################################################
         if val<=MIZthresh:
            col   = MIZcols[val+1] # colour corresponding to this value
            poly  = sf_info[n][0]
            #
            if PLOT_COMBINED==0:
               SFU.plot_poly(poly,pobj=bm,plot_holes=True,color=col,latlon=True)
               #plot_coords(bm,pts,color=col,latlon=True)
            else:
               # merge polygons
               MPlist.append(poly)
               # print('appending poly')
               # print(len(MPlist))
               # if 0:
               #    # test unary union
               #    uu = shops.unary_union(MPlist)
               #    print('unary union ok')
         #########################################################

      #########################################################
      if PLOT_COMBINED==1:
         MP = shops.unary_union(MPlist) # merge neighbouring polygons
         if not hasattr(MP,'geoms'):
            MP = shgeom.MultiPolygon(MP)
         print('len(MP)='+str(len(MP.geoms)))

         print("Doing plot now...")
         Ngrps = len(MP.geoms)
         
         cols  = ['c','b','g','m','r'] # colours corresponding to each group of MIZ ice
         for m in range(Ngrps):
            Ninner   = len(MP.geoms[m].interiors)
            print('No of interior regions')
            print([m,Ninner])
            #
            m_    = np.mod(m,len(cols))
            SFU.plot_poly(MP.geoms[m],pobj=bm,color=cols[m_],plot_holes=True,latlon=True)
         ####################################################
      else:
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
      #########################################################
