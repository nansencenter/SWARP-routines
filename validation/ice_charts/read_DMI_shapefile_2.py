import os,sys
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
from matplotlib import lines  as mlines

# plot patches
from matplotlib import patches,cm,collections

# pyshp->shapefile info:
# https://pypi.python.org/pypi/pyshp
import shapefile
import rtree

# shapely info:
# official docs: http://toblerity.org/shapely/shapely.geometry.html
# more demos: http://toblerity.org/shapely/manual.html
import shapely.geometry as shgeom
import shapely.ops      as shops

swarp = os.getenv('SWARP_ROUTINES')
pyswarp  = swarp+'/py_funs'
if pyswarp not in sys.path:
   sys.path.append(pyswarp)

import shapefile_utils  as SFU
import fns_plotting     as Mplt


def plot_patches(ax,MP,col):
   patch_list  = []

   if not hasattr(MP,'geoms'):
      ccl   = list(MP.exterior.coords)
      patch_list.append(patches.Polygon(ccl,True))
   else:
      for gg in MP.geoms:
         ccl   = list(gg.exterior.coords)
         patch_list.append(patches.Polygon(ccl,True))

   pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=.5)
   pc.set_facecolor(col)
   ax.add_collection(pc)
   return

def test_plot(MP_WTR,MP_PACK,MP_MIZ,MIZout,MPlist_MIZ=[],figname=None):
   fig      = plt.figure()
   ax       = fig.add_subplot(111)

   # plot water and pack as patches
   x0,y0,x1,y1 = MP_WTR .bounds
   X0,Y0,X1,Y1 = MP_PACK.bounds
   x0 = min(x0,X0)
   y0 = min(y0,Y0)
   x1 = max(x1,X1)
   y1 = max(y1,Y1)

   plot_patches(ax,MP_WTR,'c')
   plot_patches(ax,MP_PACK,'r')
   plot_patches(ax,MP_MIZ,'y')

   ax.set_xlim([x0,x1])
   ax.set_ylim([y0,y1])

   lstils   = ['.b','.r']
   for MIZi in MIZout:
      Lon   = MIZi['lon']
      Lat   = MIZi['lat']
      ax.plot(Lon,Lat,':')
      for flag in range(2):
         lon   = Lon[MIZi['flags']==flag]
         lat   = Lat[MIZi['flags']==flag]
         ax.plot(lon,lat,lstils[flag])

   for poly in MPlist_MIZ:
      lon,lat  = poly.exterior.coords.xy
      ax.plot(lon,lat,':k')

   if figname is not None:
      print('Saving '+figname+'...\n')
      fig.savefig(figname,bbox_inches='tight')
   else:
      plt.show(fig)

   ax.cla()
   plt.close(fig)
   return

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

# ================================================================================
# MIZ definition stuff
form_cats_MIZ   = ['FA','FB','FC']
form_vals   = SFU.DMI_form_dictionary()
MIZthresh   = 4 # in MIZ if (FA<=MIZthresh) or (FB<=MIZthresh) or (FC<=MIZthresh)
# ================================================================================


# string values (for plotting later)
MIZvals  = ['X'] # values of forms which will be in MIZ (floe size < 500m)
for n in range(MIZthresh+1):
   MIZvals.append(str(n))
MIZcols     = ['c','k','g','m','b','r'] # colours corresponding to each form categatory 

RECORDS_PRINTED   = 0

for fname in snames:
   fname_full  = indir+"/"+fname+'.shp'
   print('\nOpening '+fname_full+'...')
   sf    = shapefile.Reader(fname_full)
   print('Processing '+fname_full+'...')

   # ============================================================
   # get shapes and records from shapefile
   # as list of lists:
   # sf_info[n]   = [poly,record]
   # where poly is a shapely polygon,
   # and record is all the metadata
   sf_info  = SFU.extract_shapefile_info(sf,get_holes=True)
   # ============================================================


   # ============================================================
   # separate water, MIZ and pack:
   # find MIZ from values of FA,FB,FC
   Npolys      = len(sf_info)
   MIZ_forms   = [] # MIZ polygons
   WTR_forms   = [] # water polygons
   PACK_forms  = [] # pack ice polygons

   print('\nTesting if polygons are in the MIZ...')
   for n in range(Npolys):

      shp_info = sf_info[n]
      lut      = shp_info[1]
      if not RECORDS_PRINTED:
         print('\nExample record:')
         print(lut)
         print('\n')
         RECORDS_PRINTED   = 1

      lst      = [n]
      isMIZ    = False

      # Try to define MIZ:
      fdct  = {}
      for key in form_cats_MIZ:
         sval  = lut[key]        # string 
         val   = form_vals[sval] # convert from string to integer with form_vals dictionary
         fdct.update({key:val})
         
         # define MIZ if any of the forms are <= threshhold
         isMIZ = (isMIZ or val<=MIZthresh)

      if isMIZ:
         MIZ_forms.append([n,fdct])
      elif lut['POLY_TYPE']=='I':
         # not MIZ but ice
         PACK_forms.append(n)
      else:
         # not ice
         WTR_forms.append(n)
   # ============================================================


   # ============================================================
   # list of shapely polygons - convert to shapely multi-polygons and take their union
   MPlist_MIZ     = []
   MPlist_PACK    = []
   MPlist_WTR     = []
   MIZform_vals   = []

   # merge pack polygons (also dissolving holes)
   Npack = len(PACK_forms)
   print('Merging '+str(Npack)+' pack polygons...\n')
   for n in PACK_forms:
      poly,lut    = sf_info[n]
      MPlist_PACK.append(poly)
   MP_PACK  = shops.unary_union(MPlist_PACK) # merge neighbouring polygons (dissolve holes)

   # merge wtr polygons (also dissolving holes)
   Nwtr  = len(WTR_forms)
   print('Merging '+str(Nwtr)+' water polygons...\n')
   for n in WTR_forms:
      poly,lut    = sf_info[n]
      MPlist_WTR.append(poly)
   MP_WTR   = shops.unary_union(MPlist_WTR) # merge neighbouring polygons (dissolve holes)


   # merge MIZ polygons (also dissolving holes)
   N_MIZ = len(MIZ_forms)
   print('Merging '+str(N_MIZ)+' MIZ polygons...\n')
   for n,fdct in MIZ_forms:
      poly,lut    = sf_info[n]
      MPlist_MIZ.append(poly)
   MP_MIZ   = shops.unary_union(MPlist_MIZ) # merge neighbouring polygons (dissolve holes)

   if not hasattr(MP_MIZ,'geoms'):
      # 1 polygon
      MP_list_MIZ_merged   = [MP_MIZ]
   else:
      MP_list_MIZ_merged   = []
      for poly in MP_MIZ.geoms:
         MP_list_MIZ_merged.append(poly)

   # ==================================================================
   # classify MIZ boundary
   # TODO match to original shapes to extract correct areas
   # TODO are holes important?
   N_MIZ_merged   = len(MP_list_MIZ_merged)
   print('Sorting '+str(N_MIZ_merged)+' MIZ polygons...\n')

   MIZout   = []
   count    = 0
   for poly in MP_list_MIZ_merged:
      count += 1
      print("Progress: ",count,N_MIZ_merged)
      llc            = list(poly.exterior.coords)
      lonc,latc      = poly.exterior.coords.xy
      Nc             = len(llc)
      boundary_flags = Nc*[2] # 2 = unknown, 0 = ice edge, 1 = MIZ-pack edge

      # =================================================
      # ice edge
      buf   = 1.e-5
      if poly.exterior.intersects(MP_WTR.buffer(buf)):
         IE = poly.exterior.intersection(MP_WTR.buffer(buf)) # coords in ice edge
         rest  = poly.exterior.symmetric_difference(IE)
         Poly  = rest.union(IE)

         if not hasattr(IE,'geoms'):
            LLc   = list(IE.coords)
         else:
            LLc   = []
            for IE_ in IE.geoms:
               LLc.extend(list(IE_.coords))

         # reduce list of points to those in the original poly
         # - intersection gets filled somehow
         llc   = list(poly.exterior.coords)
         for cc in LLc:
            if cc in llc:
               idx   = llc.index(cc)
               boundary_flags[idx]  = 0
      # =================================================


      # =================================================
      # MIZ-pack edge
      if poly.exterior.intersects(MP_PACK.buffer(buf)):
         IE = poly.exterior.intersection(MP_PACK.buffer(buf)) # coords in MIZ-pack edge

         if not hasattr(IE,'geoms'):
            LLc   = list(IE.coords)
         else:
            LLc   = []
            for IE_ in IE.geoms:
               LLc.extend(list(IE_.coords))

         # reduce list of points to those in the original poly
         # - intersection gets filled somehow
         llc   = list(poly.exterior.coords)
         for cc in LLc:
            if cc in llc:
               idx   = llc.index(cc)
               boundary_flags[idx]  = 1
      # =================================================


      for flag in range(3):
         print('#'+str(flag)+': '+str(boundary_flags.count(flag)))
      print('Total number of points: '+str(Nc)+"\n")

      MIZout.append({'lon'    :np.array(lonc),
                     'lat'    :np.array(latc),
                     'flags'  :np.array(boundary_flags),
                     'record' :lut})
      # ============================================================
   # ==================================================================


   cdate = fname[:8]
   Fname = outdir+'/'+cdate+'/'+fname
   if not os.path.exists(outdir):
      os.mkdir(outdir)
   if not os.path.exists(outdir+'/'+cdate):
      os.mkdir(outdir+'/'+cdate)
   test_plot(MP_WTR,MP_PACK,MP_MIZ,MIZout,MPlist_MIZ=MPlist_MIZ,figname=Fname+'.png')


   # ============================================================
   # write text files
   tfil1 = Fname+'_MIZpolys.txt'
   # tfil2 = Fname+'_MIZpolys_info.txt'
   print('\nSaving polygons + boundary flags to '+tfil1)
   # print('Saving metadata to '+tfil2)
   print('\n')
   tf1   = open(tfil1,'w')
   # tf2   = open(tfil2,'w')
   blk   = 3*' '

   for npoly,dct in enumerate(MIZout):
      npts  = len(dct['lon'])

      # ==================================================
      # headers
      if npoly==0:
         ss = 'Polygon'+blk+\
              'lon'+blk+\
              'lat'+blk+\
              'flag\n'
         tf1.write(ss)

         # keys  = dct['record'].keys()
         # ss    = 'Polygon'+blk
         # for key in keys[:-1]:
         #    ss   +=key+blk
         # ss   += keys[-1]+'\n'
         # tf2.write(ss)
      # ==================================================

      # # ==================================================
      # # metadata for poly
      # ss = str(npoly)+blk
      # for key in keys[:-1]:
      #    ss   += str(dct['record'][key])+blk
      # ss   += str(dct['record'][keys[-1]])+'\n'
      # tf2.write(ss)
      # # ==================================================


      # ==================================================
      # all points + classification
      for npt in range(npts):
         ss0   = [npoly,\
                  dct['lon'][npt],
                  dct['lat'][npt],
                  dct['flags'][npt]]

         ss = ''
         for val in ss0[:-1]:
            ss   += str(val)+blk
         ss   += str(ss0[-1])+"\n"
         tf1.write(ss)
      # ==================================================

   tf1.close()
   # tf2.close()
   # ============================================================

   sys.exit()


   ####################################################

   Nmiz  = len(MIZ_forms)
   print('Number of good polygons: '+str(Nmiz))
   if Nmiz==0:
      print('No polygons in MIZ are found')
      break
   elif 0:
      # test:
      for m in range(Nmiz):
         n     = MIZ_forms[m][0]
         rec   = recs[n].record
         ind_n = sdct['POLY_TYPE']
         print('\nPOLY_TYPE: '+rec[ind_n])
         print('Forms:')
         print(MIZ_forms[m][1])
      #return

   ####################################################
   # plot outlines of polygons (colour-coded)
   # for key in form_cats_MIZ:
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
         n        = MIZ_forms[tval][0]
         dct      = MIZ_forms[tval][1]

         figname  = Fname+'_MIZpoly'+str(n)+'.png'

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
         bm    = Mplt.start_map(bbox,cres='f')
         print('bbox=')
         print(bbox)

      else:
         # plot all polygons in MIZ
         to_plot  = range(len(MIZ_forms))

         PLOT_COMBINED  = 1 # merge neighbouring MIZ polygons together
         if PLOT_COMBINED==0:
            figname  = Fname+'_MIZ_'+key+'.png'
         else:
            figname  = Fname+'_MIZ_'+key+'2.png'
         ttl      = 'All polygons in MIZ - using '+key

         # set plot domain from overall bbox
         bbox  = sf.bbox # [lonmin,latmin,lonmax,latmax]
         bm    = Mplt.start_map(bbox)
      ####################################################
   
      ####################################################
      Ncols = len(MIZcols)
      x_all = []
      y_all = []

      MPlist   = []
      for m in to_plot:
         n     = MIZ_forms[m][0]
         fdct  = MIZ_forms[m][1]
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
         MP = shops.unary_union(MPlist) # merge neighbouring polygons (dissolve holes)
         if not hasattr(MP,'geoms'):
            MP = shgeom.MultiPolygon(MP)
         print('len(MP)='+str(len(MP.geoms)))

         print("Doing plot now...")
         Ngrps = len(MP.geoms)
         
         ####################################################
         cols  = ['c','b','g','m','r'] # colours corresponding to each group of MIZ ice
         for m in range(Ngrps):
            Ninner   = len(MP.geoms[m].interiors)
            print('No of interior regions')
            print([m,Ninner])
            #
            m_ = np.mod(m,len(cols))
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
   
      Mplt.finish_map(bm)
      plt.title(ttl,y='1.06')
      plt.savefig(figname)
      plt.close()
      fig.clf()
      print('Saving figure to '+figname)
      #########################################################

      if PLOT_COMBINED==1:
         # save MP to shapefile
         ofil  = 'test_outputs/test.shp'
         print('\nSaving MultiPolygon to shapefile: '+ofil+'\n')
         SFU.MultiPolygon2ShapeFile(MP,ofil)

         if 1:
            #########################################################
            # test shapefile
            print('Testing write to '+ofil+'...')
            sf2      = shapefile.Reader(ofil)
            sf2_info = SFU.extract_shapefile_info(sf2,get_holes=True)
            #
            fig2  = plt.figure()
            bm2   = Mplt.start_map(sf2.bbox)
            cols  = ['c','b','g','m','r'] # colours corresponding to each group of MIZ ice

            for m in range(sf2.numRecords):
               poly  = sf2_info[m][0]
               rec   = sf2_info[m][1]

               print('Looking at polygon '+str(m))
               for key2 in rec.keys():
                  print(key2+' = '+str(rec[key2]))

               Ninner   = len(poly.interiors)
               print('No of interior regions: '+str(Ninner)+'\n')
               m_    = np.mod(m,len(cols))
               col   = cols[m_]
               SFU.plot_poly(poly,pobj=bm2,plot_holes=True,color=col,latlon=True)

            Mplt.finish_map(bm2)
            figname2 = 'test_outputs/test_shp.png'
            plt.savefig(figname2)
            plt.close()
            fig2.clf()
            print('Saving figure to '+figname2+'\n')
            #########################################################
