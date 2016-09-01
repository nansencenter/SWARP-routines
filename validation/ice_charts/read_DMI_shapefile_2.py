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

# shapely info:
# official docs: http://toblerity.org/shapely/shapely.geometry.html
# more demos: http://toblerity.org/shapely/manual.html
import shapely.geometry as shgeom
import shapely.ops      as shops

import shapefile_utils  as SFU
import fns_plotting     as Mplt
from getopt import getopt

indir          = None
outdir         = None
MIZ_criteria   = "FA_only"

opts,args   = getopt(sys.argv[1:],"",["MIZ_criteria=","indir=","outdir="])
for opt,arg in opts:
   if opt=='--indir':
      indir = arg
   if opt=='--outdir':
      outdir   = arg
   if opt=='--MIZ_criteria':
      MIZ_criteria   = arg

if indir is None:
   raise ValueError('Specify input dir with --indir=')
if outdir is None:
   raise ValueError('Specify output dir with --outdir=')

class shapes:
   def __init__(self,sf_info):
      self.polygons  = []
      self.records   = []
      for poly,lut in sf_info:
         self.polygons.append(poly)
         self.records .append(lut)
      return


def plot_patches(ax,MPlist,col):
   patch_list  = []

   for poly in MPlist:
      ccl   = list(poly.exterior.coords)
      patch_list.append(patches.Polygon(ccl,True))

   pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=.5)
   pc.set_facecolor(col)
   ax.add_collection(pc)
   return

class MIZ_from_shapefile:
   def __init__(self,fname_full,MIZ_criteria="FA_only"):
      """
      MIZshp   = MIZ_from_shapefile(fname_full,MIZ_criteria="FA_only")
      """
      sf    = shapefile.Reader(fname_full)
      print('Processing '+fname_full+'...')

      # ============================================================
      # get shapes and records from shapefile
      # as list of lists:
      # sf_info[n]   = [poly,record]
      # where poly is a shapely polygon,
      # and record is all the metadata
      sf_info     = SFU.extract_shapefile_info(sf,get_holes=True)
      self.shapes = shapes(sf_info)
      self.Npolys = len(sf_info)

      print('\nExample record:')
      print(self.shapes.records[0])
      print('\n')
      # ============================================================


      # ================================================================================
      # MIZ definition stuff
      if MIZ_criteria=="FA_only":
         form_cats_MIZ     = ['FA']
         match_all_crits   = True
         #  in MIZ if (FA<=MIZthresh)
      elif MIZ_criteria=="FB_only":
         form_cats_MIZ     = ['FB']
         match_all_crits   = True
         #  in MIZ if (FB<=MIZthresh)
      elif MIZ_criteria=="FC_only":
         form_cats_MIZ     = ['FB']
         match_all_crits   = True
         #  in MIZ if (FC<=MIZthresh)
      else:
         form_cats_MIZ  = ['FA','FB','FC']
         if MIZ_criteria=="all_Fcats_small":
            match_all_crits   = True
            #  in MIZ if (FA<=MIZthresh) and (FB<=MIZthresh) and (FC<=MIZthresh)
         elif MIZ_criteria=="any_Fcats_small":
            match_all_crits   = False
            #  in MIZ if (FA<=MIZthresh) or (FB<=MIZthresh) or (FC<=MIZthresh)
         else:
            raise ValueError("unknown option for 'MIZ_criteria'\n"+\
                              "- valid options: 'FA_only', 'FB_only', 'FC_only',\n"+\
                              " 'all_Fcats_small' or 'any_Fcats_small'")

      form_vals         = SFU.DMI_form_dictionary()
      MIZthresh         = 4
      self.MIZ_criteria = MIZ_criteria
      # ================================================================================


      # ============================================================
      # separate water, MIZ and pack:
      # find MIZ from values of FA,FB,FC
      self.MIZ_forms    = [] # MIZ polygons
      self.WTR_forms    = [] # water polygons
      self.PACK_forms   = [] # pack ice polygons

      print('\nTesting if polygons are in the MIZ...')
      for n,lut in enumerate(self.shapes.records):

         # ========================================================
         # Try to define MIZ:
         fdct  = {}
         bools = np.zeros((len(form_cats_MIZ)))
         for i,key in enumerate(form_cats_MIZ):
            sval  = lut[key]        # string 
            val   = form_vals[sval] # integer from form_vals dictionary 
            fdct.update({key:val})
            
            bools[i] = val<=MIZthresh

         if match_all_crits:
            # define MIZ if all of the forms are <= threshhold
            isMIZ = np.all(bools)
         else:
            # define MIZ if any of the forms are <= threshhold
            isMIZ = np.any(bools)
         # ========================================================

         if isMIZ:
            self.MIZ_forms.append([n,fdct])
         elif lut['POLY_TYPE']=='I':
            # not MIZ but ice
            self.PACK_forms.append(n)
         else:
            # not ice
            self.WTR_forms.append(n)
      # ============================================================


      # merge water,pack and MIZ polygons
      self.merge_polygons()

      # classify MIZ boundary
      self.classify_MIZ_bdy()

      # determine which of the original shapes
      # correspond to which of the merged polygons
      # - can then get areas of each merged one
      self.get_areas()

      return

   def write_text_file(self,tfil1,MIZtype=None):
      """
      write_text_file(self,textfil,MIZtype=None)
      MIZtype='inner','outer' or None (both in/out)
      """
      # ============================================================
      # write text files
      # tfil2 = Fname+'_MIZpolys_info.txt'
      print('\nSaving polygons + boundary flags to '+tfil1)
      # print('Saving metadata to '+tfil2)
      print('\n')
      tf1   = open(tfil1,'w')
      # tf2   = open(tfil2,'w')
      blk   = 3*' '

      Count = 0
      for npoly,dct in enumerate(self.MIZclass):
         npts  = len(dct['lon'])

         if MIZtype is not None:
            if dct['type']!=MIZtype:
               continue

         # ==================================================
         # headers
         if Count==0:
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
            ss0   = [Count,\
                     dct['lon'][npt],
                     dct['lat'][npt],
                     dct['flags'][npt]]

            ss = ''
            for val in ss0[:-1]:
               ss   += str(val)+blk
            ss   += str(ss0[-1])+"\n"
            tf1.write(ss)
         # ==================================================

         Count+= 1

      tf1.close()
      # tf2.close()
      # ============================================================
      return

   def merge_polygons(self):
      # ============================================================
      # merge pack polygons (also dissolving holes)
      MPlist_PACK = []
      Npack       = len(self.PACK_forms)

      print('Merging '+str(Npack)+' pack polygons...\n')
      for n in self.PACK_forms:
         MPlist_PACK.append(self.shapes.polygons[n])

      self.MP_PACK  = shops.unary_union(MPlist_PACK)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_PACK,'geoms'):
         # single polygon
         self.MPlist_PACK  = [self.MP_PACK]
      else:
         # many polygons
         self.MPlist_PACK  = []
         for poly in self.MP_PACK.geoms:
            self.MPlist_PACK.append(poly)
      # ============================================================


      # ============================================================
      # merge wtr polygons (also dissolving holes)
      MPlist_WTR  = []
      Nwtr        = len(self.WTR_forms)

      print('Merging '+str(Nwtr)+' water polygons...\n')
      for n in self.WTR_forms:
         MPlist_WTR.append(self.shapes.polygons[n])

      self.MP_WTR = shops.unary_union(MPlist_WTR)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_WTR,'geoms'):
         # single polygon
         self.MPlist_WTR   = [self.MP_WTR]
      else:
         # many polygons
         self.MPlist_PACK  = []
         for poly in self.MP_PACK.geoms:
            self.MPlist_PACK.append(poly)
      # ============================================================


      # ============================================================
      # merge MIZ polygons (also dissolving holes)
      MPlist_MIZ  = []
      N_MIZ       = len(self.MIZ_forms)

      print('Merging '+str(N_MIZ)+' MIZ polygons...\n')
      for n,fdct in self.MIZ_forms:
         MPlist_MIZ.append(self.shapes.polygons[n])

      self.MP_MIZ = shops.unary_union(MPlist_MIZ)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_MIZ,'geoms'):
         # 1 polygon
         self.MPlist_MIZ  = [self.MP_MIZ]
      else:
         self.MPlist_MIZ  = []
         for poly in self.MP_MIZ.geoms:
            self.MPlist_MIZ.append(poly)
      # ===============================================================
      return
   
   def classify_MIZ_bdy(self):

      # ==================================================================
      # TODO match to original shapes to extract correct areas
      # TODO are holes important?
      N_MIZ_merged   = len(self.MPlist_MIZ)
      print('Sorting '+str(N_MIZ_merged)+' MIZ polygons...\n')

      self.MIZclass  = []
      count          = 0
      for poly in self.MPlist_MIZ:
         count += 1
         print("Progress: "+str(count)+"/"+str(N_MIZ_merged))
         llc            = list(poly.exterior.coords)
         lonc,latc      = poly.exterior.coords.xy
         Nc             = len(llc)
         boundary_flags = Nc*[2] # 2 = unknown, 0 = ice edge, 1 = MIZ-pack edge

         # =================================================
         # ice edge
         buf   = 1.e-5
         if poly.exterior.intersects(self.MP_WTR.buffer(buf)):
            IE = poly.exterior.intersection(self.MP_WTR.buffer(buf))
               # coords in ice edge

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
         if poly.exterior.intersects(self.MP_PACK.buffer(buf)):
            IE = poly.exterior.intersection(self.MP_PACK.buffer(buf))
               # coords in MIZ-pack edge

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


         # ============================================================
         # finalise outputs
         # - which list to append poly to?
         # (inner MIZ or outer/proper MIZ)
         for flag in range(3):
            print('#'+str(flag)+': '+str(boundary_flags.count(flag)))
         print('Total number of points: '+str(Nc)+"\n")

         Dct   = {'lon'    :np.array(lonc),
                  'lat'    :np.array(latc),
                  'flags'  :np.array(boundary_flags)}

         if boundary_flags.count(0)>0:
            # does have an ice-water boundary
            Dct.update({'type':'outer'})
         else:
            # doesn't have an ice-water boundary
            Dct.update({'type':'inner'})

         self.MIZclass.append(Dct)
         # ============================================================
      return
      # ============================================================

   # ============================================================
   def get_areas(self):
   
      Nmiz                    = len(self.MPlist_MIZ)
      self.MIZ_areas          = np.zeros((Nmiz))
      self.MIZ_forms_sorted   = Nmiz*[[]]
      for i,Poly in enumerate(self.MPlist_MIZ):
         for n,fdct in self.MIZ_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.MIZ_areas[i]   += lut['AREA']
               self.MIZ_forms_sorted[i].append(n)


      Npack                   = len(self.MPlist_PACK)
      self.PACK_areas         = np.zeros((Npack))
      self.PACK_forms_sorted  = Npack*[[]]
      for i,Poly in enumerate(self.MPlist_PACK):
         for n in self.PACK_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.PACK_areas[i]   += lut['AREA']
               self.PACK_forms_sorted[i].append(n)


      Nwater                  = len(self.MPlist_WTR)
      self.WTR_areas         = np.zeros((Nwater))
      self.WTR_forms_sorted  = Npack*[[]]
      for i,Poly in enumerate(self.MPlist_WTR):
         for n in self.WTR_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.WTR_areas[i]   += lut['AREA']
               self.WTR_forms_sorted[i].append(n)

      return
   # ============================================================


   def test_plot(self,figname=None):
      fig      = plt.figure(figsize=(20,20))
      ax       = fig.add_subplot(111)

      # ====================================================
      # limits of figure (not automatic with patches)
      x0,y0,x1,y1 = self.MP_WTR .bounds
      X0,Y0,X1,Y1 = self.MP_PACK.bounds
      x0 = min(x0,X0)
      y0 = min(y0,Y0)
      x1 = max(x1,X1)
      y1 = max(y1,Y1)
      X0,Y0,X1,Y1 = self.MP_MIZ.bounds
      x0 = min(x0,X0)
      y0 = min(y0,Y0)
      x1 = max(x1,X1)
      y1 = max(y1,Y1)
      ax.set_xlim([x0,x1])
      ax.set_ylim([y0,y1])
      ax.set_xlabel('Longitude, $^\circ$E',fontsize=18)
      ax.set_ylabel('Latitude, $^\circ$N',fontsize=18)
      for li in ax.get_xticklabels():
         li.set_fontsize(18)
      for li in ax.get_yticklabels():
         li.set_fontsize(18)
      # ====================================================


      # ====================================================
      plot_patches(ax,self.MPlist_WTR,'c')
      plot_patches(ax,self.MPlist_PACK,'r')

      # distinguish between "inner" and "outer" MIZ
      ML1   = []
      ML2   = []
      for i,poly in enumerate(self.MPlist_MIZ):
         if self.MIZclass[i]['type']=='outer':
            ML1.append(poly)
         else:
            ML2.append(poly)

      plot_patches(ax,ML1,'y')
      plot_patches(ax,ML2,'g')
      # ====================================================


      # ====================================================
      # plot ice edges (flag=0) and inner edges (flag=1)
      lstils   = ['.b','.r']
      for MIZi in self.MIZclass:
         Lon   = MIZi['lon']
         Lat   = MIZi['lat']
         ax.plot(Lon,Lat,':')
         for flag in range(2):
            lon   = Lon[MIZi['flags']==flag]
            lat   = Lat[MIZi['flags']==flag]
            ax.plot(lon,lat,lstils[flag])
      # ====================================================


      # ====================================================
      # plot all original MIZ polygons
      for n,fdct in self.MIZ_forms:
         poly     = self.shapes.polygons[n]
         lon,lat  = poly.exterior.coords.xy
         ax.plot(lon,lat,':k')
      # ====================================================


      # ====================================================
      if figname is not None:
         print('Saving '+figname+'...\n')
         fig.savefig(figname,bbox_inches='tight')
      else:
         plt.show(fig)
      # ====================================================

      ax.cla()
      plt.close(fig)
      return

############################################################################
# indir    = 'test_inputs'
# indir    = 'test_inputs2'
fnames   = os.listdir(indir)
snames   = []

# outdir   = 'test_outputs'
if not os.path.exists(outdir):
   os.mkdir(outdir)

for fname in fnames:
   ext   = os.path.splitext(fname)
   if ext[1]=='.shp':
      snames.append(ext[0])



# # string values (for plotting later)
# MIZvals  = ['X'] # values of forms which will be in MIZ (floe size < 500m)
# for n in range(MIZthresh+1):
#    MIZvals.append(str(n))
# MIZcols     = ['c','k','g','m','b','r'] # colours corresponding to each form categatory 


for fname in snames:
   fname_full  = indir+"/"+fname+'.shp'
   print('\nOpening '+fname_full+'...')
   
   MIZshp   = MIZ_from_shapefile(fname_full,MIZ_criteria=MIZ_criteria)

   # ==================================================================
   #  make a figure
   cdate = fname[:8]
   Fname = outdir+'/'+cdate+'/'+fname
   if not os.path.exists(outdir):
      os.mkdir(outdir)
   if not os.path.exists(outdir+'/'+cdate):
      os.mkdir(outdir+'/'+cdate)
   # MIZshp.test_plot()
   MIZshp.test_plot(figname=Fname+'.png')
   # ==================================================================


   # ============================================================
   # write text files
   tfil1 = Fname+'_MIZpolys.txt'
   MIZshp.write_text_file(tfil1)

   tfil1 = Fname+'_MIZpolys_out.txt'
   MIZshp.write_text_file(tfil1,MIZtype='outer')
   # ============================================================
