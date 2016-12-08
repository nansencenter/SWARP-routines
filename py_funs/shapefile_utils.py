import numpy as np
import shapefile

############################################################################
def DMI_form_dictionary():

   # dictionary to map from strings to integers
   form_vals      = {'X':-1,'-9':np.NaN}
   form_vals_max  = 10
   for n in range(form_vals_max+1):
      form_vals.update({str(n):n})
   
   return form_vals
############################################################################


############################################################################
def AARI_form_dictionary():

   # dictionary to map from strings to integers
   form_vals      = {'X':-1,'-9':np.NaN,'99':np.NaN}
   form_vals_max  = 10
   for n in range(form_vals_max+1):
      form_vals.update({"%2.2d" %(n):n})
   
   return form_vals
############################################################################


############################################################################
def form_cols(cmapname=None):
   """
   make a dictionary mapping form to rgba color
   """
   import matplotlib.pyplot as plt
   import matplotlib.cm as mplcm
   import matplotlib.colors as colors

   NUM_COLORS = 15
   if cmapname is None:
      # see http://matplotlib.org/examples/color/colormaps_reference.html
      # for a full list of colormaps
      cmapname = 'Set1'
   cm          = plt.get_cmap(cmapname)
   cNorm       = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
   scalarMap   = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
   col_list    = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]

   i     = 0
   cols  = {np.nan:col_list[i]}
   for j in range(-1,11):
      i    += 1
      cols.update({j:col_list[i]})

   return cols
############################################################################


############################################################################
def xy2list(x,y):
   
   lst   = []
   for n in range(len(x)):
      lst.append([x[n],y[n]])

   return lst
############################################################################

############################################################################
def xy2tuple_list(x,y):
   # make list of tuples from x,y
   
   lst   = []
   for n in range(len(x)):
      lst.append((x[n],y[n]))

   return lst
############################################################################

############################################################################
def get_poly(shp,get_holes=True):
   # convert points from a shapefile shape into a shapely polygon

   import shapely.geometry as shgeom

   pts   = shp.points
   Npts  = len(pts)
   #
   parts    = shp.parts
   Nparts   = len(parts)
   parts.append(Npts)

   bdy   = pts[:parts[1]]

   #########################################################################
   if get_holes and (Nparts>1):
      # add rings (holes) (if present and wanted)
      rings = []
      for n in range(1,Nparts):
         rpts  = pts[parts[n]:parts[n+1]]
         rings.append(rpts)

      poly  = shgeom.Polygon(bdy,rings)
   else:
      poly  = shgeom.Polygon(bdy)
   #########################################################################

   # try to ensure poly is valid
   return poly.buffer(0)
############################################################################

############################################################################
def plot_poly(poly,pobj=None,plot_holes=True,**kwargs):

   import shapely.geometry as shgeom

   if pobj is None:
      from matplotlib import pyplot as plt
      pobj  = plt

   # external polygon
   q     = shgeom.Polygon(poly.exterior)
   x,y   = q.boundary.coords.xy
   pobj.plot(x,y,**kwargs)

   if plot_holes:
      # interior polygons
      qi    = poly.interiors
      for q in qi:
         x,y   = shgeom.Polygon(q).boundary.coords.xy
         pobj.plot(x,y,linestyle='--',**kwargs)

   return
############################################################################

############################################################################
def plot_Mpoly(MP,pobj=None,colors=None,plot_holes=True,**kwargs):
   # plot MultiPolygon MP
   import numpy as np
   import shapely.geometry as shgeom

   if pobj is None:
      from matplotlib import pyplot as plt
      pobj  = plt
   
   col   = None
   Ng    = len(MP.geoms)

   for m,poly in enumerate(MP.geoms):
      if colors is not None:
         m_    = np.mod(m,len(colors))
         col   = colors[m_]

      L  = len(shgeom.Polygon(poly.exterior).boundary.coords)
      plot_poly(poly,pobj=pobj,color=col,plot_holes=plot_holes,**kwargs)
      if ((not poly.is_valid) and (L>50)) or (L>50):
         print("polgon "+str(m)+' ('+str(Ng)+'): '+str(len(poly.interiors))+' interiors')
         print('Valid: '+str(poly.is_valid))
         print('Length of boundary: '+str(L) +'\n')

   return
############################################################################

############################################################################
def extract_shapefile_info(sfile,get_holes=True):

   Nrec     = sfile.numRecords
   fields   = sfile.fields[1:]
   Nfields  = len(fields)

   out   = []

   for n in range(Nrec):
      shp   = sfile.shape(n)
      rec   = sfile.record(n)

      # get field values:
      lut   = {}
      for m in range(Nfields):
         fld   = fields[m][0]
         val   = rec[m]
         lut.update({fld:val})

      # get shapely polygon representation of points
      poly  = get_poly(shp,get_holes=get_holes)

      # add lookup table and shapely polygon to out
      out.append([poly,lut])

   return out
############################################################################

############################################################################
def prepMultiPolygon4ShapeFile(MP):
   # convert shapely MultiPolygon to form usable by shapefile
   # also calculate field values
   import shapely.geometry as shgeom

   Np          = len(MP.geoms)
   parts_list  = Np*[0]
   rec_list    = Np*[0]

   for n in range(Np):
      parts = []

      ############################################################################
      # record list: [AREA_ID]
      rec_list[n] = [n]
      # TODO do stereographic projection
      #      > calc area,perimeter and width and append
      #      > eg rec_list[n].append(Area)
      ############################################################################

      ############################################################################
      poly  = MP.geoms[n]
      q     = shgeom.Polygon(poly.exterior)
      x,y   = q.boundary.coords.xy
      lst   = xy2list(x,y)

      # pyshp: if points are anti-clockwise, the points are treated as exterior ones
      # TODO check if shapely does the same
      parts = [lst]

      Nint  = len(poly.interiors)
      for m in range(Nint):
         q     = shgeom.Polygon(poly.interiors[m])
         x,y   = q.boundary.coords.xy
         lst   = xy2list(x,y)

         # pyshp: if points are clockwise, the points are treated as interior ones
         # TODO check if shapely does the same
         parts.append(lst)
      ############################################################################
      
      parts_list[n]  = parts
   return parts_list,rec_list
############################################################################

############################################################################
def MultiPolygon2ShapeFile(MP,filename):
   # save MultiPolygon to shapefile
   # 1. convert shapely MultiPolygon to form usable by shapefile
   # 2. save shapefile

   import shapefile
   w  = shapefile.Writer(shapefile.POLYGON)

   # define attributes
   w.field('AREA_ID','N','40')

   # TODO calculate these and add fields for them
   # w.field('Area','N','40') # name,type ('C'=character, 'N'=number), size (?)
   # w.field('Perimeter','N','40')
   # w.field('Width','N','40')

   parts_list,rec_list  = prepMultiPolygon4ShapeFile(MP)

   Nrec  = len(parts_list)
   for n in range(Nrec):
      # add polygons
      w.poly(parts=parts_list[n])
      rec   = rec_list[n]
      w.record(AREA_ID=rec[0])

   # save file
   w.save(filename)

   return
############################################################################

############################################################################
def SKcontours2MP(contours,x,y):

   import shapely.geometry as shgeom
   import numpy as np

   plist = []
   for cont in contours:
      ivec  = np.array(cont[:,0],dtype=int)
      jvec  = np.array(cont[:,1],dtype=int)
      #
      xvec  = x[ivec,jvec]
      yvec  = y[ivec,jvec]
      tups  = xy2tuple_list(xvec,yvec)
      #
      poly  = shgeom.Polygon(tups)
      if 1:#poly.is_valid:
         plist.append(poly)

   MP = shgeom.MultiPolygon(plist)
   return MP
############################################################################


############################################################################
class shapes:
   def __init__(self,sf_info):
      self.polygons  = []
      self.records   = []
      for poly,lut in sf_info:
         self.polygons.append(poly)
         self.records .append(lut)
      return
############################################################################

class shapefile_info:
   def __init__(self,fname_full,MIZ_criteria="FA_only",chart_source="DMI"):

      # for info
      self.shapefile    = fname_full
      self.chart_source = chart_source
      if chart_source=='DMI':
         self.form_dict = DMI_form_dictionary()
      elif chart_source=='AARI':
         self.form_dict = AARI_form_dictionary()
      
      import shapefile
      sf    = shapefile.Reader(fname_full)
      print('Reading '+fname_full+'...')

      # ============================================================
      # get shapes and records from shapefile
      # self.shapes.records[n]
      # self.shapes.polygons[n]
      # give the info for the n-th shape
      # *poly is a shapely polygon
      # *record is the metadata
      sf_info     = extract_shapefile_info(sf,get_holes=True)
      self.shapes = shapes(sf_info)
      self.Npolys = len(sf_info)

      print('\nExample record:')
      print(self.shapes.records[0])
      print('\n')
      # ============================================================


      # ============================================================
      # bounds: mostly for plotting
      for i,P in enumerate(self.shapes.polygons):
         lon0,lat0,lon1,lat1  = P.bounds
         if i==0:
            self.lon_min,self.lat_min,\
                  self.lon_max,self.lat_max  = lon0,lat0,lon1,lat1
         else:
            self.lon_min  = min(self.lon_min,lon0)
            self.lat_min  = min(self.lat_min,lat0)
            self.lon_max  = max(self.lon_max,lon1)
            self.lat_max  = max(self.lat_max,lat1)
      # ============================================================

      return
   # ============================================================


   # ============================================================
   def test_plot(self,figname=None,cmapname=None):

      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      fig      = plt.figure(figsize=(12,12))
      ax       = fig.add_subplot(111)

      # ====================================================
      # limits of figure (not automatic with patches)
      ax.set_xlim([self.lon_min,self.lon_max])
      ax.set_ylim([self.lat_min,self.lat_max])
      ax.set_xlabel('Longitude, $^\circ$E',fontsize=18)
      ax.set_ylabel('Latitude, $^\circ$N',fontsize=18)
      for li in ax.get_xticklabels():
         li.set_fontsize(18)
      for li in ax.get_yticklabels():
         li.set_fontsize(18)
      # ====================================================


      # ====================================================
      clist = form_cols(cmapname=cmapname)
      for n,poly in enumerate(self.shapes.polygons):
         rec   = self.shapes.records[n]
         if rec['POLY_TYPE']=='W':
            Fplt.plot_patches(ax,[poly],'c')
         else:
            key   = self.form_dict[rec['FA']]
            col   = clist[key]
            Fplt.plot_patches(ax,[poly],col)
      # ====================================================
      

      # ====================================================
      # legend
      import matplotlib.patches as mpatches
      Patches  = [mpatches.Patch(color='c', label='W')]
      for key in sorted(clist):
         col   = clist[key]
         if np.isnan(key):
            lab   = '-'
         elif key==-1:
            lab   = 'X'
         else:
            lab   = str(key)
         Patches.append(mpatches.Patch(color=col,label=lab))

      ax.legend(handles=Patches,loc='lower right')
      # ====================================================


      # ====================================================
      # show or save fig
      if figname is None:
         plt.show(fig)
      else:
         print('\nSaving '+figname+'...\n')
         fig.savefig(figname,bbox_inches='tight')
      # ====================================================


   # ============================================================

############################################################################
class MIZ_from_shapefile:
   def __init__(self,fname_full,MIZ_criteria="FA_only",chart_source="DMI"):
      """
      MIZshp   = MIZ_from_shapefile(fname_full,MIZ_criteria="FA_only")
      """

      # for info
      self.shapefile    = fname_full
      self.chart_source = chart_source
      self.MIZ_criteria = MIZ_criteria

      sf    = shapefile.Reader(fname_full)
      print('Processing '+fname_full+'...')

      # ============================================================
      # get shapes and records from shapefile
      # as list of lists:
      # sf_info[n]   = [poly,record]
      # where poly is a shapely polygon,
      # and record is all the metadata
      sf_info     = extract_shapefile_info(sf,get_holes=True)
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

      if chart_source=="DMI":
         form_vals         = DMI_form_dictionary()
      else:
         form_vals         = AARI_form_dictionary()
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

      if self.N_MIZ==0:
         print("No MIZ polygons - not writing text file")
         return

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
      import shapely.ops      as shops
      # ============================================================
      # merge pack polygons (also dissolving holes)
      MPlist_PACK = []
      Npack       = len(self.PACK_forms)

      print('Merging '+str(Npack)+' pack polygons...')
      bad   = []
      for n in self.PACK_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_PACK.append(P)
         else:
            bad.append(n)

      for n in bad:
         self.PACK_forms.remove(n)

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

      self.N_PACK = len(self.MPlist_PACK)
      print(str(self.N_PACK)+ ' pack polygons after merging\n')
      # ============================================================


      # ============================================================
      # merge wtr polygons (also dissolving holes)
      MPlist_WTR  = []
      Nwtr        = len(self.WTR_forms)

      print('Merging '+str(Nwtr)+' water polygons...')
      bad   = []
      for n in self.WTR_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_WTR.append(P)
         else:
            bad.append(n)

      for n in bad:
         self.WTR_forms.remove(n)

      self.MP_WTR = shops.unary_union(MPlist_WTR)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_WTR,'geoms'):
         # single polygon
         self.MPlist_WTR   = [self.MP_WTR]
      else:
         # many polygons
         self.MPlist_WTR  = []
         for poly in self.MP_WTR.geoms:
            self.MPlist_WTR.append(poly)

      self.N_WTR  = len(self.MPlist_WTR)
      print(str(self.N_WTR)+ ' water polygons after merging\n')
      # ============================================================


      # ============================================================
      # merge MIZ polygons (also dissolving holes)
      MPlist_MIZ  = []
      N_MIZ       = len(self.MIZ_forms)

      print('Merging '+str(N_MIZ)+' MIZ polygons...')
      bad   = []
      for n,fdct in self.MIZ_forms:
         P  = self.shapes.polygons[n]
         if P.is_valid:
            MPlist_MIZ.append(P)
         else:
            bad.append([n,fdct])

      for nd in bad:
         self.MIZ_forms.remove(nd)

      self.MP_MIZ = shops.unary_union(MPlist_MIZ)
         # merge neighbouring polygons (dissolve holes)

      if not hasattr(self.MP_MIZ,'geoms'):
         # 1 polygon
         self.MPlist_MIZ  = [self.MP_MIZ]
      else:
         self.MPlist_MIZ  = []
         for poly in self.MP_MIZ.geoms:
            self.MPlist_MIZ.append(poly)

      self.N_MIZ  = len(self.MPlist_MIZ)
      print(str(self.N_MIZ)+ ' MIZ polygons after merging\n')
      # ===============================================================
      return
   
   def classify_MIZ_bdy(self):

      # ==================================================================
      # TODO match to original shapes to extract correct areas
      # TODO are holes important?
      N_MIZ_merged   = len(self.MPlist_MIZ)
      print('Sorting '+str(N_MIZ_merged)+' MIZ polygons...\n')

      self.MIZclass  = []
      if N_MIZ_merged==0:
         return

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
               self.MIZ_areas[i]   += float(lut['AREA'])
               self.MIZ_forms_sorted[i].append(n)


      Npack                   = len(self.MPlist_PACK)
      self.PACK_areas         = np.zeros((Npack))
      self.PACK_forms_sorted  = Npack*[[]]
      for i,Poly in enumerate(self.MPlist_PACK):
         for n in self.PACK_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.PACK_areas[i]   += float(lut['AREA'])
               self.PACK_forms_sorted[i].append(n)


      Nwater                 = len(self.MPlist_WTR)
      self.WTR_areas         = np.zeros((Nwater))
      self.WTR_forms_sorted  = Nwater*[[]]
      for i,Poly in enumerate(self.MPlist_WTR):
         for n in self.WTR_forms:
            poly  = self.shapes.polygons[n]
            lut   = self.shapes.records[n]
            if poly.within(Poly):
               self.WTR_areas[i]   += float(lut['AREA'])
               self.WTR_forms_sorted[i].append(n)

      return
   # ============================================================


   def test_plot(self,figname=None):

      from matplotlib import pyplot as plt
      import fns_plotting as Fplt

      fig      = plt.figure(figsize=(12,12))
      ax       = fig.add_subplot(111)

      # ====================================================
      # limits of figure (not automatic with patches)
      x0,y0,x1,y1 = self.MP_WTR .bounds
      X0,Y0,X1,Y1 = self.MP_PACK.bounds
      x0 = min(x0,X0)
      y0 = min(y0,Y0)
      x1 = max(x1,X1)
      if self.N_MIZ>0:
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
      Fplt.plot_patches(ax,self.MPlist_WTR,'c')
      Fplt.plot_patches(ax,self.MPlist_PACK,'r')

      # distinguish between "inner" and "outer" MIZ
      ML1   = []
      ML2   = []
      for i,poly in enumerate(self.MPlist_MIZ):
         if self.MIZclass[i]['type']=='outer':
            ML1.append(poly)
         else:
            ML2.append(poly)

      Fplt.plot_patches(ax,ML1,'y')
      Fplt.plot_patches(ax,ML2,'g')
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
