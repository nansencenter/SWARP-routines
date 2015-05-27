import sys,os
import numpy as np
from matplotlib import pyplot as plt

pdir  = '../py_funs'
if pdir not in sys.path:
   sys.path.append(pdir)

import mod_reading   as Mrdg
import fns_plotting  as Mplt

# contours with skimage:
# skimage info: http://scikit-image.org/docs/dev/auto_examples/plot_contours.html
from skimage import measure
import shapefile_utils  as SFU
import shapely.geometry as shgeom
import shapely.ops      as shops

#######################################################################
def get_grid(grid_dir='.'):

   gfil  = grid_dir+'/regional.grid.a'
   dfil  = grid_dir+'/regional.depth.a'
   #
   plon  = Mrdg.get_array_from_HYCOM_binary(gfil,1)
   nx,ny = plon.shape
   #
   plat     = Mrdg.get_array_from_HYCOM_binary(gfil,2,dims=(nx,ny))
   depths   = Mrdg.get_array_from_HYCOM_binary(dfil,1,dims=(nx,ny))

   return plon,plat,depths,nx,ny
#######################################################################

########################################################
def get_land(depths):
   B              = np.ones(depths.shape)
   B[depths>0.]   = 0
   return B
########################################################

########################################################
def get_ice(ficem):
   B           = np.ones(ficem.shape)
   B[ficem>0.15] = 0
   return B
########################################################

########################################################
def get_ice_edge(ficem,depths,plon,plat,bm):

   # get binaries
   Bi       = get_ice(ficem)
   Bl       = get_land(depths)

   # get ice boundary/coasts
   i_contours  = measure.find_contours(Bi,0.)
   l_contours  = measure.find_contours(Bl,0.)

   # project lon/lat so discontinuities etc don't confuse shapely
   X,Y   = bm(plon,plat)
   MPi   = SFU.SKcontours2MP(i_contours,X,Y) # multipolygon
   MPl   = SFU.SKcontours2MP(l_contours,X,Y) # multipolygon

   if 0:
      # check coasts
      print('check coasts')
      plt.imshow(Bl,cmap='jet')
      for contour in l_contours:
         ivec  = np.array(contour[:,0],dtype=int)
         jvec  = np.array(contour[:,1],dtype=int)
         plt.plot(jvec,ivec)  # reverse order since using imshow
      plt.show()

   elif 0:
      # check ice edge
      print('check ice edge')
      plt.imshow(Bi,cmap='jet')
      for contour in i_contours:
         ivec  = np.array(contour[:,0],dtype=int)
         jvec  = np.array(contour[:,1],dtype=int)
         plt.plot(jvec,ivec)  # reverse order since using imshow
      plt.show()

   return MPi,MPl
########################################################


##################################################################
# read in arrays from TP4 binaries (.a)
# lon/lat from grid binary
tdir  = 'eg_TP4files'
plon,plat,depths,nx,ny  = get_grid(tdir)
# 
afil  = tdir+'/TP4archv_wav.2015_139_180000.a'
bfil  = afil[:-2]+'.b'

# make dictionary from bfil to get record number in afil
vlst  = Mrdg.get_record_numbers_HYCOM(bfil)

# get recno from dictionary created from bfil
vbl   = 'ficem'
ficem = Mrdg.get_array_from_HYCOM_binary(afil,vlst[vbl],dims=(nx,ny))
vbl   = 'dfloe'
dfloe = Mrdg.get_array_from_HYCOM_binary(afil,vlst[vbl],dims=(nx,ny))
##################################################################

fig      = plt.figure()
bm       = Mplt.start_HYCOM_map('TP4')
if 1:
   bm.pcolor(plon,plat,ficem,latlon=True,vmin=0.,vmax=1.)
   plt.colorbar()

MPi,MPl  = get_ice_edge(ficem,depths,plon,plat,bm)

# simplify ice
MP = shops.unary_union(MPi)

cols  = ['b','r','m','c','g']
if 1:
   # SFU.plot_Mpoly(MPi,pobj=bm,plot_holes=True)
   SFU.plot_Mpoly(MP,pobj=bm,plot_holes=True,colors=['m'],linewidth=1.5)
elif 0:
   SFU.plot_Mpoly(MPl,pobj=bm,colors=cols,plot_holes=True)
   # SFU.plot_Mpoly(shops.unary_union(MPl),pobj=bm,colors=cols,plot_holes=True)

plist1 = []
plist2 = []
for p in MP.geoms:
   b1    = shgeom.Polygon(p.exterior).boundary
   lst2  = []
   OK    = True
   for q in MPl.geoms:
      b2    = shgeom.Polygon(q.exterior).boundary
      if 1:
         # get intersection
         isec  = b1.intersection(b2)
         if isec.length>0:
            lst2.append(isec)
            tst   = isec

            if tst.type=='LineString':
               xx,yy = tst.coords.xy
               bm.plot(xx,yy,'k',linewidth=1.5)
            elif tst.type=='MultiLineString':
               for ll in tst.geoms:
                  xx,yy = ll.coords.xy
                  bm.plot(xx,yy,'k',linewidth=1.5)
            elif tst.type=='Polygon':
               xx,yy = shgeom.Polygon(tst.exterior).boundary.coords.xy
               bm.plot(xx,yy,'k',linewidth=1.5)
            # elif tst.type=='MultiPolygon':
            #    for pp in tst.geoms:
            #       xx,yy = shgeom.Polygon(pp.exterior).boundary.coords.xy
            #       bm.plot(xx,yy,'k',linewidth=1.5)

      else:
         # get non-intersecting part of boundary
         # NB CRASHES
         ncomm = b1.intersection(b1.symmetric_difference(b2))
         if ncomm.length>0.0:
            lst2.append(ncomm)
            print('non-common type: '+ncomm.type)
            print('non-common length: '+str(ncomm.length))

            tst   = ncomm
            # # plot 
            # if tst.type=='LineString':
            #    xx,yy = tst.coords.xy
            #    bm.plot(xx,yy,'k',linewidth=1.5)
            # elif tst.type=='MultiLineString':
            #    for ll in tst.geoms:
            #       xx,yy = ll.coords.xy
            #       bm.plot(xx,yy,'k',linewidth=1.5)
            # elif tst.type=='Polygon':
            #    xx,yy = shgeom.Polygon(tst.exterior).boundary.coords.xy
            #    bm.plot(xx,yy,'k',linewidth=1.5)
            # # elif tst.type=='MultiPolygon':
            # #    for pp in tst.geoms:
            # #       xx,yy = shgeom.Polygon(pp.exterior).boundary.coords.xy
            # #       bm.plot(xx,yy,'k',linewidth=1.5)
      
   # if isec.type!='Polygon':
   #    plist1.append(p)
   #    plist2.append(lst2)


if 0:
   Mplt.finish_map(bm)

if 1:
   plt.savefig('out/testMIZ.png')
   plt.close()
   fig.clf()
else:
   plt.show()
