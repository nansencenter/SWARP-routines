import sys,os
import numpy as np
from matplotlib import pyplot as plt

pdir  = '../py_funs'
if pdir not in sys.path:
   sys.path.append(pdir)

import mod_reading   as Mrdg
import fns_plotting  as Mplt

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

if 0:
   vbl   = 'ficem'
   Vmax  = 1.
   Vmin  = 0.
   Bc    = 0.  # look for contours of this binary array from getbin
   Bstr  = 'ice edge'
   Ccol  = 'c'

   def getbin(V):
      B        = np.zeros(V.shape)
      B[V>.15] = 1.
      return B
elif 0:
   vbl   = 'hicem'
   Vmax  = 5.
   Vmin  = 0.
   Bc    = np.nan  # look for contours of this binary array from getbin
   Ccol  = 'k'

elif 1:
   vbl   = 'dfloe'
   Vmax  = 300.
   Vmin  = 0.
   Bc    = 0.  # look for contours of this binary array from getbin
   Bstr  = 'MIZ'
   Ccol  = 'k'

   def getbin(V):
      B        = np.zeros(V.shape)
      B[np.logical_and(V>0.,V<200.)]   = 1.
      return B

print('*******************************************************')
print('Looking at '+vbl)
if not np.isnan(Bc):
   print('>>> Looking for contours of type: '+Bstr)
print('*******************************************************'+'\n')

##################################################################
# read in arrays from TP4 binaries (.a)
# lon/lat from grid binary
tdir  = 'eg_TP4files'
plon,plat,depths,nx,ny  = get_grid(tdir)
# 
afil  = tdir+'/TP4archv_wav.2015_139_180000.a'
bfil  = afil[:-2]+'.b'
if 0:
   vlst  = {'ficem':1,'hicem':2,'dfloe':3,'swh':4,'mwp':5}
else:
   # make dictionary from bfil to get record number in afil
   vlst  = Mrdg.get_record_numbers_HYCOM(bfil)

# get recno from dictionary created from bfil
recno = vlst[vbl]
V     = Mrdg.get_array_from_HYCOM_binary(afil,recno,dims=(nx,ny))
#print('recno='+str(recno))
##################################################################

if 0:
   # plot V with imshow:
   plt.imshow(V)
   plt.title(vbl)
   plt.colorbar()
   plt.show()
elif 1:
   # plot V with basemap:
   fig   = plt.figure()
   bm = Mplt.start_HYCOM_map('TP4')

   contours = []
   if not np.isnan(Bc):
      # binary of ice:
      B  = getbin(V)

      # contours with skimage:
      # skimage info: http://scikit-image.org/docs/dev/auto_examples/plot_contours.html
      from skimage import measure
      contours = measure.find_contours(B,Bc)

   if 1:
      # plot conc
      bm.pcolor(plon,plat,V,vmin=Vmin,vmax=Vmax,latlon=True)
      bm.colorbar()
   else:
      # plot binary
      bm.pcolor(plon,plat,B,vmin=0.,vmax=1.,latlon=True)
      bm.colorbar()

   Nconts   = len(contours)
   if Nconts>0:
      print('number of contours: '+str(Nconts)+'...\n')
      for n,contour in enumerate(contours):
         # print('contour: '+str(n))
         # print('number of points: '+str(len(contour))+'\n')
         ivec  = np.array(contour[:,0],dtype=int)
         jvec  = np.array(contour[:,1],dtype=int)
         bm.plot(plon[ivec,jvec],plat[ivec,jvec],Ccol,linewidth=1.5,latlon=True)
   
   plt.title(vbl)
   Mplt.finish_map(bm)

   if 0:
      plt.show()
   else:
      figname  = 'out/TP4test_'+vbl+'.png'
      print('\nsaving to '+figname+'...\n')
      plt.savefig(figname)
      plt.close()
      fig.clf()
