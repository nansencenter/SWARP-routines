import sys,os
import numpy as np
from matplotlib import pyplot as plt

pdir  = '../py_funs'
if pdir not in sys.path:
   sys.path.append(pdir)

import mod_reading   as Mrdg
import fns_plotting  as Mplt

##################################################################
# read in arrays from TP4 binaries (.a)
# lon/lat from grid binary
tdir  = 'eg_TP4files'
gfil  = tdir+'/regional.grid.a'
plon  = Mrdg.get_array_from_HYCOM_binary(gfil,1)
nx,ny = plon.shape

if 0:
   plt.imshow(plon)
   plt.colorbar()
   plt.show()
   sys.exit('L22')

plat  = Mrdg.get_array_from_HYCOM_binary(gfil,2,dims=(nx,ny))

# 
afil  = tdir+'/TP4archv_wav.2015_139_180000.a'
vlst  = {'ficem':1,'hicem':2,'dfloe':3,'swh':4,'mwp':5}
vbl   = 'ficem'
recno = vlst[vbl]
V     = Mrdg.get_array_from_HYCOM_binary(afil,recno,dims=(nx,ny))
##################################################################

if 1:
   # plot V with imshow:
   plt.imshow(V)
   plt.title(vbl)
   plt.colorbar()
   plt.show()
elif 1:
   # plot V with basemap:
   fig   = plt.figure()
   bm = Mplt.start_HYCOM_map('TP4')

   if 0:
      # plot conc
      bm.pcolor(plon,plat,V,vmin=0.,vmax=1.,latlon=True)
      bm.colorbar()
   else:
      # binary of ice:
      B           = np.zeros((nx,ny))
      B[V>0.15]   = 1.
      bm.pcolor(plon,plat,B,vmin=0.,vmax=1.,latlon=True)
      bm.colorbar()

      # skimage info: http://scikit-image.org/docs/dev/auto_examples/plot_contours.html
      from skimage import measure
      contours = measure.find_contours(B, 0.)
      for n,contour in enumerate(contours):
         print('contour: '+str(n))
         print('number of points: '+str(len(contour))+'\n')
         bm.plot(plon[contour[:,0],contour[:,1]],plat[contour[:,0],contour[:,1]],'c',linewidth=1.5,latlon=True)
   
   plt.title(vbl)
   Mplt.finish_map(bm)

   if 0:
      plt.show()
   else:
      plt.savefig('out/TP4test.png')
      plt.close()
      fig.clf()
