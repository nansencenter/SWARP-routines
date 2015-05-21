import sys,os
import numpy as np
from matplotlib import pyplot as plt

pdir  = '../py_funs'
if pdir not in sys.path:
   sys.path.append(pdir)

import mod_reading   as Mrdg
import fns_plotting  as Mplt

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

if 1:
   plt.imshow(V)
   plt.title(vbl)
   plt.colorbar()
   plt.show()
else:
   # plot V:
   fig   = plt.figure()
   bm = Mplt.start_HYCOM_map('TP4')
   bm.pcolor(plon,plat,V)
   bm.colorbar()
   plt.title(vbl)
   Mplt.finish_map(bm)

   plt.savefig('out/TP4test.png')
   fig.clf()
