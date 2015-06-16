import sys,os
import numpy as np
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
import time

rootdir  = os.getenv('SWARP_ROUTINES')
rootdir  = '/Users/gianot/SWARP-routines'
sys.path.append(rootdir+'/py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution as Leqs
import f_vals_smoother as valsm
# basemap (OSISAF grid projection)
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
      projection='stere',lat_ts=70,lat_0=90,lon_0=-45)

if 1:
   # diff between OSISAF and model
   nt = 55
else:
   # MIZ from conc
   nt = 70

npfil = 'npz/poly'+str(nt)+'.npz'
if not os.path.exists(npfil):
   poly_test   = poly_list[nt]
   xy_coords   = poly_test.xy_list
   xy_coords2  = [tuple(xyc) for xyc in xy_coords]
   fvals2      = 1*poly_test.function_vals

   # save file
   np.savez(npfil,xy=xy_coords,func_vals=fvals2)
else:
   npz         = np.load(npfil)
   xy_coords2  = [tuple(xyc) for xyc in npz['xy']]
   fvals2      = npz['func_vals']

# TODO call AOD code here (or MIZ code) (instead of loading saved files)
# TODO move "smoother" function here

t0 = time.clock()
print('\n**********************************************************************')
print('Calculating potential...\n')
fun_sol  = Leqs.dirichlet_fund_soln(xy_coords2,fvals2)#,bmap=hqm)
t1 = time.clock()
print('\nTime to get potential (s): '+str(t1-t0))
print('**********************************************************************\n')

print('\n**********************************************************************')
print('Calculating stream function...\n')
stream   = Leqs.dirichlet_stream_func(potential=fun_sol)
t2 = time.clock()
print('\nTime to get stream function (s): '+str(t2-t1))
print('**********************************************************************\n')

###########################################################################################
#add a test function to eliminate contours that don't cross from "0" to "1":
class contour_selection:

   def __init__(self,func_vals):
      self.func_vals = 1*func_vals
      return

   def selector_binary(self,i0,il):
      # assume func_vals are a binary function
      # remove contours that end on the opposite value
      # OR the "unknown" part (between 0,1)
      f0 = self.func_vals[i0]
      fl = self.func_vals[il]
      if f0==1.:
         keep  = (fl<1.)
      elif f0==0.:
         keep  = (fl>0.)
      else:
         keep  = False
      return keep

   def selector_binary_v2(self,i0,il):
      # assume func_vals are a binary function
      # remove contours that end on the same value (or group of values)
      f0 = self.func_vals[i0]
      fl = self.func_vals[il]
      if f0==1.:
         keep  = (fl<1.)
      elif f0==0.:
         keep  = (fl>0.)
      else:
         keep  = ((fl==0.) or (fl==1.))
      return keep
###########################################################################################


CS = contour_selection(fun_sol.func_vals)
if 1:
   selector_function = CS.selector_binary
   fstr              = '_v1'
elif 0:
   selector_function = CS.selector_binary_v2
   fstr              = '_v2'
else:
   selector_function = None
   fstr              = ''

outdir   = 'out/test_Lap'
if not os.path.exists(outdir):
   os.mkdir(outdir)
outdir   = outdir+'/poly'+str(nt)
if not os.path.exists(outdir):
   os.mkdir(outdir)

print('\n**********************************************************************')
print('Getting streamlines...\n')
if 0:
   # euclidean space
   AI = stream.get_contour_lengths(pobj=plt,show=False,\
                                   test_function=selector_function,\
                                   func_vals_orig=1*fun_sol.func_vals)
   figname  = outdir+'/test_Laplacian_planar'+fstr+'.png'
   print('Saving plot to figure '+figname)
else:
   # spherical stuff
   AI = stream.get_contour_lengths(pobj=plt,bmap=hqm,show=False,\
                                   test_function=selector_function,\
                                   func_vals_orig=1*fun_sol.func_vals)
   figname  = outdir+'/test_Laplacian_spherical'+fstr+'.png'
   print('Saving plot to figure '+figname)

ttl   = 'Median length (km) '+str(np.round(10.*AI.length_median/1.e3)/10.)
plt.title(ttl)
plt.savefig(figname)
plt.close()

t3 = time.clock()
print('\nTime to get streamlines (mins): '+str((t3-t2)/60.))
print(ttl)
print('**********************************************************************\n')
