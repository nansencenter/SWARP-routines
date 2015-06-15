from netCDF4 import Dataset
import sys,os
import glob
import numpy as np
import subprocess
import shutil
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from skimage import measure as msr
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

rootdir  = os.getenv('SWARP_ROUTINES')
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

fun_sol  = Leqs.dirichlet_fund_soln(xy_coords2,fvals2)#,bmap=hqm)
stream   = Leqs.dirichlet_stream_func(potential=fun_sol)
if 0:
   # euclidean space
   merged_levels  = stream.get_contour_lengths(pobj=plt)
   plt.savefig('test_Laplacian_planar.png')
else:
   # spherical stuff
   merged_levels  = stream.get_contour_lengths(pobj=plt,bmap=hqm,show=False)
   plt.savefig('test_Laplacian_spherical.png')
