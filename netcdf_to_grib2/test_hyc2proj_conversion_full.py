import os,sys
import numpy as np
from datetime import datetime
#
import pygrib
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d
#
SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import mod_reading as m_rdg
import mod_grib2_setup as m_g2s

ncgv        = m_rdg.nc_get_var
KEEP_MASK   = 1

#######################################################################
# test: all vbls's, all times
# - no plotting
#######################################################################

#######################################################################
# file inputs:
# ncfil = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc" # hyc2proj netcdf
ncfil = "test_ncfiles/SWARPwavesice_forecast_start20150723T000000Z.nc" # hyc2proj netcdf
print('\nConverting '+ncfil+'\n')

# output file:
outdir  = 'out'
if not os.path.exists(outdir):
   os.makedirs(outdir)
fil_out = outdir+'/test_hyc2proj_to_grib2.grb2'

# do conversion
ncinfo   = m_g2s.hyc2proj_to_grib2(ncfil,fil_out,KEEP_MASK=KEEP_MASK)
#######################################################################

#######################################################################
Nt             = ncinfo.number_of_time_records
vbl_list       = ncinfo.variable_list
time_indices   = range(Nt)
Nv             = len(vbl_list)
first          = 1

gr       = pygrib.open(fil_out)
grbmsgs  = gr.read() # get list of all messages:

diff_arr = []

for vindx in range(Nv):
   for tindx in range(Nt):
      time_index     = time_indices[tindx]
      vbl_name       = vbl_list    [vindx]
      #
      print('Test write to '+fil_out)
      print("Doing test of encoding of variable "+vbl_name)
      print("at time record number "+str(time_index))
      print(' ')
      data  = ncgv(ncfil,vbl_name,time_index)
      if KEEP_MASK==1:
         data_arr = data[:,:] # masked array
      else:
         data_arr = data[:,:].data # array

      grb_index   = Nv*time_index+vindx
      msg         = grbmsgs[grb_index]
      if first:
         # print extra info (grid, shape of earth etc)
         first = 0
         grb   = g2d(msg.tostring(),gribmsg=True)
         print(grb)
      print(msg)
      print('\n')

      #################################################################
      # compare data arrays
      print('Comparing data arrays between netcdf and grib2 files...\n')

      if KEEP_MASK:
         Data_arr = data_arr.data
         nmask    = data_arr.mask
         gdat     = msg.data()[0].data
         gmask    = msg.data()[0].mask

         # get min/max for plotting:
         clim     = [msg.data()[0].min(),
                     msg.data()[0].max()]
      else:
         Data_arr = data_arr
         nmask    = data[:,:].mask
         gdat     = msg.data()[0]
         gmask    = data[:,:].mask

         # get min/max for plotting:
         clim     = [gdat[not np.isnan(gdat)].min(),
                     gdat[not np.isnan(gdat)].max()]
      #################################################################

      #################################################################
      # print diff between arrays:
      Darr     = Data_arr[nmask==False]
      Gdat     = gdat    [gmask==False]
      ddat     = abs(Gdat-Darr)
      max_diff = ddat.max()
      print('>> Max diff in data arrays: '+str(max_diff))
      print('\n')
      #################################################################

      diff_arr.append([tindx,vbl_name,max_diff])

# finished test
gr.close()
#######################################################################

print('\n#########################################################')
print('Time index,variable name,max difference')
print('#########################################################')
for darr in diff_arr:
   if darr[0]==0:
      print('\n')
   print(darr)
