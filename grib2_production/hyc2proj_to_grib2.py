import os,sys
import numpy as np
from datetime import datetime
#
import pygrib
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d
#
import mod_reading as tw_rdg
import mod_grib2_setup as tw_g2s
ncgv  = tw_rdg.nc_get_var

##########################################################
# inputs - loop over these?
vbl_list       = ['fice','hice'] # TODO: change names? -> icec,icetk (also in mod_grib2_setup.py)
                                 # TODO: get vbl_list from netcdf file directly

# file inputs:
pfil  = 'inputs/proj.in' # proj.in file used by hyc2proj
ncfil = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc" # hyc2proj netcdf

# output file:
outdir  = 'out'
if not os.path.exists(outdir):
   os.makedirs(outdir)
fil_out = outdir+'/test_hycproj_to_grib2.grb2'
##########################################################

##########################################################
# read proj.in file
proj_in  = tw_rdg.read_proj_infile(pfil)

# read hyc2proj netcdf:
ncinfo         = tw_rdg.nc_getinfo(ncfil)
time_indices   = range(ncinfo.number_of_time_records)

# Open outfile before starting to encode messages
f_out = open(fil_out,'wb')

time_index  = time_indices[0]
vbl         = vbl_list[0]
data        = ncgv(ncfil,vbl,time_index)

##########################################################################################################
# INITIALISATION

# INPUT 1) discipline code (pygrib/ncepgrib2.py, grib2message)
discipline_code   = 10  #Oceanographic Products > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table0-0.shtml

# INPUT 2) identification section: pygrib/g2clib.c ("unpack1" method)
ident_sect  = tw_g2s.def_ident_sect(ncinfo.reftime,ncinfo.reftime_sig)

# CALL Grib2Encode:
grb_out     = g2e(discipline_code,ident_sect)
##########################################################################################################

##########################################################################################################
# DEFINE GRID:

# INPUT 1)
# grid definition section:
nx,ny    = ncinfo.shape
Npts     = nx*ny  # no of points
gdsinfo  = tw_g2s.def_gdsinfo(proj_in,Npts)

# INPUT 2)
# grid definition template:
# > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
gdt   = tw_g2s.def_grid_template(proj_in,ncinfo)

# CALL addgrid
grb_out.addgrid(gdsinfo,gdt)
# sys.exit('L132')
##########################################################################################################

##########################################################################################################
# add product definition template, data representation template
# and data (including bitmap which is read from data mask).

# INPUT 1):
pdtnum   = 0 # product_definition_template_number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-0.shtml

# INPUT 2):
pdtmpl   = tw_g2s.def_prod_template(ncinfo,data.name)

# INPUT 3):
# data representation template number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
drtnum  = 40 # grid point data, jpeg 2000 compression
# drtnum   = 0 # grid point data, simple packing

# INPUT 4):
# data representation template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
drtmpl   = tw_g2s.def_datarep_template(drtnum)

# INPUT 5)
data_arr = data[:,:] # masked array or array

# CALL addfield:
grb_out.addfield(pdtnum,pdtmpl,drtnum,drtmpl,data_arr)

# finalize the grib message.
grb_out.end()

print('Writing grib message to '+fil_out)
f_out.write(grb_out.msg)
##########################################################################################################

f_out.close()

if 1:
   print('Test write to '+fil_out)
   gr = pygrib.open(fil_out)

   N        = 1
   grbmsgs  = gr.read(N) # get list of 1st N messages:
   for msg in grbmsgs:
      print(msg)
      print('\n')
      grb   = g2d(msg.tostring(),gribmsg=True)
      print(grb)

      if 1:
         # compare grids:
         lon         = ncgv(ncfil,'longitude')
         lat         = ncgv(ncfil,'latitude')
         glat,glon   = msg.latlons()
         #
         dlat  = abs(glat.transpose()-lat[:,:])
         dlon  = abs(glon.transpose()-lon[:,:])
         print('max diff in lat arrays: '+str(dlat.max()))
         print('max diff in lon arrays: '+str(dlon.max()))

         # compare data arrays:
         gdat  = msg.data()[0]
         Gdat  = gdat.data[gdat.mask==False]

         Data_arr = data_arr.transpose()
         Darr     = Data_arr.data[Data_arr.mask==False]
         
         ddat  = abs(Gdat-Darr)

         # # make plots
         # from matplotlib import pyplot as plt 
         # from matplotlib import Basemap
         # fig   = plt.figure()

      # data           = msg.values
      # Z              = data.data
      # Z[data.mask]   = np.NaN
      # Zmin           = data.min()
      # Zmax           = data.max()
      # #
      # lon_0 = lon.mean()
      # lat_0 = lat.mean()
      # R1    = grb.earthRmajor

   gr.close()
