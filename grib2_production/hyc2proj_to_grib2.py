import os,sys
import numpy as np
from datetime import datetime,timedelta
#
import pygrib
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d
#
import mod_reading as tw_rdg
import mod_grib2_setup as tw_g2s

##########################################################
# read proj.in file

# grid:
pfil     = 'inputs/proj.in'
proj_in  = tw_rdg.read_proj_infile(pfil)

# read hyc2proj netcdf:
ncfil    = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc"
ncgv     = tw_rdg.nc_get_var
lon      = ncgv(ncfil,'longitude')
lat      = ncgv(ncfil,'latitude')
time_vec = ncgv(ncfil,'time')
# xx       = ncgv(ncfil,'x')
# yy       = ncgv(ncfil,'y')

data  = ncgv(ncfil,'fice',0)
# hice  = nc.variables['hice']

nx,ny = lon.shape


reftime_h   = time_vec[0] # hours since ref point
time_info   = time_vec.units.split()
time_info   = [ time_info[i] for i in [0,2] ]
time_fmt    = '%Y-%m-%dT%H:%M:%SZ'
refpoint    = datetime.strptime(time_info[1],time_fmt)
reftime     = refpoint+timedelta(hours=reftime_h)
reftime_sig = 'start of forecast'

# Outfile
fil_out = 'test_hycproj_to_grib2.grb2'
f_out    = open(fil_out,'wb')

##########################################################################################################
# INITIALISATION

# INPUT 1) discipline code (pygrib/ncepgrib2.py, grib2message)
discipline_code   = 10  #Oceanographic Products > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table0-0.shtml

# INPUT 2) identification section: pygrib/g2clib.c ("unpack1" method)
ident_sect  = tw_g2s.def_ident_sect(reftime,reftime_sig)

# CALL Grib2Encode:
grb_out     = g2e(discipline_code,ident_sect)
##########################################################################################################

##########################################################################################################
# DEFINE GRID:

# INPUT 1)
# grid definition section:
Npts     = nx*ny  # no of points
gdsinfo  = tw_g2s.def_gdsinfo(proj_in,Npts)

# INPUT 2)
# grid definition template:
gdt   = []

# grid parameters for polar stereographic projection
# > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
i_earth_shape  = 6   #spherical (assumed radius) > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-2.shtml
gdt.append(i_earth_shape)

#not applicable (parameters for spherical earth)
gdt.extend([0,0])

#not applicable (oblate spherical)
gdt.extend([0,0,0,0])

# # ??? 2 extra values compared to the table
# gdt.extend([0,0])

# no of points in each dirn
gdt.extend([nx,ny])

# lon,lat of 1st grid point:
fac   = 1e6 # degrees -> micro-degrees
lon0  = int(lon[0,0]*fac)
if lon0<0:
   lon0  = lon0+360*fac
lat0  = int(lat[0,0]*fac)
if lon0<0:
   lat0  = lat0+360*fac
gdt.extend([lat0,lon0])

# resolution and component flags > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-3.shtml
# - most not applicable (only rcf[4] - definition of reference for vector directions)
if proj_in.Rotate_vectors_to_output_grid:
   rcf[4]   = 8*'0'
   # rcf[4]=='0' => Resolved u and v components of vector quantities relative to easterly and northerly directions
else:
   rcf      = '00001000'
   # rcf[4]=='1' => Resolved u and v components of vector quantities relative to grid

#convert from binary string to integer
rcf   = int(rcf,2)
gdt.append(rcf)

# projection center
fac   = 1e6 # degrees -> micro-degrees
latc  = proj_in.Central_latitude_of_projection  # lat where resolution is specified
lonc  = proj_in.Central_longitude_of_projection # orientation = meridian which is parallel to the y-axis
latc  = int(latc*fac)
lonc  = int(lonc*fac)
if lonc<0:
   lonc  = lonc+fac*360
gdt.extend([latc,lonc])

# resolution (km)
fac   = 1e6 # km-> mm
dx    = proj_in.Projection_grid_increment_left_right  # x dirn
dy    = proj_in.Projection_grid_increment_bottom_top  # y dirn
dx    = int(fac*dx)
dy    = int(fac*dy)
gdt.extend([dx,dy])

# projection flags > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-5.shtml
pf    = 8*'0'
# pf[0]=='0' => north pole is in projection plane

#convert from binary string to integer
pf = int(pf,2)
gdt.append(pf)

# scanning mode - not applicable to stereographic grid??
scan_mode   = 8*'0'

#convert from binary string to integer
scan_mode   = int(scan_mode,2)
gdt.append(scan_mode)

# convert to numpy array
gdt   = np.array(gdt,'i')

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
pdtmpl   = [] # product definition template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp4-0.shtml

# Parameter category & Parameter number
if data.name=='fice':
   param_cat   = 2 # Discipline_code = 10: 2 -> Ice > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
   param_num   = 0   # ice cover > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
elif data.name=='hice':
   param_cat   = 2 # Discipline_code = 10: 2 -> Ice > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
   param_num   = 1   # ice thickness > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
elif data.name=='swh':
   param_cat   = 0 # Discipline_code = 10: 0 -> Waves > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
   param_num   = 4 # total significant wave height > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
elif data.name=='mwp':
   param_cat   = 0 # Discipline_code = 10: 0 -> Waves > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
   param_num   = 11# primary mean wave period > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml

pdtmpl.extend([param_cat,param_num])

# Type of generating process > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-3.shtml
type_gen_proc  = 0 # Analysis
# type_gen_proc  = 1 # Initialization
# type_gen_proc  = 2 # Forecast
pdtmpl.append(type_gen_proc)

# not applicable?
pdtmpl.extend([0,0,0,0])

# units of forecast time > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-4.shtml
ft_unit  = 0 # hours
pdtmpl.append(ft_unit)

# forecast time (time since "reftime" = start of forecast )
ft = 0
pdtmpl.append(ft)

# Type of fixed surface  > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-5.shtml
tofs  = 101 # mean sea level
pdtmpl.append(tofs)

# rest not applicable??
pdtmpl.extend(6*[0])

# convert to numpy array
pdtmpl   = np.array(pdtmpl,'i')

# INPUT 3):
# data representation template number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
drtnum  = 40 # grid point data, jpeg 2000 compression
# drtnum   = 0 # grid point data, simple packing

# INPUT 4):
# data representation template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
drtmpl   = []
if drtnum==40:
   # > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp5-40.shtml
   # drtmpl  = [1129381888,0,2,11,0,0,255] # drtnum = 40
   drtmpl      = 6*[0] # drtnum = 40
   drtmpl[0]   = 0 # lossless compression
   drtmpl[5]   = 255 # missing
else:
   # data representation template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp5-0.shtml
   drtmpl      = 5*[0]
   drtmpl[4]   = 0 # originals are floating point

   # # 10^D*y = R + 10^E*X
   # # negative values of E&D
   # # > http://meteo.ieec.uned.es:8086/PFC_JMEstepa/documentos/GRIB.pdf
   # # page: I.2-GRIB Reg. - 3 
   # # Reg 92.1.5: If applicable, negative values shall be indicated by setting the most significant bit to "1"
   # ref_val        = data.add_offset # R
   # bin_scale_fac  = np.log10(data.scale_factor) # E
   # dec_scale_fac  = 0 # D
   # drtmpl.extend([ref_val,bin_scale_fac,dec_scale_fac])

   # # bits per packed value
   # bppv  = 8
   # drtmpl.extend(bppv)

   # # original format > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-1.shtml
   # orig_fmt = 0 # floating point
   # drtmpl.extend(orig_fmt)

drtmpl   = np.array(drtmpl)

# INPUT 5)
data  = data[:,:] # masked array or array

# CALL addfield:
grb_out.addfield(pdtnum,pdtmpl,drtnum,drtmpl,data)

# finalize the grib message.
grb_out.end()

##########################################################################################################

print('Writing grib message to '+fil_out)
f_out.write(grb_out.msg)
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
      #
      # lat,lon        = msg.latlons()
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
