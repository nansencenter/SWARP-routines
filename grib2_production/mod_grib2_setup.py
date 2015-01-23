import os,sys
import numpy as np
from datetime import datetime,timedelta

####################################################################################################
def def_ident_sect(reftime,reftime_sig):
   ident_sect  = np.array(13*[0],'i')

   # info for center
   i_centre          = 255 #missing
   i_subcentre       = 255 #missing
   ident_sect[0:2]   = [i_centre,i_subcentre]

   #Table version
   i_table_version         = 2 # version of master tables
   i_local_table_version   = 0 # local tables not used
   ident_sect[2:4]         = [i_table_version,i_local_table_version]

   # ref time
   if reftime_sig=='start of forecast':
      i_sig_reftime  = 1#start of forecast > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table1-2.shtml
   rt                = reftime
   ident_sect[4:11]  = [i_sig_reftime,rt.year,rt.month,rt.day,rt.hour,rt.minute,rt.second]

   # extra info
   iProd_status      = 1#operational test product > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table1-3.shtml
   iData_type        = 1#forecast product > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table1-4.shtml
   ident_sect[11:]   = [iProd_status,iData_type]

   return ident_sect
####################################################################################################

####################################################################################################
def def_gdsinfo(ncinfo,gdtnum):
   # define grid definition section:
   # >http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_sect3.shtml 
   # inputs:
   # 1) ncinfo: object defined by mod_reading.py from netcdf file
   # 2) gdtnum: grid template number (from def_grid_template)
   
   i_Source_grid     = 0 # source of grid definition > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-0.shtml
   N_xtra_gps        = 0 # Number of octets needed for each additional grid points defn.
   N_interpret_gps   = 0 # interpretation of additional points

   gdsinfo  = [i_Source_grid,ncinfo.Npts,N_xtra_gps,N_interpret_gps,gdtnum]

   return np.array(gdsinfo,'i')
####################################################################################################

####################################################################################################
def def_grid_template(ncinfo):
   # grid definition template:
   # link from > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-1.shtml

   proj_in     = ncinfo.proj_info
   proj_name   = proj_in.grid_mapping_name

   if proj_name=='polar_stereographic':
      # stereographic > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
      gdtnum   = 20

      nx    = ncinfo.Npts_x
      ny    = ncinfo.Npts_y
      gdt   = []

      # grid parameters for polar stereographic projection
      # > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
      i_earth_shape  = 1   #spherical (specified radius) > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-2.shtml
      gdt.append(i_earth_shape)

      # parameters for spherical earth
      re_h2c         = 6378273 # radius of earth (in m) from Hyc2proj/mod_toproj.F90 (subroutine polar_stereographic)
      re_scale_fac   = 0       # radius = scaled_radius*10^(-scale_fac)
      gdt.extend([re_scale_fac,re_h2c])

      #not applicable (oblate spherical)
      gdt.extend([0,0,0,0])

      # no of points in each dirn
      gdt.extend([nx,ny])

      # lon,lat of 1st grid point:
      fac   = 1e6 # degrees -> micro-degrees
      lon0  = int(ncinfo.lon0*fac)
      lat0  = int(ncinfo.lat0*fac)
      if lon0<0:
         lon0  = lon0+360*fac
      if lon0<0:
         lat0  = lat0+360*fac
      gdt.extend([lat0,lon0])

      # resolution and component flags > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-3.shtml
      # - most not applicable (only rcf[4] - definition of reference for vector directions)
      if proj_in.false_easting==1.0:
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
      latc  = proj_in.latitude_of_projection_origin  # lat where resolution is specified
      lonc  = proj_in.longitude_of_projection_origin # orientation = meridian which is parallel to the y-axis
      latc  = int(latc*fac)
      lonc  = int(lonc*fac)
      if lonc<0:
         lonc  = lonc+fac*360
      gdt.extend([latc,lonc])

      # resolution
      fac   = 1e3 # m -> mm
      dx    = proj_in.x_resolution  # x dirn
      dy    = proj_in.y_resolution  # y dirn
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

   return gdtnum,np.array(gdt,'i')
####################################################################################################

####################################################################################################
def def_prod_template(ncinfo,data_name):
   pdtmpl   = [] # product definition template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp4-0.shtml

   # Parameter category & Parameter number
   if data_name=='fice':
      param_cat   = 2 # Discipline_code = 10: 2 -> Ice > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
      param_num   = 0   # ice cover > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
   elif data_name=='hice':
      param_cat   = 2 # Discipline_code = 10: 2 -> Ice > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
      param_num   = 1   # ice thickness > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
   elif data_name=='swh':
      param_cat   = 0 # Discipline_code = 10: 0 -> Waves > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
      param_num   = 4 # total significant wave height > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-10-2.shtml
   elif data_name=='mwp':
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
   if ncinfo.timeunits=='hour':
      ft_unit  = 1 # hours
   pdtmpl.append(ft_unit)

   # forecast time (time since "reftime" = start of forecast )
   ft = ncinfo.datatime
   pdtmpl.append(ft)

   # Type of fixed surface  > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-5.shtml
   tofs  = 1   # ground or water surface
   pdtmpl.append(tofs)

   fs_scale_fac   = 0
   fs_scale_val   = 1
   pdtmpl.extend([fs_scale_fac,fs_scale_val])

   # 2nd fixed surface: missing
   pdtmpl.extend([255,0,0])

   # convert to numpy array
   return np.array(pdtmpl,'i')
####################################################################################################

####################################################################################################
def def_datarep_template(drtnum):

   if drtnum==40:
      # jpeg compression
      # > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp5-40.shtml
      # drtmpl  = [1129381888,0,2,11,0,0,255] # drtnum = 40
      drtmpl      = 6*[0] # drtnum = 40
      drtmpl[0]   = 0 # lossless compression
      drtmpl[5]   = 255 # missing

   else:
      # simple packing
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

   return np.array(drtmpl,'i')
####################################################################################################
