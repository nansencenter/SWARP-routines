import os,sys
import numpy as np

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

def def_gdsinfo(proj_in,Npts):
   # define grid definition section:
   # >http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_sect3.shtml 
   # inputs:
   # 1) proj_in: object defined by mod_reading.py from proj.in file
   #              (used by hyc2proj)
   # 2) Npts: total number of points in the grid 
   
   i_Source_grid     = 0 # source of grid definition > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-0.shtml
   N_xtra_gps        = 0 # Number of octets needed for each additional grid points defn.
   N_interpret_gps   = 0 # interpretation of additional points

   if proj_in.projection_name=='polar_stereographic':
      i_grid_templ_no   = 20  # polar stereographic, grid template no
                              # > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-1.shtml

   gdsinfo  = [i_Source_grid,Npts,N_xtra_gps,N_interpret_gps,i_grid_templ_no]

   return np.array(gdsinfo,'i')
