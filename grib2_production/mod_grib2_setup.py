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
