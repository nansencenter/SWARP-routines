import os,sys

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import mod_grib2_setup as m_g2s

if 1:
   fildir      = 'test_ncfiles'
   # file_list   = ['SWARPwavesice_forecast_start20150723T000000Z.nc']
   file_list   = ['SWARPiceonly_forecast_start20150723T000000Z.nc']
else:
   fildir      = '../netcdf_production/output'
   file_list   = [ f for f in os.listdir(fildir)
                     if os.path.isfile(os.path.join(fildir,f)) ]

#######################################################################
for ncfil in file_list:

   # set name of output file:
   outfil   = ncfil.replace('.nc','.grb2')
   outdir  = 'out'
   if not os.path.exists(outdir):
      os.makedirs(outdir)
   fil_out = outdir+'/'+outfil

   # do conversion
   ncinfo   = m_g2s.hyc2proj_to_grib2(fildir+'/'+ncfil,fil_out)
#######################################################################
