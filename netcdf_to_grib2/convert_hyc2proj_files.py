import os
import mod_grib2_setup as m_g2s

# fildir      = 'test_ncfiles'
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
