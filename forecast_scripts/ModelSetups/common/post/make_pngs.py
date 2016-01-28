import os,sys
from getopt import getopt
import mod_reading as Mr
import matplotlib
matplotlib.use('Agg')

# ==========================================================================
# options
ncfile   = None
outdir   = '.'
FCtype   = None
scalars  = 1
pairs    = 1
vectors  = 1

opts,args   = getopt(sys.argv[1:],"",["ncfile=","outdir=","FCtype="])
for opt,arg in opts:
   if opt=='--ncfile':
      ncfile   = arg
   if opt=='--outdir':
      outdir   = arg
   if opt=='--FCtype':
      FCtype   = arg

if ncfile is None:
   raise ValueError('no netcdf file specified (use --ncfile==)')

if FCtype is None:
   raise ValueError('forecast type not specified (use --FCtype==)')
# ==========================================================================

print(' ')
print('Using netcdf file: '+ncfile)
print('Saving figures to: '+outdir)
print(' ')

nci   = Mr.nc_getinfo(ncfile)
print('Start of forecast:')
print(nci.datetimes[0])
print('\nEnd of forecast:')
print(nci.datetimes[-1])
print('\n')

if 'area' in nci.ncattrs.list():
   HYCOMreg = nci.ncattrs.area[:3]
elif 'area_name' in nci.ncattrs.list():
   HYCOMreg = nci.ncattrs.area_name[:3]
else:
   HYCOMreg = 'TP4'

if FCtype!="ice_only":
   # =============================================================
   # wave-ice FC's

   if scalars:
      # scalars
      scalars  = []
      clabs = {}
      clims = {}

      vname = 'dmax'
      clabs.update({vname:'Maximum floe size, m'})
      clims.update({vname:[0,300]})
      scalars.append(Mr.make_plot_options(vname,ice_mask=True))

      for V in scalars:
         vname    = V.name
         figdir   = outdir+'/'+vname
         if not os.path.exists(figdir):
            os.mkdir(figdir)
         nci.make_png_all(V,figdir=figdir,HYCOMreg=HYCOMreg,\
               clabel=clabs[vname],clim=clims[vname])

   if pairs:
      # pairs
      pairs = []
      clabs = {}
      clims = {}

      vname = 'swh'
      clabs.update({vname:'$H_s$, m'})
      clims.update({vname:[0,5]})
      v1    = Mr.make_plot_options(vname,wave_mask=True)
      #
      vname = 'mwd'
      v2    = Mr.make_plot_options(vname,vec_opt=5,dir_from=True,wave_mask=True)
      pairs.append([v1,v2])


      for v1,v2 in pairs:
         vname    = v1.name
         figdir   = outdir+'/'+v1.name
         if not os.path.exists(figdir):
            os.mkdir(figdir)
         nci.make_png_pair_all(v1,v2,figdir=figdir,HYCOMreg=HYCOMreg,\
               clabel=clabs[vname],clim=clims[vname])

   # if vectors:
   if 0:
      # vectors
      vectors  = []
      clabs    = {}
      clims    = {}
      if 0:
         vec_opt  = 1   # plots speed
      else:
         vec_opt  = 2   # plots speed with arrows showing direction

      vname = 'taux_wav'
      clabs.update({vname:'Wave stress, Pa'})
      clims.update({vname:[0,.25]})
      vectors.append(Mr.make_plot_options(vname,vec_opt=vec_opt,\
         ice_mask=True, wave_mask=True ))

      for V in vectors:
         vname    = V.name
         figdir   = outdir+'/'+vname
         if not os.path.exists(figdir):
            os.mkdir(figdir)
         nci.make_png_all(V,figdir=figdir,HYCOMreg=HYCOMreg,\
               clabel=clabs[vname],clim=clims[vname])
   # =============================================================

else:

   # =============================================================
   # ice-only var's

   if scalars:
      # scalars
      scalars  = []
      clabs = {}
      clims = {}

      vname = 'icec'
      clabs.update({vname:'Ice concentration'})
      clims.update({vname:[0,1]})
      scalars.append(Mr.make_plot_options(vname,ice_mask=True))

      vname = 'icetk'
      clabs.update({vname:'Ice thickness, m'})
      if HYCOMreg=='FR1' or HYCOMreg=='BS1':
         clims.update({vname:[0,2]})
      else:
         clims.update({vname:[0,5]})
      scalars.append(Mr.make_plot_options(vname,ice_mask=True))

      vname = 'sst'
      clabs.update({vname:'SST, $^\circ$C'})
      clims.update({vname:[-2,6]})
      scalars.append(Mr.make_plot_options(vname))

      vname = 'ssh'
      clabs.update({vname:'SSH, m'})
      if HYCOMreg=='FR1' or HYCOMreg=='BS1':
         clims.update({vname:[-4.5,1.5]})
      else:
         clims.update({vname:[-7,5]})
      scalars.append(Mr.make_plot_options(vname))

      vname = 'sss'
      clabs.update({vname:'SSS, ppt'})
      clims.update({vname:[5,35]})
      scalars.append(Mr.make_plot_options(vname))

      for V in scalars:
         vname    = V.name
         figdir   = outdir+'/'+vname
         if not os.path.exists(figdir):
            os.mkdir(figdir)
         nci.make_png_all(V,figdir=figdir,HYCOMreg=HYCOMreg,\
               clabel=clabs[vname],clim=clims[vname])

   if vectors:
      # vectors
      vectors  = []
      clabs    = {}
      clims    = {}
      if 0:
         vec_opt  = 1   # plots speed
      else:
         vec_opt  = 2   # plots speed with arrows showing direction

      vname = 'uice'
      clabs.update({vname:'Ice speed, km/day'})
      clims.update({vname:[0,50]})
      vectors.append(Mr.make_plot_options(vname,vec_opt=vec_opt,\
         ice_mask=True, conv_fac = 24*3600/1.e3 ))# m/s -> km/day

      vname = 'usurf'
      clabs.update({vname:'Surface speed, km/day'})
      clims.update({vname:[0,50]})
      vectors.append(Mr.make_plot_options(vname,vec_opt=vec_opt,\
         conv_fac = 24*3600/1.e3 ))# m/s -> km/day

      for V in vectors:
         vname    = V.name
         figdir   = outdir+'/'+vname
         if not os.path.exists(figdir):
            os.mkdir(figdir)
         nci.make_png_all(V,figdir=figdir,HYCOMreg=HYCOMreg,\
               clabel=clabs[vname],clim=clims[vname])
   # =============================================================
