import os,sys
from getopt import getopt
import mod_reading as Mr
import matplotlib
matplotlib.use('Agg')

ncfile   = None
outdir   = '.'
opts,args   = getopt(sys.argv[1:],"",["ncfile=","outdir="])
for opt,arg in opts:
   if opt=='--ncfile':
      ncfile   = arg
   if opt=='--outdir':
      outdir   = arg

if ncfile is None:
   raise ValueError('no netcdf file specified (use --ncfile==)')

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

if 1:
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
   clims.update({vname:[0,5]})
   scalars.append(Mr.make_plot_options(vname,ice_mask=True))

   for V in scalars:
      vname    = V.name
      figdir   = outdir+'/'+vname
      if not os.path.exists(figdir):
         os.mkdir(figdir)
      nci.make_png_all(V,figdir=figdir,\
            clabel=clabs[vname],clim=clims[vname])

if 1:
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
      nci.make_png_all(V,figdir=figdir,\
            clabel=clabs[vname],clim=clims[vname])
