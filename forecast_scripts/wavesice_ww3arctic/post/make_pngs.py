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

   vname = 'dmax'
   clabs.update({vname:'Maximum floe size, m'})
   clims.update({vname:[0,300]})
   scalars.append(Mr.make_plot_options(vname,ice_mask=True))

   for V in scalars:
      vname    = V.name
      figdir   = outdir+'/'+vname
      if not os.path.exists(figdir):
         os.mkdir(figdir)
      nci.make_png_all(V,figdir=figdir,\
            clabel=clabs[vname],clim=clims[vname])

if 1:
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
      nci.make_png_pair_all(v1,v2,figdir=figdir,\
            clabel=clabs[vname],clim=clims[vname])
