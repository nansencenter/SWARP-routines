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
   clabs = {'icec':'Ice concentration','icetk':'Ice thickness'}
   clims = {'icec':[0,1],'icetk':[0,4]}
   for vname in clabs.keys():
      figdir   = outdir+'/'+vname
      if not os.path.exists(figdir):
         os.mkdir(figdir)
      nci.make_png_all(vname,ice_mask=True,figdir=figdir,\
            clabel=clabs[vname],clim=clims[vname])

if 1:
   # vectors
   if 0:
      vec_opt  = 1   # plots speed
      vecs     = {'uice':'ice_speed','usurf':'surf_speed'}
   else:
      vec_opt  = 2   # plots speed with arrows showing direction
      vecs     = {'uice':'ice_vel','usurf':'surf_vel'}

   masks    = {'uice':True,'usurf':False}
   clabs    = {'uice':'Ice speed, km/day','usurf':'Surface speed, km/day'}
   clims    = {'uice':[0,40],'usurf':[0,40]}
   conv_fac = 24*3600/1.e3 # m/s -> km/day

   for vname in vecs.keys():
      figdir   = outdir+'/'+vecs[vname]
      if not os.path.exists(figdir):
         os.mkdir(figdir)
      nci.make_png_all(vname,figdir=figdir,\
         conv_fac=conv_fac,vec_opt=vec_opt,ice_mask=masks[vname],\
            clabel=clabs[vname],clim=clims[vname])
