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

# scalars
# TODO label colorbars
clabs = {'icec':'Ice concentration','icetk':'Ice thickness'}
for vname in clabs.keys():
   figdir   = outdir+'/'+vname
   if not os.path.exists(figdir):
      os.mkdir(figdir)
   nci.make_png_all(vname,ice_mask=True,figdir=figdir,clabel=clabs[vname])


# vectors
# TODO put directions on them
# TODO label colorbars
vecs     = {'uice':'ice_speed','usurf':'surf_speed'}
masks    = {'uice':True,'usurf':False}
clabs    = {'uice':'Ice speed, km/day','usurf':'Surface speed, km/day'}
conv_fac = 24*3600/1.e3 # m/s -> km/day

for vname in vecs.keys():
   figdir   = outdir+'/'+vecs[vname]
   if not os.path.exists(figdir):
      os.mkdir(figdir)
   nci.make_png_all(vname,figdir=figdir,\
      conv_fac=conv_fac,vec_mag=True,ice_mask=masks[vname],\
         clabel=clabs[vname])
