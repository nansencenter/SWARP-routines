import os,sys
from getopt import getopt
import mod_reading as Mr
import matplotlib
matplotlib.use('Agg')

# ==========================================================================
# options
indir    = None
outdir   = '.'
FCtype   = 'ice_only'

opts,args   = getopt(sys.argv[1:],"",["indir=","outdir="])
for opt,arg in opts:
   if opt=='--indir':
      indir = arg
   if opt=='--outdir':
      outdir   = arg

if indir  is None:
   raise ValueError('no input directory specified (use --indir==)')
# ==========================================================================

print(' ')
print('Using directory: '+indir)
print('Saving figures to: '+outdir)
print(' ')

nci   = Mr.file_list(indir,'DAILY','.a')
print('Start of forecast:')
print(nci.datetimes[0])
print('\nEnd of forecast:')
print(nci.datetimes[-1])
print('\n')

HYCOMreg = nci.HYCOM_region

if FCtype=="ice_only":
   # =============================================================
   # ice-only var's

   # scalars
   scalars  = []
   clabs = {}
   clims = {}

   Figdir   = outdir+'/ice_edge_OSISAF'
   if not os.path.exists(Figdir):
      os.mkdir(Figdir)
   nci.compare_ice_edge_OSISAF_all(figdir=Figdir)

   vo0      = Mr.make_plot_options('hice',ice_mask=True)
   Figdir   = outdir+'/hice'
   if not os.path.exists(Figdir):
      os.mkdir(Figdir)
   nci.make_png_all(vo0,figdir=Figdir,date_label=1,clabel='Thickness, m',clim=[0,5])

   vo1      = Mr.make_plot_options('temp',layer=1)
   Figdir   = outdir+'/sst'
   if not os.path.exists(Figdir):
      os.mkdir(Figdir)
   nci.make_png_all(vo1,figdir=Figdir,date_label=1,clabel='SST, $^\circ$C',clim=[-2,6])
   # sys.exit('HEY!')

   Figdir   = outdir+'/usurf'
   if not os.path.exists(Figdir):
      os.mkdir(Figdir)
   vo2   = Mr.make_plot_options('utot',layer=1,vec_opt=1,conv_fac=24*3600/1.e3)
   vo3   = Mr.make_plot_options('utot',layer=1,vec_opt=2,conv_fac=24*3600/1.e3)
   nci.make_png_pair_all(vo2,vo3,figdir=Figdir,date_label=1,\
         clabel='Surface speed, km/day',clim=[0,40])
   # =============================================================
