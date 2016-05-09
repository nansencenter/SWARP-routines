# python $pydir/analyse_DAILY.py --infile=<<DAILY average .a file>> --outdir=<<path to output dir>> --plotting=<<True/False>>
import os,sys,matplotlib
matplotlib.use('Agg')
import mod_reading   as mr
import MIZchar       as mc
from getopt import getopt

if 0:
   # plot regions ie Greenland Sea, Labrador Sea,...
   mc.plot_regions_v2()


# ==========================================================================
# options
infile   = None
outdir   = '.'
plotting = False

opts,args   = getopt(sys.argv[1:],"",["infile=","outdir=","plotting="])
for opt,arg in opts:
   if opt=='--infile':
      infile   = arg
   if opt=='--outdir':
      outdir   = arg
   if opt=='--plotting':
      from distutils.util import strtobool
      plotting = strtobool(arg)
      # if arg not in ['True','False']:
      #    raise ValueError('Use --plotting=True or --plotting=False only')
      #    plotting = (arg=='True')

if infile is None:
   raise ValueError('specify infile with --infile="..."')

hi = mr.HYCOM_binary_info(infile)
if not os.path.exists(outdir):
   os.mkdir(outdir)
# ==========================================================================
basename  = infile.split('/')[-1]


# ==========================================================================
# plot OSISAF ice_edge on fice:
if 1:
   odir     = outdir+'/../../IceEdgeOSISAF' # in analysis
   if not os.path.exists(odir):
      os.mkdir(odir)
   figname  = odir+'/'+basename[:-2]+'_vIceEdgeOSISAF.png'
   hi.compare_ice_edge_obs(obs_type='OSISAF',date_label=1,figname=figname)
# ==========================================================================



# ==========================================================================
if 1:
   odir  = outdir+'/vOSISAF'
   if not os.path.exists(odir):
      os.mkdir(odir)

   # Distance to Ice Edge
   hi.areas_of_disagreement(obs_type='OSISAF',outdir=odir,do_sort=True,plotting=plotting)
# ==========================================================================

if 0:
   # ==========================================================================
   # MIZ widths (hice,fice)
   for vname in ['fice','hice']:
      # Get MIZ width for all variables
      odir  = outdir+'/'+vname
      if not os.path.exists(odir):
         os.mkdir(odir)
      hi.MIZmap(var_name=vname,do_sort=True,outdir=odir,plotting=plotting)
   # ==========================================================================
