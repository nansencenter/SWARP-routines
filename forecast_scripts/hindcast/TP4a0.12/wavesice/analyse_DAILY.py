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

hi       = mr.HYCOM_binary_info(infile)
if not os.path.exists(outdir):
   os.mkdir(outdir)
# ==========================================================================


# ==========================================================================
# Distance to Ice Edge
odir  = outdir+'/vOSISAF'
if not os.path.exists(odir):
   os.mkdir(odir)
hi.areas_of_disagreement(obs_type='OSISAF',outdir=odir,do_sort=True,plotting=plotting)
# ==========================================================================


# ==========================================================================
# MIZ widths (hice,fice)
for vname in ['fice','hice']:
   # Get MIZ width for all variables
   odir  = outdir+'/'+vname
   if not os.path.exists(odir):
      os.mkdir(odir)
   hi.MIZmap(var_name=vname,do_sort=True,outdir=odir,plotting=plotting)
# ==========================================================================
