import os,sys
import shapefile_utils  as SFU
from getopt import getopt
import matplotlib
matplotlib.use('Agg')

indir          = None
outdir         = None
MIZ_criteria   = "FA_only"
chart_source   = "DMI"

opts,args   = getopt(sys.argv[1:],"",["MIZ_criteria=","indir=","outdir=","chart_source="])
for opt,arg in opts:
   if opt=='--indir':
      indir = arg
   if opt=='--outdir':
      outdir   = arg
   if opt=='--MIZ_criteria':
      MIZ_criteria   = arg
   if opt=='--chart_source':
      chart_source   = arg

if indir is None:
   raise ValueError('Specify input dir with --indir=')
if outdir is None:
   raise ValueError('Specify output dir with --outdir=')



############################################################################
fnames   = os.listdir(indir)
snames   = []

if not os.path.exists(outdir):
   os.mkdir(outdir)

cdates   = []
for fname in fnames:
   ext   = os.path.splitext(fname)
   if ext[1]=='.shp':
      snames.append(ext[0])
      if chart_source=='DMI':
         cdate = ext[0][:8]
      else:
         raise ValueError('unknown chart type')
      cdates.append(int(cdate))


# sort in order of increasing date
dummy    = sorted([(e,i) for i,e in enumerate(cdates)])
snames_  = [snames[i] for e,i in dummy]
snames   = 1*snames_
del(cdates)
del(snames_)
del(dummy)

if not os.path.exists(outdir):
   os.mkdir(outdir)
figdir   = outdir+'/png'
if not os.path.exists(figdir):
   os.mkdir(figdir)
txtdir   = outdir+'/MIZ_polys_classified'
if not os.path.exists(txtdir):
   os.mkdir(txtdir)


for fname in snames:
   fname_full  = indir+"/"+fname+'.shp'
   print('\nOpening '+fname_full+'...')
   
   MIZshp   = SFU.MIZ_from_shapefile(fname_full,MIZ_criteria=MIZ_criteria)

   # ==================================================================
   #  make a figure
   cdate = fname[:8]
   Fname = figdir+'/'+fname
   MIZshp.test_plot(figname=Fname+'.png')
   # ==================================================================


   # ============================================================
   # write text files
   Fname = txtdir+'/'+fname
   # tfil1 = Fname+'_MIZpolys_all.txt'
   # MIZshp.write_text_file(tfil1)

   tfil1 = Fname+'_MIZpolys.txt'
   MIZshp.write_text_file(tfil1,MIZtype='outer')
   # ============================================================

   # sys.exit()
