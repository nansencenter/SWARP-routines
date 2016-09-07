import os,sys
import shapefile_utils  as SFU
from getopt import getopt
import matplotlib
import datetime as dtm
matplotlib.use('Agg')

indir          = None
outdir         = None
MIZ_criteria   = "FA_only"
chart_source   = "DMI"
overwrite      = False

opts,args   = getopt(sys.argv[1:],"",\
      ["MIZ_criteria=","indir=","outdir=","chart_source=",\
      "overwrite="])
for opt,arg in opts:
   if opt=='--indir':
      indir = arg
   if opt=='--outdir':
      outdir   = arg
   if opt=='--MIZ_criteria':
      MIZ_criteria   = arg
   if opt=='--chart_source':
      chart_source   = arg
   if opt=='--overwrite':
      overwrite   = arg

if indir is None:
   raise ValueError('Specify input dir with --indir=')
if outdir is None:
   raise ValueError('Specify output dir with --outdir=')



############################################################################
fnames      = os.listdir(indir)
snames      = []
datetimes   = []

if not os.path.exists(outdir):
   os.mkdir(outdir)

for fname in fnames:
   ext   = os.path.splitext(fname)
   if ext[1]=='.shp':
      snames.append(ext[0])
      if chart_source=='DMI':
         cdate = ext[0][:8]
      elif chart_source=='AARI':
         cdate = ext[0][9:17]
      else:
         raise ValueError('unknown chart type')
      datetimes.append(dtm.datetime.strptime(cdate,"%Y%m%d"))


# sort in order of increasing date
dummy       = sorted([(e,i) for i,e in enumerate(datetimes)])
snames      = [snames[i] for e,i in dummy]
datetimes   = [e for e,i in dummy]
del(dummy)

if not os.path.exists(outdir):
   os.mkdir(outdir)
figdir   = outdir+'/png'
if not os.path.exists(figdir):
   os.mkdir(figdir)
txtdir   = outdir+'/MIZ_polys_classified'
if not os.path.exists(txtdir):
   os.mkdir(txtdir)


for i,fname in enumerate(snames):
   # print(fname,datetimes[i])
   dto   = datetimes[i]
   cyear = dto.strftime("%Y")

   fname_full  = indir+'/'+fname+'.shp'
   if not os.path.exists(figdir+'/'+cyear):
      os.mkdir(figdir+'/'+cyear)
   figname  = figdir+'/'+cyear+'/'+fname+'.png'
   if not os.path.exists(txtdir+'/'+cyear):
      os.mkdir(txtdir+'/'+cyear)
   txtname  = txtdir+'/'+cyear+'/'+fname+'_MIZpolys.txt'

   if os.path.exists(figname) and os.path.exists(txtname) and (not overwrite):
      print('\nSkipping '+fname_full+'...')
      continue

   print('\nOpening '+fname_full+'...')
   
   MIZshp   = SFU.MIZ_from_shapefile(fname_full,\
         MIZ_criteria=MIZ_criteria,chart_source=chart_source)

   # ==================================================================
   #  make a figure
   MIZshp.test_plot(figname=figname)
   # ==================================================================


   # ============================================================
   # write text file
   MIZshp.write_text_file(txtname,MIZtype='outer')
   # ============================================================

   # sys.exit()
