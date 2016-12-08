import os,sys
import matplotlib
matplotlib.use('Agg')
import shapefile_utils  as SFU
from getopt import getopt
import datetime as dtm
import mod_reading as mr

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
      if arg.lower() in ['true','t','yes','y','1']:
         overwrite   = True
      elif arg.lower() not in ['false','f','no','n','0']:
         raise ValueError('Unknown value for overwrite: '+arg+'- use eg T or F')

if indir is None:
   raise ValueError('Specify input dir with --indir=')
if outdir is None:
   raise ValueError('Specify output dir with --outdir=')



############################################################################
fnames      = os.listdir(indir)
snames      = []
datetimes   = []
#print(fnames)

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

   Txtdir   = txtdir+'/'+cyear+'/'+fname
   if not os.path.exists(Txtdir):
      os.mkdir(Txtdir)
   txtname  = Txtdir+'/'+fname+'_MIZpolys.txt'

   if os.path.exists(figname) and (not overwrite):
      print('\n'+figname+' exists')
      print('-skipping '+fname_full+'...')
      continue

   print('\nOpening '+fname_full+'...')
   
   if 1:
      MIZshp   = SFU.MIZ_from_shapefile(fname_full,\
            MIZ_criteria=MIZ_criteria,chart_source=chart_source)

      # ============================================================
      # write text file
      MIZshp.write_text_file(txtname,MIZtype='outer')
      # ============================================================


      # ==================================================================
      #  make a figure
      MIZshp.test_plot(figname=figname)
      # ==================================================================


   # ==================================================================
   # calc MIZ width
   chart_list  = mr.polygon_file_list(file_info={'files':[txtname],'dates':[dto]})

   verts = None
   if 1:
      # Greenland Sea as defined by Dany/Strong
      verts = [(-44,50),(-44,89),(8,89),(8,50)]
      verts.append(verts[0])
   mil   = chart_list.get_solutions(time_index=0,vertices=verts)
   # ==================================================================


   # ==================================================================
   # save summary file
   ofil  = Txtdir+'/'+fname+'_summary.txt'
   out   = mil.save_summary(ofil,date_time=dto)

   # save shapefile
   ofil  = Txtdir+'/'+fname+'.shp'
   out   = mil.save_shapefile(filename=ofil)
   # ==================================================================
