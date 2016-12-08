import os,sys
import numpy as np
import mod_reading as mr

if len(sys.argv)>1:
   indir = os.path.abspath(sys.argv[1])
else:
   indir = os.getcwd()

dirs  = os.listdir(indir)
tdir  = indir+'/time_series'
if not os.path.exists(tdir):
   os.mkdir(tdir)


# ================================================================================
# 1st loop to see how many differnet time series eg regions are present
Files    = {}
regions  = []
bad      = []
for ddir in dirs:
   if ddir=='time_series':
      continue

   Ddir     = indir+'/'+ddir
   files    = os.listdir(Ddir)

   dfiles   = {}
   for f in files:
      if 'summary.txt' in f:
         reg   = f.split('_summary.txt')[0]
         reg   = reg.split(ddir)[-1]
         if reg=='':
            reg   = 'all'
         else:
            reg   = reg[1:]
         regions.append(reg)
         regions  = list(set(regions))
         dfiles.update({reg:f})

   if len(dfiles)==0:
      print('\n'+60*'--')
      print('no summary.txt files in '+Ddir) 
      print(60*'--'+'\n')
   else:
      Files.update({Ddir:dfiles})
# ================================================================================

Ndates   = len(Files)
Nseries  = len(regions)

Data  = {}
for reg in regions:
   Data.update({reg:{}})

for Ddir in Files:
   dfiles   = Files[Ddir]
   dfiles2  = {}
   for reg in regions:
      dfiles2.update({reg:[]})
      if reg in dfiles:
         dfiles2[reg]   = dfiles[reg]
   Files[Ddir] = dfiles2

L0    = Ndates*[np.nan]
Init  = True
Dates = 1*L0
for n,Ddir in enumerate(Files):
   for reg in regions:

      f0 = Files[Ddir][reg]
      if len(f0)>0:
         # read the file
         f  = Ddir+'/'+f0
         f_info   = mr.read_MIZpoly_summary(f)
         dto      = f_info.info['datetime']
         fields   = f_info.info['fields']
         # print(f)

         # assign Data
         if Init:
            for fld in fields:
               Data[reg].update({fld:1*L0})
               Init  = False

         for fld in fields:
            Data[reg][fld][n] = getattr(f_info,fld)

   # load the date;
   Dates[n] = dto

Tdir  = indir+'/time_series'
if not os.path.exists(Tdir):
   os.mkdir(Tdir)

TS = {}
for reg in regions:
   ts = mr.time_series(Dates,Data[reg],filename=Tdir+'/time_series_'+reg+'.txt',overwrite=True)
   TS.update({reg:ts})
