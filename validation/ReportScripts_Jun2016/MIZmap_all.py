import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')

# =============================================================
if 0:
   # on compute node, can't use files from /home
   # - need to copy to /work NB should be done before going to compute node
   py_funs  = os.getenv('SWARP_ROUTINES')+'/py_funs'
   sys.path.remove(py_funs) # remove dir from PYTHONPATH
   py_funs2 = 'py_funs'
   sys.path.append(py_funs2)

print('\nPaths:')
for d in sys.path:
   if 'py_funs' in d:
      print(d)
print('\n')
# =============================================================


import mod_reading as MR
import fns_plotting as FP

mdir     = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice_erai2/'
subdir   = 'binaries/archv_wav/'

# write progress to text file
if not os.path.exists('log'):
   os.mkdir('log')
stat  = 'log/status.txt'
its   = 0

plotting = True

for year in range(2015,2017):
   cyear = str(year)
   mdir2 = mdir+cyear
   wlist = os.listdir(mdir2)

   # =============================================================
   # eliminate non-result dirs's and sort into numerical order
   Wlist = []
   Ilist = []
   for wdir in wlist:
      if cyear+'_' in wdir:
         Wlist.append(wdir)
         Ilist.append(int(wdir[-3:]))
   ii    = sorted([(I,i) for i,I in enumerate(Ilist)])
   wlist = [Wlist[cc[1]] for cc in ii]
   del(Wlist,Ilist,ii)
   print('\nResults dirs ('+cyear+'):')
   print(wlist)
   print('\n')
   # sys.exit()
   # =============================================================

   for wdir in wlist:

      its  += 1
      fid   = open(stat,'w')
      fid.write(wdir+'\n')
      ss = str(its)+' | '+str(len(wlist))
      fid.write(ss)
      fid.close()

      # place to save results
      outdir   = mdir2+'/'+wdir+'/analysis'
      if not os.path.exists(outdir):
         os.mkdir(outdir)

      outdir  += '/MIZmap'
      if not os.path.exists(outdir):
         os.mkdir(outdir)
      else:
         continue

      flist = MR.file_list(mdir2+'/'+wdir+'/'+subdir,'archv_wav','.a')
      step  = 0.5
      start_date  = flist.datetimes[2].strftime('%Y%m%dT%H%M%SZ')
      end_date    = flist.datetimes[-1].strftime('%Y%m%dT%H%M%SZ')
      # end_date    = flist.datetimes[4].strftime('%Y%m%dT%H%M%SZ')

      regions  = ['gre','bar','lab','beau']
      out      = flist.MIZmap_all(end_date=end_date,start_date=start_date,step=step,\
            outdir=outdir,regions=regions,plotting=plotting)

      if 1:
         # do 1 week at a time (can do all in parallel)
         print('\nExit python\n')
         sys.exit()

print('\nExit python\n')
