import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')

# =============================================================
if 1:
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
subdir   = 'binaries/DAILY/'

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
      else:
         continue

      outdir  += '/AODs'
      if not os.path.exists(outdir):
         os.mkdir(outdir)

      flist = MR.file_list(mdir2+'/'+wdir+'/'+subdir,'DAILY','.a')

      regions  = ['gre','bar','lab','beau']
      obs_type = 'OSISAF'
      out      = flist.AODs_all(outdir=outdir,regions=regions,obs_type=obs_type,plotting=plotting)

      # ============================================
      # weekly average of conc anomalies
      cdates   = []
      nlist    = []
      wt       = 1./7
      anom_av  = 0.

      # time series
      data  = {'Bias_both_ice':[],\
               'Bias_either_ice':[],\
               'RMSE_both_ice':[],\
               'RMSE_either_ice':[]}
      keys  = data.keys()

      # cdates   = []
      # for day in range(17,24):
      #    cdates.append('201508'+str(day))
      # for i,cdate in enumerate(cdates):

      for i,dto in enumerate(out.times_to_analyse):
         cdate    = dto.strftime('%Y%m%d')
         nfil     = outdir+'/'+cdate+'/conc_anomaly_'+obs_type+'_'+cdate+'.npz'
         fid      = np.load(nfil)
         anom_av += wt*fid['anomaly']
         if i==0:
            cdate0   = cdate
            mask     = fid['mask']
            lon      = fid['lon']
            lat      = fid['lat']
         fid.close()

         tfil  = nfil.replace('.npz','.txt')
         info  = MR.read_MIZpoly_summary(tfil)
         for key in keys:
            data[key].append(getattr(info,key))

      anom_av  = np.ma.array(anom_av,mask=mask)
      cmax     = .5
      anom_fig = outdir+'/time_series/conc_anomaly_weekly_av_'+cdate+'.png'
      FP.plot_anomaly(lon,lat,anom_av,anom_fig,\
            text=cdate,\
            clim=[-cmax,cmax],\
            clabel='Concentration anomaly')

      anom_ts  = MR.time_series(out.times_to_analyse,data,\
                     filename=outdir+'/time_series/conc_anomaly_weekly_'+cdate+'.txt')
      # ============================================

      # sys.exit()

print('\nExit python\n')
