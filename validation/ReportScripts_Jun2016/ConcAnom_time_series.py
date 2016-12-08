import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import datetime as dtm

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
refdate  = dtm.datetime(2015,3,1)
xticks   = []
xtlabs   = []
tfac     = 1./(24*3600) # seconds to days

years = range(2015,2017)
# years =[2015]
for year in years:
   for M in range(2,13,2):
      dto   = dtm.datetime(year,M,1)
      xtlabs.append(dto.strftime('%Y%m%d'))
      xticks.append((dto-refdate).total_seconds()*tfac)

labs     = ['RMSE','Bias']
cols     = {'RMSE':'r','Bias':'b'}

ODIR  = 'out/'
if not os.path.exists(ODIR):
   os.mkdir(ODIR)

PO       = MR.plot_object()
figname  = ODIR+'/ConcAnom'+'.png'
txtname  = ODIR+'/ConcAnom'+'.txt'

# write header to text file
tid   = open(txtname,'w')
blk   = 3*' '
ss    = 'Date'+blk+\
         'RMSE'+blk+'Bias\n'
tid.write(ss)

Tmin  = 1e30
Tmax  = -1e30
for year in years:
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
   if 0:
      print('\nResults dirs ('+cyear+'):')
      print(wlist)
      print('\n')
   # =============================================================

   for wdir in wlist:
      # loop over weeks
      its  += 1
      fid   = open(stat,'w')
      fid.write(wdir+'\n')
      ss = str(its)+' | '+str(len(wlist))
      fid.write(ss)
      fid.close()

      # place to get results
      outdir   = mdir2+'/'+wdir+'/analysis/AODs/time_series/'
      if not os.path.exists(outdir):
         continue

      olist    = os.listdir(outdir)
      for ofil in olist:
         if ('conc_anomaly' in ofil) and ('.txt' in ofil):
            TS  = MR.read_time_series(outdir+ofil)
            continue
            
      rmse  =  TS.data['RMSE_either_ice']
      bias  =  TS.data['Bias_either_ice']

      tdays = []
      for i,date in enumerate(TS.dates):
         dt = date-refdate
         ti = dt.total_seconds()*tfac
         tdays.append(ti)
         Tmin  = min(Tmin,ti)
         Tmax  = max(Tmax,ti)


         # =================================================================
         # write to file:
         ss = '%f   %f\n' %(rmse[i],bias[i])
         ss = dto.strftime('%Y%m%dT%H%M%SZ')+blk+ss
         tid.write(ss)
         # =================================================================
                    
      
      lines = []
      lin  ,= PO.ax.plot(tdays,rmse,color=cols['RMSE'])
      lines.append(lin)
      lin  ,= PO.ax.plot(tdays,bias,color=cols['Bias'])
      lines.append(lin)

# finished loop over years
PO.ax.legend(lines,labs,loc=4)
tid.close()
print('\nFinished writing to '+txtname)

Xticks   = []
Xtlabs   = []
for i,xt in enumerate(xticks):
   if xt>=Tmin and xt<=Tmax:
      Xticks.append(xt)
      Xtlabs.append(xtlabs[i])

PO.ax.set_xticks(Xticks)
PO.ax.set_xticklabels(Xtlabs,rotation=270)
PO.ax.set_ylabel('Concentration Error')
PO.fig.savefig(figname,bbox_inches='tight')
print('\nSaving '+figname)
PO.ax.cla()
plt.close(PO.fig)
