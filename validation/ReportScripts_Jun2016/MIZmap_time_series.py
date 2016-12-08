import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import datetime as dtm

import mod_reading as MR
import fns_plotting as FP

mdir     = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice_erai2/'
var_name = 'dmax'

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

years    = range(2015,2017)
events   = {}
events.update({'gre':\
      [dtm.datetime(2015,12,18,12,0)]})
events.update({'beau':None})
events.update({'bar':None})
events.update({'lab':None})

# years =[2015]
for year in years:
   for M in range(2,13,2):
      dto   = dtm.datetime(year,M,1)
      xtlabs.append(dto.strftime('%Y%m%d'))
      xticks.append((dto-refdate).total_seconds()*tfac)

regions  = ['gre','bar','beau','lab']
# for legend
#labs     = ['Over','Under','Bias']
#cols     = {'Over':'r','Under':'b','Bias':'g'}

ODIR  = 'out/'
if not os.path.exists(ODIR):
   os.mkdir(ODIR)

for reg in regions:
   PO       = MR.plot_object()
   figname  = ODIR+'/MIZwidth_'+var_name+'_'+reg+'.png'
   txtname  = ODIR+'/MIZwidth_'+var_name+'_'+reg+'.txt'

   # write header to text file
   tid   = open(txtname,'w')
   blk   = 3*' '
   ss    = 'Date'+blk+'Mean_intersecting_width\n'
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

         # place to get results
         outdir   = mdir2+'/'+wdir+'/analysis/MIZmap/time_series/'
         if not os.path.exists(outdir):
            continue

         olist    = os.listdir(outdir)
         for ofil in olist:
            if (reg in ofil):
               O_ts  = MR.read_time_series(outdir+ofil)
               continue

         W_miz =  O_ts.data['Mean_intersecting_width']
         W_miz[np.isnan(W_miz)]  = 0.  # should change some nans to 0?
         
         tdays = []
         for i,date in enumerate(O_ts.dates):
            dt = date-refdate
            ti = dt.total_seconds()*tfac
            tdays.append(ti)
            Tmin  = min(Tmin,ti)
            Tmax  = max(Tmax,ti)


            # =================================================================
            # write to file:
            # 'Total_area_over'+blk+'Width_over'+blk+\
            # 'Total_area_under'+blk+'Width_under'+blk+\
            # 'Bias\n'
            ss = '%f\n' %(W_miz[i])
            ss = dto.strftime('%Y%m%dT%H%M%SZ')+blk+ss
            tid.write(ss)
            # =================================================================
                       
         
         lines = []
         lin  ,= PO.ax.plot(tdays,W_miz/1.e3,color='b')
         lines.append(lin)

   # finished loop over years
   # PO.ax.legend(lines,labs,loc=4)
   if events[reg] is not None:
      ylim  = PO.ax.get_ylim()
      for ev in events[reg]:
         xev   = (ev-refdate).total_seconds()*tfac
         PO.ax.plot([xev,xev],ylim,'r',linewidth=1.5)
      PO.ax.set_ylim(ylim)

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
   PO.ax.set_ylabel('MIZ width, km')
   PO.fig.savefig(figname,bbox_inches='tight')
   print('\nSaving '+figname)
   PO.ax.cla()
   plt.close(PO.fig)
