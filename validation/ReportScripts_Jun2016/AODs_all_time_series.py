import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import datetime as dtm

import mod_reading as MR
import fns_plotting as FP

mdir     = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice_erai2/'

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

regions  = ['gre','bar','beau','lab']
labs     = ['Over','Under','Bias']
cols     = {'Over':'r','Under':'b','Bias':'g'}

ODIR  = 'out/'
if not os.path.exists(ODIR):
   os.mkdir(ODIR)

for reg in regions:
   PO       = MR.plot_object()
   figname  = ODIR+'/AOD_'+reg+'.png'
   txtname  = ODIR+'/AOD_'+reg+'.txt'

   # write header to text file
   tid   = open(txtname,'w')
   blk   = 3*' '
   ss    = 'Date'+blk+\
            'Total_area_over'+blk+'Width_over'+blk+\
            'Total_area_under'+blk+'Width_under'+blk+\
            'Bias\n'
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
            if (reg in ofil) and ('Over' in ofil):
               O_ts  = MR.read_time_series(outdir+ofil)
               continue
            if (reg in ofil) and ('Under' in ofil):
               U_ts  = MR.read_time_series(outdir+ofil)
               continue
               
         A_over   =  O_ts.data['Total_area']
         W_over   =  O_ts.data['Mean_intersecting_width']
         A_under  =  U_ts.data['Total_area']
         W_under  = -U_ts.data['Mean_intersecting_width']

         # change some nans to 0
         for i in range(len(W_under)):
            Wo = W_over[i]
            Wu = W_under[i]
            if np.isnan(Wo) and (not np.isnan(Wu)):
               W_over[i]   = 0.
               A_over[i]   = 0.
            if np.isnan(Wu) and (not np.isnan(Wo)):
               W_under[i]   = 0.
               A_under[i]   = 0.
         Bias  =  (A_over*W_over+A_under*W_under)/(A_over+A_under)
         
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
            ss = '%f   %f   %f   %f   %f\n'\
                     %(A_over[i],W_over[i],A_under[i],W_under[i],Bias[i])
            ss = dto.strftime('%Y%m%dT%H%M%SZ')+blk+ss
            tid.write(ss)
            # =================================================================
                       
         
         lines = []
         lin  ,= PO.ax.plot(tdays,W_over/1.e3,color=cols['Over'])
         lines.append(lin)
         lin  ,= PO.ax.plot(tdays,W_under/1.e3,color=cols['Under'])
         lines.append(lin)
         lin  ,= PO.ax.plot(tdays,Bias/1.e3,color=cols['Bias'])
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
   PO.ax.set_ylabel('Ice edge disagreement, km')
   PO.fig.savefig(figname,bbox_inches='tight')
   print('\nSaving '+figname)
   PO.ax.cla()
   plt.close(PO.fig)
