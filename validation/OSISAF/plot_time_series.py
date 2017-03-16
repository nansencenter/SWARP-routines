# to launch from inside python:
# exec(open('plot_time_series.py').read(),globals())

import os,sys
import mod_reading as mr
from datetime import datetime as DT
import numpy as np

subdirs  = ['FC0days' , 'FC1days' , 'FC2days' , 'FC3days' , 'FC4days' , 'FC5days',\
                        'PFC1days', 'PFC2days', 'PFC3days', 'PFC4days', 'PFC5days']
Cfields  = ['RMSE_both_ice', 'Bias_both_ice', 'RMSE_either_ice', 'Bias_either_ice']
Ufields  = ['Maximum_total_width', 'Mean_intersecting_width', 'Mean_total_width',\
            'Maximum_intersecting_width', 'Total_perimeter', 'Total_area']

Ofields  = 1*Ufields
Bfields  = Ufields[1:3]
akey     = Ufields[-1]  # area


show  = 0
if 0:
   # reduce period plotted
   dto_lims = [DT(2017,1,1),DT(2017,3,12)]
   Mstep    = 1
else:
   # full period plotted
   dto_lims = [DT(2016,1,1),DT(2017,3,12)]
   Mstep    = 2

M0       = dto_lims[0].month
xlim     = []
refdate  = DT(2010,1,1)
tfac     = 1./(24*3600) # seconds to days
for dto in dto_lims:
   xlim.append((dto-refdate).total_seconds()*tfac)


for ng in range(2):
   if ng==0:
      # plot FC's only
      psubdirs       = ['FC0days','FC2days','FC3days','FC4days']
      figdir         = 'time_series/figs/FC'
      DO_DATA_COMP   = 0
      data_comp_str  = [psubdirs[1],psubdirs[-1]]
      # data_comp_str  = [psubdirs[1],psubdirs[2]]
      # continue
   else:
      # compare the persistence FCs too
      psubdirs       = ['FC0days','FC2days','PFC1days','PFC2days']
      figdir         = 'time_series/figs/PFC'
      DO_DATA_COMP   = 0
      data_comp_str  = [psubdirs[1],psubdirs[-1]]
      # continue

   for sdir in subdirs:

      tsU   = mr.read_time_series('time_series/Model_Underest_'+sdir+'.txt')
      tsO   = mr.read_time_series('time_series/Model_Overest_'+sdir+'.txt')
      tsB   = 'time_series/Model_Bias_'+sdir+'.txt'


      # ============================================================
      # get last time series (bias) from over & under time series
      dates = tsU.dates
      areaU = tsU.data[akey] 
      areaO = tsO.data[akey] 
      areaT = areaU+areaO
      data  = {akey:areaT}
      for wkey in Bfields:
         # get weighted avg for 2 mean widths
         widU  = tsU.data[wkey] 
         widO  = tsO.data[wkey] 
         widB  = (-areaU*widU+areaO*widO)/areaT
         data.update({wkey:widB})

      tsB   = mr.time_series(dates,data,filename=tsB,overwrite=True)
      # ============================================================


   xtlabs   = []
   xticks   = []
   for year in [2016,2017]:
      for M in range(M0,13,Mstep):
         dto   = DT(year,M,1)
         if dto>=dto_lims[0] and dto<=dto_lims[1]:
            xtlabs.append(dto.strftime('%Y%m%d'))
            xticks.append((dto-refdate).total_seconds()*tfac)

   #now make the plots
   fignames    = []
   fignames.append(figdir+'/conc_anomaly_RMSE.png')
   fignames.append(figdir+'/conc_anomaly_bias.png')
   fignames.append(figdir+'/Model_Underest.png')
   fignames.append(figdir+'/Model_Overest.png')
   fignames.append(figdir+'/Model_Bias.png')


   for n,figname in enumerate(fignames):
      pobj  = mr.plot_object()
      if n==0:
         fld      = Cfields[0]
         tstr     = 'time_series/conc_anomaly_'
         yscaling = 1.
         ylabel   = 'RMSE conc. anomaly'
         ylim     = [.1,.25]
      elif n==1:
         fld      = Cfields[1]
         tstr     = 'time_series/conc_anomaly_'
         yscaling = 1.
         ylabel   = 'Bias in conc. anomaly'
         ylim     = [-.1,.1]
      elif n==2:
         fld      = Bfields[0]
         tstr     = 'time_series/Model_Underest_'
         yscaling = 1.e-3 # m to km
         ylabel   = 'Ice edge underestimation, km'
         ylim     = [0,150.]
      elif n==3:
         fld      = Bfields[0]
         tstr     = 'time_series/Model_Overest_'
         yscaling = 1.e-3 # m to km
         ylabel   = 'Ice edge overestimation, km'
         ylim     = [0,150.]
      elif n==4:
         fld      = Bfields[0]
         tstr     = 'time_series/Model_Bias_'
         yscaling = 1.e-3 # m to km
         ylabel   = 'Ice edge bias, km'
         ylim     = [-150,50.]

      lines       = []
      legs        = []
      data_comp   = []
      for sdir in psubdirs:
         ts = mr.read_time_series(tstr+sdir+'.txt')
         pobj,lin,info  = ts.plot(fld,pobj=pobj,refdate=refdate,yscaling=yscaling)
         lines.append(lin)
         legs.append(sdir[:-4])
         if DO_DATA_COMP:
            if sdir in data_comp_str:
               data_comp.append(ts.data[fld])

      pobj.ax.legend(lines,legs)
      pobj.ax.set_xticks(xticks)
      pobj.ax.set_xticklabels(xtlabs,rotation=270)
      pobj.ax.set_ylabel(ylabel)
      pobj.ax.set_ylim(ylim)
      pobj.ax.set_xlim(xlim)


      # ====================================================================
      if DO_DATA_COMP:
         # highlight good or bad time periods
         Ddiff = np.abs(data_comp[1])-np.abs(data_comp[0]) #eg PFC - FC

         if 0:#n==4:
            # plot the difference
            from matplotlib import pyplot as plt
            FIG   = plt.figure()
            AX    = FIG.add_subplot(111)
            AX.plot(Ddiff*yscaling)
            AX.set_ylabel(ylabel)
            AX.set_title(fld+' : '+data_comp_str[1]+' - '+data_comp_str[0])
            plt.show(FIG)
            AX.cla()
            plt.close(FIG)

         jref  = 0
         dt1   = dates[jref]
         dD1   = Ddiff[jref]
         x1    = (dt1-refdate).total_seconds()*tfac
         for jd,ddiff in enumerate(Ddiff):
            if jd==jref:
               continue


            if ddiff*dD1>0:
               # hasn't crossed 0
               continue

            # =================================================
            # has crossed 0 if we get here
            x2 = (dates[jd]-refdate).total_seconds()*tfac
            if dD1>0:
               # PFC worse:
               # - highlight region on plot
               pobj.ax.axvspan(x1,x2,ymin=0,ymax=1,\
                  color='red',alpha=0.3)

            # reset ref point of zero-crossing check
            dt1   = dates[jd]
            dD1   = ddiff
            x1    = x2
            # =================================================

      # ====================================================================

      print('\nSaving '+figname)
      pobj.fig.savefig(figname,bbox_inches='tight')
      if show:
         pobj.fig.show()
      pobj.ax.cla()
      pobj.fig.clear()
