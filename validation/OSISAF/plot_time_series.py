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


if 1:
   # choose what to plot
   psubdirs = ['FC0days','FC2days','FC3days','FC4days']
   figdir   = 'time_series/figs/FC'
else:
   # choose what to plot
   psubdirs = ['FC0days','FC2days','PFC1days','PFC2days']
   figdir   = 'time_series/figs/PFC'

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

refdate  = dates[0]
xtlabs   = []
xticks   = []
tfac     = 1./(24*3600) # seconds to days
dto_end  = DT(2017,3,12)
for year in [2016,2017]:
   for M in range(2,13,2):
      dto   = DT(year,M,1)
      if dto<dto_end:
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

   lines = []
   legs  = []
   for sdir in psubdirs:
      ts = mr.read_time_series(tstr+sdir+'.txt')
      pobj,lin,info  = ts.plot(fld,pobj=pobj,yscaling=yscaling)
      lines.append(lin)
      legs.append(sdir[:-4])

   pobj.ax.legend(lines,legs)
   pobj.ax.set_xticks(xticks)
   pobj.ax.set_xticklabels(xtlabs,rotation=270)
   pobj.ax.set_ylabel(ylabel)
   pobj.ax.set_ylim(ylim)

   print('\nSaving '+figname)
   pobj.fig.savefig(figname,bbox_inches='tight')
   pobj.fig.show()
   pobj.ax.cla()
   pobj.fig.clear()
