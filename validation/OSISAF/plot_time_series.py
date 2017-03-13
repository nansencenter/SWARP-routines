import os,sys
import mod_reading as mr
from datetime import datetime as DT
import numpy as np

vdir  = '/work/timill/RealTime_Models/validation/TP4a0.12/ice_only'
ddirs = os.listdir(vdir)
"""
daily directories, with contents (eg)
ls 20160106/OSISAF/
FC0days  FC1days  FC2days  FC3days  FC4days  FC5days  PFC1days  PFC2days  PFC3days  PFC4days  PFC5days  figs
"""


# sort the dates
datetimes   = [DT.strptime(ddir,'%Y%m%d') for ddir in ddirs]
DI          = sorted([(dto,i) for i,dto in enumerate(datetimes)])
ddirs       = [ddirs[i] for dto,i in DI]
datetimes   = [dto for dto,i in DI]


subdirs  = ['FC0days' , 'FC1days' , 'FC2days' , 'FC3days' , 'FC4days' , 'FC5days',\
                        'PFC1days', 'PFC2days', 'PFC3days', 'PFC4days', 'PFC5days']
Cfields  = ['RMSE_both_ice', 'Bias_both_ice', 'RMSE_either_ice', 'Bias_either_ice']
Ufields  = ['Maximum_total_width', 'Mean_intersecting_width', 'Mean_total_width',\
            'Maximum_intersecting_width', 'Total_perimeter', 'Total_area']

Ofields  = 1*Ufields
Bfields  = Ufields[1:3]
akey     = Ufields[-1]  # area


psubdirs = []
for sdir in subdirs:
   if 'PFC' not in sdir:
      # for plotting
      psubdirs.append(sdir)

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


#now make the plots
fignames    = []
fignames.append('time_series/conc_anomaly_RMSE.png')
fignames.append('time_series/conc_anomaly_bias.png')
fignames.append('time_series/Model_Underest.png')
fignames.append('time_series/Model_Overest.png')
fignames.append('time_series/Model_Bias.png')


for n,figname in enumerate(fignames):
   pobj  = mr.plot_object()
   if n==0:
      fld      = Cfields[0]
      tstr     = 'time_series/conc_anomaly_'
      yscaling = 1.
   elif n==1:
      fld      = Cfields[1]
      tstr     = 'time_series/conc_anomaly_'
      yscaling = 1.
   elif n==2:
      fld      = Bfields[0]
      tstr     = 'time_series/Model_Underest_'
      yscaling = 1.e-3 # m to km
   elif n==3:
      fld      = Bfields[0]
      tstr     = 'time_series/Model_Overest_'
      yscaling = 1.e-3 # m to km
   elif n==4:
      fld      = Bfields[0]
      tstr     = 'time_series/Model_Bias_'
      yscaling = 1.e-3 # m to km

   lines = []
   legs  = []
   for sdir in psubdirs:
      ts = mr.read_time_series(tstr+sdir+'.txt')
      pobj,lin,info  = ts.plot(fld,pobj=pobj,yscaling=yscaling)
      lines.append(lin)
      legs.append(sdir[:-4])

   pobj.ax.legend(lines,legs)
   print('\nSaving '+figname)
   pobj.fig.savefig(figname)
   # pobj.fig.show()
   pobj.ax.cla()
   pobj.fig.clear()
