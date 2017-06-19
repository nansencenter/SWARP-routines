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

# time series files and objects
tsfilsC  = {}
tseriesC = {}
tsfilsU  = {}
tseriesU = {}
tsfilsO  = {}
tseriesO = {}


for sdir in subdirs:
   tsfilsC.update({sdir:'time_series/conc_anomaly_'+sdir+'.txt'})
   tsfilsU.update({sdir:'time_series/Model_Underest_'+sdir+'.txt'})
   tsfilsO.update({sdir:'time_series/Model_Overest_'+sdir+'.txt'})

for n,ddir in enumerate(ddirs):
   dto   = datetimes[n]
   print('\n')
   print(dto)

   for sdir in subdirs:
      Sdir  = vdir+'/'+ddir+'/OSISAF/'+sdir

      # ===========================================================================
      # get summary files and read them
      Cfil  = None
      Ofil  = None
      Ufil  = None

      if os.path.exists(Sdir):
         files = os.listdir(Sdir)
         for f in files:
            if 'conc_anomaly_OSISAF_' in f and '.txt' in f:
               Cfil  = Sdir+'/'+f
            if 'Under_summary.txt' in f:
               Ufil  = Sdir+'/'+f
            if 'Over_summary.txt' in f:
               Ofil  = Sdir+'/'+f

      dates = [dto]
      if Cfil is not None:
         tmp,dataC  = mr.read_MIZpoly_summary(Cfil).get_time_series(inputs_only=True)
      else:
         dataC = {}
         for Cfld in Cfields:
            dataC.update({Cfld:[np.nan]})

      if Ufil is not None:
         tmp,dataU  = mr.read_MIZpoly_summary(Ufil).get_time_series(inputs_only=True)
      else:
         dataU = {}
         for Ufld in Ufields:
            dataU.update({Ufld:[np.nan]})

      if Ofil is not None:
         tmp,dataO  = mr.read_MIZpoly_summary(Ofil).get_time_series(inputs_only=True)
      else:
         dataO = {}
         for Ofld in Ofields:
            dataO.update({Ofld:[np.nan]})
      # ===========================================================================

      if n==0:
         tseriesC.update({sdir:mr.time_series(dates,dataC,filename=tsfilsC[sdir],overwrite=True)})
         tseriesU.update({sdir:mr.time_series(dates,dataU,filename=tsfilsU[sdir],overwrite=True)})
         tseriesO.update({sdir:mr.time_series(dates,dataO,filename=tsfilsO[sdir],overwrite=True)})
      else:
         tseriesC[sdir] = tseriesC[sdir].extend(dates,dataC,filename=tsfilsC[sdir],overwrite=True)
         tseriesU[sdir] = tseriesU[sdir].extend(dates,dataU,filename=tsfilsU[sdir],overwrite=True)
         tseriesO[sdir] = tseriesO[sdir].extend(dates,dataO,filename=tsfilsO[sdir],overwrite=True)
