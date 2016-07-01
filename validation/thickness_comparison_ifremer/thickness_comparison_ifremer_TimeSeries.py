import os,sys
if 1:
   from   mpl_toolkits.basemap import pyproj
else:
   import pyproj
import mod_reading          as mr
import numpy                as np
import geometry_sphere      as GS
import fns_plotting         as FP
from matplotlib import pyplot as plt
import datetime as dtm


###############################################
def wgs84_ellipsoid_mat():
   a,fi = 6378137,298.257223563
   # f=flattening=(a-b)/a=1/fi
   f    = 1./fi
   b    = (1-f)*a
   ecc  = np.sqrt(1-pow(b/a,2)) # eccentricity
   # print(a,fi,f,b,ecc)
   # print(a/(a-b))
   return a,ecc
###############################################


###############################################

bmap    = FP.start_HYCOM_map('Arctic')
ODLmap  = pyproj.Proj('+init=EPSG:3413')
"""
srs='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45'+\
    ' +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
"""


# directories
compname  = os.getenv('HOSTNAME')
if 'hexagon' in compname:
   #hexagon
   HEX      = True
   mdir     = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice_erai2/'
   gridpath = None


# make time series plots
figdir   = 'out/'
tso1     = mr.read_time_series(figdir+'/thickness_error_winter1.txt')
tso2     = mr.read_time_series(figdir+'/thickness_error_winter2.txt')
vlist    = ['Bias_either_ice','RMSE_either_ice']
figs     = ['thickness_error_winter1.png','thickness_error_winter2.png']
legs     = ['Bias','RMSE']

refdate  = dtm.datetime(2015,1,1)
tfac     = 1./(24*3600) # seconds to days
xticks   = []
xtlabs   = []
for i in range(1,13):
   xticks.append((dtm.datetime(2015,i,1)-refdate).total_seconds()*tfac)
   xtlabs.append(dtm.datetime(2015,i,1).strftime('%Y%m%d'))
   xticks.append((dtm.datetime(2015,i,15)-refdate).total_seconds()*tfac)
   xtlabs.append(dtm.datetime(2015,i,15).strftime('%Y%m%d'))
for i in range(1,3):
   xticks.append((dtm.datetime(2016,i,1)-refdate).total_seconds()*tfac)
   xtlabs.append(dtm.datetime(2016,i,1).strftime('%Y%m%d'))
   xticks.append((dtm.datetime(2016,i,15)-refdate).total_seconds()*tfac)
   xtlabs.append(dtm.datetime(2016,i,15).strftime('%Y%m%d'))

po = None
for i,tso in enumerate([tso1,tso2]):
   Tmin     = (tso.dates[0 ]-refdate).total_seconds()*tfac
   Tmax     = (tso.dates[-1]-refdate).total_seconds()*tfac
   print(Tmin,Tmax)
   Xticks   = []
   Xtlabs   = []
   for j,xt in enumerate(xticks):
      xtl   = xtlabs[j]
      if i==1 and np.mod(j,2)==1:
         continue
      if xt>=Tmin and xt<=Tmax:
         Xticks.append(xt)
         Xtlabs.append(xtl)

   lines    = []
   for v in vlist:
      po,lin   = tso.plot(v,pobj=po,refdate=refdate,yscaling=1.,linestyle='-',marker='^')
      lines.append(lin)

   # legend, axes labels
   po.ax.legend(lines,legs,loc=4)
   # po.ax.set_xlabel('Days since '+refdate.strftime('%Y-%m-%d %H:%M'))
   po.ax.set_xticks(Xticks)
   po.ax.set_xticklabels(Xtlabs,rotation=270)
   po.ax.set_ylabel('Thickness error, m')

   figname  = figdir+'/'+figs[i]
   print('Saving '+figname)
   po.fig.savefig(figname,bbox_inches='tight')
   po.ax.cla()

plt.close(po.fig)
