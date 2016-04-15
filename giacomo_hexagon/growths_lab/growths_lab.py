# Calling in all the needed packages
from netCDF4 import Dataset
import sys,os
import time
import datetime
import glob
import numpy as np
import subprocess
import shutil
import matplotlib
#NOTE to be used only on servers (i.e. Hexagon)
matplotlib.use('Agg')
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from skimage import measure as msr
from skimage import morphology as morph
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

sys.path.append('./utilities')
import mod_reading as Mrdg
import Laplace_eqn_solution_2tw as Leqs
import test_HYCOM_diag as Hdiag
import MIZchar as widths

class reader:
    def __init__(self,dadate,bmm):

        self.year = dadate[:4]
        self.month = dadate[4:6]
        self.day = dadate[6:8]
        gigio = datetime.datetime(int(float(self.year)),int(float(self.month)),int(float(self.day)),0,0)
        gigio = gigio.strftime('%j')
        gigio = int(float(gigio))
        self.day_j = gigio
        self.jday_start(gigio)
        self.tidx_i,self.tidx_w = self.time_frame(gigio)

        #Reading the datasets
        #self._read_osi_(dadate,bmo)
        #self._read_osi_2(dadate,bmo)
        #self._read_mdl_DAILY(self.data_date,dadate,bmm)
        self._read_mdl_ice_only(self.data_date,bmm)
        self._read_mdl_waves_ice(self.data_date,bmm)

    def jday_start(self,jday):
        if jday < 60 or jday > 274:
            print 'ERROR: dataset not in the melting season'
        else:
            std = 60
            check = 0
            while (check < 1):
               if jday < (std+8):
                  sjd = std
                  check = 5
               else:
                  std += 7
        data_date = datetime.datetime(2015,1,1)+datetime.timedelta(sjd)
        data_date = data_date.strftime('%Y%m%d')
        self.jdaystart = sjd
        self.data_date = data_date
        return()

    def time_frame(self,jday):
        dd = jday - self.jdaystart
        tf_i = 4 + ((dd -1)*8)
        tf_w = 2 + ((dd -1)*4)
        return(tf_i,tf_w)

    def _read_osi_(self,dadate,basemap):
        day = self.day
        month = self.month
        year = self.year
        # Read in OSI_SAF file
        outdir = '/work/timill/giacomo/osisaf_repo'
        ncfil = outdir+'/ice_conc_nh_polstere-100_multi_'+dadate+'1200.nc'
        clon = 'lon'
        clat = 'lat'
        cconc = 'ice_conc'
        lonO = Mrdg.nc_get_var(ncfil,clon) # lon[:,:] is a numpy array
        latO = Mrdg.nc_get_var(ncfil,clat) # lat[:,:] is a numpy array
        conc = Mrdg.nc_get_var(ncfil,cconc,time_index=0)
        xc = Mrdg.nc_get_var(ncfil,'xc')
        yc = Mrdg.nc_get_var(ncfil,'yc')

        X2,Y2 = basemap(lonO[:,:],latO[:,:],inverse=False)
        XO = np.copy(X2)
        YO = np.copy(Y2)

        Z2 = conc[:,:]
        ZO = Z2/100

        self.lonO = lonO[:,:]
        self.latO = latO[:,:]
        self.XO = XO[:,:]
        self.YO = YO[:,:]
        self.ZO = ZO

        print('OSI-SAF: lonO,latO,XO,YO,ZO')
        
        return() 
   
    def _read_osi_2(self,dadate,basemap):

        day = self.day
        month = self.month
        year = self.year
        dadate2 = datetime.date(int(year),int(month),int(day))+datetime.timedelta(1)
        dadate2 = dadate2.strftime('%Y%m%d')
        # Read in OSI_SAF file
        outdir = '/work/timill/giacomo/osisaf_repo'
        ncfil = outdir+'/ice_conc_nh_polstere-100_multi_'+dadate2+'1200.nc'
        clon = 'lon'
        clat = 'lat'
        cconc = 'ice_conc'
        lonO = Mrdg.nc_get_var(ncfil,clon) # lon[:,:] is a numpy array
        latO = Mrdg.nc_get_var(ncfil,clat) # lat[:,:] is a numpy array
        conc = Mrdg.nc_get_var(ncfil,cconc,time_index=0)
        xc = Mrdg.nc_get_var(ncfil,'xc')
        yc = Mrdg.nc_get_var(ncfil,'yc')

        X2,Y2 = basemap(lonO[:,:],latO[:,:],inverse=False)
        XO = np.copy(X2)
        YO = np.copy(Y2)

        Z2 = conc[:,:]
        ZO = Z2/100

        self.lonO2 = lonO[:,:]
        self.latO2 = latO[:,:]
        self.XO2 = XO[:,:]
        self.YO2 = YO[:,:]
        self.ZCO2 = ZO

        print('OSI-SAF plus 1 day: lonO2,latO2,XO2,YO2,ZCO2')
        
        return()

    def _read_mdl_DAILY(self,data_date,dadate,basemap):
        
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        jdaystart = self.jdaystart

        if jdaystart < 100:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                 'wavesice/2015_GOOD/2015_0'+str(jdaystart)+'/netcdf/DAILY/'
        else:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                 'wavesice/2015_GOOD/2015_'+str(jdaystart)+'/netcdf/DAILY/'

        outdir = '/work/timill/giacomo/DAILY_WIM_2015_melting'
        ncfil = outdir+'/TP4DAILY_start'+data_date+'_dump'+dadate+'.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'fice'
        sthic = 'hice'
        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        lonM = lonM[:,:]
        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        latM = latM[:,:]
        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
        self.lonM,self.latM,self.XM,self.YM = lonM,latM,X,Y
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
        conc = conc[:,:]
        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=0)
        thic = thic[:,:]

        self.ZM,self.ZT = conc,thic
        
        print('DAILY: lonM,latM,X,Y,ZC,ZT')

        return()

    def _read_mdl_ice_only(self,data_date,basemap):
        
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        jdaystart = self.jdaystart

        # TIME INDEX
        tidx = self.tidx_i

        # Read TP4arch
        if jdaystart < 100:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                   'ice_only/2015_GOOD/2015_0'+str(jdaystart)+'/final_output'
        else:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                   'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+data_date+'T000000Z.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'icec'
        sthic = 'icetk'
        shsnow = 'hsnow'
        sqtot = 'qtot'
        sqcool = 'flx_cool'
        sqother = 'flx_itop'
        sqatm = 'flx_ow'
        sconc_old = 'fi_old'
        sthic_old = 'hi_old'
        shsnow_old = 'hs_old'
        shsnow_diff = 'hs_int'
        sconc_diff = 'dfice'
        sthic_diff = 'dhice'
        staux = 'taux'
        stauy = 'tauy'
        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        lonM = lonM[:,:]
        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        latM = latM[:,:]
        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
        self.lonM,self.latM,self.X,self.Y = lonM,latM,X,Y
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
        conc = conc[:,:]
        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
        thic = thic[:,:]
        hsnow = Mrdg.nc_get_var(ncfil,shsnow,time_index=tidx)
        hsnow = hsnow[:,:]
        qtot = Mrdg.nc_get_var(ncfil,sqtot,time_index=tidx)
        qtot = qtot[:,:]
        qcool = Mrdg.nc_get_var(ncfil,sqcool,time_index=tidx)
        qcool = qcool[:,:]
        qother = Mrdg.nc_get_var(ncfil,sqother,time_index=tidx)
        qother = qother[:,:]
        qatm = Mrdg.nc_get_var(ncfil,sqatm,time_index=tidx)
        qatm = qatm[:,:]
        conc_old = Mrdg.nc_get_var(ncfil,sconc_old,time_index=tidx)
        conc_old = conc_old[:,:]
        thic_old = Mrdg.nc_get_var(ncfil,sthic_old,time_index=tidx)
        thic_old = thic_old[:,:]
        hsnow_old = Mrdg.nc_get_var(ncfil,shsnow_old,time_index=tidx)
        hsnow_old = hsnow_old[:,:]
        hsnow_diff = Mrdg.nc_get_var(ncfil,shsnow_diff,time_index=tidx)
        hsnow_diff = hsnow_diff[:,:]
        conc_diff = Mrdg.nc_get_var(ncfil,sconc_diff,time_index=tidx)
        conc_diff = conc_diff[:,:]
        thic_diff = Mrdg.nc_get_var(ncfil,sthic_diff,time_index=tidx)
        thic_diff = thic_diff[:,:]
        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
        taux = taux[:,:]
        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
        tauy = tauy[:,:]

        self.ZMi,self.ZTI,self.HSI,self.QT,self.QC,self.QO,self.QA,self.FOI,\
                self.TOI,self.HSO,self.HSD,self.FDI,self.TDI,self.TXI,self.TYI =\
                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
                hsnow_diff,conc_diff,thic_diff,taux,tauy
        
        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')

        return()

    def _read_mdl_waves_ice(self,data_date,basemap):
    
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        jdaystart = self.jdaystart

        # TIME INDEX
        tidx = self.tidx_w

        # Read TP4arch_wav
        if jdaystart < 100:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                   'wavesice/2015_GOOD/2015_0'+str(jdaystart)+'/final_output'
        else:
           outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                   'wavesice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'

        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+data_date+'T000000Z.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'icec'
        sthic = 'icetk'
        sdmax = 'dmax'
        sswh = 'swh'
        staux = 'taux_wav'
        stauy = 'tauy_wav'
        smwd = 'mwd'
        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        lonM = lonM[:,:]
        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        latM = latM[:,:]
        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
        self.lonM,self.latM,self.XM,self.YM = lonM,latM,X,Y
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
        conc = conc[:,:]
        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
        thic = thic[:,:]
        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=tidx)
        dmax = dmax[:,:]
        swh = Mrdg.nc_get_var(ncfil,sswh,time_index=tidx)
        swh = swh[:,:]
        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
        taux = taux[:,:]
        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
        tauy = tauy[:,:]
        mwd = Mrdg.nc_get_var(ncfil,smwd,time_index=tidx)
        mwd = mwd[:,:]
        self.ZCW,self.ZTW,self.ZDW,self.SWHW,self.MWDW,self.TXW,self.TYW =\
                conc,thic,dmax,swh,mwd,taux,tauy

        print('Wave Ice: lonM,latM,X,Y,ZC,ZT,ZD,swh,mwd,taux,tauy')
       
        return()

def basemap_creator(region,cres='i'):
    
    if region=='TP4':
        lonc     = -45.
        latc     = 90.
        lat_ts   = latc
        rad      = 33 # radius in deg
        width    = 2*rad*111.e3
        height   = 2*rad*111.e3
        #
        bm = Basemap(width=width,height=height,\
                     resolution=cres,projection='stere',\
                     lat_ts=lat_ts,lat_0=latc,lon_0=lonc)
    
    elif region=='BS1':
        lonc     = 48.
        latc     = 74.
        lat_ts   = latc
        rad      = 12 # radius in deg
        width    = 2*rad*111.e3
        height   = 2*rad*111.e3
        #
        bm = Basemap(width=width,height=height,\
                     resolution=cres,projection='stere',\
                     lat_ts=lat_ts,lat_0=latc,lon_0=lonc)
    
    elif region=='FR1':
        lonc     = 0.5
        latc     = 78.75
        lat_ts   = latc
        rad      = 6.0 # radius in deg
        width    = 2*rad*111.e3
        height   = 2*rad*111.e3
        #
        bm = Basemap(width=width,height=height,\
                     resolution=cres,projection='stere',\
                     lat_ts=lat_ts,lat_0=latc,lon_0=lonc)

    elif region=='OSI':
        lonc    = -45
        latc    = 90
        lat_ts  = 70
        width   = 7600000
        height  = 11200000
    
        bm  = Basemap(width=width,height=height,\
                resolution=cres,rsphere=(6378273,6356889.44891),\
                projection='stere',lat_ts=lat_ts,lat_0=latc,lon_0=lonc)
    
    return bm

dadate = sys.argv[1] 
x = sys.argv[2]
y = sys.argv[3]

bmm = basemap_creator('TP4')
data = reader(dadate,bmm)

dmax  = data.ZDW[x,y]
print dmax
f_old = data.FOI[x,y]
h_old = data.TOI[x,y]
hs_old = data.HSO[x,y]
Qatm = data.QA[x,y]
Qother = data.QO[x,y]
Qcool = data.QC[x,y]

fusi  = 3.02e8
fuss  = 1.10e8
dmin  = 20.
dt    = 400.

print(Qatm*(1-f_old),Qother*f_old)
print(Qatm*(1-f_old)*dt/fusi,Qother*f_old*dt/fusi)
print(' ')

if dmax>200:
   beta_lat    = (4*h_old*dmax)/pow(dmax,2)
   alpha_lat   = beta_lat/(1+beta_lat)
elif dmax>dmin:
   M = np.log(dmax/dmin)
   M = M/np.log(2)
   M = np.floor(M)
   for el in range(M+1):
      nm = N0*(1-f)*(f*xi**2)**el
      dm = linZD[n]/(xi**el)
      nsum = nsum + nm 
      ndsum = ndsum + (nm*dm)
      ndsum2 = ndsum2 + (nm*dm**2)
   linDmean[n] = ndsum/nsum
   linD2mean[n] = ndsum2/nsum
   linSlat[n] = N0*4*ndsum*linthic_old[n]
   linSbot[n] = N0*ndsum2
   linAlpha[n] = linSlat[n]/(linSlat[n]+linSbot[n])
   linBeta[n] = linSlat[n]/linSbot[n]

qdist = (Qcool+(1-f_old)*Qatm)*dt
qlat  = alpha_lat*qdist
qvrt  = (1-alpha_lat)*qdist+f_old*Qother*dt

df       = -qlat/(fusi*h_old+fuss*hs_old)
dV_lat   = h_old*df
dh       = -qvrt/fusi
dV_vrt   = f_old*dh
Glat     = dV_lat/dt
Gvrt     = dV_vrt/dt
Glath    = data.TOI[x,y]*data.FDI[x,y]/dt
Gvrth    = data.COI[x,y]*data.TDI[x,y]/dt

print(qlat,qvrt)
print(df,dh)
print(Glat,Gvrt)

fac   = 24*3600*1e3 #m/s ->mm/day
print('\nmm/day WIM')
print(Glat*fac,Gvrt*fac)
print('\nmm/day HYCOM')
print(Glat*fac,Gvrt*fac)


