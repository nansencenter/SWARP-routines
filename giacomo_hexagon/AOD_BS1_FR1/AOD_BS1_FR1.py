############################################################################
# THIS script will:
# 1) Read in WIM ice_only and waves_ice
# 2) (optional) read in OSI
# 3) (optional) read in TP4 reanalysis
# 4) Run an analysis of the MIZ
# 5) Run an analysis of growths and stresses for the WIM
# 6) (optional) run an analysis of the AOD
############################################################################
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

# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#CLASSES
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Data reader
class reader:
    def __init__(self,dadate,bmm,bmo):

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
        self._read_osi_(dadate,bmo)
        #self._read_osi_2(dadate,bmo)
        self._read_mdl_DAILY(self.data_date,dadate,bmm)
        #self._read_mdl_ice_only(self.data_date,bmm)
        #self._read_mdl_waves_ice(self.data_date,bmm)

    def jday_start(self,jday):
        if jday < 60 or jday > 274:
            print 'ERROR: dataset not in the melting season'
        else:
           std = 60
           check = 0
           while check < 1:
              if jday < (std+8):
                 stj = std
                 check = 5
              else:
                 std += 7
        data_date = datetime.datetime(2015,1,1)+datetime.timedelta(stj)
        data_date = data_date.strftime('%Y%m%d')
        print data_date
        self.jdaystart = stj
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

        #if jdaystart < 100:
        #   outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
        #         'ice_only/2015_GOOD/2015_0'+str(jdaystart)+'/netcdf/DAILY/'
        #else:
        #   outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
        #         'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/netcdf/DAILY/'

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

# From WIMchar
class WIM_growth_stress:
    def __init__(self,obasemap,mbasemap,XMo,YMo,ZM,XO,YO,ZO,TXI,TYI,TXW,TYW,FOI,TOI,ZDW,QA,QC,QO,TDI,FDI,HSD,region):

        ZMo = reproj_mdl2osi(XMo,YMo,ZM,XO,YO)

        mrg = np.load('./geo_info/bar_gre_bal_reg.npy')
        org = np.load('./geo_info/bar_gre_bal_reg_OSI.npy')
        mms = np.load('./geo_info/bar_gre_bal_reg_mask.npy')
        oms = np.load('./geo_info/bar_gre_bal_reg_mask_OSI.npy')

        ZMo,ZO,TXI,TYI,TXW,TYW,FOI,TOI,ZDW,QA,QC,QO,TDI,FDI,HSD = \
                mask_region(ZMo,region,mrg,org,mms,oms),mask_region(ZO,region,mrg,org,mms,oms),\
                mask_region(TXI,region,mrg,org,mms,oms),mask_region(TYI,region,mrg,org,mms,oms),\
                mask_region(TXW,region,mrg,org,mms,oms),mask_region(TYW,region,mrg,org,mms,oms),\
                mask_region(FOI,region,mrg,org,mms,oms),mask_region(TOI,region,mrg,org,mms,oms),\
                mask_region(ZDW,region,mrg,org,mms,oms),mask_region(QA,region,mrg,org,mms,oms),\
                mask_region(QC,region,mrg,org,mms,oms),mask_region(QO,region,mrg,org,mms,oms),\
                mask_region(TDI,region,mrg,org,mms,oms),mask_region(FDI,region,mrg,org,mms,oms),\
                mask_region(HSD,region,mrg,org,mms,oms)

        self.ZMo = ZMo

        # Useful variables
        self.grid_cell_area = 1
        self.delta_t = 400

        self.poly_ovr_und(ZMo,ZO)

        # Recalculation of fluxes and growths
        #self.OSI_ice_growth(self.ZO,self.ZCO2)

        self.stress_h,self.stress_h_ovr,self.stress_h_und = self.stresses(TXI,TYI,obasemap)
        self.stress_w,self.stress_w_ovr,self.stress_w_und = self.stresses(TXW,TYW,obasemap)
        self.Dmean,self.D2mean,self.Slat,self.Sbot,self.Alpha,self.Beta = self.floes_analysis(FOI,TOI,ZDW)
        self.qlat,self.qvrt = self.heat_flux(QA,QC,QO,FOI)
        self.glat_HYCOM,self.gvrt_HYCOM,self.glat_WIM,self.gvrt_WIM,self.CNW,self.TNW =\
        self.WIM_ice_growth(TOI,TDI,FOI,FDI,HSD)

        # Save the data
        self.save_data(region)

        return
    
    def poly_ovr_und(self,conc_m,conc_o):
        mdl = np.ma.copy(conc_m)
        osi = np.ma.copy(conc_o)
        mdl[mdl<.15] = 0
        osi[osi<.15] = 0
        mdl[mdl>=.15] = 1
        osi[osi>=.15] = 1
        polys = osi - mdl
        ovr = np.ma.copy(polys)
        und = np.ma.copy(polys)
        ovr[ovr==1] = 0
        und[und==-1] = 0
        ovr = abs(ovr)
        und = abs(und)
        ovr = closing(ovr)
        und = closing(und)
        self.ovr_poly = ovr
        self.und_poly = und
        
        return(ovr,und)
   
    def floes_analysis(self,conc_old,thic_old,ZD,mode=1):

        # COSTANTS/VARIABLES related
        # In how many pieces will the floe break? That is xi^2.
        xi = 2
        
        # Starting number of floes
        N0 = 1
        
        # Fragility of the ice. Having it variable would be awesome.
        f = 0.9
        
        # Extreme of the Distribution
        dmax = 200
        dmin = 20
        
        # Grid Cell
        Asq = self.grid_cell_area

        # New arrays (using ZD as "template")
        Alpha = np.ma.zeros(ZD.shape)
        Alpha = np.ma.masked_where(np.ma.getmask(ZD),Alpha)
        Beta= np.ma.zeros(ZD.shape)
        Beta = np.ma.masked_where(np.ma.getmask(ZD),Beta)
        Slat = np.ma.zeros(ZD.shape)
        Slat = np.ma.masked_where(np.ma.getmask(ZD),Slat)
        Sbot = np.ma.zeros(ZD.shape)
        Sbot = np.ma.masked_where(np.ma.getmask(ZD),Sbot)
        Dmean = np.ma.zeros(ZD.shape)
        Dmean = np.ma.masked_where(np.ma.getmask(ZD),Dmean)
        D2mean = np.ma.zeros(ZD.shape)
        D2mean = np.ma.masked_where(np.ma.getmask(ZD),D2mean)
        
        # Reshaping for cycle
        linconc_old = np.reshape(conc_old,conc_old.size)
        linthic_old = np.reshape(thic_old,thic_old.size)
        linZD = np.reshape(ZD,ZD.size)
        linAlpha = np.reshape(Alpha,Alpha.size)
        linBeta = np.reshape(Beta,Beta.size)
        linSlat = np.reshape(Slat,Slat.size)
        linSbot = np.reshape(Sbot,Sbot.size)
        linDmean = np.reshape(Dmean,Dmean.size)
        linD2mean = np.reshape(D2mean,D2mean.size)
      

        if mode == 1:
            # fsd with continuous PDF (unlike RG method)
            # * PDF function is A*D^{-(1+alpha)}/(D_min^{-alpha}-D_max^{-alpha}) if D\in[Dmin,Dmax] else 0
            # * cumulative probability function is P(d>D) = (D^{-alpha}-D_max^{-alpha})/(D_min^{-alpha}-D_max^{-alpha})
            
            alpha = np.log(xi*xi*f)/np.log(xi)

            for n,en in enumerate(linZD):
                if en != 0 and en > 20:
                    A = alpha/(pow(dmin,-alpha)-pow(en,-alpha))
                    linDmean[n] = A/(alpha-1)*(pow(dmin,1-alpha)-pow(en,1-alpha))
                    linD2mean[n] = A/(alpha-2)*(pow(dmin,2-alpha)-pow(en,2-alpha))
                    linBeta[n] = 4*linthic_old[n]*linDmean[n]/linD2mean[n]
                    linAlpha[n] = linBeta[n]/(1+linBeta[n])
                    linSbot[n] = np.ma.copy(linconc_old[n])*Asq
                    linSlat[n] = linBeta[n]*linSbot[n]
                #elif en >= 200:
                #    linDmean[n] = 200
                #    linD2mean[n] = 40000
                #    linBeta[n] = 4*linthic_old[n]*linDmean[n]/linD2mean[n]
                #    linAlpha[n] = linBeta[n]/(1+linBeta[n])
                #    linSbot[n] = np.ma.copy(linconc_old[n])*Asq
                #    linSlat[n] = linBeta[n]*linSbot[n]
                #elif en == 0:
                #elif en <= 20 and linconc_old[n] >= .15:
                #    en = 21
                #    A = alpha/(pow(dmin,-alpha)-pow(en,-alpha))
                #    linDmean[n] = A/(alpha-1)*(pow(dmin,1-alpha)-pow(en,1-alpha))
                #    linD2mean[n] = A/(alpha-2)*(pow(dmin,2-alpha)-pow(en,2-alpha))
                #    linBeta[n] = 4*linthic_old[n]*linDmean[n]/linD2mean[n]
                #    linAlpha[n] = linBeta[n]/(1+linBeta[n])
                #    linSbot[n] = np.ma.copy(linconc_old[n])*Asq
                #    linSlat[n] = linBeta[n]*linSbot[n]
                else:
                    linDmean[n] = 0
                    linD2mean[n] = 0
                    linBeta[n] = 0
                    linAlpha[n] = 0
                    linSbot[n] = 0
                    linSlat[n] = 0

            Dmean = np.reshape(linDmean,Dmean.shape)
            self.D2mean = np.reshape(linD2mean,D2mean.shape)
            Slat = np.reshape(linSlat,Slat.shape)
            Sbot = np.reshape(linSbot,Sbot.shape)
            Alpha = np.reshape(linAlpha,Alpha.shape)
            Beta = np.reshape(linBeta,Beta.shape)

            self.Dmean = Dmean
            self.D2mean = D2mean
            self.Slat = Slat
            self.Sbot = Sbot
            self.Alpha = Alpha
            self.Beta = Beta
           
        elif mode == 2:
            # Number of iterations
            ZD = np.copy(ZD)
            ZD[ZD<=dmin] = dmin
            nM = np.log(ZD/dmin)
            M = nM/np.log(2)
            M[M>(np.log(dmax/dmin)/np.log(2))] = 999
            M[M<0] = 0
            M = np.floor(M)
            M = np.ma.masked_where(np.ma.getmask(ZD),M)
            self.frag_iter = M

            linM = np.reshape(M,M.size)
            
            for n,en in enumerate(linM):
                nsum = 0
                ndsum = 0
                ndsum2= 0
                if en > 0 and en != 999 and en is not None:
                    nen = np.int(en)
                    for el in range(nen+1):
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
                elif en == 999 and linZD[n] != 0 and linthic_old[n] != 0:
                    linDmean[n] = linZD[n]
                    linD2mean[n] = linZD[n]**2
                    linSlat[n] = 4*linthic_old[n]*N0*linZD[n]
                    linSbot[n] = N0*linD2mean[n]
                    linAlpha[n] = linSlat[n]/(linSlat[n]+linSbot[n])
                    linBeta[n] = linSlat[n]/linSbot[n]
                else:
                    linDmean[n] = 0
                    linD2mean[n] = 0
                    linSbot[n] = 0
                    linSlat[n] = 0
                    linAlpha[n] = 0
                    linBeta[n] = 0
                
            Dmean = np.reshape(linDmean,Dmean.shape)
            self.D2mean = np.reshape(linD2mean,D2mean.shape)
            Slat = np.reshape(linSlat,Slat.shape)
            Sbot = np.reshape(linSbot,Sbot.shape)
            Alpha = np.reshape(linAlpha,Alpha.shape)
            Beta = np.reshape(linBeta,Beta.shape)
            
            Sbot = np.ma.copy(conc_old)*Asq
            Slat = Beta*Sbot
            
            self.Dmean = Dmean
            self.D2mean = D2mean
            self.Slat = Slat
            self.Sbot = Sbot
            self.Alpha = Alpha
            self.Beta = Beta
                
        return(Dmean,D2mean,Slat,Sbot,Alpha,Beta)
        
    def heat_flux(self,qatm,qcool,qother,f_old):
        
        Alpha = self.Alpha
        dt = self.delta_t
        Asq = self.grid_cell_area

        qdist = (qcool + (qatm*(1-f_old)))*Asq*dt
        qlat = Alpha*qdist
        qvrt = qother*f_old*Asq*dt + (1-Alpha)*qdist
        self.qlat = qlat
        self.qvrt = qvrt

        return(qlat,qvrt)

    def WIM_ice_growth(self,thic_old,thic_new,conc_old,conc_new,snow_new):
        
        # Needed variables
        qlat = self.qlat
        qvrt = self.qvrt
        alpha = self.Alpha
        beta = self.Beta
        Asq = self.grid_cell_area
        dt = self.delta_t

        # Latent heats per volume
        Lsice = 3.02*10**8
        Lsnow = 1.10*10**8
        fmax = .995

        # Linearizerz!
        linqlat = np.reshape(qlat,qlat.size)
        linqvrt = np.reshape(qvrt,qvrt.size)
        linbeta = np.reshape(beta,beta.size)
        lintoi = np.reshape(thic_old,thic_old.size)
        lintni = np.reshape(thic_new,thic_new.size)
        lincoi = np.reshape(conc_old,conc_old.size)
        lincni = np.reshape(conc_new,conc_new.size)
        lintsn = np.reshape(snow_new,snow_new.size)
        linqlat = np.reshape(qlat,qlat.size)
        linqvrt = np.reshape(qvrt,qvrt.size)

        # Creating the new conc and thic differences
        CNW = np.ma.zeros(conc_old.shape)
        CNW = np.ma.masked_where(np.ma.getmask(conc_old),CNW)
        lincnw = np.reshape(CNW,CNW.size)

        TNW = np.ma.zeros(thic_old.shape)
        TNW = np.ma.masked_where(np.ma.getmask(thic_old),TNW)
        lintnw = np.reshape(TNW,TNW.size)

        for n,en in enumerate(linqlat):
            if en < 0 and en is not None:
                lincnw[n] = -linqlat[n]/float(lintoi[n]*Asq*Lsice)
                f_tmp = lincnw[n]+lincoi[n]
                teff = 0
                if  f_tmp > fmax:
                    lincnw[n] = fmax - lincoi[n]
                    teff = (lintoi[n]/float(fmax))*(f_tmp - fmax)
                lintnw[n] = -(linqvrt[n]/float(lincoi[n]*Asq*Lsice)) + teff
            elif en > 0 and en is not None:
                lincnw[n] = -linqlat[n]/float((Lsice*lintoi[n]+Lsnow*lintsn[n])*Asq)
                lintnw[n] = -(linqvrt[n]/float(lincoi[n]*Asq*Lsice))
            elif en == 0 and en is not None:
                lincnw[n] = 0
                if lincoi[n] != 0:
                    lintnw[n] = -(linqvrt[n]/float(lincoi[n]*Asq*Lsice))
                else:
                    lintnw[n] = 0

        glath = lintoi*lincni/dt
        gvrth = lincoi*lintni/dt
        glatw = lintoi*lincnw/dt
        gvrtw = lincoi*lintnw/dt

        CNW = np.reshape(lincnw,CNW.shape)
        TNW = np.reshape(lintnw,TNW.shape)

        glat_HYCOM = np.reshape(glath,qlat.shape)
        gvrt_HYCOM = np.reshape(gvrth,qvrt.shape)
        glat_WIM = np.reshape(glatw,qlat.shape)
        gvrt_WIM = np.reshape(gvrtw,qvrt.shape)

        self.CNW = CNW
        self.TNW = TNW
        self.glat_HYCOM = glat_HYCOM
        self.gvrt_HYCOM = gvrt_HYCOM
        self.glat_WIM = glat_WIM
        self.gvrt_WIM = gvrt_WIM

        return(glat_HYCOM,gvrt_HYCOM,glat_WIM,gvrt_WIM,CNW,TNW)

    def OSI_ice_growth(self,conc_old,conc_new):
        C1 = np.ma.copy(conc_old)
        C2 = np.ma.copy(conc_new)
        DC = C2 - C1
        #GLat = DC.data*10^5

        self.diff_OSI = DC
        #self.glat_OSI = Glat
        return(DC)#,Glat)

    def stresses(self,tx,ty,basemap):
        tau = np.sqrt(tx**2+ty**2)
        ntau = reproj_mdl2osi(XMo,YMo,tau,XO,YO)
        ovr_tau = ntau * self.ovr_poly
        und_tau = ntau * self.und_poly
        return(tau,ovr_tau,und_tau)

    def save_data(self,region):

        # REGION
        reg_repo = './outputs/WIM/'+str(region)
        if not os.path.exists(reg_repo):
            os.mkdir(reg_repo)

        stressdir = reg_repo+'/stress'
        if not os.path.exists(stressdir):
            os.mkdir(stressdir)
        
        stressh = linearizer(self.stress_h)
        stressh_ovr = linearizer(self.stress_h)
        stressh_und = linearizer(self.stress_h)
        stressw = linearizer(self.stress_w)
        stressw_ovr = linearizer(self.stress_w)
        stressw_und = linearizer(self.stress_w)

        filname = str(stressdir)+'/'+str(dadate)+'_h.txt'
        with open(filname,'a') as f:
            for l,el in enumerate(stressh):
                f.write("%s %s %s" % (str(el), str(stressh_ovr[l]), str(stressh_und[l])))
                f.write('\n')
            f.close()

        filname = str(stressdir)+'/'+str(dadate)+'_w.txt'
        with open(filname,'a') as f:
            for l,el in enumerate(stressw):
                f.write("%s %s %s" % (str(el), str(stressw_ovr[l]), str(stressw_und[l])))
                f.write('\n')
            f.close()

        glath,gvrth = linearizer(self.glat_HYCOM),linearizer(self.gvrt_HYCOM)
        glatw,gvrtw = linearizer(self.glat_WIM),linearizer(self.gvrt_WIM)

        growthdir = reg_repo+'/growth'
        if not os.path.exists(growthdir):
            os.mkdir(growthdir)

        filname = str(growthdir)+'/'+str(dadate)+'_h.txt'
        with open(filname,'a') as f:
            for l,el in enumerate(glath):
                f.write("%s %s" % (str(el),str(gvrth[l])))
                f.write('\n')
            f.close()

        filname = str(growthdir)+'/'+str(dadate)+'_w.txt'
        with open(filname,'a') as f:
            for l,el in enumerate(glatw):
                f.write("%s %s" % (str(el),str(gvrtw[l])))
                f.write('\n')
            f.close()

        # STRESS STATS
        h_stress_max = np.nanmax(self.stress_h)
        w_stress_max = np.nanmax(self.stress_w)
        h_stress_min = np.nanmin(self.stress_h)
        w_stress_min = np.nanmin(self.stress_w)
        h_stress_avg = np.nanmean(self.stress_h)
        w_stress_avg = np.nanmean(self.stress_w)

        # GROWTH STATS
        h_gvrt_max = np.nanmax(self.gvrt_HYCOM)
        w_gvrt_max = np.nanmax(self.gvrt_WIM)
        h_gvrt_min = np.nanmin(self.gvrt_HYCOM)
        w_gvrt_min = np.nanmin(self.gvrt_WIM)
        h_glat_max = np.nanmax(self.glat_HYCOM)
        w_glat_max = np.nanmax(self.glat_WIM)
        h_glat_min = np.nanmin(self.glat_HYCOM)
        w_glat_min = np.nanmin(self.glat_WIM)
        h_gvrt_avg = np.nanmean(self.gvrt_HYCOM)
        w_gvrt_avg = np.nanmean(self.gvrt_WIM)
        h_glat_avg = np.nanmean(self.glat_HYCOM)
        w_glat_avg = np.nanmean(self.glat_WIM)
    
        filname = reg_repo+'/stress_stats.txt'
        with open(filname,'a') as f:
            row1 = [dadate,h_stress_max,h_stress_min,h_stress_avg,w_stress_max,w_stress_min,w_stress_avg]
            str1 = ' '.join(map(str,row1))
            f.write(str1+'\n')
            f.close()

        filname = reg_repo+'/growth_stats.txt'
        with open(filname,'a') as f:
            row1 = [dadate,h_gvrt_max,h_gvrt_min,h_gvrt_avg,w_gvrt_max,w_gvrt_min,w_gvrt_avg,\
                    h_glat_max,h_glat_min,h_glat_avg,w_glat_max,w_glat_min,w_glat_avg]
            str1 = ' '.join(map(str,row1))
            f.write(str1+'\n')
            f.close()

        return()

# From MIZwidth
class MIZwidth:
    def __init__(self,X,Y,Z,basemap,region):
        ZM = mask_region(Z,region)
        if np.nanmax(ZM) < 1:
            MIZraw,ice_ext,pack_ice = self.binary_diff(ZM,.15,.80)
            self.MIZraw,self.ice_ext,self.pack_ice = MIZraw,ice_ext,pack_ice
        else:
            MIZraw,ice_ext,pack_ice = self.binary_diff(ZM,.1,300)
            self.MIZraw,self.ice_ext,self.pack_ice = MIZraw,ice_ext,pack_ice

        MIZcont = self.poly_maker(self.MIZraw)
        self.MIZcont = MIZcont

        MIZ = self.mapper(self.MIZraw,Z)
        self.MIZ = MIZ

        # Analysing the polygons one by one, dividing them into their respective regions
        poly_list = []
        area_list = []
        perim_list = []
        clonlat_list = []
        for n,cont in enumerate(MIZcont):
            polygon,area,perim,clonlat = self.poly_stat(cont,pack_ice,ice_ext,MIZ,X,Y,n,region,basemap=basemap)
            poly_list.append(polygon)
            area_list.append(area)
            perim_list.append(perim)
            clonlat_list.append(clonlat)
        self.poly_list = poly_list
        self.regional_save(Z,poly_list,area_list,perim_list,clonlat_list,basemap,region)
        return
    
    def binary_diff(self,data,threshm,threshM):
        
        # generating the binary file, get the difference, divide into overprediction and
        # underprediction, apply the closing on both and finally re-map the land
        ice_ext = np.ma.copy(data)
        pack_ice = np.ma.copy(data)

        # OLD method where ice-ice = ocean-ocean = 0
        ice_ext[ice_ext<threshm] = 0
        pack_ice[pack_ice<threshM] = 0
        ice_ext[ice_ext>=threshm] = 1  
        pack_ice[pack_ice>=threshM] = 1  

        # Difference MODEL - OSISAF
        ddata = ice_ext - pack_ice
        # Cut out the nans for both
        thenans = np.isnan(ddata)
        ddata[thenans] = 0
        # Closing operation (see def closing)
        ddata = closing(ddata)
        # Add the under to the over to get total diff.
        return(ddata,ice_ext,pack_ice)
    
    def poly_maker(self,ddata):
        # finding contours from the difference data map
        MIZcont = msr.find_contours(ddata,.5)
        MIZcont = sorted(MIZcont, key=len)
        return(MIZcont)

    def mapper(self,DIFF,Z):
        DN = np.ma.copy(DIFF)
        DN = np.ma.masked_where(np.ma.getmask(Z),DN)
        return(DN)

    def poly_stat(self,cont,pack,ext,miz,X,Y,num,region,basemap):

            # this function find the different polygons

            # (NOTE i and j are inverted -> i = [:,1], j = [:,0])

            # NOTE if a point has non integer i coordinate is going to be a vertical edge,
            # if has non integer j coordinate is going to be a horizontal edge hence the use
            # of different arounds depending on the indexes of the point

            ext = ext
            pack = pack
            in_cont = []
            out_cont = []
            unk_cont = []
            func_vals= []
            func_mod = 0 # value of func_vals at model ice edge
            func_osi = 1 # value of func_vals at OSISAF ice edge
            func_unk = 2 # value of func_vals at other parts of contour
            
            for n,el in enumerate(cont):
                #getting all the neighbours
                around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) 
                around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1]))
                check_cont = 0
                if el[0]/int(el[0]) == 1:
                    for h,v in around1:
                        if pack[h][v] == ext[h][v] == 1:
                        #if pack[h][v] == ext[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif pack[h][v] == ext[h][v] == 0:
                        #elif pack[h][v] == ext[h][v] == 1:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
                else:
                    for h,v in around2:
                        if pack[h][v] == ext[h][v] == 1:
                        #if pack[h][v] == ext[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif pack[h][v] == ext[h][v] == 0:
                        #elif pack[h][v] == ext[h][v] == 1:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
            
            in_cont = np.array(in_cont)
            out_cont = np.array(out_cont)
            unk_cont = np.array(unk_cont)
            f_vals = np.array(func_vals)

            # changes indexes to x and y 
            x = np.ma.copy(cont)
            y = np.ma.copy(cont)
            xvec = range(len(cont[:,0]))
            for n,en in enumerate(xvec):
                en = X[cont[n,0]][cont[n,1]]
                x[n] = en
            yvec = range(len(cont[:,1]))
            for n,en in enumerate(yvec):
                en = Y[cont[n,0]][cont[n,1]]
                y[n] = en
            xy_list = zip(x,y)
            xy_list = np.array(xy_list)
        
            # Centroid - Location in lon/lat of the centroid of the polygon
            DX = 0
            DY = 0
            B = 0
            for n in range(len(cont)-1):
                DX += ((cont[n][0]+cont[n+1][0])*\
                      ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1])))  
                DY += ((cont[n][1]+cont[n+1][1])*\
                      ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1])))  
                B += ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1]))
            A = (0.5)*B 
            CX = (1/float(6*A))*DX
            CY = (1/float(6*A))*DY
            self.centroid = [CX,CY]
            self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
            clonlat = [self.centroid_longitude,self.centroid_latitude]
            lon_list,lat_list=basemap(cont[:,0],cont[:,1],inverse=True)
            poly = []
            for n, en in enumerate(f_vals):
                poly.append([num,lon_list[n],lat_list[n],en])

            a = 0
            x0,y0 = cont[0]
            for [x1,y1] in cont[1:]:
                dx = x1-x0
                dy = y1-y0
                a += 0.5*(y0*dx - x0*dy)
                x0 = x1 
                y0 = y1 
            E_area = abs(a*100)
            
            # Calculating perimeter in xy coordinates (unit = 10km)
            perim = 0
            for n in range(len(cont)-1):
                perim += np.sqrt(pow(cont[n+1][0]-cont[n][0],2)+pow(cont[n+1][1]-cont[n][1],2))
            E_perim = abs(perim*10)

            return(poly,E_area,E_perim,clonlat)
    
    def regional_save(self,data,poly_list,area,perimeter,clonlat,basemap,region):

        # Prep the lists for polygon.txt
        head_poly_list = ('Polygon_number','Longitude','Latitude','Func_value')
    
        if data.shape == (1120,760):
            # REGION
            reg_repo = './outputs/MIZ/OSI/'+str(region)
            if not os.path.exists(reg_repo):
                os.mkdir(reg_repo)
        else:
            # REGION
            reg_repo = './outputs/MIZ/WIM/'+str(region)
            if not os.path.exists(reg_repo):
                os.mkdir(reg_repo)

        # POLYGONS
        polydir = reg_repo+'/polygons'
        if not os.path.exists(polydir):
            os.mkdir(polydir)

        filname = str(polydir)+'/'+str(dadate)+'.txt'
        with open(filname,'w+') as f:
            f.write(str(head_poly_list)+'\n')
            for n,en in enumerate(poly_list):
                for l,el in enumerate(en):
                    string = ' '.join(map(str,el))
                    for item in string:
                        f.write(item)
                    f.write('\n')
            f.close()

        # WIDTHS
        widthsdir = reg_repo+'/widths'
        if not os.path.exists(widthsdir):
            os.mkdir(widthsdir)

        reg_width_list = []
        reg_widths = widths.single_file(filname,basemap)
        for w,wd in enumerate(reg_widths):
            if wd is not None:
                reg_width_list.append(wd.int_widths)
        filname = str(widthsdir)+'/'+str(dadate)+'.txt'
        with open(filname,'w+') as f:
            for l,el in enumerate(reg_width_list):
                string = ' '.join(map(str,el))
                for item in string:
                    f.write(item)
                f.write('\n')
            f.close()

        # AVG WIDTH
        uglyway = open(filname,'r+')
        uglierway = uglyway.readlines()
        w_num = 0
        w_den = 0
        for el in uglierway:
            pw = np.array(map(float,el.split()))
            w_num += sum(pw)
            w_den += len(pw)
        avg_width = w_num/float(w_den)

        # MIZ LOCATION
        num_lon = 0
        num_lat = 0
        den = 0
        for n,en in enumerate(clonlat):
            num_lon += en[0]*area[n]
            num_lat += en[1]*area[n]
            den += area[0]
        mean_lon = num_lon/float(den)
        mean_lat = num_lat/float(den)

        # Regional Lists
        # NOTE list is [ice extent, ice area, miz extent, miz area]
        lst = []
    
        # Which mode?
        if np.nanmax(data) < 1:
            lindata = np.reshape(data,data.size)
            for pt in lindata:
                if pt > .15 and pt < .80:
                    lst.append([1,pt*100,1,pt*100])
                elif pt >= .80:
                    lst.append([1,pt*100,0,0])
        else:
            lindata = np.reshape(data,data.size)
            for pt in lindata:
                if pt > .1 and pt < 300:
                    lst.append([1,pt*100,1,pt*100])
                elif pt == 300:
                    lst.append([1,pt*100,0,0])

        stats = map(sum,zip(*lst))
    
        filname = reg_repo+'/stats.txt'
        with open(filname,'a') as f:
            row1 = [dadate,stats[0],stats[1],stats[2],stats[3],avg_width,mean_lon,mean_lat]
            str1 = ' '.join(map(str,row1))
            f.write(str1+'\n')
            f.close()
        return()

    def stat_chart(self,number,save=False):
    
        vs = cont[number]
    
        a = 0
        x0,y0 = vs[0]
        for [x1,y1] in vs[1:]:
            dx = x1-x0
            dy = y1-y0
            a += 0.5*(y0*dx - x0*dy)
            x0 = x1 
            y0 = y1 
        self.E_area = abs(a*100)
    
        # Calculating perimeter in xy coordinates (unit = 10km)
        perim = 0
        for n in range(len(vs)-1):
            perim += np.sqrt(pow(vs[n+1][0]-vs[n][0],2)+pow(vs[n+1][1]-vs[n][1],2))
        self.E_perim = abs(perim*10)
    
        # class definition
        lcont = len(vs)
        if lcont > 200:
            self.polygon_class = 'H'
        elif lcont > 100 and lcont <= 200:
            self.polygon_class = 'B'
        elif lcont > 30 and lcont <= 100:
            self.polygon_class = 'M'
        elif lcont <= 30:
            self.polygon_class = 'S'
    
        # The Statistical Chart is an output that tries to compress as many informations
        # as possible in a single figure.
        # The figure is:
        # 1) Top left arctic map with studied polygon highlighted
        # 2) Top right a close up of the polygon with contours and min. distances
        # 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
        # 4) Bottom right is a legend for the small image
    
        pname = 'Polygon_'+str(number)
        self.name = pname
        pclass = self.polygon_class
        region = self.region
        print ''
        print 'Statistic Chart for ',pname
        print 'Class and Region ',pclass,region
        DN = self.ip_miz
        ij = self.ij_list
        inside_contour = self.in_ice_edge
        outside_contour = self.out_ice_edge
        unknown_contour = self.unknown_edge
        clon = '%1.2f' %self.centroid_longitude
        clat = '%1.2f' %self.centroid_latitude
        clonlat = '{0}/{1}'.format(clon,clat)
        E_area = self.E_area
        area = '%1.4e' %(E_area) 
        E_perim = self.E_perim
        perim = '%1.4e' %(E_perim)
        
        # Setting up the plot (2x2) and subplots
        fig = plt.figure(figsize=(15,10))
        gs = gridspec.GridSpec(2,2,width_ratios=[2,1],height_ratios=[4,2])
        plt.suptitle(pname+', class '+pclass+', '+region,fontsize=18)
        main = plt.subplot(gs[0,0])
        polyf = plt.subplot(gs[0,1])
        tab = plt.subplot(gs[1,0])
        leg = plt.subplot(gs[1,1])
        tab.set_xticks([])
        leg.set_xticks([])
        tab.set_yticks([])
        leg.set_yticks([])
        tab.set_frame_on(False)
        leg.set_frame_on(False)
        
        # Main image on the top left
        main.imshow(DN[::-1],interpolation='nearest',cmap='winter')
        x1,x2,y1,y2 = np.min(ij[:,1])-10,np.max(ij[:,1])+10,np.min(ij[:,0])-10,np.max(ij[:,0])+10
        main.axvspan(x1,x2,ymin=1-((y1-320)/float(len(DN)-320)),\
            ymax=1-((y2-320)/float(len(DN)-320)),color='red',alpha=0.3)
        main.axis([0,760,0,800])
        
        # Polygon image on the top right
        polyf.imshow(DN,interpolation='nearest',cmap='winter')
        polyf.axis([x1,x2,y2,y1])
        if len(inside_contour) != 0:
            polyf.plot(inside_contour[:,1],inside_contour[:,0],'ro',markersize=4)
        if len(outside_contour) != 0:
            polyf.plot(outside_contour[:,1],outside_contour[:,0],'yo',markersize=4)
        if len(unknown_contour) != 0:
            polyf.plot(unknown_contour[:,1],unknown_contour[:,0],'go',markersize=4)
        
        # Legend on the bottom right
        mc = mlines.Line2D([],[],color='red',marker='o')
        oc = mlines.Line2D([],[],color='yellow',marker='o')
        uc = mlines.Line2D([],[],color='green',marker='o')
        leg.legend([mc,oc,uc],('15% Cont.','80% Cont.','Unknown Cont.'))
                 
        # Statistics text on the bottom left
        txt =  '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
               '2) Area = '+str(area)+' km^2\n'+ \
               '3) Perimeter = '+str(perim)+' km\n'+ \
               '4) Avg. Width = '+'\n'+ \
               '5) Widths SD = '
        tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
        if save:
            outdir = str(out_dir)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            valid_class = outdir+'/'+dadate
            if not os.path.exists(valid_class):
                os.mkdir(valid_class)
            if not os.path.exists(valid_class+'/'+region):
                os.mkdir(valid_class+'/'+region)
            fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
            plt.close()
        else:
            plt.show(False)
        print 'Statistic chart done for '+str(pname)
        print ''
        return()

# From ice_edge_validation
class AODwidth:
    def __init__(self,XMo,YMo,ZM,XO,YO,ZO,basemap,region):
        if ZM.shape != ZO.shape:
            ZMo = reproj_mdl2osi(XMo,YMo,ZM,XO,YO)

        ZMo,ZO = mask_region(ZMo,region),mask_region(ZO,region)

        DN,mdl_bin,osi_bin,ovr_bin,und_bin = self.poly_ovr_und(ZMo,ZO)
        ovr,und = self.poly_maker(DN)
        DN = self.mapper(DN,ZMo,ZO)

        self.DN,self.mdl_bin,self.osi_bin,self.ovr_bin,self.und_bin =\
                DN,mdl_bin,osi_bin,ovr_bin,und_bin

        ovr_list = []
        ovr_area = []
        ovr_perim = []
        ovr_clonlat = []
        und_list = []
        und_area = []
        und_perim = []
        und_clonlat = []

        # Polygon analysis
        for n,el in enumerate(ovr):
            # classification of overestimation polygons
            poly_list,E_area,E_perim,clonlat = self.poly_stat(el,'1',osi_bin,mdl_bin,DN,XO,YO,n,basemap=basemap)
            ovr_list.append(poly_list)
            ovr_area.append(E_area)
            ovr_perim.append(E_perim)
            ovr_clonlat.append(clonlat)
        try:
            n
        except NameError:
            n = -1 
        for n2,el in enumerate(und):
            # classification of negative (n+1 for good enumeration)
            n += 1 
            poly_list,E_area,E_perim,clonlat = self.poly_stat(el,'0',osi_bin,mdl_bin,DN,XO,YO,n,basemap=basemap)
            und_list.append(poly_list)
            und_area.append(E_area)
            und_perim.append(E_perim)
            und_clonlat.append(clonlat)
        self.ovr_list = ovr_list
        self.und_list = und_list
        self.regional_save(ZO,ZMo,ovr_list,und_list,ovr_area,und_area,ovr_clonlat,und_clonlat,region,basemap)
        return
    
    def poly_ovr_und(self,conc_m,conc_o):
    
        mdl = np.ma.copy(conc_m)
        osi = np.ma.copy(conc_o)
        osi[osi.mask==True] = 999
        mdl[osi==999] = 999
        mdl[mdl<.15] = 0
        osi[osi<.15] = 0
        mdl[mdl>=.15] = 1
        osi[osi>=.15] = 1
        aod = osi - mdl
        ovr = np.ma.copy(aod)
        und = np.ma.copy(aod)
        ovr[ovr==1] = 0
        und[und==-1] = 0
        ovr = abs(ovr)
        und = abs(und)
        ovr = closing(ovr)
        und = closing(und)
        self.ovr_poly = ovr
        self.und_poly = und
        
        return(aod,mdl,osi,ovr,und)

    def poly_maker(self,ddata):
        # finding contours from the difference data map
        ovr = msr.find_contours(ddata,.5)
        ovr = sorted(ovr, key=len)
        ovr = ovr[::-1]
        und = msr.find_contours(ddata,-.5)
        und = sorted(und, key=len)
        und = und[::-1]
        return(ovr,und)
    
    def mapper(self,DN,ZMo,ZO):
        DN = np.ma.copy(DN)
        mask1 = np.ma.getmask(ZMo)
        mask2 = np.ma.mask_or(np.ma.getmask(ZO),mask1)
        DN = np.ma.masked_where(mask2,DN)
        return(DN)

    def poly_stat(self,cont,poly_status,osi_bin,mdl_bin,aod,X,Y,number,basemap=None):

        # NOTE if a point has non integer i coordinate is going to be a vertical edge,
        # if has non integer j coordinate is going to be a horizontal edge hence the use
        # of different arounds depending on the indexes of the point

        # NOTE if the polygon is an OVERESTIMATION the contour finding is inverted

        in_cont = []
        out_cont = []
        unk_cont = []
        func_vals= []
        func_mod = 0 # value of func_vals at model ice edge
        func_osi = 1 # value of func_vals at OSISAF ice edge
        func_unk = 2 # value of func_vals at other parts of contour
        
        if poly_status == 0:
            for n,el in enumerate(cont):
                #getting all the neighbours
                around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
                around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
                check_cont = 0
                if el[0]/int(el[0]) == 1:
                    for h,v in around1:
                        if osi_bin[h][v] == mdl_bin[h][v] == 1:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif osi_bin[h][v] == mdl_bin[h][v] == 0:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
                else:
                    for h,v in around2:
                        if osi_bin[h][v] == mdl_bin[h][v] == 1:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif osi_bin[h][v] == mdl_bin[h][v] == 0:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
        else:
            for n,el in enumerate(cont):
                #getting all the neighbours
                around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
                around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
                check_cont = 0
                if el[0]/int(el[0]) == 1:
                    for h,v in around1:
                        if osi_bin[h][v] == mdl_bin[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif osi_bin[h][v] == mdl_bin[h][v] == 1:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
                else:
                    for h,v in around2:
                        if osi_bin[h][v] == mdl_bin[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif osi_bin[h][v] == mdl_bin[h][v] == 1:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
        
        in_cont = np.array(in_cont)
        out_cont = np.array(out_cont)
        unk_cont = np.array(unk_cont)
        f_vals = np.array(func_vals)

        a = 0
        x0,y0 = cont[0]
        for [x1,y1] in cont[1:]:
            dx = x1-x0
            dy = y1-y0
            a += 0.5*(y0*dx - x0*dy)
            x0 = x1 
            y0 = y1 
        E_area = abs(a*100)
        # Calculating perimeter in xy coordinates (unit = 10km)
        perim = 0
        for n in range(len(cont)-1):
            perim += np.sqrt(pow(cont[n+1][0]-cont[n][0],2)+pow(cont[n+1][1]-cont[n][1],2))
        E_perim = abs(perim*10)

        # changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
        x = []
        y = []
        xvec = range(len(cont[:,0]))
        for n,en in enumerate(xvec):
            en = X[cont[n,0]][cont[n,1]]
            x.append(en)
        yvec = range(len(cont[:,1]))
        for n,en in enumerate(yvec):
            en = Y[cont[n,0]][cont[n,1]]
            y.append(en)
        xy_list = zip(x,y)
        xy_list = np.array(xy_list)
    
        ## Centroid - Location in lon/lat of the centroid of the polygon
        #DX = 0
        #DY = 0
        #B = 0
        #for n in range(len(xy_list)-1):
        #    DX += ((xy_list[n][0]+xy_list[n+1][0])*\
        #          ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
        #    DY += ((xy_list[n][1]+xy_list[n+1][1])*\
        #          ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
        #    B += ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1]))
        #A = (0.5)*B 
        #CX = (1/float(6*A))*DX
        #CY = (1/float(6*A))*DY
        #centroid_longitude,centroid_latitude = basemap(CX,CY,inverse=True)
        #clonlat = [centroid_longitude,centroid_latitude]
        lon,lat=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
        clon = np.mean(lon)
        clat = np.mean(lat)
        clonlat = [clon,clat]

        poly_list = []
        for n, en in enumerate(f_vals):
            poly_list.append([number,lon[n],lat[n],en])

        return(poly_list,E_area,E_perim,clonlat)

    def regional_save(self,ZO,ZMo,ovr_list,und_list,ovr_area,und_area,ovr_clonlat,und_clonlat,region,basemap):
    
        # NOTE data1 is MODEL, data2 is OSISAF
        # This gives the difference as:
        # -1 = OSI PRES, MDL MISS -> UNDERESTIMATION
        #  0 = OSI PRES, MDL PRES (OR OSI MISS, MDL MISS) -> HIT (NO-HIT)
        # +1 = OSI MISS, MDL PRES -> OVERESTIMATION
        
        # NOTE outputs are difference, overpredictions and underpredictions (CLOSED & MAPPED)
        
        # generating the binary file, get the difference, divide into overprediction and
        # underprediction, apply the closing on both and finally re-map the land
        ndata1 = np.ma.copy(ZMo)
        ndata2 = np.ma.copy(ZO)
        
        # NEW method where ice-ice = 2 ; ocean-ocean = -2
        ndata1[ndata1>=.15] = 1
        ndata1[ndata1<.15] = 2
        ndata2[ndata2<.15] = 0
        ndata2[ndata2>=.15] = 3
        
        # Difference MODEL - OSISAF
        ddata = ndata1 - ndata2
        # Cut out the nans for both
        thenans = np.isnan(ddata)
        ddata[thenans] = 0
        ddata[ndata1.mask==True] = 0
        ddata[ndata2.mask==True] = 0
        # Work with over and under (clean them up)
        ovr_data = np.copy(ddata)
        und_data = np.copy(ddata)
        ovr_data[ovr_data==-1] = 0
        und_data[und_data==1] = 0
        und_data = abs(und_data)
        # Closing operation (see def closing)
        o_data = closing(ovr_data)
        u_data = closing(und_data)
        # Add the under to the over to get total diff.
        o_data[u_data==1] = -1
        ddata = np.copy(o_data)
    
        tot_hits = ddata[ddata==2].size
        tot_ovr = ddata[ddata==1].size
        tot_und = ddata[ddata==-1].size
        tot_neg_hits = ddata[ddata==-2].size
        tot_unk = ddata[ddata==0].size
        total = ddata.size - tot_unk
    
        # Prep the lists for polygon.txt
        head_poly_list = 'Polygon_number','Longitude','Latitude','Func_value'

        ovr_width_list = []
        und_width_list = []
        
        # Output directory
        outdir = './outputs/AOD'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        reg_repo = outdir+'/'+str(region)
        if not os.path.exists(reg_repo):
            os.mkdir(reg_repo)
    
        # Writing the overestimation polygons list
        if len(ovr_list) == 0:
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_polygons.txt'
           with open(filname,'w+') as f:
               f.write(str(head_poly_list)+'\n')
               string = '0','0','0','0'
               f.write(item)
               f.write('\n')
               f.close()
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_widths.txt'
           with open(filname,'w+') as f:
               string = '0','0','0','0'
               f.write(item)
               f.write('\n')
               f.close()
           mean_ovr_lon = 0
           mean_ovr_lat = 0
           avg_ovr = 0
        else:
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_polygons.txt'
           with open(filname,'w+') as f:
               f.write(str(head_poly_list)+'\n')
               for n,en in enumerate(ovr_list):
                   for l,el in enumerate(en):
                       string = ' '.join(map(str,el))
                       for item in string:
                           f.write(item)
                       f.write('\n')
               f.close()
    
           # Width calculation for OVERESTIMATION
           ovr_widths = widths.single_file(filname,basemap)
           for w,wd in enumerate(ovr_widths):
               if wd is not None:
                   ovr_width_list.append(wd.int_widths)
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_widths.txt'
           with open(filname,'w+') as f:
               for l,el in enumerate(ovr_width_list):
                   string = ' '.join(map(str,el))
                   for item in string:
                       f.write(item)
                   f.write('\n')
               f.close()
           tot = 0
           n = 0
           for en in ovr_width_list:
               for el in en:
                   if el != 0:
                       tot += el
                       n += 1
           avg_ovr = tot/float(n)

           # OVER LOCATION
           num_lon = 0
           num_lat = 0
           den = 0
           for n,en in enumerate(ovr_clonlat):
               num_lon += en[0]*ovr_area[n]
               num_lat += en[1]*ovr_area[n]
               den += ovr_area[n]
           mean_ovr_lon = num_lon/float(den)
           mean_ovr_lat = num_lat/float(den)
          
        # Writing the underestimation polygons list
        if len(und_list) == 0:
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_polygons.txt'
           with open(filname,'w+') as f:
               f.write(str(head_poly_list)+'\n')
               string = '0','0','0','0'
               f.write(item)
               f.write('\n')
               f.close()
           filname = str(reg_repo)+'/'+str(dadate)+'_ovr_widths.txt'
           with open(filname,'w+') as f:
               string = '0','0','0','0'
               f.write(item)
               f.write('\n')
               f.close()
           mean_und_lat = 0
           mean_und_lon = 0
           avg_und = 0
        else:
           filname = str(reg_repo)+'/'+str(dadate)+'_und_polygons.txt'
           with open(filname,'w+') as f:
               f.write(str(head_poly_list)+'\n')
               for n,en in enumerate(und_list):
                   for l,el in enumerate(en):
                       string = ' '.join(map(str,el))
                       for item in string:
                           f.write(item)
                       f.write('\n')
               f.close()
    
           # Width calculation for UNDERESTIMATION
           und_widths = widths.single_file(filname,basemap)
           for w,wd in enumerate(und_widths):
               if wd is not None:
                   und_width_list.append(wd.int_widths)   
           filname = str(reg_repo)+'/'+str(dadate)+'_und_widths.txt'
           with open(filname,'w+') as f:
               for l,el in enumerate(und_width_list):
                   string = ' '.join(map(str,el))
                   for item in string:
                       f.write(item)
                   f.write('\n')
               f.close()
           tot = 0
           n = 0
           for en in und_width_list:
               for el in en:
                   if el != 0:
                       tot += el
                       n += 1
           avg_und = tot/float(n)
   
           # UNDER LOCATION
           num_lon = 0
           num_lat = 0
           den = 0
           for n,en in enumerate(und_clonlat):
               num_lon += en[0]*und_area[n]
               num_lat += en[1]*und_area[n]
               den += und_area[n]
           mean_und_lon = num_lon/float(den)
           mean_und_lat = num_lat/float(den)
   
        filname = str(reg_repo)+'/AOD_stats.txt'
        with open(filname,'a') as f:
            row1 = [total,tot_hits,tot_ovr,tot_und,tot_neg_hits,tot_unk,\
                    avg_ovr,avg_und,mean_ovr_lon,mean_ovr_lat,mean_und_lon,mean_und_lat]
            str1 = ' '.join(map(str,row1))
            f.write(str1)
            f.write('\n')
            f.close()
    
        return()

############################################################################
# UTILITIES
############################################################################
def reproj_mdl2osi(Xold,Yold,Zold,Xnew,Ynew):

    # getting ready for reprojection
    XI = Xold.reshape(Xold.size)
    YI = Yold.reshape(Yold.size)
    ZI = Zold.reshape(Zold.size)
    C = [XI,YI]
    C = np.array(C)
    C = C.T

    # Interpolation can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
    ZNn = grd(C,ZI,(Xnew,Ynew),method='nearest')
    
    return(ZNn)
# Functions that linearize arrays and viceversa
def linearizer(arr):
    linarr = np.reshape(arr,arr.size)
    return(linarr)
def arrayizer(lin,arr=None):
    if arr is not None:
        narr = np.reshape(lin,arr.shape)
    else:
        arr = str(lin)
        arr = arr[3::]
        narr = np.reshape(lin,arr.shape)
    return(narr)
# Function on basemap creations
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
# Function that masks regions (only BS1 and FR1) will be used
def mask_region(data,region,mreg=None,oreg=None,mmask=None,omask=None,mask=False):
    if mreg is None or oreg is None:
        mdl_reg = np.load('./geo_info/bar_gre_bal_reg.npy')
        osi_reg = np.load('./geo_info/bar_gre_bal_reg_OSI.npy')
        mdl_mask = np.load('./geo_info/bar_gre_bal_reg_mask.npy')
        osi_mask = np.load('./geo_info/bar_gre_bal_reg_mask_OSI.npy')
    else:
        mdl_reg = mreg
        osi_reg = oreg
        mdl_mask = mmask
        osi_mask = omask
    if region == 'gre':
        mdl_reg[mdl_reg!=1] = 0
        osi_reg[osi_reg!=1] = 0
    elif region == 'bar':
        mdl_reg[mdl_reg!=2] = 0
        mdl_reg[mdl_reg!=0] = 1
        osi_reg[osi_reg!=2] = 0
        osi_reg[osi_reg!=0] = 1
    elif region == 'bal':
        mdl_reg[mdl_reg!=3] = 0
        mdl_reg[mdl_reg!=0] = 1
        osi_reg[osi_reg!=3] = 0
        osi_reg[osi_reg!=0] = 1
    else:
        print('no region selected, three masks applied: 1=gre,2=bar,3=bal')
    if data.shape == mdl_reg.shape:
        ndata = data*mdl_reg
    else:
        ndata = data*osi_reg
    if mask:
        if data.mask.shape == mdl_mask.shape and data.shape != mdl_reg.shape:
            ndata = np.ma.mask_or(data.mask,mdl_mask)
        else:
            ndata = np.ma.mask_or(data.mask,osi_mask)
    return(ndata)
# Function with visual operation of closing
def closing(data):
    # apply closing to avoid small polynias and clean up a little
    kernel = np.ones((3,3),np.uint8)
    data_cl = morph.closing(data,kernel)
    return(data_cl)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Beginning of the script
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# MISSION and DATE
dadate = sys.argv[1] 
region = sys.argv[2] 

# Start the count and analyse the datasets
time0 = time.time()
print('Working on '+str(dadate))
print('')

# Let's go, go, GO!
time1 = time.time()
print('Basemap creation...')
#bmm = basemap_creator('TP4')
bmo = basemap_creator('OSI')
#if region == 'bar':
#   bmb = basemap_creator('BS1')
#   print('Barents!')
#elif region == 'gre':
#   bmf = basemap_creator('FR1')
#   print('Greenland!')
elapsedtime = time.time() - time1
print 'Basemap created in ',elapsedtime
print('')

mrg = np.load('./geo_info/bar_gre_bal_reg.npy')
org = np.load('./geo_info/bar_gre_bal_reg_OSI.npy')
mms = np.load('./geo_info/bar_gre_bal_reg_mask.npy')
oms = np.load('./geo_info/bar_gre_bal_reg_mask_OSI.npy')

print('Reading datasets...')
data = reader(dadate,bmo,bmo)
print('DONE')
print('')

if region == 'bar':
   # Barents
   print('Studying the BARENTS region...')
   print('')
   
   XMo,YMo,latM,lonM,XO,YO,latO,lonO = data.XM,data.YM,data.latM,data.lonM,\
           data.XO,data.YO,data.latO,data.lonO
   
   ZM,ZO = data.ZM,data.ZO

   print('Analysing AOD...')
   bar_AOD = AODwidth(XMo,YMo,ZM,XO,YO,ZO,bmo,'bar')
   print('DONE')
   print('')

elif region == 'gre':
   print('Studying the GREENLAND region...')
   print('')
   
   XMo,YMo,latM,lonM,XO,YO,latO,lonO = data.XM,data.YM,data.latM,data.lonM,\
           data.XO,data.YO,data.latO,data.lonO

   ZM,ZO = data.ZM,data.ZO
   
   print('Analysing AOD...')
   gre_AOD = AODwidth(XMo,YMo,ZM,XO,YO,ZO,bmo,'gre')
   print('DONE')
   print('')
else:
   print('Studying the BARENTS region...')
   print('')
   
   XMo,YMo,latM,lonM,XO,YO,latO,lonO = data.XM,data.YM,data.latM,data.lonM,\
           data.XO,data.YO,data.latO,data.lonO
   
   ZM,ZO = data.ZM,data.ZO

   print('Analysing AOD...')
   bar_AOD = AODwidth(XMo,YMo,ZM,XO,YO,ZO,bmo,'bar')
   print('DONE')
   print('')

   print('Studying the GREENLAND region...')
   print('')
   
   print('Analysing AOD...')
   gre_AOD = AODwidth(XMo,YMo,ZM,XO,YO,ZO,bmo,'gre')
   print('DONE')
   print('')


elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime
