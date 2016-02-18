############################################################################
# THIS script will:
# 1) Read in ice_only and waves_ice
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
#matplotlib.use('Agg')
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

# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#CLASSES
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

############################################################################
# WIM thermo and stress analyzer
class WIM:
    def __init__(self,dadate,obasemap,mbasemap):
        self.year = dadate[:4]
        self.month = dadate[4:6]
        self.day = dadate[6:8]
        gigio = datetime.datetime(int(float(self.year)),int(float(self.month)),int(float(self.day)),0,0)
        gigio = gigio.strftime('%j')
        gigio = int(float(gigio))
        self.day_j = gigio
        self.jday_start(gigio)
        self.WIM_timeframes = self.time_frame(gigio)

        # Useful variables
        self.grid_cell_area = 1
        self.delta_t = 400

        #Reading the datasets
        self._read_osi_(dadate,obasemap)
        #self._read_osi_2(dadate,obasemap)
        self._read_mdl_ice_only(dadate,mbasemap)
        self._read_mdl_waves_ice(dadate,mbasemap)

        self.poly_over_under(self.ZCI,self.ZCO)

        # Recalculation of fluxes and growths
        #self.OSI_ice_growth(self.ZCO,self.ZCO2)

        self.stress_h,self.stress_w = self.stresses(self.TXI,self.TYI),\
                self.stresses(self.TXW,self.TYW)
        self.Dmean,self.D2mean,self.Slat,self.Sbot,self.Alpha,self.Beta =\
        self.floes_analysis(self.FOI,self.TOI,self.ZDW)
        self.qlat,self.qvrt = self.heat_flux(self.QA,self.QC,self.QO,self.FOI)
        self.glat_HYCOM,self.gvrt_HYCOM,self.glat_WIM,self.gvrt_WIM,self.CNW,self.TNW =\
        self.WIM_ice_growth(self.TOI,self.TDI,self.FOI,self.FDI,self.HSD)

        # Ice conc BIAS
        #tt1 = self.reproj_wim2osi(self.ZCI)
        #self.BIAS = self.ZCO - tt1

        return
    
    def jday_start(self,jday):
        if jday < 60 or jday > 274:
            print 'ERROR: dataset not in the melting season'
        elif jday > 60 and jday < 68:
            jdaystart = 60
        elif jday > 67 and jday < 75:
            jdaystart = 67
        elif jday > 74 and jday < 82:
            jdaystart = 74
        elif jday > 81 and jday < 89:
            jdaystart = 81
        elif jday > 88 and jday < 96:
            jdaystart = 88
        elif jday > 95 and jday < 103:
            jdaystart = 95
        elif jday > 102 and jday < 110:
            jdaystart = 102
        elif jday > 109 and jday < 117:
            jdaystart = 109
        elif jday > 116 and jday < 124:
            jdaystart = 116
        elif jday > 123 and jday < 131:
            jdaystart = 123
        elif jday > 130 and jday < 138:
            jdaystart = 130
        elif jday > 137 and jday < 145:
            jdaystart = 137
        elif jday > 144 and jday < 159:
            jdaystart = 144
        elif jday > 158 and jday < 166:
            jdaystart = 158
        elif jday > 165 and jday < 173:
            jdaystart = 165
        elif jday > 172 and jday < 180:
            jdaystart = 172
        elif jday > 179 and jday < 187:
            jdaystart = 179
        elif jday > 186 and jday < 194:
            jdaystart = 186
        elif jday > 193 and jday < 201:
            jdaystart = 193
        elif jday > 200 and jday < 208:
            jdaystart = 200
        elif jday > 207 and jday < 215:
            jdaystart = 207
        elif jday > 214 and jday < 222:
            jdaystart = 214
        elif jday > 221 and jday < 236:
            jdaystart = 221
        elif jday > 235 and jday < 243:
            jdaystart = 235
        elif jday > 242 and jday < 250:
            jdaystart = 242
        elif jday > 249 and jday < 257:
            jdaystart = 249
        elif jday > 256 and jday < 264:
            jdaystart = 256
        elif jday > 263 and jday < 271:
            jdaystart = 263
        elif jday > 270 and jday < 275:
            jdaystart = 270
        self.jdaystart = jdaystart
        return()

    def time_frame(self,jday):
        # Time frames for SWARP.nc files
        dd = jday - self.jdaystart
        tf_12 = 4 + ((dd -1)*8)
        if tf_12 == 52:
            self.jdaystart2 = int(self.jdaystart) + 7
            tf_18 = 54
            tf_24 = 0
            tf_30 = 2
            tf_36 = 4
        else:
            self.jdaystart2 = 0
            tf_18 = tf_12 + 2
            tf_24 = tf_12 + 4
            tf_30 = tf_12 + 6
            tf_36 = tf_12 + 8
        return(tf_12,tf_18,tf_24,tf_30,tf_36)

    def _read_osi_(self,dadate,basemap):
        day = self.day
        month = self.month
        year = self.year
        # Read in OSI_SAF file
        #outdir = '/work/timill/giacomo/osisaf_repo'
        outdir = '/home/charlie/Documents/ncdata'
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
        self.ZCO = ZO

        print('OSI-SAF: lonO,latO,XO,YO,ZCO')
        
        return() 
   
    def _read_osi_2(self,dadate,basemap):

        day = self.day
        month = self.month
        year = self.year
        dadate2 = datetime.date(int(year),int(month),int(day))+datetime.timedelta(1)
        dadate2 = dadate2.strftime('%Y%m%d')
        # Read in OSI_SAF file
        #outdir = '/work/timill/giacomo/osisaf_repo'
        outdir = '/home/charlie/Documents/ncdata'
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

    def _read_mdl_ice_only(self,dadate,basemap):
        
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        jdaystart = self.jdaystart

        # TIME INDEX
        tidx = self.WIM_timeframes[0]

        # Read TP4arch
        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/'+\
                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
        outdir = '/home/charlie/Documents/ncdata'
        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
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
        self.lonMI,self.latMI,self.XI,self.YI = lonM,latM,X,Y
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

        self.ZCI,self.ZTI,self.HSI,self.QT,self.QC,self.QO,self.QA,self.FOI,\
                self.TOI,self.HSO,self.HSD,self.FDI,self.TDI,self.TXI,self.TYI =\
                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
                hsnow_diff,conc_diff,thic_diff,taux,tauy
        
        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')

        return()

    def _read_mdl_waves_ice(self,dadate,basemap):
    
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        jdaystart = self.jdaystart
        jday = jdaystart + 1
        gigio = datetime.datetime.strptime( str(year)+str(jday), '%Y%j')
        stdate = gigio.strftime('%Y%m%d')

        # TIME INDEX
        tidx = self.WIM_timeframes[0]

        # Read TP4arch_wav
        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/\
                waves_ice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
        outdir = '/home/charlie/Documents/ncdata'
        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+dadate+'T000000Z.nc'
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
        self.lonMW,self.latMW,self.XW,self.YW = lonM,latM,X,Y
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
    
    def poly_over_under(self,conc_m,conc_o):
    
        mdl = np.ma.copy(conc_m)
        osi = np.ma.copy(conc_o)
        if mdl.size != osi.size:
            mdl = self.reproj_wim2osi(mdl)
        osi[osi.mask==True] = 999
        mdl[osi==999] = 999
        mdl[mdl<.15] = 0
        osi[osi<.15] = 0
        mdl[mdl>=.15] = 1
        osi[osi>=.15] = 1
        polys = osi - mdl
        over = np.ma.copy(polys)
        under = np.ma.copy(polys)
        over[over==1] = 0
        under[under==-1] = 0
        over = abs(over)
        under = abs(under)
        over = self.closing(over)
        under = self.closing(under)
        self.over_poly = over
        self.under_poly = under
        
        return(over,under)
   
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
            print('Calculating FSD: Alpha, Beta, Slat, Sbot')
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
                elif en <= 20 and linconc_old[n] >= .15:
                    en = 21
                    A = alpha/(pow(dmin,-alpha)-pow(en,-alpha))
                    linDmean[n] = A/(alpha-1)*(pow(dmin,1-alpha)-pow(en,1-alpha))
                    linD2mean[n] = A/(alpha-2)*(pow(dmin,2-alpha)-pow(en,2-alpha))
                    linBeta[n] = 4*linthic_old[n]*linDmean[n]/linD2mean[n]
                    linAlpha[n] = linBeta[n]/(1+linBeta[n])
                    linSbot[n] = np.ma.copy(linconc_old[n])*Asq
                    linSlat[n] = linBeta[n]*linSbot[n]
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
            print('Calculating FSD: Alpha, Beta, Slat, Sbot')
            
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
        
        print('Calculating heat fluxes: qlat,qvrt')

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
        
        print('Calculating ice growths: Glat/vrt_HYCOM,Glat/vrt_WIM')

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

    def stresses(self,tx,ty):
        tau = np.sqrt(tx**2+ty**2)
        return(tau)

############################################################################
# UTILITIES
############################################################################
# Function that reprojects model into observational grid
    def reproj_wim2osi(self,Zold):
    
        # getting the correct Xold and Yold from OSI basemap
        bm = basemap_creator('OSI')
        Xold,Yold = bm(self.lonMI,self.latMI,inverse=False)
        Xnew,Ynew = self.XO,self.YO
    
        # getting ready for reprojection
        X3 = Xold.reshape(Xold.size)
        Y3 = Yold.reshape(Yold.size)
        Z3 = Zold.reshape(Zold.size)
        C = [X3,Y3]
        C = np.array(C)
        C = C.T
    
        # Interpolation can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
        ZNn = grd(C,Z3,(Xnew,Ynew),method='nearest')
        
        return(ZNn)
# Functions that linearize arrays and viceversa
    def linearizer(self,arr):
        linarr = np.reshape(arr,arr.size)
        return(linarr)
    def arrayizer(self,lin,arr=None):
        if arr is not None:
            narr = np.reshape(lin,arr.shape)
        else:
            arr = str(lin)
            arr = arr[3::]
            narr = np.reshape(lin,arr.shape)
        return(narr)
# Reproject,linearize and plot
    def tbp(self,arr1,arr2):
        if arr1.size != arr2.size:
            arr1_good = self.reproj_wim2osi(arr1)
        else:
            arr1_good = arr1
        linarr1 = np.reshape(arr1_good,arr1_good.size)
        linarr2 = np.reshape(arr2,arr2.size)
        tbp1 = []
        tbp2 = []
        for n,en in enumerate(linarr1):
            if linarr1.mask[n] != True or linarr1.data[n] != 0 or linarr2.mask[n] != True or linarr2.data[n] != 0:
                tbp1.append(linarr1[n])
                tbp2.append(linarr2[n])
        return(tbp1,tbp2)
# Closing, image processing for excluding little polygons
    def closing(self,DN):
        # apply closing to avoid small polynias and clean up a little
        kernel = np.ones((3,3),np.uint8)
        DN_cl = morph.closing(DN,kernel)
        return(DN_cl)

##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# FUNCTIONS
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Baffin Bay reader
def baffin_bay():
    fil = ('./geo_info/baffin_bay.txt')
    bbay = open(fil)
    bblon = [] 
    bblat = [] 
    bblone = [] 
    bblate = [] 
    bblonw = [] 
    bblatw = [] 
    for n,en in enumerate(bbay):
        nen = en.split(';')
        lon = float(nen[1])
        lat = float(nen[2])
        bblon.append(lon)
        bblat.append(lat)
    for n,en in enumerate(bblon):
        if n < 1266:
            bblone.append(bblon[n])
            bblate.append(bblat[n])
        else:
            bblonw.append(bblon[n])
            bblatw.append(bblat[n])
    return(bblone,bblate,bblonw,bblatw)

############################################################################
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

############################################################################
# Function that gives regional info about hits, over, under and neg-hits
# As inputs model data and mode (ic, fsd)
def regional_stats(data,widths,lons,lats,mode):
    
    # widths
    bar_avg_widths = widths[0]
    les_avg_widths = widths[1]
    ncb_avg_widths = widths[2]
    lab_avg_widths = widths[3]
    gre_avg_widths = widths[4]
 
    # lons
    bar_avg_lons = lons[0]
    les_avg_lons = lons[1]
    ncb_avg_lons = lons[2]
    lab_avg_lons = lons[3]
    gre_avg_lons = lons[4]
 
    # lats
    bar_avg_lats = lats[0]
    les_avg_lats = lats[1]
    ncb_avg_lats = lats[2]
    lab_avg_lats = lats[3]
    gre_avg_lats = lats[4]
 
    # we want sea ice extent(area), MIZ extent(area), navigability area?

    # Regional file again, to be sure
    regions = np.load('./geo_info/regions.npy')

    # Regional Lists
    # NOTE list is [ice extent, ice area, miz extent, miz area]
    bar_lst = []
    gre_lst = []
    lab_lst = []
    ncb_lst = []
    les_lst = []

    # Which mode?
    if mode == 'ic':
        # Going through every point to assess location
        for i in range(len(data[:,0])):
            for j in range(len(data[0,:])):
                if regions[i,j] == 1:
                    if data[i,j] > .15 and data[i,j] < .80:
                        bar_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                    elif data[i,j] >= .80:
                        bar_lst.append([1,data[i,j]*100,0,0])
                elif regions[i,j] == 2:
                    if data[i,j] > .15 and data[i,j] < .80:
                        les_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                    elif data[i,j] >= .80:
                        les_lst.append([1,data[i,j]*100,0,0])
                elif regions[i,j] == 3:
                    if data[i,j] > .15 and data[i,j] < .80:
                        ncb_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                    elif data[i,j] >= .80:
                        ncb_lst.append([1,data[i,j]*100,0,0])
                elif regions[i,j] == 4:
                    if data[i,j] > .15 and data[i,j] < .80:
                        lab_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                    elif data[i,j] >= .80:
                        lab_lst.append([1,data[i,j]*100,0,0])
                elif regions[i,j] == 5:
                    if data[i,j] > .15 and data[i,j] < .80:
                        gre_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                    elif data[i,j] >= .80:
                        gre_lst.append([1,data[i,j]*100,0,0])
    elif mode == 'fsd':
       # Going through every point to assess location
       for i in range(len(data[:,0])):
           for j in range(len(data[0,:])):
               if regions[i,j] == 1:
                   if data[i,j] > .1 and data[i,j] < 300:
                       bar_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                   elif data[i,j] == 300:
                       bar_lst.append([1,data[i,j]*100,0,0])
               elif regions[i,j] == 2:
                   if data[i,j] > .1 and data[i,j] < 300:
                       les_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                   elif data[i,j] == 300:
                       les_lst.append([1,data[i,j]*100,0,0])
               elif regions[i,j] == 3:
                   if data[i,j] > .1 and data[i,j] < 300:
                       ncb_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                   elif data[i,j] == 300:
                       ncb_lst.append([1,data[i,j]*100,0,0])
               elif regions[i,j] == 4:
                   if data[i,j] > .1 and data[i,j] < 300:
                       lab_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                   elif data[i,j] == 300:
                       lab_lst.append([1,data[i,j]*100,0,0])
               elif regions[i,j] == 5:
                   if data[i,j] > .1 and data[i,j] < 300:
                       gre_lst.append([1,data[i,j]*100,1,data[i,j]*100])
                   elif data[i,j] == 300:
                       gre_lst.append([1,data[i,j]*100,0,0])
  
    #bar_lst = np.array(bar_lst)
    #les_lst = np.array(les_lst)
    #ncb_lst = np.array(ncb_lst)
    #lab_lst = np.array(lab_lst)
    #gre_lst = np.array(gre_lst)

    bar_stats = map(sum,zip(*bar_lst))
    les_stats = map(sum,zip(*les_lst))
    ncb_stats = map(sum,zip(*ncb_lst))
    lab_stats = map(sum,zip(*lab_lst))
    gre_stats = map(sum,zip(*gre_lst))

    filname = str(out_dir)+'/bar_stats.txt'
    with open(filname,'a') as f:
        row1 = [dadate,bar_stats[0],bar_stats[1],bar_stats[2],bar_stats[3],bar_avg_widths,bar_avg_lons,bar_avg_lats]
        str1 = ' '.join(map(str,row1))
        f.write(str1+'\n')
        f.close()

    filname = str(out_dir)+'/les_stats.txt'
    with open(filname,'a') as f:
        row1 = [dadate,les_stats[0],les_stats[1],les_stats[2],les_stats[3],les_avg_widths,les_avg_lons,les_avg_lats]
        str1 = ' '.join(map(str,row1))
        f.write(str1+'\n')
        f.close()

    filname = str(out_dir)+'/ncb_stats.txt'
    with open(filname,'a') as f:
        row1 = [dadate,ncb_stats[0],ncb_stats[1],ncb_stats[2],ncb_stats[3],ncb_avg_widths,ncb_avg_lons,ncb_avg_lats]
        str1 = ' '.join(map(str,row1))
        f.write(str1+'\n')
        f.close()

    filname = str(out_dir)+'/lab_stats.txt'
    with open(filname,'a') as f:
        row1 = [dadate,lab_stats[0],lab_stats[1],lab_stats[2],lab_stats[3],lab_avg_widths,lab_avg_lons,lab_avg_lats]
        str1 = ' '.join(map(str,row1))
        f.write(str1+'\n')
        f.close()

    filname = str(out_dir)+'/gre_stats.txt'
    with open(filname,'a') as f:
        row1 = [dadate,gre_stats[0],gre_stats[1],gre_stats[2],gre_stats[3],gre_avg_widths,gre_avg_lons,gre_avg_lats]
        str1 = ' '.join(map(str,row1))
        f.write(str1+'\n')
        f.close()

    return()

############################################################################
# Function that gives statistics from widths
def widths_stats(lst):
    tot = 0
    n = 0
    for en in lst:
        for el in en:
            if el != 0:
                tot += el
                n += 1
    avg = tot/float(n)
    return(avg)


############################################################################
# Function that masks everything but the MIZ - data is the array
# that will be masked, conc  is the reference ice_concentration
def mask_as_MIZ(data,dmax):
    mask1 = np.ma.make_mask(dmax<20)
    mask2 = np.ma.make_mask(dmax>200)
    mask3 = np.ma.mask_or(mask1,mask2)
    miz_data = np.ma.copy(data)
    miz_data = np.ma.masked_where(mask3,miz_data)
    return(miz_data)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Beginning of the script
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# MISSION and DATE
dadate = sys.argv[1] 
out_dir = sys.argv[2]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Start the count and analyse the datasets
time0 = time.time()

# Call in the basemap
# TODO still having problems with model's basemap
bmm = basemap_creator('TP4')
bmo = basemap_creator('OSI')
#bmb = basemap_creator('BS1')
#bmf = basemap_creator('FR1')

# Call in the regions
# LABELS -> [BAR = 1, LES = 2, NCB = 3, LAB = 4, GRE = 5]
regions = np.load('./geo_info/regions.npy')

data = WIM(dadate,bmo,bmm)
np.save(str(out_dir)+'/'+str(dadate)+'_glat_H',data.reproj_wim2osi(data.glat_HYCOM))
np.save(str(out_dir)+'/'+str(dadate)+'_gvrt_H',data.reproj_wim2osi(data.gvrt_HYCOM))
np.save(str(out_dir)+'/'+str(dadate)+'_glat_W',data.reproj_wim2osi(data.glat_WIM))
np.save(str(out_dir)+'/'+str(dadate)+'_gvrt_W',data.reproj_wim2osi(data.gvrt_WIM))
np.save(str(out_dir)+'/'+str(dadate)+'_stress_H',data.reproj_wim2osi(data.stress_h))
np.save(str(out_dir)+'/'+str(dadate)+'_stress_W',data.reproj_wim2osi(data.stress_w))
np.save(str(out_dir)+'/'+str(dadate)+'_over',data.reproj_wim2osi(data.over_poly))
np.save(str(out_dir)+'/'+str(dadate)+'_under',data.reproj_wim2osi(data.under_poly))


elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime


# OLD OLD OLD OLD OLD OLD OLD OLD
#
#    def _read_mdl_ice_only_1(self,dadate,basemap):
#        
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[1]
#
#        # Read TP4arch
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/'+\
#                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        shsnow = 'hsnow'
#        sqtot = 'qtot'
#        sqcool = 'flx_cool'
#        sqother = 'flx_itop'
#        sqatm = 'flx_ow'
#        sconc_old = 'fi_old'
#        sthic_old = 'hi_old'
#        shsnow_old = 'hs_old'
#        shsnow_diff = 'hs_int'
#        sconc_diff = 'dfice'
#        sthic_diff = 'dhice'
#        staux = 'taux'
#        stauy = 'tauy'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMI,self.latMI,self.XI,self.YI = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        hsnow = Mrdg.nc_get_var(ncfil,shsnow,time_index=tidx)
#        hsnow = hsnow[:,:]
#        qtot = Mrdg.nc_get_var(ncfil,sqtot,time_index=tidx)
#        qtot = qtot[:,:]
#        qcool = Mrdg.nc_get_var(ncfil,sqcool,time_index=tidx)
#        qcool = qcool[:,:]
#        qother = Mrdg.nc_get_var(ncfil,sqother,time_index=tidx)
#        qother = qother[:,:]
#        qatm = Mrdg.nc_get_var(ncfil,sqatm,time_index=tidx)
#        qatm = qatm[:,:]
#        conc_old = Mrdg.nc_get_var(ncfil,sconc_old,time_index=tidx)
#        conc_old = conc_old[:,:]
#        thic_old = Mrdg.nc_get_var(ncfil,sthic_old,time_index=tidx)
#        thic_old = thic_old[:,:]
#        hsnow_old = Mrdg.nc_get_var(ncfil,shsnow_old,time_index=tidx)
#        hsnow_old = hsnow_old[:,:]
#        hsnow_diff = Mrdg.nc_get_var(ncfil,shsnow_diff,time_index=tidx)
#        hsnow_diff = hsnow_diff[:,:]
#        conc_diff = Mrdg.nc_get_var(ncfil,sconc_diff,time_index=tidx)
#        conc_diff = conc_diff[:,:]
#        thic_diff = Mrdg.nc_get_var(ncfil,sthic_diff,time_index=tidx)
#        thic_diff = thic_diff[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#
#        self.ZCI1,self.ZTI1,self.HSI1,self.QT1,self.QC1,self.QO1,self.QA1,self.FOI1,\
#                self.TOI1,self.HSO1,self.HSD1,self.FDI1,self.TDI1,self.TXI1,self.TYI1 =\
#                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
#                hsnow_diff,conc_diff,thic_diff,taux,tauy
#        
#        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
#                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')
#
#        return()
#
#    def _read_mdl_ice_only_2(self,dadate,basemap):
#        
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[2]
#
#        # Read TP4arch
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/'+\
#                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        shsnow = 'hsnow'
#        sqtot = 'qtot'
#        sqcool = 'flx_cool'
#        sqother = 'flx_itop'
#        sqatm = 'flx_ow'
#        sconc_old = 'fi_old'
#        sthic_old = 'hi_old'
#        shsnow_old = 'hs_old'
#        shsnow_diff = 'hs_int'
#        sconc_diff = 'dfice'
#        sthic_diff = 'dhice'
#        staux = 'taux'
#        stauy = 'tauy'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMI,self.latMI,self.XI,self.YI = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        hsnow = Mrdg.nc_get_var(ncfil,shsnow,time_index=tidx)
#        hsnow = hsnow[:,:]
#        qtot = Mrdg.nc_get_var(ncfil,sqtot,time_index=tidx)
#        qtot = qtot[:,:]
#        qcool = Mrdg.nc_get_var(ncfil,sqcool,time_index=tidx)
#        qcool = qcool[:,:]
#        qother = Mrdg.nc_get_var(ncfil,sqother,time_index=tidx)
#        qother = qother[:,:]
#        qatm = Mrdg.nc_get_var(ncfil,sqatm,time_index=tidx)
#        qatm = qatm[:,:]
#        conc_old = Mrdg.nc_get_var(ncfil,sconc_old,time_index=tidx)
#        conc_old = conc_old[:,:]
#        thic_old = Mrdg.nc_get_var(ncfil,sthic_old,time_index=tidx)
#        thic_old = thic_old[:,:]
#        hsnow_old = Mrdg.nc_get_var(ncfil,shsnow_old,time_index=tidx)
#        hsnow_old = hsnow_old[:,:]
#        hsnow_diff = Mrdg.nc_get_var(ncfil,shsnow_diff,time_index=tidx)
#        hsnow_diff = hsnow_diff[:,:]
#        conc_diff = Mrdg.nc_get_var(ncfil,sconc_diff,time_index=tidx)
#        conc_diff = conc_diff[:,:]
#        thic_diff = Mrdg.nc_get_var(ncfil,sthic_diff,time_index=tidx)
#        thic_diff = thic_diff[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#
#        self.ZCI2,self.ZTI2,self.HSI2,self.QT2,self.QC2,self.QO2,self.QA2,self.FOI2,\
#                self.TOI2,self.HSO2,self.HSD2,self.FDI2,self.TDI2,self.TXI2,self.TYI2 =\
#                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
#                hsnow_diff,conc_diff,thic_diff,taux,tauy
#        
#        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
#                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')
#
#        return()
#
#    def _read_mdl_ice_only_3(self,dadate,basemap):
#        
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[3]
#
#        # Read TP4arch
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/'+\
#                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        shsnow = 'hsnow'
#        sqtot = 'qtot'
#        sqcool = 'flx_cool'
#        sqother = 'flx_itop'
#        sqatm = 'flx_ow'
#        sconc_old = 'fi_old'
#        sthic_old = 'hi_old'
#        shsnow_old = 'hs_old'
#        shsnow_diff = 'hs_int'
#        sconc_diff = 'dfice'
#        sthic_diff = 'dhice'
#        staux = 'taux'
#        stauy = 'tauy'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMI,self.latMI,self.XI,self.YI = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        hsnow = Mrdg.nc_get_var(ncfil,shsnow,time_index=tidx)
#        hsnow = hsnow[:,:]
#        qtot = Mrdg.nc_get_var(ncfil,sqtot,time_index=tidx)
#        qtot = qtot[:,:]
#        qcool = Mrdg.nc_get_var(ncfil,sqcool,time_index=tidx)
#        qcool = qcool[:,:]
#        qother = Mrdg.nc_get_var(ncfil,sqother,time_index=tidx)
#        qother = qother[:,:]
#        qatm = Mrdg.nc_get_var(ncfil,sqatm,time_index=tidx)
#        qatm = qatm[:,:]
#        conc_old = Mrdg.nc_get_var(ncfil,sconc_old,time_index=tidx)
#        conc_old = conc_old[:,:]
#        thic_old = Mrdg.nc_get_var(ncfil,sthic_old,time_index=tidx)
#        thic_old = thic_old[:,:]
#        hsnow_old = Mrdg.nc_get_var(ncfil,shsnow_old,time_index=tidx)
#        hsnow_old = hsnow_old[:,:]
#        hsnow_diff = Mrdg.nc_get_var(ncfil,shsnow_diff,time_index=tidx)
#        hsnow_diff = hsnow_diff[:,:]
#        conc_diff = Mrdg.nc_get_var(ncfil,sconc_diff,time_index=tidx)
#        conc_diff = conc_diff[:,:]
#        thic_diff = Mrdg.nc_get_var(ncfil,sthic_diff,time_index=tidx)
#        thic_diff = thic_diff[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#
#        self.ZCI3,self.ZTI3,self.HSI3,self.QT3,self.QC3,self.QO3,self.QA3,self.FOI3,\
#                self.TOI3,self.HSO3,self.HSD3,self.FDI3,self.TDI3,self.TXI3,self.TYI3 =\
#                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
#                hsnow_diff,conc_diff,thic_diff,taux,tauy
#        
#        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
#                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')
#
#        return()
#
#    def _read_mdl_ice_only_4(self,dadate,basemap):
#        
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[4]
#
#        # Read TP4arch
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/'+\
#                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        shsnow = 'hsnow'
#        sqtot = 'qtot'
#        sqcool = 'flx_cool'
#        sqother = 'flx_itop'
#        sqatm = 'flx_ow'
#        sconc_old = 'fi_old'
#        sthic_old = 'hi_old'
#        shsnow_old = 'hs_old'
#        shsnow_diff = 'hs_int'
#        sconc_diff = 'dfice'
#        sthic_diff = 'dhice'
#        staux = 'taux'
#        stauy = 'tauy'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMI,self.latMI,self.XI,self.YI = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        hsnow = Mrdg.nc_get_var(ncfil,shsnow,time_index=tidx)
#        hsnow = hsnow[:,:]
#        qtot = Mrdg.nc_get_var(ncfil,sqtot,time_index=tidx)
#        qtot = qtot[:,:]
#        qcool = Mrdg.nc_get_var(ncfil,sqcool,time_index=tidx)
#        qcool = qcool[:,:]
#        qother = Mrdg.nc_get_var(ncfil,sqother,time_index=tidx)
#        qother = qother[:,:]
#        qatm = Mrdg.nc_get_var(ncfil,sqatm,time_index=tidx)
#        qatm = qatm[:,:]
#        conc_old = Mrdg.nc_get_var(ncfil,sconc_old,time_index=tidx)
#        conc_old = conc_old[:,:]
#        thic_old = Mrdg.nc_get_var(ncfil,sthic_old,time_index=tidx)
#        thic_old = thic_old[:,:]
#        hsnow_old = Mrdg.nc_get_var(ncfil,shsnow_old,time_index=tidx)
#        hsnow_old = hsnow_old[:,:]
#        hsnow_diff = Mrdg.nc_get_var(ncfil,shsnow_diff,time_index=tidx)
#        hsnow_diff = hsnow_diff[:,:]
#        conc_diff = Mrdg.nc_get_var(ncfil,sconc_diff,time_index=tidx)
#        conc_diff = conc_diff[:,:]
#        thic_diff = Mrdg.nc_get_var(ncfil,sthic_diff,time_index=tidx)
#        thic_diff = thic_diff[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#
#        self.ZCI4,self.ZTI4,self.HSI4,self.QT4,self.QC4,self.QO4,self.QA4,self.FOI4,\
#                self.TOI4,self.HSO4,self.HSD4,self.FDI4,self.TDI4,self.TXI4,self.TYI4 =\
#                conc,thic,shsnow,qtot,qcool,qother,qatm,conc_old,thic_old,hsnow_old,\
#                hsnow_diff,conc_diff,thic_diff,taux,tauy
#        
#        print('Ice Only: lonM,latM,X,Y,ZC,ZT,hsnow,qtot,qcool,qother,qatm,conc_old,thic_old,'+\
#                'hsnow_old,hsnow_diff,conc_diff,thic_diff,taux,tauy')
#
#        return()
#
#    def _read_mdl_waves_ice_1(self,dadate,basemap):
#    
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#        jday = jdaystart + 1
#        gigio = datetime.datetime.strptime( str(year)+str(jday), '%Y%j')
#        stdate = gigio.strftime('%Y%m%d')
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[1]
#
#        # Read TP4arch_wav
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/\
#                waves_ice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        sdmax = 'dmax'
#        sswh = 'swh'
#        staux = 'taux_wav'
#        stauy = 'tauy_wav'
#        smwd = 'mwd'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMW,self.latMW,self.XW,self.YW = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=tidx)
#        dmax = dmax[:,:]
#        swh = Mrdg.nc_get_var(ncfil,sswh,time_index=tidx)
#        swh = swh[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#        mwd = Mrdg.nc_get_var(ncfil,smwd,time_index=tidx)
#        mwd = mwd[:,:]
#        self.ZCW1,self.ZTW1,self.ZDW1,self.SWHW1,self.MWDW1,self.TXW1,self.TYW1 =\
#                conc,thic,dmax,swh,mwd,taux,tauy
#
#        print('Wave Ice: lonM,latM,X,Y,ZC,ZT,ZD,swh,mwd,taux,tauy')
#       
#        return()
#
#    def _read_mdl_waves_ice_2(self,dadate,basemap):
#    
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#        jday = jdaystart + 1
#        gigio = datetime.datetime.strptime( str(year)+str(jday), '%Y%j')
#        stdate = gigio.strftime('%Y%m%d')
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[2]
#
#        # Read TP4arch_wav
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/\
#                waves_ice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        sdmax = 'dmax'
#        sswh = 'swh'
#        staux = 'taux_wav'
#        stauy = 'tauy_wav'
#        smwd = 'mwd'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMW,self.latMW,self.XW,self.YW = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=tidx)
#        dmax = dmax[:,:]
#        swh = Mrdg.nc_get_var(ncfil,sswh,time_index=tidx)
#        swh = swh[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#        mwd = Mrdg.nc_get_var(ncfil,smwd,time_index=tidx)
#        mwd = mwd[:,:]
#        self.ZCW2,self.ZTW2,self.ZDW2,self.SWHW2,self.MWDW2,self.TXW2,self.TYW2 =\
#                conc,thic,dmax,swh,mwd,taux,tauy
#
#        print('Wave Ice: lonM,latM,X,Y,ZC,ZT,ZD,swh,mwd,taux,tauy')
#       
#        return()
#
#    def _read_mdl_waves_ice_3(self,dadate,basemap):
#    
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#        jday = jdaystart + 1
#        gigio = datetime.datetime.strptime( str(year)+str(jday), '%Y%j')
#        stdate = gigio.strftime('%Y%m%d')
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[3]
#
#        # Read TP4arch_wav
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/\
#                waves_ice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        sdmax = 'dmax'
#        sswh = 'swh'
#        staux = 'taux_wav'
#        stauy = 'tauy_wav'
#        smwd = 'mwd'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMW,self.latMW,self.XW,self.YW = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=tidx)
#        dmax = dmax[:,:]
#        swh = Mrdg.nc_get_var(ncfil,sswh,time_index=tidx)
#        swh = swh[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#        mwd = Mrdg.nc_get_var(ncfil,smwd,time_index=tidx)
#        mwd = mwd[:,:]
#        self.ZCW3,self.ZTW3,self.ZDW3,self.SWHW3,self.MWDW3,self.TXW3,self.TYW3 =\
#                conc,thic,dmax,swh,mwd,taux,tauy
#
#        print('Wave Ice: lonM,latM,X,Y,ZC,ZT,ZD,swh,mwd,taux,tauy')
#       
#        return()
#
#    def _read_mdl_waves_ice_4(self,dadate,basemap):
#    
#        # NetCDF reader 
#        day = self.day
#        month = self.month
#        year = self.year
#        jdaystart = self.jdaystart
#        jday = jdaystart + 1
#        gigio = datetime.datetime.strptime( str(year)+str(jday), '%Y%j')
#        stdate = gigio.strftime('%Y%m%d')
#
#        # TIME INDEX
#        tidx = self.WIM_timeframes[4]
#
#        # Read TP4arch_wav
#        outdir = '/work/timill/RealTime_Models/results_hindicast/TP4a0.12/\
#                waves_ice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
#        outdir = '/home/charlie/Documents/ncdata'
#        ncfil = outdir+'/SWARP_hindcast_wavesice_start'+dadate+'T000000Z.nc'
#        slon = 'longitude'
#        slat = 'latitude'
#        sconc = 'icec'
#        sthic = 'icetk'
#        sdmax = 'dmax'
#        sswh = 'swh'
#        staux = 'taux_wav'
#        stauy = 'tauy_wav'
#        smwd = 'mwd'
#        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
#        lonM = lonM[:,:]
#        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
#        latM = latM[:,:]
#        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
#        self.lonMW,self.latMW,self.XW,self.YW = lonM,latM,X,Y
#        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=tidx)
#        conc = conc[:,:]
#        thic = Mrdg.nc_get_var(ncfil,sthic,time_index=tidx)
#        thic = thic[:,:]
#        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=tidx)
#        dmax = dmax[:,:]
#        swh = Mrdg.nc_get_var(ncfil,sswh,time_index=tidx)
#        swh = swh[:,:]
#        taux = Mrdg.nc_get_var(ncfil,staux,time_index=tidx)
#        taux = taux[:,:]
#        tauy = Mrdg.nc_get_var(ncfil,stauy,time_index=tidx)
#        tauy = tauy[:,:]
#        mwd = Mrdg.nc_get_var(ncfil,smwd,time_index=tidx)
#        mwd = mwd[:,:]
#        self.ZCW4,self.ZTW4,self.ZDW4,self.SWHW4,self.MWDW4,self.TXW4,self.TYW4 =\
#                conc,thic,dmax,swh,mwd,taux,tauy
#
#        print('Wave Ice: lonM,latM,X,Y,ZC,ZT,ZD,swh,mwd,taux,tauy')
#       
#        return()
 
# Read in the data
#if run == 'OSI':
#    osisaf = reader('Osisaf',dadate,iqm)
#    lonO,latO,X,Y,Z = osisaf.lonO,osisaf.latO,osisaf.X,osisaf.Y,osisaf.Z
#    regional_daily(X,Y,Z,basemap=lqm,run=run)
#elif run == 'MIC':
#    model = reader('Model',dadate,iqm)
#    lonM,latM,XMC,YMC,ZMM,ZMD = model.lonM,model.latM,model.X,model.Y,model.ZC,model.ZD
#    # Reproject model data
#    Z = reproj(XMC,YMC,ZMM,XO,YO)
#    X = np.copy(XO)
#    Y = np.copy(YO)
#    regional_daily(X,Y,Z,basemap=lqm,run=run)
#elif run == 'MFD':
#    model = reader('Model',dadate,iqm)
#    lonM,latM,XMF,YMF,ZMM,ZMD = model.lonM,model.latM,model.X,model.Y,model.ZC,model.ZD
#    # Reproject model data
#    Z = reproj(XC,YC,ZMD,XO,YO)
#    X = np.copy(XO)
#    Y = np.copy(YO)
#    regional_daily(X,Y,Z,basemap=lqm,run=run)
#else:
#    print('Manual Mode')

##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## RUNS
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
############################################################################
## REGIONAL DAILY RUN 
#def regional_daily(X,Y,Z,basemap=None,run=None):
#    
#    import MIZchar as widths
#
#    # We need the contours and 4 state map (IP=2,MIZ=1,Ocean=0,Land=NaN)
#    if run == 'OSI' or run == 'MIC':
#        mode = 'ic'
#        MIZ_p = MIZ_poly(X,Y,Z,mode)
#    elif run == 'MFD':
#        mode = 'fsd'
#        MIZ_p = MIZ_poly(X,Y,Z,mode)
#    else:
#        print('Unknown Run')
#    
#    IPMIZ,B1,B2,MIZ_cont = MIZ_p.MIZ, MIZ_p.B1, MIZ_p.B2, MIZ_p.MIZcont
#
#    # Prep the lists for polygon.txt
#    global bar_poly_list
#    bar_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
#    global les_poly_list
#    les_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
#    global ncb_poly_list
#    ncb_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
#    global lab_poly_list
#    lab_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
#    global gre_poly_list
#    gre_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
#
#    # This list of classes may be useful
#    poly_list=[]
#
#    # Analysing the polygons one by one, dividing them into their respective regions
#    for n,el in enumerate(MIZ_cont):
#        polygons=poly_stat(el,B2,B1,IPMIZ,X,Y,n,basemap=basemap)
#        poly_list.append(polygons)
#
#    # BAR REGION
#    reg_repo = str(out_dir)+'/bar'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#    # DAILY
#    dailydir = reg_repo+'/daily'
#    if not os.path.exists(dailydir):
#        os.mkdir(dailydir)
#    # need variables for weighted location
#    num_lon = 0
#    num_lat = 0
#    den = 0
#    for poly in poly_list:
#        if poly.region == 'bar':
#            filname = str(dailydir)+'/'+str(dadate)+'.txt'
#            with open(filname,'a') as f:
#                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
#                    poly.centroid_longitude,poly.centroid_latitude]]
#                # Keep some info aside
#                num_lon += poly.centroid_longitude*poly.E_area 
#                num_lat += poly.centroid_latitude*poly.E_area 
#                den += poly.E_area
#                str1 = ' '.join(map(str,row1))
#                f.write(str1+'\n')
#                f.close
#    mean_lon_bar = num_lon/float(den)
#    mean_lat_bar = num_lat/float(den)
#    # POLYGONS
#    polydir = reg_repo+'/polygons'
#    if not os.path.exists(polydir):
#        os.mkdir(polydir)
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(bar_poly_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    # WIDTHS
#    widthsdir = reg_repo+'/widths'
#    if not os.path.exists(widthsdir):
#        os.mkdir(widthsdir)
#    reg_width_list = []
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    reg_widths = widths.single_file(filname,basemap)
#    for w,wd in enumerate(reg_widths):
#        if wd is not None:
#            reg_width_list.append(wd.int_widths)
#    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(reg_width_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    uglyway = open(filname,'r+')
#    uglierway = uglyway.readlines()
#    w_num = 0
#    w_den = 0
#    for el in uglierway:
#        pw = np.array(map(float,el.split()))
#        w_num += sum(pw)
#        w_den += len(pw)
#    avg_width_bar = w_num/float(w_den)
#
#    # LES REGION
#    reg_repo = str(out_dir)+'/les'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#    # DAILY
#    dailydir = reg_repo+'/daily'
#    if not os.path.exists(dailydir):
#        os.mkdir(dailydir)
#    # need variables for weighted location
#    num_lon = 0
#    num_lat = 0
#    den = 0
#    for poly in poly_list:
#        if poly.region == 'les':
#            filname = str(dailydir)+'/'+str(dadate)+'.txt'
#            with open(filname,'a') as f:
#                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
#                    poly.centroid_longitude,poly.centroid_latitude]]
#                # Keep some info aside
#                num_lon += poly.centroid_longitude*poly.E_area 
#                num_lat += poly.centroid_latitude*poly.E_area 
#                den += poly.E_area
#                str1 = ' '.join(map(str,row1))
#                f.write(str1+'\n')
#                f.close
#    mean_lon_les = num_lon/float(den)
#    mean_lat_les = num_lat/float(den)
#    # POLYGONS
#    polydir = reg_repo+'/polygons'
#    if not os.path.exists(polydir):
#        os.mkdir(polydir)
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(les_poly_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    # WIDTHS
#    widthsdir = reg_repo+'/widths'
#    if not os.path.exists(widthsdir):
#        os.mkdir(widthsdir)
#    reg_width_list = []
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    reg_widths = widths.single_file(filname,basemap)
#    for w,wd in enumerate(reg_widths):
#        if wd is not None:
#            reg_width_list.append(wd.int_widths)
#    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(reg_width_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    uglyway = open(filname,'r+')
#    uglierway = uglyway.readlines()
#    w_num = 0
#    w_den = 0
#    for el in uglierway:
#        pw = np.array(map(float,el.split()))
#        w_num += sum(pw)
#        w_den += len(pw)
#    avg_width_les = w_num/float(w_den)
#
#    # NCB REGION
#    reg_repo = str(out_dir)+'/ncb'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#    # DAILY
#    dailydir = reg_repo+'/daily'
#    if not os.path.exists(dailydir):
#        os.mkdir(dailydir)
#    # need variables for weighted location
#    num_lon = 0
#    num_lat = 0
#    den = 0
#    for poly in poly_list:
#        if poly.region == 'ncb':
#            filname = str(dailydir)+'/'+str(dadate)+'.txt'
#            with open(filname,'a') as f:
#                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
#                    poly.centroid_longitude,poly.centroid_latitude]]
#                # Keep some info aside
#                num_lon += poly.centroid_longitude*poly.E_area 
#                num_lat += poly.centroid_latitude*poly.E_area 
#                den += poly.E_area
#                str1 = ' '.join(map(str,row1))
#                f.write(str1+'\n')
#                f.close
#    mean_lon_ncb = num_lon/float(den)
#    mean_lat_ncb = num_lat/float(den)
#    # POLYGONS
#    polydir = reg_repo+'/polygons'
#    if not os.path.exists(polydir):
#        os.mkdir(polydir)
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(ncb_poly_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    # WIDTHS
#    widthsdir = reg_repo+'/widths'
#    if not os.path.exists(widthsdir):
#        os.mkdir(widthsdir)
#    reg_width_list = []
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    reg_widths = widths.single_file(filname,basemap)
#    for w,wd in enumerate(reg_widths):
#        if wd is not None:
#            reg_width_list.append(wd.int_widths)
#    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(reg_width_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    uglyway = open(filname,'r+')
#    uglierway = uglyway.readlines()
#    w_num = 0
#    w_den = 0
#    for el in uglierway:
#        pw = np.array(map(float,el.split()))
#        w_num += sum(pw)
#        w_den += len(pw)
#    avg_width_ncb = w_num/float(w_den)
#
#    # LAB REGION
#    reg_repo = str(out_dir)+'/lab'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#    # DAILY
#    dailydir = reg_repo+'/daily'
#    if not os.path.exists(dailydir):
#        os.mkdir(dailydir)
#    # need variables for weighted location
#    num_lon = 0
#    num_lat = 0
#    den = 0
#    for poly in poly_list:
#        if poly.region == 'lab':
#            filname = str(dailydir)+'/'+str(dadate)+'.txt'
#            with open(filname,'a') as f:
#                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
#                    poly.centroid_longitude,poly.centroid_latitude]]
#                # Keep some info aside
#                num_lon += poly.centroid_longitude*poly.E_area 
#                num_lat += poly.centroid_latitude*poly.E_area 
#                den += poly.E_area
#                str1 = ' '.join(map(str,row1))
#                f.write(str1+'\n')
#                f.close
#    mean_lon_lab = num_lon/float(den)
#    mean_lat_lab = num_lat/float(den)
#    # POLYGONS
#    polydir = reg_repo+'/polygons'
#    if not os.path.exists(polydir):
#        os.mkdir(polydir)
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(lab_poly_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    # WIDTHS
#    widthsdir = reg_repo+'/widths'
#    if not os.path.exists(widthsdir):
#        os.mkdir(widthsdir)
#    reg_width_list = []
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    reg_widths = widths.single_file(filname,basemap)
#    for w,wd in enumerate(reg_widths):
#        if wd is not None:
#            reg_width_list.append(wd.int_widths)
#    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(reg_width_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    uglyway = open(filname,'r+')
#    uglierway = uglyway.readlines()
#    w_num = 0
#    w_den = 0
#    for el in uglierway:
#        pw = np.array(map(float,el.split()))
#        w_num += sum(pw)
#        w_den += len(pw)
#    avg_width_lab = w_num/float(w_den)
#
#    # GRE REGION
#    reg_repo = str(out_dir)+'/gre'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#    # DAILY
#    dailydir = reg_repo+'/daily'
#    if not os.path.exists(dailydir):
#        os.mkdir(dailydir)
#    # need variables for weighted location
#    num_lon = 0
#    num_lat = 0
#    den = 0
#    for poly in poly_list:
#        if poly.region == 'gre':
#            filname = str(dailydir)+'/'+str(dadate)+'.txt'
#            with open(filname,'a') as f:
#                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
#                    poly.centroid_longitude,poly.centroid_latitude]]
#                # Keep some info aside
#                num_lon += poly.centroid_longitude*poly.E_area 
#                num_lat += poly.centroid_latitude*poly.E_area 
#                den += poly.E_area
#                str1 = ' '.join(map(str,row1))
#                f.write(str1+'\n')
#                f.close
#    mean_lon_gre = num_lon/float(den)
#    mean_lat_gre = num_lat/float(den)
#    # POLYGONS
#    polydir = reg_repo+'/polygons'
#    if not os.path.exists(polydir):
#        os.mkdir(polydir)
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(gre_poly_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    # WIDTHS
#    widthsdir = reg_repo+'/widths'
#    if not os.path.exists(widthsdir):
#        os.mkdir(widthsdir)
#    reg_width_list = []
#    filname = str(polydir)+'/'+str(dadate)+'.txt'
#    reg_widths = widths.single_file(filname,basemap)
#    for w,wd in enumerate(reg_widths):
#        if wd is not None:
#            reg_width_list.append(wd.int_widths)
#    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
#    with open(filname,'a') as f:
#        for l,el in enumerate(reg_width_list):
#            string = ' '.join(map(str,el))
#            for item in string:
#                f.write(item)
#            f.write('\n')
#        f.close()
#    uglyway = open(filname,'r+')
#    uglierway = uglyway.readlines()
#    w_num = 0
#    w_den = 0
#    for el in uglierway:
#        pw = np.array(map(float,el.split()))
#        w_num += sum(pw)
#        w_den += len(pw)
#    avg_width_gre = w_num/float(w_den)
#
#    # daily regional gen_stats
#    # Widhts, lons and lats lists
#    widths = [avg_width_bar,avg_width_les,avg_width_ncb,avg_width_lab,avg_width_gre]
#    lons = [avg_width_bar,avg_width_les,avg_width_ncb,avg_width_lab,avg_width_gre]
#    lats = [avg_width_bar,avg_width_les,avg_width_ncb,avg_width_lab,avg_width_gre]
#    
#    regional_stats(Z,widths,lons,lats,mode)
#
#    plt.imshow(IPMIZ)
#    plt.colorbar()
#    plt.savefig(str(out_dir)+'/'+str(dadate)+'.png',bbox_inches='tight')
#    plt.close()
#   


#############################################################################
## OLD OLD OLD OLD OLD OLD but useful?
#class MIZ_poly:
#    def __init__(self,X1,Y1,Z1,typo):
#        ZM = np.copy(Z1)
#        if typo == 'ic':
#            self.MIZraw,self.B1,self.B2 = self.binary_diff(ZM,.15,.80)
#        elif typo == 'fsd':
#            self.MIZraw,self.B1,self.B2 = self.binary_diff(ZM,.1,300)
#        self.MIZcont = self.poly_maker(self.MIZraw)
#        self.MIZ = self.mapper(self.MIZraw,ZM,typo)
#        return
#    
#    def binary_diff(self,data,threshm,threshM):
#        def closing(data):
#            # apply closing to avoid small polynias and clean up a little
#            kernel = np.ones((3,3),np.uint8)
#            data_cl = morph.closing(data,kernel)
#            return(data_cl)
#        
#        # generating the binary file, get the difference, divide into overprediction and
#        # underprediction, apply the closing on both and finally re-map the land
#        ndata1 = np.copy(data)
#        ndata2 = np.copy(data)
#
#        # OLD method where ice-ice = ocean-ocean = 0
#        ndata1[ndata1<threshm] = 0
#        ndata2[ndata2<threshM] = 0
#        ndata1[ndata1>=threshm] = 1  
#        ndata2[ndata2>=threshM] = 1  
#
#        # Difference MODEL - OSISAF
#        ddata = ndata1 - ndata2
#        # Cut out the nans for both
#        thenans = np.isnan(ddata)
#        ddata[thenans] = 0
#        # Closing operation (see def closing)
#        ddata = closing(ddata)
#        # Add the under to the over to get total diff.
#        return(ddata,ndata1,ndata2)
#    
#    def poly_maker(self,ddata):
#        # finding contours from the difference data map
#        MIZcont = msr.find_contours(ddata,.5)
#        MIZcont = sorted(MIZcont, key=len)
#        MIZcont = MIZcont[::-1]
#        return(MIZcont)
#    
#    def mapper(self,DIFF,Z1,typo):
#        DN = np.copy(DIFF)
#        IP = np.copy(Z1)
#        ZM = np.copy(Z1)
#        thenan = np.isnan(ZM)
#        if typo == 'ic':
#            IP[IP>=.80] = 2
#            IP[IP<.80] = 0
#        elif typo == 'fsd':
#            IP[IP==300] = 2
#            IP[IP<300] = 0
#        #getting back the ice pack
#        DN = IP + DN
#        DN[DN>2] = 1
#        # getting back the terrain lost during binary_mod (for contour finding reasons)
#        DN[thenan] = None
#        return(DN)
#
#class poly_stat:
#    def __init__(self,cont,OUT,INS,IPMIZ,X,Y,number,basemap=None):
#        self.number = number #get the number
#        self.ij_list = cont #get the contour
#        self.class_def(cont) #calculate class from contour
#        xy_list = self.ij2xy(cont,X,Y) #get xy from ij (cont)
#        self.ip_miz = IPMIZ #get ip_miz
#        f_vals = self.split_cont(OUT,INS) #
#        lon,lat = self.lon_lat(xy_list,basemap)
#        self.area_perimeter()
#        self.add2list(self.region,number,lon,lat,f_vals)
#
#    # class definition
#    def class_def(self,contour):
#        lcont = len(contour)
#        if lcont > 200:
#            self.polygon_class = 'H'
#        elif lcont > 100 and lcont <= 200:
#            self.polygon_class = 'B'
#        elif lcont > 30 and lcont <= 100:
#            self.polygon_class = 'M'
#        elif lcont <= 30:
#            self.polygon_class = 'S'
#        return()
#
#        #self.stat_chart(save=True)
#
#    # output: xy_coords
#    # NOTE weird error displayed (non-integer index?) it works but I need to check this
#    def ij2xy(self,cont,X,Y):
#        # changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
#        x = []
#        y = []
#        xvec = range(len(cont[:,0]))
#        for n,en in enumerate(xvec):
#            en = X[cont[n,0]][cont[n,1]]
#            x.append(en)
#        yvec = range(len(cont[:,1]))
#        for n,en in enumerate(yvec):
#            en = Y[cont[n,0]][cont[n,1]]
#            y.append(en)
#        xy_list = zip(x,y)
#        xy_list = np.array(xy_list)
#        self.xy_list = xy_list
#        return(xy_list)
#    
#    # output: f_vals
#    def split_cont(self,OUT,INS):
#        # this function find the different contours
#        # NOTE if a point has non integer i coordinate is going to be a vertical edge,
#        # if has non integer j coordinate is going to be a horizontal edge hence the use
#        # of different arounds depending on the indexes of the point
#        # NOTE if the polygon is an OVERESTIMATION the contour finding is inverted
#        vs = self.ij_list # list of (i,j) pixel indices
#        in_cont = []
#        out_cont = []
#        unk_cont = []
#        func_vals= []
#        func_mod = 0 # value of func_vals at model ice edge
#        func_osi = 1 # value of func_vals at OSISAF ice edge
#        func_unk = 2 # value of func_vals at other parts of contour
#        
#        for n,el in enumerate(vs):
#            #getting all the neighbours
#            around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
#            around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
#            check_cont = 0
#            if el[0]/int(el[0]) == 1:
#                for h,v in around1:
#                    if OUT[h][v] == INS[h][v] == 1:
#                        in_cont.append(el)
#                        func_val=func_mod
#                        check_cont = 1
#                    elif OUT[h][v] == INS[h][v] == 0:
#                        out_cont.append(el)
#                        func_val=func_osi
#                        check_cont = 1
#                if check_cont == 0:
#                    unk_cont.append(el)
#                    func_val=func_unk
#                func_vals.append(func_val)
#            else:
#                for h,v in around2:
#                    if OUT[h][v] == INS[h][v] == 1:
#                        in_cont.append(el)
#                        func_val=func_mod
#                        check_cont = 1
#                    elif OUT[h][v] == INS[h][v] == 0:
#                        out_cont.append(el)
#                        func_val=func_osi
#                        check_cont = 1
#                if check_cont == 0:
#                    unk_cont.append(el)
#                    func_val=func_unk
#                func_vals.append(func_val)
#
#        # let's try the above method, if it doesn't work, use bot
#        #for n,el in enumerate(vs):
#        #    #getting all the neighbours
#        #    around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
#        #    around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
#        #    check_cont = 0
#        #    if el[0]/int(el[0]) == 1:
#        #        for h,v in around1:
#        #            if OUT[h][v] == INS[h][v] == 0:
#        #                in_cont.append(el)
#        #                func_val=func_mod
#        #                check_cont = 1
#        #            elif OUT[h][v] == INS[h][v] == 1:
#        #                out_cont.append(el)
#        #                func_val=func_osi
#        #                check_cont = 1
#        #        if check_cont == 0:
#        #            unk_cont.append(el)
#        #            func_val=func_unk
#        #        func_vals.append(func_val)
#        #    else:
#        #        for h,v in around2:
#        #            if OUT[h][v] == INS[h][v] == 0:
#        #                in_cont.append(el)
#        #                func_val=func_mod
#        #                check_cont = 1
#        #            elif OUT[h][v] == INS[h][v] == 1:
#        #                out_cont.append(el)
#        #                func_val=func_osi
#        #                check_cont = 1
#        #        if check_cont == 0:
#        #            unk_cont.append(el)
#        #            func_val=func_unk
#        #        func_vals.append(func_val)
#        
#        in_cont = np.array(in_cont)
#        out_cont = np.array(out_cont)
#        unk_cont = np.array(unk_cont)
#        f_vals = np.array(func_vals)
#        self.in_ice_edge = in_cont
#        self.out_ice_edge = out_cont
#        self.unknown_edge = unk_cont
#        self.f_vals = f_vals
#        return(f_vals)
#    
#    # output: lon,lat,region of polygon
#    # NOTE cont_list might be XY or IJ (no LONLAT)
#    def lon_lat(self,cont_list,basemap,typo='xy'):
#
#        # Centroid - Location in lon/lat of the centroid of the polygon
#        DX = 0
#        DY = 0
#        B = 0
#        for n in range(len(cont_list)-1):
#            DX += ((cont_list[n][0]+cont_list[n+1][0])*\
#                  ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1])))  
#            DY += ((cont_list[n][1]+cont_list[n+1][1])*\
#                  ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1])))  
#            B += ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1]))
#        A = (0.5)*B 
#        CX = (1/float(6*A))*DX
#        CY = (1/float(6*A))*DY
#        self.centroid = [CX,CY]
#        if typo == 'xy':
#            centroid_longitude,centroid_latitude = basemap(CX,CY,inverse=True)
#            self.centroid_longitude = centroid_longitude
#            self.centroid_latitude = centroid_latitude
#            lon_list,lat_list=basemap(cont_list[:,0],cont_list[:,1],inverse=True)
#            bblone,bblate,bblonw,bblatw = baffin_bay()
#            # Defining the region
#            if centroid_longitude < 90 and centroid_longitude > 16:
#                self.region = 'bar'
#            elif centroid_longitude <= 16 and centroid_longitude > -44:
#                self.region = 'gre'
#            elif centroid_longitude <= -44 and centroid_longitude > -90:
#                n = 0
#                N = len(bblatw[:-1])
#                if centroid_latitude <= bblatw[-1] and centroid_longitude >= bblonw[-1]:
#                    self.region = 'lab'
#                else:
#                    while n < N:
#                        if centroid_latitude >= bblatw[n+1] and centroid_latitude <= bblatw[n]:
#                            if centroid_longitude >= bblonw[n]: 
#                                self.region = 'lab'
#                                break
#                        n += 1
#                    else:
#                        self.region = 'ncb'
#            elif centroid_longitude < 180 and centroid_longitude > 90:
#                self.region = 'les'
#            else:
#                self.region = 'ncb'
#        elif typo == 'ij':
#            regions = np.load('./geo_info/regions.npy')
#            if regions[CX,CY] == 1:
#                self.region = 'bar'
#            elif regions[CX,CY] == 2:
#                self.region = 'les'
#            elif regions[CX,CY] == 3:
#                self.region = 'ncb'
#            elif regions[CX,CY] == 4:
#                self.region = 'lab'
#            elif regions[CX,CY] == 5:
#                self.region = 'gre'
#       
#        return(lon_list,lat_list)
#
#    # output: euclidean area (E_area) and perimeter (E_perim)
#    def area_perimeter(self):
#        # Calculating the area of irregular polygon (from perimeter)
#        vs = self.ij_list
#        a = 0
#        x0,y0 = vs[0]
#        for [x1,y1] in vs[1:]:
#            dx = x1-x0
#            dy = y1-y0
#            a += 0.5*(y0*dx - x0*dy)
#            x0 = x1 
#            y0 = y1 
#        self.E_area = abs(a*100)
#        # Calculating perimeter in xy coordinates (unit = 10km)
#        perim = 0
#        for n in range(len(vs)-1):
#            perim += np.sqrt(pow(vs[n+1][0]-vs[n][0],2)+pow(vs[n+1][1]-vs[n][1],2))
#        self.E_perim = abs(perim*10)
#        return
#
#    # output: add the polygon points to the regional list
#    def add2list(self,region,number,lon,lat,fvals):
#
#        if region == 'bar':
#            for n, en in enumerate(fvals):
#                bar_poly_list.append([number,lon[n],lat[n],en])
#        elif region == 'les':
#            for n, en in enumerate(fvals):
#                les_poly_list.append([number,lon[n],lat[n],en])
#        elif region == 'ncb':
#            for n, en in enumerate(fvals):
#                ncb_poly_list.append([number,lon[n],lat[n],en])
#        elif region == 'lab':
#            for n, en in enumerate(fvals):
#                lab_poly_list.append([number,lon[n],lat[n],en])
#        elif region == 'gre':
#            for n, en in enumerate(fvals):
#                gre_poly_list.append([number,lon[n],lat[n],en])
#
#    # Function stat chart
#    def stat_chart(self,save=False):
#        # The Statistical Chart is an output that tries to compress as many informations
#        # as possible in a single figure.
#        # The figure is:
#        # 1) Top left arctic map with studied polygon highlighted
#        # 2) Top right a close up of the polygon with contours and min. distances
#        # 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
#        # 4) Bottom right is a legend for the small image
#        
#        nmbr = self.number
#        pname = 'Polygon_'+str(nmbr)
#        self.name = pname
#        pclass = self.polygon_class
#        region = self.region
#        print ''
#        print 'Statistic Chart for ',pname
#        print 'Class and Region ',pclass,region
#        DN = self.ip_miz
#        ij = self.ij_list
#        inside_contour = self.in_ice_edge
#        outside_contour = self.out_ice_edge
#        unknown_contour = self.unknown_edge
#        clon = '%1.2f' %self.centroid_longitude
#        clat = '%1.2f' %self.centroid_latitude
#        clonlat = '{0}/{1}'.format(clon,clat)
#        E_area = self.E_area
#        area = '%1.4e' %(E_area) 
#        E_perim = self.E_perim
#        perim = '%1.4e' %(E_perim)
#        
#        # Setting up the plot (2x2) and subplots
#        fig = plt.figure(figsize=(15,10))
#        gs = gridspec.GridSpec(2,2,width_ratios=[2,1],height_ratios=[4,2])
#        plt.suptitle(pname+', class '+pclass+', '+region,fontsize=18)
#        main = plt.subplot(gs[0,0])
#        polyf = plt.subplot(gs[0,1])
#        tab = plt.subplot(gs[1,0])
#        leg = plt.subplot(gs[1,1])
#        tab.set_xticks([])
#        leg.set_xticks([])
#        tab.set_yticks([])
#        leg.set_yticks([])
#        tab.set_frame_on(False)
#        leg.set_frame_on(False)
#        
#        # Main image on the top left
#        main.imshow(DN[::-1],interpolation='nearest',cmap='winter')
#        x1,x2,y1,y2 = np.min(ij[:,1])-10,np.max(ij[:,1])+10,np.min(ij[:,0])-10,np.max(ij[:,0])+10
#        main.axvspan(x1,x2,ymin=1-((y1-320)/float(len(DN)-320)),\
#            ymax=1-((y2-320)/float(len(DN)-320)),color='red',alpha=0.3)
#        main.axis([0,760,0,800])
#        
#        # Polygon image on the top right
#        polyf.imshow(DN,interpolation='nearest',cmap='winter')
#        polyf.axis([x1,x2,y2,y1])
#        if len(inside_contour) != 0:
#            polyf.plot(inside_contour[:,1],inside_contour[:,0],'ro',markersize=4)
#        if len(outside_contour) != 0:
#            polyf.plot(outside_contour[:,1],outside_contour[:,0],'yo',markersize=4)
#        if len(unknown_contour) != 0:
#            polyf.plot(unknown_contour[:,1],unknown_contour[:,0],'go',markersize=4)
#        
#        # Legend on the bottom right
#        mc = mlines.Line2D([],[],color='red',marker='o')
#        oc = mlines.Line2D([],[],color='yellow',marker='o')
#        uc = mlines.Line2D([],[],color='green',marker='o')
#        leg.legend([mc,oc,uc],('15% Cont.','80% Cont.','Unknown Cont.'))
#                 
#        # Statistics text on the bottom left
#        txt =  '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
#               '2) Area = '+str(area)+' km^2\n'+ \
#               '3) Perimeter = '+str(perim)+' km\n'+ \
#               '4) Avg. Width = '+'\n'+ \
#               '5) Widths SD = '
#        tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
#        if save:
#            outdir = str(out_dir)
#            if not os.path.exists(outdir):
#                os.mkdir(outdir)
#            valid_class = outdir+'/'+dadate
#            if not os.path.exists(valid_class):
#                os.mkdir(valid_class)
#            if not os.path.exists(valid_class+'/'+region):
#                os.mkdir(valid_class+'/'+region)
#            fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
#            plt.close()
#        else:
#            plt.show(False)
#        print 'Statistic chart done for '+str(pname)
#        print ''
#        return
#   
        #self._read_mdl_ice_only_1(dadate,mbasemap)
        #self._read_mdl_ice_only_2(dadate,mbasemap)
        #self._read_mdl_ice_only_3(dadate,mbasemap)
        #self._read_mdl_ice_only_4(dadate,mbasemap)
        #self._read_mdl_waves_ice_1(dadate,mbasemap)
        #self._read_mdl_waves_ice_2(dadate,mbasemap)
        #self._read_mdl_waves_ice_3(dadate,mbasemap)
        #self._read_mdl_waves_ice_4(dadate,mbasemap)
       
        #self.Dmean1,self.D2mean1,self.Slat1,self.Sbot1,self.Alpha1,self.Beta1 = self.floes_analysis(self.FOI1,self.TOI1,self.ZDW1)
        #self.qlat1,self.qvrt1 = self.heat_flux(self.QA1,self.QC1,self.QO1,self.FOI1)
        #self.glat_HYCOM1,self.gvrt_HYCOM1,self.glat_WIM1,self.gvrt_WIM1,self.CNW1,self.TNW1 =\
        #self.WIM_ice_growth(self.TOI1,self.TDI1,self.FOI1,self.FDI1,self.HSD1)

        #self.Dmean2,self.D2mean2,self.Slat2,self.Sbot2,self.Alpha2,self.Beta2 = self.floes_analysis(self.FOI2,self.TOI2,self.ZDW2)
        #self.qlat2,self.qvrt2 = self.heat_flux(self.QA2,self.QC2,self.QO2,self.FOI2)
        #self.glat_HYCOM2,self.gvrt_HYCOM2,self.glat_WIM2,self.gvrt_WIM2,self.CNW2,self.TNW2 =\
        #self.WIM_ice_growth(self.TOI2,self.TDI2,self.FOI2,self.FDI2,self.HSD2)
       
        #self.Dmean3,self.D2mean3,self.Slat3,self.Sbot3,self.Alpha3,self.Beta3 = self.floes_analysis(self.FOI3,self.TOI3,self.ZDW3)
        #self.qlat3,self.qvrt3 = self.heat_flux(self.QA3,self.QC3,self.QO3,self.FOI3)
        #self.glat_HYCOM3,self.gvrt_HYCOM3,self.glat_WIM3,self.gvrt_WIM3,self.CNW3,self.TNW3 =\
        #self.WIM_ice_growth(self.TOI3,self.TDI3,self.FOI3,self.FDI3,self.HSD3)
       
        #self.Dmean4,self.D2mean4,self.Slat4,self.Sbot4,self.Alpha4,self.Beta4 = self.floes_analysis(self.FOI4,self.TOI4,self.ZDW4)
        #self.qlat4,self.qvrt4 = self.heat_flux(self.QA4,self.QC4,self.QO4,self.FOI4)
        #self.glat_HYCOM4,self.gvrt_HYCOM4,self.glat_WIM4,self.gvrt_WIM4,self.CNW4,self.TNW4 =\
        #self.WIM_ice_growth(self.TOI4,self.TDI4,self.FOI4,self.FDI4,self.HSD4)

        # OVER_UNDER mask
        #self.wind_stress_x_over = self.over_poly * self.reproj_wim2osi(self.TXI)
        #self.wind_stress_y_over = self.over_poly * self.reproj_wim2osi(self.TYI)
        #self.wave_stress_x_over = self.over_poly * self.reproj_wim2osi(self.TXW)
        #self.wave_stress_y_over = self.over_poly * self.reproj_wim2osi(self.TYW)
        #self.glat_HYCOM_over = self.over_poly * self.reproj_wim2osi(self.glat_HYCOM)
        #self.gvrt_HYCOM_over = self.over_poly * self.reproj_wim2osi(self.gvrt_HYCOM)
        #self.glat_WIM_over = self.over_poly * self.reproj_wim2osi(self.glat_WIM)
        #self.gvrt_WIM_over = self.over_poly * self.reproj_wim2osi(self.gvrt_WIM)

        #self.wind_stress_x_under = self.under_poly * self.reproj_wim2osi(self.TXI)
        #self.wind_stress_y_under = self.under_poly * self.reproj_wim2osi(self.TYI)
        #self.wave_stress_x_under = self.under_poly * self.reproj_wim2osi(self.TXW)
        #self.wave_stress_y_under = self.under_poly * self.reproj_wim2osi(self.TYW)
        #self.glat_HYCOM_under = self.under_poly * self.reproj_wim2osi(self.glat_HYCOM)
        #self.gvrt_HYCOM_under = self.under_poly * self.reproj_wim2osi(self.gvrt_HYCOM)
        #self.glat_WIM_under = self.under_poly * self.reproj_wim2osi(self.glat_WIM)
        #self.gvrt_WIM_under = self.under_poly * self.reproj_wim2osi(self.gvrt_WIM)


