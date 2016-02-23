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
            jdaystart = '060'
        elif jday > 67 and jday < 75:
            jdaystart = '067'
        elif jday > 74 and jday < 82:
            jdaystart = '074'
        elif jday > 81 and jday < 89:
            jdaystart = '081'
        elif jday > 88 and jday < 96:
            jdaystart = '088'
        elif jday > 95 and jday < 103:
            jdaystart = '095'
        elif jday > 102 and jday < 110:
            jdaystart = '102'
        elif jday > 109 and jday < 117:
            jdaystart = '109'
        elif jday > 116 and jday < 124:
            jdaystart = '116'
        elif jday > 123 and jday < 131:
            jdaystart = '123'
        elif jday > 130 and jday < 138:
            jdaystart = '130'
        elif jday > 137 and jday < 145:
            jdaystart = '137'
        elif jday > 144 and jday < 159:
            jdaystart = '144'
        elif jday > 158 and jday < 166:
            jdaystart = '158'
        elif jday > 165 and jday < 173:
            jdaystart = '165'
        elif jday > 172 and jday < 180:
            jdaystart = '172'
        elif jday > 179 and jday < 187:
            jdaystart = '179'
        elif jday > 186 and jday < 194:
            jdaystart = '186'
        elif jday > 193 and jday < 201:
            jdaystart = '193'
        elif jday > 200 and jday < 208:
            jdaystart = '200'
        elif jday > 207 and jday < 215:
            jdaystart = '207'
        elif jday > 214 and jday < 222:
            jdaystart = '214'
        elif jday > 221 and jday < 236:
            jdaystart = '221'
        elif jday > 235 and jday < 243:
            jdaystart = '235'
        elif jday > 242 and jday < 250:
            jdaystart = '242'
        elif jday > 249 and jday < 257:
            jdaystart = '249'
        elif jday > 256 and jday < 264:
            jdaystart = '256'
        elif jday > 263 and jday < 271:
            jdaystart = '263'
        elif jday > 270 and jday < 275:
            jdaystart = '270'
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
        outdir = '/work/timill/giacomo/osisaf_repo'
        #outdir = '/home/charlie/Documents/ncdata'
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
        tidx = 0

        # Read TP4arch
        outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
                'ice_only/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
        #outdir = '/home/charlie/Documents/ncdata'
        ncfil = outdir+'/SWARP_hindcast_ice_only_start'+dadate+'T000000Z.nc'
        print outdir
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

        # TIME INDEX
        tidx = 0

        # Read TP4arch_wav
        outdir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/'+\
              'wavesice/2015_GOOD/2015_'+str(jdaystart)+'/final_output'
        #outdir = '/home/charlie/Documents/ncdata'
        print outdir
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
np.savetxt(str(out_dir)+'/'+str(dadate)+'_glat_H',data.reproj_wim2osi(data.glat_HYCOM))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_gvrt_H',data.reproj_wim2osi(data.gvrt_HYCOM))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_glat_W',data.reproj_wim2osi(data.glat_WIM))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_gvrt_W',data.reproj_wim2osi(data.gvrt_WIM))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_stress_H',data.reproj_wim2osi(data.stress_h))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_stress_W',data.reproj_wim2osi(data.stress_w))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_over',data.reproj_wim2osi(data.over_poly))
np.savetxt(str(out_dir)+'/'+str(dadate)+'_under',data.reproj_wim2osi(data.under_poly))


elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime

