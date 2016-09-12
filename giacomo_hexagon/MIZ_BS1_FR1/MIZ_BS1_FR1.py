# vv. 1.3 29/02/2016
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

sys.path.append('../utilities')
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


# From MIZwidth
class MIZwidth:
    def __init__(self,X,Y,Z,basemap,region,save=False):
        ZM = mask_region(Z,region)
        MIZraw,ice_ext,pack_ice = self.binary_diff(ZM,.15,.80)
        self.MIZraw,self.ice_ext,self.pack_ice = MIZraw,ice_ext,pack_ice

        MIZcont = self.poly_maker(self.MIZraw)
        self.MIZcont = MIZcont

        MIZ = self.mapper(self.MIZraw,ZM)
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
        if save:
            self.regional_save(ZM,poly_list,area_list,perim_list,clonlat_list,basemap,region)
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
                        in_cont.append(el)
                        func_val=func_mod
                        check_cont = 1
                    elif pack[h][v] == ext[h][v] == 0:
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
                        in_cont.append(el)
                        func_val=func_mod
                        check_cont = 1
                    elif pack[h][v] == ext[h][v] == 0:
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
        #for n in range(len(cont)-1):
        #    DX += ((cont[n][0]+cont[n+1][0])*\
        #          ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1])))  
        #    DY += ((cont[n][1]+cont[n+1][1])*\
        #          ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1])))  
        #    B += ((cont[n][0]*cont[n+1][1])-(cont[n+1][0]*cont[n][1]))
        #A = (0.5)*B 
        #CX = (1/float(6*A))*DX
        #CY = (1/float(6*A))*DY
        #self.centroid = [CX,CY]
        #self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
        #clonlat = [self.centroid_longitude,self.centroid_latitude]
        lon_list,lat_list=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
        clon = np.mean(lon_list)
        clat = np.mean(lat_list)
        clonlat = [clon,clat]
        print clonlat

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
            reg_repo = './outputs/MIZ/MDL/'+str(region)
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
            den += area[n]
        mean_lon = num_lon/float(den)
        mean_lat = num_lat/float(den)

        # Regional Lists
        # NOTE list is [ice extent, ice area, miz extent, miz area]
        lst = []
    
        # Which mode?
        if np.nanmax(data) <= 1:
            lindata = np.reshape(data,data.size)
            for pt in lindata:
                if pt > .15 and pt < .80:
                    lst.append([1,1])
                elif pt >= .80:
                    lst.append([1,0])
        else:
            lindata = np.reshape(data,data.size)
            for pt in lindata:
                if pt > .1 and pt < 300:
                    lst.append([1,1])
                elif pt == 300:
                    lst.append([1,0])

        stats = map(sum,zip(*lst))
    
        filname = reg_repo+'/stats.txt'
        with open(filname,'a') as f:
            row1 = [dadate,stats[0],stats[1],avg_width,mean_lon,mean_lat]
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
        mdl_reg = np.load('../geo_info/bar_gre_bal_reg.npy')
        osi_reg = np.load('../geo_info/bar_gre_bal_reg_OSI.npy')
        mdl_mask = np.load('../geo_info/bar_gre_bal_reg_mask.npy')
        osi_mask = np.load('../geo_info/bar_gre_bal_reg_mask_OSI.npy')
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

# Let's go, go, GO!
time1 = time.time()
print('Basemap creation...')
bmm = basemap_creator('TP4')
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

mrg = np.load('../geo_info/bar_gre_bal_reg.npy')
org = np.load('../geo_info/bar_gre_bal_reg_OSI.npy')
mms = np.load('../geo_info/bar_gre_bal_reg_mask.npy')
oms = np.load('../geo_info/bar_gre_bal_reg_mask_OSI.npy')

print('Reading datasets...')
data = reader(dadate,bmm,bmo)
print('DONE')
print('')

if region == 'bar':
   # Barents
   print('Studying the BARENTS region...')
   print('')
   
   ZM,XM,YM,latM,lonM,ZO,XO,YO,latO,lonO = data.ZM,data.XM,data.YM,data.latM,data.lonM,\
           data.ZO,data.XO,data.YO,data.latO,data.lonO
   
   #ZZM = mask_region(ZM,region)
   #ZZO = mask_region(ZO,region)

   #miz = ZZM[(ZZM>.15)&(ZZM<.8)].size
   #sie = ZZM[.15<ZZM].size
   #nans = np.isnan(ZZM)
   #tot = mrg[mrg==2].size

   #omiz = ZZO[(ZZO>.15)&(ZZO<.8)].size
   #osie = ZZO[.15<ZZO].size
   #onans = np.isnan(ZZO)
   #otot = org[org==2].size

   #filname = './'+str(region)+'_ext.txt'
   #with open(filname,'a') as f:
   #   row1 = [dadate,tot,sie,miz,otot,osie,omiz]
   #   str1 = ' '.join(map(str,row1))
   #   f.write(str1+'\n')
   #   f.close()

   print('Analysing model MIZ (widths)...')
   bar_mdl_MIZ = MIZwidth(XM,YM,ZM,bmm,'bar')
   print('DONE')
   print('')
   
   print('Analysing osisaf MIZ (widths)...')
   bar_osi_MIZ = MIZwidth(XO,YO,ZO,bmo,'bar')
   print('DONE')
   print('')
   
elif region == 'gre':
   print('Studying the GREENLAND region...')
   print('')
   
   ZM,XM,YM,latM,lonM,ZO,XO,YO,latO,lonO = data.ZM,data.XM,data.YM,data.latM,data.lonM,\
           data.ZO,data.XO,data.YO,data.latO,data.lonO
   
   #ZZM = mask_region(ZM,region)
   #ZZO = mask_region(ZO,region)

   #miz = ZZM[(ZZM>.15)&(ZZM<.8)].size
   #sie = ZZM[.15<ZZM].size
   #nans = np.isnan(ZZM)
   #tot = mrg[mrg==1].size

   #omiz = ZZO[(ZZO>.15)&(ZZO<.8)].size
   #osie = ZZO[.15<ZZO].size
   #onans = np.isnan(ZZO)
   #otot = org[org==1].size

   #filname = './'+str(region)+'_ext.txt'
   #with open(filname,'a') as f:
   #   row1 = [dadate,tot,sie,miz,otot,osie,omiz]
   #   str1 = ' '.join(map(str,row1))
   #   f.write(str1+'\n')
   #   f.close()

   print('Analysing model MIZ (widths)...')
   gre_mdl_MIZ = MIZwidth(XM,YM,ZM,bmm,'gre')
   print('DONE')
   print('')
   
   print('Analysing osisaf MIZ (widths)...')
   gre_osi_MIZ = MIZwidth(XO,YO,ZO,bmo,'gre')
   print('DONE')
   print('')
   
elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime
