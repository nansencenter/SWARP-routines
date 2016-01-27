############################################################################
# THIS SCRIPT WILL GIVE ARTIC,REGIONAL or FULL INFO ABOUT ICE EDGES
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

sys.path.append('../../py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution_2tw as Leqs
import MIZchar as widths

# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')

############################################################################
#CLASSES
############################################################################
# This Class will read in and prepare every file for polygon detection
class reader:
    def __init__(self,name,dadate,basemap):
        self.year = dadate[:4]
        self.month = dadate[4:6]
        self.day = dadate[6:8]
        gigio = datetime.datetime(int(float(self.year)),int(float(self.month)),int(float(self.day)),0,0)
        gigio = gigio.strftime('%j')
        gigio = int(float(gigio))
        self.julian = gigio - 1
        self.filname = name+'_'+dadate
        if name == 'Osisaf':
            self.lonO,self.latO,self.X,self.Y,self.Z = self._read_osi_(dadate,basemap) 
            return
        elif name == 'Model':
            self.lonM,self.latM,self.X,self.Y,self.ZC,self.ZD = self._read_mdl_(dadate,basemap)
            return
        elif name == 'Aari':
            self.ad = self._read_aari_(dadate,basemap)
    
    def _read_osi_(self,dadate,basemap):
        day = self.day
        month = self.month
        year = self.year
        # Read in OSI_SAF file
        #outdir = '/work/shared/nersc/msc/OSI-SAF/'+str(year)+'_nh_polstere/'
        outdir = './tmp/OSI'
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
        Z2 = conc[:,:].data
        mask2 = conc[:,:].mask
        Z2[mask2] = np.NaN
        ZO = Z2/100
        return(lonO,latO,XO,YO,ZO) 
   
    def _read_mdl_(self,dadate,basemap):
    
        #TODO
        # work on the binary reader, for now *.nc will be used 
    
        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        # Read TP4arch_wav
        #outdir = '/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/'+dadate+'/netcdf/'
        outdir = './tmp/MDL'
        ncfil = outdir+'/TP4archv_wav_start'+str(dadate)+'_000000Z_dump'+str(dadate)+'_120000Z.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'fice'
        sdmax = 'dmax'
        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=0)
        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
        ZD = dmax[:,:].data
        mask = dmax[:,:].mask
        ZD[mask] = np.NaN
        ZC = conc[:,:].data
        mask = conc[:,:].mask
        ZC[mask] = np.NaN
        
        return(lonM,latM,X,Y,ZC,ZD)

############################################################################
# This Class will find AOD polygons for 2 sets with different grid.
class AOD_poly:
    def __init__(self,X1,Y1,Z1,X2,Y2,Z2):
        ZM = Z1
        ZO = Z2
        # NOTE THRESHOLD FOR ICE EDGE VALIDATION
        self.DN,self.B1,self.B2,oDN,uDN = self.binary_diff(ZM,ZO,.15)
        self.over,self.under = self.poly_maker(self.DN)
        self.DN = self.mapper(self.DN,ZM,ZO)
        return
    
    def binary_diff(self,data1,data2,thresh):
        # NOTE data1 is MODEL, data2 is OSISAF
        # This gives the difference as:
        # -1 = OSI PRES, MDL MISS -> UNDERESTIMATION
        #  0 = OSI PRES, MDL PRES (OR OSI MISS, MDL MISS) -> HIT (NO-HIT)
        # +1 = OSI MISS, MDL PRES -> OVERESTIMATION
        def closing(DN):
            # apply closing to avoid small polynias and clean up a little
            kernel = np.ones((3,3),np.uint8)
            DN_cl = morph.closing(DN,kernel)
            return(DN_cl)
        
        # NOTE outputs are difference, overpredictions and underpredictions (CLOSED & MAPPED)
        
        # generating the binary file, get the difference, divide into overprediction and
        # underprediction, apply the closing on both and finally re-map the land
        ndata1 = np.copy(data1)
        ndata2 = np.copy(data2)

        # OLD method where ice-ice = ocean-ocean = 0
        ndata1[ndata1<thresh] = 0
        ndata2[ndata2<thresh] = 0
        ndata1[ndata1>=thresh] = 1  
        ndata2[ndata2>=thresh] = 1  

        # Difference MODEL - OSISAF
        ddata = ndata1 - ndata2
        # Cut out the nans for both
        thenans = np.isnan(ddata)
        ddata[thenans] = 0
        # Work with over and under (clean them up)
        over_data = np.copy(ddata)
        under_data = np.copy(ddata)
        over_data[over_data==-1] = 0
        under_data[under_data==1] = 0
        under_data = abs(under_data)
        # Closing operation (see def closing)
        o_data = closing(over_data)
        u_data = closing(under_data)
        # Add the under to the over to get total diff.
        o_data[u_data==1] = -1
        ddata = np.copy(o_data)
        return(ddata,ndata1,ndata2,o_data,u_data)
    
    def poly_maker(self,ddata):
        # finding contours from the difference data map
        over = msr.find_contours(ddata,.5)
        over = sorted(over, key=len)
        over = over[::-1]
        under = msr.find_contours(ddata,-.5)
        under = sorted(under, key=len)
        under = under[::-1]
        return(over,under)
    
    def mapper(self,DIFF,INS,OUT):
        # getting back the terrain lost during binary_mod (for contour finding reasons)
        DN = DIFF
        ZO = OUT
        ZI = INS
        thenan = np.isnan(ZO)
        thenan2 = np.isnan(ZI)
        DN[thenan] = None
        DN[thenan2] = None
        return(DN)

############################################################################
# This Class will analyse contours and adds polygon's info to the daily regional file
class poly_stat:
    def __init__(self,cont,poly_status,OUT,INS,DIFF,X,Y,number,basemap=None):
        self.number = number
        self.polygon_status = poly_status
        self.ij_list = cont
        xy_list = self.ij2xy(cont,X,Y)
        self.difference = DIFF
        f_vals = self.split_cont(OUT,INS,poly_status)
        lon,lat = self.lon_lat(xy_list,basemap)
        self.area_perimeter()
        if run == 'arctic':
            self.add2list_arctic(number,lon,lat,f_vals,poly_status)
        else:
            self.add2list_regional(number,lon,lat,f_vals)

        ## OLD CLASSI DEFINITION, not needed for now
        ## class definition
        #if len(cont) > 200:
        #    self.polygon_class = 'H'
        #elif len(cont) > 100 and len(cont) <= 200:
        #    self.polygon_class = 'B'
        #elif len(cont) > 30 and len(cont) <= 100:
        #    self.polygon_class = 'M'
        #elif len(cont) <= 30:
        #    self.polygon_class = 'S'

    # output: xy_coords
    # NOTE weird error displayed (non-integer index?) it works but I need to check this
    def ij2xy(self,cont,X,Y):
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
        self.xy_list = xy_list
        return(xy_list)
    
    # output: f_vals
    def split_cont(self,OUT,INS,polygon_status):
        # this function find the different contours
        # NOTE if a point has non integer i coordinate is going to be a vertical edge,
        # if has non integer j coordinate is going to be a horizontal edge hence the use
        # of different arounds depending on the indexes of the point
        # NOTE if the polygon is an OVERESTIMATION the contour finding is inverted
        vs = self.ij_list # list of (i,j) pixel indices
        in_cont = []
        out_cont = []
        unk_cont = []
        func_vals= []
        func_mod = 0 # value of func_vals at model ice edge
        func_osi = 1 # value of func_vals at OSISAF ice edge
        func_unk = 2 # value of func_vals at other parts of contour
        
        if polygon_status == 0:
            for n,el in enumerate(vs):
                #getting all the neighbours
                around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
                around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
                check_cont = 0
                if el[0]/int(el[0]) == 1:
                    for h,v in around1:
                        if OUT[h][v] == INS[h][v] == 1:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif OUT[h][v] == INS[h][v] == 0:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
                else:
                    for h,v in around2:
                        if OUT[h][v] == INS[h][v] == 1:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif OUT[h][v] == INS[h][v] == 0:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
        else:
            for n,el in enumerate(vs):
                #getting all the neighbours
                around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
                around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
                check_cont = 0
                if el[0]/int(el[0]) == 1:
                    for h,v in around1:
                        if OUT[h][v] == INS[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif OUT[h][v] == INS[h][v] == 1:
                            out_cont.append(el)
                            func_val=func_osi
                            check_cont = 1
                    if check_cont == 0:
                        unk_cont.append(el)
                        func_val=func_unk
                    func_vals.append(func_val)
                else:
                    for h,v in around2:
                        if OUT[h][v] == INS[h][v] == 0:
                            in_cont.append(el)
                            func_val=func_mod
                            check_cont = 1
                        elif OUT[h][v] == INS[h][v] == 1:
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
        self.in_ice_edge = in_cont
        self.out_ice_edge = out_cont
        self.unknown_edge = unk_cont
        self.f_vals = f_vals
        return(f_vals)
    
    # output: lon,lat,region of polygon
    def lon_lat(self,xy_list,basemap):
        # Centroid - Location in lon/lat of the centroid of the polygon
        DX = 0
        DY = 0
        B = 0
        for n in range(len(xy_list)-1):
            DX += ((xy_list[n][0]+xy_list[n+1][0])*\
                  ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
            DY += ((xy_list[n][1]+xy_list[n+1][1])*\
                  ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
            B += ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1]))
        A = (0.5)*B 
        CX = (1/float(6*A))*DX
        CY = (1/float(6*A))*DY
        self.centroid = [CX,CY]
        centroid_longitude,centroid_latitude = basemap(CX,CY,inverse=True)
        self.centroid_longitude = centroid_longitude
        self.centroid_latitude = centroid_latitude
        lon_list,lat_list=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
        
        # Defining the region
        if centroid_longitude < 90 and centroid_longitude > 16:
            self.region = 'Barents_Kara_sea'
        elif centroid_longitude <= 16 and centroid_longitude > -44:
            self.region = 'Greenland_sea'
        elif centroid_longitude <= -44 and centroid_longitude > -90:
            n = 0
            N = len(bblatw[:-1])
            if centroid_latitude <= bblatw[-1] and centroid_longitude >= bblonw[-1]:
                self.region = 'Labrador_sea'
            else:
                while n < N:
                    if centroid_latitude >= bblatw[n+1] and centroid_latitude <= bblatw[n]:
                        if centroid_longitude >= bblonw[n]: 
                            self.region = 'Labrador_sea'
                            break
                    n += 1
                else:
                    self.region = 'North_Canada_Beaufort_sea'
        elif centroid_longitude < 180 and centroid_longitude > 90:
            self.region = 'Laptev_East_Siberian_sea'
        else:
            self.region = 'North_Canada_Beaufort_sea'
       
        return(lon_list,lat_list)

    # output: add the polygon points to the regional list
    def add2list_arctic(self,number,lon,lat,fvals,status):
        if status == '0':
            for n, en in enumerate(fvals):
                under_arctic_list.append([number,lon[n],lat[n],en])
        else:
            for n, en in enumerate(fvals):
                over_arctic_list.append([number,lon[n],lat[n],en])

    # output: add the polygon points to the regional list
    def add2list_regional(self,number,lon,lat,fvals):
        region = self.region
        if region == 'Barents_Kara_sea':
            for n, en in enumerate(fvals):
                bar_poly_stat.append([number,lon[n],lat[n],en])
        elif region == 'Greenland_sea':
            for n, en in enumerate(fvals):
                gre_poly_stat.append([number,lon[n],lat[n],en])
        elif region == 'Labrador_sea':
            for n, en in enumerate(fvals):
                lab_poly_stat.append([number,lon[n],lat[n],en])
        elif region == 'North_Canada_Beaufort_sea':
            for n, en in enumerate(fvals):
                ncb_poly_stat.append([number,lon[n],lat[n],en])
        elif region == 'Laptev_East_Siberian_sea':
            for n, en in enumerate(fvals):
                les_poly_stat.append([number,lon[n],lat[n],en])

    # output: euclidean area (E_area) and perimeter (E_perim)
    def area_perimeter(self):
        # Calculating the area of irregular polygon (from perimeter)
        vs = self.ij_list
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
        return
    
############################################################################
# FUNCTIONS
############################################################################
# Function that gives arctic info about hits, over, under and neg-hits
def arctic_stats(data1,data2,thresh=.15):

    # NOTE data1 is MODEL, data2 is OSISAF
    # This gives the difference as:
    # -1 = OSI PRES, MDL MISS -> UNDERESTIMATION
    #  0 = OSI PRES, MDL PRES (OR OSI MISS, MDL MISS) -> HIT (NO-HIT)
    # +1 = OSI MISS, MDL PRES -> OVERESTIMATION
    def closing(DN):
        # apply closing to avoid small polynias and clean up a little
        kernel = np.ones((3,3),np.uint8)
        DN_cl = morph.closing(DN,kernel)
        return(DN_cl)
    
    # NOTE outputs are difference, overpredictions and underpredictions (CLOSED & MAPPED)
    
    # generating the binary file, get the difference, divide into overprediction and
    # underprediction, apply the closing on both and finally re-map the land
    ndata1 = np.copy(data1)
    ndata2 = np.copy(data2)
    
    # NEW method where ice-ice = 2 ; ocean-ocean = -2
    ndata1[ndata1>=thresh] = 1
    ndata1[ndata1<thresh] = 2
    ndata2[ndata2<thresh] = 0
    ndata2[ndata2>=thresh] = 3
    
    # Difference MODEL - OSISAF
    ddata = ndata1 - ndata2
    # Cut out the nans for both
    thenans = np.isnan(ddata)
    ddata[thenans] = 0
    # Work with over and under (clean them up)
    over_data = np.copy(ddata)
    under_data = np.copy(ddata)
    over_data[over_data==-1] = 0
    under_data[under_data==1] = 0
    under_data = abs(under_data)
    # Closing operation (see def closing)
    o_data = closing(over_data)
    u_data = closing(under_data)
    # Add the under to the over to get total diff.
    o_data[u_data==1] = -1
    ddata = np.copy(o_data)

    total = ddata.size
    tot_hits = ddata[ddata==2].size
    tot_over = ddata[ddata==1].size
    tot_under = ddata[ddata==-1].size
    tot_neg_hits = ddata[ddata==-2].size
    tot_unk = ddata[ddata==0].size

    return(total,tot_hits,tot_over,tot_under,tot_neg_hits,tot_unk)

############################################################################
# Function that gives regional info about hits, over, under and neg-hits
def regional_stats(X,Y,DN,ZM,ZO,bm):
    lon,lat = bm(X[:,:],Y[:,:],inverse=True)
    
    # Regional Lists
    bar_lst = []
    gre_lst = []
    lab_lst = []
    ncb_lst = []
    les_lst = []


    # Going through every point to assess status and location
    for i in range(len(DN[:,0])):
        for j in range(len(DN[0,:])):
            ptlon = lon[i,j]
            ptlat = lat[i,j]
            if ptlon < 90 and ptlon > 16:
                region = 'Barents_Kara_sea'
            elif ptlon < 16 and ptlon > -44:
                region = 'Greenland_sea'
            elif ptlon < -44 and ptlon > -90:
                n = 0
                N = len(bblatw[:-1])
                if ptlat <= bblatw[-1] and ptlon >= bblonw[-1]:
                    region = 'Labrador_sea'
                else:
                    while n < N:
                        if ptlat >= bblatw[n+1] and ptlat <= bblatw[n]:
                            if ptlon >= bblonw[n]: 
                                region = 'Labrador_sea'
                                break
                        n += 1
                    else:
                        region = 'North_Canada_Beaufort_sea'
            elif ptlon < 180 and ptlon > 90:
                region = 'Laptev_East_Siberian_sea'
            else:
                region = 'North_Canada_Beaufort_sea'

            # LIST as [ hit, over, under, neg-hit, unkwn ]
            if region == 'Barents_Kara_sea':
                if DN[i,j] == 2:
                    bar_lst.append([1,0,0,0,0])
                elif DN[i,j] == 1:
                    bar_lst.append([0,1,0,0,0])
                elif DN[i,j] == -1:
                    bar_lst.append([0,0,1,0,0])
                elif DN[i,j] == -2:
                    bar_lst.append([0,0,0,1,0])
                elif DN[i,j] == 0:
                    bar_lst.append([0,0,0,0,1])
            elif region == 'Greenland_sea':
                if DN[i,j] == 2:
                    gre_lst.append([1,0,0,0,0])
                elif DN[i,j] == 1:
                    gre_lst.append([0,1,0,0,0])
                elif DN[i,j] == -1:
                    gre_lst.append([0,0,1,0,0])
                elif DN[i,j] == -2:
                    gre_lst.append([0,0,0,1,0])
                elif DN[i,j] == 0:
                    gre_lst.append([0,0,0,0,1])
            elif region == 'Labrador_sea':
                if DN[i,j] == 2:
                    lab_lst.append([1,0,0,0,0])
                elif DN[i,j] == 1:
                    lab_lst.append([0,1,0,0,0])
                elif DN[i,j] == -1:
                    lab_lst.append([0,0,1,0,0])
                elif DN[i,j] == -2:
                    lab_lst.append([0,0,0,1,0])
                elif DN[i,j] == 0:
                    lab_lst.append([0,0,0,0,1])
            elif region == 'North_Canada_Beaufort_sea':
                if DN[i,j] == 2:
                    ncb_lst.append([1,0,0,0,0])
                elif DN[i,j] == 1:
                    ncb_lst.append([0,1,0,0,0])
                elif DN[i,j] == -1:
                    ncb_lst.append([0,0,1,0,0])
                elif DN[i,j] == -2:
                    ncb_lst.append([0,0,0,1,0])
                elif DN[i,j] == 0:
                    ncb_lst.append([0,0,0,0,1])
            elif region == 'Laptev_East_Siberian_sea':
                if DN[i,j] == 2:
                    les_lst.append([1,0,0,0,0])
                elif DN[i,j] == 1:
                    les_lst.append([0,1,0,0,0])
                elif DN[i,j] == -1:
                    les_lst.append([0,0,1,0,0])
                elif DN[i,j] == -2:
                    les_lst.append([0,0,0,1,0])
                elif DN[i,j] == 0:
                    les_lst.append([0,0,0,0,1])      
    
    return(bar_lst,gre_lst,lab_lst,ncb_lst,les_lst)

############################################################################
# Function stat chart
def stat_chart(self,save=False):
    # The Statistical Chart is an output that tries to compress as many informations
    # as possible in a single figure.
    # The figure is:
    # 1) Top left arctic map with studied polygon highlighted
    # 2) Top right a close up of the polygon with contours and min. distances
    # 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
    # 4) Bottom right is a legend for the small image
    
    nmbr = self.number
    pname = 'Polygon_'+str(nmbr)
    self.name = pname
    pclass = self.polygon_class
    if self.polygon_status == 1:
        pstat = 'Overestimation'
    else:
        pstat = 'Underestimation'
    region = self.region
    print ''
    print 'Statistic Chart for ',pname
    print 'Class and Region ',pclass,region
    print 'Status ',pstat
    DN = self.difference
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
    pos_p = mpatches.Patch(color='lightgreen')
    neg_p = mpatches.Patch(color='royalblue')
    leg.legend([mc,oc,uc,pos_p,neg_p],(\
          'Model Cont.','Observation Cont.','Unknown Cont.', \
          'Model Underestimate','Model Overestimate'),loc='center')
             
    # Statistics text on the bottom left
    txt = '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
          '2) Area = '+str(area)+' km^2\n'+ \
          '3) Perimeter = '+str(perim)+' km\n'+ \
          '4) Status = '+str(pstat)
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
    return

############################################################################
# Function that reprojects model into observational grid
# NOTE the basemap used in the reprojector is based on OSISAF's grid 
# Theoretically any stereographical map could be reprojected into this grid.
def reproj(X1,Y1,Z1,X2,Y2):

    # getting ready for reprojection
    X3 = X1.reshape(X1.size)
    Y3 = Y1.reshape(Y1.size)
    Z3 = Z1.reshape(Z1.size)
    C = [X3,Y3]
    C = np.array(C)
    C = C.T
    
    # Interpolation can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
    ZN = grd(C,Z3,(X2,Y2),method='linear')
    
    return(ZN)

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
    return(tot,avg)

############################################################################
# GEO INFO
############################################################################
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

###########################################################################
# DIFFERENT RUNS
###########################################################################
# ARTIC RUN 
def arctic(XO,YO,ZO,ZM,DN,over,under,basemap):

    # daily regional gen_stats
    gen_arctic_stats = arctic_stats(ZM,ZO)
    
    # Prep the lists for polygon.txt
    global over_arctic_list
    over_arctic_list = [['Polygon_number','Longitude','Latitude','Func_value']]
    global under_arctic_list
    under_arctic_list = [['Polygon_number','Longitude','Latitude','Func_value']]

    poly_list = []
    over_width_list = []
    under_width_list = []
    
    # Outpud directory
    outdir = str(out_dir)+'/'+str(dadate)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    reg_repo = outdir+'/arctic'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)

    # Polygon analysis
    for n,el in enumerate(over):
        # classification of overestimation polygons
        aod=poly_stat(el,'1',BO,BM,DN,XO,YO,n,basemap=basemap)
        poly_list.append(aod)

    # Writing the overestimation polygons list
    filname = str(reg_repo)+'/over_polygons.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(over_arctic_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    # Width calculation for OVERESTIMATION
    over_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(over_widths):
        if wd is not None:
            over_width_list.append(wd.int_widths)

    filname = str(reg_repo)+'/over_widths.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(over_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    try:
        n
    except NameError:
        n = -1 
    for n2,el in enumerate(under):
        # classification of negative (n+1 for good enumeration)
        n += 1 
        aod=poly_stat(el,'0',BO,BM,DN,XO,YO,n,basemap=basemap)
        poly_list.append(aod)
        
    # Writing the overestimation polygons list
    filname = str(reg_repo)+'/under_polygons.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(under_arctic_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    # Width calculation for UNDERESTIMATION
    under_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(under_widths):
        if wd is not None:
            under_width_list.append(wd.int_widths)   

    filname = str(reg_repo)+'/under_widths.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(under_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    outdir = str(out_dir)+'/'+str(dadate)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    reg_repo = outdir+'/arctic'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)

    filname = str(reg_repo)+'/daily_stats.txt'
    with open(filname,'a') as f:
        row1 = ['Tot.Points','IceHits','Over','Under','SeaHits','Unknown','Avg.OverWidths','Avg.UnderWidths']
        str1 = ' '.join(map(str,row1))
        row2 = [gen[0],gen[1],gen[2],gen[3],gen[4],gen[5],oavg,uavg]
        str2 = ' '.join(map(str,row2))
        f.write(str1+'\n'+str2)
        f.close()

    return(gen_arctic_stats,over_width_list,under_width_list)
        
###########################################################################
# REGIONAL RUN
def regional(XO,YO,ZO,ZM,DN,over,under,basemap):

    # daily regional gen_stats
    bar_lst,gre_lst,lab_lst,ncb_lst,les_lst = gen_region_stats(XO,YO,DN,ZM,ZO,basemap)
    
    # Prep the lists for polygon.txt
    global bar_poly_stat
    bar_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    global gre_poly_stat
    gre_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    global lab_poly_stat
    lab_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    global les_poly_stat
    les_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    global ncb_poly_stat
    ncb_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]

    poly_list=[]
    
    STCHidx = None

    for n,el in enumerate(over):
        # classification of positive polygons
        aod=poly_stat(el,'1',BO,BM,DN,XO,YO,n,basemap=hqm,STCH=STCHidx)
        poly_list.append(aod)
    try:
        n
    except NameError:
        n = -1 
    for n2,el in enumerate(under):
        # classification of negative (n+1 for good enumeration)
        n += 1 
        aod=poly_stat(el,'0',BO,BM,DN,XO,YO,n,basemap=hqm,STCH=STCHidx)
        poly_list.append(aod)
    
    outdir = str(out_dir)+'/'+str(dadate)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    reg_repo = outdir+'/Barents_Kara_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(bar_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    filname = str(reg_repo)+'/polygons.txt'
    results = aod_widths.single_file(filname,basemap)
    

    reg_repo = outdir+'/Greenland_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(gre_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    
    reg_repo = outdir+'/Labrador_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(lab_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    reg_repo = outdir+'/North_Canada_Beaufort_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(ncb_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    reg_repo = outdir+'/Laptev_East_Siberian_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(les_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

###########################################################################
# FULL RUN
def full(XO,YO,ZO,ZM,DN,over,under,basemap):

    # daily regional gen_stats
    bar_lst,gre_lst,lab_lst,ncb_lst,les_lst = gen_region_stats(XO,YO,DN,ZM,ZO,basemap)
    
    # Prep the lists for polygon.txt
    bar_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    gre_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    lab_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    les_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]
    ncb_poly_stat = [['Polygon_number','Longitude','Latitude','Func_value']]

    poly_list=[]
    
    STCHidx = None

    for n,el in enumerate(over):
        # classification of positive polygons
        aod=poly_stat(el,'1',BO,BM,DN,XO,YO,n,basemap=hqm,STCH=STCHidx)
        poly_list.append(aod)
    try:
        n
    except NameError:
        n = -1 
    for n2,el in enumerate(under):
        # classification of negative (n+1 for good enumeration)
        n += 1 
        aod=poly_stat(el,'0',BO,BM,DN,XO,YO,n,basemap=hqm,STCH=STCHidx)
        poly_list.append(aod)
    
    outdir = str(out_dir)+'/'+str(dadate)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    reg_repo = outdir+'/Barents_Kara_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(bar_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    filname = str(reg_repo)+'/polygons.txt'
    results = aod_widths.single_file(filname,basemap)
    

    reg_repo = outdir+'/Greenland_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(gre_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    
    reg_repo = outdir+'/Labrador_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(lab_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    reg_repo = outdir+'/North_Canada_Beaufort_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(ncb_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

    reg_repo = outdir+'/Laptev_East_Siberian_sea'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    
    with open(str(reg_repo)+'/polygons.txt','a') as f:
        for l,el in enumerate(les_poly_stat):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()

###########################################################################
# Beginning of the script
###########################################################################

###########################################################################
# Defining the basemap
# NOTE Based on OSISAF's grid
# low quality map
lqm = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),\
        projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
# better quality map
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
        projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
###########################################################################

# MISSION and DATE
run = sys.argv[1]
dadate = sys.argv[2] 
out_dir = sys.argv[3]

# Start the count and analyse the datasets
time0 = time.time()

# Read in the data
osisaf = reader('Osisaf',dadate,hqm)
model = reader('Model',dadate,hqm)
lonO,latO,XO,YO,ZO = osisaf.lonO,osisaf.latO,osisaf.X,osisaf.Y,osisaf.Z
lonM,latM,XM,YM,ZMM = model.lonM,model.latM,model.X,model.Y,model.ZC

# Reproject model data
ZM = reproj(XM,YM,ZMM,XO,YO)

# Analyse the AODs
AOD = AOD_poly(XM,YM,ZM,XO,YO,ZO)
over,under,DN,BM,BO = AOD.over,AOD.under,AOD.DN,AOD.B1,AOD.B2

# Get custom geo info
bblone,bblate,bblonw,bblatw = baffin_bay()
       
if run == 'arctic':
    gen,o_list,u_list = arctic(XO,YO,ZO,ZM,DN,over,under,lqm)
    otot,oavg = widths_stats(o_list)
    utot,uavg = widths_stats(u_list)

elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime
