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

#sys.path.append('./py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution_2tw as Leqs

# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#CLASSES
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

############################################################################
# This Class will read in and prepare every file for polygon detection
class reader:
    def __init__(self,name,dadate,basemap,jday):
        self.year = dadate[:4]
        self.month = dadate[4:6]
        self.day = dadate[6:8]
        self.jday = jday
        gigio = datetime.datetime(int(float(self.year)),int(float(self.month)),int(float(self.day)),0,0)
        gigio = gigio.strftime('%j')
        gigio = int(float(gigio))
        self.julian = gigio - 1
        self.filname = name+'_'+dadate
        if name == 'Osisaf':
            self.lonO,self.latO,self.X,self.Y,self.Z = self._read_osi_(dadate,basemap) 
            return
        elif name == 'Model':
            #self.lonM,self.latM,self.X,self.Y,self.ZC,self.ZD = self._read_mdl_(dadate,basemap)
            self.lonM,self.latM,self.X,self.Y,self.ZC = self._read_mdl_(dadate,basemap)
            return
        elif name == 'Aari':
            self.ad = self._read_aari_(dadate,basemap)
    
    def _read_osi_(self,dadate,basemap):
        day = self.day
        month = self.month
        year = self.year
        # Read in OSI_SAF file
        #outdir = '/work/shared/nersc/msc/OSI-SAF/'+str(year)+'_nh_polstere/'
        outdir = '../osisaf_repo'
        ncfil = outdir+'/ice_conc_nh_polstere-100_reproc_'+dadate+'1200.nc'
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
        jday = self.jday
        # Read TP4arch_wav
        outdir = '../dailync'
        ncfil = outdir+'/TP4DAILY_'+year+'_'+jday+'.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'fice00'
        #sdmax = 'dmax'
        lonM = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        latM = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
        #dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=0)
        X,Y = basemap(lonM[:,:],latM[:,:],inverse=False)
        #ZD = dmax[:,:].data
        #mask = dmax[:,:].mask
        #ZD[mask] = np.NaN
        ZC = conc[:,:].data
        mask = conc[:,:].mask
        ZC[mask] = np.NaN
        
        return(lonM,latM,X,Y,ZC) #,ZD)

############################################################################
# This Class will find MIZ contours
class MIZ_poly:
    def __init__(self,X1,Y1,Z1,typo):
        ZM = np.copy(Z1)
        if typo == 'ic':
            self.MIZraw,self.B1,self.B2 = self.binary_diff(ZM,.15,.80)
        elif typo == 'fsd':
            self.MIZraw,self.B1,self.B2 = self.binary_diff(ZM,.1,300)
        self.MIZcont = self.poly_maker(self.MIZraw)
        self.MIZ = self.mapper(self.MIZraw,ZM,typo)
        return
    
    def binary_diff(self,data,threshm,threshM):
        def closing(data):
            # apply closing to avoid small polynias and clean up a little
            kernel = np.ones((3,3),np.uint8)
            data_cl = morph.closing(data,kernel)
            return(data_cl)
        
        # generating the binary file, get the difference, divide into overprediction and
        # underprediction, apply the closing on both and finally re-map the land
        ndata1 = np.copy(data)
        ndata2 = np.copy(data)

        # OLD method where ice-ice = ocean-ocean = 0
        ndata1[ndata1<threshm] = 0
        ndata2[ndata2<threshM] = 0
        ndata1[ndata1>=threshm] = 1  
        ndata2[ndata2>=threshM] = 1  

        # Difference MODEL - OSISAF
        ddata = ndata1 - ndata2
        # Cut out the nans for both
        thenans = np.isnan(ddata)
        ddata[thenans] = 0
        # Closing operation (see def closing)
        ddata = closing(ddata)
        # Add the under to the over to get total diff.
        return(ddata,ndata1,ndata2)
    
    def poly_maker(self,ddata):
        # finding contours from the difference data map
        MIZcont = msr.find_contours(ddata,.5)
        MIZcont = sorted(MIZcont, key=len)
        MIZcont = MIZcont[::-1]
        return(MIZcont)
    
    def mapper(self,DIFF,Z1,typo):
        DN = np.copy(DIFF)
        IP = np.copy(Z1)
        ZM = np.copy(Z1)
        thenan = np.isnan(ZM)
        if typo == 'ic':
            IP[IP>=.80] = 2
            IP[IP<.80] = 0
        elif typo == 'fsd':
            IP[IP==300] = 2
            IP[IP<300] = 0
        #getting back the ice pack
        DN = IP + DN
        DN[DN>2] = 1
        # getting back the terrain lost during binary_mod (for contour finding reasons)
        DN[thenan] = None
        return(DN)

############################################################################
# This Class will analyse contours and adds polygon's info to the daily regional file
class poly_stat:
    def __init__(self,cont,OUT,INS,IPMIZ,X,Y,number,basemap=None):
        self.number = number #get the number
        self.ij_list = cont #get the contour
        self.class_def(cont) #calculate class from contour
        xy_list = self.ij2xy(cont,X,Y) #get xy from ij (cont)
        self.ip_miz = IPMIZ #get ip_miz
        f_vals = self.split_cont(OUT,INS) #
        lon,lat = self.lon_lat(xy_list,basemap)
        self.area_perimeter()
        self.add2list(self.region,number,lon,lat,f_vals)

    # class definition
    def class_def(self,contour):
        lcont = len(contour)
        if lcont > 200:
            self.polygon_class = 'H'
        elif lcont > 100 and lcont <= 200:
            self.polygon_class = 'B'
        elif lcont > 30 and lcont <= 100:
            self.polygon_class = 'M'
        elif lcont <= 30:
            self.polygon_class = 'S'
        return()

        #self.stat_chart(save=True)

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
    def split_cont(self,OUT,INS):
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

        # let's try the above method, if it doesn't work, use bot
        #for n,el in enumerate(vs):
        #    #getting all the neighbours
        #    around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1])) # vertical boundaries - OK!
        #    around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5)) # horizontal boundaries
        #    check_cont = 0
        #    if el[0]/int(el[0]) == 1:
        #        for h,v in around1:
        #            if OUT[h][v] == INS[h][v] == 0:
        #                in_cont.append(el)
        #                func_val=func_mod
        #                check_cont = 1
        #            elif OUT[h][v] == INS[h][v] == 1:
        #                out_cont.append(el)
        #                func_val=func_osi
        #                check_cont = 1
        #        if check_cont == 0:
        #            unk_cont.append(el)
        #            func_val=func_unk
        #        func_vals.append(func_val)
        #    else:
        #        for h,v in around2:
        #            if OUT[h][v] == INS[h][v] == 0:
        #                in_cont.append(el)
        #                func_val=func_mod
        #                check_cont = 1
        #            elif OUT[h][v] == INS[h][v] == 1:
        #                out_cont.append(el)
        #                func_val=func_osi
        #                check_cont = 1
        #        if check_cont == 0:
        #            unk_cont.append(el)
        #            func_val=func_unk
        #        func_vals.append(func_val)
        
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
    # NOTE cont_list might be XY or IJ (no LONLAT)
    def lon_lat(self,cont_list,basemap,typo='xy'):

        # Centroid - Location in lon/lat of the centroid of the polygon
        DX = 0
        DY = 0
        B = 0
        for n in range(len(cont_list)-1):
            DX += ((cont_list[n][0]+cont_list[n+1][0])*\
                  ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1])))  
            DY += ((cont_list[n][1]+cont_list[n+1][1])*\
                  ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1])))  
            B += ((cont_list[n][0]*cont_list[n+1][1])-(cont_list[n+1][0]*cont_list[n][1]))
        A = (0.5)*B 
        CX = (1/float(6*A))*DX
        CY = (1/float(6*A))*DY
        self.centroid = [CX,CY]
        if typo == 'xy':
            centroid_longitude,centroid_latitude = basemap(CX,CY,inverse=True)
            self.centroid_longitude = centroid_longitude
            self.centroid_latitude = centroid_latitude
            lon_list,lat_list=basemap(cont_list[:,0],cont_list[:,1],inverse=True)
            bblone,bblate,bblonw,bblatw = baffin_bay()
            # Defining the region
            if centroid_longitude < 90 and centroid_longitude > 16:
                self.region = 'bar'
            elif centroid_longitude <= 16 and centroid_longitude > -44:
                self.region = 'gre'
            elif centroid_longitude <= -44 and centroid_longitude > -90:
                n = 0
                N = len(bblatw[:-1])
                if centroid_latitude <= bblatw[-1] and centroid_longitude >= bblonw[-1]:
                    self.region = 'lab'
                else:
                    while n < N:
                        if centroid_latitude >= bblatw[n+1] and centroid_latitude <= bblatw[n]:
                            if centroid_longitude >= bblonw[n]: 
                                self.region = 'lab'
                                break
                        n += 1
                    else:
                        self.region = 'ncb'
            elif centroid_longitude < 180 and centroid_longitude > 90:
                self.region = 'les'
            else:
                self.region = 'ncb'
        elif typo == 'ij':
            regions = np.load('./geo_info/regions.npy')
            if regions[CX,CY] == 1:
                self.region = 'bar'
            elif regions[CX,CY] == 2:
                self.region = 'les'
            elif regions[CX,CY] == 3:
                self.region = 'ncb'
            elif regions[CX,CY] == 4:
                self.region = 'lab'
            elif regions[CX,CY] == 5:
                self.region = 'gre'
       
        return(lon_list,lat_list)

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

    # output: add the polygon points to the regional list
    def add2list(self,region,number,lon,lat,fvals):

        if region == 'bar':
            for n, en in enumerate(fvals):
                bar_poly_list.append([number,lon[n],lat[n],en])
        elif region == 'les':
            for n, en in enumerate(fvals):
                les_poly_list.append([number,lon[n],lat[n],en])
        elif region == 'ncb':
            for n, en in enumerate(fvals):
                ncb_poly_list.append([number,lon[n],lat[n],en])
        elif region == 'lab':
            for n, en in enumerate(fvals):
                lab_poly_list.append([number,lon[n],lat[n],en])
        elif region == 'gre':
            for n, en in enumerate(fvals):
                gre_poly_list.append([number,lon[n],lat[n],en])

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
        return
   
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
# Function on basemap creations as lqm (low qual) and iqm (medium qual)
def basemap_creator():
    # NOTE Based on OSISAF's grid
    # low quality map
    lqm = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),\
            projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    # intermediate quality map
    iqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
            projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    # high quality map
    hqm = Basemap(width=7600000,height=11200000,resolution='h',rsphere=(6378273,6356889.44891),\
            projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
   
    #TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO
    # NOTE Based on WIM's grid
    #lqmM = Basemap(width=6090000,height=8810000,resolution='l',rsphere=(6378273,6356889.44891),\
    #        projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    ## better quality map
    #iqmM = Basemap(width=6090000,height=8810000,resolution='i',rsphere=(6378273,6356889.44891),\
    #        projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    return(lqm,iqm,hqm)

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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# RUNS
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

###########################################################################
# REGIONAL DAILY RUN 
def regional_daily(X,Y,Z,basemap=None,run=None):
    
    import MIZchar as widths

    # We need the contours and 4 state map (IP=2,MIZ=1,Ocean=0,Land=NaN)
    if run == 'OSI' or run == 'MIC':
        mode = 'ic'
        MIZ_p = MIZ_poly(X,Y,Z,mode)
    elif run == 'MFD':
        mode = 'fsd'
        MIZ_p = MIZ_poly(X,Y,Z,mode)
    else:
        print('Unknown Run')
    
    IPMIZ,B1,B2,MIZ_cont = MIZ_p.MIZ, MIZ_p.B1, MIZ_p.B2, MIZ_p.MIZcont

    # Prep the lists for polygon.txt
    global bar_poly_list
    bar_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
    global les_poly_list
    les_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
    global ncb_poly_list
    ncb_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
    global lab_poly_list
    lab_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]
    global gre_poly_list
    gre_poly_list = [['Polygon_number','Longitude','Latitude','Func_value']]

    # This list of classes may be useful
    poly_list=[]

    # Analysing the polygons one by one, dividing them into their respective regions
    for n,el in enumerate(MIZ_cont):
        polygons=poly_stat(el,B2,B1,IPMIZ,X,Y,n,basemap=basemap)
        poly_list.append(polygons)

    # BAR REGION
    reg_repo = str(out_dir)+'/bar'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    # DAILY
    dailydir = reg_repo+'/daily'
    if not os.path.exists(dailydir):
        os.mkdir(dailydir)
    # need variables for weighted location
    num_lon = 0
    num_lat = 0
    den = 0
    for poly in poly_list:
        if poly.region == 'bar':
            filname = str(dailydir)+'/'+str(dadate)+'.txt'
            with open(filname,'a') as f:
                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
                    poly.centroid_longitude,poly.centroid_latitude]]
                # Keep some info aside
                num_lon += poly.centroid_longitude*poly.E_area 
                num_lat += poly.centroid_latitude*poly.E_area 
                den += poly.E_area
                str1 = ' '.join(map(str,row1))
                f.write(str1+'\n')
                f.close
    mean_lon_bar = num_lon/float(den)
    mean_lat_bar = num_lat/float(den)
    # POLYGONS
    polydir = reg_repo+'/polygons'
    if not os.path.exists(polydir):
        os.mkdir(polydir)
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(bar_poly_list):
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
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    reg_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(reg_widths):
        if wd is not None:
            reg_width_list.append(wd.int_widths)
    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(reg_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    uglyway = open(filname,'r+')
    uglierway = uglyway.readlines()
    w_num = 0
    w_den = 0
    for el in uglierway:
        pw = np.array(map(float,el.split()))
        w_num += sum(pw)
        w_den += len(pw)
    avg_width_bar = w_num/float(w_den)

    # LES REGION
    reg_repo = str(out_dir)+'/les'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    # DAILY
    dailydir = reg_repo+'/daily'
    if not os.path.exists(dailydir):
        os.mkdir(dailydir)
    # need variables for weighted location
    num_lon = 0
    num_lat = 0
    den = 0
    for poly in poly_list:
        if poly.region == 'les':
            filname = str(dailydir)+'/'+str(dadate)+'.txt'
            with open(filname,'a') as f:
                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
                    poly.centroid_longitude,poly.centroid_latitude]]
                # Keep some info aside
                num_lon += poly.centroid_longitude*poly.E_area 
                num_lat += poly.centroid_latitude*poly.E_area 
                den += poly.E_area
                str1 = ' '.join(map(str,row1))
                f.write(str1+'\n')
                f.close
    mean_lon_les = num_lon/float(den)
    mean_lat_les = num_lat/float(den)
    # POLYGONS
    polydir = reg_repo+'/polygons'
    if not os.path.exists(polydir):
        os.mkdir(polydir)
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(les_poly_list):
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
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    reg_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(reg_widths):
        if wd is not None:
            reg_width_list.append(wd.int_widths)
    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(reg_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    uglyway = open(filname,'r+')
    uglierway = uglyway.readlines()
    w_num = 0
    w_den = 0
    for el in uglierway:
        pw = np.array(map(float,el.split()))
        w_num += sum(pw)
        w_den += len(pw)
    avg_width_les = w_num/float(w_den)

    # NCB REGION
    reg_repo = str(out_dir)+'/ncb'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    # DAILY
    dailydir = reg_repo+'/daily'
    if not os.path.exists(dailydir):
        os.mkdir(dailydir)
    # need variables for weighted location
    num_lon = 0
    num_lat = 0
    den = 0
    for poly in poly_list:
        if poly.region == 'ncb':
            filname = str(dailydir)+'/'+str(dadate)+'.txt'
            with open(filname,'a') as f:
                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
                    poly.centroid_longitude,poly.centroid_latitude]]
                # Keep some info aside
                num_lon += poly.centroid_longitude*poly.E_area 
                num_lat += poly.centroid_latitude*poly.E_area 
                den += poly.E_area
                str1 = ' '.join(map(str,row1))
                f.write(str1+'\n')
                f.close
    mean_lon_ncb = num_lon/float(den)
    mean_lat_ncb = num_lat/float(den)
    # POLYGONS
    polydir = reg_repo+'/polygons'
    if not os.path.exists(polydir):
        os.mkdir(polydir)
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(ncb_poly_list):
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
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    reg_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(reg_widths):
        if wd is not None:
            reg_width_list.append(wd.int_widths)
    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(reg_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    uglyway = open(filname,'r+')
    uglierway = uglyway.readlines()
    w_num = 0
    w_den = 0
    for el in uglierway:
        pw = np.array(map(float,el.split()))
        w_num += sum(pw)
        w_den += len(pw)
    avg_width_ncb = w_num/float(w_den)

    # LAB REGION
    reg_repo = str(out_dir)+'/lab'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    # DAILY
    dailydir = reg_repo+'/daily'
    if not os.path.exists(dailydir):
        os.mkdir(dailydir)
    # need variables for weighted location
    num_lon = 0
    num_lat = 0
    den = 0
    for poly in poly_list:
        if poly.region == 'lab':
            filname = str(dailydir)+'/'+str(dadate)+'.txt'
            with open(filname,'a') as f:
                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
                    poly.centroid_longitude,poly.centroid_latitude]]
                # Keep some info aside
                num_lon += poly.centroid_longitude*poly.E_area 
                num_lat += poly.centroid_latitude*poly.E_area 
                den += poly.E_area
                str1 = ' '.join(map(str,row1))
                f.write(str1+'\n')
                f.close
    mean_lon_lab = num_lon/float(den)
    mean_lat_lab = num_lat/float(den)
    # POLYGONS
    polydir = reg_repo+'/polygons'
    if not os.path.exists(polydir):
        os.mkdir(polydir)
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(lab_poly_list):
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
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    reg_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(reg_widths):
        if wd is not None:
            reg_width_list.append(wd.int_widths)
    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(reg_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    uglyway = open(filname,'r+')
    uglierway = uglyway.readlines()
    w_num = 0
    w_den = 0
    for el in uglierway:
        pw = np.array(map(float,el.split()))
        w_num += sum(pw)
        w_den += len(pw)
    avg_width_lab = w_num/float(w_den)

    # GRE REGION
    reg_repo = str(out_dir)+'/gre'
    if not os.path.exists(reg_repo):
        os.mkdir(reg_repo)
    # DAILY
    dailydir = reg_repo+'/daily'
    if not os.path.exists(dailydir):
        os.mkdir(dailydir)
    # need variables for weighted location
    num_lon = 0
    num_lat = 0
    den = 0
    for poly in poly_list:
        if poly.region == 'gre':
            filname = str(dailydir)+'/'+str(dadate)+'.txt'
            with open(filname,'a') as f:
                row1 = [[poly.number,poly.polygon_class,poly.E_area,poly.E_perim, \
                    poly.centroid_longitude,poly.centroid_latitude]]
                # Keep some info aside
                num_lon += poly.centroid_longitude*poly.E_area 
                num_lat += poly.centroid_latitude*poly.E_area 
                den += poly.E_area
                str1 = ' '.join(map(str,row1))
                f.write(str1+'\n')
                f.close
    mean_lon_gre = num_lon/float(den)
    mean_lat_gre = num_lat/float(den)
    # POLYGONS
    polydir = reg_repo+'/polygons'
    if not os.path.exists(polydir):
        os.mkdir(polydir)
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(gre_poly_list):
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
    filname = str(polydir)+'/'+str(dadate)+'.txt'
    reg_widths = widths.single_file(filname,basemap)
    for w,wd in enumerate(reg_widths):
        if wd is not None:
            reg_width_list.append(wd.int_widths)
    filname = str(widthsdir)+'/'+str(dadate)+'.txt'
    with open(filname,'a') as f:
        for l,el in enumerate(reg_width_list):
            string = ' '.join(map(str,el))
            for item in string:
                f.write(item)
            f.write('\n')
        f.close()
    uglyway = open(filname,'r+')
    uglierway = uglyway.readlines()
    w_num = 0
    w_den = 0
    for el in uglierway:
        pw = np.array(map(float,el.split()))
        w_num += sum(pw)
        w_den += len(pw)
    avg_width_gre = w_num/float(w_den)

    # daily regional gen_stats
    # Widhts, lons and lats lists
    widths = [avg_width_bar,avg_width_les,avg_width_ncb,avg_width_lab,avg_width_gre]
    lons = [mean_lon_bar,mean_lon_les,mean_lon_ncb,mean_lon_lab,mean_lon_gre]
    lats = [mean_lat_bar,mean_lat_les,mean_lat_ncb,mean_lon_lat,mean_lat_gre]
    
    regional_stats(Z,widths,lons,lats,mode)

    plt.imshow(IPMIZ)
    plt.colorbar()
    plt.savefig(str(out_dir)+'/'+str(dadate)+'.png',bbox_inches='tight')
    plt.close()
   
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Beginning of the script
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# MISSION and DATE
run = sys.argv[1]
dadate = sys.argv[2] 
out_dir = sys.argv[3]
jday = sys.argv[4]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Start the count and analyse the datasets
time0 = time.time()

# Call in the basemap (for now just OSI's) TODO the model's
lqm,iqm,hqm = basemap_creator()

# Call in the regions
# LABELS -> [BAR = 1, LES = 2, NCB = 3, LAB = 4, GRE = 5]
regions = np.load('./geo_info/regions.npy')

# Loading X_OSISAF and Y_OSISAF
XO = np.load('./geo_info/X_OSI.npy')
YO = np.load('./geo_info/Y_OSI.npy')

# Read in the data
if run == 'OSI':
    osisaf = reader('Osisaf',dadate,iqm,jday)
    lonO,latO,X,Y,Z = osisaf.lonO,osisaf.latO,osisaf.X,osisaf.Y,osisaf.Z
elif run == 'MIC':
    model = reader('Model',dadate,iqm,jday)
    #lonM,latM,XMC,YMC,ZMM,ZMD = model.lonM,model.latM,model.X,model.Y,model.ZC,model.ZD
    lonM,latM,XMC,YMC,ZMM = model.lonM,model.latM,model.X,model.Y,model.ZC
    # Reproject model data
    Z = reproj(XMC,YMC,ZMM,XO,YO)
    X = np.copy(XO)
    Y = np.copy(YO)
elif run == 'MFD':
    model = reader('Model',dadate,iqm)
    lonM,latM,XMF,YMF,ZMM,ZMD = model.lonM,model.latM,model.X,model.Y,model.ZC,model.ZD
    # Reproject model data
    Z = reproj(XC,YC,ZMD,XO,YO)
    X = np.copy(XO)
    Y = np.copy(YO)

regional_daily(X,Y,Z,basemap=lqm,run=run)

elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime
