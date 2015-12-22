############################################################################
# THIS SCRIPT WILL COMPARE MODEL IC with OSISAF IC
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
# FUNCTIONS
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

    ## find the regions
    #for i in range(len(DN[:,0])):
    #    for j in range(len(DN[0,:])):
    #        lon = lonO[i,j]
    #        lat = latO[i,j]
    #        if lon >= 16 and lon < 90:
    #            region[i,j] = 1
    #        elif lon >= 90 and lon < 180:
    #            region[i,j] = 2
    #        elif lon >= -180 and lon < -90:
    #            region[i,j] = 3
    #        elif lon >= -90 and lon < -44:
    #            if lat > bblatw[0]:
    #                region[i,j] = 3
    #            elif lat < bblatw[0] and lat > bblatw[-1]:
    #                n = 0
    #                N = len(bblatw[:-1])
    #                while n < N:
    #                    if lat >= bblatw[n+1] and lat < bblatw[n]:
    #                        if lon >= bblonw[n]:
    #                            region[i,j] = 4
    #                            break
    #                        else:
    #                            region[i,j] = 3
    #                            break
    #                    else:
    #                        n += 1
    #            else:
    #                if lon >= bblonw[-1]:
    #                    region[i,j] = 4
    #                else:
    #                    region[i,j] = 3
    #        elif (lon >= -44 and lon < 0) or (lon >= 0 and lon <= 16):
    #            region[i,j] = 5
    #        else:
    #            region[i,j] = 6

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

    return(gen_arctic_stats,over_width_list,under_width_list)

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
dadate = sys.argv[1] 
out_dir = sys.argv[2]

# Start the count and analyse the datasets
time0 = time.time()

# Read in the data
osisaf = reader('Osisaf',dadate,hqm)
model = reader('Model',dadate,hqm)
lonO,latO,XO,YO,ZO = osisaf.lonO,osisaf.latO,osisaf.X,osisaf.Y,osisaf.Z
lonM,latM,XM,YM,ZMM = model.lonM,model.latM,model.X,model.Y,model.ZC

# Reproject model data
ZM = reproj(XM,YM,ZMM,XO,YO)

# Difference
DN = ZM-ZO

# Regions [BAR = 1, LES = 2, NCB = 3, LAB = 4, GRE = 5]
regions = np.load('regions.npy')

# Arctic Stats
bias = np.nansum(DN)/len(DN)
rmsd = np.sqrt(np.nansum(np.power(DN,2)))/len(DN)
mdl_extent = len(ZM[ZM>.15])*100
obs_extent = len(ZO[ZO>.15])*100
mdl_iarea = np.copy(ZM)*100
osi_iarea = np.copy(ZO)*100

# Regional Stats
bar_mvals = []
les_mvals = []
ncb_mvals = []
lab_mvals = []
gre_mvals = []

bar_ovals = []
les_ovals = []
ncb_ovals = []
lab_ovals = []
gre_ovals = []

bar_diff = []
les_diff = []
ncb_diff = []
lab_diff = []
gre_diff = []

for i in range(len(DN[:,0])):
    for j in range(len(DN[0,:])):
        if regions[i,j] == 1:
            bar_mvals.append(ZM[i,j])
            bar_ovals.append(ZO[i,j])
            bar_diff.append(DN[i,j])
        elif regions[i,j] == 2:
            les_mvals.append(ZM[i,j])
            les_ovals.append(ZO[i,j])
            les_diff.append(DN[i,j])
        elif regions[i,j] == 3:
            ncb_mvals.append(ZM[i,j])
            ncb_ovals.append(ZO[i,j])
            ncb_diff.append(DN[i,j])
        elif regions[i,j] == 4:
            lab_mvals.append(ZM[i,j])
            lab_ovals.append(ZO[i,j])
            lab_diff.append(DN[i,j])
        elif regions[i,j] == 5:
            gre_mvals.append(ZM[i,j])
            gre_ovals.append(ZO[i,j])
            gre_diff.append(DN[i,j])

bar_bias = np.nansum(bar_diff)/len(bar_diff)
les_bias = np.nansum(les_diff)/len(les_diff)
ncb_bias = np.nansum(ncb_diff)/len(ncb_diff)
lab_bias = np.nansum(lab_diff)/len(lab_diff)
gre_bias = np.nansum(gre_diff)/len(gre_diff)
bar_rmsd = np.sqrt(np.nansum(np.power(bar_diff,2)))/len(bar_diff)
les_rmsd = np.sqrt(np.nansum(np.power(les_diff,2)))/len(les_diff)
ncb_rmsd = np.sqrt(np.nansum(np.power(ncb_diff,2)))/len(ncb_diff)
lab_rmsd = np.sqrt(np.nansum(np.power(lab_diff,2)))/len(lab_diff)
gre_rmsd = np.sqrt(np.nansum(np.power(gre_diff,2)))/len(gre_diff)

bar_mextent = len(bar_mvals[bar_mvals>.15])*100
les_mextent = len(les_mvals[les_mvals>.15])*100
ncb_mextent = len(ncb_mvals[ncb_mvals>.15])*100
lab_mextent = len(lab_mvals[lab_mvals>.15])*100
gre_mextent = len(gre_mvals[gre_mvals>.15])*100
bar_marea = np.copy(bar_mvals)*100
les_marea = np.copy(les_mvals)*100
ncb_marea = np.copy(ncb_mvals)*100
lab_marea = np.copy(lab_mvals)*100
gre_marea = np.copy(gre_mvals)*100

bar_oextent = len(bar_ovals[bar_ovals>.15])*100
les_oextent = len(les_ovals[les_ovals>.15])*100
ncb_oextent = len(ncb_ovals[ncb_ovals>.15])*100
lab_oextent = len(lab_ovals[lab_ovals>.15])*100
gre_oextent = len(gre_ovals[gre_ovals>.15])*100
bar_oarea = np.copy(bar_ovals)*100
les_oarea = np.copy(les_ovals)*100
ncb_oarea = np.copy(ncb_ovals)*100
lab_oarea = np.copy(lab_ovals)*100
gre_oarea = np.copy(gre_ovals)*100

outdir = str(out_dir)
if not os.path.exists(outdir):
    os.mkdir(outdir)

# BIAS FILE
filname = str(outdir)+'/bias.txt'
with open(filname,'a') as f:
    row1 = [dadate,bias,bar_bias,les_bias,ncb_bias,lab_bias,gre_bias,'artbarlesncblabgre']
    str1 = ' '.join(map(str,row1))
    f.write(str1)
    f.write('\n')
    f.close()

# RMSD FILE
filname = str(outdir)+'/rmsd.txt'
with open(filname,'a') as f:
    row1 = [dadate,rmsd,bar_rmsd,les_rmsd,ncb_rmsd,lab_rmsd,gre_rmsd,'artbarlesncblabgre']
    str1 = ' '.join(map(str,row1))
    f.write(str1)
    f.write('\n')
    f.close()

# SEA ICE EXTENT and AREA FILE
filname = str(outdir)+'/sea_ice_extent_area.txt'
with open(filname,'a') as f:
    row1 = [dadate,mdl_extent,osi_extent,mdl_iarea,osi_iarea,'mdlext,osiext,mdlarea,osiarea']
    str1 = ' '.join(map(str,row1))
    f.write(str1)
    f.write('\n')
    f.close()

#       
#if run == 'arctic':
#    gen,o_list,u_list = arctic(XO,YO,ZO,ZM,DN,over,under,lqm)
#    otot,oavg = widths_stats(o_list)
#    utot,uavg = widths_stats(u_list)
#
#    outdir = str(out_dir)+'/'+str(dadate)
#    if not os.path.exists(outdir):
#        os.mkdir(outdir)
#    reg_repo = outdir+'/arctic'
#    if not os.path.exists(reg_repo):
#        os.mkdir(reg_repo)
#
#    filname = str(reg_repo)+'/daily_stats.txt'
#    with open(filname,'a') as f:
#        row1 = ['Tot.Points','IceHits','Over','Under','SeaHits','Unknown','Avg.OverWidths','Avg.UnderWidths']
#        str1 = ' '.join(map(str,row1))
#        row2 = [gen[0],gen[1],gen[2],gen[3],gen[4],gen[5],oavg,uavg]
#        str2 = ' '.join(map(str,row2))
#        f.write(str1+'\n'+str2)
#        f.close()

elapsedtime = time.time() - time0
print str(dadate)+' done in ',elapsedtime
