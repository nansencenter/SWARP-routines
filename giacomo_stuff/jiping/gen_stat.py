############################################################################
# SCRIPT FOR ICE EDGE VALIDATION - TO BE USED WITH ITS BASH LAUNCHER
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
from scipy.stats import pearsonr as prs_idx
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

sys.path.append('../py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution_2tw as Leqs
############################################################################

############################################################################
# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')
############################################################################

############################################################################
# Functions and classes
############################################################################
def binary_mod(data1,thresh1):
  odata = np.copy(data2)
  odata[odata<thresh1] = 0
  odata[odata>=thresh1] = 1 
  thenans = np.isnan(odata)
  odata[thenans] = 0
  return(odata)
############################################################################
# read in and prepare every file for polygon detection
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
    		self.X,self.Y,self.Z = self._read_osi_(dadate,basemap) 
    		return
    	elif name == 'Model':
    		self.X,self.Y,self.ZC,self.ZD = self._read_mdl_(dadate,basemap)
    		return
    	elif name == 'Aari':
    		self.ad = self._read_aari_(dadate,basemap)
    
    def _read_osi_(self,dadate,basemap):
    	day = self.day
    	month = self.month
    	year = self.year
    	# Read in OSI_SAF file
    	#outdir = '/work/shared/nersc/msc/OSI-SAF/'+str(year)+'_nh_polstere/'
    	outdir = './data/OSI'
    	ncfil = outdir+'/ice_conc_nh_polstere-100_multi_'+dadate+'1200.nc'
    	clon = 'lon'
    	clat = 'lat'
    	cconc = 'ice_conc'
    	lon2 = Mrdg.nc_get_var(ncfil,clon) # lon[:,:] is a numpy array
    	lat2 = Mrdg.nc_get_var(ncfil,clat) # lat[:,:] is a numpy array
    	conc = Mrdg.nc_get_var(ncfil,cconc,time_index=0)
    	xc = Mrdg.nc_get_var(ncfil,'xc')
    	yc = Mrdg.nc_get_var(ncfil,'yc')
    	X2,Y2 = basemap(lon2[:,:],lat2[:,:],inverse=False)
    	XO = np.copy(X2)
    	YO = np.copy(Y2)
    	Z2 = conc[:,:].data
    	mask2 = conc[:,:].mask
    	Z2[mask2] = np.NaN
    	ZO = Z2/100
        return(XO,YO,ZO) 
    
    def _read_mdl_(self,dadate,basemap):

        # NetCDF reader 
        day = self.day
        month = self.month
        year = self.year
        # Read TP4arch_wav
        #outdir = '/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/'+dadate+'/netcdf/'
        outdir = './data/MDL'
        ncfil = outdir+'/TP4archv_wav_start'+str(dadate)+'_000000Z_dump'+str(dadate)+'_120000Z.nc'
        slon = 'longitude'
        slat = 'latitude'
        sconc = 'fice'
        sdmax = 'dmax'
        lon = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
        lat = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
        conc = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
        dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=0)
        X,Y = basemap(lon[:,:],lat[:,:],inverse=False)
        ZD = dmax[:,:].data
        mask = dmax[:,:].mask
        ZD[mask] = np.NaN
        ZC = conc[:,:].data
        mask = conc[:,:].mask
        ZC[mask] = np.NaN
        
        return(X,Y,ZC,ZD)
		
    def _read_aari_(self,dadate,basemap):
    	# Read in AARI charts for Barents Sea
    	# NOTE this reader's output is a dictionary with a list of polygons and
    	# indexes for contours identification
    	def universal_reader(fil,reg):
            icb = open(fil)
            icblist = []
            for n,en in enumerate(icb):
            	nen = en.split(';')
            	icblist.append(nen)
                icb_cont = np.array(icblist)
                # we want to cut first row (general info) and first column (number of points)
                icb_cont = icb_cont[1:,1:]
                # it is unknown if the first line will be the poly index or the contour
            	# the finder is the first line and will be used to find the in/out
            	finder = icb_cont[0]
            	def is_number(s):
            	    try:
            	        float(s)
            	        return True
            	    except ValueError:
            	        return False
            	idx = 0
            	for el in finder:
            		if is_number(el):
            			idx += 1
            		else:
            			break
             	# change 'in' as '1.' and 'out' as '0.'
             	for el in icb_cont:
            		if el[idx] == 'in':
            		 	el[idx] = 1
            		elif el[idx] == 'out':
            		 	el[idx] = 0
            	# cut out extra columns
             	icb_cont = icb_cont[:,idx:]
             	# reads as string, we want floats
             	icb_cont = icb_cont.astype(np.float)
             	# find out number of polygons
            	if idx == 1:
            		bpolyn = np.int(icb_cont[:,0].max())
            	else:
            		bpolyn = 1
             	# split the array into polygons - creation of dict
             	ad = {}
             	for n in range(bpolyn):
             		 poly = []
             		 for m,em in enumerate(icb_cont):
             				if em[0] == n+1:
             					 poly.append(em)
             		 ad[str(reg)+'_'+str(n+1)] = np.array(poly)
            	return(ad)
            
            empty = 0
            outdir = './ice_charts/AARI/' 
            ncfil = outdir+'aari_bar_'+dadate+'.txt'
            if os.path.isfile(ncfil):
             	print('Ice Chart for Barents sea in '+dadate+'\n')
             	region = 'Barents_Kara_sea'
            	adb = universal_reader(ncfil,region)
            	b_check = 1
            else:
             	print('No Ice Chart for Barents sea in '+dadate+'\n')
            	empty += 1
            	adb = 0
            
            # Read in AARI charts for Greenland Sea
            ncfil2 = outdir+'aari_gre_'+dadate+'.txt'
            if os.path.isfile(ncfil2):
             	print('Ice Chart for Greenland sea in '+dadate+'\n')
             	region = 'Greenland_sea'
            	adg = universal_reader(ncfil2,region)
            	g_check = 1
            else:
             	print('No Ice Chart for Greenland sea in '+dadate+'\n')
            	empty += 1
            	adg = 0
            
            # Read in AARI charts for both Barents and Greenland
            # sometimes the charts are united in case the some polygons belong to both regions
            if empty == 2:
             	ncfil3 = ''.join( glob.glob('*'+dadate+'.txt'))
             	if ncfil3 == '':
            		print('No Common (Barents+Greenland) Ice Chart in '+dadate+'\n')
            		print ''
            		print 'No Aari charts available'
            		aari_presence = 0
            		adbg = 0
             	else:
            		region = 'Barents_Greenland_seas'
            		adbg = universal_reader(ncfil3,region)
            else:
            	print('No Common (Barents+Greenland) Ice Chart in '+dadate+'\n')
            	adbg = 0
            
            return(adb,adg,adbg)

############################################################################
def _reprojector_(self,X1,Y1,Z1,X2,Y2,Z2):
    time0 = time.time()
    # low quality map
    self.lqm = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),\
          projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    # better quality map
    self.hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
          projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    
    # getting ready for reprojection
    X3 = X1.reshape(X1.size)
    Y3 = Y1.reshape(Y1.size)
    Z3 = Z1.reshape(Z1.size)
    C = [X3,Y3]
    C = np.array(C)
    C = C.T
    
    # Interpolation can be done with other methods ('linear','cubic'<--doesn't work for our data)
    ZN = grd(C,Z3,(X2,Y2),method='nearest')
    
    elapsedtime = time.time() - time0
    print 'Reprojection done in ',elapsedtime
    print ''
    return(ZN)
############################################################################
# TODO modify this to make a 3 value map (OpenWater,MIZ,PackIce)
class aari_stat:
    # gets single contour and sorts it
    def __init__(self,data,region,basemap=None,PLOT=None,STCH=None):
    	self.name = name
    	self.flllist = data
    	self.region = region
    	self.lonlat_list = data[:,1:]
    	self.fval = data[:,0]
    	self.typo = 'AARI'
    	if basemap is not None:
    	    self.ll2xy(data[:,1:],basemap)
    	    self.split_cont(self.fxylist)
    	else:
    	    self.split_cont(self.flllist)
    	self.centroid(basemap=basemap)
    	self.dist_edges()
    	self.laplacian_solution()
    	if STCH:             		
    	    self.stat_chart(save=True)
    	if PLOT:             		
    	    self.poly_contour_plot()
    	if self.region == 'Barents_Kara_sea':
    	    bar_m_widths.append([self.med_width])
    	    bar_area.append([self.p_area])
    	    bar_perim.append([self.p_perim])
    	elif self.region == 'Greenland_sea':
    	    gre_m_widths.append([self.med_width])
    	    gre_area.append([self.p_area])
    	    gre_perim.append([self.p_perim])
    	else:
    	    bar_gre_widths.append([self.med_width])
    	    bar_gre_area.append([self.p_area])
    	    bar_gre_perim.append([self.p_perim])
    
    def ll2xy(self,cont,basemap):
    	fvals = self.fval
    	lonlat = self.lonlat_list
    	lon = lonlat[:,0]
    	lat = lonlat[:,1]
    	x,y = basemap(lon,lat,inverse=False)
    	xy_list = zip(fvals,x,y)
    	xy_list = np.array(xy_list)
    	self.fxylist = xy_list
    	print xy_list
    	return
    
    def split_cont(self,trilist):
     	# this function find the different contours
     	# NOTE if a point has non integer i coordinate is going to be a vertical edge,
     	# if has non integer j coordinate is going to be a horizontal edge hence the use
     	# of different arounds depending on the indexes of the point
     	# NOTE if the polygon is an OVERESTIMATION the contour finding is inverted
    	cont = trilist
     	in_cont = []
     	out_cont = []
    	for el in cont:
    		if el[0] == 1:
    			in_cont.append([el[1],el[2]])
    		else:
    			out_cont.append([el[1],el[2]])
     	in_cont = np.array(in_cont)
     	out_cont = np.array(out_cont)
     	self.in_ice_edge = in_cont
     	self.out_ice_edge = out_cont
     	return
    
    def centroid(self,basemap=None):
    	xy_list = self.fxylist[:,1:]
     	# Centroid - Location in lon/lat of the centroid of the polygon
     	DX	= 0
     	DY	= 0
     	B		= 0
     	for n in range(len(xy_list)-1):
     		 DX	+= ((xy_list[n][0]+xy_list[n+1][0])*\
     		       ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
     		 DY	+= ((xy_list[n][1]+xy_list[n+1][1])*\
     		       ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
     		 B	+= ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1]))
     	A = (0.5)*B 
     	CX = (1/float(6*A))*DX
     	CY = (1/float(6*A))*DY
    	centroid = [CX,CY]
    	self.centroid = centroid
    	if basemap is not None:
    		self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
     	
    def area_perimeter(self):
    	# Calculating the area of irregular polygon (from perimeter)
    	vs = self.fxylist[:,1:]
    	a = 0
    	x0,y0 = vs[0]
    	for [x1,y1] in vs[1:]:
    		dx = x1-x0
    		dy = y1-y0
    		a += 0.5*(y0*dx - x0*dy)
    		x0 = x1 
    		y0 = y1 
    	self.area_euclidean = abs(a*100)
    	# Calculating perimeter in xy coordinates (unit = 10km)
    	perim = 0
    	for n in range(len(vs)-1):
    		perim += np.sqrt(pow(vs[n+1][0]-vs[n][0],2)+pow(vs[n+1][1]-vs[n][1],2))
    	self.perimeter_euclidean = abs(perim*10)
    	return
    
    def dist_edges(self):
    	# Calculating the min distance from a model point to any osisaf point
    	# and vice-versa (osisaf to model)
    	# same script applied for the Model2Model(in&out) products - see legends
    	inside_contour = self.in_ice_edge
    	outside_contour	= self.out_ice_edge
    	tcont = self.fxylist[:,1:]
    	UKW   = 'Unknown contours < 20%'
    	DMW   = 'Contour difference < 40%'
    	width_in2out = []
    	width_out2in = []
    	if len(inside_contour) == 0 or len(outside_contour) == 0:
    		width_in2out	= [0,0,0,0,0]
    		width_out2in	= [0,0,0,0,0]
    	else:
    		dmo = 100*(abs(len(inside_contour)-len(outside_contour))/float(len(tcont)))
    		if dmo >= 40:
    			DMW = 'WARNING - contours difference > 40%'
    		for n,en in enumerate(inside_contour):
    			dist_pt = []
    			for m,em in enumerate(outside_contour):
    				dist1	= np.sqrt(pow(en[0]-em[0],2)+pow(en[1]-em[1],2))
    				dist_pt.append([dist1,en[0],en[1],em[0],em[1]])
    			idx,value = min(enumerate(dist_pt),key=itg(1))
    			width_in2out.append(value)
    		for n,en in enumerate(outside_contour):
    			dist_pt = []
    			for m,em in enumerate(inside_contour):
    				dist1	= np.sqrt(pow(en[0]-em[0],2)+pow(en[1]-em[1],2))
    				dist_pt.append([dist1,en[0],en[1],em[0],em[1]])
    			idx,value = min(enumerate(dist_pt),key=itg(1))
    			width_out2in.append(value)
    	self.widths_in2out = np.array(width_in2out)
    	self.widths_out2in = np.array(width_out2in)
    	self.DMW = DMW
    	return
    
    def laplacian_solution(self):
    	lon = self.flllist[:,1]
    	lat = self.flllist[:,2]
    	xy = self.fxylist[:,1:]
     	fval = self.fval
     	region = self.region
    	typo = self.typo
     	daname = 'Aari_'+str(self.region)
     	f_out = './outputs/'+str(typo)+'/'+str(dadate)
     	basemap = hqm
     	results = Leqs.get_MIZ_widths(lon,lat,fval,name=daname,region=region,fig_outdir=f_out,basemap=basemap,xy_coords2=xy)
     	self.AI = results[0]
     	self.fun_sol = results[1]
     	self.stream = results[2]
    	# If the result from the Laplacian solution is NaN (is different from itself)
    	# the basic distance will be calculated and used instead. In case both methods
    	# fail a nominal 10km (min resolution) width is associated to the poly 
    	if results[3] != results[3]:
    		in2out = self.widths_in2out
    		out2in = self.widths_out2in
    		if in2out.tolist() != [0,0,0,0,0] and out2in.tolist() != [0,0,0,0,0]: 
    			in2out_m = np.mean(in2out[:,0])
    			out2in_m = np.mean(out2in[:,0])
    			self.med_width = (in2out_m+out2in_m)/float(2)
    			#self.med_width = (np.mean(self.widths_in2out[:,0])+np.mean(self.widths_out2in[:,0]))/float(2)
    			print 'M.WIDTH(E) : ',self.med_width
    			self.width_status = 'Euc.'
    		else:
    			self.med_width = 10
    			print 'M.WIDTH(MIN) : ',self.med_width
    			self.width_status = 'Min.'
    	else:
    		self.med_width = results[3]
    		print 'M.WIDTH(L) : ',self.med_width
    		self.width_status = 'Lap.'
    	# similar to above, if the area is less than 100km^2 we use the Euclidean area and perimeter as
    	# a comparison. Both results + difference are printed out 
    	if results[4] < 100:
    		self.p_area = results[4]
    		self.p_perim = results[5]
    		self.area_perimeter()
    		self.p_area_e = self.area_euclidean
    		self.p_perim_e = self.perimeter_euclidean
    		self.ap_status = 'LE'
    		print ''
    		print 'CHECK OUT THE DIFFERENCES'
    		print 'AREAS(L.,E.,DIF.): ',self.p_area,self.p_area_e,abs(self.p_area - self.p_area_e)
    		print 'PERIM.(L.,E.,DIF.): ',self.p_perim,self.p_perim_e,abs(self.p_perim - self.p_perim_e)
    	else:
    		self.p_area = results[4]
    		self.p_perim = results[5]
    		self.ap_status = 'L'
    		print 'AREA: ',self.p_area
    		print 'PERIM.: ',self.p_perim
     	return
    
    def poly_contour_plot(self):
    	# used to plot contours in the total arctic basemap for visual observations
    	inside_contour = self.in_ice_edge
    	outside_contour	= self.out_ice_edge
    	unknown_contour = self.unknown_edge
    	if len(inside_contour) != 0:
    		plt.plot(inside_contour[:,1],inside_contour[:,0],'ro',markersize=5)
    	if len(outside_contour) != 0:
    		plt.plot(outside_contour[:,1],outside_contour[:,0],'yo',markersize=5)
    	if len(unknown_contour) != 0:
    		plt.plot(unknown_contour[:,1],unknown_contour[:,0],'go',markersize=5)
    	if dist_in2out != [0,0,0,0,0]:
    		for n,en in enumerate(dist_in2out):
    			plt.plot([en[2],en[4]],[en[1],en[3]],color='black',alpha=0.1)
    	if dist_out2in != [0,0,0,0,0]:
    		for n,en in enumerate(dist_out2in):
    			plt.plot([en[2],en[4]],[en[1],en[3]],color='magenta',alpha=0.1)
    	return
    
    def stat_chart(self,save=False):
    	# The Statistical Chart is an output that tries to compress as many informations
    	# as possible in a single figure.
    	# The figure is:
    	# 1) Top left arctic map with studied polygon highlighted
    	#	2) Top right a close up of the polygon with contours and min. distances
    	# 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
    	# 4) Bottom right is a legend for the small image
    	typo = self.typo 
    	region = self.region
    	pname = 'Aari_'+str(region)
    	print ''
    	print 'Statistic Chart for ',pname
    	print ''
    	xy_list = self.fxylist[:,1:]
    	inside_contour = self.in_ice_edge
    	outside_contour	= self.out_ice_edge
    	if basemap is not None:
    		clon = '%1.2f' %self.centroid_longitude
    		clat = '%1.2f' %self.centroid_latitude
    		clonlat = '{0}/{1}'.format(clon,clat)
    	else:
    		clonlat = 'Basemap missing!'
    	area = self.p_area
    	area = '%1.4e' %area + ' ' + self.ap_status
    	perim = self.p_perim
    	perim = '%1.4e' %perim + ' ' + self.ap_status
    	dist = self.med_width
    	dist = '%1.2f' %dist_in + ' ' + self.width_status
    
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
    	main.imshow(DN,interpolation='nearest',cmap='winter')
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
    	if dist_in2out.tolist() != [0,0,0,0,0]:
    		for n,en in enumerate(dist_in2out):
    			polyf.plot([en[2],en[4]],[en[1],en[3]],color='black',alpha=0.1)
    	if dist_out2in.tolist() != [0,0,0,0,0]:
    		for n,en in enumerate(dist_out2in):
    			polyf.plot([en[2],en[4]],[en[1],en[3]],color='magenta',alpha=0.1)
    
    	if typo == 'DFP':
    		# Legend on the bottom right
    		mc = mlines.Line2D([],[],color='red',marker='o')
    		oc = mlines.Line2D([],[],color='yellow',marker='o')
    		uc = mlines.Line2D([],[],color='green',marker='o')
    		md = mlines.Line2D([],[],color='grey')
    		od = mlines.Line2D([],[],color='magenta')
    		leg.legend([mc,oc,uc,md,od],(\
    		      'Open Water Cont.','Ice Pack Cont.','Unknown Cont.','Dist. Water to Pack', \
    		      'Dist. Pack to Water'),loc='center')
    		           
    		# Statistics text on the bottom left
    		txt = '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
    		      '2) Area = '+str(area)+' km^2\n'+ \
    		      '3) Perimeter = '+str(perim)+' km\n'+ \
    		      '4) Mean W-P Width = '+str(dist_in)+' km\n'+ \
    		      '5) Mean P-W Width = '+str(dist_out)+' km\n'+ \
    		      '6) '+str(DMW)+'\n'+ \
    		      '7) '+str(UKW)
    		tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
    		if save:
    			valid_class = './outputs/'+str(typo)+'/'+str(dadate)
    			if not os.path.exists(valid_class):
    				os.mkdir(valid_class)
    			if not os.path.exists(valid_class+'/'+region):
    				os.mkdir(valid_class+'/'+region)
    			fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
    			plt.close()
    		else:
    			plt.show(False)
    		print 'Statistic chart done for '+str(pname)
    	elif typo == 'ICP':
    		# Legend on the bottom right
    		mc = mlines.Line2D([],[],color='red',marker='o')
    		oc = mlines.Line2D([],[],color='yellow',marker='o')
    		uc = mlines.Line2D([],[],color='green',marker='o')
    		md = mlines.Line2D([],[],color='grey')
    		od = mlines.Line2D([],[],color='magenta')
    		leg.legend([mc,oc,uc,md,od],(\
    		      '15 % Cont.','80 % Cont.','Unknown Cont.','Dist. 15 to 80', \
    		      'Dist. 80 to 15'),loc='center')
    		           
    		# Statistics text on the bottom left
    		txt = '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
    		      '2) Area = '+str(area)+' km^2\n'+ \
    		      '3) Perimeter = '+str(perim)+' km\n'+ \
    		      '4) Mean 15-80 Width = '+str(dist_in)+' km\n'+ \
    		      '5) Mean 80-15 Width = '+str(dist_out)+' km\n'+ \
    		      '6) '+str(DMW)+'\n'+ \
    		      '7) '+str(UKW)
    		tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
    		if save:
    			valid_class = './outputs/'+str(typo)+'/'+str(dadate)
    			if not os.path.exists(valid_class):
    				os.mkdir(valid_class)
    			if not os.path.exists(valid_class+'/'+region):
    				os.mkdir(valid_class+'/'+region)
    			fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
    			plt.close()
    		else:
    			plt.show(False)
    		print 'Statistic chart done for '+str(pname)
    	else:
    		# Legend on the bottom right
    		mc = mlines.Line2D([],[],color='red',marker='o')
    		oc = mlines.Line2D([],[],color='yellow',marker='o')
    		uc = mlines.Line2D([],[],color='green',marker='o')
    		md = mlines.Line2D([],[],color='grey')
    		od = mlines.Line2D([],[],color='magenta')
    		pos_p = mpatches.Patch(color='lightgreen')
    		neg_p = mpatches.Patch(color='royalblue')
    		leg.legend([mc,oc,uc,md,od,pos_p,neg_p],(\
    		      'Model Cont.','Observation Cont.','Unknown Cont.','Dist. Mdl to Osi', \
    		      'Dist. Osi to Mdl','Model Underestimate','Model Overestimate'),loc='center')
    	           
    		# Statistics text on the bottom left
    		txt = '1) Center Lon/Lat = '+str(clonlat)+' degrees\n'+ \
    		      '2) Area = '+str(area)+' km^2\n'+ \
    		      '3) Perimeter = '+str(perim)+' km\n'+ \
    		      '4) Mean M-O Width = '+str(dist_in)+' km\n'+ \
    		      '5) Mean O-M Width = '+str(dist_out)+' km\n'+ \
    		      '6) Status = '+str(pstat)+'\n'+ \
    		      '7) '+str(DMW)+'\n'+ \
    		      '8) '+str(UKW)
    		tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
    		if save:
    			outdir = './outputs/'+str(typo)
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
    	return
###########################################################################
class ds_stats:
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
    		self.X,self.Y,self.Z = self._read_osi_(dadate,basemap) 
    		return
    	elif name == 'Model':
    		self.X,self.Y,self.ZC,self.ZD = self._read_mdl_(dadate,basemap)
    		return
    	elif name == 'Aari':
    		self.ad = self._read_aari_(dadate,basemap)

#    def correlation(self,dset1,dset2):
#        if len(dset1) != len(dset2):
#            print('DIFFERENT DATASETS LENGTHS!')
#            quit
#        else
#            pearson = prs_idx(dset1,dset2)

###########################################################################
# Beginning of the script
###########################################################################

###########################################################################
# Defining the basemap
# low quality map
lqm = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),\
      projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
# better quality map
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
      projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
###########################################################################

mission = sys.argv[1]

# AOD RUN
if mission == 'AOD':
    time0 = time.time()
    dadate = sys.argv[2] 
    osisaf = reader('Osisaf',dadate,hqm)
    model = reader('Model',dadate,hqm)
    XO,YO,ZO = osisaf.X,osisaf.Y,osisaf.Z
    XM,YM,ZM = model.X,model.Y,model.ZC
    AOD = AOD_poly(XM,YM,ZM,XO,YO,ZO)
    over,under,DN,BM,BO = AOD.over,AOD.under,AOD.DN,AOD.B1,AOD.B2
    
    bar_poly_stat = []
    gre_poly_stat = []
    lab_poly_stat = []
    les_poly_stat = []
    ncb_poly_stat = []
    
    poly_list=[]
    for n,el in enumerate(over):
    	# classification of positive polygons
    	aod=poly_stat(el,'1',BO,BM,DN,XO,YO,n,'AOD',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(aod)
    try:
    	n
    except NameError:
    	n = -1 
    for n2,el in enumerate(under):
    	# classification of negative (n+1 for good enumeration)
    	n += 1 
    	aod=poly_stat(el,'0',BO,BM,DN,XO,YO,n,'AOD',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(aod)
    
    outdir = './outputs/AOD/'+str(dadate)
    if not os.path.exists(outdir):
    	os.mkdir(outdir)
    
    reg_repo = outdir+'/Barents_Kara_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/AOD/'+str(dadate)+'/Barents_Kara_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), bar_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Greenland_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/AOD/'+str(dadate)+'/Greenland_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), gre_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Labrador_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/AOD/'+str(dadate)+'/Labrador_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), lab_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Laptev_East_Siberian_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/AOD/'+str(dadate)+'/Laptev_East_Siberian_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), les_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/North_Canada_Beaufort_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/AOD/'+str(dadate)+'/North_Canada_Beaufort_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), ncb_poly_stat)))
    f.close()
    
    elapsedtime = time.time() - time0
    print str(dadate)+' done in ',elapsedtime

# ICP model test run
if mission == 'ICP':
    time0 = time.time()
    dadate = sys.argv[2] 
    model = reader('Model',dadate,hqm)
    XM,YM,ZM = model.X,model.Y,model.ZC
    SDA = SDA_poly(XM,YM,ZM)
    over,under,DN,BM1,BM2 = SDA.over,SDA.under,SDA.DN,SDA.B1,SDA.B2
    
    bar_poly_stat = []
    gre_poly_stat = []
    lab_poly_stat = []
    les_poly_stat = []
    ncb_poly_stat = []
    
    poly_list=[]
    for n,el in enumerate(over):
    	# classification of positive polygons
    	sda=poly_stat(el,'1',BM2,BM1,DN,XM,YM,n,'ICP',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(sda)
    try:
    	n
    except NameError:
    	n = -1 
    for n2,el in enumerate(under):
    	# classification of negative (n+1 for good enumeration)
    	n += 1 
    	sda=poly_stat(el,'0',BM2,BM1,DN,XM,YM,n,'ICP',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(sda)
    
    outdir = './outputs/ICP/'+str(dadate)
    if not os.path.exists(outdir):
    	os.mkdir(outdir)
    
    reg_repo = outdir+'/Barents_Kara_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/ICP/'+str(dadate)+'/Barents_Kara_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), bar_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Greenland_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/ICP/'+str(dadate)+'/Greenland_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), gre_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Labrador_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/ICP/'+str(dadate)+'/Labrador_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), lab_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Laptev_East_Siberian_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/ICP/'+str(dadate)+'/Laptev_East_Siberian_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), les_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/North_Canada_Beaufort_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/ICP/'+str(dadate)+'/North_Canada_Beaufort_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), ncb_poly_stat)))
    f.close()
    
    elapsedtime = time.time() - time0
    print str(dadate)+' done in ',elapsedtime

# DFP model
if mission == 'DFP':
    time0 = time.time()
    dadate = sys.argv[2] 
    model = reader('Model',dadate,hqm)
    XM,YM,ZM = model.X,model.Y,model.ZD
    SDA = SDA_poly(XM,YM,ZM)
    over,under,DN,BM1,BM2 = SDA.over,SDA.under,SDA.DN,SDA.B1,SDA.B2
    	
    bar_poly_stat = []
    gre_poly_stat = []
    lab_poly_stat = []
    les_poly_stat = []
    ncb_poly_stat = []
    
    poly_list=[]
    for n,el in enumerate(over):
    	# classification of positive polygons
    	sda=poly_stat(el,'1',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(sda)
    try:
    	n
    except NameError:
    	n = -1 
    for n2,el in enumerate(under):
    	# classification of negative (n+1 for good enumeration)
    	n += 1 
    	sda=poly_stat(el,'0',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
    	poly_list.append(sda)
    
    outdir = './outputs/DFP/'+str(dadate)
    if not os.path.exists(outdir):
    	os.mkdir(outdir)
    
    reg_repo = outdir+'/Barents_Kara_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/DFP/'+str(dadate)+'/Barents_Kara_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), bar_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Greenland_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/DFP/'+str(dadate)+'/Greenland_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), gre_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Labrador_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/DFP/'+str(dadate)+'/Labrador_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), lab_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/Laptev_East_Siberian_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/DFP/'+str(dadate)+'/Laptev_East_Siberian_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), les_poly_stat)))
    f.close()
    
    reg_repo = outdir+'/North_Canada_Beaufort_sea'
    if not os.path.exists(reg_repo):
    	os.mkdir(reg_repo)
    
    f = open('./outputs/DFP/'+str(dadate)+'/North_Canada_Beaufort_sea/analysis.txt', 'w')
    f.write('\n'.join(map(lambda x: str(x), ncb_poly_stat)))
    f.close()
    
    elapsedtime = time.time() - time0
    print str(dadate)+' done in ',elapsedtime

