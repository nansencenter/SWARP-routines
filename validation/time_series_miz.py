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

sys.path.append('../py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution_2 as Leqs
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
def baffin_bay():
  fil = ('./data/baffin_bay.txt')
  bbay = open(fil)
  bblon = []
  bblat = []
  bblone = []
  bblate = []
  bblonw = []
  bblatw = []
  for n,en in enumerate(bbay):
    nen = en.split(';')
    lon = float(nen[0])
    lat = float(nen[1])
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
def north_canada():
  fil = ('./data/north_canada.txt')
  ncan = open(fil)
  nclon = []
  nclat = []
  for n,en in enumerate(ncan):
    nen = en.split(';')
    lon = float(nen[0])
    lat = float(nen[1])
    nclon.append(lon)
    nclat.append(lat)
  return(nclon,nclat)
############################################################################
def open_close(self):
   kernel = np.ones((3,3),np.uint8)
   DN_op = np.copy(DN)
   DN_cl = np.copy(DN)
   DN_op = cv2.morphologyEx(DN_op,cv2.MORPH_OPEN,kernel)
   DN_cl = cv2.morphologyEx(DN_cl,cv2.MORPH_CLOSE,kernel)
   return(DN_op,DN_cl)
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
		outdir = './data/'
		ncfil = outdir+'ice_conc_nh_polstere-100_multi_'+dadate+'1200.nc'
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

#	#TODO I couldn't make the binary version work properly
#	# work on the binary reader, for now *.nc will be used 
#		# Binary reader version
#		def get_grid(grid_dir='.'):
#			gfil  = grid_dir+'/regional.grid.a'
#			dfil  = grid_dir+'/regional.depth.a'
#			#   
#			plon  = Mrdg.get_array_from_HYCOM_binary(gfil,1)
#			nx,ny = plon.shape
#			#   
#			plat     = Mrdg.get_array_from_HYCOM_binary(gfil,2,dims=(nx,ny))
#			depths   = Mrdg.get_array_from_HYCOM_binary(dfil,1,dims=(nx,ny))
#			return plon,plat,depths,nx,ny
#
#		icp		= 'fice' #daily average
#		icpma = 1.
#		icpmi = 0.
#		dfp   = 'dfloe'
#		dfpma = 300.
#		dfpmi = 0.
#		Bc    = 0.  # look for contours of this binary array from getbin
#		Ccol  = 'k'
#
#		##################################################################
#		# read in arrays from TP4 binaries (.a)
#		# lon/lat from grid binary
#		tdir  = './data/MDL'
#		gridir = './data'
#		plon,plat,depths,nx,ny  = get_grid(gridir)
#		# 
#		afil  = tdir+'/TP4archv_wav.'+str(self.year)+'_'+str(self.julian)+'_120000.a'
#		bfil  = afil[:-2]+'.b'
#		
#		# make dictionary from bfil to get record number in afil
#		vlst  = Mrdg.get_record_numbers_HYCOM(bfil)
#		
#		# get recno from dictionary created from bfil
#		icprecno = vlst[icp]
#		dfprecno = vlst[dfp]
#		ZC     = Mrdg.get_array_from_HYCOM_binary(afil,icprecno,dims=(nx,ny))
#		ZD     = Mrdg.get_array_from_HYCOM_binary(afil,dfprecno,dims=(nx,ny))
#		##################################################################
#
#	 	X,Y = basemap(plon[:,:],plat[:,:],inverse=False)

		# NetCDF reader 
		day = self.day
		month = self.month
		year = self.year
		# Read TP4arch_wav
		#outdir = '/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/'+dadate+'/netcdf/'
		outdir = './data/'
	 	ncfil = outdir+'TP4archv_wav_start'+str(dadate)+'_000000Z_dump'+str(dadate)+'_120000Z.nc'
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
				else:
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
# This class will find AOD polygons for 2 sets with different grid.
# NOTE the basemap used in the reprojector is thought for OSISAF database but 
# may fit any stereographical map, the reprojection should be fine as well.
class SDA_poly:
	def __init__(self,X1,Y1,Z1):
		ZM = Z1
		if len(ZM[ZM>1]) > 0:
			self.DN,self.B1,self.B2,oDN,uDN = self.binary_diff(ZM,.01,299.99)
		else:
			self.DN,self.B1,self.B2,oDN,uDN = self.binary_diff(ZM,.15,.8)
		self.over,self.under = self.poly_maker(self.DN)
		self.DN = self.mapper(self.DN,ZM)
		return

	def binary_diff(self,data1,thresh1,thresh2):

		def closing(DN):
			# apply closing to avoid small polynias and clean up a little
			kernel = np.ones((3,3),np.uint8)
			DN_cl = morph.closing(DN,kernel)
			return(DN_cl)

		# NOTE outputs are difference, overpredicitions and underpredictions (CLOSED & MAPPED)

	 	# generating the binary file, get the difference, divide into overprediction and
		# underprediction, apply the closing on both and finally re-map the land
	 	ndata1 = np.copy(data1)
		ndata2 = np.copy(data1)
	 	ndata1[ndata1<thresh1] = 0
	 	ndata2[ndata2<thresh2] = 0
	 	ndata1[ndata1>=thresh1] = 1	
	 	ndata2[ndata2>=thresh2] = 1	
		ddata = ndata2 - ndata1
		thenans = np.isnan(ddata)
		ddata[thenans] = 0
		over_data = np.copy(ddata)
		under_data = np.copy(ddata)
		over_data[over_data==-1] = 0
		under_data[under_data==1] = 0
		under_data = abs(under_data)
		o_data = closing(over_data)
		u_data = closing(under_data)
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
 	
	def mapper(self,DIFF,BIN):
	 	# getting back the terrain lost during binary_mod (for contour finding reasons)
	 	DN = DIFF
	 	ZO = BIN
	 	thenan = np.isnan(ZO)
	 	DN[thenan] = None
	 	return(DN)

############################################################################
# This class will find AOD polygons for 2 sets with different grid.
# NOTE the basemap used in the reprojector is thought for OSISAF database but 
# may fit any stereographical map, the reprojection should be fine as well.
class AOD_poly:
  def __init__(self,X1,Y1,Z1,X2,Y2,Z2):
    ZM = self._reprojector_(X1,Y1,Z1,X2,Y2,Z2)
    ZO = Z2
    self.DN,self.B1,self.B2,oDN,uDN = self.binary_diff(ZM,ZO,.15)
    self.over,self.under = self.poly_maker(self.DN)
    self.DN = self.mapper(self.DN,ZM,ZO)
    return

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

  def binary_diff(self,data1,data2,thresh):

    def closing(DN):
      # apply closing to avoid small polynias and clean up a little
      kernel = np.ones((3,3),np.uint8)
      DN_cl = morph.closing(DN,kernel)
      return(DN_cl)


    # NOTE outputs are difference, overpredicitions and underpredictions (CLOSED & MAPPED)

    # generating the binary file, get the difference, divide into overprediction and
    # underprediction, apply the closing on both and finally re-map the land
    ndata1 = np.copy(data1)
    ndata2 = np.copy(data2)
    ndata1[ndata1<thresh] = 0
    ndata2[ndata2<thresh] = 0
    ndata1[ndata1>=thresh] = 1  
    ndata2[ndata2>=thresh] = 1  
    ddata = ndata2 - ndata1
    thenans = np.isnan(ddata)
    ddata[thenans] = 0
    over_data = np.copy(ddata)
    under_data = np.copy(ddata)
    over_data[over_data==-1] = 0
    under_data[under_data==1] = 0
    under_data = abs(under_data)
    o_data = closing(over_data)
    u_data = closing(under_data)
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

############################################################################

class poly_stat:
	# gets single contour and sorts it
	def __init__(self,cont,poly_status,OUT,INS,DIFF,X,Y,number,typo,basemap=None,PLOT=None,STCH=None):
		# class definition
		if len(cont) > 200:
			self.polygon_class = 'H'
		elif len(cont) > 100 and len(cont) <= 200:
			self.polygon_class = 'B'
		elif len(cont) > 30 and len(cont) <= 100:
			self.polygon_class = 'M'
		elif len(cont) <= 30:
			self.polygon_class = 'S'
		self.polygon_status = poly_status
		self.typo = typo
		self.number = number
		self.ij_list = cont
		self.ij2xy(cont,X,Y)
		self.difference = DIFF
		self.split_cont(OUT,INS,poly_status)
		self.lon_lat(basemap)
		self.dist_edges()
		self.area_perimeter()
		self.laplacian_solution()
		if STCH:             		
			self.stat_chart(save=True)
		if PLOT:             		
			self.poly_contour_plot()
		if self.region == 'Barents_Kara_sea':
			bar_poly_stat.append([str(self.name), str(self.polygon_status), str(self.polygon_class), str(self.centroid_longitude), \
												str(self.centroid_latitude), str(self.L_width), str(self.E_width), str(self.L_area), str(self.E_area), \
												str(self.L_perim), str(self.E_perim)])
		elif self.region == 'Greenland_sea':
			gre_poly_stat.append([str(self.name), str(self.polygon_status), str(self.polygon_class), str(self.centroid_longitude), \
												str(self.centroid_latitude), str(self.L_width), str(self.E_width), str(self.L_area), str(self.E_area), \
												str(self.L_perim), str(self.E_perim)])
		elif self.region == 'Labrador_sea':
			lab_poly_stat.append([str(self.name), str(self.polygon_status), str(self.polygon_class), str(self.centroid_longitude), \
												str(self.centroid_latitude), str(self.L_width), str(self.E_width), str(self.L_area), str(self.E_area), \
												str(self.L_perim), str(self.E_perim)])
		elif self.region == 'Laptev_East_Siberian_sea':
			les_poly_stat.append([str(self.name), str(self.polygon_status), str(self.polygon_class), str(self.centroid_longitude), \
												str(self.centroid_latitude), str(self.L_width), str(self.E_width), str(self.L_area), str(self.E_area), \
												str(self.L_perim), str(self.E_perim)])
		elif self.region == 'North_Canada_Beaufort_sea':
			ncb_poly_stat.append([str(self.name), str(self.polygon_status), str(self.polygon_class), str(self.centroid_longitude), \
												str(self.centroid_latitude), str(self.L_width), str(self.E_width), str(self.L_area), str(self.E_area), \
												str(self.L_perim), str(self.E_perim)])

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
		return
	
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
	 	func_mod = 0	# value of func_vals at model ice edge
	 	func_osi = 1	# value of func_vals at OSISAF ice edge
	 	func_unk = 2	# value of func_vals at other parts of contour
	 
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
	 	func_vals = np.array(func_vals)
	 	self.in_ice_edge = in_cont
	 	self.out_ice_edge = out_cont
	 	self.unknown_edge = unk_cont
	 	self.f_vals = func_vals
	 	return
	
	def lon_lat(self,basemap):
	 	xy_list = self.xy_list
	 	
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
		self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
		self.lon_list,self.lat_list=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
		
		# defining the region
		centroid_longitude = self.centroid_longitude
		centroid_latitude = self.centroid_latitude
		bblone,bblate,bblonw,bblatw = baffin_bay()
		if centroid_longitude < 90 and centroid_longitude > 8:
			self.region = 'Barents_Kara_sea'
		elif centroid_longitude < 8 and centroid_longitude > -44:
			self.region = 'Greenland_sea'
		elif centroid_longitude < -44 and centroid_longitude > -90:
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
		
		# Calculating lon/lat for model,osisaf,unknown contours and distances
		inside_contour = self.in_ice_edge
		outside_contour	= self.out_ice_edge
		unknown_contour	= self.unknown_edge
		if len(inside_contour) != 0:
			 self.inside_contour_lon,self.inside_contour_lat = basemap(inside_contour[:,0],inside_contour[:,1],inverse=True)
		if len(outside_contour) != 0:
			 self.outside_contour_lon,self.outside_contour_lat = basemap(outside_contour[:,0],outside_contour[:,1],inverse=True)
		if len(unknown_contour) != 0:
			 self.unknown_contour_lon,self.unknown_contour_lat = basemap(unknown_contour[:,0],unknown_contour[:,1],inverse=True)
	 	return
	 	
	def laplacian_solution(self):
	 	lon = self.lon_list
	 	lat = self.lat_list
	 	xy = self.xy_list
	 	fval = self.f_vals
	 	region = self.region
		typo = self.typo
	 	daname = 'Polygon_'+str(self.number)
	 	f_out = './outputs/'+str(typo)+'/'+str(dadate)
	 	basemap = hqm
	 	results = Leqs.get_MIZ_widths(lon,lat,fval,name=daname,region=region,fig_outdir=f_out,basemap=basemap,xy_coords2=xy)
	 	self.AI = results[0]
	 	self.fun_sol = results[1]
	 	self.stream = results[2]
		self.L_width = results[3]
		self.L_area = results[4]
		self.L_perim = results[5]
		if results[3] != results[3]:
			in2out = self.widths_in2out
			out2in = self.widths_out2in
			if in2out.tolist() != [0,0,0,0,0] and out2in.tolist() != [0,0,0,0,0]: 
				self.med_width = self.E_width
				print 'M.WIDTH(E) : ',self.med_width
			else:
				self.med_width = 10
				print 'M.WIDTH(MIN) : ',self.med_width
		else:
			self.med_width = self.L_width
			print 'M.WIDTH(L) : ',self.med_width
		if results[4] < 100:
			print ''
			print 'CHECK OUT THE DIFFERENCES'
			print 'AREAS(L.,E.,DIF.): ',self.L_area,self.E_area,abs(self.L_area - self.E_area)
			print 'PERIM.(L.,E.,DIF.): ',self.L_perim,self.E_perim,abs(self.L_perim - self.E_perim)
		else:
			print 'AREA: ',self.L_area
			print 'PERIM.: ',self.L_perim
	 	return

	def poly_contour_plot(self):
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

	def dist_edges(self):
		# Calculating the min distance from a model point to any osisaf point
		# and vice-versa (osisaf to model)
		# same script applied for the Model2Model(in&out) products - see legends
		inside_contour = self.in_ice_edge
		outside_contour	= self.out_ice_edge
		unknown_contour = self.unknown_edge
		tcont = self.ij_list
		UKW   = 'Unknown contours < 20%'
		DMW   = 'Contour difference < 40%'
		width_in2out = []
		width_out2in = []
		if len(inside_contour) == 0 or len(outside_contour) == 0:
			width_in2out	= [0,0,0,0,0]
			width_out2in	= [0,0,0,0,0]
			self.E_width = 0
		else:
			unk = 100*(len(unknown_contour)/float(len(tcont)))
			dmo = 100*(abs(len(inside_contour)-len(outside_contour))/float(len(tcont)))
			if unk >= 20:
				UKW	= 'WARNING - unknown contours > 20%'
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
			gigio = np.array(width_in2out)
			topo = np.array(width_out2in)
			in2out_m = np.mean(gigio[:,0])
			out2in_m = np.mean(topo[:,0])
			self.E_width = (in2out_m+out2in_m)/float(2)
		self.widths_in2out = np.array(width_in2out)
		self.widths_out2in = np.array(width_out2in)
		self.UKW = UKW
 		self.DMW = DMW
 		return

	def stat_chart(self,save=False):
		# The Statistical Chart is an output that tries to compress as many informations
		# as possible in a single figure.
		# The figure is:
		# 1) Top left arctic map with studied polygon highlighted
		#	2) Top right a close up of the polygon with contours and min. distances
		# 3) Bottom left a recap of statistics about the polygon (only Euclidean and min dist for now)
		# 4) Bottom right is a legend for the small image

		nmbr = self.number
		pname = 'Polygon_'+str(nmbr)
		self.name = pname
		typo = self.typo
		pclass = self.polygon_class
		if self.polygon_status == 0:
			pstat = 'Overestimation'
		else:
			pstat = 'Underestimation'
		region = self.region
		print ''
		print 'Statistic Chart for ',pname
		print ''
		print 'Class and Region ',pclass,region
		print ''
		print 'Status ',pstat
		print ''
		DN = self.difference
		DMW = self.DMW
		UKW = self.UKW
		ij = self.ij_list
		inside_contour = self.in_ice_edge
		outside_contour	= self.out_ice_edge
		unknown_contour = self.unknown_edge
		dist_in2out = self.widths_in2out
		dist_out2in = self.widths_out2in
		clon = '%1.2f' %self.centroid_longitude
		clat = '%1.2f' %self.centroid_latitude
		clonlat = '{0}/{1}'.format(clon,clat)
		L_area = self.L_area
		E_area = self.E_area
		area = '%1.4e %1.4e' %(L_area,E_area) 
		L_perim = self.L_perim
		E_perim = self.E_perim
		perim = '%1.4e %1.4e' %(L_perim,E_perim)
		if dist_in2out.tolist() != [0,0,0,0,0]:
			# changing units from decakilometer to kilometer
			dist_in = np.median(dist_in2out[:,0])*10
			dist_in = '%1.2f' %dist_in
		else:
			dist_in = 'NaN'
		if dist_out2in.tolist() != [0,0,0,0,0]:
			# changing units from decakilometer to kilometer
			dist_out = np.median(dist_out2in[:,0])*10
			dist_out = '%1.2f' %dist_out
		else:
			dist_out = 'NaN'
		pstat = self.polygon_status

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
class reg_stats:
  def _init_(self,bar,gre,lab,les,ncb,typo):
    bar = self.bar
    gre = self.gre
    lab = self.lab
    les = self.les
    ncb = self.ncb
    typo = self.typo
    bar_plot = self._plotter_(bar,'Barents_Kara',typo)
    gre_plot = self._plotter_(gre,'Greenland',typo)
    lab_plot = self._plotter_(lab,'Labrador',typo)
    les_plot = self._plotter_(les,'Laptev_EastSea',typo)
    ncb_plot = self._plotter_(ncb,'NorthCanada_Beaufort',typo)

  def _plotter_(self,data,name,typo):
    m_width = data[:,0]
    area = data[:,1]
    perim = data[:,2]
    dadate = data[:,3]
    x = [date2num(dadate) for (m_width,area,perim,dadate) in data]
    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].plot(x,m_width,'r-o')
    axarr[0].set_title('Mean Width')
    axarr[1].plot(x,area,'r-o')
    axarr[1].set_title('Area')
    axarr[2].plot(x,perim,'r-o')
    axarr[2].set_title('Perim')
    fig.savefig(str(typo)+str(dadate[0])+'_'+str(dadate[-1])+'_'+str(name)+'.png')
    plt.close(fig) 
    return
###########################################################################

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
	dadate = '20150512'
	model = reader('Model','20150512',hqm)
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
		aod=poly_stat(el,'1',BM2,BM1,DN,XM,YM,n,'ICP',basemap=hqm,PLOT=None,STCH=True)
		poly_list.append(aod)
	try:
		n
	except NameError:
		n = -1 
	for n2,el in enumerate(under):
		# classification of negative (n+1 for good enumeration)
		n += 1 
		aod=poly_stat(el,'0',BM2,BM1,DN,XM,YM,n,'ICP',basemap=hqm,PLOT=None,STCH=True)
		poly_list.append(aod)
	
	bar_m_widths = sum(np.array(bar_m_widths))
	bar_area = sum(np.array(bar_area))
	bar_perim = sum(np.array(bar_perim))
	gre_m_widths = sum(np.array(gre_m_widths))
	gre_area = sum(np.array(gre_area))
	gre_perim = sum(np.array(gre_perim))
	lab_m_widths = sum(np.array(lab_m_widths))
	lab_area = sum(np.array(lab_area))
	lab_perim = sum(np.array(lab_perim))
	les_m_widths = sum(np.array(les_m_widths))
	les_area = sum(np.array(les_area))
	les_perim = sum(np.array(les_perim))
	ncb_m_widths = sum(np.array(ncb_m_widths))
	ncb_area = sum(np.array(ncb_area))
	ncb_perim = sum(np.array(ncb_perim))			
	bar_stats = [bar_m_widths,bar_area,bar_perim,dadate]
	gre_stats = [gre_m_widths,gre_area,gre_perim,dadate]
	lab_stats = [lab_m_widths,lab_area,lab_perim,dadate]
	les_stats = [les_m_widths,les_area,les_perim,dadate]
	ncb_stats = [ncb_m_widths,ncb_area,ncb_perim,dadate]
	
	elapsedtime = time.time() - time0
	print str(dadate)+' done in ',elapsedtime

# DFP model
if mission == 'DFP':
	time0 = time.time()
	dadate = '20150512'
	model = reader('Model','20150512',hqm)
	XM,YM,ZM = model.X,model.Y,model.ZD
	SDA = SDA_poly(XM,YM,ZM)
	over,under,DN,BM1,BM2 = SDA.over,SDA.under,SDA.DN,SDA.B1,SDA.B2
	
	bar_m_widths = []
	bar_area = []
	bar_perim = []
	gre_m_widths = []
	gre_area = []
	gre_perim = []
	lab_m_widths = []
	lab_area = []
	lab_perim = []
	les_m_widths = []
	les_area = []
	les_perim = []
	ncb_m_widths = []
	ncb_area = []
	ncb_perim = []
	
	poly_list=[]
	for n,el in enumerate(over):
		# classification of positive polygons
		aod=poly_stat(el,'1',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
		poly_list.append(aod)
	try:
		n
	except NameError:
		n = -1 
	for n2,el in enumerate(under):
		# classification of negative (n+1 for good enumeration)
		n += 1 
		aod=poly_stat(el,'0',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
		poly_list.append(aod)
	
	bar_m_widths = sum(np.array(bar_m_widths))
	bar_area = sum(np.array(bar_area))
	bar_perim = sum(np.array(bar_perim))
	gre_m_widths = sum(np.array(gre_m_widths))
	gre_area = sum(np.array(gre_area))
	gre_perim = sum(np.array(gre_perim))
	lab_m_widths = sum(np.array(lab_m_widths))
	lab_area = sum(np.array(lab_area))
	lab_perim = sum(np.array(lab_perim))
	les_m_widths = sum(np.array(les_m_widths))
	les_area = sum(np.array(les_area))
	les_perim = sum(np.array(les_perim))
	ncb_m_widths = sum(np.array(ncb_m_widths))
	ncb_area = sum(np.array(ncb_area))
	ncb_perim = sum(np.array(ncb_perim))			
	bar_stats = [bar_m_widths,bar_area,bar_perim,dadate]
	gre_stats = [gre_m_widths,gre_area,gre_perim,dadate]
	lab_stats = [lab_m_widths,lab_area,lab_perim,dadate]
	les_stats = [les_m_widths,les_area,les_perim,dadate]
	ncb_stats = [ncb_m_widths,ncb_area,ncb_perim,dadate]
	
	elapsedtime = time.time() - time0
	print str(dadate)+' done in ',elapsedtime

# Aari polys
if 0:
	time0 = time.time()
	dadate = '20150113'
	aari = reader('Aari',dadate,hqm)
	# Now we have the aari dictionary with us, let's split it!
	adb,adg,adbg = aari.ad
#	SDA = SDA_poly(XM,YM,ZM)
#	over,under,DN,BM1,BM2 = SDA.over,SDA.under,SDA.DN,SDA.B1,SDA.B2
#	
#	bar_m_widths = []
#	bar_area = []
#	bar_perim = []
#	gre_m_widths = []
#	gre_area = []
#	gre_perim = []
#	lab_m_widths = []
#	lab_area = []
#	lab_perim = []
#	les_m_widths = []
#	les_area = []
#	les_perim = []
#	ncb_m_widths = []
#	ncb_area = []
#	ncb_perim = []
#	
#	poly_list=[]
#	for n,el in enumerate(over):
#		# classification of positive polygons
#		aod=poly_stat(el,'1',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
#		poly_list.append(aod)
#	try:
#		n
#	except NameError:
#		n = -1 
#	for n2,el in enumerate(under):
#		# classification of negative (n+1 for good enumeration)
#		n += 1 
#		aod=poly_stat(el,'0',BM2,BM1,DN,XM,YM,n,'DFP',basemap=hqm,PLOT=None,STCH=True)
#		poly_list.append(aod)
#	
#	bar_m_widths = sum(np.array(bar_m_widths))
#	bar_area = sum(np.array(bar_area))
#	bar_perim = sum(np.array(bar_perim))
#	gre_m_widths = sum(np.array(gre_m_widths))
#	gre_area = sum(np.array(gre_area))
#	gre_perim = sum(np.array(gre_perim))
#	lab_m_widths = sum(np.array(lab_m_widths))
#	lab_area = sum(np.array(lab_area))
#	lab_perim = sum(np.array(lab_perim))
#	les_m_widths = sum(np.array(les_m_widths))
#	les_area = sum(np.array(les_area))
#	les_perim = sum(np.array(les_perim))
#	ncb_m_widths = sum(np.array(ncb_m_widths))
#	ncb_area = sum(np.array(ncb_area))
#	ncb_perim = sum(np.array(ncb_perim))			
#	bar_stats = [bar_m_widths,bar_area,bar_perim,dadate]
#	gre_stats = [gre_m_widths,gre_area,gre_perim,dadate]
#	lab_stats = [lab_m_widths,lab_area,lab_perim,dadate]
#	les_stats = [les_m_widths,les_area,les_perim,dadate]
#	ncb_stats = [ncb_m_widths,ncb_area,ncb_perim,dadate]
#	
#	elapsedtime = time.time() - time0
#	print str(dadate)+' done in ',elapsedtime

