############################################################################
# Calling in all the needed packages
from netCDF4 import Dataset
import sys,os
import time
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
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

sys.path.append('../py_funs')
import mod_reading as Mrdg
import Laplace_eqn_solution as Leqs
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
			self.filname = name+'_'+dadate
			if name == 'Osisaf':
				 self.XO,self.YO,self.ZO = self._read_osi_(dadate,basemap) 
			elif name == 'Model':
				 self.X,self.Y,self.ZC,self.ZD = self._read_mdl_(dadate,basemap)
			elif name == 'Aari':
				 self._read_aari_(dadate,basemap)
	 	 return
	 
	 def _read_osi_(self,dadate,basemap):
	 	 # Read in OSI_SAF file
	 	 outdir = #TODO
	 	 ncfil = outdir+'ice_conc_nh_polstere-100_multi_'+dadate+'.nc'
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
			# Read TP4arch_wav
			outdir = #TODO
			ncfil = outdir+'TP4archv'+dadate+'.nc'
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

	 def binary_mod(data,thresh):
			# binary maker
			ndata = np.copy(data)
			ndata[ndata<thresh] = 0
			ndata[ndata>=thresh] = 1	
			thenans = np.isnan(ndata)
			ndata[thenans] = 0
			return(ndata)
			
	 def closing(self,DN):
			# apply closing to avoid small polynias and clean up a little
			kernel = np.ones((3,3),np.uint8)
			DN_cl = closing(DN,kernel)
			return(DN_cl)
			
	 def mapper(self,DIFF,OUT,INS):
			# getting back the terrain lost during binary_mod (for contour finding reasons)
			DN = DIFF
			ZO = OUT
			ZI = INS
			thenan = np.isnan(ZO)
			thenan2 = np.isnan(ZI)
			DN[thenan] = None
			DN[thenan2] = None
			self.difference = DN
			return
 
	 def _poly_maker_(self):
			if name == 'Osisaf':
				 ZO = self.ZO
				 X,Y = self.XO,self.YO
				 ic80 = self.binary_mod(ZO,.8)
				 ic80m = self.closing(ic80)
				 ic15 = self.binary_mod(ZO,.15)
				 ic15m = self.closing(ic15)
				 icp = ic15 - ic180
				 icpm = ic15m - ic80m
				 diff = icp
				 diffm = icpm
				 cont = msr.find_contours(diff,.5)
				 for n,en in enumerate(cont):
						aod

			if name == 'Model':
				 ZC = self.ZC
				 ZD = self.ZD
				 X,Y = self.X,self.Y
				 ic80 = self.binary_mod(ZC,.8)
				 ic80m = self.closing(ic80)
				 ic15 = self.binary_mod(ZC,.15)
				 ic15m = self.closing(ic15)
				 icp = ic80 - ic15
				 icpm = ic80m - ic15m
				 diff = icp
				 diffm = icpm
				 fsIP = self.binary_mod(ZC,300)
				 fsIPm = self.closing(fsIP)
				 fsOW = self.binary_mod(ZC,.01)
				 fsOWm = self.closing(fsOW)
				 fsp = fsIP - fsOW
				 fspm = fsIPm - fsOWm
				 fsdiff = fsp
				 fsdiffm = fspm
#			if name == 'Aari':
#				 # TODO
#				 ZO = self.osisaf
#				 ic80 = self.binary_mod(ZO,.8)
#				 ic80m = self.closing(ic80)
#				 ic15 = self.binary_mod(ZO,.15)
#				 ic15m = self.closing(ic15)
#				 icp = ic80 - ic15
#				 icpm = ic80m - ic15m
#				 diff = icp
#				 diffm = icpm

			# finding contours from the difference data map
			pos = msr.find_contours(DN,.5)
			neg = msr.find_contours(DN,-.5)
			# here happens the classification of the polygons - see aod_poly class
			poly_list=[]
			for n,el in enumerate(pos):
				 # classification of positive polygons
				 aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=0)
				 poly_list.append(aod)
			try:
				 n
			except NameError:
				 n = -1
			for n2,el in enumerate(neg):
				 # classification of negative (n+1 for good enumeration)
				 n += 1
				 aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=1)
				 poly_list.append(aod)

	 def _read_aari_(dadate,basemap):
	 # Read in AARI charts for Barents Sea
	 outdir = #TODO
	 ncfil = outdir+'aari_bar_'+dadate+'.txt'))
	 if ncfil == '':
			print('No Ice Chart for Barents sea in '+dadate+'\n')
	 else:
			region = 'Barents'
			icb = open(ncfil)
			icblist = []
			for n,en in enumerate(icb):
				 nen = en.split(';')
				 icblist.append(nen)
			icb_cont = np.array(icblist)
			# we want to cut first row (general info) and first column (number of points)
			icb_cont = icb_cont[1:,1:]
			# change 'in' as '1.' and 'out' as '0.'
			for el in icb_cont:
				 if el[1] == 'in':
						el[1] = 1
				 else:
						el[1] = 0
			# reads as string, we want floats
			icb_cont = icb_cont.astype(np.float)
			# find out number of polygons
			bpolyn = np.int(icb_cont[:,0].max())
			# split the array into polygons - creation of dict
			ad = {}
			for n in range(bpolyn):
				 poly = []
				 for m,em in enumerate(icb_cont):
						if em[0] == n+1:
							 poly.append(em)
				 ad["apoly"+str(n+1)] = np.array(poly)
	 
	 # Read in AARI charts for Greenland Sea
	 ncfil2 = 'aari_gre_'+dadate+'.txt'
	 if ncfil2 == '':
			print('No Ice Chart for Greenland sea in '+dadate+'\n')
	 else:
			region = 'Greenland'
			icg = open(ncfil2)
			icglist = []
			for n,en in enumerate(icg):
				 nen = en.split(';')
				 icglist.append(nen)
			icg_cont = np.array(icglist)
			# we want to cut first row (general info) and first column (number of points)
			icg_cont = icg_cont[1:,1:]
			# change 'in' as '1.' and 'out' as '0.'
			for el in icg_cont:
				 if el[1] == 'in':
						el[1] = 1
				 else:
						el[1] = 0
			# reads as string, we want floats
			icg_cont = icg_cont.astype(np.float)
			# find out number of polygons
			gpolyn = np.int(icg_cont[:,0].max())
			# split the array into polygons - creation of dict
			if ncfil == '':
				 ad = {}
			for n in range(gpolyn):
				 poly = []
				 for m,em in enumerate(icg_cont):
						if em[0] == n+1:
							 poly.append(em)
				 if ncfil == '':
						ad["apoly"+str(n+1)] = np.array(poly)
				 else:
						ad["apoly"+str(n+bpolyn+1)] = np.array(poly)
	 
	 # Read in AARI charts for both Barents and Greenland
	 # sometimes the charts are united in case the some polygons belong to both regions
	 if ncfil == '' and ncfil2 == '':
			ncfil3 = ''.join( glob.glob('*'+dadate+'.txt'))
			if ncfil3 == '':
				 print('No Common (Barents+Greenland) Ice Chart in '+dadate+'\n')
			else:
				 region = 'Barents/Greenland'
				 icbg = open(ncfil3)
				 icbglist = []
				 for n,en in enumerate(icbg):
						nen = en.split(';')
						icbglist.append(nen)
				 icbg_cont = np.array(icbglist)
				 # we want to cut first row (general info) and first column (number of points)
				 icbg_cont = icbg_cont[1:,1:]
				 # change 'in' as '1.' and 'out' as '0.'
				 for el in icbg_cont:
						if el[1] == 'in':
							 el[1] = 1
						else:
							 el[1] = 0
				 # reads as string, we want floats
				 icbg_cont = icbg_cont.astype(np.float)
				 # find out number of polygons
				 bgpolyn = np.int(icbg_cont[:,0].max())
				 ad = {}
				 for n in range(bgpolyn):
						poly = []
						for m,em in enumerate(icbg_cont):
							 if em[0] == n+1:
									poly.append(em)
						ad["apoly"+str(n+1)] = np.array(poly)
	 return(ad)

	 for i in ad:
			apoly = aari_poly(i,region,ad[i][:,1],ad[i][:,2],ad[i][:,3],X2,Y2,lqm)
			apoly_list.append(apoly)

############################################################################
# MIZ width
# we want to analyze the polygon and only if is class H get stats:
# 1) Width
# 2) Area
# NOTE if everything goes as planned we should hav only one class H per
# region hence the actual MIZ

class MIZ_poly:

	 # gets single contour and sorts it
	 def __init__(self,cont,OUT,INS,DIFF,X,Y,polygon_status=1):
			self.polygon_status = polygon_status
			# class definition
			if len(cont) > 200:
				 self.polygon_class = 'H'
				 self.ij_list = cont
				 self._ij2xy_(cont,X,Y)
				 self._split_cont(OUT,INS)
			elif len(cont) > 100 and len(cont) <= 200:
				 self.polygon_class = 'B'
				 self.ij_list = cont
				 self._ij2xy_(cont,X,Y)
				 self._split_cont(OUT,INS)
			elif len(cont) > 30 and len(cont) <= 100:
				 self.polygon_class = 'M'
				 self.ij_list = cont
				 self._ij2xy_(cont,X,Y)
				 self._split_cont(OUT,INS)
			elif len(cont) <= 30:
				 self.polygon_class = 'S'
				 self.ij_list = cont
				 self._ij2xy_(cont,X,Y)
				 self._split_cont(OUT,INS)
			return
	
	def _ij2xy_(self,cont,X,Y):
			# changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
			x = []
			y = []
			xvec =	range(len(cont[:,0]))
			for n,en in enumerate(xvec):
				 en	=	X[cont[n,0]][cont[n,1]]
				 x.append(en)
			yvec =	range(len(cont[:,1]))
			for n,en in enumerate(yvec):
				 en	=	Y[cont[n,0]][cont[n,1]]
				 y.append(en)
			xy_list = zip(x,y)
			xy_list = np.array(xy_list)
			self.xy_list = xy_list
			return
	
	 def _split_cont(self,OUT,INS,basemap):
			# this function find the different contours
			# NOTE if a point has non integer i coordinate is going to be a vertical edge,
			# if has non integer j coordinate is going to be a horizontal edge hence the use
			# of different arounds depending on the indexes of the point
			# NOTE if the polygon is an OVERESTIMATION the contour finding is inverted
			vs = self.ij_list # list of (i,j) pixel indices
			npoly = self.polygon_name
			in_cont = []
			out_cont = []
			ukn_cont = []
			func_vals= []
			func_mod = 0	# value of func_vals at model ice edge
			func_osi = 1	# value of func_vals at OSISAF ice edge
			func_unk = 2	# value of func_vals at other parts of contour
		
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
							 ukn_cont.append(el)
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
							 ukn_cont.append(el)
							 func_val=func_unk
						func_vals.append(func_val)

			in_cont = np.array(in_cont)
			out_cont = np.array(out_cont)
			ukn_cont = np.array(ukn_cont)
			func_vals = np.array(func_vals)
			self.in_ice_edge = in_cont
			self.out_ice_edge = out_cont
			self.unknown_edge = ukn_cont
			self.function_vals = func_vals
			return

	 def _lon_lat_(self):
			xy_list = self.xy_list
			self.lon_list,self.lat_list=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
			
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
			self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
			
			# defining the region
			centroid_longitude = self.centroid_longitude
			centroid_latitude = self.centroid_latitude
			if centroid_longitude < 80 and centroid_longitude > 8:
				self.polygon_region = 'Barents_Kara_sea'
			elif centroid_longitude < 8 and centroid_longitude > -44:
				self.polygon_region = 'Greenland_sea'
			elif centroid_longitude < -44 and centroid_longitude > -65:
				self.polygon_region = 'Labrador_sea'
			else:
				self.polygon_region = 'Misc.'
			
			# Calculating lon/lat for model,osisaf,unknown contours and distances
			inside_contour = self.inside_contour
			outside_contour	= self.outside_contour
			unknown_contour	= self.unknown_contour
			if len(inside_contour) != 0:
				 self.inside_contour_lon,self.inside_contour_lat = basemap(inside_contour[:,0],inside_contour[:,1],inverse=True)
			if len(outside_contour) != 0:
				 self.outside_contour_lon,self.outside_contour_lat = basemap(outside_contour[:,0],outside_contour[:,1],inverse=True)
			if len(unknown_contour) != 0:
				 self.unknown_contour_lon,self.unknown_contour_lat = basemap(unknown_contour[:,0],unknown_contour[:,1],inverse=True)
			return
			
	 def _laplacian_solution_(self,basemap):
			name = self.name
			lon = self.lon_list
			lat = self.lat_list
			xy = self.xy_list
			fval = self.f_vals
			results = Leqs.get_MIZ_widths(lon,fval,name=name,fig_outdir=dadate,basemap=basemap,xy_coords2=xy)
			self.AI = results[0]
			self.fun_sol = results[1]
			self.stream = results[2]
			return
############################################################################
# In this class every polygon already classified will be analyzed and stats
# will be saved and printed with the stat_chart() function
class aod_stats:
	def __init__(self,aod,dadate,basemap=None):
		self.ij_list = aod.ij_list
		self.xy_list = aod.xy_list
		self.number_of_points = len(aod.ij_list)
		self.polygon_status = aod.status
		self.name = aod.polygon_name
		self.pclass = aod.polygon_class
		if MODEL2MODEL:
			self.inside_contour = aod.inside_ice_edge
			self.outside_contour = aod.outside_ice_edge
		else:
			self.model_contour = aod.model_ice_edge
			self.osisaf_contour = aod.osisaf_ice_edge
		self.unknown_contour = aod.unknown_edge
		self.f_vals = aod.function_vals

		#get euclidean area
		#TODO get area on sphere?
		self.funtion_poly_area()
		
		#get euclidean perimeter
		#(geometry_planar.curve_info could be used)
		self.function_poly_perimeter()
		
		# distances between contours set self.widths
		self.function_dist_edges()

		if basemap is not None:
			xy_list = self.xy_list
			self.lon_list,self.lat_list=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
			
			# Perimeter mean - used to check diff with centroid 
			xmean	=	np.mean(xy_list[:,0])
			ymean	=	np.mean(xy_list[:,1])
			self.mean_perimeter_lon,self.mean_perimeter_lat	=	basemap(xmean,ymean,inverse=True)
			
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
			A	= (0.5)*B 
			CX	= (1/float(6*A))*DX
			CY	= (1/float(6*A))*DY
			self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
			# defining the region
			centroid_longitude = self.centroid_longitude
			centroid_latitude = self.centroid_latitude
			if centroid_longitude < 80 and centroid_longitude > 8:
				self.polygon_region = 'Barents_Kara_sea'
			elif centroid_longitude < 8 and centroid_longitude > -44:
				self.polygon_region = 'Greenland_sea'
      elif centroid_longitude < -44 and (centroid_longitude > -65 and centroid_latitude < 68):
				self.polygon_region = 'Labrador_sea'
			else:
				self.polygon_region = 'Misc.'

			# Calculating lon/lat for model,osisaf,unknown contours and distances
			if MODEL2MODEL:
				inside_contour = self.inside_contour
				outside_contour	= self.outside_contour
				unknown_contour	= self.unknown_contour
				dist_in2out = self.widths_in2out
				dist_out2in = self.widths_out2in
				if len(inside_contour) != 0:
					self.inside_contour_lon,self.inside_contour_lat = basemap(inside_contour[:,0],inside_contour[:,1],inverse=True)
				if len(outside_contour) != 0:
					self.outside_contour_lon,self.outside_contour_lat = basemap(outside_contour[:,0],outside_contour[:,1],inverse=True)
				if len(unknown_contour) != 0:
					self.unknown_contour_lon,self.unknown_contour_lat = basemap(unknown_contour[:,0],unknown_contour[:,1],inverse=True)
				distlonlat	= []
				if dist_in2out != [0,0,0,0,0]:
					for n,en in enumerate(dist_in2out):
						x1,y1 = basemap(en[1],en[2],inverse=True)
						x2,y2 = basemap(en[3],en[4],inverse=True)
						distlonlat.append([en[0],x1,y1,x2,y2])
					self.width_in2out_lonlat = distlonlat
				if dist_out2in != [0,0,0,0,0]:
					for n,en in enumerate(dist_out2in):
						x1,y1 = basemap(en[1],en[2],inverse=True)
						x2,y2 = basemap(en[3],en[4],inverse=True)
						distlonlat.append([en[0],x1,y1,x2,y2])
					self.width_out2in_lonlat = distlonlat
			else:
				inside_contour	= self.model_contour
				outside_contour = self.osisaf_contour
				unknown_contour	= self.unknown_contour
				dist_mdl2osi = self.widths_mdl2osi
				dist_osi2mdl = self.widths_osi2mdl
				if len(inside_contour) != 0:
					self.model_contour_lon,self.model_contour_lat = basemap(inside_contour[:,0],inside_contour[:,1],inverse=True)
				if len(outside_contour) != 0:
					self.osisaf_contour_lon,self.osisaf_contour_lat = basemap(outside_contour[:,0],outside_contour[:,1],inverse=True)
				if len(unknown_contour) != 0:
					self.unknown_contour_lon,self.unknown_contour_lat = basemap(unknown_contour[:,0],unknown_contour[:,1],inverse=True)
				distlonlat	= []
				if dist_mdl2osi != [0,0,0,0,0]:
					for n,en in enumerate(dist_mdl2osi):
						x1,y1 = basemap(en[1],en[2],inverse=True)
						x2,y2 = basemap(en[3],en[4],inverse=True)
						distlonlat.append([en[0],x1,y1,x2,y2])
					self.width_mdl2osi_lonlat = distlonlat
				if dist_osi2mdl != [0,0,0,0,0]:
					for n,en in enumerate(dist_osi2mdl):
						x1,y1 = basemap(en[1],en[2],inverse=True)
						x2,y2 = basemap(en[3],en[4],inverse=True)
						distlonlat.append([en[0],x1,y1,x2,y2])
					self.width_osi2mdl_lonlat = distlonlat
		return

	def funtion_poly_area(self):
	  # Calculating the area of irregular polygon (from perimeter)
		vs = self.ij_list
		a	= 0
		x0,y0 = vs[0]
		for [x1,y1] in vs[1:]:
			dx = x1-x0
			dy = y1-y0
			a += 0.5*(y0*dx - x0*dy)
			x0 = x1
			y0 = y1
		self.area_euclidean	=	a*100
		return

	def function_poly_perimeter(self):
		# Calculating perimeter in xy coordinates (unit = 10km)
		xy		= self.ij_list
		perim	= 0
		for n in range(len(xy)-1):
			perim += np.sqrt(pow(xy[n+1][0]-xy[n][0],2)+pow(xy[n+1][1]-xy[n][1],2))
		self.perimeter	= perim*10
		return

	def function_dist_edges(self):
		# Calculating the min distance from a model point to any osisaf point
		# and vice-versa (osisaf to model)
		# same script applied for the Model2Model(in&out) products - see legends
		if MODEL2MODEL:
			inside_contour = self.inside_contour
			outside_contour	= self.outside_contour
		else:
			inside_contour	= self.model_contour
			outside_contour = self.osisaf_contour
		unknown_contour = self.unknown_contour
		tcont = self.ij_list
		UKW   = 'Unknown contours < 20%'
		DMW   = 'Contour difference < 40%'
		width_in2out = []
		width_out2in = []
		if len(inside_contour) == 0:
			DMW = 'Only OSISAF data, MODEL not present'
			if MODEL2MODEL:
				if FLOES:
					DMW = 'Only Ice Pack, Open Water not present'
				else:
					DMW = 'Only 80% edge, 15% not present'
			self.widths_in2out	= [0,0,0,0,0]
			self.widths_out2in	= [0,0,0,0,0]
		elif len(outside_contour) == 0:
			DMW = 'Only MODEL data, OSISAF not present'
			if MODEL2MODEL:
				if FLOES:
					DMW = 'Only Open Water, Ice Pack not present'
				else:
					DMW = 'Only 15% edge, 80% not present'
			self.widths_in2out	= [0,0,0,0,0]
			self.widths_out2in	= [0,0,0,0,0]
		else:
			ukn = 100*(len(unknown_contour)/float(len(tcont)))
			dmo = 100*(abs(len(inside_contour)-len(outside_contour))/float(len(tcont)))
			if ukn >= 20:
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
			if MODEL2MODEL:
				self.widths_in2out = width_in2out
				self.widths_out2in = width_out2in
			else:
				self.widths_mdl2osi = width_in2out
				self.widths_osi2mdl = width_out2in
		self.UKW		= UKW
 		self.DMW		= DMW
 		return

	def function_laplacian_solution(self,dadate,basemap):
		print('This function will find the Laplacian and Stream solutions of the polygon...\n')
		name = self.name
		lon = self.lon_list
		lat = self.lat_list
		xy = self.xy_list
		fval = self.f_vals
		results = Leqs.get_MIZ_widths(lon,fval,name=name,fig_outdir=dadate,basemap=basemap,xy_coords2=xy)
		self.AI = results[0]
		self.fun_sol = results[1]
		self.stream = results[2]
		return

	def quick_figure(self):
		DN = aod.difference
		DMW = self.DMW
		UKW = self.UKW
		ij = self.ij_list
		if MODEL2MODEL:
			inside_contour = self.inside_contour
			outside_contour	= self.outside_contour
		else:
			inside_contour	= self.model_contour
			outside_contour = self.osisaf_contour
		unknown_contour = self.unknown_contour
		if MODEL2MODEL:
			dist_in2out = self.widths_in2out
			dist_out2in = self.widths_out2in
		else:
			dist_in2out = self.widths_mdl2osi
			dist_out2in = self.widths_osi2mdl
		x1,x2,y1,y2 = np.min(ij[:,1])-10,np.max(ij[:,1])+10,np.min(ij[:,0])-10,np.max(ij[:,0])+10
		plt.imshow(DN,interpolation='nearest',cmap='winter')
		plt.axis([x1,x2,y2,y1])
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
		plt.show(False)
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

###########################################################################
# finding contours from the difference data map
pos = msr.find_contours(DN,.5)
neg = msr.find_contours(DN,-.5)
###########################################################################
# NOTE here happens the classification of the polygons - see aod_poly class
poly_list=[]
for n,el in enumerate(pos):
	 # classification of positive polygons
	 aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=0)
	 poly_list.append(aod)
try:
	 n
except NameError:
	 n = -1
for n2,el in enumerate(neg):
	 # classification of negative (n+1 for good enumeration)
	 n += 1
	 aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=1)
	 poly_list.append(aod)

## changing name for every polygon, easier to work with
# NOTE decided to move all the info to the stat class, keep this silent
#d = {}
#for x in range(len(poly_list)):
#	 d[poly_list[x].polygon_name+'ID']=poly_list[x]
#for key,value in sorted(d.items()):
#	 globals()[key] = value

elapsedtime = time.time() - time0
print 'Polygon identification done in ',elapsedtime
print ''
###########################################################################

###########################################################################
# Fancy stats & aod_stats_production
# NOTE the following stats are just an overview of the number of polygons per region
# it has no effect whatsoever on the final result of the script but it's nice to see
# what is happening and how are the polygons distributed on the regions
# NOTE inside the fancy stats there is the stats production - see class aod_stats
print 'Calculating statistics for every polygon...'
print ''
scount = 0
mcount = 0
bcount = 0
hcount = 0
sbk_sea = 0
sgr_sea = 0
slb_sea = 0
smiscel = 0
mbk_sea = 0
mgr_sea = 0
mlb_sea = 0
mmiscel = 0
bbk_sea = 0
bgr_sea = 0
blb_sea = 0
bmiscel = 0
hbk_sea = 0
hgr_sea = 0
hlb_sea = 0
hmiscel = 0
poly_stat_list = []
for el in poly_list:
	 if el.polygon_class == 'S':
			scount += 1
			print('Name - '+str(el.polygon_name))
			print('Class - '+str(el.polygon_class))
			poly = aod_stats(el,dadate,lqm)
			if poly.polygon_region == 'Barents_Kara_sea':
				 sbk_sea += 1
			elif poly.polygon_region == 'Greenland_sea':
				 sgr_sea += 1
			elif poly.polygon_region == 'Labrador_sea':
				 slb_sea += 1
			elif poly.polygon_region == 'Misc.':
				 smiscel += 1
			print('Region - '+str(poly.polygon_region))
			print ''
	 elif el.polygon_class == 'M':
			mcount += 1
			print('Name - '+str(el.polygon_name))
			print('Class - '+str(el.polygon_class))
			poly = aod_stats(el,dadate,lqm)
			if poly.polygon_region == 'Barents_Kara_sea':
				 mbk_sea += 1
			elif poly.polygon_region == 'Greenland_sea':
				 mgr_sea += 1
			elif poly.polygon_region == 'Labrador_sea':
				 mlb_sea += 1
			elif poly.polygon_region == 'Misc.':
				 mmiscel += 1
			print('Region - '+str(poly.polygon_region))
			print ''
	 elif el.polygon_class == 'B':
			print('Name - '+str(el.polygon_name))
			print('Class - '+str(el.polygon_class))
			bcount += 1
			poly = aod_stats(el,dadate,lqm)
			if poly.polygon_region == 'Barents_Kara_sea':
	 			bbk_sea += 1
	 		elif poly.polygon_region == 'Greenland_sea':
	 			bgr_sea += 1
	 		elif poly.polygon_region == 'Labrador_sea':
	 			blb_sea += 1
	 		elif poly.polygon_region == 'Misc.':
	 			bmiscel += 1
	 		print('Region - '+str(poly.polygon_region))
	 		print ''
	 elif el.polygon_class == 'H':
	 		print('Name - '+str(el.polygon_name))
	 		print('Class - '+str(el.polygon_class))
	 		hcount += 1
	 		poly = aod_stats(el,dadate,lqm)
	 		if poly.polygon_region == 'Barents_Kara_sea':
	 			hbk_sea += 1
	 		elif poly.polygon_region == 'Greenland_sea':
	 			hgr_sea += 1
	 		elif poly.polygon_region == 'Labrador_sea':
	 			hlb_sea += 1
	 		elif poly.polygon_region == 'Misc.':
	 			hmiscel += 1
	 		print('Region - '+str(poly.polygon_region))
	 		print ''
	 poly_stat_list.append(poly)

# changing name to every polygon_stat, easier to work with
d = {}
for x in range(len(poly_stat_list)):
	d[poly_stat_list[x].name]=poly_stat_list[x]
for key,value in sorted(d.items()):
	globals()[key] = value

# printing out the localization stats produced above
print '..statistics done.'
print ''
print 'Total number of Polygons	=',len(poly_list)
print ''
print 'Small Polygons (30 points or less) = ',scount
print 'Of which:	'
print sbk_sea,' in the Barents_Kara sea'
print sgr_sea,' in the Greenland sea'
print slb_sea,' in the Labrador sea'
print smiscel,' in the other seas'
print ''
print 'Medium Polygons (100 points or less) = ',mcount
print 'Of which:	'
print mbk_sea,' in the Barents_Kara sea'
print mgr_sea,' in the Greenland sea'
print mlb_sea,' in the Labrador sea'
print mmiscel,' in the other seas'
print ''
print 'Big Polygons (200 points or less) = ',bcount
print 'Of which:	'
print bbk_sea,' in the Barents_Kara sea'
print bgr_sea,' in the Greenland sea'
print blb_sea,' in the Labrador sea'
print bmiscel,' in the other seas'
print ''
print 'Huge Polygons (more than 200 points) = ',hcount
print 'Of which:	'
print hbk_sea,' in the Barents_Kara sea'
print hgr_sea,' in the Greenland sea'
print hlb_sea,' in the Labrador sea'
print hmiscel,' in the other seas'
print ''

elapsedtime = time.time() - time0
print 'Job done in ',elapsedtime

###########################################################################
# NOTE the user get to choose if the stat_charts are saved
# the script will automatically save the charts in regional directories
# all of which will be moved in the directory with the date of the products
# ./outputs/aod/YYYYMMDD/Barents&Kara,Greenland,Labrador,Miscellaneous
if FIGURE:
	print 'Saving figures...'
	print ''
	time1 = time.time()
	#figure_save(X2,Y2,DN,'DN'+dadate,hqm)
	for poly in poly_stat_list:
		poly.stat_chart(save=True)
	# moving all the charts in the proper dir
	fin_dir = './outputs/aod/'+str(dadate)+'/'
	if not os.path.exists(fin_dir):
		os.makedirs(fin_dir)
	for filname in glob.glob(os.path.join('*.png')):
		shutil.move(filname,fin_dir)
	if os.path.exists('Model_Osisaf'):
		shutil.move('Model_Osisaf',fin_dir)
	if os.path.exists('Model2Model_Floe_Size'):
		shutil.move('Model2Model_Floe_Size',fin_dir)
	if os.path.exists('Model2Model_Ice_Conc'):
		shutil.move('Model2Model_Ice_Conc',fin_dir)
	elapsedtime1 = time.time() - time1
	elapsedtime = time.time() - time0
	print 'Figures saved and moved in ',elapsedtime1
	print ''
	print 'Time needed for the job ',elapsedtime

###########################################################################
# TODO if the user choose to save, it can save the polygons H & B as npz 
#	if MODEL2MODEL:
#		if FLOES:
#			outdir='/./outputs/aod/'+str(dadate)+'/Model2Model_Floe_Size/npz/'		
#			if not os.path.exists(outdir):             		# save a poly for testing
#				os.mkdir(outdir)
#			for poly in poly_stat_list:
#				if poly.pclass == 'H' or poly.pclass == 'B':
#					npfil = outdir+poly.pname+'.npz'
#					print('saving to '+npfil)
#					#
#					xy_coords = poly.xy_list
#					xy_coords2 = [tuple(xyc) for xyc in xy_coords]
#					fvals2 = 1*poly.function_vals
#					# save file
#					np.savez(npfil,xy=xy_coords,func_vals=fvals2)
#		else:
#			outdir='/./outputs/aod/'+str(dadate)+'/Model2Model_Ice_Conc/npz/'		
#			if not os.path.exists(outdir):             		# save a poly for testing
#				os.mkdir(outdir)
#			for poly in poly_stat_list:
#				if poly.pclass == 'H' or poly.pclass == 'B':
#					npfil = outdir+poly.pname+'.npz'
#					print('saving to '+npfil)
#					#
#					xy_coords = poly.xy_list
#					xy_coords2 = [tuple(xyc) for xyc in xy_coords]
#					fvals2 = 1*poly.function_vals
#					# save file
#					np.savez(npfil,xy=xy_coords,func_vals=fvals2)	
#	else:
#			outdir='/./outputs/aod/'+str(dadate)+'/Model_Osisaf/npz/'		
#			if not os.path.exists(outdir):             		# save a poly for testing
#				os.mkdir(outdir)
#			for poly in poly_stat_list:
#				if poly.pclass == 'H' or poly.pclass == 'B':
#					npfil = outdir+poly.pname+'.npz'
#					print('saving to '+npfil)
#					#
#					xy_coords = poly.xy_list
#					xy_coords2 = [tuple(xyc) for xyc in xy_coords]
#					fvals2 = 1*poly.function_vals
#					# save file
#					np.savez(npfil,xy=xy_coords,func_vals=fvals2)	
###########################################################################

###########################################################################
# This part of the script has the purpose of saving npz polygons for testing
if 0:
	outdir='../python_tutorial/npz'
	if not os.path.exists(outdir):
		os.mkdir(outdir)

	if MODEL2MODEL:
		# save a poly for testing
		nt = 70
		npfil = outdir+'/poly'+str(nt)+'.npz'
		print('saving to '+npfil)
		#
		poly_test = poly_list[nt]
		xy_coords = poly_test.xy_list
		xy_coords2 = [tuple(xyc) for xyc in xy_coords]
		fvals2 = 1*poly_test.function_vals
		
		# save file
		np.savez(npfil,xy=xy_coords,func_vals=fvals2)
	else:
		# save a poly for testing
		nt = 55
		npfil = outdir+'/poly'+str(nt)+'.npz'
		print('saving to '+npfil)
		#
		poly_test = poly_list[nt]
		xy_coords = poly_test.xy_list
		xy_coords2 = [tuple(xyc) for xyc in xy_coords]
		fvals2 = 1*poly_test.function_vals
		
		# save file
		np.savez(npfil,xy=xy_coords,func_vals=fvals2)
			
