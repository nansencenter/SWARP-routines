############################################################################
# AREA OF DISAGREEMENT - for any info/questions - timothy.williams@nersc.no
#                                               - gcmdnt90@gmail.com
############################################################################

############################################################################
# The following script will:
# 1) Reproject any dataset (as long as Stereographical) into the OSISAF
#			arctic sea ice dataset grid (http://www.osi-saf.org/)
# 2) Depending on the user's choice will study areas of disagreement between:
#		a) Ice concentrations of Model and Observational data - Ice cover except 
#				sparse ice (default 15% threshold) 
#		b) Ice concentration contours between the model - concentration dependent 
#				MIZ (default 15% to 80%) 
#		c) Floe size distribution contours between - Floe size dependent MIZ 
#				(default Dmax < 300 m)
# 3) Identify polygons with relative contours, localize and classify them
# 4) Running stats for every polygon - producing stat charts
# 5) TODO run Laplace script for every class H polygon
############################################################################

############################################################################
# INPUTS
# The script will automatically find any .nc file in the data directory and 
# distinguish between Model and Observational data
# NOTE the NetCDF files will have to be located into /data directory,
# to download the files from Hexagon repositories check ./data/fetch_data.sh
# NOTE This script is strictly bounded to the presence of ../py_funs/  and
# ./data/ directories, copy them properly in case of copying of the script
############################################################################

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
from skimage.morphology import opening, closing
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
# function used to make pcolor plots fancier
def finish_map(m):
	# *m is a basemap
	m.drawcoastlines()
	m.fillcontinents(color='gray')
	# draw parallels and meridians.
	m.drawparallels(np.arange(60.,91.,10.),\
	        labels=[True,False,True,True]) # labels = [left,right,top,bottom]
	m.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
	        labels=[True,False,False,True])
	m.drawmapboundary() # fill_color='aqua')
	return
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
def open_close(DD):
	kernel = np.ones((3,3),np.uint8)
	DN_op = np.copy(DD)
	DN_cl = np.copy(DD)
	DN_op = opening(DN_op,kernel)
	DN_cl = closing(DN_cl,kernel)
	return(DN_op,DN_cl)
############################################################################
# function that calculates eigenvectors and uses them to find sensible
# starting point for "extended" Laplacian width analisis
def start_point(data):
	i_vec = data[:,0]
	j_vec = data[:,1]
	i_vec = i_vec - np.mean(i_vec)
	j_vec = j_vec - np.mean(j_vec)
	ndata = [i_vec,j_vec]
	ndata = np.matrix(ndata)
	A = ndata * ndata.T
	L,D,R = np.linalg.svd(A)
	u = np.array(L[:,0]*1000)
	v = np.array(L[:,1]*1000)
	dist_pt = []
# TODO I have to decide which point to use to calculate the min distance!
# NOTE I could use the in_cont the out_cont or a new set o median pts!
	if DMW == 'Contour difference < 40%':
		for m,em in enumerate(			):
			dist1	= np.sqrt(pow(u[0]-em[0],2)+pow(u[1]-em[1],2))
			dist_pt.append([dist1,em[0],em[1]])
		idx,value = min(enumerate(dist_pt),key=itg(1))
		start_pt = [value[1],value[2]]
	else:
		for m,em in enumerate(			):
			dist1	= np.sqrt(pow(v[0]-em[0],2)+pow(v[1]-em[1],2))
			dist_pt.append([dist1,em[0],em[1]])
		idx,value = min(enumerate(dist_pt),key=itg(1))
		start_pt = [value[1],value[2]]
	return(start_pt)
############################################################################
# function used to create binary datasets - used in Model2Osisaf
def binary_mod(data1,data2,thresh):
	mdata = np.copy(data1)
	mdata[mdata<thresh] = 0
	mdata[mdata>=thresh] = 1
	odata = np.copy(data2)
	odata[odata<thresh] = 0
	odata[odata>=thresh] = 1	
	ddata = odata - mdata
	thenans = np.isnan(ddata)
	ddata[thenans] = 0
	return(mdata,odata,ddata)
############################################################################
# function used to create binary datasets - used in Model2Model_ice_conc
def binary_mod_2(data1,data2,thresh1,thresh2):
	mdata = np.copy(data1)
	mdata[mdata<thresh1] = 0
	mdata[mdata>=thresh1] = 1
	odata = np.copy(data2)
	odata[odata<thresh2] = 0
	odata[odata>=thresh2] = 1	
	ddata = odata - mdata
	thenans = np.isnan(ddata)
	ddata[thenans] = 0
	return(mdata,odata,ddata)
############################################################################
# function used to create binary datasets - used in Model2Model_floe_size
def binary_mod_3(data1,thresh):
	mdata = np.copy(data1)
	mdata[mdata>.01] = 1
	mdata[mdata<=.01] = 0
	odata = np.copy(data1)
	odata[odata<thresh] = 0
	odata[odata>=thresh] = 1
	ddata = odata - mdata
	thenans = np.isnan(ddata)
	ddata[thenans] = 0
	return(mdata,odata,ddata)
############################################################################
# function used to save a fancy version of the Difference Dataset
def figure_save(X,Y,Z,name,m):
	f = plt.figure()
	m.pcolor(X,Y,Z,cmap='winter')
	finish_map(m)
	fname = name+'.png'
	plt.colorbar()
	plt.title(name)
	plt.savefig(fname,format='png',dpi=1000)
	plt.close()
	f.clf()
	return
###########################################################################
# function to print initial informations about the script outputs, pretty
# much useless if not for major data errors. Why? Because Science.
def get_stats(data,name):
	stat0	= (data == 0).sum()
	stat1	= (data == 1).sum()
	stat2	= (data == -1).sum()
	stat	= data.size
	phit	= (stat0 / float(stat)) * 100
	pone	= (stat1 / float(stat)) * 100
	pmone	= (stat2 / float(stat)) * 100
	print ' ' 
	print 'Number of elements:	%d' %stat
	print 'Hit:			%d - %1.2f %%' %(stat0,phit)
	print 'Underprediction(+1):	%d - %1.2f %%' %(stat1,pone)
	print 'Overprediction(-1):	%d - %1.2f %%' %(stat2,pmone)
	print ' ' 
	return()

###########################################################################
# GETTING THE POLYGONS from AARI-Icecharts text files
class aari_poly:
	# gets single contour and sorts it
	def __init__(self,name,region,fvals,lon,lat,X,Y,basemap=None):
		self.lon = np.array(lon)
		self.lat = np.array(lat)
		self.polygon_name = name
		self.polygon_region = region
		# setting up lonlat_list as array
		lonlat = zip(self.lon,lat)
		lonlat_list = np.array(lonlat)
		self.lonlat_list = lonlat_list
		self.fvals = fvals
		# class definition
		if len(lonlat_list) <= 30:
			self.polygon_class = 'S'
		elif len(lonlat_list) > 30 and len(lonlat_list) <= 100:
			self.polygon_class = 'M'
		elif len(lonlat_list) > 100 and len(lonlat_list) <= 200:
			self.polygon_class = 'B'
		elif len(lonlat_list) > 200:
			self.polygon_class = 'H'
		if basemap is not None:
			x_list,y_list = basemap(lonlat_list[:,0],lonlat_list[:,1])
			xy_list = zip(x_list,y_list)
			self.xy_list = np.array(xy_list)
			#self._xy2ij_(X,Y)
			self._split_cont_()
		return

  # This can't be done because of different grid!
	def _xy2ij_(self,X,Y):
		xy_list = self.xy_list
		i = []
		j = []
		xh,xv = X.shape
		for n,en in enumerate(xy_list[:,0]):
			for m in range(xh):
				for l in range(xv):
					if en == X[m,l]:
						i.append(m)
		yh,yv = Y.shape
		for n,en in enumerate(xy_list[:,1]):
			for m in range(yh):
				for l in range(yv):
					if en == X[m,l]:
						j.append(l)
		ij_list = zip(i,j)
		self.ij_list = np.array(ij_list)
		return
	
	def _split_cont_(self):
		# this function find the different contours
		#ij = self.ij_list # list of (i,j) pixel indices
		xy = self.xy_list
		lonlat = self.lonlat_list
		fvals = self.fvals
		npoly = self.polygon_name
		#ij_in = []
		#ij_out = []
		xy_in = []
		xy_out = []
		lonlat_in = []
		lonlat_out = []

		for n,en in enumerate(fvals):
			if en == 1:
				#ij_in.append(ij_list[n,:])
				xy_in.append(xy[n,:])
				lonlat_in.append(lonlat[n,:])
			else:
				#ij_out.append(ij_list[n,:])
				xy_out.append(xy[n,:])
				lonlat_out.append(lonlat[n,:])
		#self.ij_in_cont = ij_in
		#self.ij_out_cont = ij_out
		self.xy_in_cont = np.array(xy_in)
		self.xy_out_cont = np.array(xy_out)
		self.lonlat_in_cont = np.array(lonlat_in)
		self.lonlat_out_cont = np.array(lonlat_out)
		return

	def function_dist_edges(self):
		# Calculating the min distance from in point to any out point
		# and vice-versa (out to in)
		inside_contour = self.xy_in_cont
		outside_contour = self.xy_out_cont
		tcont = self.xy_list
		UKW   = 'Unknown contours < 20%'
		DMW   = 'Contour difference < 40%'
		width_in2out	= [0,0,0,0,0]
		width_out2in	= [0,0,0,0,0]
		if len(inside_contour) == 0:
			DMW = 'Only outside contour, inside condtour not present'
		elif len(outside_contour) == 0:
			DMW = 'Only inside contour, outside not present'
		else:
			width_in2out = []
			width_out2in = []
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
		self.widths_in2out = width_in2out
		self.widths_out2in = width_out2in
		self.UKW = UKW
 		self.DMW = DMW
 		return

	def function_laplacian_solution(self):
		print('This function will find the Laplacian and Stream solutions of the polygon...\n')
		name = self.name
		bmp = basemap
		ddt = dadate
		lon = self.lon
		lat = self.lat
		xy = self.xy_list
		fval = self.fvals
		results = Leqs.get_MIZ_widths(lon,fval,name=name,fig_outdir=ddt,basemap=bmp,xy_coords2=xy)
		self.AI = results[0]
		self.fun_sol = results[1]
		self.stream = results[2]
		return

	def get_stats(self):
		xy_list = self.xy_list
		iowidths = self.widths_in2out
		oiwidths = self.widths_out2in
		iomedian = np.median(iowidths)
		oimedian = np.median(oiwidths)
		# Centroid - Location in lon/lat of the centroid of the polygon
		DX = 0
		DY = 0
		B = 0
		for n in range(len(xy_list)-1):
			DX += ((xy_list[n][0]+xy_list[n+1][0])*\
				((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
			DY += ((xy_list[n][1]+xy_list[n+1][1])*\
				((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1])))  
			B	+= ((xy_list[n][0]*xy_list[n+1][1])-(xy_list[n+1][0]*xy_list[n][1]))
			A	= (0.5)*B 
			CX	= (1/float(6*A))*DX
			CY	= (1/float(6*A))*DY
		self.centroid_x,self.centroid_y = CX,CY
		self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
		# defining the region
		centroid_longitude = self.centroid_longitude
		bblone,bblate,bblonw,bblatw = baffin_bay()
		if centroid_longitude < 80 and centroid_longitude > 8:
		 	self.polygon_region = 'Barents_Kara_sea'
		elif centroid_longitude < 8 and centroid_longitude > -44:
		 	self.polygon_region = 'Greenland_sea'
		elif centroid_longitude < -44 and centroid_longitude > -90:
			for n,en in enumerate(bblatw[1:len(bblatw)]):
				if (centroid_latitude >= bblatw[n+1] and centroid_latitude <= bblatw[n]):
					if centroid_longitude >= bblonw[n]: 
						self.polygon_region = 'Labrador_sea'
						break
				else:
					self.polygon_region = 'Misc.'
		else:
		 	self.polygon_region = 'Misc.'
		return
	
############################################################################
# GETTING THE POLYGONS
# In this class we define polygons from the contour findings
# What does it do?
# 1) Classification of the Polygon
# 2) Storage of all the points in:
#			a) indexes
#			b) x&y coords 
#			c) Lon/Lat
# 3) Definition of different contours
class aod_poly:
	# gets single contour and sorts it
	def __init__(self,cont,OUT,INS,DIFF,X,Y,number,polygon_status=1):
		self.polygon_status = polygon_status
		self.polygon_name = 'polygon'+str(number)
		# class definition
		if len(cont) <= 30:
			self.polygon_class = 'S'
		elif len(cont) > 30 and len(cont) <= 100:
			self.polygon_class = 'M'
		elif len(cont) > 100 and len(cont) <= 200:
			self.polygon_class = 'B'
		elif len(cont) > 200:
			self.polygon_class = 'H'
		self.ij_list = cont
		self._diff_with_terrain_(DIFF,OUT,INS) 
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
	
	def _diff_with_terrain_(self,DIFF,OUT,INS):
		# getting back the terrain lost during binary_mod (for contour finding reasons)
		DN = DIFF
		ZO = OUT
		ZM = INS
		thenan = np.isnan(ZO)
		thenan2 = np.isnan(ZM)
		DN[thenan] = None
		DN[thenan2] = None
		self.difference = DN
		return
	
	def _split_cont(self,OUT,INS):
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
		
		if self.polygon_status==1:
			# if polygon is negative - overestimation of ice
			self.status='Underestimate'
			for n,el in enumerate(vs):
				#getting all the neighbours
				around2 = ((el[0]+.5,el[1]),(el[0]-.5,el[1]))
				around1 = ((el[0],el[1]+.5),(el[0],el[1]-.5))
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
						ukn_cont.append(el)
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
						ukn_cont.append(el)
						func_val=func_unk
					func_vals.append(func_val)
		else:
			# if polygon is positive - underestimation of ice
			self.status='Overestimate'
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
		func_vals	= np.array(func_vals)
		# include sorted contours as port of object
		if MODEL2MODEL:
			self.inside_ice_edge = in_cont
			self.outside_ice_edge = out_cont
		else:
			self.model_ice_edge	= in_cont
			self.osisaf_ice_edge = out_cont
		self.unknown_edge = ukn_cont
		self.function_vals = func_vals
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
			self.centroid_x,self.centroid_y = CX,CY
			self.centroid_longitude,self.centroid_latitude = basemap(CX,CY,inverse=True)
			# defining the region
			centroid_longitude = self.centroid_longitude
			centroid_latitude = self.centroid_latitude
			bblone,bblate,bblonw,bblatw = baffin_bay()
			if centroid_longitude < 80 and centroid_longitude > 8:
				self.polygon_region = 'Barents_Kara_sea'
			elif centroid_longitude < 8 and centroid_longitude > -44:
				self.polygon_region = 'Greenland_sea'
			elif centroid_longitude < -44 and centroid_longitude > -90:
				for n,en in enumerate(bblatw[1:len(bblatw)]):
					if (centroid_latitude >= bblatw[n+1] and centroid_latitude <= bblatw[n]):
						if centroid_longitude >= bblonw[n]: 
							self.polygon_region = 'Labrador_sea'
							break
					else:
						self.polygon_region = 'Misc.'
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
		if len(inside_contour) == 0:
			DMW = 'Only OSISAF data, MODEL not present'
			if MODEL2MODEL:
				if FLOES:
					DMW = 'Only Ice Pack, Open Water not present'
				else:
					DMW = 'Only 80% edge, 15% not present'
			width_in2out = [0,0,0,0,0]
			width_out2in = [0,0,0,0,0]
		elif len(outside_contour) == 0:
			DMW = 'Only MODEL data, OSISAF not present'
			if MODEL2MODEL:
				if FLOES:
					DMW = 'Only Open Water, Ice Pack not present'
				else:
					DMW = 'Only 15% edge, 80% not present'
			width_in2out = [0,0,0,0,0]
			width_out2in = [0,0,0,0,0]
		else:
			width_in2out = []
			width_out2in = []
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

	def function_laplacian_solution(self):
		print('This function will find the Laplacian and Stream solutions of the polygon...\n')
		name = self.name
		f_out = dadate
		basemap = hqm
		lon = self.lon_list
		lat = self.lat_list
		xy = self.xy_list
		fval = self.f_vals
		results = Leqs.get_MIZ_widths(lon,lat,fval,name=name,fig_outdir=f_out,basemap=basemap,xy_coords2=xy)
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

	def stat_chart(self,save=False):
		pname = self.name
		pclass = self.pclass
		pstat = self.polygon_status
		region = self.polygon_region
		print ''
		print 'Statistic Chart for '+str(pname)
		print ''
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
		clon = '%1.2f' %self.centroid_longitude
		clat = '%1.2f' %self.centroid_latitude
		clonlat = '{0}/{1}'.format(clon,clat)
		area = self.area_euclidean
		area = '%1.4e' %area 
		perim = self.perimeter
		perim = '%1.4e' %perim
		if dist_in2out != [0,0,0,0,0]:
			dist_in = np.array(dist_in2out)
			# changing units from decakilometer to kilometer
			dist_in = dist_in[:,0]*10
			dist_in = np.median(dist_in)
			dist_in = '%1.2f' %dist_in
		else:
			dist_in = 'NaN'
		if dist_out2in != [0,0,0,0,0]:
			dist_out = np.array(dist_out2in)
			# changing units from decakilometer to kilometer
			dist_out = dist_out[:,0]*10
			dist_out = np.median(dist_out)
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
		if dist_in2out != [0,0,0,0,0]:
			for n,en in enumerate(dist_in2out):
				polyf.plot([en[2],en[4]],[en[1],en[3]],color='black',alpha=0.1)
		if dist_out2in != [0,0,0,0,0]:
			for n,en in enumerate(dist_out2in):
				polyf.plot([en[2],en[4]],[en[1],en[3]],color='magenta',alpha=0.1)

		if MODEL2MODEL:
			if FLOES:
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
					valid_class = 'Model2Model_Floe_Size'
					if not os.path.exists(valid_class):
						os.mkdir(valid_class)
					if not os.path.exists(valid_class+'/'+region):
						os.mkdir(valid_class+'/'+region)
					fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
				print 'Statistic chart done for '+str(pname)
			else:
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
					valid_class = 'Model2Model_Ice_Conc'
					if not os.path.exists(valid_class):
						os.mkdir(valid_class)
					if not os.path.exists(valid_class+'/'+region):
						os.mkdir(valid_class+'/'+region)
					fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
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
			      'Model Cont.','Osisaf Cont.','Unknown Cont.','Dist. Mdl to Osi', \
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
				valid_class = 'Model_Osisaf'
				if not os.path.exists(valid_class):
					os.mkdir(valid_class)
				if not os.path.exists(valid_class+'/'+region):
					os.mkdir(valid_class+'/'+region)
				fig.savefig(valid_class+'/'+region+'/'+pclass+'_'+pname+'.png',bbox_inches='tight')
			print 'Statistic chart done for '+str(pname)

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
# Getting a nice blank space before user's input 
if 0:
   print ''
   # User's inputs
   FIGURE = raw_input('[1] to plot&save all figures, [Enter] to run without plot&save	') 
   print ''
   MODEL2MODEL = raw_input('[1] for model2model, [ENTER] to model2osisaf	') 
   print ''
   if MODEL2MODEL:
           FLOES = raw_input('[1] for Floe Size Distribution, [ENTER] for Ice concentration	') 
           print ''
   else:
           FLOES = []
else:
   FIGURE   = ''
   MODEL2MODEL = 1
   FLOES = 0

time0 = time.time()
###########################################################################

###########################################################################
# NOTE apparently daily files for the wavesice forecast lose information about maximum floe size.
# TP4archv file from 1200h will be used instead, the following "if" allow the user to analyze a TP4DAILY from the ice_only
if 0:
	# Read TP4 Daily
	print ''
	ncfil = ''.join( glob.glob('./data/TP4DAILY*.nc'))
	print('TP4DAILY ice_only file = ' +ncfil+'\n')
	slon		= 'longitude'
	slat		= 'latitude'
	sconc		= 'fice'
	lon		= Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
	lat		= Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
	conc		= Mrdg.nc_get_var(ncfil,sconc,time_index=0)
	X,Y      = lqm(lon[:,:],lat[:,:],inverse=False)
	Z        = conc[:,:].data
	mask     = conc[:,:].mask
	Z[mask]  = np.NaN
	
	# Date
	dadate		= ncfil[-11:-3]

else:
	# Read TP4arch_wav
	print ''
	ncfil = ''.join( glob.glob('./data/TP4archv*.nc'))
	print('TP4archv_wav file = ' +ncfil+'\n')
	slon = 'longitude'
	slat = 'latitude'
	sconc = 'fice'
	sdmax = 'dmax'
	lon = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
	lat = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
	conc = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
	dmax = Mrdg.nc_get_var(ncfil,sdmax,time_index=0)
	X,Y = lqm(lon[:,:],lat[:,:],inverse=False)
	if FLOES:
		Z = dmax[:,:].data
		mask = dmax[:,:].mask
	else:
		Z = conc[:,:].data
		mask = conc[:,:].mask
	Z[mask] = np.NaN

	# Date
	dadate		= ncfil[-19:-11]
# maybe we don't need this but i'm just going to leave it here...
#	globals()[dadate] = dadate
###########################################################################

###########################################################################
# Read in OSI_SAF file
ncfil2 = ''.join( glob.glob('./data/ice_conc_nh_polstere-100_multi_*.nc'))
print('OSISAF file = '+ncfil2+'\n')
clon     = 'lon'
clat     = 'lat'
cconc    = 'ice_conc'
lon2     = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
lat2     = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
conc		 = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
xc			 = Mrdg.nc_get_var(ncfil2,'xc')
yc			 = Mrdg.nc_get_var(ncfil2,'yc')
X2,Y2		 = lqm(lon2[:,:],lat2[:,:],inverse=False)
XO			 = np.copy(X2)
YO			 = np.copy(Y2)
Z2       = conc[:,:].data
mask2    = conc[:,:].mask
Z2[mask2] = np.NaN
ZO				= Z2/100

# getting ready for reprojection
X3 = X.reshape(X.size)
Y3 = Y.reshape(Y.size)
Z3 = Z.reshape(Z.size)
C = [X3,Y3]
C = np.array(C)
C = C.T

###########################################################################
# Reprojection
###########################################################################
# NOTE it's possible to study the areas between different thresholds of the same dataset (i.e. model)
if MODEL2MODEL:
	# this is a model2model analysis
	ZN = grd(C,Z3,(X2,Y2),method='nearest')
	if FLOES:
		# model2model_floe_size
		print 'Floe Size Distribution'
		print ''
		# Binary between model
		BN,BO,DN = binary_mod_3(ZN,300)
	else:
		# model2model_ice_conc
		# Binary between model
		BN,BO,DN = binary_mod_2(ZN,ZN,.15,.80)
else:
  # This is a model2osisaf analysis
	# Interpolation can be done with other methods ('linear','cubic'<--doesn't work for our data)
	ZN = grd(C,Z3,(X2,Y2),method='nearest')
	# BINARY
	BN,BO,DN = binary_mod(ZN,ZO,.15)
	# binary contour <-not used in this analysis
	#NDN = binary_cont(XO,YO,DN,BO,BN)
	elapsedtime = time.time() - time0
	print 'Reprojection done in ',elapsedtime
	print ''
###########################################################################

###########################################################################
# AARI charts section
if MODEL2MODEL and FLOES:
	apoly_list = []
	
	# Read in AARI charts for Barents Sea
	ncfil3 = ''.join( glob.glob('./ice_charts/AARI/*_bar_'+dadate+'.txt'))
	if ncfil3 == '':
		print('No Ice Chart for Barents sea in '+dadate+'\n')
	else:
		print('AARI Barents ice chart = '+ncfil3+'\n')
		region = 'Barents'
		icb = open(ncfil3)
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
		for i in ad:
			apoly = aari_poly(i,region,ad[i][:,1],ad[i][:,2],ad[i][:,3],X2,Y2,lqm)
			apoly_list.append(apoly)
	
	# Read in AARI charts for Greenland Sea
	ncfil4 = ''.join( glob.glob('./ice_charts/AARI/*_gre_'+dadate+'.txt'))
	if ncfil4 == '':
		print('No Ice Chart for Greenland sea in '+dadate+'\n')
	else:
		print('AARI Greenland ice chart = '+ncfil4+'\n')
		region = 'Greenland'
		icg = open(ncfil4)
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
		if ncfil3 == '':
			ad = {}
		for n in range(gpolyn):
			poly = []
			for m,em in enumerate(icg_cont):
				if em[0] == n+1:
					poly.append(em)
			if ncfil3 == '':
				ad["apoly"+str(n+1)] = np.array(poly)
			else:
				ad["apoly"+str(n+bpolyn+1)] = np.array(poly)
		for i in ad:
			apoly = aari_poly(i,region,ad[i][:,1],ad[i][:,2],ad[i][:,3],X2,Y2,lqm)
			apoly_list.append(apoly)
	
	# Read in AARI charts for both Barents and Greenland
	# sometimes the charts are united in case the some polygons belong to both regions
	if ncfil3 == '' and ncfil4 == '':
		ncfil5 = ''.join( glob.glob('./ice_charts/AARI/*'+dadate+'.txt'))
		if ncfil5 == '':
			print('No Common (Barents+Greenland) Ice Chart in '+dadate+'\n')
		else:
			print('AARI Common Barents and Greenland sea ice chart = '+ncfil5+'\n')
			region = 'Barents/Greenland'
			icbg = open(ncfil5)
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
			for i in ad:
				apoly = aari_poly(i,region,ad[i][:,1],ad[i][:,2],ad[i][:,3],X2,Y2,lqm)
				apoly_list.append(apoly)
	
	if apoly_list != []:
		ad = {}
		for x in range(len(apoly_list)):
			ad[apoly_list[x].polygon_name]=apoly_list[x]
		for key,value in sorted(ad.items()):
			globals()[key] = value
	
		elapsedtime = time.time() - time0
		print 'AARI Polygon identification done in ',elapsedtime
###########################################################################

###########################################################################
# print some stats for the "Difference" Dataset
get_stats(DN,dadate)
###########################################################################

###########################################################################
# Cutting out closed bays/seas
# USING THE BORDERS FROM arcGIS -> see baffin_bay function
if 0:
	# NOTE i,j are inverted for full maps (i.e. ZO,ZN,BO,BN,DN etc...)
	a = DN.shape
	
	# North Canada mess of islands
	for i in range(a[1]):
	 	for j in range(a[0]):
	 		 if i < 250 and 617 < j < 800 and (DN[j][i] == 1 or DN[j][i] == -1):
	 				DN[j][i] = 0
	 		 elif i < 300 and 650 < j < 800 and (DN[j][i] == 1 or DN[j][i] == -1):
	 				DN[j][i] = 0	 
	 		 elif i < 350 and 350 < j < 750 and (DN[j][i] == 1 or DN[j][i] == -1):
	 				DN[j][i] = 0	 
	
	# Hudson Bay & Northwestern Passages
	for i in range(a[1]):
	 	for j in range(a[0]):
	 		 if i < 270 and j > 760 and (DN[j][i] == 1 or DN[j][i] == -1):
	 				DN[j][i] = 0
	 		 elif i < 276 and j > 875 and (DN[j][i] == 1 or DN[j][i] == -1):
	 				DN[j][i] = 0

###########################################################################
# finding contours from the difference data map
pos = msr.find_contours(DN,.5)
pos = sorted(pos, key=len)
pos = pos[::-1]
neg = msr.find_contours(DN,-.5)
neg = sorted(neg, key=len)
neg = neg[::-1]
###########################################################################

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

	# save a poly for testing
	nt = 1
	npfil = outdir+'/poly'+str(nt)+'.npz'
	print('saving to '+npfil)
	#
	poly_test = poly_list[nt]
	xy_coords = poly_test.xy_list
	xy_coords2 = [tuple(xyc) for xyc in xy_coords]
	fvals2 = 1*poly_test.function_vals
	
	# save file
	np.savez(npfil,xy=xy_coords,func_vals=fvals2)

	# save a poly for testing
	nt = 2
	npfil = outdir+'/poly'+str(nt)+'.npz'
	print('saving to '+npfil)
	#
	poly_test = poly_list[nt]
	xy_coords = poly_test.xy_list
	xy_coords2 = [tuple(xyc) for xyc in xy_coords]
	fvals2 = 1*poly_test.function_vals
	
	# save file
	np.savez(npfil,xy=xy_coords,func_vals=fvals2)

	# save a poly for testing
	nt = 3
	npfil = outdir+'/poly'+str(nt)+'.npz'
	print('saving to '+npfil)
	#
	poly_test = poly_list[nt]
	xy_coords = poly_test.xy_list
	xy_coords2 = [tuple(xyc) for xyc in xy_coords]
	fvals2 = 1*poly_test.function_vals
	
	# save file
	np.savez(npfil,xy=xy_coords,func_vals=fvals2)
				
