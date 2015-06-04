# Python script for reprojecting model product into OSI-SAF grid

from netCDF4 import Dataset
import sys,os
import glob
import numpy as np
import subprocess
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from skimage import measure as msr
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

sys.path.append('../py_funs')
import mod_reading as Mrdg

############################################################################
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
############################################################################
def binary_mod(data,thresh):
	ndata									= np.copy(data)
	ndata[ndata<thresh]		= 0
	ndata[ndata>=thresh]	= 1
	#thenans								= np.isnan(ndata)
	#ndata[thenans]				= 0
	return(ndata)
############################################################################
############################################################################
def binary_diff(data1,data2,thresh):
	mdata									= np.copy(data1)
	mdata[mdata<thresh]		= 0
	mdata[mdata>=thresh]	= 1
	odata									= np.copy(data2)
	odata[odata<thresh]		= 0
	odata[odata>=thresh]	= 1	
	ddata									= odata - mdata
	thenans								= np.isnan(ddata)
	ddata[thenans]				= 0
	return(ddata)
###########################################################################
###########################################################################
def figure_save(X,Y,Z,name,m):
	f = plt.figure()
	m.pcolor(X,Y,Z,cmap='YlGnBu',vmin=-1,vmax=1)
	finish_map(m)
	fname = './outputs/aod/'+name+'.png'
	plt.colorbar()
	plt.title(name)
	plt.savefig(fname,format='png',dpi=1000)
	plt.close()
	f.clf()
###########################################################################
###########################################################################
def get_stats(data,name):
	stat0	= (data == 0).sum()
	stat1	= (data == 1).sum()
	stat2	= (data == -1).sum()
	stat	= data.size
	pone	= (stat1 / float(stat)) * 100
	pmone	= (stat2 / float(stat)) * 100
	print "Number of elements:	",stat
	print "Hit:			",stat0
	print "Underprediction(+1):	",stat1,pone
	print "Overprediction(-1):	",stat2,pmone
	print " " 
	f = open('./outputs/aod/'+name+'.txt','w')
	f.write("Number of elements:	%s \n" % stat)
	f.write("Hit:			%s \n" % stat0)
	f.write("Underprediction(+1):	%s - %s \n" % (stat1,pone))
	f.write("Overprediction(-1):	%s - %s \n" % (stat2,pmone))
	f.write(" ")
	f.close
	return()
###########################################################################
###########################################################################
# FINDS THE CONT OF THE MDL(2) AND THE CONT OF OSI(-2)
def binary_cont(X,Y,D,O,M):
	ND = np.copy(D)
	for m,el in enumerate(X2):
		for n,el in enumerate(X2[m]):
			if D[m][n] == 1:
				around = ((m,n+1),(m,n-1),(m+1,n),(m-1,n),(m+1,n+1),(m-1,n+1),(m+1,n-1),(m-1,n-1))
				for hor,ver in around:
					if O[hor][ver] == M[hor][ver] == 1:
						ND[hor][ver] = 2
					elif O[hor][ver] == M[hor][ver] == 0:
						ND[hor][ver] = -2
					elif np.isnan(O[hor][ver]) or np.isnan(M[hor][ver]):
						ND[hor][ver] = 1
			elif D[m][n] == -1:
				around = ((m,n+1),(m,n-1),(m+1,n),(m-1,n),(m+1,n+1),(m-1,n+1),(m+1,n-1),(m-1,n-1))
				for hor,ver in around:
					if O[hor][ver] == M[hor][ver] == 1:
						ND[hor][ver] = -2
					elif O[hor][ver] == M[hor][ver] == 0:
						ND[hor][ver] = 2
					elif np.isnan(O[hor][ver]) or np.isnan(M[hor][ver]):
						ND[hor][ver] = -1
	return(ND)
############################################################################
# GETTING THE POLYGONS
class aod_poly:
	# gets single contour and sorts it
	def __init__(self,cont,OSI,MOD,X,Y,polygon_status=1):
		self.polygon_status		= polygon_status
		self.ij_list					= cont
		self._ij2xy_(cont,X,Y)
		self._split_cont(OSI,MOD)
		return
	
	def _ij2xy_(self,cont,X,Y):
		x			= []
		y			= []
		xvec	=	range(len(cont[:,0]))
		for n,en in enumerate(xvec):
			en	=	X[cont[n,0]][cont[n,1]]
			y.append(en)
		yvec	=	range(len(cont[:,1]))
		for n,en in enumerate(yvec):
			en	=	Y[cont[n,0]][cont[n,1]]
			x.append(en)
		xy_list = zip(x,y)
		xy_list = np.array(xy_list)
		self.xy_list = xy_list
		return

	def _split_cont(self,OSI,MOD):
		vs=self.ij_list # list of (i,j) pixel indices
		mdl_cont = []
		osi_cont = []
		ukn_cont = []
		func_vals= []
		func_mod = 0	# value of func_vals at model ice edge
		func_osi = 1	# value of func_vals at OSISAF ice edge
		func_unk = 2	# value of func_vals at other parts of contour
	
		if self.polygon_status==0:
			# if polygon is negative - overestimation of ice
			self.model_status="Model OVERestimate"

			for n,el in enumerate(vs):
				#getting all the neighbours
				around = ((el[0],el[1]+1),(el[0],el[1]-1),\
									(el[0]+1,el[1]),(el[0]-1,el[1]),\
									(el[0]+1,el[1]+1),(el[0]-1,el[1]+1),\
									(el[0]+1,el[1]-1),(el[0]-1,el[1]-1))
				check_cont = 0
				for h,v in around:
					if OSI[h][v] == MOD[h][v] == 1:
						mdl_cont.append(el)
						func_vals.append(func_mod)
						check_cont = 1
					elif OSI[h][v] == MOD[h][v] == 0:
						osi_cont.append(el)
						func_vals.append(func_osi)
						check_cont = 1
				if check_cont == 0:
					ukn_cont.append(el)
				  func_vals.append(func_unk)
		else:
			# if polygon is positive - underestimation of ice
			self.model_status="Model UNDERestimate"
			for n,el in enumerate(vs):
				around = ((el[0],el[1]+1),(el[0],el[1]-1),\
									(el[0]+1,el[1]),(el[0]-1,el[1]),\
									(el[0]+1,el[1]+1),(el[0]-1,el[1]+1),\
									(el[0]+1,el[1]-1),(el[0]-1,el[1]-1))
				check_cont = 0
				for h,v in around:
					if OSI[h][v] == MOD[h][v] == 0:
						mdl_cont.append(el)
						func_vals.append(func_mod)
						check_cont = 1
					elif OSI[h][v] == MOD[h][v] == 1:
						osi_cont.append(el)	
						func_vals.append(func_osi)
						check_cont = 1
				if check_cont == 0:
					ukn_cont.append(el)
					func_vals.append(func_unk)
		mdl_cont = np.array(mdl_cont)
		osi_cont = np.array(osi_cont)
		ukn_cont = np.array(ukn_cont)
		func_vals = np.array(func_vals)
		# include sorted contours as port of object
		self.model_ice_edge		=mdl_cont
		self.osisaf_ice_edge	=osi_cont
		self.unknown_edge=ukn_cont
		self.func_vals=func_vals
		return

############################################################################
class aod_stats:
	def __init__(self,aod,basemap=None):
		self.ij_list					= aod.ij_list
		self.xy_list					= aod.xy_list
		self.number_of_points	=	len(aod.ij_list)
		self.polygon_status		= aod.polygon_status
		self.mcont						= aod.model_ice_edge
		self.ocont						= aod.osisaf_ice_edge
		self.ucont						= aod.unknown_edge

		#get euclidean area
		#TODO get area on sphere?
		self.poly_area()

		#get euclidean perimeter
		#(geometry_planar.curve_info could be used)
		#TODO self.poly_perimeter() sets self.perimeter

		if basemap is not None:
			xy_list = self.xy_list
			self.list_lon,self.list_lat=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
			xmean	=	np.mean(xy_list[:,0])
			ymean	=	np.mean(xy_list[:,1])
			self.mlon,self.mlat	=	basemap(xmean,ymean,inverse=True)
			
			DX = 0
			DY = 0
			for n in range(len(xy_list)-1):
				DX += ((xy_list[n][0]+xy_list[n+1][0])-((xy_list[n][0]*xy_list[n+1])-(xy_list[n+1][0]*xy_list[n][1])))  
				DY += ((xy_list[n][1]+xy_list[n+1][1])-((xy_list[n][0]*xy_list[n+1])-(xy_list[n+1][0]*xy_list[n][1])))  
			CX = (1/float(6*self.area))*DX
			CY = (1/float(6*self.area))*DY
			self.mclon,self.mclat = basemap(CX,CY,inverse=True)

		# distances between contours set self.widths
		self.dist_edges()
		return

	def poly_area(self):
		vs=self.ij_list
		a = 0
		x0,y0 = vs[0]
		for [x1,y1] in vs[1:]:
			dx = x1-x0
			dy = y1-y0
			a += 0.5*(y0*dx - x0*dy)
			x0 = x1
			y0 = y1
		
		self.area	=	a
		return

	def dist_edges(self):
		mcont = self.mcont
		ocont = self.ocont
		ucont = self.ucont
		tcont = self.ij_list
		width = []
		ukn		= 100*(len(ucont)/float(len(tcont)))
		dmo		= 100*(abs(len(mcont)-len(ocont))/float(len(tcont)))
		if ukn >= 20:
			print "WARNING - unknown contours over 20%"
			UKW	= 1
		if dmo >= 20:
			print "WARNING - difference between Model Cont. and Osi Cont. over 20%"
			DMW	= 1
		for n,en in enumerate(mcont):
			dist_pt = []
			for m,em in enumerate(ocont):
				dist1	= np.sqrt(pow(en[0]-em[0],2)+pow(en[1]-em[1],2))
				dist_pt.append([dist1,en[0],en[1],em[0],em[1]])
			idx,value = min(enumerate(dist_pt),key=itg(1))
			width.append(value)

		self.widths = width
		return

#	def poly_fig(self):
#		fig,ax = plt.subplots()
#		ax.imshow(DN)
#		axins.imshow(
			

############################################################################
###########################################################################
#def find_poly(D,O,M):
#	ND		= np.copy(D)
#	pos		= msr.find_contours(D,.9)
#	neg		= msr.find_contours(D,-.9)
#	poly	=[]
##	for n,el in enumerate(neg):
##		if len(neg[n]) > 5:
##			poly.append(neg[n])
##	for n,el in enumerate(pos):
##		if len(pos[n]) > 5:
##			poly.append(pos[n])	
#	pos				=	np.array(pos)
#	neg				= np.array(neg)
#	mdl_cont	= []
#	osi_cont	= []
## now we need to distinguish the 2 borders
#	for num,pol in enumerate(pos):
#		for num2,pts in enumerate(pol):
#			pos_mdl	= []
#			pos_osi	= []
#			around = ((pts[0],pts[1]+1),(pts[0],pts[1]-1),(pts[0]+1,pts[1]),(pts[0]-1,pts[1]),(pts[0]+1,pts[1]+1),(pts[0]-1,pts[1]+1),(pts[0]+1,pts[1]-1),(pts[0]-1,pts[1]-1))
#			for hor,ver in around:
#				if O[hor][ver] == M[hor][ver] == 1:
#					pos_mdl.append(pts)	
#				elif O[hor][ver] == M[hor][ver] == 0:
#					pos_osi.append(pts)
#	for num,pol in enumerate(neg):
#		for num2,pts in enumerate(pol):
#			neg_mdl	= []
#			neg_osi	= []
#			around = ((pts[0],pts[1]+1),(pts[0],pts[1]-1),(pts[0]+1,pts[1]),(pts[0]-1,pts[1]),(pts[0]+1,pts[1]+1),(pts[0]-1,pts[1]+1),(pts[0]+1,pts[1]-1),(pts[0]-1,pts[1]-1))
#			for hor,ver in around:
#				if O[hor][ver] == M[hor][ver] == 1:
#					neg_mdl.append(pts)	
#				elif O[hor][ver] == M[hor][ver] == 0:
#					neg_osi.append(pts)
#	return(mdl_cont,osi_cont)
#############################################################################
###########################################################################


# DEFINING THE BASEMAP
m = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)

# CHOOSE
FIGURE = raw_input('1 for figure, [Enter] for not:	')

# READ TP4 DAILY
ncfil = ''.join( glob.glob('./data/TP4DAILY*.nc'))
print('TP4DAILY ice_only file = ' +ncfil+'\n')
slon		= 'longitude'
slat		= 'latitude'
sconc		= 'fice'
lon			= Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
lat			= Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
conc		= Mrdg.nc_get_var(ncfil,sconc,time_index=0)
X,Y      = m(lon[:,:],lat[:,:],inverse=False)
Z        = conc[:,:].data
mask     = conc[:,:].mask
Z[mask]  = np.NaN

# DATE
dadate		= ncfil[-11:-3]

# READ IN OSI-SAF FILE
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
X2,Y2		 = m(lon2[:,:],lat2[:,:],inverse=False)
XO			 = np.copy(X2)
YO			 = np.copy(Y2)
Z2       = conc[:,:].data
mask2    = conc[:,:].mask
Z2[mask2] = np.NaN
ZO				= Z2/100

# REPROJECTION
X3 = X.reshape(X.size)
Y3 = Y.reshape(Y.size)
Z3 = Z.reshape(Z.size)
C = [X3,Y3]
C = np.array(C)
C = C.T

# INTERPOLATION CAN BE DONE WITH ANOTHER METHOD ('linear')
ZN = grd(C,Z3,(X2,Y2),method='nearest')

# BINARY
BO	= binary_mod(ZO,.15)
BN	= binary_mod(ZN,.15)

# DIFFERENCES
DN = binary_diff(ZN,ZO,.15)

# CALL GENERAL STATS
#print "DN stats"
#get_stats(DN,'DNstats'+dadate)
#print "DL stats"
#get_stats(DL,'DLstats'+dadate)

NDN = binary_cont(XO,YO,DN,BO,BN)


		
pos = msr.find_contours(DN,.9)
neg = msr.find_contours(DN,-.9)
aod_polys=[]
for n,el in enumerate(pos):
	aod=aod_poly(el,BO,BN,XO,YO,polygon_status=1)
	aod_polys.append(aod)
for n,el in enumerate(neg):
	aod=aod_poly(el,BO,BN,XO,YO,polygon_status=0)
	aod_polys.append(aod)

if FIGURE:
	figure_save(X2,Y2,DN,'DN'+dadate,hqm)
	figure_save(X2,Y2,DL,'DL'+dadate,hqm)

# PLOT RESULTS
#plt.figure(0) 
#plt.imshow(ZO)
#plt.title('OSI')
#
#plt.figure(1)
#plt.imshow(ZN)
#plt.title('NRS')
#
#plt.figure(2)
#plt.imshow(ZL)
#plt.title('LIN')
#
#plt.figure(3)
#plt.imshow(ZC)
#plt.title('CUB')
#
#plt.show()


