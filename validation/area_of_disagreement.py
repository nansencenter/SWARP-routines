# Python script for reprojecting model product into OSI-SAF grid

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
# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')
############################################################################

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
############################################################################
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
############################################################################
def figure_save(X,Y,Z,name,m):
	f = plt.figure()
	m.pcolor(X,Y,Z,cmap='winter',vmin=-1,vmax=1)
	finish_map(m)
	fname = name+'.png'
	plt.colorbar()
	plt.title(name)
	plt.savefig(fname,format='png',dpi=1000)
	plt.close()
	f.clf()
	return
###########################################################################
###########################################################################
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
	#f = open('./outputs/aod/'+name+'.txt','w')
	#f.write('Number of elements:	%s \n' % stat)
	#f.write('Hit:			%s \n' % stat0)
	#f.write('Underprediction(+1):	%s - %s \n' % (stat1,pone))
	#f.write('Overprediction(-1):	%s - %s \n' % (stat2,pmone))
	#f.write(' ')
	#f.close
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
	def __init__(self,cont,OSI,MOD,DIFF,X,Y,number,polygon_status=1):
		self.polygon_status		= polygon_status
		self.polygon_name			= 'polygon'+str(number)
		# class definition
		if len(cont) <= 30:
			self.polygon_class = 'S'
		elif len(cont) > 30 and len(cont) <= 100:
			self.polygon_class = 'M'
		elif len(cont) > 100 and len(cont) <= 200:
			self.polygon_class = 'B'
		elif len(cont) > 200:
			self.polygon_class = 'H'
		self.ij_list					= cont
		self._diff_with_terrain_(DIFF,OSI,MOD) 
		self._ij2xy_(cont,X,Y)
		self._split_cont(OSI,MOD)
		return
	
	def _ij2xy_(self,cont,X,Y):
		# changes indexes to x and y (NOTE i and j are inverted -> i = [:,1], j = [:,0])
		x			= []
		y			= []
		xvec	=	range(len(cont[:,0]))
		for n,en in enumerate(xvec):
		        en	=	X[cont[n,0]][cont[n,1]]
		        x.append(en)
		yvec	=	range(len(cont[:,1]))
		for n,en in enumerate(yvec):
		        en	=	Y[cont[n,0]][cont[n,1]]
		        y.append(en)
		xy_list = zip(x,y)
		xy_list = np.array(xy_list)
		self.xy_list = xy_list
		return
	
	def _diff_with_terrain_(self,DIFF,OSI,MOD):
		# getting back the terrain lost during binary_mod (for contour finding reasons)
		DN			= DIFF
		ZO			= OSI
		ZM			= MOD
		thenan	= np.isnan(ZO)
		thenan2	= np.isnan(ZM)
		DN[thenan]	= None
		DN[thenan2]	= None
		self.difference = DN
		return
	
	def _split_cont(self,OSI,MOD):
		vs = self.ij_list # list of (i,j) pixel indices
		npoly = self.polygon_name
		mdl_cont = []
		osi_cont = []
		ukn_cont = []
		func_vals= []
		func_mod = 0	# value of func_vals at model ice edge
		func_osi = 1	# value of func_vals at OSISAF ice edge
		func_unk = 2	# value of func_vals at other parts of contour
		
		if self.polygon_status==1:
			# if polygon is negative - overestimation of ice
			self.model_status='Model OVERestimate'
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
						func_val=func_mod
						check_cont = 1
					elif OSI[h][v] == MOD[h][v] == 0:
						osi_cont.append(el)
						func_val=func_osi
						check_cont = 1
				if check_cont == 0:
					ukn_cont.append(el)
					func_val=func_unk
				func_vals.append(func_val)
		else:
			# if polygon is positive - underestimation of ice
			self.model_status='Model UNDERestimate'
			for n,el in enumerate(vs):
				around = ((el[0],el[1]+1),(el[0],el[1]-1),\
				          (el[0]+1,el[1]),(el[0]-1,el[1]),\
				          (el[0]+1,el[1]+1),(el[0]-1,el[1]+1),\
				          (el[0]+1,el[1]-1),(el[0]-1,el[1]-1))
				check_cont = 0
				for h,v in around:
					if OSI[h][v] == MOD[h][v] == 0:
						mdl_cont.append(el)
						func_val=func_mod
						check_cont = 1
					elif OSI[h][v] == MOD[h][v] == 1:
						osi_cont.append(el)	
						func_val=func_osi
						check_cont = 1
				if check_cont == 0:
					ukn_cont.append(el)
					func_val=func_unk
				func_vals.append(func_val)
		mdl_cont = np.array(mdl_cont)
		osi_cont = np.array(osi_cont)
		ukn_cont = np.array(ukn_cont)
		func_vals	= np.array(func_vals)
		# include sorted contours as port of object
		self.model_ice_edge 	= mdl_cont
		self.osisaf_ice_edge = osi_cont
		self.unknown_edge = ukn_cont
		self.function_vals = func_vals
		return

############################################################################
class aod_stats:
	def __init__(self,aod,basemap=None):
		self.ij_list = aod.ij_list
		self.xy_list = aod.xy_list
		self.number_of_points = len(aod.ij_list)
		self.polygon_status = aod.polygon_status
		self.name = aod.polygon_name
		self.pclass = aod.polygon_class
		self.patch_sign = aod.model_status
		self.mcont = aod.model_ice_edge
		self.ocont = aod.osisaf_ice_edge
		self.ucont = aod.unknown_edge
		self.diff = aod.difference
		
		#get euclidean area
		#TODO get area on sphere?
		self.poly_area()
		
		#get euclidean perimeter
		#(geometry_planar.curve_info could be used)
		self.poly_perimeter()
		
		# distances between contours set self.widths
		self.dist_edges()

		if basemap is not None:
			xy_list = self.xy_list
			self.list_lon,self.list_lat=basemap(xy_list[:,0],xy_list[:,1],inverse=True)
			
			# Perimeter mean - used to check diff with centroid 
			xmean	=	np.mean(xy_list[:,0])
			ymean	=	np.mean(xy_list[:,1])
			self.mlon,self.mlat	=	basemap(xmean,ymean,inverse=True)
			
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
			self.mclon,self.mclat = basemap(CX,CY,inverse=True)
			# defining the region
			mclon = self.mclon
			if mclon < 80 and mclon > 8:
				self.polygon_region = 'Barents_Kara_sea'
			elif mclon < 8 and mclon > -44:
				self.polygon_region = 'Greenland_sea'
			elif mclon < -44 and mclon > -70:
				self.polygon_region = 'Labrador_sea'
			else:
				self.polygon_region = 'Misc.'

			# Calculating lon/lat for model,osisaf,unknown contours and distances
			mcont		= self.mcont
			ocont		= self.ocont
			ucont		= self.ucont
			dist		= self.widths
			if len(mcont) != 0:
				self.mcontlon,self.mcontlat = basemap(mcont[:,0],mcont[:,1],inverse=True)
			if len(ocont) != 0:
				self.ocontlon,self.ocontlat = basemap(ocont[:,0],ocont[:,1],inverse=True)
			if len(ucont) != 0:
				self.ucontlon,self.ucontlat = basemap(ucont[:,0],ucont[:,1],inverse=True)
			distlonlat	= []
			if dist != [0,0,0,0,0]:
				for n,en in enumerate(dist):
					x1,y1 = basemap(en[1],en[2],inverse=True)
					x2,y2 = basemap(en[3],en[4],inverse=True)
					distlonlat.append([en[0],x1,y1,x2,y2])
				self.widthlonlat = distlonlat
		return

	def poly_area(self):
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
		self.area	=	a*100
		return

	def poly_perimeter(self):
		# Calculating perimeter in xy coordinates (unit = 10km)
		xy		= self.ij_list
		perim	= 0
		for n in range(len(xy)-1):
			perim += np.sqrt(pow(xy[n+1][0]-xy[n][0],2)+pow(xy[n+1][1]-xy[n][1],2))
		self.perimeter	= perim*10
		return

	def dist_edges(self):
		# Calculating the min distance from a model point to any osisaf point
		mcont = self.mcont
		ocont = self.ocont
		ucont = self.ucont
		tcont = self.ij_list
		UKW   = 'Unknown contours < 20%'
		DMW   = 'Contour difference < 40%'
		width = []
		osi_width = []
		if len(mcont) == 0:
			DMW = 'Only OSISAF data, MODEL not present'
			if MODEL2MODEL:
				DMW = 'Only 80% edge, 15% not present'
			self.widths	= [0,0,0,0,0]
			self.osi_widths	= [0,0,0,0,0]
		elif len(ocont) == 0:
			DMW = 'Only MODEL data, OSISAF not present'
			if MODEL2MODEL:
				DMW = 'Only 15% edge, 80% not present'
			self.widths	= [0,0,0,0,0]
			self.osi_widths	= [0,0,0,0,0]
		else:
			ukn = 100*(len(ucont)/float(len(tcont)))
			dmo = 100*(abs(len(mcont)-len(ocont))/float(len(tcont)))
			if ukn >= 20:
				UKW	= 'WARNING - unknown contours > 20%'
			if dmo >= 40:
				DMW = 'WARNING - contours difference > 40%'
			for n,en in enumerate(mcont):
				dist_pt = []
				for m,em in enumerate(ocont):
					dist1	= np.sqrt(pow(en[0]-em[0],2)+pow(en[1]-em[1],2))
					dist_pt.append([dist1,en[0],en[1],em[0],em[1]])
				idx,value = min(enumerate(dist_pt),key=itg(1))
				width.append(value)
			for n,en in enumerate(ocont):
				dist_pt = []
				for m,em in enumerate(mcont):
					dist1	= np.sqrt(pow(en[0]-em[0],2)+pow(en[1]-em[1],2))
					dist_pt.append([dist1,en[0],en[1],em[0],em[1]])
				idx,value = min(enumerate(dist_pt),key=itg(1))
				osi_width.append(value)
			self.widths = width
			self.osi_widths = osi_width
		self.UKW		= UKW
 		self.DMW		= DMW
 		return

	def stat_chart(self,save=False):
		pname = self.name
		pclass = self.pclass
		region = self.polygon_region
		print ''
		print 'Statistic Chart for '+str(pname)
		print ''
		DN = self.diff
		DMW = self.DMW
		UKW = self.UKW
		ij = self.ij_list
		mcont = self.mcont
		ocont = self.ocont
		ucont = self.ucont
		dist = self.widths
		osi_dist = self.osi_widths
		clon = '%1.2f' %self.mclon
		clat = '%1.2f' %self.mclat
		clonlat = '{0}/{1}'.format(clon,clat)
		area = self.area
		area = '%1.4e' %area 
		perim = self.perimeter
		perim = '%1.4e' %perim
		if dist != [0,0,0,0,0]:
			mdist = np.array(self.widths)
			# changing units from decakilometer to kilometer
			mdist = mdist[:,0]*10
			mdist = np.median(mdist)
			mdist = '%1.2f' %mdist
		else:
			mdist = 'NaN'
		if osi_dist != [0,0,0,0,0]:
			odist = np.array(self.osi_widths)
			# changing units from decakilometer to kilometer
			odist = odist[:,0]*10
			odist = np.median(odist)
			odist = '%1.2f' %odist
		else:
			odist = 'NaN'
		pstat = self.polygon_status
		if pstat == 1:
			status = 'Underestimation'
		else:
			status = 'Overestimation'

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
		main.imshow(DN[::-1],cmap='winter')
		x1,x2,y1,y2 = np.min(ij[:,1])-10,np.max(ij[:,1])+10,np.min(ij[:,0])-10,np.max(ij[:,0])+10
		main.axvspan(x1,x2,ymin=1-((y1-320)/float(len(DN)-320)),\
		      ymax=1-((y2-320)/float(len(DN)-320)),color='red',alpha=0.3)
		main.axis([0,760,0,800])

		# Polygon image on the top right
		polyf.imshow(DN,cmap='winter')
		polyf.axis([x1,x2,y2,y1])
		if len(mcont) != 0:
			polyf.plot(mcont[:,1],mcont[:,0],'ro',markersize=4)
		if len(ocont) != 0:
			polyf.plot(ocont[:,1],ocont[:,0],'yo',markersize=4)
		if len(ucont) != 0:
			polyf.plot(ucont[:,1],ucont[:,0],'go',markersize=4)
		if dist != [0,0,0,0,0]:
			for n,en in enumerate(dist):
				polyf.plot([en[2],en[4]],[en[1],en[3]],color='black',alpha=0.1)
		if osi_dist != [0,0,0,0,0]:
			for n,en in enumerate(osi_dist):
				polyf.plot([en[2],en[4]],[en[1],en[3]],color='magenta',alpha=0.1)

		if MODEL2MODEL:
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
			      '4) Mean 15-80 Width = '+str(mdist)+' km\n'+ \
			      '5) Mean 80-15 Width = '+str(odist)+' km\n'+ \
			      '6) '+str(DMW)+'\n'+ \
			      '7) '+str(UKW)
			tab.text(.2,.2,txt,fontsize=15,bbox=dict(boxstyle='round',facecolor='white',alpha=1))
			if save:
				valid_class = 'Model2Model'
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
			       '4) Mean M-O Width = '+str(mdist)+' km\n'+ \
			       '5) Mean O-M Width = '+str(odist)+' km\n'+ \
			       '6) Status = '+str(status)+'\n'+ \
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

############################################################################
###########################################################################

# DEFINING THE BASEMAP
# low quality map
m = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),\
      projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
# better quality map
hqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
      projection='stere',lat_ts=70,lat_0=90,lon_0=-45)

# Getting a nice blank space before info's printing
print ''
FIGURE = raw_input('[1] to plot&save all figures, [Enter] to plot without save	') 
print ''
MODEL2MODEL = raw_input('[1] for model2model, [ENTER] to model2osisaf	') 
print ''

time0 = time.time()

# READ TP4 DAILY
ncfil = ''.join( glob.glob('./data/TP4DAILY*.nc'))
print('TP4DAILY ice_only file = ' +ncfil+'\n')
slon		= 'longitude'
slat		= 'latitude'
sconc		= 'fice'
lon		= Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
lat		= Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
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

# getting ready for reprojection
X3 = X.reshape(X.size)
Y3 = Y.reshape(Y.size)
Z3 = Z.reshape(Z.size)
C = [X3,Y3]
C = np.array(C)
C = C.T


# NOTE it's possible to study the areas between different thresholds of the same dataset (i.e. model)
if MODEL2MODEL:
	ZN = grd(C,Z3,(X2,Y2),method='nearest')
	# Binary between model
	BN,BO,DN = binary_mod_2(ZN,ZN,.15,.80)
else:
	# INTERPOLATION CAN BE DONE WITH OTHER METHODS ('linear','cubic'<--doesn't work for our data)
	ZN = grd(C,Z3,(X2,Y2),method='nearest')
	# BINARY
	BN,BO,DN = binary_mod(ZN,ZO,.15)
	# binary contour <-not used in this analysis
	#NDN = binary_cont(XO,YO,DN,BO,BN)
	elapsedtime = time.time() - time0
	print 'Reprojection done in ',elapsedtime
	print ''
get_stats(DN,dadate)
		
# finding contours from the difference data map
pos = msr.find_contours(DN,.9)
neg = msr.find_contours(DN,-.9)
poly_list=[]
for n,el in enumerate(pos):
	# classification of positive polygons
	aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=1)
	poly_list.append(aod)
try:
	n
except NameError:
	n = -1
for n2,el in enumerate(neg):
	# classification of negative (n+1 for good enumeration)
	n += 1
	aod=aod_poly(el,BO,BN,DN,XO,YO,n,polygon_status=0)
	poly_list.append(aod)

# changing name for every polygon, easier to work with
d = {}
for x in range(len(poly_list)):
	d['poly{0}'.format(x)]=poly_list[x]
for key,value in sorted(d.items()):
	globals()[key] = value

elapsedtime = time.time() - time0
print 'Polygon identification done in ',elapsedtime
print ''

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
		poly = aod_stats(el,m)
		if poly.polygon_region == 'Barents_Kara_sea':
			sbk_sea += 1
		elif poly.polygon_region == 'Greenland_sea':
			sgr_sea += 1
		elif poly.polygon_region == 'Labrador_sea':
			slb_sea += 1
		elif poly.polygon_region == 'Misc.':
			smiscel += 1
	elif el.polygon_class == 'M':
		mcount += 1
		poly = aod_stats(el,m)
		if poly.polygon_region == 'Barents_Kara_sea':
			mbk_sea += 1
		elif poly.polygon_region == 'Greenland_sea':
			mgr_sea += 1
		elif poly.polygon_region == 'Labrador_sea':
			mlb_sea += 1
		elif poly.polygon_region == 'Misc.':
			mmiscel += 1
	elif el.polygon_class == 'B':
		bcount += 1
		poly = aod_stats(el,m)
		if poly.polygon_region == 'Barents_Kara_sea':
			bbk_sea += 1
		elif poly.polygon_region == 'Greenland_sea':
			bgr_sea += 1
		elif poly.polygon_region == 'Labrador_sea':
			blb_sea += 1
		elif poly.polygon_region == 'Misc.':
			bmiscel += 1
	elif el.polygon_class == 'H':
		hcount += 1
		poly = aod_stats(el,m)
		if poly.polygon_region == 'Barents_Kara_sea':
			hbk_sea += 1
		elif poly.polygon_region == 'Greenland_sea':
			hgr_sea += 1
		elif poly.polygon_region == 'Labrador_sea':
			hlb_sea += 1
		elif poly.polygon_region == 'Misc.':
			hmiscel += 1
	poly_stat_list.append(poly)

# changing name to every polygon_stat, easier to work with
d = {}
for x in range(len(poly_stat_list)):
	d['poly{0}stats'.format(x)]=poly_stat_list[x]
for key,value in sorted(d.items()):
	globals()[key] = value

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
	if os.path.exists('Model2Model'):
		shutil.move('Model2Model',fin_dir)
	elapsedtime1 = time.time() - time1
	elapsedtime = time.time() - time0
	print 'Figures saved and moved in ',elapsedtime1
	print ''
	print 'Time needed for the job ',elapsedtime

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
			
