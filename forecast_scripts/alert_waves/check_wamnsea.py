# Set up email alert for when WAMNSEA waves become very large
# - possibly limit search to near ice?

from netCDF4 import Dataset
import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure as msr

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import mod_reading as Mrdg

############################################################################
# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')
############################################################################

SEND_EMAIL  = 1
CHECK_NC    = 1
odir        = 'out' # where to put temporary outputs
if not os.path.exists(odir):
   # not in forecast_scripts directory, go to work_py
   odir  = '/work/timill/work_py/out' # where to put temporary outputs
   if not os.path.exists(odir):
      os.mkdir(odir)

#print 'Select the contour method:   '
#print 'Type 0 for basemap contour package'
#print 'Type 1 for Threshold method'
#chc = raw_input()
chc = '1'

############################################################################
def finish_map(bm):
   # finish_map(bm)
   # *bm is a basemap
   bm.drawcoastlines()
   bm.fillcontinents(color='gray')

   # draw parallels and meridians.
   bm.drawparallels(np.arange(60.,91.,10.),\
         labels=[True,False,True,True]) # labels = [left,right,top,bottom]
   bm.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
         labels=[True,False,False,True])
   bm.drawmapboundary() # fill_color='aqua')

   return
############################################################################
def dist_cont2cont(xice,yice,xwav,ywav,bm):
  n1 = len(xice)
  dist2 = np.zeros(n1)
  for n in range(n1):
    x1n = xice[n]
    y1n = yice[n]
    dist1 = np.sqrt(pow(xwav-x1n,2)+pow(ywav-y1n,2))
    dist2[n] = dist1.min()
  imin = dist2.argmin()
  dist = dist2[imin] # want xmin,ymin on ice corresponding to dist
  xmin = xice[imin]
  ymin = yice[imin]
  lon,lat = bm(xmin,ymin,inverse=True)
  return dist,lon,lat
############################################################################
def binary_cont(data,thresh):
   Z = np.copy(data)
   Z[Z<thresh] = 0
   Z[Z>=thresh] = 1
   return Z
############################################################################
def binary_diff(data,thresh1,thresh2):
	Z1 = np.copy(data)
	Z2 = np.copy(data)
	Z1[Z1<thresh1] = 0
	Z1[Z1>=thresh1] = 1
	Z2[Z2<thresh2] = 0
	Z2[Z2>=thresh2] = 1
	D = Z2 - Z1
	return D,Z1,Z2
############################################################################
# setup stereographic basemap
# - for plotting and also for limiting search area
# lat_ts is latitude of true scale.
# (lon_0,lat_0) is central point -> (x,y)=(0,0)
rad   = 18.          # approx radius of image (degrees)
xmax  = rad*111.e3   # half width of image [m]
ymax  = rad*111.e3   # half height of image [m]
cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
#
lat_ts   = 74. # deg N
lon_0    = 50. # deg E
lat_0    = 74. # deg N
#
bm = Basemap(width=2*xmax,height=2*ymax,\
             resolution=cres,projection='stere',\
             lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)

if CHECK_NC:
	# define netcdf file  
	import time
	from datetime import date, timedelta
	tday = date.today()
	yday = date.today() - timedelta(1)
	yday2 = date.today() - timedelta(2)
	cday = tday.strftime('%Y%m%d')
	cyear = tday.strftime('%Y')
	pday = yday.strftime('%Y%m%d')
	pday2 = yday2.strftime('%Y%m%d')
#	cday  = '20150707'
#	cyear = '2015'
#	pday  = '20150706'
#	pday2 = '20150705'
#	tday  = '2015-07-07'
#	yday2 = '2015-07-05'
	wmsc = '/work/shared/nersc/msc/WAMNSEA/'+ cyear + '/forecasts/'
	ncfil = wmsc + 'wam_nsea.fc.' + cday + '.nc' # should be determined from today's date use "fc"
#	ncfil = 'wam_nsea.fc.'+cday+'.nc'
	print('WAMNSEA file = ' + ncfil+'\n')
	
	# get info about nc file
	ncinfo_wav = Mrdg.nc_getinfo(ncfil)
	times = ncinfo_wav.timevalues # hours from 1st time in file
	Ntimes = len(times)
	print('Time values (h):')
	print(times)
	print(' ')
	
	# get lon/lat and restrict to relevant area
	slon = 'longitude'
	slat = 'latitude'
	sswh = 'significant_wave_height'
	fnday = 'max_time'
	lon = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
	lat = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
	X,Y = bm(lon[:,:],lat[:,:],inverse=False)
	in_area = np.logical_and(abs(X)<xmax,abs(Y)<ymax)
	
	# read in OSISAF conc file
	# plot conc + ice edge (15% bm.pcontour? )
	# get lon/lat and restrict to relevant area
	osisaf = '/work/shared/nersc/msc/OSI-SAF/' + str(yday.year) + '_nh_polstere'
	ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday + '1200.nc'
	if not os.path.exists(ncfil2):
		ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday2 + '1200.nc'
#	ncfil2 = 'ice_conc_nh_polstere-100_multi_201507061200.nc'
	print('OSISAF file = '+ncfil2+'\n')
	clon = 'lon'
	clat = 'lat'
	cconc = 'ice_conc'
	edge_level = 15
	lon2 = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
	lat2 = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
	X2,Y2 = bm(lon2[:,:],lat2[:,:],inverse=False)
	#in_area = np.logical_and(abs(X)<xmax,abs(Y)<ymax)
	conc = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
	Z2 = conc[:,:].data
	mask2 = conc[:,:].mask
	Z2[mask2] = np.NaN
	
	
	#check_list=0
	#check_list=range(Ntimes)
	check_list=[Ntimes-1]
	for loop_i in check_list:
		swh = Mrdg.nc_get_var(ncfil,sswh,time_index=loop_i)
		# swh[:,:] is a numpy masked array
		Z = swh[:,:].data
		mask = swh[:,:].mask
		Zmax = swh[in_area].max() # restrict max calc to relevant area
		Zmin = swh[in_area].min() # restrict max calc to relevant area
		
		# make plot
		Z[mask] = np.NaN
		# if mask is 'True' (data is missing or invalid) set to NaN

##############################################################################
		
		if 0:
			# make test plot showing ice edge contour with bm.contour
			g = plt.figure()
			bm.pcolor(X2,Y2,Z2,vmin=0,vmax=100)
			# plot ice edge
			# (get 15% conc contour)
			# http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines
			cs1   = bm.contour(X2,Y2,Z2,[edge_level])
			print(cs1)
			coll1 = cs1.collections
			nlev1 = len(coll1)
			for nl in range(nlev1):
			  p = coll1[nl].get_paths() # only one conc contour so use 0 
			  nseg = len(p)
			  for ns in range(nseg):
					# loop over segments
					v = p[ns].vertices
					x = v[:,0]
					y = v[:,1]
					bm.plot(x,y,'m',linewidth=1)
			
			finish_map(bm)
			figname  = odir+'/bm_contour.png'
			plt.savefig(figname)
			plt.close()
			g.clf()
		
			# make test plot showing ice edge contour with threshold 
			g = plt.figure()
			Z4 = np.copy(Z2)
			Z4[Z4<15]=0
			Z4[Z4>=15]=1
			thenans = np.isnan(Z4)
			Z4[thenans] = 0
			bm.contour(X2,Y2,Z4,1)
			finish_map(bm)
			figname  = odir+'/threshold_contour.png'
			plt.savefig(figname)
			plt.close()
			g.clf()
		
################################################################################
		
		if 0:
			# make test plot showing ice edge contour with skimage.measure
			# apparently this is still not working properly, probably because
			# skimage works with indexes and the basemap module wants lon/lat
			# TODO from indexes to lon/lat
			g = plt.figure()
			bm.pcolor(X2,Y2,Z2,vmin=0,vmax=100)
			# we first have to create a binary dataset for the 15%
			Z15 = binary_cont(Z2,15)
			# get contour
			cont = msr.find_contours(Z15,.5)
			for n,en in enumerate(cont):
				bm.plot(en[:,0],en[:,1],'c-',linewidth=2)

			finish_map(bm)
			figname  = odir+'/bm_contour.png'
			plt.savefig(figname)
			plt.close()
			g.clf()
		
################################################################################
		
		if 1:
			################################################################################
			# SENSIBLE POINTS WITH BM.CONTOUR
			################################################################################
			
			# Plotting waves
			f  = plt.figure()
			bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
			print('range in '+sswh+' (m):') 
			print(Zmin,Zmax)
			print(' ')
			
			cb = plt.colorbar()
			cb.set_label("Significant Wave Height [m]",rotation=270)
			
			##############################################################################
			# plot ice edge
			# (get 15% conc contour)
			##############################################################################
			
			cs1   = bm.contour(X2,Y2,Z2,[edge_level])
			coll1 = cs1.collections
			nlev1 = len(coll1)
			for nl in range(nlev1):
			  p     = coll1[nl].get_paths() # only one conc contour so use 0 
			  nseg  = len(p)
			  for ns in range(nseg):
					# loop over segments
					v     = p[ns].vertices
					x     = v[:,0]
					y     = v[:,1]
					bm.plot(x,y,'m',linewidth=1.5)
			
			##############################################################################
			
			# working on the waves threshold
			Hthresh=3
			Hmax=np.ceil(Zmax)
			Hlev=np.arange(Hthresh,Hmax,.5)
			
			if len(Hlev)==0:
				#no waves over threshhold
				out_list = []
			else:
				#some waves over threshhold
				cs0   = bm.contour(X,Y,Z,Hlev)
				coll0 = cs0.collections
				nlev0 = len(coll0)
				dist_thresh=50.e3
				out_list=[]
				for nl in range(nlev0):
				  swh0=cs0.levels[nl]
				  p     = coll0[nl].get_paths() # only one conc contour so use 0 
				  nseg  = len(p)
				  for ns in range(nseg):
						# loop over segments
						v = p[ns].vertices
						x = v[:,0]
						y = v[:,1]
						
						#loop over all ice edges
						for nl1 in range(nlev1):
						  p1 = coll1[nl1].get_paths() # only one conc contour so use 0 
						  nseg1 = len(p1)
						  for ns1 in range(nseg1):
								# loop over segments
								v1 = p1[ns1].vertices
								x1 = v1[:,0]
								y1 = v1[:,1]
								dist,lon_min,lat_min=dist_cont2cont(x1,y1,x,y,bm) # output (lon,lat) for ice
								if dist<dist_thresh:
								  dist_list=[dist,lon_min,lat_min,swh0]
								  out_list.append(dist_list)
			
			##############################################################################
			
			if 1:
				# add test point to plot
				# - to check if SAR image is ordered in the right place
				# - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
				#   then use trial and error
				nout = len(out_list)
				for mm in range(nout):
					list0=out_list[mm]
					if list0[3] >= 3 and list0[0] <= 5:
						lon_plot = list0[1]
						lat_plot = list0[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'og',markersize=5)
					elif list0[3] >= 4 and list0[0] <= 20:
						lon_plot = list0[1]
						lat_plot = list0[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'om',markersize=5)
					elif list0[3] >= 5 and list0[0] <= 50:
						lon_plot = list0[1]
						lat_plot = list0[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'ow',markersize=5)
			
			finish_map(bm)
			
			# date+time to title and file name
			# label '$H_s$, m' to colorbar
			wav = Dataset(ncfil)
			fnday = getattr(wav,'max_time')
			fnday = fnday.replace(" ", "_")
			fnday = fnday.replace(":", "")
			fnday = fnday.replace("-", "")
			figname  = odir+'/img/'+fnday+'.png'
                        if not os.path.exists(odir+'/img'):
                           os.mkdir(odir+'/img')
			plt.savefig(figname)
			print('saving figure:')
			print(figname+'\n')
			plt.close()
			f.clf()
			# bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function
			
			
			# swh to find when there are large waves (>4m) in the vicinity of the ice.
			# write and send an email to warn when this happens (so we can order some SAR images)
			
			# filename of text file to form contents of email message 
                        if not os.path.exists(odir+'/lst'):
                           os.mkdir(odir+'/lst')
			textfile = odir+'/lst/'+fnday+'_list.txt'
			nout=len(out_list)
                        SEND_EMAIL2=(nout>0)
			if SEND_EMAIL2:
			  tf=open(textfile,'w')
			
			  for mm in range(nout):
			    list0=out_list[mm]
			    line='\n' #get info from out_list
			    if list0[3] >= 3 and list0[0] <= 5:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list0[ii])+line
			    elif list0[3] >= 4 and list0[0] <= 20:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list0[ii])+line
			    elif list0[3] >= 5 and list0[0] <= 50:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list0[ii])+line
			
			    tf.write(line)
			
			  tf.close()
			
			
			##############################################################################
			# SENSIBLE POINTS WITH THRESHOLD
			##############################################################################
			
			# Plotting waves
			f  = plt.figure()
			bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
			print('range in '+sswh+' (m):') 
			print(Zmin,Zmax)
			print(' ')
			
			cb = plt.colorbar()
			cb.set_label('SWH [m]',rotation=270)
			
			##############################################################################
			# plot ice edge
			# (get 15% conc contour)
			##############################################################################
			
			# ice pack
			Z3 = np.copy(Z2)
			Zval = np.ma.array(Z3)
			Z3 = np.ma.masked_where(Zval < 15, Zval)
			bm.pcolor(X2,Y2,Z3,cmap='RdPu')
			
			# ice edge contour with threshold 
			Z4 = np.copy(Z2)
			Z4[Z4<15]=0
			Z4[Z4>=15]=1
			thenans = np.isnan(Z4)
			Z4[thenans] = 0
			cs4 = bm.contour(X2,Y2,Z4,1)
			coll4 = cs4.collections
			nlev4 = len(coll4)
			for nl4 in range(nlev4):
			  p4 = coll4[nl4].get_paths() # only one conc contour so use 0 
			  nseg4 = len(p4)
			  for ns4 in range(nseg4):
					v4 = p4[ns4].vertices
					x4 = v4[:,0]
					y4 = v4[:,1]
					bm.plot(x4,y4,'m',linewidth=1.5)
					
					# working on the waves threshold
					Hthresh=3
					Hmax=np.ceil(Zmax)
					Hlev=np.arange(Hthresh,Hmax,.5)
					
					if len(Hlev)==0:
						#no waves over threshhold
						out_list = []
					else:
						#some waves over threshhold
						cs3 = bm.contour(X,Y,Z,Hlev)
						coll3 = cs3.collections
						nlev3 = len(coll3)
						dist_thresh=50.e3
						out_list=[]
						for nl3 in range(nlev3):
						  swh3=cs3.levels[nl3]
						  p3 = coll3[nl3].get_paths() # only one conc contour so use 0 
						  nseg3  = len(p3)
						  for ns3 in range(nseg3):
								# loop over segments
								v3 = p3[ns3].vertices
								x3 = v[:,0]
								y3 = v[:,1]
								
								#loop over all ice edges
								for nl4 in range(nlev4):
								  p4 = coll4[nl4].get_paths() # only one conc contour so use 0 
								  nseg4 = len(p4)
								  for ns4 in range(nseg4):
										# loop over segments
										v4 = p4[ns4].vertices
										x4 = v4[:,0]
										y4 = v4[:,1]
										dist,lon_min,lat_min=dist_cont2cont(x4,y4,x3,y3,bm) # output (lon,lat) for ice
										if dist<dist_thresh:
										  dist_list=[dist,lon_min,lat_min,swh3]
										  out_list.append(dist_list)
			
##############################################################################
			
			if 1:
				# add test point to plot
				# - to check if SAR image is ordered in the right place
				# - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
				#   then use trial and error
				nout = len(out_list)
				for mm in range(nout):
					list1=out_list[mm]
					if list1[3] >= 3 and list1[0] <= 5:
						lon_plot = list1[1]
						lat_plot = list1[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'og',markersize=5)
					elif list1[3] >= 4 and list1[0] <= 20:
						lon_plot = list1[1]
						lat_plot = list1[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'om',markersize=5)
					elif list1[3] >= 5 and list1[0] <= 50:
						lon_plot = list1[1]
						lat_plot = list1[2]
						print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
						x_plot,y_plot  = bm(lon_plot,lat_plot)
						bm.plot(x_plot,y_plot,'ow',markersize=5)
			
			finish_map(bm)
			
			# date+time to title and file name
			# label '$H_s$, m' to colorbar
			wav = Dataset(ncfil)
			fnday = getattr(wav,'max_time')
			fnday = fnday.replace(" ", "_")
			fnday = fnday.replace(":", "")
			fnday = fnday.replace("-", "")
			figname  = odir+'/img/'+fnday+'_threshold.png'
			plt.savefig(figname)
			print('saving figure:')
			print(figname+'\n')
			plt.close()
			f.clf()
			# bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function
			
			
			# swh to find when there are large waves (>4m) in the vicinity of the ice.
			# write and send an email to warn when this happens (so we can order some SAR images)
			
			# filename of text file to form contents of email message 
			textfile = odir+'/lst/'+fnday+'_threshold_list.txt'
			nout=len(out_list)
                        SEND_EMAIL3 = (nout>0)
			if SEND_EMAIL3:
			  tf=open(textfile,'w')
			
			  for mm in range(nout):
			    list1=out_list[mm]
			    line='\n' #get info from out_list
			    if list1[3] >= 3 and list1[0] <= 5:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list1[ii])+line
			    elif list1[3] >= 4 and list1[0] <= 20:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list1[ii])+line
			    elif list1[3] >= 5 and list1[0] <= 50:
						for ii in range(3,-1,-1):
						   line=3*' '+str(list1[ii])+line
			
			    tf.write(line)
			
			  tf.close()

################################################################
# EMAIL SYSTEM
if SEND_EMAIL and (SEND_EMAIL2 or SEND_EMAIL3):
	import subprocess
	subprocess.call(["chmod +x out/img/*"])
	subprocess.call(["chmod +x out/lst/*"])
	subprocess.check_call(['/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/alert_waves/waves_alert.sh', fnday])

################################################################
