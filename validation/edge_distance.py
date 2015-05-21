# Calculation of the minimum distance between model product and OSI-SAF daily average

from netCDF4 import Dataset
import sys,os
import glob
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure as ms

sys.path.append('../py_funs')
import mod_reading as Mrdg

# Functions to be called

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

############################################################################
def dist_edges(mdl,osi,bm):
	nmdl	= zip(*mdl_c)
	nosi	= zip(*osi_c)
	xm		= nmdl[0]
	ym		= nmdl[1]
	xf		= nosi[0]
	yf		= nosi[1]
	dist = []
	for n, el in enumerate(xm):
		dist2 = np.zeros(shape=0)
		dist4	= []
		for m, em in enumerate(xf):
			dist1	=	np.sqrt(pow(xm[n]-xf[m],2)+pow(ym[n]-yf[m],2))
			dist2	= np.append(dist2,dist1)
		dist3	= np.amin(dist2)
		lon,lat	= bm(xm[n],ym[n],inverse=True)
		dist4	= [dist3,lon,lat]
		dist.append(dist4)
	return dist
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

# GETTING THE NC FILES

#   # define product file
#   import time
#   from datetime import date, timedelta
#   tday  = date.today()
#   yday  = date.today() - timedelta(1)
#   yday2 = date.today() - timedelta(2)
#   cday  = tday.strftime('%Y%m%d')
#   cyear = tday.strftime('%Y')
#   pday  = yday.strftime('%Y%m%d')
#   pday2 = yday2.strftime('%Y%m%d')
print("Data files:   ")
print glob.glob("./data/*.nc")
#   date = raw_input('Insert date [YYYYMMDD]: ')
date='20150414'

####################################################################################
# DEFINING EDGE LEVELS
edge_levels = [.15]   # conc thresholds
####################################################################################

# READ TP4 DAILY
ncfil = ''.join( glob.glob('./data/TP4DAILY_start*_dump'+date+'.nc'))
print('TP4DAILY ice_only file = ' +ncfil+'\n')
slon     = 'longitude'
slat     = 'latitude'
sconc    = 'fice'
lon      = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
lat      = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
X,Y      = bm(lon[:,:],lat[:,:],inverse=False)
conc     = Mrdg.nc_get_var(ncfil,sconc,time_index=0)
Z        = conc[:,:].data
mask     = conc[:,:].mask
Z[mask]  = np.NaN

# READ IN OSI-SAF FILE
ncfil2 = './data/ice_conc_nh_polstere-100_multi_'+date+'1200.nc'
print('OSISAF file = '+ncfil2+'\n')
clon     = 'lon'
clat     = 'lat'
cconc    = 'ice_conc'
lon2     = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
lat2     = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
X2,Y2    = bm(lon2[:,:],lat2[:,:],inverse=False)
conc     = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
Z2       = conc[:,:].data
mask2    = conc[:,:].mask
Z2[mask2] = np.NaN

# MIN DIST BETWEEN EDGES (SKIMAGE)

MIN_SKI	=	1

# MIN DIST BETWEEN EDGES (BMCONTOUR)
# it's not really working, probably I haven't really understood how the bm.contour command works
# TODO

MIN_BM	= 1

if MIN_BM:
	f	= plt.figure()
	# LOOP THROUGH THE EDGE LEVELS
	for edge in edge_levels:
		g				= plt.figure()
		cs1			= bm.contour(X,Y,Z,[edge])
		coll1		= cs1.collections
		nlev1		= len(coll1)
		cs2			= bm.contour(X2,Y2,Z2,[edge])
		coll2		= cs2.collections
		nlev2		= len(coll2)
		mdl_c		= []
		osi_c		= []
		ptlist	= []
		out_list	= []
		for nl in range(nlev1):
			p     = coll1[nl].get_paths() # only one conc contour so use 0
			nseg  = len(p)
			for ns in range(nseg):
				# loop over segments
				v     = p[ns].vertices
				x     = v[:,0]
				y     = v[:,1]
				for el in v:
					mdl_c.append(el)
		for nl in range(nlev2):
			p2    = coll2[nl].get_paths()
			nseg2 = len(p2)
			for ns in range(nseg2):
				v2    = p2[ns].vertices
				x2    = v2[:,0]
				y2    = v2[:,1]
				for el in v2:
					osi_c.append(el)

		ptlist	= dist_edges(mdl_c,osi_c,bm)
							
		finish_map(bm)
		figname  = 'outputs/'+date+'_'+str(edge)+'_map.png'
		plt.savefig(figname)
		plt.close()
		g.clf()
		pts = []
		for n, el in enumerate(ptlist):
			pts.append(ptlist[n][0])
		points = np.asarray(pts)
		plt.plot(range(len(points)),points)
		textfile = 'outputs/'+date+'_'+str(edge)+'_dist_list.txt'
		np.savetxt(textfile,points)

	figname  = 'outputs/'+date+'_dist.png'
	plt.savefig(figname)
	plt.close()
	f.clf()

#      ##############################################################################
#
#         # make test plot showing ice edge contour with threshold
#         g = plt.figure()
#         Z4 = np.copy(Z2)
#         Z4[Z4<15]=0
#         Z4[Z4>=15]=1
#         thenans = np.isnan(Z4)
#         Z4[thenans] = 0
#         bm.contour(X2,Y2,Z4,1)
#         finish_map(bm)
#         figname  = odir+'/threshold_contour.png'
#         plt.savefig(figname)
#         plt.close()
#         g.clf()
#
#      ################################################################################
#
#      ################################################################################
#      # SENSIBLE POINTS WITH BM.CONTOUR
#      ################################################################################
#
#      # Plotting waves
#      f  = plt.figure()
#      bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
#      print('range in '+sswh+' (m):')
#      print(Zmin,Zmax)
#      print(' ')
#
#      cb = plt.colorbar()
#      cb.set_label("Significant Wave Height [m]",rotation=270)
#
#
#
#      ##############################################################################
#      # plot ice edge
#      # (get 15% conc contour)
#      ##############################################################################
#
#      cs1   = bm.contour(X2,Y2,Z2,[edge_level])#,'k',linewidth=2)
#      coll1 = cs1.collections
#      nlev1 = len(coll1)
#      for nl in range(nlev1):
#        p     = coll1[nl].get_paths() # only one conc contour so use 0
#        nseg  = len(p)
#        for ns in range(nseg):
#           # loop over segments
#           v     = p[ns].vertices
#           x     = v[:,0]
#           y     = v[:,1]
#           bm.plot(x,y,'m',linewidth=1.5)
#
#      ##############################################################################
#
#      # working on the waves threshold
#      Hthresh=3
#      Hmax=np.ceil(Zmax)
#      Hlev=np.arange(Hthresh,Hmax,.5)
#
#      if len(Hlev)==0:
#         #no waves over threshhold
#         out_list = []
#      else:
#         #some waves over threshhold
#         cs0   = bm.contour(X,Y,Z,Hlev)#,'k',linewidth=2)
#         print(cs0)
#         coll0 = cs0.collections
#         nlev0 = len(coll0)
#         dist_thresh=50.e3
#         out_list=[]
#         for nl in range(nlev0):
#           swh0=cs0.levels[nl]
#           p     = coll0[nl].get_paths() # only one conc contour so use 0
#           nseg  = len(p)
#           for ns in range(nseg):
#              # loop over segments
#              v     = p[ns].vertices
#              x     = v[:,0]
#              y     = v[:,1]
#
#              #loop over all ice edges
#              for nl2 in range(nlev1):
#                p2    = coll0[nl2].get_paths() # only one conc contour so use 0
#                nseg2 = len(p2)
#                for ns2 in range(nseg2):
#                   # loop over segments
#                   v2    = p2[ns2].vertices
#                   x2    = v2[:,0]
#                   y2    = v2[:,1]
#                   dist,lon_min,lat_min=dist_cont2cont(x2,y2,x,y,bm) # output (lon,lat) for ice
#                   if dist<dist_thresh:
#                     dist_list=[dist,lon_min,lat_min,swh0]
#                     out_list.append(dist_list)
#
#      ##############################################################################
#
#      if 1:
#         # add test point to plot
#         # - to check if SAR image is ordered in the right place
#         # - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
#         #   then use trial and error
#         nout = len(out_list)
#         for mm in range(nout):
#          list0=out_list[mm]
#          if list0[3] >= 3 and list0[0] <= 5:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#          elif list0[3] >= 4 and list0[0] <= 20:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#          elif list0[3] >= 5 and list0[0] <= 50:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#
#      finish_map(bm)
#
#      # date+time to title and file name
#      # label '$H_s$, m' to colorbar
#      wav = Dataset(ncfil)
#      fnday = getattr(wav,'max_time')
#      fnday = fnday.replace(" ", "_")
#      fnday = fnday.replace(":", "")
#      fnday = fnday.replace("-", "")
#      figname  = odir+'/img/'+fnday+'.png'
#      plt.savefig(figname)
#      print('saving figure:')
#      print(figname+'\n')
#      plt.close()
#      f.clf()
#      # bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function
#
#
#      # swh to find when there are large waves (>4m) in the vicinity of the ice.
#      # write and send an email to warn when this happens (so we can order some SAR images)
#
#      # filename of text file to form contents of email message
#      textfile = odir+'/lst/'+fnday+'_list.txt'
#      nout=len(out_list)
#      if nout>0:
#        SEND_EMAIL2=1
#        tf=open(textfile,'w')
#
#        for mm in range(nout):
#          list0=out_list[mm]
#          line='\n' #get info from out_list
#          if list0[3] >= 3 and list0[0] <= 5:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#          elif list0[3] >= 4 and list0[0] <= 20:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#          elif list0[3] >= 5 and list0[0] <= 50:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#
#          tf.write(line)
#
#        tf.close()
#      else:
#        SEND_EMAIL2=0
#      # f  = open(textfile,'w')
#      # f.write('Hei Tim!')
#      # f.close()
#
#      #############################################################################
#
#      ##############################################################################
#      # SENSIBLE POINTS WITH THRESHOLD
#      ##############################################################################
#
#      # Plotting waves
#      f  = plt.figure()
#      bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
#      print('range in '+sswh+' (m):')
#      print(Zmin,Zmax)
#      print(' ')
#
#      cb = plt.colorbar()
#      cb.set_label("Significant Wave Height [m]",rotation=270)
#
#      ##############################################################################
#      # plot ice edge
#      # (get 15% conc contour)
#      ##############################################################################
#
#      # ice pack
#      Z3 = np.copy(Z2)
#      Zval = np.ma.array(Z3)
#      Z3 = np.ma.masked_where(Zval < 15, Zval)
#      bm.pcolor(X2,Y2,Z3,cmap='Greys')
#
#      # ice edge contour with threshold
#      Z4 = np.copy(Z2)
#      Z4[Z4<15]=0
#      Z4[Z4>=15]=1
#      thenans = np.isnan(Z4)
#      Z4[thenans] = 0
#      cs4 = bm.contour(X2,Y2,Z4,1)
#      coll4 = cs4.collections
#      nlev4 = len(coll4)
#      print nlev4
#      for nl in range(nlev4):
#        p4    = coll4[nl].get_paths() # only one conc contour so use 0
#        nseg4 = len(p4)
#        print nseg4
#        for ns in range(nseg4):
#           v4 = p4[ns].vertices
#           x4 = v4[:,0]
#           y4 = v4[:,1]
#           bm.plot(x4,y4,'k',linewidth=1.5)
#
#           ##############################################################################
#
#           # working on the waves threshold
#           Hthresh=3
#           Hmax=np.ceil(Zmax)
#           Hlev=np.arange(Hthresh,Hmax,.5)
#
#           if len(Hlev)==0:
#              #no waves over threshhold
#              out_list = []
#           else:
#              #some waves over threshhold
#              cs0   = bm.contour(X,Y,Z,Hlev)
#              print(cs0)
#              coll0 = cs0.collections
#              nlev0 = len(coll0)
#              dist_thresh=50.e3
#              out_list=[]
#              for nl in range(nlev0):
#                swh0=cs0.levels[nl]
#                p     = coll0[nl].get_paths() # only one conc contour so use 0
#                nseg  = len(p)
#                for ns in range(nseg):
#                   # loop over segments
#                   v     = p[ns].vertices
#                   x     = v[:,0]
#                   y     = v[:,1]
#
#                   #loop over all ice edges
#                   for nl2 in range(nlev1):
#                     p2    = coll0[nl2].get_paths() # only one conc contour so use 0
#                     nseg2 = len(p2)
#                     for ns2 in range(nseg2):
#                        # loop over segments
#                        v2    = p2[ns2].vertices
#                        x2    = v2[:,0]
#                        y2    = v2[:,1]
#                        dist,lon_min,lat_min=dist_cont2cont(x2,y2,x,y,bm) # output (lon,lat) for ice
#                        if dist<dist_thresh:
#                          dist_list=[dist,lon_min,lat_min,swh0]
#                          out_list.append(dist_list)
#
#         ##############################################################################
#
#      if 1:
#         # add test point to plot
#         # - to check if SAR image is ordered in the right place
#         # - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
#         #   then use trial and error
#         nout = len(out_list)
#         for mm in range(nout):
#          list0=out_list[mm]
#          if list0[3] >= 3 and list0[0] <= 5:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#          elif list0[3] >= 4 and list0[0] <= 20:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#          elif list0[3] >= 5 and list0[0] <= 50:
#             lon_plot = list0[1]
#             lat_plot = list0[2]
#             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
#             x_plot,y_plot  = bm(lon_plot,lat_plot)
#             bm.plot(x_plot,y_plot,'og',markersize=5)
#
#      finish_map(bm)
#
#      # date+time to title and file name
#      # label '$H_s$, m' to colorbar
#      wav = Dataset(ncfil)
#      fnday = getattr(wav,'max_time')
#      fnday = fnday.replace(" ", "_")
#      fnday = fnday.replace(":", "")
#      fnday = fnday.replace("-", "")
#      figname  = odir+'/img/'+fnday+'_threshold.png'
#      plt.savefig(figname)
#      print('saving figure:')
#      print(figname+'\n')
#      plt.close()
#      f.clf()
#      # bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function
#
#
#      # swh to find when there are large waves (>4m) in the vicinity of the ice.
#      # write and send an email to warn when this happens (so we can order some SAR images)
#
#      # filename of text file to form contents of email message
#      textfile = odir+'/lst/'+fnday+'_threshold_list.txt'
#      nout=len(out_list)
#      if nout>0:
#        SEND_EMAIL2=1
#        tf=open(textfile,'w')
#
#        for mm in range(nout):
#          list0=out_list[mm]
#          line='\n' #get info from out_list
#          if list0[3] >= 3 and list0[0] <= 5:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#          elif list0[3] >= 4 and list0[0] <= 20:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#          elif list0[3] >= 5 and list0[0] <= 50:
#             for ii in range(3,-1,-1):
#                line=3*' '+str(list0[ii])+line
#
#          tf.write(line)
#
#        tf.close()
#      else:
#        SEND_EMAIL2=0
#      # f  = open(textfile,'w')
#      # f.write('Hei Tim!')
#      # f.close()
##   else:
##      print 'ERROR: Please type 0 for bm.contour or 1 for threshold'
##      exit
#
#################################################################
## EMAIL SYSTEM
#if SEND_EMAIL and SEND_EMAIL2:
#   # import smtplib for the actual sending function
#   import smtplib
#
#   # Import the email modules we'll need
#   from email.mime.text import MIMEText
#
#   if CHECK_NC==0:
#      textfile = 'message.txt'
#      f  = open(textfile,'w')
#      f.write('Hei Tim!')
#      f.close()
#
#   # Open a plain text file for reading.  For this example, assume that
#   # the text file contains only ASCII characters.
#   fp = open(textfile, 'rb')
#   # Create a text/plain message
#   msg = MIMEText(fp.read())
#   fp.close()
#
#   COMMASPACE  = ', '
#   sender      = 'timill@hpc.uib.no'            # sender's email
#   receivers   = ['timothy.williams@nersc.no']  # list of receivers' emails
#
#   msg['Subject'] = 'The contents of %s' % textfile
#   msg['From']    = sender
#   msg['To']      = COMMASPACE.join(receivers)
#
#   # Send the message via our own SMTP server, but don't include the
#   # envelope header.
#   s = smtplib.SMTP('localhost')
#   s.sendmail(sender, receivers, msg.as_string())
#   s.quit()
#################################################################
