# Set up email alert for when WAMNSEA waves become very large
# - possibly limit search to near ice?

from netCDF4 import Dataset
import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure

sys.path.append('../../py_funs')
import mod_reading as Mrdg

SEND_EMAIL  = 0
CHECK_NC    = 1
odir        = 'out' # where to put temporary outputs

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
def dist_cont2cont(xice,yice,xwav,ywav,bm):

  n1=len(xice)
  dist2=np.zeros(n1)
  for n in range(n1):
    x1n=xice[n]
    y1n=yice[n]
    dist1=np.sqrt(pow(xwav-x1n,2)+pow(ywav-y1n,2))
    dist2[n]=dist1.min()
    
  imin=dist2.argmin()
  dist=dist2[imin] # want xmin,ymin on ice corresponding to dist
  xmin = xice[imin]
  ymin = yice[imin]
  lon,lat    = bm(xmin,ymin,inverse=True)
  return dist,lon,lat
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
   tday  = date.today()
   yday  = date.today() - timedelta(1)
   yday2 = date.today() - timedelta(2)
   cday  = tday.strftime('%Y%m%d')
   pday  = yday.strftime('%Y%m%d')
   pday2 = yday2.strftime('%Y%m%d')
   wmsc  = '/work/shared/nersc/msc/WAMNSEA/'
   ncfil = wmsc + 'wam_nsea.fc.' + cday + '.nc' # should be determined from today's date use "fc"
   #ncfil = 'wam_nsea.fc.20150413.nc'
   print('WAMNSEA file = ' + ncfil+'\n')

   # get info about nc file
   ncinfo_wav  = Mrdg.nc_getinfo(ncfil)
   times       = ncinfo_wav.timevalues # hours from 1st time in file
   Ntimes      = len(times)
   print('Time values (h):')
   print(times)
   print(' ')

   # get lon/lat and restrict to relevant area
   slon     = 'longitude'
   slat     = 'latitude'
   sswh     = 'significant_wave_height'
   fnday    = 'max_time'
   lon      = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
   lat      = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
   X,Y      = bm(lon[:,:],lat[:,:],inverse=False)
   in_area  = np.logical_and(abs(X)<xmax,abs(Y)<ymax)

   # read in OSISAF conc file
   # plot conc + ice edge (15% bm.pcontour? )
   # get lon/lat and restrict to relevant area
   osisaf = '/work/shared/nersc/msc/OSI-SAF/' + str(yday.year) + '_nh_polstere'
   ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday + '1200.nc'
   if not os.path.exists(ncfil2):
      ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday2 + '1200.nc'
   #ncfil2 = 'ice_conc_nh_polstere-100_multi_201504121200.nc'
   print('OSISAF file = '+ncfil2+'\n')
   clon     = 'lon'
   clat     = 'lat'
   cconc    = 'ice_conc'
   edge_level = 15
   lon2     = Mrdg.nc_get_var(ncfil2,clon) # lon[:,:] is a numpy array
   lat2     = Mrdg.nc_get_var(ncfil2,clat) # lat[:,:] is a numpy array
   X2,Y2    = bm(lon2[:,:],lat2[:,:],inverse=False)
   #in_area  = np.logical_and(abs(X)<xmax,abs(Y)<ymax)
   conc     = Mrdg.nc_get_var(ncfil2,cconc,time_index=0)
   Z2       = conc[:,:].data
   mask2    = conc[:,:].mask
   Z2[mask2] = np.NaN


   #check_list=0
   #check_list=range(Ntimes)
   check_list=[Ntimes-1]
   for loop_i in check_list:
      swh   = Mrdg.nc_get_var(ncfil,sswh,time_index=loop_i)
         # swh[:,:] is a numpy masked array
      Z        = swh[:,:].data
      mask     = swh[:,:].mask
      Zmax     = swh[in_area].max() # restrict max calc to relevant area
      Zmin     = swh[in_area].min() # restrict max calc to relevant area

      ##############################################################################
      # make plot
      Z[mask]  = np.NaN
         # if mask is 'True' (data is missing or invalid) set to NaN

      ##############################################################################
      if 1:
        # make test plot showing ice edge contour with bm.contour
        g  = plt.figure()
        bm.pcolor(X2,Y2,Z2,vmin=0,vmax=100)
        # plot ice edge
        # (get 15% conc contour)
        # http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines
        cs1   = bm.contour(X2,Y2,Z2,[edge_level])#,'k',linewidth=2)
        print(cs1)
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
             bm.plot(x,y,'k',linewidth=1)

        # make test plot showing ice edge contour with skimage.measure.find_contours
        conc = measure.find_contours(Z2,edge_level)
        for el in conc:
          ell = np.array(el).tolist()
          for n in range(len(ell)):
            if np.str(ell[n][0]) == 'nan' or np.str(ell[n][1]) == 'nan':
              print 'nan'
            else:
              bm.plot([ell[n][0]],[ell[n][1]],'r',linewidth=1 ) 

        finish_map(bm)
        figname  = odir+'/test2.png'
        plt.savefig(figname)
        plt.close()
        g.clf()
      ##############################################################################

      f  = plt.figure()
      bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
      print('range in '+sswh+' (m):') 
      print(Zmin,Zmax)
      print(' ')

      bm.colorbar()
      ##############################################################################

      ##############################################################################
      # plot ice edge
      # (get 15% conc contour)
      # http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines
      cs1   = bm.contour(X2,Y2,Z2,[edge_level])#,'k',linewidth=2)
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
           bm.plot(x,y,'--m',linewidth=1.5)
      ##############################################################################

      ##############################################################################
      Hthresh=3
      Hmax=np.ceil(Zmax)
      Hlev=np.arange(Hthresh,Hmax,.5)

      if len(Hlev)==0:
         #no waves over threshhold
         out_list = []
      else:
         #some waves over threshhold
         cs0   = bm.contour(X,Y,Z,Hlev)#,'k',linewidth=2)
         print(cs0)
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
              v     = p[ns].vertices
              x     = v[:,0]
              y     = v[:,1]

              #loop over all ice edges
              for nl2 in range(nlev1):
                p2    = coll0[nl2].get_paths() # only one conc contour so use 0 
                nseg2 = len(p2)
                for ns2 in range(nseg2):
                   # loop over segments
                   v2    = p2[ns2].vertices
                   x2    = v2[:,0]
                   y2    = v2[:,1]
                   dist,lon_min,lat_min=dist_cont2cont(x2,y2,x,y,bm) # output (lon,lat) for ice
                   if dist<dist_thresh:
                     dist_list=[dist,lon_min,lat_min,swh0]
                     out_list.append(dist_list)
              #bm.plot(x,y,'--g',linewidth=1.5)

      ##############################################################################
      #
      if 1:
         # add test point to plot
         # - to check if SAR image is ordered in the right place
         # - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
         #   then use trial and error
         # TODO add automatic estimate for SAR image? (plot points on ice edge that have high waves near them?)
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
             bm.plot(x_plot,y_plot,'og',markersize=5)
          elif list0[3] >= 5 and list0[0] <= 50:
             lon_plot = list0[1]
             lat_plot = list0[2]
             print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
             x_plot,y_plot  = bm(lon_plot,lat_plot)
             bm.plot(x_plot,y_plot,'og',markersize=5)

      finish_map(bm)

      #TODO add date+time to title and file name
      #TODO add label '$H_s$, m' to colorbar
      wav = Dataset(ncfil)
      fnday = getattr(wav,'max_time')
      fnday = fnday.replace(" ", "_")
      fnday = fnday.replace(":", "")
      fnday = fnday.replace("-", "")
      figname  = odir+'/'+fnday+'.png'
      plt.savefig(figname)
      print('saving figure:')
      print(figname+'\n')
      plt.close()
      f.clf()
      # bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function


      # swh to find when there are large waves (>4m) in the vicinity of the ice.
      # write and send an email to warn when this happens (so we can order some SAR images)

   # filename of text file to form contents of email message 
   textfile = odir+'/'+fnday+'_list.txt'
   nout=len(out_list)
   if nout>0:
     SEND_EMAIL2=1
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
   else:
     SEND_EMAIL2=0
   # f  = open(textfile,'w')
   # f.write('Hei Tim!')
   # f.close()

################################################################
if SEND_EMAIL and SEND_EMAIL2:
   # import smtplib for the actual sending function
   import smtplib

   # Import the email modules we'll need
   from email.mime.text import MIMEText

   if CHECK_NC==0:
      textfile = 'message.txt'
      f  = open(textfile,'w')
      f.write('Hei Tim!')
      f.close()

   # Open a plain text file for reading.  For this example, assume that
   # the text file contains only ASCII characters.
   fp = open(textfile, 'rb')
   # Create a text/plain message
   msg = MIMEText(fp.read())
   fp.close()

   COMMASPACE  = ', '
   sender      = 'timill@hpc.uib.no'            # sender's email
   receivers   = ['timothy.williams@nersc.no']  # list of receivers' emails

   msg['Subject'] = 'The contents of %s' % textfile
   msg['From']    = sender
   msg['To']      = COMMASPACE.join(receivers)

   # Send the message via our own SMTP server, but don't include the
   # envelope header.
   s = smtplib.SMTP('localhost')
   s.sendmail(sender, receivers, msg.as_string())
   s.quit()
################################################################
