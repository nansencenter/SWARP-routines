# Set up email alert for when WAMNSEA waves become very large
# - possibly limit search to near ice?

from netCDF4 import Dataset
import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

sys.path.append('../py_funs')
import mod_reading as Mrdg
# import basemap_gridlines as bmg #TODO - put Tim's basemap_gridlines function in py_funs

SEND_EMAIL  = 0
CHECK_NC    = 1

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
   cday = tday.strftime('%Y%m%d')
   pday = yday.strftime('%Y%m%d')
   wmsc  = '/work/shared/nersc/msc/WAMNSEA/'
   ncfil = wmsc + 'wam_nsea.fc.' + cday + '.nc' # should be determined from today's date use "fc"
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
   lon      = Mrdg.nc_get_var(ncfil,slon) # lon[:,:] is a numpy array
   lat      = Mrdg.nc_get_var(ncfil,slat) # lat[:,:] is a numpy array
   X,Y      = bm(lon[:,:],lat[:,:],inverse=False)
   in_area  = np.logical_and(abs(X)<xmax,abs(Y)<ymax)

   # TODO read in OSISAF conc file
   # plot conc + ice edge (15% bm.pcontour? )
   # get lon/lat and restrict to relevant area
   osisaf = '/work/shared/nersc/msc/OSI-SAF/' + str(yday.year) + '_nh_polstere'
   ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday + '1200.nc'
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


   # for loop_i in range(Ntimes):
   for loop_i in [0]:
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
      print(cs1)
      coll1 = cs1.collections
      p     = coll1[0].get_paths() # only one conc contour so use 0 
      nseg  = len(p)
      for ns in range(nseg):
         # loop over segments
         v     = p[ns].vertices
         x     = v[:,0]
         y     = v[:,1]
         bm.plot(x,y,'--m',linewidth=1.5)
      ##############################################################################

      ##############################################################################
      #
      bm.drawcoastlines()
      bm.fillcontinents(color='gray')

      # draw parallels and meridians.
      bm.drawparallels(np.arange(60.,91.,10.),\
            labels=[True,False,True,True]) # labels = [left,right,top,bottom]
      bm.drawmeridians(np.arange(-180.,181.,20.),latmax=90.,\
            labels=[True,False,False,True])
      bm.drawmapboundary() # fill_color='aqua')
      #
      figname  = 'test.png'
      plt.savefig(figname)
      print('saving figure:')
      print(figname+'\n')
      plt.close()
      f.clf()
      # bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function


      # TODO search in swh to find when there are large waves (>4m) in the vicinity of the ice.
      # TODO write and send an email to warn when this happens (so we can order some SAR images)

   # filename of text file to form contents of email message 
   textfile = 'message.txt'
   # f  = open(textfile,'w')
   # f.write('Hei Tim!')
   # f.close()

################################################################
if SEND_EMAIL:
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
