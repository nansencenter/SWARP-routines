# Set up email alert for when WAMNSEA waves become very large
# - possibly limit search to near ice?

from netCDF4 import Dataset
import sys,os

sys.path.append('../netcdf_to_grib2')
import mod_reading as mr

if 0:
   # define netcdf file  
   wmsc  = '/work/shared/nersc/msc/WAMNSEA/'
   ncfil = wmsc+'wam_nsea.an.20141101.nc'#TODO should be determined from today's date

   # get lon/lat and swh
   lon   = mr.nc_get_var(ncfil,'longitude') # lon[:,:] is a numpy array
   lat   = mr.nc_get_var(ncfil,'latitude')  # lat[:,:] is a numpy array
   swh   = mr.nc_get_var(ncfil,'significant_wave_height') # swh[i,:,:] is a numpy masked array

   # TODO search in swh to find when there are large waves (>4m) in the vicinity of the ice.
   # TODO write and send an email to warn when this happens (so we can order some SAR images)

   # filename of text file to form contents of email message 
   textfile = 'message.txt'
   if 1:
      f  = open(textfile,'w')
      f.write('Hei Tim!')
      f.close()
   else:
      # write a proper message to warn there are big waves near the ice
      print(' ')
else:
   textfile = 'test.txt'

################################################################
# import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.mime.text import MIMEText

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
