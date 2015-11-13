# Set up email alert for when WAMNSEA waves become very large
# - possibly limit search to near ice?

from netCDF4 import Dataset
import os,sys
import numpy as np
import time
from datetime import date, timedelta

############################################################################
WANT2SHOW   = 0 # show the figure - make sure you don't use this if you need to run inside crontab
if not WANT2SHOW:
   # don't need to show figure
   # - this lets the script be run through crontab or on compute node
   import matplotlib
   matplotlib.use('Agg')
############################################################################

from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure as msr

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import mod_reading as Mrdg
import fns_plotting as Fplt

############################################################################
# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')
############################################################################

SEND_EMAIL  = 1     # send email alert if there are big waves close to ice
SEND_EMAIL3 = 0     # change to 1 if threshhold method says there are big waves close to ice

tday  = date.today()
cday  = tday.strftime('%Y%m%d')
#
yday  = date.today() - timedelta(1)
pyear = tday.strftime('%Y')
pday  = yday.strftime('%Y%m%d')
#
yday2 = date.today() - timedelta(2)
pday2 = yday2.strftime('%Y%m%d')

# where to put outputs
odir  = '/work/timill/RealTime_Models/check_ww3arctic/'
if not os.path.exists(odir):
   os.mkdir(odir)
# daily sub-dir
odir  = odir+pday
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

############################################################################
def binary_cont(data,thresh):
   Z              = np.zeros(data.shape,dtype=int)
   Z[Z>=thresh]   = 1
   return Z
############################################################################

############################################################################
def binary_diff(data,thresh1,thresh2):
   Z1 = np.zeros(data.shape,dtype=int)
   Z2 = np.zeros(data.shape,dtype=int)
   #
   Z1[Z1>=thresh1]   = 1
   Z2[Z2>=thresh2]   = 1
   #
   D = Z2 - Z1
   return D,Z1,Z2
############################################################################

############################################################################
def Hs_filter():
   # filter some contours out, depending on wave height and distance to ice edge
   # - distance to threshhold
   thrdic      = {5  :50}
   thrdic.update({4.5:35})
   thrdic.update({4  :20})
   thrdic.update({3.5:10})
   thrdic.update({3  :5 })

   # - symbol
   symdic      = {5  :'+g'}
   symdic.update({4.5:'^g'})
   symdic.update({4  :'vg'})
   symdic.update({3.5:'og'})
   symdic.update({3  :'*g'})

   # - marker size
   ms_dic      = {5  :7}
   ms_dic.update({4.5:5})
   ms_dic.update({4  :5})
   ms_dic.update({3.5:5})
   ms_dic.update({3  :7})
   return thrdic,symdic,ms_dic
############################################################################


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
############################################################################


#####################################################################################
# define netcdf file  
#####################################################################################

wmsc = '/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/'+ pyear + '/forecast/'
ncfil = wmsc + 'SWARP_WW3_ARCTIC-12K_'+pday+'.fc.nc'
print('\nWW3 Arctic file = ' + ncfil+'\n')
	
# get info about nc file
nci      = Mrdg.nc_getinfo(ncfil)
times    = nci.timevalues # hours from 1st time in file
Ntimes   = len(times)
irec     = Ntimes-1
print('Time values (h):')
print(times)
print(' ')
	
# get lon/lat and restrict to relevant area
sswh     = 'hs'
lon,lat  = nci.get_lonlat()
X,Y      = bm(lon,lat,inverse=False)
#
in_area  = np.logical_and(X>bm.xmin,X<bm.xmax)
in_area  = np.logical_and(in_area,Y<bm.ymax)
in_area  = np.logical_and(in_area,Y>bm.ymin)
	
# read in yesterday's OSISAF conc file
# plot conc + ice edge (15% bm.pcontour? )
# get lon/lat and restrict to relevant area
osisaf = '/work/shared/nersc/msc/OSI-SAF/' + str(yday.year) + '_nh_polstere'
ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday + '1200.nc'

# use day before yesterday's conc if yesterday's not present
if not os.path.exists(ncfil2):
   osisaf = '/work/shared/nersc/msc/OSI-SAF/' + str(yday2.year) + '_nh_polstere'
   ncfil2 = osisaf + '/ice_conc_nh_polstere-100_multi_' + pday2 + '1200.nc'
   # ncfil2 = 'ice_conc_nh_polstere-100_multi_201507061200.nc'

print('OSISAF file = '+ncfil2+'\n')
nci2        = Mrdg.nc_getinfo(ncfil2)
clon        = 'lon'
clat        = 'lat'
cconc       = 'ice_conc'
edge_level  = 15
lon2,lat2   = nci2.get_lonlat()
X2,Y2       = bm(lon2,lat2,inverse=False)
#in_area = np.logical_and(abs(X)<xmax,abs(Y)<ymax)
conc        = nci2.get_var(cconc,time_index=0)
Z2          = conc.values.data
mask2       = conc.values.mask
Z2[mask2]   = np.NaN
	
	
#check_list=0
#check_list=range(Ntimes)
check_list  = [irec]
for loop_i in check_list:
   swh = nci.get_var(sswh,time_index=loop_i)
   # swh.values is a numpy masked array

   Z     = swh.values.data
   mask  = swh.values.mask
   Zmax  = swh.values[in_area].max() # restrict max calc to relevant area
   Zmin  = swh.values[in_area].min() # restrict max calc to relevant area
   print(Zmin,Zmax)
   
   # make plot
   Z[mask] = np.NaN
   # if mask is 'True' (data is missing or invalid) set to NaN

##############################################################################
		
		
################################################################################
		
   #################################################################################
   if 1:
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
      Z3    = np.copy(Z2)
      Zval  = np.ma.array(Z3)
      Z3    = np.ma.masked_where(Zval < 15, Zval)
      bm.pcolor(X2,Y2,Z3,cmap='RdPu')
      
      # ice edge contour with threshold 
      Z4          = np.copy(Z2)
      Z4[Z4<15]   = 0
      Z4[Z4>=15]  = 1
      thenans     = np.isnan(Z4)
      Z4[thenans] = 0
      cs4         = bm.contour(X2,Y2,Z4,1)
      coll4       = cs4.collections
      nlev4       = len(coll4)
      for nl4 in range(nlev4):
         p4     = coll4[nl4].get_paths() # only one conc contour so use 0 
         nseg4  = len(p4)
         for ns4 in range(nseg4):
            v4 = p4[ns4].vertices
            x4 = v4[:,0]
            y4 = v4[:,1]
            bm.plot(x4,y4,'m',linewidth=1.5)
            
            # working on the waves threshold
            # Creating array with .5 interval from 3 meters to Max Wave Height
            # IF MWH < 3 meters --> empty array
            Hthresh  = 3
            Hmax     = np.ceil(Zmax)
            Hlev     = np.arange(Hthresh,Hmax,.5)
            
            if len(Hlev)==0:
               #no waves over threshhold
               out_list = []
            elif 1:
               #some waves over threshhold
               # - use basemap.contour
               cs3         = bm.contour(X,Y,Z,Hlev)
               coll3       = cs3.collections
               nlev3       = len(coll3)
               dist_thresh = 50.e3
               out_list    = []

               # loop over Hs values
               for nl3 in range(nlev3):
                  swh3  = cs3.levels[nl3]
                  p3    = coll3[nl3].get_paths() # only one conc contour so use 0 
                  nseg3 = len(p3)
                  for ns3 in range(nseg3):
                     # loop over segments
                     v3 = p3[ns3].vertices
                     x3 = v3[:,0]
                     y3 = v3[:,1]
                     
                     #loop over all ice edges
                     for nl4 in range(nlev4):
                        p4    = coll4[nl4].get_paths() # only one conc contour so use 0 
                        nseg4 = len(p4)
                        for ns4 in range(nseg4):
                           # loop over segments
                           v4 = p4[ns4].vertices
                           x4 = v4[:,0]
                           y4 = v4[:,1]
                           dist,lon_min,lat_min=dist_cont2cont(x4,y4,x3,y3,bm) # output (lon,lat) for ice
                           if dist<dist_thresh:
                              dist_list   = [dist,lon_min,lat_min,swh3]
                              out_list.append(dist_list)

               # add labels - this shouldn't affect above results since bm.contour has been called already
               #              and the results processed
               texth = plt.clabel(cs3,inline=1,fontsize=10,fmt='%1.1f')

            # elif 1:
            #    #some waves over threshhold
            #    # - use skimage
            #    import skimage.measure as msr

            #    # loop over Hs values
            #    for Hs in Hlev:
            #       B1       = 0.*Z
            #       B1[Z>Hs] = 1.
            #       conts    = msr.find_contours(B1,.5)

            #       for v3 in conts:
            #          x3 = v3[:,0]
            #          y3 = v3[:,1]
            #          
            #          #loop over all ice edges
            #          for nl4 in range(nlev4):
            #             p4    = coll4[nl4].get_paths() # only one conc contour so use 0 
            #             nseg4 = len(p4)
            #             for ns4 in range(nseg4):
            #                # loop over segments
            #                v4 = p4[ns4].vertices
            #                x4 = v4[:,0]
            #                y4 = v4[:,1]
            #                dist,lon_min,lat_min=dist_cont2cont(x4,y4,x3,y3,bm) # output (lon,lat) for ice
            #                if dist<dist_thresh:
            #                   dist_list   = [dist,lon_min,lat_min,swh3]
            #                   out_list.append(dist_list)

            #    # plot with bm.contour (plot numbers on figure without changing the v3 conts)
            #    CS = bm.contour(X,Y,Z,Hlev)
      ##############################################################################
			
      ##########################################################################################
      if 1:
         # plot nearest points on ice edge
         nout = len(out_list)
         sym_skip = [] # skip some symbols (plots can get crowded)
         # sym_skip = ['*'] # skip some symbols (plots can get crowded)
         if nout==0:
            print('No large waves close to ice\n')
         else:
            for mm in range(nout):
               dist_list                  = 1*out_list[mm]
               dist_list[0]               = dist_list[0]/1.e3 #km
               dist,lon_plot,lat_plot,Hs  = dist_list
            
               # filter some out, depending on wave height and distance to ice edge
               thrdic,symdic,ms_dic = Hs_filter()
               for Hsc in  thrdic.keys():
                  if Hs >= Hsc and dist <= thrdic[Hsc]:
                     if symdic[Hsc][0] not in sym_skip:
                        print('Adding test point ('+str(lon_plot)+'E,'+str(lat_plot)+'N)\n')
                        bm.plot(lon_plot,lat_plot,symdic[Hsc],markersize=ms_dic[Hsc],latlon=True)
      ##########################################################################################

      ##########################################################################################
      if 0:
         # add manual test point(s) to plot
         # - to check if SAR image is ordered in the right place
         # - get initial estimate from ncview (use OSISAF file not wamnsea - lon/lat are weird in those files),
         #   then use trial and error
         print('Adding manual points to plot:')

         man_list = [[40.,82.5]] # list of lon/lat points
         for lonm,latm in man_list:
            print('('+str(lonm)+' E, '+str(latm)+' N)')
            bm.plot(lonm,latm,'^g',markersize=7,latlon=True)
         print('\n')
      ##########################################################################################

      finish_map(bm)

      # date+time to title and file name
      # label '$H_s$, m' to colorbar
      dt       = nci.timeval_to_datetime(nci.timevalues[loop_i])
      fnday    = dt.strftime("%Y%m%dT%H%M%SZ")
      figname  = odir+'/img/'+fnday+'.png'
      if not os.path.exists(odir+'/img'):
         os.mkdir(odir+'/img')
      print('saving figure:')
      print(figname+'\n')
      plt.savefig(figname)

      if WANT2SHOW:
         # show figure to enable zooming etc
         print('Showing figure - close it to continue\n')
         plt.show()

      plt.close()
      f.clf()
      # bmg.latlon_grid(bm,10.,10.) #TODO - get Tim's basemap_gridlines function
      
      
      # swh to find when there are large waves (>4m) in the vicinity of the ice.
      # write and send an email to warn when this happens (so we can order some SAR images)
      
      # filename of text file to form contents of email message 
      textfile = odir+'/lst/'+fnday+'_list.txt'
      if not os.path.exists(odir+'/lst'):
         os.mkdir(odir+'/lst')
      nout=len(out_list)
      SEND_EMAIL3 = (nout>0)
      if SEND_EMAIL3:
         tf    = open(textfile,'w')
         hdr   = ['Dist (km)','lon (deg N)','lat (deg E)','Hs (m)']
         blk   = 3*' '
         line  = ''
         ordr  = [3,0,1,2]
         for ii in ordr:
            line  = line+blk+hdr[ii]
         tf.write(line+'\n')
         print(line)
         
         thrdic   = Hs_filter()[0]
         for mm in range(nout):
            #get info from out_list
            dist_list         = 1*out_list[mm]
            dist,lon,lat,Hs   = dist_list
            dist              = dist/1.e3 # dist now in km
            dist_list         = [dist,lon,lat,Hs]
            line              = ''
            if 0:
               # just print everything
               for ii in ordr:
                  line  = line+blk+str(dist_list[ii])
               print(line)
               tf.write(line+'\n')
            else:
               # filter some out, depending on wave height and distance to ice edge
               for Hsc in  thrdic.keys():
                  if Hs >= Hsc and dist <= thrdic[Hsc]:
                     for ii in ordr:
                        line  = line+blk+str(dist_list[ii])
                     print(line)
                     tf.write(line+'\n')

         tf.close()
   #################################################################################


# ################################################################
# # EMAIL SYSTEM
# if SEND_EMAIL and SEND_EMAIL3:
#    import subprocess
#    subprocess.call(["chmod",'-R',"+rw",odir+"/img"])
#    subprocess.call(["chmod",'-R',"+rw",odir+"/lst"])
#    #
#    awdir = SR+'/forecast_scripts/alert_waves_ww3arctic'
#    subprocess.check_call([awdir+'/waves_alert.sh', pday])
# elif SEND_EMAIL:
#    import subprocess
#    # just send pic
#    subprocess.call(["chmod",'-R',"+rw",odir+"/img"])
#    #
#    awdir = SR+'/forecast_scripts/alert_waves_ww3arctic'
#    subprocess.check_call([awdir+'/waves_pic.sh', pday])
# ################################################################
