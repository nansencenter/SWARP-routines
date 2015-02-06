import os,sys
import numpy as np
from datetime import datetime
#
import pygrib
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d
#
import mod_reading as m_rdg
import mod_grib2_setup as m_g2s

ncgv        = m_rdg.nc_get_var
KEEP_MASK   = 1
DO_TEST     = 0

#######################################################################
if DO_TEST==1:
   # test options
   tindx    = 0#time index to test
   vindx    = 0#variable index to test
   #
   PLOT_OPT = 0 # 0: no plot
                # 1: simple plot with imshow
                # 2: more complicated plot with basemap
#######################################################################

#######################################################################
# file inputs:
ncfil = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc" # hyc2proj netcdf

# output file:
outdir  = 'out'
if not os.path.exists(outdir):
   os.makedirs(outdir)
fil_out = outdir+'/test_hyc2proj_to_grib2.grb2'

# do conversion
ncinfo   = m_g2s.hyc2proj_to_grib2(ncfil,fil_out,KEEP_MASK=KEEP_MASK)

#######################################################################
if DO_TEST==1:
   Nt             = ncinfo.number_of_time_records
   vbl_list       = ncinfo.variable_list
   time_indices   = range(Nt)
   Nv             = len(vbl_list)
   time_index     = time_indices[tindx]
   vbl_name       = vbl_list    [vindx]
   #
   print('Test write to '+fil_out)
   print("Doing test of encoding of variable "+vbl_name)
   print("at time record number "+str(time_index))
   print(' ')
   data  = ncgv(ncfil,vbl_name,time_index)
   if KEEP_MASK==1:
      data_arr = data[:,:] # masked array
   else:
      data_arr = data[:,:].data # array

   grb_index   = 1+Nv*time_index+vindx
   gr          = pygrib.open(fil_out)
   grbmsgs     = gr.read(grb_index) # get list of 1st "grb_index" messages:

   for msg in [grbmsgs[-1]]:
      #last message
      print(msg)
      print('\n')
      grb   = g2d(msg.tostring(),gribmsg=True)
      print(grb)

      # lon,lat from grib2
      glat,glon   = msg.latlons()

      # lat,lon from nc file
      lon   = ncgv(ncfil,'longitude')
      lat   = ncgv(ncfil,'latitude')
      nlon  = lon[:,:]
      nlat  = lat[:,:]

      #################################################################
      # calc difference between grid arrays:
      print('Comparing grids between netcdf and grib2 files...\n')

      # small errors in lon can make 1 set be around -180
      # while the other can be around +180
      dlat  = abs(glat-nlat)
      dlon  = abs(glon-nlon)
      glon[dlon>1]   = glon[dlon>1]+360.
      dlon  = abs(glon-nlon)
      dlon0 = abs(glon[dlon<1]-nlon[dlon<1])

      # points that still have big differences in lon
      # - this is just around the north pole
      print('> "Bad" grid points:')
      print('>> Latitudes (degrees): '+str(glat[dlon>1])+' '+str(nlat[dlon>1]))
      print('>> Longitudes (degrees): '+str(glon[dlon>1])+' '+str(nlon[dlon>1]))
      print('\n')

      # rest of points - lon/lat should be close
      print('> Comparing rest of grid points...')
      print('>> Max diff in lat arrays (degrees): '+str(dlat.max()))
      print('>> Max diff in lon arrays (degrees): '+str(dlon0.max()))
      print('\n')
      glon[glon>180] = glon[glon>180]-360
      #################################################################

      #################################################################
      # compare data arrays
      print('Comparing data arrays between netcdf and grib2 files...\n')

      if KEEP_MASK:
         Data_arr = data_arr.data
         nmask    = data_arr.mask
         gdat     = msg.data()[0].data
         gmask    = msg.data()[0].mask

         # get min/max for plotting:
         clim     = [msg.data()[0].min(),
                     msg.data()[0].max()]
      else:
         Data_arr = data_arr
         nmask    = data[:,:].mask
         gdat     = msg.data()[0]
         gmask    = data[:,:].mask

         # get min/max for plotting:
         clim     = [gdat[not np.isnan(gdat)].min(),
                     gdat[not np.isnan(gdat)].max()]
      #################################################################

      #################################################################
      # print diff between arrays:
      Darr  = Data_arr[nmask==False]
      Gdat  = gdat    [gmask==False]
      ddat  = abs(Gdat-Darr)
      print('>> Max diff in data arrays: '+str(ddat.max()))
      print('\n')
      #################################################################

      #################################################################
      if PLOT_OPT>0:
         # make plots:
         figdir   = 'figs/'
         if not os.path.exists(figdir):
            os.makedirs(figdir)

         # plot titles:
         ttl   = data.standard_name
         ttl   = ttl.replace('_',' ')
 
         ######################################################################
         def do_test_plot(Z,grb=None,**kwargs):

            fig   = plt.figure()

            if grb is None:
               # simple plot with imshow (reverse i index)
               nx,ny = Z.shape
               Z2    = [Z[i,:] for i in range(nx-1,-1,-1)]
               plt.imshow(Z2,vmin=clim[0],vmax=clim[1])
            else:
               # more complicated plot with basemap (TODO fix this)
               from mpl_toolkits.basemap import Basemap

               nx = grb.points_in_x_direction
               ny = grb.points_in_y_direction
               dx = grb.gridlength_in_x_direction
               dy = grb.gridlength_in_y_direction
               #
               width    = dx*nx
               height   = dy*ny
               m        = Basemap(  projection=grb.proj4_proj,\
                                    lon_0=grb.proj4_lon_0,\
                                    lat_0=grb.proj4_lat_0,\
                                    lat_ts=grb.proj4_lat_ts,\
                                    width=width,height=height,\
                                    rsphere=grb.earthRmajor,\
                                    resolution='i',area_thresh=10000)

               # plot variable:
               m.pcolor(nlon,nlat,Z,vmin=clim[0],vmax=clim[1],latlon=True)
               m.drawcoastlines()
               m.fillcontinents()

               # draw parallels & meridians:
               parallels = np.arange(0.,90,10.)
               m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
               meridians = np.arange(-180.,180.,10.)
               m.drawmeridians(meridians,labels=[0,0,0,1],
                               fontsize=10,latmax=90)

            plt.colorbar()
            return fig
         ######################################################################

         ######################################################################
         from matplotlib import pyplot as plt 
         print('Making nc test plot...')
         Z        = Data_arr
         Z[nmask] = np.nan

         if PLOT_OPT==1:
            # simple plot with imshow (reverse i index)
            fig   = do_test_plot(Z,vmin=clim[0],vmax=clim[1])
         else:
            # make basemap figure:
            # TODO: not working at the moment
            fig   = do_test_plot(Z,grb=grb,vmin=clim[0],vmax=clim[1])

         # add title and save:
         plt.title(ttl)
         fnam  = figdir+'test_nc_'+vbl_name+'.png'
         print('Saving to: '+fnam+'\n')
         plt.savefig(fnam)
         plt.close()
         fig.clf()
         ######################################################################

         ######################################################################
         print('Making grb2 test plot...')
         Z        = gdat
         Z[gmask] = np.nan
         if PLOT_OPT==1:
            # simple plot with imshow (reverse i index)
            do_test_plot(Z,vmin=clim[0],vmax=clim[1])
         else:
            # make basemap figure:
            # TODO: not working at the moment
            do_test_plot(Z,grb=grb,vmin=clim[0],vmax=clim[1])

         # add title and save:
         plt.title(ttl)
         fnam  = figdir+'test_grb2_'+vbl_name+'.png'
         print('Saving to: '+fnam+'\n')
         plt.savefig(fnam)
         plt.close()
         fig.clf()
         #################################################################

      # finished with message
      ####################################################################

   # finished test
   gr.close()
   #######################################################################
