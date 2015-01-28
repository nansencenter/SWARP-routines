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
ncgv  = m_rdg.nc_get_var

##########################################################
# inputs - loop over these?
DO_TEST  = 1

# file inputs:
pfil  = 'inputs/proj.in' # proj.in file used by hyc2proj
ncfil = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc" # hyc2proj netcdf

# output file:
outdir  = 'out'
if not os.path.exists(outdir):
   os.makedirs(outdir)
fil_out = outdir+'/test_hyc2proj_to_grib2.grb2'
##########################################################

##########################################################
# # read proj.in file
# proj_in  = m_rdg.read_proj_infile(pfil)

# read hyc2proj netcdf:
ncinfo         = m_rdg.nc_getinfo(ncfil)
time_indices   = range(ncinfo.number_of_time_records)
vbl_list       = ncinfo.variable_list

# Open outfile before starting to encode messages
f_out = open(fil_out,'wb')

################################################################################################################
if DO_TEST==1:
   # just encode and check one variable and one time record:
   time_indices   = [time_indices[0]]
   vbl_list       = [vbl_list[0]]
   #
   print("Doing test...")
   print("Encoding variable "+vbl_list[0])
   print("at time record number "+str(time_indices[0]))
   print(' ')
################################################################################################################

################################################################################################################
for time_index in time_indices:
   for vbl_name in vbl_list:
      # encode all variables and all times as one message

      data  = ncgv(ncfil,vbl_name,time_index)

      ##########################################################################################################
      # INITIALISATION

      # INPUT 1) discipline code (pygrib/ncepgrib2.py, grib2message)
      discipline_code   = 10  #Oceanographic Products > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table0-0.shtml

      # INPUT 2) identification section: pygrib/g2clib.c ("unpack1" method)
      ident_sect  = m_g2s.def_ident_sect(ncinfo.reftime,ncinfo.reftime_sig)

      # CALL Grib2Encode:
      grb_out  = g2e(discipline_code,ident_sect)
      ##########################################################################################################

      ##########################################################################################################
      # DEFINE GRID:

      # INPUT 2)
      # grid definition template:
      # > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
      gdtnum,gdt  = m_g2s.def_grid_template(ncinfo)

      # INPUT 1)
      # grid definition section:
      gdsinfo  = m_g2s.def_gdsinfo(ncinfo,gdtnum)

      # CALL addgrid
      grb_out.addgrid(gdsinfo,gdt)
      ##########################################################################################################

      ##########################################################################################################
      # add product definition template, data representation template
      # and data (including bitmap which is read from data mask).

      # INPUT 1):
      pdtnum   = 0 # product_definition_template_number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-0.shtml

      # INPUT 2):
      pdtmpl   = m_g2s.def_prod_template(ncinfo,vbl_name)

      # INPUT 3):
      # data representation template number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
      # drtnum  = 40 # grid point data, jpeg 2000 compression
      drtnum   = 0 # grid point data, simple packing

      # INPUT 4):
      # data representation template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
      drtmpl   = m_g2s.def_datarep_template(drtnum)

      # INPUT 5)
      KEEP_MASK   = 1
      if KEEP_MASK:
         data_arr = data[:,:] # masked array
      else:
         data_arr = data[:,:].data # array

      # CALL addfield:
      grb_out.addfield(pdtnum,pdtmpl,drtnum,drtmpl,data_arr)

      # finalize the grib message.
      grb_out.end()

      print('Writing grib message to '+fil_out+'\n')
      f_out.write(grb_out.msg)
      ##########################################################################################################

f_out.close()

if DO_TEST:
   print('Test write to '+fil_out)

   gr       = pygrib.open(fil_out)
   N        = 1
   grbmsgs  = gr.read(N) # get list of 1st N messages:

   for msg in grbmsgs:
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
      if 1:
         # make plots:
         figdir   = 'figs/'
         if not os.path.exists(figdir):
            os.makedirs(figdir)

         # plot titles:
         ttl   = data.standard_name
         ttl   = ttl.replace('_',' ')
 
         PLOT_OPT = 1 # 1: simple plot with imshow
                      # else: more complicated plot with basemap (TODO fix this)

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
               m.pcolor(nlon,nlat,Z,vmin=clim[0],vmax=clim[1])
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
