import os,sys
import numpy as np
from datetime import datetime
#
import pygrib
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d
#
import mod_reading as tw_rdg
import mod_grib2_setup as tw_g2s
ncgv  = tw_rdg.nc_get_var

##########################################################
# inputs - loop over these?

# file inputs:
pfil  = 'inputs/proj.in' # proj.in file used by hyc2proj
ncfil = "test_ncfiles/TP4DAILY_start20120723_dump20120723.nc" # hyc2proj netcdf

# output file:
outdir  = 'out'
if not os.path.exists(outdir):
   os.makedirs(outdir)
fil_out = outdir+'/test_hycproj_to_grib2.grb2'
##########################################################

##########################################################
# # read proj.in file
# proj_in  = tw_rdg.read_proj_infile(pfil)

# read hyc2proj netcdf:
ncinfo         = tw_rdg.nc_getinfo(ncfil)
time_indices   = range(ncinfo.number_of_time_records)
vbl_list       = ncinfo.variable_list

# Open outfile before starting to encode messages
f_out = open(fil_out,'wb')

time_index  = time_indices[0]
vbl_name    = vbl_list[0]
data        = ncgv(ncfil,vbl_name,time_index)

##########################################################################################################
# INITIALISATION

# INPUT 1) discipline code (pygrib/ncepgrib2.py, grib2message)
discipline_code   = 10  #Oceanographic Products > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table0-0.shtml

# INPUT 2) identification section: pygrib/g2clib.c ("unpack1" method)
ident_sect  = tw_g2s.def_ident_sect(ncinfo.reftime,ncinfo.reftime_sig)

# CALL Grib2Encode:
grb_out     = g2e(discipline_code,ident_sect)
##########################################################################################################

##########################################################################################################
# DEFINE GRID:

# INPUT 2)
# grid definition template:
# > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-20.shtml
gdtnum,gdt  = tw_g2s.def_grid_template(ncinfo)

# INPUT 1)
# grid definition section:
gdsinfo  = tw_g2s.def_gdsinfo(ncinfo,gdtnum)

# CALL addgrid
grb_out.addgrid(gdsinfo,gdt)
# sys.exit('L132')
##########################################################################################################

##########################################################################################################
# add product definition template, data representation template
# and data (including bitmap which is read from data mask).

# INPUT 1):
pdtnum   = 0 # product_definition_template_number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-0.shtml

# INPUT 2):
pdtmpl   = tw_g2s.def_prod_template(ncinfo,vbl_name)

# INPUT 3):
# data representation template number > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
# drtnum  = 40 # grid point data, jpeg 2000 compression
drtnum   = 0 # grid point data, simple packing

# INPUT 4):
# data representation template > http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table5-0.shtml
drtmpl   = tw_g2s.def_datarep_template(drtnum)

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

print('Writing grib message to '+fil_out)
f_out.write(grb_out.msg)
##########################################################################################################

f_out.close()

if 1:
   print('Test write to '+fil_out)
   gr = pygrib.open(fil_out)

   N        = 1
   grbmsgs  = gr.read(N) # get list of 1st N messages:
   for msg in grbmsgs:
      print(msg)
      print('\n')
      grb   = g2d(msg.tostring(),gribmsg=True)
      print(grb)

      if 1:
         # lon,lat from grib2
         glat,glon   = msg.latlons()

         # lat,lon from nc file
         lon   = ncgv(ncfil,'longitude')
         lat   = ncgv(ncfil,'latitude')
         nlon  = lon[:,:]
         nlat  = lat[:,:]

         if 1:
            # calc difference between arrays:
            dlat  = abs(glat-nlat)
            dlon  = abs(glon-nlon)
            glon[dlon>1]   = glon[dlon>1]+360.
            dlon  = abs(glon-nlon)
            dlon0 = abs(glon[dlon<1]-nlon[dlon<1])
            #
            print('"bad" grid points:')
            print(glon[dlon>1],nlon[dlon>1])
            print(glat[dlon>1],nlat[dlon>1])
            print('\n')
            #
            print('compare rest of grid points:')
            print('max diff in lat arrays (degrees): '+str(dlat.max()))
            print('max diff in lon arrays (degrees): '+str(dlon0.max()))
            print('\n')
            glon[glon>180] = glon[glon>180]-360

         if 1:
            # compare data arrays:
            if KEEP_MASK:
               Data_arr = data_arr.data
               nmask    = data_arr.mask
               gdat     = msg.data()[0].data
               gmask    = msg.data()[0].mask
            else:
               Data_arr = data_arr
               nmask    = data[:,:].mask
               gdat     = msg.data()[0]
               gmask    = data[:,:].mask

            if 1:
               # print diff between arrays:
               Darr  = Data_arr[nmask==False]
               Gdat  = gdat    [gmask==False]
               ddat  = abs(Gdat-Darr)
               print('max diff in data arrays: '+str(ddat.max()))
               print('\n')

            if 1:
               # make plots:
               figdir   = 'figs/'
               if not os.path.exists(figdir):
                  os.makedirs(figdir)

               # plot titles:
               ttl   = data.standard_name
               ttl   = ttl.replace('_',' ')
               clim  = [0,1]

               # make basemap:
               from matplotlib import pyplot as plt 
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

               if 1:
                  print('Making nc test plot...')
                  Z        = Data_arr
                  Z[nmask] = np.nan

                  if 1:
                     # simple plot with imshow (reverse i index)
                     Z2    = [Z[i,:] for i in range(nx-1,-1,-1)]
                     fig   = plt.figure()
                     plt.imshow(Z2,vmin=clim[0],vmax=clim[1])
                  elif 0:
                     # use basemap...
                     # plot is empty for some reason?????
                     fig   = plt.figure()
                     m.pcolor(nlon,nlat,Z,vmin=clim[0],vmax=clim[1])
                     m.drawcoastlines()
                     m.fillcontinents()

                     # draw parallels.
                     parallels = np.arange(0.,90,10.)
                     m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
                     meridians = np.arange(-180.,180.,10.)
                     m.drawmeridians(meridians,labels=[0,0,0,1],
                                     fontsize=10,latmax=90)

                  plt.colorbar()
                  plt.title(ttl)

                  fnam  = figdir+'test_nc.png'
                  print('Saving to: '+fnam+'\n')
                  plt.savefig(fnam)
                  plt.close()
                  fig.clf()

               if 1:
                  print('Making grb2 test plot...')
                  Z        = gdat
                  Z[gmask] = np.nan
                  if 1:
                     # simple plot with imshow (reverse i index)
                     Z2 = [Z[i,:] for i in range(nx-1,-1,-1)]
                     fig   = plt.figure()
                     plt.imshow(Z2,vmin=clim[0],vmax=clim[1])
                  elif 0:
                     # use basemap...
                     # plot is empty for some reason?????
                     fig   = plt.figure()
                     m.pcolor(nlon,nlat,Z,vmin=clim[0],vmax=clim[1])
                     m.drawcoastlines()
                     m.fillcontinents()

                     # draw parallels.
                     parallels = np.arange(0.,90,10.)
                     m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
                     meridians = np.arange(-180.,180.,10.)
                     m.drawmeridians(meridians,labels=[0,0,0,1],
                                     fontsize=10,latmax=90)

                  plt.colorbar()
                  plt.title(ttl)

                  fnam  = figdir+'test_grb2.png'
                  print('Saving to: '+fnam+'\n')
                  plt.savefig(fnam)
                  plt.close()
                  fig.clf()

   gr.close()
