import os,sys
import pygrib
import numpy as np
from ncepgrib2 import Grib2Encode as g2e
from ncepgrib2 import Grib2Decode as g2d

#plotting
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
import basemap_gridlines as bmg

TEST_REWRITE   = 0
if TEST_REWRITE:
   fil_out  = 'out/test_write_multi_1.grb2'

# fil1  = 'eg_grib2/multi_1.at_10m.tp.200911.grb2'
fil1  = 'out/SWARPwavesice_forecast_start20150723T000000Z.grb2'
print('reading '+fil1+'\n')
gr       = pygrib.open(fil1)

if 1:
   N        = 1
   grbmsgs  = gr.read(N) # get list of 1st N messages:
   # grbmsgs  = gr.read()
   for msg in grbmsgs:
      print(msg)
      print('\n')
      grb   = g2d(msg.tostring(),gribmsg=True)
      print(grb)
      #
      lat,lon        = msg.latlons()
      data           = msg.values
      Z              = data.data
      Z[data.mask]   = np.NaN
      Zmin           = data.min()
      Zmax           = data.max()
      #
      lon_0 = lon.mean()
      lat_0 = lat.mean()
      R1    = grb.earthRmajor

      ###############################################################################
      if 0:
         # plot variable
         bm = Basemap(projection='stere',lon_0=lon_0,lat_0=lat_0,lat_ts=lat_0,\
                     llcrnrlat=lat.min(),urcrnrlat=lat.max(),\
                     llcrnrlon=lon.min(),urcrnrlon=lon.max(),\
                     rsphere=R1,resolution='i')
         #
         figname  = fil1.split('.grb2')[0]+'.png'
         f        = plt.figure()
         X,Y      = bm(lon,lat)
         bm.pcolor(X,Y,Z,vmin=Zmin,vmax=Zmax)
         bm.colorbar()
         bmg.latlon_grid(bm,10.,10.)
         
         nm    = str(msg.name)
         fct   = str(msg.forecastTime)
         fcu   = msg.fcstimeunits
         ttl   = nm+'. Forecast time: '+fct+fcu
         plt.title(str(msg.name))
         plt.savefig(figname)
         plt.close()
         f.clf()
      ###############################################################################

gr.close()

if not TEST_REWRITE:
   sys.exit('Not testing rewrite')

# re-write the grib message to a new file.
f_out    = open(fil_out,'wb')

# pygrib doc for
# https://github.com/erdc-cm/pygrib/blob/master/g2clib_src/grib2.h
disc_code   = grb.discipline_code # integer: cf 
isec        = grb.identification_section
   # cf http://www.cosmo-model.org/content/consortium/generalMeetings/general2014/wg6-pompa/grib2/grib/grib2keys_1.htm
   # numpy array - rank=1, len=13, dtype=int32
   # isec   = np.array([i_centre,
   #                    i_subcentre,i3,i4,
   #                    i_sig_reftime,       #significance of ref time: eg 1=start of forecast
   #                    year,month,day,
   #                    hour,minute,second,
   #                    i_ProdStatus,        #productionStatusOfProcessedData eg 0=operational product
   #                    i_DataType])         #typeOfProcessedData eg 1=forecast product

grb_out = g2e(disc_code,isec)

gdsinfo  = grb.grid_definition_info
"""
  - gdsinfo[0] = Source of grid definition (see Code
    U{Table 3.0
    <http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-0.shtml>})
  - gdsinfo[1] = Number of grid points in the defined grid.
  - gdsinfo[2] = Number of octets needed for each additional grid points defn.
    Used to define number of points in each row for non-reg grids (=0 for
    regular grid).
  - gdsinfo[3] = Interp. of list for optional points defn (Code
    U{Table 3.11
    <http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-11.shtml>})
  - gdsinfo[4] = Grid Definition Template Number (Code
    U{Table 3.1
    <http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-1.shtml>})
"""

grb_out.addgrid(gdsinfo,grb.grid_definition_template)
# add product definition template, data representation template
# and data (including bitmap which is read from data mask).
grb_out.addfield(grb.product_definition_template_number,
                 grb.product_definition_template,
                 grb.data_representation_template_number,
                 grb.data_representation_template,data)
# finalize the grib message.
grb_out.end()

# write it to the file.
f_out.write(grb_out.msg)
# close the output file
f_out.close()
