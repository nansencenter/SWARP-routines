import numpy as np
from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen

##########################################################
class proj_obj:
   # projection info in this object
   def __init__(self,att_names,att_vals=None):

      # rest of attributes:
      if att_vals is not None:
         for n in range(len(att_names)):
            setattr(self,att_names[n],att_vals[n])
      else:
         for n in range(len(att_names)):
            setattr(self,att_names[n],0)
##########################################################

##########################################################
def nc_get_var(ncfil,vblname,time_index=None):
   #get basic info from 2d variable

   ########################################################
   class def_vbl:

      #####################################################
      def __init__(self,vbl0,time_index=None):

         # get the netcdf attributes
         ncatts   = vbl0.ncattrs()
         for att in ncatts:
            attval   = getattr(vbl0,att)
            setattr(self,att,attval)

         if hasattr(vbl0,'dimensions'):
            dims              = vbl0.dimensions
            self.dimensions   = dims

         ##################################################
         # some attributes that depend on rank
         if vbl0.ndim==1:
            vals  = vbl0[:]
         elif vbl0.ndim==2:
            vals  = vbl0[:,:]
         elif vbl0.ndim==3:
            if time_index is None:
               vals  = vbl0[:,:,:]
            else:
               vals  = vbl0[time_index,:,:]
               if hasattr(vbl0,'dimensions'):
                  self.dimensions = dims[1:]# drop time dimension
         ##################################################

         self.values = vals
         self.shape  = vals.shape
         self.ndim   = vals.ndim
         self.__getitem__  = self.values.__getitem__
      #####################################################

   ########################################################

   nc       = ncopen(ncfil)
   vbl0     = nc.variables[vblname]
   vbl      = def_vbl(vbl0,time_index)
   nc.close()

   return vbl

def nc_getinfo(ncfil,time_index=None):

   ########################################################
   class nc_info:

      #####################################################
      def __init__(self):
         self.ncattrs         = 0 # global netcdf attributes
         self.proj_info       = 0 # info about projection
         self.variable_list   = 0 # list of variable names

         self.lon0         = 0 # 1st longitude value
         self.lat0         = 0 # 1st latitude value
         self.shape        = 0 # shape of grid
         self.Npts_x       = 0 # No of points in x dirn
         self.Npts_y       = 0 # No of points in y dirn
         self.Npts         = 0 # Total no of points

         self.reftime      = 0
         # added here manually
         # - could possibly be determined
         #   from netcdf metadata though
         self.reftime_sig  = 'start of forecast'

         self.timevalues = []
         self.timeunits = 'hour'
         self.datatime  = 0

         self.number_of_time_records   = 1
   ########################################################

   ncinfo   = nc_info()
   nc       = ncopen(ncfil)

   # get global netcdf attributes
   ncinfo.ncattrs  = nc.ncattrs()

   dkeys = nc.dimensions.keys()
   vkeys = nc.variables.keys()
   Nkeys = len(vkeys)

   ########################################################
   # time info:
   time        = nc.variables['time']
   Nt          = len(time[:])
   reftime_u   = time[0] # hours since refpoint
   time_info   = time.units.split()
   #
   time_info[0]   = time_info[0].strip('s') # 1st word gives units
   if time_info[0]=='econd':
      time_info[0]   = 'second'
   #
   tu    = time.units
   i0    = tu.index('-')
   date0 = tu[i0-4:i0+6]
   i1    = tu.index(':')
   date1 = tu[i1-2:i1+6]
   #
   date2    = date0+'T'+date1+'Z'
   time_fmt = '%Y-%m-%dT%H:%M:%SZ' # eg 1950-1-1T12:00:00Z
   refpoint = datetime.strptime(date2,time_fmt)
   #
   if time_info[0]=='second':
      ncinfo.reftime = refpoint+timedelta(seconds=reftime_u)
   elif time_info[0]=='hour':
      ncinfo.reftime = refpoint+timedelta(hours=reftime_u)

   if time_info[0]=='hour':
      ncinfo.timeunits  = time_info[0]
      ncinfo.timevalues = [time[i]-time[0] for i in range(Nt)]
   elif time_info[0]=='second':
      # convert time units to hours for readability of the messages:
      ncinfo.timeunits  = 'hour'
      ncinfo.timevalues = [int((time[i]-time[0])/3600.) for i in range(Nt)]

   if time_index is not None:
      ncinfo.datatime   = ncinfo.timevalues[time_index]

   ncinfo.number_of_time_records = Nt

   ########################################################

   ########################################################
   # grid info:
   ncinfo.lon0    = nc.variables['longitude'][0,0]
   ncinfo.lat0    = nc.variables['latitude'] [0,0]
   ncinfo.shape   = nc.variables['latitude'] [:,:].shape
   ny,nx          = ncinfo.shape
   ncinfo.Npts_x  = nx # No of points in x dirn
   ncinfo.Npts_y  = ny # No of points in y dirn
   ncinfo.Npts    = nx*ny # Total no of points

   ncinfo.shape   = nc.variables['latitude'] [:,:].shape
   ########################################################

   ########################################################
   # projection info:
   proj_list   = ['stereographic','projection_3'] # could also have mercator or regular lon-lat
   for proj_name in proj_list: 
      if proj_name in vkeys:
         proj     = nc.variables[proj_name]
         att_list = proj.ncattrs()
         break

   # object with the netcdf attributes of projection variable
   # + some extra proj-dependent info 
   att_list_full  = [att_list[i] for i in range(len(att_list))]
   if proj_name=='stereographic':
      att_list_full.extend(['x_resolution','y_resolution'])

   ncinfo.proj_info  = proj_obj(att_list_full)
   for att in att_list:
      attval   = proj.getncattr(att)
      setattr(ncinfo.proj_info,att,attval)

   if proj_name=='stereographic':
      # add x,y resolution to ncinfo.proj_info
      xx    = nc.variables['x'] [0:2]
      yy    = nc.variables['y'] [0:2]
      dx    = xx[1]-xx[0] # units of 100km
      dy    = yy[1]-yy[0] # units of 100km
      fac   = 1e5 # convert from 100km to m
      #
      ncinfo.proj_info.x_resolution = dx*fac
      ncinfo.proj_info.y_resolution = dy*fac
   ########################################################
   
   ########################################################
   # variable list
   vlist = []
   bkeys = [proj_name,'longitude','latitude','model_depth']
   bkeys.extend(dkeys)
   for key in vkeys:
      if key not in bkeys:
         vlist.append(key)

   ncinfo.variable_list = vlist
   ########################################################

   nc.close()
   return ncinfo
###########################################################

###########################################################
def print_grib_messages(grb2fil,N=None):

   import pygrib
   from ncepgrib2 import Grib2Encode as g2e
   from ncepgrib2 import Grib2Decode as g2d

   gr       = pygrib.open(grb2fil)

   if N is None:
      # print all messages:
      grbmsgs  = gr.read()
   else:
      # print 1st N messages:
      grbmsgs  = gr.read(N)

   for msg in grbmsgs:
      print(msg)
      print('\n')
      grb   = g2d(msg.tostring(),gribmsg=True)
      print(grb)
      print('\n')

   gr.close()
###########################################################

##############################################################
def get_array_from_binary(fid,nx,ny,fmt_size=4,order='fortran'):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)
   recs     = nx*ny
   rec_size = recs*fmt_size
   #
   data  = fid.read(rec_size)
   if fmt_size==4:
      fmt_py   = 'f' # python string for single
   else:
      fmt_py   = 'd' # python string for double

   fld   = struct.unpack(recs*fmt_py,data)
   fld   = np.array(fld)

   if order!='fortran':
      fld   = fld.reshape((nx,ny))  # array order follows python/c convention
                                    # (index increases across array)
   else:
      fld   = fld.reshape((ny,nx)).transpose()  # need to transpose because of differences between
                                                # python/c and fortran/matlab 

   return fld
##############################################################

##############################################################
def get_array_from_HYCOM_binary(afile,recno,dims=None,afile_grid='regional.grid.a'):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)

   ######################################################################
   # get size of grid
   if dims is None:
      # check regional.grid.b file for size of grid
      bfile = afile_grid[:-1]+'.b'
      bid   = open(bfile,'r')
      nx    = int( bid.readline().split()[0] )
      ny    = int( bid.readline().split()[0] )
      bid.close()
   else:
      nx = dims[0]
      ny = dims[1]
   ######################################################################

   ######################################################################
   # set record size, skip to record number
   fmt_size = 8      # HYCOM files are double precision
   if fmt_size==4:
      fmt_py   = 'f' # python string for single
   else:
      fmt_py   = 'd' # python string for double

   n0       = 4096   # HYCOM stores records in multiples of 4096
   Nhyc     = (1+(nx*ny)/n0)*n0
   rec_size = Nhyc*fmt_size
   #
   aid   = open(afile,'rb')
   for n in range(1,recno):
      aid.seek(rec_size)
   ######################################################################

   # read data and close file
   data  = aid.read(rec_size)
   aid.close()

   # rearrange into correctly sized array
   fld   = struct.unpack(recs*fmt_py,data)
   fld   = np.array(fld[0:nx*ny]) # select the 1st nx,ny - rest of the Nhyc record is rubbish
   fld   = fld.reshape((ny,nx)).transpose()  # need to transpose because of differences between
                                             # python/c and fortran/matlab 

   return fld
##############################################################

