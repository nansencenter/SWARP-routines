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
            setattr(self,att_names[n],0)
##########################################################

# ##########################################################
# def read_proj_infile(pfil):
#    # extracts information from proj.in file used by hyc2proj
#    # when producing netcdf files from binaries
# 
# 
#    pp = open(pfil,'r')
#    pl = pp.readlines()
#    pp.close()
# 
#    proj_name   = pl[0].strip()
#    Nattr       = len(pl)-1
#    att_names   = []
#    att_vals    = []
# 
#    ############################################
#    # make lists of attribute names and values
#    for n in range(1,Nattr+1):
#       pln      = pl[n].split('#')
#       attval   = pln[0]
#       attname  = pln[1].split()
# 
#       # attval can be numeric or logical
#       if 'T' in attval:
#          attval   = True
#       elif 'F' in attval:
#          attval   = False
#       else:
#          attval   = float(attval)
# 
#       # assign names:
#       ss       = attname
#       attname  = ss[0]
#       for j in range(1,len(ss)):
#          attname  = attname+'_'+ss[j]
# 
#       # remove troublesome characters:
#       attname  = attname.replace('(','')
#       attname  = attname.replace(')','')
#       attname  = attname.replace('>','')
#       attname  = attname.replace('-','_')
#       if '=' in attname:
#          jn = attname.index('=')
#          attname  = attname[:jn-2]
# 
#       # append to lists:
#       # print(attname)
#       # print(attval)
#       att_names.append(attname)
#       att_vals.append(attval)
#    ############################################
# 
#    proj_in  = proj_obj(proj_name,att_names,att_vals)
#    return proj_in
# ##########################################################

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
   reftime_h   = time[0] # hours since refpoint
   time_info   = time.units.split()
   time_info   = [ time_info[i] for i in [0,2] ]

   time_fmt       = '%Y-%m-%dT%H:%M:%SZ' # eg 1950-1-1T12:00:00Z
   refpoint       = datetime.strptime(time_info[1],time_fmt)
   ncinfo.reftime = refpoint+timedelta(hours=reftime_h)

   ncinfo.timeunits  = time_info[0]
   if time_index is not None:
      ncinfo.datatime   = time[time_index]-reftime_h

   ncinfo.number_of_time_records = len(time[:])
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
   proj_list   = ['stereographic'] # could also have mercator or regular lon-lat
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
