import numpy as np
from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen

##########################################################
def read_proj_infile(pfil):
   # extracts information from proj.in file used by hyc2proj
   # when producing netcdf files from binaries

   ############################################
   class proj_obj:
      def __init__(self,proj_name,att_names,att_vals):

         # projection name:
         setattr(self,'projection_name',proj_name)

         # rest of attributes:
         for n in range(len(att_names)):
            setattr(self,att_names[n],att_vals[n])
   ############################################

   pp = open(pfil,'r')
   pl = pp.readlines()
   pp.close()

   proj_name   = pl[0].strip()
   Nattr       = len(pl)-1
   att_names   = []
   att_vals    = []

   ############################################
   # make lists of attribute names and values
   for n in range(1,Nattr+1):
      pln      = pl[n].split('#')
      attval   = pln[0]
      attname  = pln[1].split()

      # attval can be numeric or logical
      if 'T' in attval:
         attval   = True
      elif 'F' in attval:
         attval   = False
      else:
         attval   = float(attval)

      # assign names:
      ss       = attname
      attname  = ss[0]
      for j in range(1,len(ss)):
         attname  = attname+'_'+ss[j]

      # remove troublesome characters:
      attname  = attname.replace('(','')
      attname  = attname.replace(')','')
      attname  = attname.replace('>','')
      attname  = attname.replace('-','_')
      if '=' in attname:
         jn = attname.index('=')
         attname  = attname[:jn-2]

      # append to lists:
      # print(attname)
      # print(attval)
      att_names.append(attname)
      att_vals.append(attval)
   ############################################

   proj_in  = proj_obj(proj_name,att_names,att_vals)
   return proj_in
##########################################################

##########################################################
def nc_get_var(ncfil,vblname,time_index=None):
   #get basic info from 2d variable

   ########################################################
   class def_vbl:

      #####################################################
      def __init__(self,vbl0,time_index=None):
         if hasattr(vbl0,'name'):
            self.name   = vbl0.name
         if hasattr(vbl0,'units'):
            self.units  = vbl0.units
         if hasattr(vbl0,'standard_name'):
            self.standard_name   = vbl0.standard_name
         if hasattr(vbl0,'long_name'):
            self.long_name = vbl0.long_name
         if hasattr(vbl0,'add_offset'):
            self.add_offset   = vbl0.add_offset
         if hasattr(vbl0,'scale_factor'):
            self.scale_factor   = vbl0.scale_factor
         if hasattr(vbl0,'_FillValue'):
            self._FillValue   = vbl0._FillValue
         if hasattr(vbl0,'missing_value'):
            self.missing_value   = vbl0.missing_value
         if hasattr(vbl0,'cell_methods'):
            self.cell_methods = vbl0.cell_methods
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
                  self.dimensions = dims[1:]
         ##################################################

         self.values = vals
         self.shape  = vals.shape
         self.ndim   = vals.ndim
         self.__getitem__  = self.values.__getitem__
      #####################################################

   ########################################################

   nc    = ncopen(ncfil)
   vbl0  = nc.variables[vblname]
   vbl   = def_vbl(vbl0,time_index)
   nc.close()

   return vbl

def nc_getinfo(ncfil,time_index=None):

   ########################################################
   class nc_info:

      #####################################################
      def __init__(self):
         self.lon0         = 0 # 1st longitude value
         self.lat0         = 0 # 1st latitude value
         self.shape        = 0 # shape of grid

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
   ########################################################

   nc.close()
   return ncinfo
