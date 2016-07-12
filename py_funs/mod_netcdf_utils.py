import numpy as np
from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen
import os,sys
import mod_reading as MR


##########################################################
def lonlat_names(ncfil):
   nc = ncopen(ncfil)
   for vbl in nc.variables:
      if 'lon' in vbl or 'Lon' in vbl:
         lon   = vbl
      if 'lat' in vbl or 'Lat' in vbl:
         lat   = vbl
   nc.close()
   return lon,lat
##########################################################


##########################################################
def nc_get_var(ncfil,vblname,time_index=None):
   """
   vbl=nc_get_var(ncfil,vblname,time_index=None)
   *ncfil is string (filename)
   *vname is string (variable name)
   *time_index is record number to get
   *vbl is a mod_reading.var_object instance
   """

   nc    = ncopen(ncfil)
   vbl0  = nc.variables[vblname]

   # get the netcdf attributes
   attlist   = vbl0.ncattrs()
   attvals  = []
   for att in attlist:
      attval   = getattr(vbl0,att)
      attvals.append(attval)

   dims  = vbl0.dimensions
   shape = vbl0.shape

   ##################################################
   # some attributes that depend on rank
   if vbl0.ndim==1:
      vals  = vbl0[:]
   elif vbl0.ndim==2:
      vals  = vbl0[:,:]
   elif vbl0.ndim==3:
      if time_index is None:
         if shape[0]==1:
            vals  = vbl0[0,:,:]
            dims  = dims[1:]
         else:
            vals  = vbl0[:,:,:]
      else:
         vals  = vbl0[time_index,:,:]
         dims  = dims[1:]
   ##################################################

   nc.close()

   attlist.append('dimensions')
   attvals.append(dims)
   vbl   = MR.var_object(vals,extra_atts=[attlist,attvals])

   return vbl
########################################################


##########################################################
def nc_get_dim(ncfil,vblname):
   """
   vbl=nc_get_var(ncfil,vblname,time_index=None)
   *ncfil is string (filename)
   *vname is string (variable name)
   *vbl   is a mod_reading.var_object instance
   """

   nc    = ncopen(ncfil)

   if vblname in nc.variables:
      vbl0  = nc.variables[vblname]

      # get the netcdf attributes
      attlist   = vbl0.ncattrs()
      attvals  = []
      for att in attlist:
         attval   = getattr(vbl0,att)
         attvals.append(attval)
      Xatts = [attlist,attvals]

      vals  = vbl0[:]
      nc.close()

      return MR.var_object(vals,extra_atts=Xatts)
   else:
      raise ValueError(vbl_name+' not given as an array')
########################################################


########################################################
class nc_getinfo:

   import os
   #####################################################
   def __init__(self,ncfil,time_index=None,lonlat_file=None):

      ##################################################
      self.filename  = ncfil
      if ncfil[0]=='/':
         bn             = os.path.basename(ncfil)
         self.basedir   = ncfil.strip(bn)
      else:
         bn             = ncfil
         self.basedir   = os.getcwd()+'/'

      self.basename     = os.path.splitext(bn)[0]
      self.filetype     = 'netcdf'
      self.object_type  = 'netcdf'

      # things to work with plotting stuff in mod_reading.py
      self.HYCOM_region    = None
      self.get_fixed_lonlat = self.get_lonlat
      ##################################################

      # added here manually
      # - TODO could possibly be determined
      #   from netcdf metadata though
      # - could also be an input
      self.reftime_sig  = 'start of forecast'


      # open the file
      nc    = ncopen(ncfil)
      dkeys = nc.dimensions.keys()
      vkeys = nc.variables.keys()

      # remove dimensions from variables
      self.dimensions   = dkeys
      for key in dkeys:
         if key in vkeys:
            vkeys.remove(key)
      Nkeys = len(vkeys)

      # is time a dimension?
      time_names     = ['time','time_counter'] # NEMO outputs call time "time_counter"
      self.time_dim  = False
      for time_name in time_names:
         if time_name in self.dimensions:
            self.time_dim  = True
            self.time_name = time_name
            break


      # get global netcdf attributes
      class ncatts:
         def __init__(self,nc):
            for att in nc.ncattrs():
               attval   = getattr(nc,att)
               setattr(self,att,attval)
            return
         def atts2list(self):
            return vars(self).keys()
      self.ncattrs   = ncatts(nc)

      ########################################################
      # time info:
      if self.time_dim:

         time        = nc.variables[self.time_name]
         Nt          = len(time[:])
         reftime_u   = time[0] # hours since refpoint
         time_info   = time.units.split()

         time_info[0]   = time_info[0].strip('s') # 1st word gives units
         if time_info[0]=='econd':
            time_info[0]   = 'second'

         tu    = time.units
         if ('T' in tu) and ('Z' in tu):
            # using the T...Z format for time
            # eg hyc2proj (this is the standard)
            split1   = tu.split('T')
            ctime    = split1[1].split('Z')[0]
            cdate    = split1[0].split()[2]
         else:
            # eg WAMNSEA product from met.no
            split1   = tu.split()
            cdate    = split1[2]
            ctime    = split1[3]

         if '-' in cdate:
            # remove '-'
            # - otherwise assume YYYYMMDD format
            split2   = cdate.split('-')
            for loop_i in range(1,3):
               if len(split2[loop_i])==1:
                  split2[loop_i] = '0'+split2[loop_i]
            cdate = split2[0]+split2[1]+split2[2] # should be YYYYMMDD now

         if len(cdate)<8:
            cdate = (8-len(cdate))*'0'+cdate

         if ':' in ctime:
            # remove ':'
            # - otherwise assume HHMMSS format
            split2   = ctime.split(':')
            for loop_i in range(0,3):
               if (split2[loop_i])==1:
                  split2[loop_i] = '0'+split2[loop_i]
            ctime = split2[0]+split2[1]+split2[2] # should be HHMMSS now

         year0    = int(cdate[:4])
         mon0     = int(cdate[4:6])
         day0     = int(cdate[6:8])
         hr0      = int(ctime[:2])
         min0     = int(ctime[2:4])
         sec0     = int(float(ctime[4:]))
         refpoint = datetime(year0,mon0,day0,hr0,min0,sec0)

         # check format of time
         i32   = np.array([0],dtype='int32')
         if type(i32[0])==type(reftime_u):
            reftime_u   = int(reftime_u)

         if time_info[0]=='second':
            self.reftime = refpoint+timedelta(seconds=reftime_u)
         elif time_info[0]=='hour':
            self.reftime = refpoint+timedelta(hours=reftime_u)
         elif time_info[0]=='day':
            self.reftime = refpoint+timedelta(reftime_u) #NB works for fraction of days also

         if time_info[0]=='second':
            # convert time units to hours for readability of the messages:
            self.timeunits  = 'hour'
            self.timevalues = [int((time[i]-time[0])/3600.) for i in range(Nt)]
         else:
            self.timeunits  = time_info[0]
            self.timevalues = [time[i]-time[0] for i in range(Nt)]

         self.number_of_time_records = Nt
         self.datetimes              = []
         for tval in self.timevalues:
            self.datetimes.append(self.timeval_to_datetime(tval))
      else:
         self.datetimes = None
      ########################################################


      ########################################################
      # grid info:
      if lonlat_file is None:
         lonlat_file = ncfil
      self.lonlat_file  = lonlat_file

      # are lon,lat dimensions?
      self.lonname,self.latname  = lonlat_names(self.lonlat_file)
      self.lonlat_dim            = (self.lonname in self.dimensions)

      ##############################################################
      # basic lon-lat info
      nc2   = ncopen(self.lonlat_file)
      lon   = nc2.variables[self.lonname]
      lat   = nc2.variables[self.latname]
      if self.lonlat_dim:
         self.lon0    = lon[0]
         self.lat0    = lat[0]

         # get example variable:
         # - for eg plotting, need to make the lon/lat matrices
         # (converted from vectors)
         # have the same shape as the variables
         vbl_dims = nc.variables[vkeys[0]].dimensions
         for dkey in vbl_dims:
            if dkey==self.lonname:
               self.lon_first = True
               self.shape     = (len(lon),len(lat))
               break
            elif dkey==self.latname:
               self.lon_first = False
               self.shape     = (len(lat),len(lon))
               break
      else:
         self.lon0    = lon[0,0]
         self.lat0    = lat[0,0]
         self.shape   = lon.shape
      nc2.close()
      ##############################################################
      

      ##############################################################
      ny,nx        = self.shape
      self.Npts_x  = nx    # No of points in x dirn
      self.Npts_y  = ny    # No of points in y dirn
      self.Npts    = nx*ny # Total no of points
      ########################################################


      ########################################################
      # projection info:
      proj_list   = ['stereographic','projection_3'] # could also have mercator or regular lon-lat
      HAVE_PROJ   = 0   # if 0 assume HYCOM native grid
      for proj_name in proj_list: 
         if proj_name in vkeys:
            proj        = nc.variables[proj_name]
            att_list    = proj.ncattrs()
            HAVE_PROJ   = 1
            break

      if HAVE_PROJ:
         # object with the netcdf attributes of projection variable
         # + some extra proj-dependent info 
         att_list_full  = [att_list[i] for i in range(len(att_list))]
         att_vals_full  = []
         for att in att_list:
            att_val  = proj.getncattr(att)
            att_vals_full.append(att_val)

         # specific to stereographic
         if proj_name=='stereographic':
            # add x,y resolution to ncinfo.proj_info
            att_list_full.extend(['x_resolution','y_resolution'])

            xx = nc.variables['x'][0:2]
            yy = nc.variables['y'][0:2]
            dx = xx[1]-xx[0]
            dy = yy[1]-yy[0]

            #convert to m
            xunits   = nc.variables['x'].units.split()
            fac      = 1.
            if len(xunits)==2:
               fac   = float(xunits[0])
               xunits.remove(xunits[0])

            if xunits[0]=='km':
               fac   = fac*1.e3
            #
            att_vals_full.extend([dx*fac,dy*fac])

         self.proj_info = MR.proj_obj(att_list_full,att_vals_full)
      else:
         self.proj_info = []
      ########################################################
      
      ########################################################
      # variable list
      # - remove some other variables from vkeys
      # - eg projection,lon,lat
      # - TODO model_depth?
      bkeys = [proj_name,self.lonname,self.latname]
      # bkeys.append('model_depth')
      for key in bkeys:
         if key in vkeys:
            vkeys.remove(key)

      self.variables       = vkeys
      self.variables3d     = None  #TODO enable treatment of 3d fields
      self.all_variables   = vkeys
      ########################################################

      nc.close()
      return
   ###########################################################

   
   ###########################################################
   def nearestDate(self, pivot):
      """
      dto,time_index = self.nearestDate(dto0)
      dto0  = datetime.datetime objects
      dto   = datetime.datetime objects - nearest value in self.datetimes to dto0
      time_index: dto=self.datetimes[time_index]
      """
      dto         = min(self.datetimes, key=lambda x: abs(x - pivot))
      time_index  = self.datetimes.index(dto)
      return dto,time_index
   ###########################################################



   ###########################################################
   def timeval_to_datetime(self,timeval):

      # check format of time
      i32   = np.array([0],dtype='int32')
      if type(i32[0])==type(timeval):
         timeval  = int(timeval)

      if self.timeunits=='second':
         dt = self.reftime +timedelta(seconds=timeval)
      elif self.timeunits=='hour':
         dt = self.reftime +timedelta(hours=timeval)
      elif self.timeunits=='day':
         dt = self.reftime +timedelta(timeval) #NB works for fraction of days also
      return dt
   ###########################################################

   ###########################################################
   def get_lonlat(self,vec2mat=True):

      nc    = ncopen(self.lonlat_file)
      lono  = nc.variables[self.lonname]
      lato  = nc.variables[self.latname]

      if lono.ndim==2:
         lon   = lono[:,:]
         lat   = lato[:,:]
      else:
         lon   = lono[:]
         lat   = lato[:]
         if vec2mat:
            if self.lon_first:
               # lon in cols, lat in rows
               lon,lat  = np.meshgrid(lon,lat,indexing='ij')
            else:
               # lon in rows, lat in cols
               lon,lat  = np.meshgrid(lon,lat,indexing='xy')
      nc.close()

      return lon,lat
   ###########################################################


   ###########################################################
   def get_var(self,vname,time_index=None):
      """
      Call: self.get_var(vname,time_index=None)
      Inputs:
      vname = string (name of variable)
      time_index = integer: if time_index 
      Returns: mod_reading.var_object instance
      """

      vname = MR.check_names(vname,self.variables)
      vbl   = nc_get_var(self.filename,vname,time_index=time_index)

      return vbl
   ###########################################################


   ###########################################################
   def get_dim(self,dname):
      """
      call: self.get_dim(dname)
      inputs:
      dname = string (name of dimension)
      returns: mod_reading.var_object instance
      """

      vname = MR.check_names(dname,self.dimensions)
      vbl   = nc_get_dim(self.filename,dname)

      return vbl
   ###########################################################


   #######################################################################
   def imshow(self,var_opts,**kwargs):
      """
      pobj   = self.imshow(var_opts,time_index=0,pobj=None,\
           clim=None,add_cbar=True,clabel=None,show=True,\
           test_ijs=None)
      """

      return MR.imshow(self,var_opts,**kwargs)
   #######################################################################


   #######################################################################
   def plot_var(self,var_opts,**kwargs):
      """
      pobj,bmap = self.plot_var(var_opts,time_index=0,\
         pobj=None,bmap=None,HYCOMreg='TP4',\
         clim=None,add_cbar=True,clabel=None,show=True,\
         test_lonlats=None,date_label=0):
      """

      return MR.plot_var(self,var_opts,**kwargs)
   #######################################################################


   #######################################################################
   def plot_var_pair(self,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs):
      """
      pobj,bmap=self.plot_var_pair(var_opts1,var_opts2,pobj=None,bmap=None,**kwargs)
      """
      return MR.plot_var_pair(self,var_opts1,var_opts2,**kwargs)
   #######################################################################


   ###########################################################
   def make_png(self,var_opts,**kwargs):
      """
      pobj,bmap=self.make_png(var_opts,pobj=None,bmap=None,figdir='.',time_index=0,date_label=2,**kwargs)
      """
      return MR.make_png(self,var_opts,**kwargs)
   ###########################################################


   ###########################################################
   def make_png_pair(self,var_opts1,var_opts2,**kwargs):
      """
      pobj,bmap = self.make_png_pair(var_opts1,var_opts2,\
         pobj=None,bmap=None,figdir='.',date_label=2,**kwargs)
      """
      return MR.make_png_pair(self,var_opts1,var_opts2,**kwargs)
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs(self,**kwargs):
      """
      pobj,bmap,obsfil = self.compare_ice_edge_obs(pobj=None,bmap=None,time_index=0,\
         obs_type='OSISAF',date_label=1,figname=None,**kwargs)
      """
      return MR.compare_ice_edge_obs(self,**kwargs)
   ###########################################################


   ###########################################################
   def MIZmap(self,**kwargs):
      """
      Call  : self.MIZmap(var_name='dmax',do_sort=False,EastOnly=True,plotting=True,**kwargs)
      Inputs:
         var_name is variable to find MIZ from
         **kwargs to be passed onto MIZchar.get_MIZ_poly:
            outdir='.',do_sort=True
      Returns: MIZchar.MIZpoly object
      """
      return MR.MIZmap(self,**kwargs)
   ###########################################################


   ###########################################################
   def areas_of_disagreement(self,**kwargs):
      """
      MPdict,tfiles,Pdict = self.areas_of_disagreement(obs_type='OSISAF',\
            time_index=0,do_sort=True,EastOnly=True,\
            plotting=True,HYCOMreg='Arctic',**kwargs)
      kwargs: outdir='.',do_sort=True
      """
      return MR.areas_of_disagreement(self,**kwargs)
   ###########################################################


   ###########################################################
   def make_png_all(self,var_opts,**kwargs):
      """
      self.make_png_all(var_opts,HYCOMreg=None,figdir='.')
      """

      MR.make_png_all(self,var_opts,**kwargs)
      return
   ###########################################################


   ###########################################################
   def make_png_pair_all(self,var_opts1,var_opts2,**kwargs):
      """
      self.make_png_pair_all(var_opts1,var_opts2,HYCOMreg=None,figdir='.')
      """

      MR.make_png_pair_all(self,var_opts1,var_opts2,**kwargs)
      return
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs_all(self,**kwargs):
      out   = MR.compare_ice_edge_obs(self,**kwargs)
      return out
   ###########################################################

###########################################################
