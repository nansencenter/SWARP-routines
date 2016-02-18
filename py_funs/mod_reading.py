import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen
import fns_plotting as Fplt
from scipy.interpolate import griddata as grd

##########################################################
def basemap_OSISAF():
    bmap = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
               projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    return bmap
##########################################################


##########################################################
class read_MIZpoly_summary:

   #######################################################
   def __init__(self,tfil,cdate=None,ctime=None):
      

      self.filename     = tfil
      self.datetime     = None
      self.time_in_days = None
      self.datetime_ref = datetime(1901,1,1)

      if (date is not None) and (time is not None):
         self.datetime  = datetime.strptime(cdate+ctime,'%Y%m%d%H%M%S')
      elif (date is not None):
         self.datetime  = datetime.strptime(cdate,'%Y%m%d')

      fid   = open(tfil,'r')
      lins  = fid.readlines()
      fid.close()

      for lin in lins:
         sp    = lin.split(':')
         att   = sp[0].strip()
         if att=='Date':
            val   = sp[1].strip()
            if len(val)==8:
               dtm   = datetime.strptime(val,'%Y%m%d')
            elif len(val)==16:
               dtm   = datetime.strptime(val,'%Y%m%dT%H%M%SZ')
            else:
               raise ValueError('date in '+tfil+' has wrong format')

            if self.datetime is not None:
               if self.datetime != dtm:
                  raise ValueError('date in '+tfil+' not consistent with inputs cdate,ctime')
            else:
               self.datetime  = dtm
         else:
            val   = float(sp[1].strip())
            setattr(self,att,val)

      if self.datetime is not None:
         dt                = self.datetime-self.datetime_ref
         self.time_in_days = dt.total_seconds()/float(24*60*60)

      return
   #######################################################

##########################################################


############################################################################
# Function that reprojects model into observational grid
def reproj_mod2obs(X1,Y1,Z1,X2,Y2,mask=None):

   # getting ready for reprojection
   X1vec = X1.reshape(X1.size)
   Y1vec = Y1.reshape(Y1.size)
   Z1vec = Z1.reshape(Z1.size)
   C = [X1vec,Y1vec]
   C = np.array(C)
   C = C.T # input coords from X1,Y1; Z1 is array to interp; X2,Y2 are output matrices

   # Interpolation can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
   Z2   = grd(C,Z1vec,(X2,Y2),method='linear')
   if mask is not None:
      # get good values from Z2
      Arr         = np.zeros(Z2.shape)
      good        = np.isfinite(Z2)
      nans        = np.logical_not(good)
      Arr[good]   = Z2[good]

      # apply union of mask and model nans
      Z2 = np.ma.array(Arr,mask=np.logical_or(mask,nans))
    
   return(Z2)
##########################################################


##########################################################
class plot_object:
   """
   create object with:
   pobj  = plot_object(fig=None,ax=None,cbar=None,axpos=None)
   fig  is a pyplot.figure instance
   ax   is a subplot axis of fig
   cbar is a colorbar associated with fig and a plot on ax
   - used in the plotting routines of mod_reading
   TODO move to fns_plotting
   """

   def __init__(self,fig=None,ax=None,cbar=None,axpos=None):

      if fig is None:
         self.fig   = plt.figure()
      else:
         self.fig   = fig

      if ax is None:
         self.ax = self.fig.add_subplot(1,1,1)
      else:
         self.ax  = ax

      self.cbar   = cbar
      self.axpos  = axpos

      return

   def get(self):
      return self.fig,self.ax,self.cbar

   def renew(self,axpos=None):
      
      pobj  = plot_object(fig=self.fig,ax=self.ax,cbar=self.cbar,axpos=axpos)

      return pobj
##########################################################


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
      return
##########################################################


###########################################################
def check_names(vname,variables):

   if vname in variables:
      return vname

   lists = []

   # ice conc alt names
   lists.append(['ficem','fice','ice_conc','icec,'\
                  'concentration','sea_ice_concentration'])

   # ice thick alt names
   lists.append(['hicem','hice','ice_thick','icetk',\
                  'sea_ice_thickness','thickness','sea_ice_concentration'])

   # floe size alt names
   lists.append(['dfloe','dmax'])

   for names in lists:
      if vname in names:
         for vbl in names:
            if vbl in variables:
               return vbl

   raise ValueError(vname+'not in variable list')
   return
###########################################################


###########################################################
def check_var_opts(var_opts,variables):
   """
   var_opts = check_var_opts(var_opts,variables)
   *var_opts can be a string with variable name
    or a mod_reading.make_plot_options object
   *variables is a list of the variables in a file 
   - error raised if the variable name is not in this list
   (there is also a list of synonyms
   eg 'hice'='sea_ice_thickness'='icetk')
   """

   ###########################################################
   class new_var_opts:
      def __init__(self,var_opts,vname):
         atts  = vars(var_opts)
         for att in atts.keys():
            setattr(self,att,atts[att])

         self.name   = vname
         return
   ###########################################################

   if type(var_opts)==type('hi'):
      # if only a string is passed in
      print("Converting string '"+var_opts+"' to mod_reading.make_plot_options object")
      print("- create and pass in such an object directly")
      print("to specify more complicated plot options:")
      print("ie")
      print("var_opts=mod_reading.make_plot_options('"+var_opts+"',")
      print("   vec_opt=0,layer=0,conv_fac=1,wave_mask=False"\
               +",ice_mask=False,dir_from=True)")
      var_opts = make_plot_options(var_opts)
   elif type(var_opts)==type([]):
      # if only a list [vname,layer] is passed in
      vname,layer = var_opts
      print("Converting string '"+vname+"'(layer="+str(layer)+\
            ") to mod_reading.make_plot_options object")
      print("- create and pass in such an object directly")
      print("to specify more complicated plot options:")
      print("ie")
      print("var_opts=mod_reading.make_plot_options('"+vname+"',")
      print("   vec_opt=0,layer="+str(layer)+",conv_fac=1,wave_mask=False"\
               +",ice_mask=False,dir_from=True)")
      var_opts = make_plot_options(vname,layer=layer)

   vname       = check_names(var_opts.name,variables)
   var_opts2   = new_var_opts(var_opts,vname)

   return var_opts2
###########################################################


###########################################################
class make_plot_options:
   """
   var_opts=make_plot_options(vname,layer=0,vec_opt=0,conv_fac=1,\
      wave_mask=False,ice_mask=False,dir_from=True)
   *vname is a string with variable name
   *layer is vertical layer
   *conv_fac is a scale factor to convert units
    eg use 0.01 to convert ice concentration from percentage to fraction
   *wave_mask: if true, mask field where H_s<.01m
   *ice_mask:  if true, mask field where fice<.15
   *vec_opt:
    - vec_opt=0   - plot field as is
    - 0<vec_opt<5 - variable needs to start with u (velocity) or taux (stress)
    - vec_opt=1   - plot vector magnitude
    - vec_opt=2   - plot vector magnitude with directions as unit vectors on top
    - vec_opt=3   - quiver plot of vector with length proportional to magnitude
    - vec_opt=4   - plot direction as scalar (0=North, 90=East)
    - vec_opt=5   - input is direction - convert to unit vectors
   *dir_from: determines if direction is from or to 
   *lower_limit: mask variable where lower than this
   *upper_limit: mask variable where greater than this
   """

   def __init__(self,vname,layer=0,vec_opt=0,conv_fac=1,\
      wave_mask=False,ice_mask=False,dir_from=True,\
      lower_limit=None,upper_limit=None):

      self.name         = vname
      self.layer        = layer
      self.vec_opt      = vec_opt
      self.conv_fac     = conv_fac
      self.ice_mask     = ice_mask
      self.wave_mask    = wave_mask
      self.dir_from     = dir_from
      self.lower_limit  = lower_limit
      self.upper_limit  = upper_limit

      return
###########################################################


###########################################################
def check_pair(var_opts1,var_opts2):
   """
   Check pairs of mod_reading.make_plot_options instances are compatible
   """

   # ====================================================================
   # check options
   if var_opts1.vec_opt in [2,3,5]:
      print('Variable : '+var_opts1.name)
      print('vec_opt  : '+var_opts1.vec_opt)
      raise ValueError('Invalid plot option for pcolor plot')

   if var_opts2.vec_opt in [0,1,4]:
      print('Variable : '+var_opts1.name)
      print('vec_opt  : '+var_opts1.vec_opt)
      raise ValueError('Invalid plot option for quiver plot')
   # ====================================================================
   return
###########################################################

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
class var_object:
   """
   vbl=var_object(vals,mask_in=None,extra_atts=None):
   *vals is an array or masked array
   *mask_in is a bool array
   *xtra_atts=[attlist,attvals] - attlist and attvals are lists
   of attribute names and values for output
   *vbl.values is a masked array 
    - also has shape and ndim att's, and min,max methods
   """
   def __init__(self,vals,mask_in=None,extra_atts=None):

      if extra_atts is not None:
         attlist,attvals = extra_atts
         for i,att in enumerate(attlist):
            setattr(self,att,attvals[i])

      self.shape  = vals.shape
      self.ndim   = vals.ndim

      ####################################################
      # create masked array
      if hasattr(vals,'mask'):
         # already a masked array
         if mask_in is None:
            # don't need to change mask
            self.values = vals
         else:
            # if additional mask is passed in,
            # take union of masks
            mask        = np.logical_or(vals.mask,mask_in)
            self.values = np.ma.array(vals.data,mask=mask)
      else:
         # convert to masked array
         # - check for NaNs
         mask  = np.logical_not(np.isfinite(vals))

         if mask_in is not None:
            # if additional mask is passed in,
            # take union of masks
            mask        = np.logical_or(mask,mask_in)

         self.values = np.ma.array(vals,mask=mask)
      ####################################################

      # activate self[:,:] 
      self.__getitem__  = self.values.__getitem__
      return
   #######################################################


   #######################################################
   def min(self):
      return self.values.min()

   def max(self):
      return self.values.max()
##########################################################


##########################################################
def nc_get_var(ncfil,vblname,time_index=None):
   """
   vbl=nc_get_var(ncfil,vblname,time_index=None):
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
   vbl   = var_object(vals,extra_atts=[attlist,attvals])

   return vbl
########################################################


########################################################
class nc_getinfo:

   #####################################################
   def __init__(self,ncfil,time_index=None,lonlat_file=None):

      ##################################################
      self.filename  = ncfil
      if ncfil[0]=='/':
         self.basedir   = '/'
      else:
         import os
         self.basedir   = os.getcwd()+'/'

      ss = ncfil.split('/')
      for i in range(len(ss)-1):
         self.basedir   = self.basedir+'/'

      self.basename  = ss[-1].strip('.nc')
      ##################################################

      # added here manually
      # - TODO could possibly be determined
      #   from netcdf metadata though
      # - could also be an input
      self.reftime_sig  = 'start of forecast'


      # open the file
      nc                = ncopen(ncfil)
      self.dimensions   = nc.dimensions.keys()

      # is time a dimension?
      self.time_dim     = ('time' in self.dimensions)

      # get global netcdf attributes
      class ncatts:
         def __init__(self,nc):
            for att in nc.ncattrs():
               attval   = getattr(nc,att)
               setattr(self,att,attval)
            return
         def list(self):
            return vars(self).keys()

      self.ncattrs   = ncatts(nc)

      dkeys = nc.dimensions.keys()
      vkeys = nc.variables.keys()
      Nkeys = len(vkeys)

      ########################################################
      # time info:
      if self.time_dim:

         time        = nc.variables['time']
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
               if (split2[loop_i])==1:
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
         sec0     = int(ctime[4:6])
         refpoint = datetime(year0,mon0,day0,hr0,min0,sec0)
         #
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
         for dkey in nc.dimensions.keys():
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

         self.proj_info = proj_obj(att_list_full,att_vals_full)
      else:
         self.proj_info = []
      ########################################################
      
      ########################################################
      # variable list
      vlist = []
      bkeys = [proj_name,self.lonname,self.latname,'model_depth']
      bkeys.extend(dkeys)
      for key in vkeys:
         if key not in bkeys:
            vlist.append(key)

      self.variable_list   = vlist
      self.variables       = vlist
      ########################################################

      nc.close()
      return
   ###########################################################


   ###########################################################
   def timeval_to_datetime(self,timeval):
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
      if not self.lonlat_dim:
         lon   = nc.variables[self.lonname][:,:]
         lat   = nc.variables[self.latname][:,:]
      else:
         lon      = nc.variables[self.lonname][:]
         lat      = nc.variables[self.latname][:]
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

      # conc can have multiple names
      vlist    = self.variable_list
      cnames   = ['fice','ficem','icec','ice_conc']
      if vname in cnames:
         for vi in cnames:
            if vi in vlist:
               vname = vi

      # thickness can have multiple names
      hnames   = ['hice','hicem','icetk']
      if vname in hnames:
         for vi in hnames:
            if vi in vlist:
               vname = vi

      vbl   = nc_get_var(self.filename,vname,time_index=time_index)

      return vbl
   ###########################################################


   ###########################################################
   def plot_var(self,var_opts,time_index=0,\
         pobj=None,bmap=None,HYCOMreg='TP4',\
         clim=None,add_cbar=True,clabel=None,show=True,\
         test_lonlats=None):

      from mpl_toolkits.basemap import Basemap, cm

      if type(var_opts)==type('hi'):
         # if only a string is passed in
         print("Converting string '"+var_opts+"' to mod_reading.make_plot_options object")
         print("- create and pass in such an object directly")
         print("to specify more complicated plot options:")
         print("ie")
         print("var_opts=mod_reading.make_plot_options('"+var_opts+"',")
         print("   vec_opt=0,conv_fac=1,wave_mask=False,ice_mask=False,dir_from=True)")
         var_opts = make_plot_options(var_opts)

      var_opts    = check_var_opts(var_opts,self.variable_list)
      vname       = var_opts.name
      vec_opt     = var_opts.vec_opt
      conv_fac    = var_opts.conv_fac
      ice_mask    = var_opts.ice_mask
      wave_mask   = var_opts.wave_mask
      dir_from    = var_opts.dir_from

      if pobj is None:
         pobj  = plot_object()

      fig,ax,cbar = pobj.get()

      if bmap is None:
         # make basemap
         bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      if clim is not None:
         vmin,vmax   = clim
      else:
         vmin  = None
         vmax  = None

      lon,lat  = self.get_lonlat()
      vbl      = self.get_var(vname,time_index=time_index)
      mask     = vbl.values.mask


      ################################################################## 
      # add ice or wave masks
      if ice_mask and wave_mask:
         fice  = self.get_var('fice',time_index=time_index)
         mask1 = 1-np.logical_and(1-mask,fice[:,:]>.01) #0 if finite,non-low conc
         #
         Hs    = self.get_var('swh',time_index=time_index)
         mask2 = 1-np.logical_and(1-mask,Hs[:,:]>.01) #0 if finite,non-low waves
         #
         mask  = np.logical_or(mask1,mask2)

      elif ice_mask:
         fice  = self.get_var('fice',time_index=time_index)
         mask  = 1-np.logical_and(1-mask,fice[:,:]>.01) #0 if finite,non-low conc

      elif wave_mask:
         Hs    = self.get_var('swh',time_index=time_index)
         mask  = 1-np.logical_and(1-mask,Hs[:,:]>.01) #0 if finite,non-low waves
      ################################################################## 


      ################################################################## 
      if vec_opt==0:

         # just plot scalar
         U     = None # no quiver plot
         Marr  = vbl.values.data
         Marr  = np.ma.array(conv_fac*Marr,mask=mask) #masked array


      elif vec_opt==1:

         # plot vector magnitude
         U  = None # no quiver plot
         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         vbl2  = self.get_var(vname2,time_index=time_index)
         U     = vbl.values.data
         V     = vbl2.values.data

         # speed
         spd   = np.hypot(U,V)
         Marr  = np.ma.array(conv_fac*spd,mask=mask) #masked array

      elif vec_opt==2:

         # plot vector magnitude + direction (unit vectors)
         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         vbl2  = self.get_var(vname2,time_index=time_index)
         U     = conv_fac*vbl.values.data
         V     = conv_fac*vbl2.values.data

         # speed
         spd   = np.hypot(U,V)
         Marr  = np.ma.array(spd,mask=mask) #masked array

         # rotate vectors
         U,V   = bmap.rotate_vector(U,V,lon,lat)

         #unit vectors
         U     = np.ma.array(U/spd,mask=mask)
         V     = np.ma.array(V/spd,mask=mask)

      elif vec_opt==3:

         # plot vector direction only
         # (no pcolor, but vector length is proportional to magnitude)
         Marr  = None

         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         vbl2  = self.get_var(vname2,time_index=time_index)
         U     = conv_fac*vbl.values.data
         V     = conv_fac*vbl2.values.data

         # speed
         spd   = np.hypot(U,V)
         avg   = np.mean(np.ma.array(spd,mask=mask))
         print('avg speed: '+str(avg))

         # rotate vectors
         U,V   = bmap.rotate_vector(U,V,lon,lat)

         # scale by the average speed
         # TODO: add key
         U  = np.ma.array(U/avg,mask=mask)
         V  = np.ma.array(V/avg,mask=mask)

      elif vec_opt==4:

         # plot direction as scalar
         U  = None

         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         vbl2  = self.get_var(vname2,time_index=time_index)
         dir   = 180/np.pi*np.arctan2(vbl2.values.data,vbl.values.data)#dir-to in degrees (<180,>-180)
         dir   = 90-dir #north is 0,angle clockwise
         if dir_from:
            # direction-from
            dir[dir>0]  = dir[dir>0]-360
            Marr        = np.ma.array(dir+180,mask=np.logical_or(mask,1-np.isfinite(dir)))
         else:
            # direction-to
            dir[dir>180]   = dir[dir>180]-360
            Marr           = np.ma.array(dir,mask=np.logical_or(mask,1-np.isfinite(dir)))

      elif vec_opt==5:
         #vbl is a direction - convert to vector
         Marr  = None
         dir   = 90-vbl.values.data
         if dir_from:
            dir   = np.pi/180*(dir+180)
         else:
            dir   = np.pi/180*dir

         # rotate unit vectors
         U,V   = bmap.rotate_vector(np.cos(dir),np.sin(dir),lon,lat)
         U     = np.ma.array(U,mask=mask)
         V     = np.ma.array(V,mask=mask)
      ################################################################## 


      ################################################################## 
      # pcolor plot
      if Marr is not None:

         #########################################################################
         # add additional masking (too low or too high)
         if (var_opts.lower_limit is not None) or (upper_limit is not None):
            mask  = 1*Marr.mask
            data  = Marr.data
            good  = np.logical_not(mask)
            if (var_opts.lower_limit is not None) and (var_opts.upper_limit is not None):
               mask[good]  = np.logical_or(data[good]<var_opts.lower_limit,data[good]>var_opts.upper_limit)
            elif (var_opts.lower_limit is not None):
               mask[good]  = (data[good]<var_opts.lower_limit)
            elif (var_opts.upper_limit is not None):
               mask[good]  = (data[good]>var_opts.upper_limit)
            Marr  = np.ma.array(data,mask=mask)
         #########################################################################

         PC = bmap.pcolor(lon,lat,Marr,latlon=True,ax=ax,vmin=vmin,vmax=vmax)

         if add_cbar:

            if cbar is None:
               cbar  = fig.colorbar(PC)
            else:
               cbar  = fig.colorbar(PC,cax=cbar.ax)

            pobj  = plot_object(fig=fig,ax=ax,cbar=cbar,axpos=pobj.axpos)
            if clabel is not None:
               cbar.set_label(clabel,rotation=270,labelpad=20,fontsize=16)
      ################################################################## 


      ################################################################## 
      if pobj.axpos is not None:
         # try to make sure axes don't move round
         pobj.ax.set_position(pobj.axpos)
      ################################################################## 


      ################################################################## 
      # quiver plot
      if U is not None:
         dens  = 10   # density of arrows
         scale = 50
         QP    = bmap.quiver(lon[::dens,::dens],lat[::dens,::dens],\
                              U[::dens,::dens],V[::dens,::dens],\
                              latlon=True,scale=scale,ax=ax)

      if test_lonlats is not None:
         for lont,latt in test_lonlats:
            bmap.plot(lont,latt,'^m',markersize=5,latlon=True,ax=ax)

      Fplt.finish_map(bmap,ax=ax)
      if show:
         fig.show()

      return pobj,bmap
   ###########################################################


   ###########################################################
   def plot_var_pair(self,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs):

      # ====================================================================
      # check names
      var_opts1   = check_var_opts(var_opts1,self.variable_list)
      var_opts2   = check_var_opts(var_opts2,self.variable_list)

      # check options
      check_pair(var_opts1,var_opts2)
      # ====================================================================

      pobj,bmap   = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,**kwargs)
      self.plot_var(var_opts2,pobj=pobj,bmap=bmap,**kwargs)

      return pobj,bmap
   ###########################################################

   ###########################################################
   def make_png(self,var_opts,pobj=None,bmap=None,figdir='.',time_index=0,date_label=True,**kwargs):

      var_opts    = check_var_opts(var_opts,self.variable_list)

      new_fig  = (pobj is None)
      if new_fig:
         pobj  = plot_object()

      pobj,bmap   = self.plot_var(var_opts,pobj=pobj,bmap=bmap,time_index=time_index,**kwargs)

      dtmo     = self.datetimes[time_index]
      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')

      if date_label:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         pobj.ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)

      if pobj.axpos is not None:
         # try to make sure axes don't move round
         pobj.ax.set_position(pobj.axpos)

      vname    = var_opts.name
      Fname    = vname
      vec_opt  = var_opts.vec_opt

      if vec_opt==1:
         #magnitude only
         if vname in ['u','usurf']:
            Fname = 'surf_speed'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_speed'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_magnitude'

      elif vec_opt==2 or vec_opt==3:
         #quiver plots on top of magnitude or by itself
         if vname in ['u','usurf']:
            Fname = 'surf_vel'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress'

      elif vec_opt==4:
         #direction as scalar
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      elif vec_opt==5:
         #direction -> vector
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      Fname    = Fname.strip('_')
      figname  = figdir+'/'+self.basename+'_'+Fname+datestr+'.png'

      print('Saving to '+figname) 
      pobj.fig.savefig(figname)

      if new_fig:
         pobj.ax.cla()
         pobj.fig.clear()
         plt.close(pobj.fig)

      return pobj,bmap
   ###########################################################


   ###########################################################
   def make_png_pair(self,var_opts1,var_opts2,time_index=0,\
         pobj=None,bmap=None,figdir='.',date_label=True,**kwargs):

      # ====================================================================
      # check names
      var_opts1   = check_var_opts(var_opts1,self.variable_list)
      var_opts2   = check_var_opts(var_opts2,self.variable_list)

      # check options
      check_pair(var_opts1,var_opts2)
      # ====================================================================

      new_fig  = (pobj is None)
      if new_fig:
         pobj  = plot_object()

      pobj,bmap   = self.plot_var_pair(var_opts1,var_opts2,time_index=time_index,\
                     pobj=pobj,bmap=bmap,**kwargs)
      fig,ax,cbar = pobj.get()

      dtmo     = self.datetimes[time_index]
      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')

      if date_label:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)

      # set name with 1st variable only
      Fname    = var_opts1.name
      vec_opt  = var_opts1.vec_opt
      if vec_opt==1:
         #magnitude only
         if vname in ['u','usurf']:
            Fname = 'surf_speed'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_speed'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_magnitude'

      elif vec_opt==4:
         #direction as scalar
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      Fname    = Fname.strip('_')
      figname  = figdir+'/'+self.basename+'_'+Fname+datestr+'.png'

      print('Saving to '+figname) 
      fig.savefig(figname)

      if new_fig:
         ax.cla()
         fig.clear()
         plt.close(fig)

      return pobj,bmap
   ###########################################################


   ###########################################################
   def make_png_all(self,var_opts,HYCOMreg='TP4',figdir='.',**kwargs):

      # check names
      var_opts    = check_var_opts(var_opts,self.variable_list)
      pobj        = plot_object()
      fig,ax,cbar = pobj.get()

      bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      N  = len(self.timevalues)
      for i in range(N):

         pobj,bmap   = self.make_png(var_opts,\
                           bmap=bmap,time_index=i,\
                           figdir=figdir,show=False,**kwargs)

         ax.cla()
         if pobj.cbar is not None:
            pobj.cbar.ax.clear()   # cbar.ax.clear()

         print('\n'+str(i+1)+' records done out of '+str(N))

      plt.close(fig)
      return
   ###########################################################


   ###########################################################
   def make_png_pair_all(self,var_opts1,var_opts2,HYCOMreg='TP4',figdir='.',**kwargs):

      # ====================================================================
      # check names
      var_opts1   = check_var_opts(var_opts1,self.variable_list)
      var_opts2   = check_var_opts(var_opts2,self.variable_list)

      # check options
      check_pair(var_opts1,var_opts2)
      # ====================================================================

      pobj        = plot_object()
      fig,ax,cbar = pobj.get()
      bmap        = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      N  = len(self.timevalues)
      for i in range(N):

         pobj,bmap   = self.make_png_pair(var_opts1,var_opts2,\
                        pobj=pobj,bmap=bmap,time_index=i,\
                        figdir=figdir,show=False,**kwargs)

         if i==0:
            # Fix axes position to stop it moving round
            pobj  = pobj.renew(axpos=pobj.ax.get_position())

         ax.cla()
         if pobj.cbar is not None:
            pobj.cbar.ax.clear()   # cbar.ax.clear()

         print('\n'+str(i+1)+' records done out of '+str(N))

      plt.close(fig)
      return
   ###########################################################


   ###########################################################
   def MIZmap(self,var_name='dmax',time_index=0,do_sort=False,EastOnly=True,\
         plotting=True,HYCOM_region='Arctic',**kwargs):
      """
      Call  : self.MIZmap(var_name='dmax',do_sort=False,EastOnly=True,plotting=True,**kwargs):
      Inputs:
         var_name is variable to find MIZ from
         **kwargs to be passed onto MIZchar.get_MIZ_poly:
            outdir='.',do_sort=True
      Returns: MIZchar.MIZpoly object
      """

      import MIZchar as mc
      vname = check_names(var_name,self.variables)

      if var_name == 'dmax':
         # FSD MIZarray(1-
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,300]# for plotting
         lower_limit = .1     # for plotting
      elif var_name == 'fice':
         # conc MIZ
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,1]  # for plotting
         lower_limit = .15    # for plotting
      elif var_name == 'hice':
         # thin ice areas
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,2.] # for plotting
         lower_limit = .01    # for plotting
      else:
         raise ValueError('Wrong selection variable for MIZmap')

      print("MIZchar.get_MIZ_poly\n")
      lon,lat  = self.get_lonlat()
      MPdict   = {}
      tfiles   = {}

      if do_sort:
         # possible regions are:
         regions  = ['gre','bar','beau','lab','balt','les','can']

         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt' )
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')

         # for reg in ['gre']:
         for reg in regions:
            mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,region=reg)
            MPdict.update({reg:mp})

            fname0   = self.basename+'_'+var_name +'_'+reg
            tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
            if 'all' in tfile.keys():
               tfiles.update({reg:tfile['all']})

         if 0:
            MPdict['gre'].show_maps()
            return MPdict

      else:
         reg   = 'all'
         mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name)
         MPdict.update({reg:mp})
         #
         fname0   = self.basename+'_'+var_name
         tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
         if 'all' in tfile.keys():
            tfiles.update({reg:tfile['all']})

      Pdict    = {}
      PLOTTING = False
      for reg in tfiles.keys():

         ##########################################################
         # filenames
         tfil     = tfiles[reg]                          # text file with polygon outlines characterized
         figname  = tfil.replace('.txt','.png')          # plot of polygons
         shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
         sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
         ##########################################################


         ##########################################################
         if do_sort:
            mapreg   = reg
         else:
            mapreg   = HYCOM_region
         ##########################################################


         ##########################################################
         # process each text file to get MIZ width etc
         print("MIZchar.single_file: "+tfil+"\n")
         bmap     = Fplt.start_HYCOM_map(mapreg)
         Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
         Pdict.update({reg:Psolns})
         
         # Save summary & shapefile
         mc.save_summary  (Psolns,sumname)
         mc.save_shapefile(Psolns,filename=shpname)
         ##########################################################

         
         if plotting:
            ##########################################################
            # Make plot
            var_opts = make_plot_options(vname,lower_limit=lower_limit)
            pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=clim)[0]
            fig      = pobj.fig
            ax       = pobj.ax
            PLOTTING = True

            for MIZi in Psolns:
               # plot outlines of polygons
               lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
               bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)

               Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
               if Wavg>26:
                  MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)

                  # add text with mean width
                  xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
                  xav                  = (xmin+xmax)/2.
                  ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
                     color='k',fontsize=16,horizontalalignment='right',\
                     verticalalignment='top')

            Fplt.finish_map(bmap)
            print('Saving '+figname)
            fig.savefig(figname)
            # plt.show(fig)
            ax.cla()
            fig.clear()
            # finished region
            ##########################################################

      if PLOTTING:
         plt.close(fig)
      return mp,Pdict,tfiles
   ###########################################################


   ###########################################################
   def areas_of_disagreement(self,obs_type='OSISAF',time_index=0,do_sort=True,EastOnly=True,plotting=True,**kwargs):
      # kwargs: outdir='.',do_sort=True

      import MIZchar as mc

      if obs_type == 'OSISAF':
         var_name    = 'fice'
         lower_limit = .15
         bmap        = basemap_OSISAF()
         #
         cyear = self.datetimes[time_index].strftime('%Y')
         cdate = self.datetimes[time_index].strftime('%Y%m%d')
         obsfil   = '/work/shared/nersc/msc/OSI-SAF/'+cyear+\
                     '_nh_polstere/ice_conc_nh_polstere-100_multi_'+\
                     cdate+'1200.nc'
      else:
         raise ValueError('Wrong selection variable for areas_of_disagreement')

      vname = check_names(var_name,self.variables)

      # observation grid & compared quantity
      nci         = nc_getinfo(obsfil)
      lon2,lat2   = nci.get_lonlat()
      Xobs,Yobs   = bmap(lon2,lat2)

      vname2   = check_names(var_name,nci.variables)
      Zobs     = nci.get_var(vname2)

      # model grid & compared quantity
      Zmod        = self.get_var(vname,time_index=time_index)
      lon,lat     = self.get_lonlat()
      Xmod,Ymod   = bmap(lon,lat)

      if 1:
         #Zref,Zint should be np.ma.array
         lon_ref,lat_ref   = lon2,lat2
         Xref,Yref,Zref    = Xobs,Yobs,Zobs.values  # obs grid is reference;                 
         Xint,Yint,Zint    = Xmod,Ymod,Zmod.values  # to be interped from model grid onto obs grid;  Zint is np.ma.array

      # add the mask for the ref to Arr
      Arr   = reproj_mod2obs(Xint,Yint,Zint,Xref,Yref,mask=1*Zref.mask)

      # add the mask for Arr to Zref
      Zref  = np.ma.array(Zref.data,mask=Arr.mask)

      MPdict   = {'Over':{},'Under':{}}
      tfiles   = {'Over':{},'Under':{}}

      if 0:
         # test interpolation and matching of masks
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,2,1)
         ax1.imshow(Arr)
         ax2   = fig.add_subplot(1,2,2)
         ax2.imshow(Zref)
         plt.show(fig)
         return

      if do_sort:
         # possible regions are:
         regions  = ['gre','bar','beau','lab','balt','les','can']

         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt' )
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')

         # for reg in ['bar']:
         for reg in regions:

            # Arr,Zref are np.ma.array objects
            Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref,region=reg)
            MPdict['Over'] .update({reg:Over})
            MPdict['Under'].update({reg:Under})

            for OU in ['Over','Under']:

               fname0   = self.basename+'_v'+obs_type +'_'+OU+'_'+reg
               tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
               if 'all' in tfile.keys():
                  tfiles[OU].update({reg:tfile['all']})

         if 0:
            MPdict['Over'] [reg].show_maps()
            MPdict['Under'][reg].show_maps()
            return MPdict
      else:
         reg         = 'all'
         Over,Under  = mc.get_AOD_polys(Arr.values,Zref.values,lon_ref,lat_ref)
         MPdict['Over'] .update({reg:Over})
         MPdict['Under'].update({reg:Under})

         for OU in ['Over','Under']:

            fname0   = self.basename+'_v'+obs_type+'_'+OU+'_'+reg
            tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
            if 'all' in tfile.keys():
               tfiles[OU].update({reg:tfile['all']})

      print(tfiles)
      print(MPdict)
      Pdict = {'Over':{},'Under':{}}
      for OU in ['Over','Under']:
         PLOTTING = False
         for reg in tfiles[OU].keys():

            ##########################################################
            # filenames
            tfil     = tfiles[OU][reg]                          # text file with polygon outlines characterized
            figname  = tfil.replace('.txt','.png')          # plot of polygons
            shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
            sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
            ##########################################################


            ##########################################################
            if do_sort:
               mapreg   = reg
            else:
               mapreg   = self.HYCOM_region
            ##########################################################


            ##########################################################
            # process each text file to get MIZ width etc
            print("MIZchar.single_file: "+tfil+"\n")
            bmap     = Fplt.start_HYCOM_map(mapreg)
            Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
            Pdict[OU].update({reg:Psolns})
            
            # Save summary & shapefile
            mc.save_summary  (Psolns,sumname)
            mc.save_shapefile(Psolns,filename=shpname)
            ##########################################################

            
            if plotting:
               ##########################################################
               # Make plot
               var_opts = make_plot_options(vname,lower_limit=lower_limit)
               pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=[0,1])[0]
               fig      = pobj.fig
               ax       = pobj.ax
               PLOTTING = True

               for MIZi in Psolns:
                  # plot outlines of polygons
                  lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
                  bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)

                  Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
                  if Wavg>26:
                     MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)

                     # add text with mean width
                     xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
                     xav                  = (xmin+xmax)/2.
                     ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
                        color='k',fontsize=16,horizontalalignment='right',\
                        verticalalignment='top')

               Fplt.finish_map(bmap)
               print('Saving '+figname)
               fig.savefig(figname)
               # plt.show(fig)
               ax.cla()
               fig.clear()
               # finished region
               ##########################################################

         if PLOTTING:
            plt.close(fig)

      return MPdict,tfiles,Pdict
   ###########################################################

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

   import struct

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
def get_array_from_HYCOM_binary(afile,recno,dims=None,grid_dir='.'):
   # routine to get the array from the .a (binary) file
   # * recno=1 is 1st record 
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)

   import struct

   ######################################################################
   # get size of grid
   if dims is None:
      if 'regional.grid' in afile:
         # afile is a grid file
         # - check .b file for size of grid
         bfile = afile[:-2]+'.b'
      else:
         # check regional.grid.b file for size of grid
         bfile = grid_dir+'/regional.grid.b'
         if os.path.exists(bfile):
            sys.exit('Grid file not present: '+bfile)

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
   fmt_size = 4      # HYCOM files are single precision
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
      aid.seek(rec_size,1) # seek in bytes (1: reference is current position)
   ######################################################################

   # read data and close file
   data  = aid.read(rec_size)
   aid.close()

   # rearrange into correctly sized array
   fld   = struct.unpack('>'+Nhyc*fmt_py,data) # NB BIG-ENDIAN so need '>'
   fld   = np.array(fld[0:nx*ny]) # select the 1st nx,ny - rest of the Nhyc record is rubbish
   fld   = fld.reshape((ny,nx)).transpose()  # need to transpose because of differences between
                                             # python/c and fortran/matlab 

   land_thresh          = 1.e30# on land 1.2677e30 
   fld[fld>land_thresh] = np.nan

   return fld
##############################################################

##############################################################
def get_record_numbers_HYCOM(bfile):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)


   bid   = open(bfile,'r')
   word  = bid.readline().split()[0] # 1st word in line 
   while word!='field':
      word  = bid.readline().split()[0] # 1st word in line 

   # have found table title
   n     = 0
   lut2d = {}
   min2d = {}
   max2d = {}
   lut3d = {}
   min3d = {}
   max3d = {}
   lin   = bid.readline()
   EOF   = (lin=='')
   while not EOF:
      n     = n+1
      word  = lin.split()[0]        # 1st word in line 
      layer = int(lin.split()[4])   # layer number (5th entry)
      xmin  = float(lin.split()[6]) # min where defined
      xmax  = float(lin.split()[7]) # max where defined

      if layer==0:
         # surface/2D var
         lut2d.update({word:n})
         min2d.update({word:xmin})
         max2d.update({word:xmax})
      else:
         # 3D var
         LUT   = {layer:n}
         MIN   = {layer:xmin}
         MAX   = {layer:xmax}
         if word not in lut3d.keys():
            lut3d.update({word:LUT})
            min3d.update({word:MIN})
            max3d.update({word:MAX})
         else:
            lut3d[word].update(LUT)
            min3d[word].update(MIN)
            max3d[word].update(MAX)
      #
      lin   = bid.readline()
      EOF   = (lin=='')

   bid.close()

   return [lut2d,min2d,max2d],[lut3d,min3d,max3d]
   ######################################################################


######################################################################
class HYCOM_binary_info:
   def __init__(self,fname,gridpath=None):
      from datetime import datetime,timedelta

      ss = fname.split('.')
      if ss[-1]=='a':
         self.afile = fname
         self.bfile = fname[:-1]+'b'
      elif ss[-1]=='b':
         self.bfile = fname
         self.afile = fname[:-1]+'a'
      else:
         raise ValueError('HYCOM binaries should have extensions .a or .b')

      basename                = fname.split('/')[-1]
      self.basename           = basename[:-2]
      self.HYCOM_region       = basename[:3]

      # info from bfile
      lut2d,lut3d = get_record_numbers_HYCOM(self.bfile)
      #
      self.record_numbers  ,self.minvals2d,self.maxvals2d   = lut2d
      self.record_numbers3d,self.minvals3d,self.maxvals3d   = lut3d

      self.variables          = self.record_numbers.keys()
      self.variables3d        = self.record_numbers3d.keys() 
      self.all_variables      = 1*self.variables
      self.all_variables.extend(1*self.variables3d)

      #######################################################################
      #path to regional.grid.[ab] files
      if gridpath is not None:
         self.gridpath = gridpath
      else:
         wsn            = '/work/shared/nersc/msc/ModelInput'
         gridpath_lut   = {'FR1':wsn+'/FramStrait_Hyc2.2.12/FR1a0.03-clean//topo',\
                           'BS1':wsn+'/BS1a0.045-clean/topo',\
                           'TP4':wsn+'/../REANALYSIS/topo'}
         self.gridpath = gridpath_lut[self.HYCOM_region]

      # get grid size
      bid   = open(self.gridpath+'/regional.grid.b','r')
      line  = bid.readline()
      while 'idm' not in line:
         line  = bid.readline()
      nx    = int( line.split()[0] )

      line  = bid.readline()
      while 'jdm' not in line:
         line  = bid.readline()
      ny    = int( line.split()[0] )

      self.dims   = [nx,ny]
      self.Nx     = nx
      self.Ny     = ny
      #######################################################################


      #######################################################################
      # date
      bid   = open(self.bfile)
      line  = bid.readline()
      its   = 0
      while 'model day' not in line and its<1200:
         its   = its+1
         line  = bid.readline()

      line              = bid.readline()
      self.time_value   = float(line.split()[-5]) # model time (days)
      bid.close()

      self.reference_date  = datetime(1900,12,31)
      self.datetime        = self.reference_date+timedelta(self.time_value)
      self.datetimes       = [self.datetime]
      self.time_values     = [self.time_value]
      self.time_units      = 'days'
      #######################################################################

      return # __init__
   #######################################################################


   #######################################################################
   def get_lonlat(self):

      gfil  = self.gridpath+'/regional.grid.a'
      plon  = get_array_from_HYCOM_binary(gfil,1,\
                  dims=self.dims,grid_dir=self.gridpath)
      plat  = get_array_from_HYCOM_binary(gfil,2,\
                  dims=self.dims,grid_dir=self.gridpath)

      return plon,plat
   #######################################################################


   #######################################################################
   def get_fixed_lonlat(self,bmap):
      gfil     = self.gridpath+'/regional.grid.a'

      if 1:
         # get ulon,vlat
         lon   = get_array_from_HYCOM_binary(gfil,5,\
                     dims=self.dims,grid_dir=self.gridpath)
         lat   = get_array_from_HYCOM_binary(gfil,8,\
                     dims=self.dims,grid_dir=self.gridpath)

      else:
         # try to fix plon,plat with scpx,scpy
         lon,lat  = self.get_lonlat()
         X,Y      = bmap(lon,lat)

         scpx  = get_array_from_HYCOM_binary(gfil,10,\
                     dims=self.dims,grid_dir=self.gridpath)
         scpy  = get_array_from_HYCOM_binary(gfil,11,\
                     dims=self.dims,grid_dir=self.gridpath)
         
         nx,ny = self.dims
         X2    = np.zeros((nx+1,ny+1))
         Y2    = np.zeros((nx+1,ny+1))
         #
         X2[1:,1:]   = 1*X    +.5*scpx
         X2[0,1:]    = X[0,:] -.5*scpx[0,:]
         X2[1:,0]    = X[:,0] -.5*scpx[:,0]
         X2[0,0]     = X[0,0] -.5*scpx[0,0]
         #
         Y2[1:,1:]   = 1*Y    +.5*scpy
         Y2[0,1:]    = Y[0,:] -.5*scpy[0,:]
         Y2[1:,0]    = Y[:,0] -.5*scpy[:,0]
         Y2[0,0]     = Y[0,0] -.5*scpy[0,0]

         # return new ones
         lon,lat  = bmap(X2,Y2,inverse=True)

      return lon,lat
   #######################################################################


   #######################################################################
   def get_depths(self):

      dfil     = self.gridpath+'/regional.depth.a'
      depths   = get_array_from_HYCOM_binary(dfil,1,\
                     dims=self.dims,grid_dir=self.gridpath)

      return depths
   #######################################################################


   #######################################################################
   def get_var(self,vname):

      if type(vname)!=type([]):
         # 2d var
         layer = 0
         vname = check_names(vname,self.variables)
         recno = self.record_numbers[vname]
         xmin  = self.minvals2d[vname]
         xmax  = self.maxvals2d[vname]
         #
         vbl   = get_array_from_HYCOM_binary(self.afile,recno,\
                     dims=self.dims)
         mask  = np.array(1-np.isfinite(vbl),dtype='bool')
      else:
         # 3d var
         vname,layer = vname
         vname       = check_names(vname,self.variables3d)
         recno       = self.record_numbers3d[vname][layer]
         xmin        = self.minvals3d       [vname][layer]
         xmax        = self.maxvals3d       [vname][layer]
         vbl         = get_array_from_HYCOM_binary(self.afile,recno,\
                           dims=self.dims)

      extra_atts  = [['dimensions'],['i','j']]
      vbl         = var_object(vbl,extra_atts=extra_atts)

      ########################################################
      # consistency check between afile and bfile
      # TODO debug this
      if 0:
         if xmax==0:
            ddx   = 1.e-8
         else:
            ddx   = 1e-8*abs(xmax)

         if abs(vbl.max()-xmax)>ddx:
            ss = 'Maximum of '+vname+'('+str(layer)+') '+\
                 'inconsistent with '+self.bfile
            print(ss)
            print(self.afile+': '+str(vbl.max()) )
            print(self.bfile+': '+str(xmax)      )
            raise ValueError()

         if xmin==0:
            ddx   = 1.e-8
         else:
            ddx   = 1e-8*abs(xmin)

         if abs(vbl.min()-xmin)>ddx:
            ss = 'Minimum of '+vname+'('+str(layer)+') '+\
                 'inconsistent with '+self.bfile
            print(ss)
            print(self.afile+': '+str(vbl.min()) )
            print(self.bfile+': '+str(xmin)      )
            raise ValueError()
      ########################################################

      return vbl
   #######################################################################


   #######################################################################
   def plot_var(self,var_opts,\
         pobj=None,bmap=None,HYCOMreg=None,\
         clim=None,add_cbar=True,clabel=None,show=True,\
         test_lonlats=None):

      from mpl_toolkits.basemap import Basemap
      from matplotlib import cm

      var_opts    = check_var_opts(var_opts,self.all_variables)
      vname       = var_opts.name
      layer       = var_opts.layer
      vec_opt     = var_opts.vec_opt
      conv_fac    = var_opts.conv_fac
      ice_mask    = var_opts.ice_mask
      wave_mask   = var_opts.wave_mask
      dir_from    = var_opts.dir_from
      if vname not in self.all_variables:
         raise ValueError('Variable '+vname+'not in '+self.afile)

      if pobj is None:
         pobj  = plot_object()
      fig,ax,cbar = pobj.get()

      if HYCOMreg is None:
         HYCOMreg = self.HYCOM_region
      if bmap is None:
         # make basemap
         bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      # lon,lat  = self.get_lonlat()
      lon,lat  = self.get_fixed_lonlat(bmap)

      if clim is not None:
         vmin,vmax   = clim
      else:
         vmin  = None
         vmax  = None

      if layer==0:
         vbl   = self.get_var(vname)
      else:
         vbl   = self.get_var([vname,layer])

      mask     = vbl.values.mask

      ################################################################## 
      # add ice or wave masks
      if ice_mask and wave_mask:
         fice        = self.get_var('ficem')
         good        = np.array(1-fice.values.mask  ,dtype='bool')
         mask1       = np.zeros(fice.shape         ,dtype='bool')
         mask1[good] = (fice[good]<.01)          # water
         mask1       = np.logical_or(mask,mask1) # water or NaN
         #
         Hs          = self.get_var('swh')
         good        = np.array(1-Hs.values.mask,dtype='bool')
         mask2       = np.zeros(Hs.shape        ,dtype='bool')
         mask2[good] = (Hs[good]<.01)              # no waves
         mask2       = np.logical_or(mask,mask2)   # no waves or NaN
         #
         mask  = np.logical_or(mask1,mask2)

      elif ice_mask:
         fice        = self.get_var('ficem')
         good        = np.array(1-fice.values.mask ,dtype='bool')
         mask1       = np.zeros(fice.shape         ,dtype='bool')
         mask1[good] = (fice[good]<.01)          # water
         mask        = np.logical_or(mask,mask1) # water or NaN

      elif wave_mask:
         Hs          = self.get_var('swh')
         good        = np.array(1-Hs.values.mask,dtype='bool')
         mask2       = np.zeros(Hs.shape        ,dtype='bool')
         mask2[good] = (Hs[good]<.01)              # no waves
         mask        = np.logical_or(mask,mask2)   # no waves or NaN
      ################################################################## 


      ###############################################################
      if vec_opt==0:

         # just plot scalar
         U     = None # no quiver plot
         Marr  = np.ma.array(conv_fac*vbl.values.data,mask=mask)

      elif vec_opt==1:

         # vector magnitude of vel or stress
         U  = None # no quiver plot
         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         if layer==0:
            vbl2  = self.get_var(vname2)
         else:
            vbl2  = self.get_var([vname2,layer])

         Marr  = np.hypot(vbl.values.data,vbl2.values.data)
         Marr  = np.ma.array(conv_fac*Marr,mask=mask)

      elif vec_opt==2:

         # plot vector magnitude + direction (unit vectors)
         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         if layer==0:
            vbl2  = self.get_var(vname2)
         else:
            vbl2  = self.get_var([vname2,layer])

         U  = conv_fac*vbl.values.data
         V  = conv_fac*vbl2.values.data

         # speed (to be plotted)
         spd         = np.hypot(U,V)
         Marr        = np.ma.array(spd,mask=mask) #masked array
         good        = np.logical_not(mask)
         pos         = np.zeros(spd.shape,dtype='bool')
         pos[good]   = (spd[good]>0)
         npos        = np.logical_not(pos)

         # rotate vectors
         U,V      = bmap.rotate_vector(U,V,lon,lat)

         #unit vectors
         U[pos]   = U[pos]/spd[pos]
         V[pos]   = V[pos]/spd[pos]
         U[npos]  = 0.
         V[npos]  = 0.

         # add masks
         U  = np.ma.array(U,mask=mask)
         V  = np.ma.array(V,mask=mask)

      elif vec_opt==3:

         # plot vector direction only
         # (no pcolor, but vector length is proportional to magnitude)
         # TODO: add key
         Marr  = None

         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         if layer==0:
            vbl2  = self.get_var(vname2)
         else:
            vbl2  = self.get_var([vname2,layer])
         U     = conv_fac*vbl.values.data
         V     = conv_fac*vbl2.values.data

         # speed
         spd   = np.hypot(U,V)
         avg   = np.mean(np.ma.array(spd,mask=mask))
         print('avg speed: '+str(avg))

         # rotate vectors
         U,V   = bmap.rotate_vector(U,V,lon,lat)

         # scale by the average speed
         U  = np.ma.array(U/avg,mask=mask)
         V  = np.ma.array(V/avg,mask=mask)

      elif vec_opt==4:

         # plot direction as scalar
         U  = None

         if vname[0]=='u':
            vname2   = 'v'+vname[1:]
         elif vname[:4]=='taux':
            vname2   = 'tauy'+vname[4:]

         if layer==0:
            vbl2  = self.get_var(vname2)
         else:
            vbl2  = self.get_var([vname2,layer])
         dir   = 180/np.pi*np.arctan2(vbl2.values.data,vbl.values.data)#dir-to in degrees (<180,>-180)
         dir   = 90-dir #north is 0,angle clockwise
         if dir_from:
            # direction-from
            dir[dir>0]  = dir[dir>0]-360
            Marr        = np.ma.array(dir+180,mask=np.logical_or(mask,1-np.isfinite(dir)))
         else:
            # direction-to
            dir[dir>180]   = dir[dir>180]-360
            Marr           = np.ma.array(dir,mask=np.logical_or(mask,1-np.isfinite(dir)))

      elif vec_opt==5:
         #vbl is a direction - convert to vector
         Marr  = None
         dir   = 90-vbl.values.data
         if dir_from:
            dir   = np.pi/180*(dir+180)
         else:
            dir   = np.pi/180*dir

         # rotate unit vectors
         U,V   = bmap.rotate_vector(np.cos(dir),np.sin(dir),lon,lat)
         U     = np.ma.array(U,mask=mask)
         V     = np.ma.array(V,mask=mask)
      ################################################################## 


      ################################################################## 
      # pcolor plot
      if Marr is not None:

         #########################################################################
         # add additional masking
         if (var_opts.lower_limit is not None) or (var_opts.upper_limit is not None):
            mask  = 1*Marr.mask
            data  = Marr.data
            good  = np.logical_not(mask)
            if (var_opts.lower_limit is not None) and (var_opts.upper_limit is not None):
               mask[good]  = np.logical_or(data[good]<var_opts.lower_limit,data[good]>var_opts.upper_limit)
            elif (var_opts.lower_limit is not None):
               mask[good]  = (data[good]<var_opts.lower_limit)
            elif (var_opts.upper_limit is not None):
               mask[good]  = (data[good]>var_opts.upper_limit)
            Marr  = np.ma.array(data,mask=mask)
         #########################################################################

         PC = bmap.pcolor(lon,lat,Marr,latlon=True,ax=ax,vmin=vmin,vmax=vmax)
         if add_cbar:

            if cbar is None:
               cbar  = fig.colorbar(PC)
            else:
               cbar  = fig.colorbar(PC,cax=cbar.ax)

            pobj  = plot_object(fig=fig,ax=ax,cbar=cbar,axpos=pobj.axpos)
            if clabel is not None:
               cbar.set_label(clabel,rotation=270,labelpad=20,fontsize=16)
      ################################################################## 


      ################################################################## 
      if pobj.axpos is not None:
         # try to make sure axes don't move round
         pobj.ax.set_position(pobj.axpos)
      ################################################################## 


      ################################################################## 
      # quiver plot
      if U is not None:
         dens  = 10   # density of arrows
         scale = 50
         QP    = bmap.quiver(lon[::dens,::dens],lat[::dens,::dens],\
                              U[::dens,::dens],V[::dens,::dens],\
                              latlon=True,scale=scale,ax=ax)

      if test_lonlats is not None:
         for lont,latt in test_lonlats:
            bmap.plot(lont,latt,'^m',markersize=5,latlon=True,ax=ax)

      Fplt.finish_map(bmap,ax=ax)
      if show:
         fig.show()

      return pobj,bmap
   ###########################################################


   ###########################################################
   def MIZmap(self,var_name='dmax',do_sort=False,EastOnly=True,plotting=True,**kwargs):
      """
      Call  : self.MIZmap(var_name='dmax',**kwargs):
      Inputs:
         var_name is variable to find MIZ from
         **kwargs to be passed onto MIZchar.get_MIZ_poly:
            outdir='.',do_sort=True
      Returns: MIZchar.MIZpoly object
      """

      import MIZchar as mc

      vname = check_names(var_name,self.variables)
      if var_name == 'dmax':
         # FSD MIZarray(1-
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,300]# for plotting
         lower_limit = .1     # for plotting
      elif var_name == 'fice':
         # conc MIZ
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,1]  # for plotting
         lower_limit = .15    # for plotting
      elif var_name == 'hice':
         # thin ice areas
         Arr         = self.get_var(vname,time_index=time_index)
         clim        = [0,2.] # for plotting
         lower_limit = .01    # for plotting
      else:
         raise ValueError('Wrong selection variable for MIZmap')

      print("MIZchar.get_MIZ_poly\n")
      lon,lat  = self.get_lonlat()
      MPdict   = {}
      tfiles   = {}

      if do_sort:
         # possible regions are:
         regions  = ['gre','bar','beau','lab','balt','les','can']

         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt' )
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')

         # for reg in ['gre']:
         for reg in regions:
            mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,region=reg)
            MPdict.update({reg:mp})

            fname0   = self.basename+'_'+var_name +'_'+reg
            tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
            if 'all' in tfile.keys():
               tfiles.update({reg:tfile['all']})

         if 0:
            MPdict['gre'].show_maps()
            return MPdict

      else:
         reg   = 'all'
         mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name)
         MPdict.update({reg:mp})
         #
         fname0   = self.basename+'_'+var_name
         tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
         if 'all' in tfile.keys():
            tfiles.update({reg:tfile['all']})

      Pdict    = {}
      PLOTTING = False
      for reg in tfiles.keys():

         ##########################################################
         # filenames
         tfil     = tfiles[reg]                          # text file with polygon outlines characterized
         figname  = tfil.replace('.txt','.png')          # plot of polygons
         shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
         sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
         ##########################################################


         ##########################################################
         if do_sort:
            mapreg   = reg
         else:
            mapreg   = self.HYCOM_region
         ##########################################################


         ##########################################################
         # process each text file to get MIZ width etc
         print("MIZchar.single_file: "+tfil+"\n")
         bmap     = Fplt.start_HYCOM_map(mapreg)
         Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
         Pdict.update({reg:Psolns})
         
         # Save summary & shapefile
         mc.save_summary  (Psolns,sumname)
         mc.save_shapefile(Psolns,filename=shpname)
         ##########################################################

         
         if plotting:
            ##########################################################
            # Make plot
            var_opts = make_plot_options(vname,lower_limit=lower_limit)
            pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=clim)[0]
            fig      = pobj.fig
            ax       = pobj.ax
            PLOTTING = True

            for MIZi in Psolns:
               # plot outlines of polygons
               lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
               bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)

               Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
               if Wavg>26:
                  MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)

                  # add text with mean width
                  xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
                  xav                  = (xmin+xmax)/2.
                  ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
                     color='k',fontsize=16,horizontalalignment='right',\
                     verticalalignment='top')

            Fplt.finish_map(bmap)
            print('Saving '+figname)
            fig.savefig(figname)
            # plt.show(fig)
            ax.cla()
            fig.clear()
            # finished region
            ##########################################################

      if PLOTTING:
         plt.close(fig)
      return mp,Pdict,tfiles
   ###########################################################
   
   
   ###########################################################
   def areas_of_disagreement(self,obs_type='OSISAF',do_sort=True,EastOnly=True,plotting=True,**kwargs):
      # kwargs: outdir='.',do_sort=True

      import MIZchar as mc

      if obs_type == 'OSISAF':
         var_name = 'fice'
         bmap     = basemap_OSISAF()
         obsfil   = '/work/shared/nersc/msc/OSI-SAF/'+\
               self.datetime.strftime('%Y')+'_nh_polstere/'+\
               'ice_conc_nh_polstere-100_multi_'+\
               self.datetime.strftime('%Y%m%d')+'1200.nc'
      else:
         raise ValueError('Wrong selection variable for areas_of_disagreement')

      # observation grid & compared quantity
      nci         = nc_getinfo(obsfil)
      lon2,lat2   = nci.get_lonlat()
      Xobs,Yobs   = bmap(lon2,lat2)
      Zobs        = nci.get_var(var_name)

      # model grid & compared quantity
      Zmod        = self.get_var(var_name)
      lon,lat     = self.get_lonlat()
      Xmod,Ymod   = bmap(lon,lat)

      if 1:
         #Zref,Zint should be np.ma.array
         lon_ref,lat_ref   = lon2,lat2
         Xref,Yref,Zref    = Xobs,Yobs,Zobs.values  # obs grid is reference;                 
         Xint,Yint,Zint    = Xmod,Ymod,Zmod.values  # to be interped from model grid onto obs grid;  Zint is np.ma.array

      # add the mask for the ref to Arr
      Arr   = reproj_mod2obs(Xint,Yint,Zint,Xref,Yref,mask=1*Zref.mask)

      # add the mask for Arr to Zref
      Zref  = np.ma.array(Zref.data,mask=Arr.mask)

      MPdict   = {'Over':{},'Under':{}}
      tfiles   = {'Over':{},'Under':{}}

      if 0:
         # test interpolation and matching of masks
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,2,1)
         ax1.imshow(Arr)
         ax2   = fig.add_subplot(1,2,2)
         ax2.imshow(Zref)
         plt.show(fig)
         return

      if do_sort:
         # possible regions are:
         regions  = ['gre','bar','beau','lab','balt','les','can']

         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt' )
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')

         # for reg in ['bar']:
         for reg in regions:

            # Arr,Zref are np.ma.array objects
            Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref,region=reg)
            MPdict['Over'] .update({reg:Over})
            MPdict['Under'].update({reg:Under})

            for OU in ['Over','Under']:

               fname0   = self.basename+'_v'+obs_type +'_'+OU+'_'+reg
               tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
               if 'all' in tfile.keys():
                  tfiles[OU].update({reg:tfile['all']})

         if 0:
            MPdict['Over'] [reg].show_maps()
            MPdict['Under'][reg].show_maps()
            return MPdict
      else:
         reg         = 'all'
         Over,Under  = mc.get_AOD_polys(Arr.values,Zref.values,lon_ref,lat_ref)
         MPdict['Over'] .update({reg:Over})
         MPdict['Under'].update({reg:Under})

         for OU in ['Over','Under']:

            fname0   = self.basename+'_v'+obs_type+'_'+OU+'_'+reg
            tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
            if 'all' in tfile.keys():
               tfiles[OU].update({reg:tfile['all']})

      print(tfiles)
      print(MPdict)
      Pdict = {'Over':{},'Under':{}}
      for OU in ['Over','Under']:
         PLOTTING = False
         for reg in tfiles[OU].keys():

            ##########################################################
            # filenames
            tfil     = tfiles[OU][reg]                          # text file with polygon outlines characterized
            figname  = tfil.replace('.txt','.png')          # plot of polygons
            shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
            sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
            ##########################################################


            ##########################################################
            if do_sort:
               mapreg   = reg
            else:
               mapreg   = self.HYCOM_region
            ##########################################################


            ##########################################################
            # process each text file to get MIZ width etc
            print("MIZchar.single_file: "+tfil+"\n")
            bmap     = Fplt.start_HYCOM_map(mapreg)
            Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
            Pdict[OU].update({reg:Psolns})
            
            # Save summary & shapefile
            mc.save_summary  (Psolns,sumname)
            mc.save_shapefile(Psolns,filename=shpname)
            ##########################################################

            
            if plotting:
               ##########################################################
               # Make plot
               var_opts = make_plot_options('fice',ice_mask=True)
               pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=[0,1])[0]
               fig      = pobj.fig
               ax       = pobj.ax
               PLOTTING = True

               for MIZi in Psolns:
                  # plot outlines of polygons
                  lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
                  bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)

                  Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
                  if Wavg>26:
                     MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)

                     # add text with mean width
                     xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
                     xav                  = (xmin+xmax)/2.
                     ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
                        color='k',fontsize=16,horizontalalignment='right',\
                        verticalalignment='top')

               Fplt.finish_map(bmap)
               print('Saving '+figname)
               fig.savefig(figname)
               # plt.show(fig)
               ax.cla()
               fig.clear()
               # finished region
               ##########################################################

         if PLOTTING:
            plt.close(fig)

      return MPdict,tfiles,Pdict
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs(self,pobj=None,bmap=None,\
         obs_type='OSISAF',date_label=1,figname=None,**kwargs):

      var_opts1   = make_plot_options('ficem',\
         vec_opt=0,conv_fac=1,wave_mask=False,ice_mask=True,dir_from=True)
      var_opts1   = check_var_opts(var_opts1,self.variables)

      if 'show' in kwargs:
         show           = kwargs['show']
         kwargs['show'] = False
         pobj,bmap      = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,**kwargs)
      else:
         show        = True
         pobj,bmap   = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,show=False,**kwargs)

      fig,ax,cbar = pobj.get()

      if obs_type=='OSISAF':
         obsfil   = '/work/shared/nersc/msc/OSI-SAF/'+\
               self.datetime.strftime('%Y')+'_nh_polstere/'+\
               'ice_conc_nh_polstere-100_multi_'+\
               self.datetime.strftime('%Y%m%d')+'1200.nc'
      else:
         raise ValueError('invalid value of obs_type: '+obs_type)

      print(obsfil)
      obs      = nc_getinfo(obsfil)
      fice     = obs.get_var('ice_conc')
      lon,lat  = obs.get_lonlat()
      bmap.contour(lon,lat,fice.values[:,:],[15],colors='g',\
            linewidths=2,ax=ax,latlon=True)

      dtmo     = self.datetimes[0]
      if self.HYCOM_region=='TP4':
         xyann = (0.05,.925)
      else:
         xyann = (0.4,.925)

      if date_label==1:
         tlabel   = dtmo.strftime('%d %b %Y')
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
      elif date_label==2:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         pobj.ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)

      if figname is not None:
         fig.savefig(figname)

      if show:
         # fig.show()
         plt.show(fig)

      return pobj,bmap
   ###########################################################


   ###########################################################
   def plot_var_pair(self,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs):

      # ====================================================================
      # check options
      var_opts1   = check_var_opts(var_opts1,self.all_variables)
      var_opts2   = check_var_opts(var_opts2,self.all_variables)
      check_pair(var_opts1,var_opts2)
      # ====================================================================

      pobj,bmap   = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,**kwargs)
      self.plot_var(var_opts2,pobj=pobj,bmap=bmap,**kwargs)

      return pobj,bmap
   ###########################################################


   ###########################################################
   def make_png(self,var_opts,pobj=None,bmap=None,figdir='.',date_label=2,**kwargs):

      var_opts = check_var_opts(var_opts,self.all_variables)

      new_fig  = (pobj is None)
      if new_fig:
         pobj  = plot_object()

      if 'show' in kwargs:
         show           = kwargs['show']
         kwargs['show'] = False
         pobj,bmap      = self.plot_var(var_opts,pobj=pobj,bmap=bmap,**kwargs)
      else:
         show        = False
         pobj,bmap   = self.plot_var(var_opts,pobj=pobj,bmap=bmap,show=False,**kwargs)

      dtmo     = self.datetimes[0]
      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')

      if self.HYCOM_region=='TP4':
         xyann = (0.05,.925)
      else:
         xyann = (0.4,.925)

      if date_label==1:
         tlabel   = dtmo.strftime('%d %b %Y')
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
      elif date_label==2:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)

      if pobj.axpos is not None:
         # try to make sure axes don't move round
         pobj.ax.set_position(pobj.axpos)

      vname    = var_opts.name
      Fname    = vname
      vec_opt  = var_opts.vec_opt

      if vec_opt==1:
         #magnitude only
         if vname in ['u','usurf']:
            Fname = 'surf_speed'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_speed'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_magnitude'

      elif vec_opt==2 or vec_opt==3:
         #quiver plots on top of magnitude or by itself
         if vname in ['u','usurf']:
            Fname = 'surf_vel'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress'

      elif vec_opt==4:
         #direction as scalar
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      elif vec_opt==5:
         #direction -> vector
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      Fname    = Fname.strip('_')
      figname  = figdir+'/'+self.basename+'_'+Fname+datestr+'.png'

      print('Saving to '+figname) 
      pobj.fig.savefig(figname)

      if new_fig:
         pobj.ax.cla()
         pobj.fig.clear()
         plt.close(pobj.fig)

      return pobj,bmap
   ###########################################################


   ###########################################################
   def make_png_pair(self,var_opts1,var_opts2,\
         pobj=None,bmap=None,figdir='.',date_label=2,**kwargs):

      # ====================================================================
      # check options
      var_opts1   = check_var_opts(var_opts1,self.all_variables)
      var_opts2   = check_var_opts(var_opts2,self.all_variables)
      check_pair(var_opts1,var_opts2)
      # ====================================================================

      new_fig  = (pobj is None)
      if new_fig:
         pobj  = plot_object()

      if 'show' in kwargs:
         show           = kwargs['show']
         kwargs['show'] = False
         pobj,bmap      = self.plot_var_pair(var_opts1,var_opts2,\
               pobj=pobj,bmap=bmap,**kwargs)
      else:
         show        = False
         pobj,bmap   = self.plot_var_pair(var_opts1,var_opts2,\
               pobj=pobj,bmap=bmap,show=False,**kwargs)

      fig,ax,cbar = pobj.get()

      dtmo     = self.datetimes[0]
      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
      if self.HYCOM_region=='TP4':
         xyann = (0.05,.925)
      else:
         xyann = (0.4,.925)

      if date_label==1:
         tlabel   = dtmo.strftime('%d %b %Y')
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
      elif date_label==2:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)

      # set name with 1st variable only
      vname    = var_opts1.name
      Fname    = var_opts1.name
      vec_opt  = var_opts1.vec_opt
      if vec_opt==1:
         #magnitude only
         if vname in ['u','usurf']:
            Fname = 'surf_speed'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_speed'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_magnitude'

      elif vec_opt==4:
         #direction as scalar
         if vname in ['u','usurf']:
            Fname = 'surf_current_dirn'
         elif 'u' in vname:
            Fname = vname.strip('u')+'_vel_dirn'
         elif 'taux' in vname:
            Fname = vname.strip('taux')+'_stress_dirn'

      Fname    = Fname.strip('_')
      figname  = figdir+'/'+self.basename+'_'+Fname+datestr+'.png'

      print('Saving to '+figname) 
      fig.savefig(figname)

      if new_fig:
         ax.cla()
         fig.clear()
         plt.close(fig)

      return pobj,bmap
   ###########################################################


class file_list:
   def __init__(self,directory,pattern,extension,**kwargs):
      import os

      self.directory = directory
      self.extension = extension

      print(directory,pattern,extension)

      lst            = os.listdir(directory)
      self.file_list = []
      for fil in lst:

         fname,fext  = os.path.splitext(fil)

         # check pattern
         if (fext==extension) and (pattern in fname):
            self.file_list.append(fil)

      wsn            = '/work/shared/nersc/msc/ModelInput'
      gridpath_lut   = {'FR1':wsn+'/FramStrait_Hyc2.2.12/FR1a0.03-clean//topo',\
                        'BS1':wsn+'/BS1a0.045-clean/topo',\
                        'TP4':wsn+'/../REANALYSIS/topo'}
      HB = False
      if extension=='.a':
         self.getinfo      = HYCOM_binary_info
         HB                = True
         self.HYCOM_region = self.file_list[0][:3]
         if 'gridpath' not in kwargs:
            kwargs.update({'gridpath':gridpath_lut[self.HYCOM_region]})
         
      elif extension=='.nc':
         self.getinfo      = nc_getinfo
         self.HYCOM_region = 'TP4' #TODO pass in HYCOMreg?

      # get main objects
      objects     = []
      datetimes   = []
      time_values = []
      for fil in self.file_list:
         obj   = self.getinfo(self.directory+'/'+fil,**kwargs)
         objects.append(obj)
         datetimes.append(obj.datetimes[0])
         time_values.append(obj.time_values[0])

      self.reference_date  = obj.reference_date
      self.time_units      = obj.time_units

      # sort the time values (1st record is earliest)
      TV                   = sorted([(e,i) for i,e in enumerate(time_values)])
      self.time_values,ii  = np.array(TV).transpose()
      self.objects         = [objects  [int(i)] for i in ii]
      self.datetimes       = [datetimes[int(i)] for i in ii]

      # add some methods inherited from individual objects
      self.get_lonlat   = self.objects[0].get_lonlat
      if HB:
         self.get_depths   = self.objects[0].get_depths

      return
   ###########################################################


   ###########################################################
   def make_png_all(self,var_opts,HYCOMreg=None,figdir='.',**kwargs):

      pobj        = plot_object()
      fig,ax,cbar = pobj.get()

      if HYCOMreg is None:
         HYCOMreg = self.HYCOM_region
      bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      N  = len(self.objects)
      for i,obj in enumerate(self.objects):

         pobj,bmap   = obj.make_png(var_opts,\
                           bmap=bmap,\
                           figdir=figdir,show=False,**kwargs)

         ax.cla()
         if pobj.cbar is not None:
            pobj.cbar.ax.clear()   # cbar.ax.clear()

         print('\n'+str(i+1)+' records done out of '+str(N))

      plt.close(fig)
      return
   ###########################################################


   ###########################################################
   def make_png_pair_all(self,var_opts1,var_opts2,HYCOMreg=None,figdir='.',**kwargs):

      pobj        = plot_object()
      fig,ax,cbar = pobj.get()

      if HYCOMreg is None:
         HYCOMreg = self.HYCOM_region
      bmap        = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      N  = len(self.objects)
      for i,obj in enumerate(self.objects):

         print('\n'+str(i)+' records done out of '+str(N))

         pobj,bmap   = obj.make_png_pair(var_opts1,var_opts2,\
                        pobj=pobj,bmap=bmap,\
                        figdir=figdir,show=False,**kwargs)

         if i==0:
            # Fix axes position to stop it moving round
            pobj  = pobj.renew(axpos=pobj.ax.get_position())

         ax.cla()
         if pobj.cbar is not None:
            pobj.cbar.ax.clear()   # cbar.ax.clear()

      plt.close(fig)
      return
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs_all(self,HYCOMreg=None,figdir='.',**kwargs):

      pobj        = plot_object()
      fig,ax,cbar = pobj.get()

      if HYCOMreg is None:
         HYCOMreg = self.HYCOM_region
      bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

      if 'show' in kwargs:
         kwargs['show'] = False
      else:
         kwargs.update({'show':False})

      if 'obs_type' not in kwargs:
         obs_type = 'OSISAF'
         kwargs.update({'obs_type':obs_type})

      N  = len(self.objects)
      for i,obj in enumerate(self.objects):

         dtmo     = obj.datetime
         datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
         figname  = figdir+'/'+obj.basename+'_v'+obs_type+'_'+datestr+'.png'
         pobj,bmap   = obj.compare_ice_edge_obs(pobj=pobj,bmap=bmap,\
               figname=figname,**kwargs)

         ax.cla()
         if pobj.cbar is not None:
            pobj.cbar.ax.clear()   # cbar.ax.clear()

         print('\n'+str(i+1)+' records done out of '+str(N))

      plt.close(fig)
      return
   ###########################################################

######################################################################
