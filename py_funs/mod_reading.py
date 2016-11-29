import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen
import fns_plotting as Fplt
from scipy.interpolate import griddata as grd
import os,sys
import shapely.geometry as shg
import geometry_sphere as GS


##########################################################
def basemap_OSISAF():
    bmap = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891),\
               projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
    return bmap
##########################################################


##########################################################
class AOD_output:
  def __init__(self,summary_files,tfiles,shapefiles,types,regions,dto,anomaly_file=None):
     self.summary_files      = summary_files
     self.shapefiles         = shapefiles
     self.text_files         = tfiles
     self.regions_analysed   = regions
     self.datetime           = dto
     self.types              = types
     
     if anomaly_file is not None:
        anom_fil  = os.path.splitext(anomaly_file)[0]
        #
        self.anomaly_file_npz = anom_fil+'.npz'
        self.anomaly_file_txt = anom_fil+'.txt'

     return
##########################################################


###############################################
class time_series:

   ############################################
   def __init__(self,dates,data,units=None,filename=None,overwrite=False):
      self.dates                 = dates
      self.data                  = data
      self.units                 = units
      self.variables             = data.keys()
      self.number_of_dates       = len(dates)
      self.number_of_variables   = len(data.keys())
      
      if not os.path.exists(filename) or overwrite:
         print('Saving time series to '+filename)
         self.write_file(filename)

      self.filename  = os.path.abspath(filename)
      return
   ############################################


   ###########################################################
   def nearestDate(self, pivot):
      """
      dto,time_index = self.nearestDate(dto0)
      dto0  = datetime.datetime objects
      dto   = datetime.datetime objects - nearest value in self.datetimes to dto0
      time_index: dto=self.dates[time_index]
      """
      dto         = min(self.dates, key=lambda x: abs(x - pivot))
      time_index  = self.dates.index(dto)
      return dto,time_index
   ###########################################################


   ###########################################################
   def max(self, vname):
      """
      vmax,time_index = self.max(vname)
      *vname is in self.data.keys():
      * vmax=self.data[vname].max()
      * self.data[vname][time_index]=vmax
      """
      vmax        = self.data[vname].max()
      time_index  = list(self.data[vname]).index(vmax)
      return vmax,time_index
   ###########################################################


   ###########################################################
   def min(self, vname):
      """
      vmin,time_index = self.min(vname)
      *vname is in self.data.keys():
      * vmin=self.data[vname].min()
      * self.data[vname][time_index]=vmin
      """
      vmin        = self.data[vname].min()
      time_index  = list(self.data[vname]).index(vmin)
      return vmin,time_index
   ###########################################################


   ############################################
   def plot(self,var_name,refdate=None,time_units='days',yscaling=1.,pobj=None,**kwargs):
      if pobj is None:
         pobj  = plot_object()

      if time_units=='days':
         xfac  = 24*3600. # seconds in 1 day
      elif time_units=='hours':
         xfac  = 3600. # seconds in 1h
      elif time_units=='minutes':
         xfac  = 60. # seconds in 1min
      else:
         xfac  = 1. # seconds in 1min
      if refdate is None:
         refdate  = self.dates[0]

      x  = np.array([(dt-refdate).total_seconds()/xfac for dt in self.dates])
      y  = self.data[var_name]*yscaling

      lin  ,= pobj.ax.plot(x,y,**kwargs)

      info  = {'time_data':x,'refdate':refdate,'time_units':'days'}
      return pobj,lin,info
   ############################################


   ############################################
   def write_file(self,filename):
      fid   = open(filename,'w')
      blk   = 3*' '

      # =======================================
      # header
      ss          = 'Date'+blk
      add_units   = (self.units is not None)
      for i,vbl in enumerate(self.variables):
         su = ''
         if add_units:
            unit  = self.units[vbl]
            if unit!=su:
               su = ','+unit
         if i<self.number_of_variables-1:
            ss   += vbl+su+blk
         else:
            ss   += vbl+su+'\n'
      fid.write(ss)
      # =======================================


      # ===================================================
      # loop over dates and write variables for each date
      for n,dt in enumerate(self.dates):
         ss = dt.strftime('%Y%m%dT%H%M%SZ')+blk
         for i,vbl in enumerate(self.variables):
            if i<self.number_of_variables-1:
               ss   += str(self.data[vbl][n])+blk
            else:
               if n<self.number_of_dates-1:
                  ss   += str(self.data[vbl][n])+'\n'
               else:
                  ss   += str(self.data[vbl][n])
         fid.write(ss)
      # ===================================================

      fid.close()
      return
   #########################################

############################################

############################################
def read_time_series(tfil):

   fid   = open(tfil)
   lines = fid.readlines()
   fid.close()
   
   #########################################
   # read header
   lin      = lines[0]
   Vnames   = lin.split()[1:] #variable names : eg MIZ_width,m (1st col is date)
   lines.remove(lin)
   
   # get variable names and units
   data     = {}
   units    = {}
   vnames   = []
   for v in Vnames:
      sp    = v.split(',')
      vname = sp[0]
      vnames.append(vname)
      data.update({vname:[]})

      # add units
      unit  = ''
      if len(sp)>1:
         unit  = sp[1]
      units.update({vname:unit})
   #########################################

   #########################################
   # read dates and data
   dates = []
   for lin in lines:
      ss    = lin.split()
      cdate = ss[0]

      if len(cdate)==8:
         cdate+= 'T120000Z'

      dates.append(datetime.strptime(cdate,'%Y%m%dT%H%M%SZ'))

      for i,vname in enumerate(vnames):
         data[vname].append(float(ss[i+1]))
   #########################################

   #########################################
   # return time_series object
   for vname in vnames:
      data[vname] = np.array(data[vname])
   if len(units)==0:
      units = None
   return time_series(dates,data,units=units,filename=tfil)
############################################

##########################################################
class read_MIZpoly_summary:

   #######################################################
   def __init__(self,tfil,cdate=None,ctime=None):
      
      ########################################################
      self.info   = {'filename'    :tfil,\
                     'datetime'    :None,\
                     'time_in_days':None,\
                     'datetime_ref':datetime(1901,1,1)}
      if (cdate is not None) and (ctime is not None):
         self.info['datetime']   = datetime.strptime(cdate+' '+ctime,'%Y%m%d %H%M%S')
      elif (cdate is not None):
         self.info['datetime']   = datetime.strptime(cdate,'%Y%m%d')
      ########################################################

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

            if self.info.datetime is not None:
               if self.info['datetime'] != dtm:
                  raise ValueError('date in '+tfil+' not consistent with inputs cdate,ctime')
            else:
               self.info['datetime'] = dtm
         else:
            val   = float(sp[1].strip())
            setattr(self,att,val)

      if self.info['datetime'] is not None:
         dt                         = self.info['datetime']-self.info['datetime_ref']
         self.info['time_in_days']  = dt.total_seconds()/float(24*60*60)

      return
   #######################################################

##########################################################


############################################################################
# Function that reprojects model into observational grid
def reproj_mod2obs(X1,Y1,Z1,X2,Y2,method='linear',mask=None):
   # input coords from X1,Y1; Z1 is array to interp; X2,Y2 are output matrices


   # ========================================
   # determine resolution of source grid
   X1d   = X1.reshape(X1.size)
   Y1d   = Y1.reshape(Y1.size)
   Z1d   = Z1.reshape(Z1.size).data.astype(float)
   dx    = abs(X1d[1]-X1d[0])
   dy    = abs(Y1d[1]-Y1d[0])
   res   = np.sqrt(dx**2+dy**2)
   Z1d[Z1.reshape(Z1.size).mask] = np.nan
   # ========================================


   # ========================================
   # reduce size of source grid
   xmin     = np.min(X2)-2*res
   ymin     = np.min(Y2)-2*res
   xmax     = np.max(X2)+2*res
   ymax     = np.max(Y2)+2*res
   needed   = np.logical_and(X1d>=xmin,X1d<=xmax)
   needed   = np.logical_and(needed,Y1d<=ymax)
   needed   = np.logical_and(needed,Y1d>=ymin)
   X1d      = X1d[needed]
   Y1d      = Y1d[needed]
   Z1d      = Z1d[needed]
   # ========================================


   # ========================================
   # Interpolation 
   # - can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
   C     = np.array([X1d,Y1d]).T
   Z2    = grd(C,Z1d,(X2,Y2),method=method)
   # ========================================

   # check output and add mask
   mask2 = np.isnan(Z2)
   if mask is not None:
      # apply union of mask and model nans
      mask2=np.logical_or(mask2,mask)
    
   Z2 = np.ma.array(Z2,mask=mask2)
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
   lists.append(['ficem','fice','ice_conc','icec',\
                  'concentration','sea_ice_concentration',\
                  'ice fraction'])

   # ice thick alt names
   lists.append(['hicem','hice','ice_thick','icetk',\
                  'sea_ice_thickness','thickness','sea_ice_concentration'])

   # floe size alt names
   lists.append(['dfloe','dmax','Dfloe','Dmax'])

   # H_s alt names
   lists.append(['swh','Hs','hs'])

   for names in lists:
      if vname in names:
         for vbl in names:
            if vbl in variables:
               return vbl

   raise ValueError(vname+' not in variable list')
   return
###########################################################


###########################################################
def check_var_opts(var_opts,variables=None):
   """
   var_opts = check_var_opts(var_opts,variables=None)
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

   if variables is not None:
      vname       = check_names(var_opts.name,variables)
      var_opts    = new_var_opts(var_opts,vname)

   return var_opts
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
def nc_getinfo(ncfil,time_index=None,lonlat_file=None):
   import mod_netcdf_utils as MNU
   nci   = MNU.nc_getinfo(ncfil,time_index=time_index,\
            lonlat_file=lonlat_file)
   return nci
##########################################################


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


######################################################################
def HYCOM_binary_info(fname,gridpath=None):
   import mod_HYCOM_utils as MHU
   hbi   = MHU.HYCOM_binary_info(fname,gridpath=gridpath)
   return hbi
######################################################################


######################################################################
def GetVar(fobj,vname,layer=0,time_index=0):

   if fobj.filetype=='HYCOM_binary':
      if layer==0:
         vbl   = fobj.get_var(vname)
      else:
         vbl   = fobj.get_var([vname,layer])
   elif fobj.filetype=='netcdf':
      vbl   = fobj.get_var(vname,time_index=time_index)

   return vbl
######################################################################

#######################################################################
def imshow(fobj,var_opts,pobj=None,\
      clim=None,add_cbar=True,clabel=None,show=True,\
      test_ijs=None,time_index=0):
   """
   pobj   = imshow(fobj,var_opts,time_index=0,pobj=None,\
        clim=None,add_cbar=True,clabel=None,show=True,\
        test_ijs=None)
   """

   from mpl_toolkits.basemap import Basemap
   from matplotlib import cm

   var_opts    = check_var_opts(var_opts,fobj.all_variables)
   vname       = var_opts.name
   layer       = var_opts.layer
   vec_opt     = var_opts.vec_opt
   conv_fac    = var_opts.conv_fac
   ice_mask    = var_opts.ice_mask
   wave_mask   = var_opts.wave_mask
   dir_from    = var_opts.dir_from
   if vname not in fobj.all_variables:
      raise ValueError('Variable '+vname+'not in '+fobj.afile)

   if pobj is None:
      pobj  = plot_object()
   fig,ax,cbar = pobj.get()

   if clim is not None:
      vmin,vmax   = clim
   else:
      vmin  = None
      vmax  = None

   vbl   = GetVar(fobj,vname,layer=layer,time_index=time_index)
   # if fobj.filetype=='HYCOM_binary':
   #    if layer==0:
   #       vbl   = fobj.get_var(vname)
   #    else:
   #       vbl   = fobj.get_var([vname,layer])
   # elif fobj.filetype=='netcdf':
   #    vbl   = fobj.get_var(vname,time_index=time_index)

   mask     = vbl.values.mask

   ################################################################## 
   # add ice or wave masks
   if ice_mask and wave_mask:
      fice        = GetVar(fobj,'fice',layer=0,time_index=time_index)
      # fice        = fobj.get_var('ficem')
      good        = np.array(1-fice.values.mask  ,dtype='bool')
      mask1       = np.zeros(fice.shape         ,dtype='bool')
      mask1[good] = (fice[good]<.01)          # water
      mask1       = np.logical_or(mask,mask1) # water or NaN
      #
      Hs          = GetVar(fobj,'swh',layer=0,time_index=time_index)
      # Hs          = fobj.get_var('swh')
      good        = np.array(1-Hs.values.mask,dtype='bool')
      mask2       = np.zeros(Hs.shape        ,dtype='bool')
      mask2[good] = (Hs[good]<.01)              # no waves
      mask2       = np.logical_or(mask,mask2)   # no waves or NaN
      #
      mask  = np.logical_or(mask1,mask2)

   elif ice_mask:
      fice        = GetVar(fobj,'fice',layer=0,time_index=time_index)
      # fice        = fobj.get_var('ficem')
      good        = np.array(1-fice.values.mask ,dtype='bool')
      mask1       = np.zeros(fice.shape         ,dtype='bool')
      mask1[good] = (fice[good]<.01)          # water
      mask        = np.logical_or(mask,mask1) # water or NaN

   elif wave_mask:
      Hs          = GetVar(fobj,'swh',layer=0,time_index=time_index)
      # Hs          = fobj.get_var('swh')
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

      vbl2  = GetVar(fobj,vname2,layer=layer,time_index=time_index)
      # if layer==0:
      #    vbl2  = fobj.get_var(vname2)
      # else:
      #    vbl2  = fobj.get_var([vname2,layer])

      Marr  = np.hypot(vbl.values.data,vbl2.values.data)
      Marr  = np.ma.array(conv_fac*Marr,mask=mask)

   elif (vec_opt==2) or (vec_opt==3):
      # 2: plot vector magnitude + direction (unit vectors)
      # 3: plot vector direction only
      raise ValueError('vec_opt==2,3 disabled for imshow')

   elif vec_opt==4:

      # plot direction as scalar
      U  = None

      if vname[0]=='u':
         vname2   = 'v'+vname[1:]
      elif vname[:4]=='taux':
         vname2   = 'tauy'+vname[4:]

      vbl2  = GetVar(fobj,vname2,layer=layer,time_index=time_index)
      # if layer==0:
      #    vbl2  = fobj.get_var(vname2)
      # else:
      #    vbl2  = fobj.get_var([vname2,layer])
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

   if fobj.filetype=='netcdf':
      PC = ax.imshow(Marr,origin='lower',vmin=vmin,vmax=vmax)
   else:
      PC = ax.imshow(Marr.transpose(),origin='lower',vmin=vmin,vmax=vmax)

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
   if test_ijs is not None:
      for itst,jtst in test_ijs:
         ax.plot(jtst,itst,'^m',markersize=5)
   ################################################################## 

   if show:
      # fig.show()
      plt.show(fig)

   return pobj
###########################################################


###########################################################
def interp2points(fobj,varname,target_lonlats,time_index=0,mapping=None,**kwargs):

   lon,lat  = fobj.get_lonlat()

   # ===============================================================================
   # do interpolation in stereographic projection
   if mapping is None:
      import pyproj
      souths   = lat[lat<0]
      lat_0    = 90.
      lon_0    = -45.
      if len(souths)==len(lat):
         # only use south pole projection if all points in southern hemisphere
         lat_0    = -90.

      srs   =  '+proj=stere'     \
               +' +lon_0='+str(lon_0)  \
               +' +lat_0='+str(lat_0)  \
               +' +lat_ts='+str(lat_0) \
               +' +ellps=WGS84'

      mapping  = pyproj.Proj(srs)
   # ===============================================================================

   #source
   X,Y   = mapping(lon,lat)
   Z     = fobj.get_var(varname,time_index=time_index).values #numpy masked array

   #target
   lons,lats   = target_lonlats
   x,y         = mapping(lons,lats)

   # do interpolation
   fvals = reproj_mod2obs(X,Y,Z,x,y,**kwargs) #numpy masked array

   return fvals
#######################################################################


#######################################################################
def plot_var(fobj,var_opts,time_index=0,\
      pobj=None,bmap=None,HYCOMreg=None,\
      date_label=0,date_color='k',\
      smoothing=0,\
      clim=None,add_cbar=True,clabel=None,show=True,\
      test_lonlats=None):

   from mpl_toolkits.basemap import Basemap
   from matplotlib import cm

   var_opts    = check_var_opts(var_opts,fobj.all_variables)
   vname       = var_opts.name
   layer       = var_opts.layer
   vec_opt     = var_opts.vec_opt
   conv_fac    = var_opts.conv_fac
   ice_mask    = var_opts.ice_mask
   wave_mask   = var_opts.wave_mask
   dir_from    = var_opts.dir_from
   if vname not in fobj.all_variables:
      raise ValueError('Variable '+vname+'not in '+fobj.afile)

   if pobj is None:
      pobj  = plot_object()
   fig,ax,cbar = pobj.get()

   if HYCOMreg is None:
      HYCOMreg = fobj.HYCOM_region
      if HYCOMreg is None:
         HYCOMreg = 'TP4'

   if bmap is None:
      # make basemap
      bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')

   # lon,lat  = fobj.get_lonlat()
   lon,lat  = fobj.get_fixed_lonlat(bmap)

   if clim is not None:
      vmin,vmax   = clim
   else:
      vmin  = None
      vmax  = None

   vbl   = GetVar(fobj,vname,layer=layer,time_index=time_index)
   mask  = vbl.values.mask

   ################################################################## 
   # add ice or wave masks
   if ice_mask and wave_mask:
      # fice        = fobj.get_var('ficem')
      fice        = GetVar(fobj,'fice',layer=0,time_index=time_index)
      good        = np.array(1-fice.values.mask  ,dtype='bool')
      mask1       = np.zeros(fice.shape         ,dtype='bool')
      mask1[good] = (fice[good]<.01)          # water
      mask1       = np.logical_or(mask,mask1) # water or NaN
      #
      # Hs          = fobj.get_var('swh')
      Hs          = GetVar(fobj,'swh',layer=0,time_index=time_index)
      good        = np.array(1-Hs.values.mask,dtype='bool')
      mask2       = np.zeros(Hs.shape        ,dtype='bool')
      mask2[good] = (Hs[good]<.01)              # no waves
      mask2       = np.logical_or(mask,mask2)   # no waves or NaN
      #
      mask  = np.logical_or(mask1,mask2)

   elif ice_mask:
      fice        = GetVar(fobj,'fice',layer=0,time_index=time_index)
      # fice        = fobj.get_var('ficem')
      good        = np.array(1-fice.values.mask ,dtype='bool')
      mask1       = np.zeros(fice.shape         ,dtype='bool')
      mask1[good] = (fice[good]<.01)          # water
      mask        = np.logical_or(mask,mask1) # water or NaN

   elif wave_mask:
      # Hs          = fobj.get_var('swh')
      Hs          = GetVar(fobj,'swh',layer=0,time_index=time_index)
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

      vbl2  = GetVar(fobj,vname,layer=layer,time_index=time_index)
      Marr  = np.hypot(vbl.values.data,vbl2.values.data)
      Marr  = np.ma.array(conv_fac*Marr,mask=mask)

   elif vec_opt==2:

      # plot vector magnitude + direction (unit vectors)
      if vname[0]=='u':
         vname2   = 'v'+vname[1:]
      elif vname[:4]=='taux':
         vname2   = 'tauy'+vname[4:]

      vbl2  = GetVar(fobj,vname,layer=layer,time_index=time_index)
      U     = conv_fac*vbl.values.data
      V     = conv_fac*vbl2.values.data

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

      vbl2  = GetVar(fobj,vname,layer=layer,time_index=time_index)
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

      vbl2  = GetVar(fobj,vname,layer=layer,time_index=time_index)
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
      if smoothing>0:
         import scipy.ndimage as NDI

         # smooth by avg over neighbours
         if smoothing==1:
            kernel   = np.array([[1, 1, 1],
                                 [1, 0, 1],
                                 [1, 1, 1]],dtype='float')
         elif smoothing==2:
            kernel   = np.array([[1, 1, 1, 1],
                                 [1, 0, 0, 1],
                                 [1, 0, 0, 1],
                                 [1, 1, 1, 1]],dtype='float')

         DT                = NDI.convolve(Marr.data, kernel/kernel.sum(), mode='constant', cval=0.0)
         kernel[kernel==0] = 1.
         MSK               = NDI.convolve(np.array(Marr.mask,dtype='float'), kernel/kernel.sum(), mode='constant', cval=1.0) 
         MSK[MSK>0]        = 1.

         Marr  = np.ma.array(DT,mask=np.array(MSK,dtype='bool'))
         del DT,MSK

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

         if 1:
            # old way - works for multiple plots
            if cbar is None:
               cbar  = fig.colorbar(PC)
            else:
               cbar  = fig.colorbar(PC,cax=cbar.ax)
         else:
            # new way - doesn't work for multiple plots
            # - worked for side-by-side plots
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider  = make_axes_locatable(ax)
            cax      = divider.append_axes("right", size="5%", pad=0.15)
            cbar     = fig.colorbar(PC,cax=cax)

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
         bmap.plot(lont,latt,'^m',markersize=10,latlon=True,ax=ax)

   Fplt.finish_map(bmap,ax=ax)
   
   # date label
   if (fobj.HYCOM_region=='TP4'):
      xyann = (0.05,.925)
   else:
      xyann = (0.4,.925)
   # xyann = (0.60,.932)

   if type(date_label)==type('hi'):
      tlabel   = date_label
      pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',\
            fontsize=18,color=date_color)
            # horizontalalignment='center',\
            # fontsize=18,color=date_color)

   elif fobj.datetimes is not None:
      dtmo     = fobj.datetimes[time_index]
      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
      tlabel   = None
      if date_label==1:
         tlabel   = dtmo.strftime('%d %b %Y')
      elif date_label==2:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')

      if tlabel is not None:
         pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',\
            fontsize=18,color=date_color)
            # horizontalalignment='center',\
            # fontsize=28,color=date_color)

   if show:
      # fig.show()
      plt.show(fig)

   return pobj,bmap
###########################################################



###########################################################
def plot_var_pair(fobj,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs):

   # ====================================================================
   # check options
   var_opts1   = check_var_opts(var_opts1,fobj.all_variables)
   var_opts2   = check_var_opts(var_opts2,fobj.all_variables)
   check_pair(var_opts1,var_opts2)
   # ====================================================================

   if 'show' in kwargs:
      show  = kwargs['show']
      del kwargs['show']
   else:
      # default is show:
      show  = True
   pobj,bmap   = plot_var(fobj,var_opts1,pobj=pobj,bmap=bmap,show=False,**kwargs)

   if 'date_label' in kwargs:
      del kwargs['date_label']
   plot_var(fobj,var_opts2,pobj=pobj,bmap=bmap,show=show,**kwargs)

   return pobj,bmap
###########################################################


###########################################################
def make_png(fobj,var_opts,pobj=None,bmap=None,figdir='.',date_label=2,**kwargs):

   var_opts = check_var_opts(var_opts,fobj.all_variables)

   new_fig  = (pobj is None)
   if new_fig:
      pobj  = plot_object()

   if 'time_index' not in kwargs:
      time_index  = 0
      kwargs.update({'time_index':time_index})
   else:
      time_index  = kwargs['time_index']

   # pass on date_label to plot_var
   kwargs.update({'date_label':date_label})

   if 'show' in kwargs:
      show           = kwargs['show']
      kwargs['show'] = False
      pobj,bmap      = fobj.plot_var(var_opts,pobj=pobj,bmap=bmap,**kwargs)
   else:
      show        = False
      pobj,bmap   = fobj.plot_var(var_opts,pobj=pobj,bmap=bmap,show=False,**kwargs)

   if pobj.axpos is not None:
      # try to make sure axes don't move round
      pobj.ax.set_position(pobj.axpos)

   vname    = var_opts.name
   Fname    = vname
   vec_opt  = var_opts.vec_opt

   # =======================================================
   # set figure name
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
   if fobj.datetimes is None:
      figname  = figdir+'/'+fobj.basename+'_'+Fname+'.png'
   else:
      datestr  = fobj.datetimes[time_index].strftime('%Y%m%dT%H%M%SZ')
      figname  = figdir+'/'+fobj.basename+'_'+Fname+datestr+'.png'

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   print('Saving to '+figname) 
   pobj.fig.savefig(figname)
   # =======================================================

   if new_fig:
      pobj.ax.cla()
      pobj.fig.clear()
      plt.close(pobj.fig)

   return pobj,bmap
###########################################################


###########################################################
def make_png_pair(fobj,var_opts1,var_opts2,\
      pobj=None,bmap=None,figdir='.',date_label=2,**kwargs):

   # ====================================================================
   # check options
   var_opts1   = check_var_opts(var_opts1,fobj.all_variables)
   var_opts2   = check_var_opts(var_opts2,fobj.all_variables)
   check_pair(var_opts1,var_opts2)
   # ====================================================================

   new_fig  = (pobj is None)
   if new_fig:
      pobj  = plot_object()

   if 'time_index' not in kwargs:
      time_index  = 0
      kwargs.update({'time_index':time_index})
   else:
      time_index  = kwargs['time_index']

   # pass on date_label to plot_var
   kwargs.update({'date_label':date_label})

   if 'show' in kwargs:
      show           = kwargs['show']
      kwargs['show'] = False
      pobj,bmap      = plot_var_pair(fobj,var_opts1,var_opts2,\
            pobj=pobj,bmap=bmap,**kwargs)
   else:
      show        = False
      pobj,bmap   = plot_var_pair(fobj,var_opts1,var_opts2,\
            pobj=pobj,bmap=bmap,show=False,**kwargs)

   fig,ax,cbar = pobj.get()


   # ==================================================================
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
   if fobj.datetimes is None:
      figname  = figdir+'/'+fobj.basename+'_'+Fname+'.png'
   else:
      datestr  = fobj.datetimes[time_index].strftime('%Y%m%dT%H%M%SZ')
      figname  = figdir+'/'+fobj.basename+'_'+Fname+datestr+'.png'

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   print('Saving to '+figname) 
   fig.savefig(figname)
   # ==================================================================


   if new_fig:
      ax.cla()
      fig.clear()
      plt.close(fig)

   return pobj,bmap
   ###########################################################

###########################################################


###########################################################
def compare_ice_edge_obs(fobj,pobj=None,bmap=None,time_index=0,\
      obs_type='OSISAF',obs_path=None,obs_option='multi',\
      date_label=1,figname=None,**kwargs):
   """
   compare_ice_edge_obs(fobj,pobj=None,bmap=None,time_index=0,\
      obs_type='OSISAF',obs_path=None,obs_option='multi',\
      date_label=1,figname=None,**kwargs)
   """

   var_opts1   = make_plot_options('ficem',\
      vec_opt=0,conv_fac=1,wave_mask=False,ice_mask=True,dir_from=True)
   var_opts1   = check_var_opts(var_opts1,fobj.variables)

   if 'show' in kwargs:
      show           = kwargs['show']
      kwargs['show'] = False
      pobj,bmap      = fobj.plot_var(var_opts1,pobj=pobj,bmap=bmap,**kwargs)
   else:
      show        = True
      pobj,bmap   = fobj.plot_var(var_opts1,pobj=pobj,bmap=bmap,show=False,**kwargs)

   fig,ax,cbar = pobj.get()

   dtmo  = fobj.datetimes[time_index]
   if obs_type=='OSISAF':
      if obs_path is None:
      	 obs_path   = '/work/shared/nersc/msc/OSI-SAF/'+\
            dtmo.strftime('%Y')+'_nh_polstere/'
      if obs_option is None:
	 obs_option='multi'
      obsfil   = obs_path+\
            '/ice_conc_nh_polstere-100_'+obs_option+'_'+\
            dtmo.strftime('%Y%m%d')+'1200.nc'
   else:
      raise ValueError('invalid value of obs_type: '+obs_type)

   print(obsfil)
   obs      = nc_getinfo(obsfil)
   fice     = obs.get_var('ice_conc')
   lon,lat  = obs.get_lonlat()
   bmap.contour(lon,lat,fice.values[:,:],[15],colors='g',\
         linewidths=2,ax=ax,latlon=True)

   if 'HYCOMreg' in kwargs:
      reg   = kwargs['HYCOMreg']
   else:
      reg   = fobj.HYCOM_region
      if reg is None:
         reg   = 'TP4'
      kwargs.update({'HYCOMreg':reg})


   #############################################################
   if reg=='TP4':
      xyann = (0.05,.925)
   else:
      xyann = (0.4,.925)


   if date_label==1:
      # date only
      tlabel   = dtmo.strftime('%d %b %Y')
      pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
   elif date_label==2:
      # date + time
      tlabel   = dtmo.strftime('%d %b %Y %H:%M')
      pobj.ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)
   elif type(date_label)==type('hi'):
      # manual label
      pobj.ax.annotate(date_label,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)
   #############################################################


   if figname is not None:
      print('Saving to '+figname)
      fig.savefig(figname)

   if show:
      # fig.show()
      plt.show(fig)

   return pobj,bmap,obsfil
###########################################################


###########################################################
def compare_ice_edge_obs_all(fobj,HYCOMreg=None,figdir='.',**kwargs):

   pobj        = plot_object()
   fig,ax,cbar = pobj.get()

   # ================================================
   # set basemap for plotting
   if HYCOMreg is None:
      if hasattr(fobj,'HYCOM_region'):
         HYCOMreg = fobj.HYCOM_region
      else:
         HYCOMreg = 'Arctic'
   bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
   # ================================================

   if 'show' in kwargs:
      kwargs['show'] = False
   else:
      kwargs.update({'show':False})

   if 'obs_type' not in kwargs:
      obs_type = 'OSISAF'
      kwargs.update({'obs_type':obs_type})

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   N  = len(fobj.objects)
   for i,obj in enumerate(fobj.objects):

      dtmo              = obj.datetime
      # datestr           = dtmo.strftime('%Y%m%dT%H%M%SZ')
      datestr           = dtmo.strftime('%Y%m%d')
      figname           = figdir+'/'+obj.basename+'_v'+obs_type+'_'+datestr+'.png'
      pobj,bmap,obsfil  = obj.compare_ice_edge_obs(pobj=pobj,bmap=bmap,\
                           figname=figname,**kwargs)

      ax.cla()
      if pobj.cbar is not None:
         pobj.cbar.ax.clear()   # cbar.ax.clear()

      print('\n'+str(i+1)+' records done out of '+str(N))

   plt.close(fig)
   return
###########################################################


###########################################################
def MIZmap(fobj,var_name='dmax',time_index=0,\
      vertices=None,regions=None,no_width=False,\
      do_sort=False,EastOnly=True,plotting=True,**kwargs):
   """
   Call  : fobj.MIZmap(var_name='dmax',time_index=0,
               vertices=None,,regions=None,
               do_sort=False,EastOnly=True,plotting=True,**kwargs)

   Inputs:
      var_name is variable to find MIZ from
      - 'dmax' (<250m), 'fice' (<.8), 'swh' (>5cm), 'hice' (<45cm)
      - MIZ is determined from above condition AND (fice>.15)
      **kwargs to be passed onto MIZchar.get_MIZ_poly.write_poly_stats:
         outdir='.' - where to save results

   Returns:
      mod_reading.AOD_output object,
      which contains locations of output files, and some other info
   """

   import MIZchar as mc

   vname = check_names(var_name,fobj.variables)
   fname = check_names('fice',fobj.variables)   # need fice for ice mask
   if var_name == 'dmax':
      # FSD MIZarray(1-
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,300]# for plotting
      lower_limit = .1     # for plotting
      fice        = GetVar(fobj,fname,time_index=time_index)
   elif var_name == 'fice':
      # conc MIZ
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,1]  # for plotting
      lower_limit = .15    # for plotting
      fice        = Arr
   elif var_name == 'hice':
      # thin ice areas
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,2.] # for plotting
      lower_limit = .01    # for plotting
      fice        = GetVar(fobj,fname,time_index=time_index)
   elif var_name == 'swh':
      # waves-in-ice areas
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,4.] # for plotting
      lower_limit = .01    # for plotting
      fice        = GetVar(fobj,fname,time_index=time_index)
   else:
      raise ValueError('Wrong selection variable for MIZmap')

   print("MIZchar.get_MIZ_poly\n")
   lon,lat        = fobj.get_lonlat()
   MPdict         = {}
   tfiles         = {}
   summary_files  = {}
   shapefiles     = {}

   if (vertices is not None) and (regions is not None):
      raise ValueError('Cannot pass in both vertices and regions')
   elif vertices is not None:
      do_sort  = False
   elif regions is not None:
      do_sort  = True
   
   if do_sort:

      # ==================================================================
      # possible regions are:
      def_regions = ['gre','bar','beau','lab','balt','les','can','Antarctic']
      if regions is not None:
         # check regions are OK
         for reg in regions:
            if reg not in def_regions:
               raise ValueError('Unknown region: '+reg)
      else:

         regions  = 1*def_regions
         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt')
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')
         else:
            # forget Baltic Sea
            regions.remove('balt')
            regions.remove('Antarctic')
      # ==================================================================

      for reg in regions:
         mp = mc.get_MIZ_poly(Arr.values,lon,lat,fice.values,var_name=var_name,region=reg)
         MPdict.update({reg:mp})

         fname0   = fobj.basename+'_'+var_name +'_'+reg
         tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
         tfiles.update({reg:tfile['all']})

      if 0:
         MPdict['gre'].show_maps()
         return MPdict

   else:
      # not do_sort:
      # - do whole Arctic
      reg   = 'Arctic'
      mp    = mc.get_MIZ_poly(Arr.values,lon,lat,fice.values,var_name=var_name,vertices=vertices,region=reg)
      MPdict.update({reg:mp})
      #
      if vertices is None:
         fname0   = fobj.basename+'_'+var_name +'_'+reg
      else:
         fname0   = fobj.basename+'_'+var_name

      tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
      tfiles.update({reg:tfile['all']})

   if no_width:
      return tfiles

   Pdict    = {}
   PLOTTING = False
   for reg in tfiles.keys():

      ##########################################################
      # filenames
      tfil     = tfiles[reg]                          # text file with polygon outlines characterized
      figname  = tfil.replace('.txt','.png')          # plot of polygons
      shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
      sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
      summary_files.update({reg:sumname})
      shapefiles   .update({reg:shpname})
      ##########################################################


      ##########################################################
      # basemap for plotting
      if vertices is None:
         if do_sort:
            mapreg   = reg
         elif fobj.HYCOM_region is not None:
            mapreg   = fobj.HYCOM_region
         else:
            mapreg   = 'Arctic'
         bmap     = Fplt.start_HYCOM_map(mapreg)
      else:
         vlons,vlats = np.array(vertices).transpose()
         vx,vy       = GS.polar_stereographic_simple(vlons,vlats,NH=True,inverse=False)
         xy_verts    = [(vx[i],vy[i]) for i in range(len(vx))]

         # get approximate centre and limits for basemap
         vP          = shg.Polygon(xy_verts)
         width       = 5*np.sqrt(vP.area)
         height      = 5*np.sqrt(vP.area)
         xcen,ycen   = vP.centroid.coords.xy
         lonc,latc   = GS.polar_stereographic_simple(np.array([xcen]),np.array([ycen]),\
                        NH=True,inverse=True)
         # make basemap
         bmap        = Basemap(projection='stere',lon_0=lonc,lat_0=latc,\
                                 width=width,height=height,\
                                 resolution='i')
      ##########################################################

      # process each text file to get MIZ width etc
      print("MIZchar.single_file: "+tfil+"\n")
      # Psolns   = mc.single_file(tfil,MK_PLOT=False,METH=5)
      tfo      = mc.single_file(tfil)
      Psolns   = tfo.get_solutions(METH=5)

      Pdict.update({reg:Psolns})
      
      # Save summary & shapefile
      mc.save_summary  (Psolns,sumname)
      mc.save_shapefile(Psolns,filename=shpname)

      if vertices is not None:
         # add total area to sumname
         loncnr,latcnr  = np.array(vertices).transpose()
         tot_area       = GS.area_polygon_ellipsoid(loncnr,latcnr)

         # append to file
         fid   = open(sumname,'a')
         fid.write('total_area_of_rectangle : '+str(tot_area))
         fid.close()
      ##########################################################

      
      if plotting:
         ##########################################################
         # Make plot
         var_opts = make_plot_options(vname,lower_limit=lower_limit)
         pobj     = fobj.plot_var(var_opts,bmap=bmap,show=False,clim=clim)[0]
         fig      = pobj.fig
         ax       = pobj.ax
         PLOTTING = True

         if vertices is not None:
            # add corners of rectangle
            bmap.plot(loncnr,latcnr,latlon=True,ax=ax,color='g',linewidth=2.5)

         for MIZi in Psolns.MIZ_info_objects:
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
         plt.show(fig)
         ax.cla()
         fig.clear()
         # finished region
         ##########################################################

   if PLOTTING:
      plt.close(fig)

   # define outputs
   return AOD_output(summary_files,\
                  tfiles,\
                  shapefiles,\
                  ['MIZ_'+var_name],\
                  regions,\
                  fobj.datetimes[time_index])
###########################################################


###########################################################
def areas_of_disagreement(fobj,time_index=0,\
      obs_type='OSISAF',obs_path=None,obs_option='multi',\
      vertices=None,regions=None,\
      do_sort=True,EastOnly=True,\
      forecast_day=None,\
      plotting=2,**kwargs):
   """
   areas_of_disagreement(fobj,time_index=0,\
      obs_type='OSISAF',obs_path=None,obs_option='multi',\
      vertices=None,regions=None,\
      do_sort=True,EastOnly=True,\
      forecast_day=None,\
      plotting=2,**kwargs)

      **kwargs to be passed onto MIZchar.get_MIZ_poly:
         outdir='.' - where to save results
         mask_corners=None - can mask out smaller region
                           - corners = [lons,lats], where lons=[DL,DR,UR,UL]
   """

   import MIZchar as mc
   PRINT_INFO  = 1

   dtmo  = fobj.datetimes[time_index]
   cdate = dtmo.strftime('%Y%m%d')

   if forecast_day is None:
      DTlabel  = 1
   else:
      dtm_start   = dtmo-timedelta(forecast_day)
      DTlabel     = dtm_start.strftime('%d %b %Y')+\
         ' +%i days' %(forecast_day)

   if obs_type == 'OSISAF':
      var_name = 'fice'
      bmap     = basemap_OSISAF()
      if obs_path is None:
      	 obs_path   = '/work/shared/nersc/msc/OSI-SAF/'+\
            dtmo.strftime('%Y')+'_nh_polstere/'
      if obs_option is None:
	 obs_option='multi'
      obsfil   = obs_path+'/ice_conc_nh_polstere-100_'+\
                  obs_option+'_'+cdate+'1200.nc'
   else:
      raise ValueError('Wrong selection variable for areas_of_disagreement')

   # observation grid & compared quantity
   if PRINT_INFO:
      print('Observation file:')
      print(obsfil)
      print('\n')

   nci         = nc_getinfo(obsfil)
   lon2,lat2   = nci.get_lonlat()
   Xobs,Yobs   = bmap(lon2,lat2)
   Zobs        = GetVar(nci,var_name,time_index=0)

   # model grid & compared quantity
   Zmod        = GetVar(fobj,var_name,time_index=time_index)
   lon,lat     = fobj.get_lonlat()
   Xmod,Ymod   = bmap(lon,lat)

   if '%' in Zobs.units:
      conv_fac = .01
   else:
      conv_fac = 1

   if 1:
      #Zref,Zint should be np.ma.array
      lon_ref,lat_ref   = lon2,lat2
      Xref,Yref,Zref    = Xobs,Yobs,conv_fac*Zobs.values # obs grid is reference;                 
      Xin,Yin,Zin       = Xmod,Ymod,Zmod.values          # to be interped from model grid onto obs grid;  Zint is np.ma.array

   # add the mask for the obs (ref) to interpolated model results (Zout)
   if PRINT_INFO:
      print('Reprojecting model...')
   Zout   = reproj_mod2obs(Xin,Yin,Zin,Xref,Yref,mask=1*Zref.mask)

   # add the mask for Zout to Zref
   Zref  = np.ma.array(Zref.data,mask=Zout.mask)


   # ==================================================================
   # calc RMSE and bias
   cdiff       = Zout-Zref
   good        = np.logical_not(Zout.mask)

   # 1st estimate (lower bound)
   # - only consider pixels with ice in both datasets
   both_ice       = np.copy(good)
   both_ice[good] = np.logical_and(Zout[good]>0.,Zref[good]>0.)
   RMSEb          = np.sqrt(np.mean(cdiff[both_ice]**2))
   BIASb          = np.mean(cdiff[both_ice])

   # 2nd estimate (upper bound)
   # - consider pixels with ice in one of the datasets
   either_ice        = np.copy(good)
   either_ice[good]  = np.logical_or(Zout[good]>0.,Zref[good]>0.)
   RMSEe             = np.sqrt(np.mean(cdiff[either_ice]**2))
   BIASe             = np.mean(cdiff[either_ice])

   if 'outdir' in kwargs:
      outdir   = kwargs['outdir']
   else:
      outdir   = '.'
   if not os.path.exists(outdir):
      os.mkdir(outdir)

   anom_fil = outdir+'/conc_anomaly_'+obs_type+'_'+cdate+'.npz'
   print('Saving '+anom_fil+'\n')
   np.savez(anom_fil,lon=lon2,lat=lat2,\
         anomaly=cdiff.data,mask=cdiff.mask)

   # write RMSE and BIAS (global variables - not regional) to summary files
   sumname  = anom_fil.replace('.npz','.txt')
   fid      = open(sumname,'w')
   fid.write('RMSE_both_ice   : '+str(RMSEb)+'\n')
   fid.write('Bias_both_ice   : '+str(BIASb)+'\n')
   fid.write('RMSE_either_ice : '+str(RMSEe)+'\n')
   fid.write('Bias_either_ice : '+str(BIASe))
   fid.close()

   if plotting>0:
      if fobj.HYCOM_region=='TP4':
         HYCreg   = 'Arctic'
      else:
         HYCreg   = fobj.HYCOM_region

      # =================================================================
      # plot model + ice edge
      cdate    = fobj.datetimes[time_index].strftime('%Y%m%dT%H%M%SZ')
      figname  = outdir+'/'+fobj.basename+'_IceEdge_'+obs_type+'_'+cdate+'.png'
      fobj.compare_ice_edge_obs(obs_type=obs_type,obs_path=obs_path,obs_option=obs_option,\
            HYCOMreg=HYCreg,figname=figname,\
            date_label=DTlabel,clabel='Sea ice concentration')

      # plot obs
      var_opts = make_plot_options('fice',ice_mask=True)
      nci.make_png(var_opts,figdir=outdir,HYCOMreg=HYCreg,\
            date_label=1,clabel='Sea ice concentration')
      # =================================================================

         
      # =================================================================
      # plot anomoly
      anom_fig = anom_fil.replace('.npz','.png')
      cmax     = .5
      clabel   = 'Concentration anomaly'

      # mask "not good" region
      cdiff = np.ma.array(cdiff.data,mask=np.logical_not(either_ice))
      Fplt.plot_anomaly(lon2,lat2,cdiff,anom_fig,\
            # text=cdate,\
            text=DTlabel,\
            HYCOM_region=HYCreg,\
            clim=[-cmax,cmax],clabel=clabel)
      # =================================================================

   # ==================================================================


   MPdict         = {'Over':{},'Under':{}}
   tfiles         = {'Over':{},'Under':{}}
   summary_files  = {'Over':{},'Under':{}}
   shapefiles     = {'Over':{},'Under':{}}

   if 0:
      # test interpolation and matching of masks
      fig   = plt.figure()
      ax1   = fig.add_subplot(1,2,1)
      ax1.imshow(Zout.transpose(),origin='upper')
      ax2   = fig.add_subplot(1,2,2)
      ax2.imshow(Zref.transpose(),origin='upper')
      plt.show(fig)
      return Xin,Yin,Zin,Xref,Yref,Zref

   if (vertices is not None) and (regions is not None):
      raise ValueError('Cannot pass in both vertices and regions')
   elif vertices is not None:
      do_sort  = False
   elif regions is not None:
      do_sort  = True

   if do_sort:

      if regions is None:
         # possible regions are:
         regions  = ['gre','bar','beau','lab','balt','les','can']

         if EastOnly:
            # concentrate on the eastern Arctic
            # (and forget Baltic Sea)
            regions.remove('balt' )
            regions.remove('les' )
            regions.remove('can' )
            regions.remove('beau')
         else:
            regions.remove('balt' )

      # for reg in ['bar']:
      for reg in regions:

         # Zout,Zref are np.ma.array objects
         if PRINT_INFO:
            print('Locating AODs for region '+reg+'...')
         Over,Under  = mc.get_AOD_polys(Zout,Zref,lon_ref,lat_ref,region=reg)
         MPdict['Over'] .update({reg:Over})
         MPdict['Under'].update({reg:Under})

         for OU in ['Over','Under']:

            fname0   = fobj.basename+'_v'+obs_type +'_'+OU+'_'+reg
            tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
            if 'all' in tfile.keys():
               tfiles[OU].update({reg:tfile['all']})

      if 0:
         # Show the 4 maps for last region (for each over/under)
         MPdict['Over'] [reg].show_maps()
         MPdict['Under'][reg].show_maps()
         return MPdict
   else:
      # not do_sort:
      regions  = ['all']
      reg      = 'all'
      if PRINT_INFO:
         print('Locating AODs for all regions...')
      Over,Under  = mc.get_AOD_polys(Zout,Zref,lon_ref,lat_ref,\
                                       vertices=vertices)
      MPdict['Over'] .update({reg:Over})
      MPdict['Under'].update({reg:Under})

      for OU in ['Over','Under']:

         fname0   = fobj.basename+'_v'+obs_type+'_'+OU
         tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
         if 'all' in tfile.keys():
            tfiles[OU].update({reg:tfile['all']})

      if 0:
         # Show the 4 maps (for over/under)
         # MPdict['Over'] ['all'].show_maps()
         MPdict['Under']['all'].show_maps()
         return MPdict

   print('\n------------------------------------------------------')
   print('Polygon info printed to files:')
   for OU in ['Over','Under']:
      for key in tfiles[OU].keys():
         print(OU+', '+key+': '+tfiles[OU][key])
      print('\n')
   print('------------------------------------------------------\n')
   # print(MPdict)

   Pdict = {'Over':{},'Under':{}}
   for OU in ['Over','Under']:
      PLOTTING = False
      for reg in tfiles[OU].keys():

         ##########################################################
         # filenames
         tfil     = tfiles[OU][reg]                      # text file with polygon outlines characterized
         figname  = tfil.replace('.txt','.png')          # plot of polygons
         shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
         sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file

         # outputs
         summary_files  [OU].update({reg:sumname})
         shapefiles     [OU].update({reg:shpname})
         ##########################################################


         ##########################################################
         # basemap for plotting
         if vertices is None:
            if do_sort:
               mapreg   = reg
            else:
               mapreg   = fobj.HYCOM_region
            bmap     = Fplt.start_HYCOM_map(mapreg)
         else:
            vlons,vlats = np.array(vertices).transpose()
            vx,vy       = GS.polar_stereographic_simple(vlons,vlats,NH=True,inverse=False)
            xy_verts    = [(vx[i],vy[i]) for i in range(len(vx))]

            # get approximate centre and limits for basemap
            vP          = shg.Polygon(xy_verts)
            width       = 5*np.sqrt(vP.area)
            height      = 5*np.sqrt(vP.area)
            xcen,ycen   = vP.centroid.coords.xy
            lonc,latc   = GS.polar_stereographic_simple(np.array([xcen]),np.array([ycen]),\
                           NH=True,inverse=True)
            # make basemap
            bmap        = Basemap(projection='stere',lon_0=lonc,lat_0=latc,\
                                    width=width,height=height,\
                                    resolution='i')
         ##########################################################


         ##########################################################
         # process each text file to get MIZ width etc
         print("MIZchar.single_file: "+tfil+"\n")
         tfo      = mc.single_file(tfil)
         Psolns   = tfo.get_solutions(METH=5)
         Pdict[OU].update({reg:Psolns})
         
         # Save summary & shapefile
         mc.save_summary  (Psolns,sumname)
         mc.save_shapefile(Psolns,filename=shpname)

         if vertices is not None:
            # add total area to sumname
            loncnr,latcnr  = np.array(vertices).transpose()
            tot_area       = GS.area_polygon_ellipsoid(loncnr,latcnr)

            # append to file
            fid   = open(sumname,'a')
            fid.write('total_area_of_rectangle : '+str(tot_area))
            fid.close()
         ##########################################################

         
         if plotting==2:
            ##########################################################
            # Make plot
            var_opts = make_plot_options('fice',ice_mask=True)
            pobj     = fobj.plot_var(var_opts,bmap=bmap,show=False,clim=[0,1])[0]
            fig      = pobj.fig
            ax       = pobj.ax
            PLOTTING = True

            if 1:
               # add observation ice edge
               bmap.contour(lon2,lat2,conv_fac*Zobs.values,[.15],lat_lon=True,\
                              colors='g',linewidths=2,ax=ax,latlon=True)

            if vertices is not None:
               # add corners of rectangle
               bmap.plot(loncnr,latcnr,latlon=True,ax=ax,color='g',linewidth=2.5)

            for MIZi in Psolns.MIZ_info_objects:
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


   # define outputs
   return AOD_output(summary_files,\
                  tfiles,\
                  shapefiles,\
                  ['Over','Under'],\
                  regions,\
                  fobj.datetimes[time_index],\
                  anomaly_file=anom_fil)
###########################################################


###########################################################
def time_average(fobj,varname,start_date=None,end_date=None,**kwargs):


      # ==================================================
      # set dates to analyse
      if start_date is not None:
         dto0  = datetime.strptime(start_date,'%Y%m%dT%H%M%SZ')
      else:
         dto0  = fobj.datetimes[0]

      if end_date is not None:
         dto1  = datetime.strptime(end_date,'%Y%m%dT%H%M%SZ')
      else:
         dto1  = fobj.datetimes[-1]
      # ==================================================


      # ==================================================
      N     = 0
      Vav   = 0
      for dt in fobj.datetimes:
         if dt>=dto0 and dt<=dto1:
            idx   = fobj.datetimes.index(dt)
            V     = fobj.get_var(varname,time_index=idx,**kwargs)
            Vav  += V.values
            N    += 1
      # ==================================================

      Vav   = (1./N)*Vav
      return Vav
###########################################################


###########################################################
def smooth(lon,lat,V,radius,mask=None,resolution=None):
    import time

    nx,ny   = lon.shape
    Vav     = np.zeros((nx,ny))

    if resolution is None:
        # determine approx resolution
        lons        = np.array([lon[0,0],lon[0,1],lon[1,1],lon[1,0]])
        lats        = np.array([lat[0,0],lat[0,1],lat[1,1],lat[1,0]])
        area        = GS.area_polygon_ellipsoid(lons,lats)
        resolution  = np.sqrt(area)

    if mask is None:
        mask    = np.logical_not(np.isnan(V))

    irng    = int(np.ceil(1.5*radius/resolution))

    start = time.time()
    print('Start smoothing (radius(km),res(km),i-j range)')
    print(radius/1.e3,resolution/1.e3,irng)

    for i in range(nx):
        for j in range(ny):
            v       = 0.
            N       = 0
            lon1    = lon[i,j]
            lat1    = lat[i,j]
            if not mask[i,j]:

                # preselect region to search over
                i0  = max(0,i-irng)
                i1  = min(nx-1,i+irng+1)
                j0  = max(0,j-irng)
                j1  = min(ny-1,j+irng+1)
                #print(i0)
                #print(i1)
                #print(j0)
                #print(j1)

                # search the area to see if distance < radius
                for i_ in range(i0,i1):
                    for j_ in range(j0,j1):
                        lon2    = lon[i_,j_]
                        lat2    = lat[i_,j_]
                        dist    = GS.greatcircledist(lat1, lon1, lat2, lon2,radians=False)
                        if dist<=radius:
                            v  += V[i_,j_]
                            N  += 1
                        # print(i,j,dist,radius,v,N)
                # sys.exit()
                # print(i,j,v,N)

                # take average
                Vav[i,j]    = v/float(N)
            else:
                Vav[i,j]    = np.nan

    end = time.time()
    print ('Smoothing time (mins):')
    print ((end - start)/60.)

    Vav = np.ma.array(Vav,mask=np.isnan(Vav))

    return Vav
###########################################################


###########################################################
class MIZmap_all:
   def __init__(self,fobj,outdir='.',start_date=None,end_date=None,step=None,dir_info=None,**kwargs):
      """
      mod_reading.MIZmap_all(fobj,outdir='.',start_date=None,end_date=None,step=1,dir_info=None,**kwargs)
         fobj a file list object or multi-record netcdf object
      """

      if not os.path.exists(outdir):
         os.mkdir(outdir)

      if step is None:
         Nrec  = fobj.number_of_time_records
         if Nrec>1:
            # time step in days
            step  = (fobj.datetimes[1]-fobj.datetimes[0]).total_seconds()/(24.*3600)
         else:
            step  = 1. #not used

      # ==================================================
      # set dates to analyse, check for missing dates
      if start_date is not None:
         dto0  = datetime.strptime(start_date,'%Y%m%dT%H%M%SZ')
         DTO   = datetime.strptime(start_date,'%Y%m%dT%H%M%SZ')
      else:
         dto0  = fobj.datetimes[0]
         DTO   = fobj.datetimes[0]

      if end_date is not None:
         dto1  = datetime.strptime(end_date,'%Y%m%dT%H%M%SZ')
      else:
         dto1  = fobj.datetimes[-1]


      self.times_to_analyse  = [dto0]
      self.missing_times     = []
      while DTO<dto1:
         DTO   = DTO+timedelta(step)
         self.times_to_analyse.append(DTO)
         if DTO not in fobj.datetimes:
            self.missing_times.append(DTO)

      Ntimes   = len(self.times_to_analyse)
      self.number_of_times_analysed = Ntimes
      self.number_of_results        = 0
      # ==================================================

     
      # ==================================================
      # loop over times:
      Init  = True
      for it,dto in enumerate(self.times_to_analyse):
         
         # restrict analysis dates
         if dto in self.missing_times:
            continue
         
         idx      = fobj.datetimes.index(dto)
         cdate    = dto.strftime('%Y%m%dT%H%M%SZ')
         outdir2  = outdir+'/'+cdate
         if not os.path.exists(outdir2):
            os.mkdir(outdir2)

         if dir_info is None:
            # run MIZ analysis
            print('Running MIZ analysis for '+cdate)
            out   = fobj.MIZmap(outdir=outdir2,time_index=idx,**kwargs)
         else:
            print('Getting summary info from '+outdir2)
            # MIZ analysis already done
            # - get summary files
            out   = summary_info(outdir2,dir_info)

         # ===================================================================
         # check if any answers present
         empty = True
         regs  = out.summary_files.keys()
         if regs is not None:
            empty    = False

            if Init:
               reg      = regs[0]
               sumfile  = out.summary_files[reg] # full path
               print('Initialising variable list from '+sumfile+'...')

               sfo      = read_MIZpoly_summary(sumfile)
               # EG OF SUMMARY FILE
               # Total_perimeter :          1146733
               # Total_area :      16769683265
               # Mean_intersecting_width :            41668
               # Mean_total_width :            41828
               # Maximum_intersecting_width :           107207
               # Maximum_total_width :           107207

               v  = vars(sfo)
               del(v['info'])
               var_names   = v.keys()

               print('\nVariables')
               for v in var_names:
                  print(v)

               print('\n')

         if empty:
            # move to next date
            continue
         else:
            self.number_of_results += 1
         # ===================================================================


         # ===================================================================
         if Init:
            Init                    = False # don't need to do this again
            self.types              = out.types
            self.regions_analysed   = out.regions_analysed

            data  = {} # data will be data[OU][reg][variable]
            for reg in out.regions_analysed:
               data.update({reg:{}})

               var_names2  = {}
               for vbl in var_names:
                  # v  = vbl+'_AOD' # change name variable is stored under
                  v  = vbl # keep name of variable
                  var_names2.update({vbl:v}) # map to the new name
                  data[reg].update({v:np.nan*np.zeros((Ntimes,))})
         # ===================================================================


         # ===================================================================
         # read all summary files to add to time series
         for reg in out.regions_analysed:

            if reg in out.summary_files:
               sumfile  = out.summary_files[reg]
               sfo      = read_MIZpoly_summary(sumfile)
               # EG OF SUMMARY FILE
               # Total_perimeter :          1146733
               # Total_area :      16769683265
               # Mean_intersecting_width :            41668
               # Mean_total_width :            41828
               # Maximum_intersecting_width :           107207
               # Maximum_total_width :           107207

               for vbl in var_names:
                  v  = var_names2[vbl]
                  data[reg][v][it]  = getattr(sfo,vbl)
      # ===================================================================


      # ===================================================================
      # convert data to time_series object
      # - gives plotting options etc
      self.time_series  = {}
      outdir3  = outdir+'/time_series'
      if not os.path.exists(outdir3):
         os.mkdir(outdir3)

      ctype = out.types[0]
      for reg in out.regions_analysed:
         ofil  = outdir3+'/time_series_'+ctype+'_'+reg+'.txt'
         self.time_series.update({reg:[]})
         self.time_series[reg]  = time_series(self.times_to_analyse,data[reg],filename=ofil)
      # ===================================================================

      return
###########################################################


###########################################################
class AODs_all:
   def __init__(self,fobj,outdir='.',start_date=None,end_date=None,step=1,dir_info=None,**kwargs):

      if not os.path.exists(outdir):
         os.mkdir(outdir)

      # ==================================================
      # set dates to analyse, check for missing dates
      if start_date is not None:
         dto0  = datetime.strptime(start_date+'T120000Z','%Y%m%dT%H%M%SZ')
         DTO   = datetime.strptime(start_date+'T120000Z','%Y%m%dT%H%M%SZ')
      else:
         dto0  = fobj.datetimes[0]
         DTO   = fobj.datetimes[0]

      if end_date is not None:
         dto1  = datetime.strptime(end_date+'T120000Z','%Y%m%dT%H%M%SZ')
      else:
         dto1  = fobj.datetimes[-1]


      self.times_to_analyse  = [dto0]
      self.missing_times     = []
      while DTO<dto1:
         DTO   = DTO+timedelta(step)
         self.times_to_analyse.append(DTO)
         if DTO not in fobj.datetimes:
            self.missing_times.append(DTO)

      Ntimes   = len(self.times_to_analyse)
      self.number_of_times_analysed = Ntimes
      self.number_of_results        = 0
      # ==================================================

     
      # ==================================================
      # loop over times:
      Init  = True
      for it,dto in enumerate(self.times_to_analyse):
         
         # restrict analysis dates
         if dto in self.missing_times:
            continue
         
         idx      = fobj.datetimes.index(dto)
         cdate    = dto.strftime('%Y%m%d')
         outdir2  = outdir+'/'+cdate
         if not os.path.exists(outdir2):
            os.mkdir(outdir2)

         if dir_info is None:
            # run AOD analysis
            print('Running AOD analysis for '+cdate)
            out   = fobj.areas_of_disagreement(outdir=outdir2,time_index=idx,**kwargs)
         else:
            print('Getting summary info from '+outdir2)
            # AOD analysis already done
            # - get summary files
            out   = summary_info(outdir2,dir_info)

         # ===================================================================
         # check if any answers present
         empty = True
         for OU in out.types:
            regs  = out.summary_files[OU].keys()
            if regs is not None:
               empty    = False

               if Init:
                  reg      = regs[0]
                  sumfile  = out.summary_files[OU][reg] # full path
                  print('Initialising variable list from '+sumfile+'...')

                  sfo      = read_MIZpoly_summary(sumfile)
                  # EG OF SUMMARY FILE
                  # Total_perimeter :          1146733
                  # Total_area :      16769683265
                  # Mean_intersecting_width :            41668
                  # Mean_total_width :            41828
                  # Maximum_intersecting_width :           107207
                  # Maximum_total_width :           107207

                  v  = vars(sfo)
                  del(v['info'])
                  var_names   = v.keys()

                  print('\nVariables')
                  for v in var_names:
                     print(v)

                  print('\n')
               break

         if empty:
            # move to next date
            continue
         else:
            self.number_of_results += 1
         # ===================================================================


         # ===================================================================
         if Init:
            Init                    = False # don't need to do this again
            self.types              = out.types
            self.regions_analysed   = out.regions_analysed

            data  = {} # data will be data[OU][reg][variable]
            for OU in out.types:
               data.update({OU:{}})

               for reg in out.regions_analysed:
                  data[OU].update({reg:{}})

                  var_names2  = {}
                  for vbl in var_names:
                     # v  = vbl+'_AOD' # change name variable is stored under
                     v  = vbl # keep name of variable
                     var_names2.update({vbl:v}) # map to the new name
                     data[OU][reg].update({v:np.nan*np.zeros((Ntimes,))})
         # ===================================================================


         # ===================================================================
         # read all summary files to add to time series
         for OU in out.types:
            for reg in out.regions_analysed:

               if reg in out.summary_files[OU]:
                  sumfile  = out.summary_files[OU][reg]
                  sfo      = read_MIZpoly_summary(sumfile)
                  # EG OF SUMMARY FILE
                  # Total_perimeter :          1146733
                  # Total_area :      16769683265
                  # Mean_intersecting_width :            41668
                  # Mean_total_width :            41828
                  # Maximum_intersecting_width :           107207
                  # Maximum_total_width :           107207

                  for vbl in var_names:
                     v  = var_names2[vbl]
                     data[OU][reg][v][it]  = getattr(sfo,vbl)
         # ===================================================================


      # ===================================================================
      # convert data to time_series object
      # - gives plotting options etc
      self.time_series  = {}
      outdir3  = outdir+'/time_series'
      if not os.path.exists(outdir3):
         os.mkdir(outdir3)

      for OU in self.types:
         self.time_series.update({OU:{}})
         for reg in self.regions_analysed:
            ofil  = outdir3+'/AOD_time_series_'+OU+'_'+reg+'.txt'
            self.time_series[OU].update({reg:[]})
            self.time_series[OU][reg]  = time_series(self.times_to_analyse,data[OU][reg],filename=ofil)
      # ===================================================================

      return
###########################################################


###########################################################
class summary_info:
   """
   mod_reading.summary_info(outdir,dir_info)
   get info about summary files in directory outdir
   sorts them according to dir_info
   - dictionary with keys 'types','regions'
   """

   def __init__(self,outdir,dir_info):

      self.types              = dir_info['types']
      self.regions_analysed   = dir_info['regions']

      flist = os.listdir(outdir)
      slist = []
      for f in flist:
         if 'summary.txt' in f:
            slist.append(f)


      self.summary_files   = {}
      for OU in self.types:
         self.summary_files.update({OU:{}})

         for reg in self.regions_analysed:
            self.summary_files[OU].update({reg:{}})

            # Check for files present that correspond to type and region
            no_file  = True
            for f in slist:
               if ('_'+OU+'_' in f) and ('_'+reg+'_' in f):
                  self.summary_files[OU][reg]   = os.path.abspath(outdir+'/'+f)
                  no_file  = False
                  break

            # delete key if no file present
            if no_file:
               del self.summary_files[OU][reg]

      # print(self.summary_files)
      return
###########################################################

def get_range_animation(fobj,var_opts,bmap,percentile_min=0,percentile_max=95):

   # boundaries of bmap
   lon0     = bmap.boundarylons.min()
   lon1     = bmap.boundarylons.max()
   lat0     = bmap.boundarylats.min()
   lat1     = bmap.boundarylats.max()
   print('\nlon/lat range for plotting')
   print(lon0,lon1,lat0,lat1)
   print('\n')

   # points inside bmap
   lon,lat  = fobj.get_lonlat()
   inside   = np.logical_and(lon>lon0,lon<lon1)
   inside   = np.logical_and(lat>lat0,inside)
   inside   = np.logical_and(lat<lat1,inside)
   del lon,lat

   N  = fobj.number_of_time_records
   vmin  = 1e30
   vmax  = -1e30
   clim  = None
   for i in range(N):
      V     = fobj.get_var(var_opts.name,time_index=i)
      good  = np.logical_not(V.values.mask)
      OK    = np.logical_and(good,inside)
      if np.any(OK):
         Vok   = V[OK]
         vmin  = min(vmin,np.percentile(V[OK],percentile_min))
         vmax  = max(vmax,np.percentile(V[OK],percentile_max))
         clim  = vmin,vmax

   return clim

###########################################################
def make_png_all(fobj,var_opts,HYCOMreg=None,figdir='.',\
      percentile_min=0,percentile_max=95,**kwargs):

   var_opts    = check_var_opts(var_opts)
   pobj        = plot_object()
   fig,ax,cbar = pobj.get()

   if HYCOMreg is None:
      HYCOMreg = fobj.HYCOM_region
      if fobj.HYCOM_region is None:
         HYCOMreg = 'TP4'

   bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
   N     = fobj.number_of_time_records

   # =====================================================
   if 'clim' not in kwargs:
      clim  = get_range_animation(fobj,var_opts,bmap,\
            percentile_min=percentile_min,percentile_max=percentile_max)
      if clim is not None:
         print('\nSetting variable range:')
         print(clim)
         print('\n')
         kwargs.update({'clim':clim})
   # =====================================================


   # =====================================================
   # loop over time records
   for i in range(N):
      pobj,bmap   = fobj.make_png(var_opts,time_index=i,pobj=pobj,\
            HYCOMreg=HYCOMreg,figdir=figdir,**kwargs)

      if i==0:
         # Fix axes position to stop it moving round
         pobj  = pobj.renew(axpos=pobj.ax.get_position())

      pobj.ax.cla()
      if pobj.cbar is not None:
         pobj.cbar.ax.clear()   # cbar.ax.clear()

      print('\n'+str(i+1)+' records done out of '+str(N))
   # =====================================================

   plt.close(pobj.fig)
   return
###########################################################


###########################################################
def make_png_pair_all(fobj,var_opts1,var_opts2,\
      HYCOMreg=None,figdir='.',\
      percentile_min=0,percentile_max=95,**kwargs):

   pobj        = plot_object()
   fig,ax,cbar = pobj.get()

   if HYCOMreg is None:
      HYCOMreg = fobj.HYCOM_region
      if fobj.HYCOM_region is None:
         HYCOMreg = 'TP4'

   bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
   N     = fobj.number_of_time_records

   # =====================================================
   if 'clim' not in kwargs:
      clim  = get_range_animation(fobj,var_opts,bmap,\
            percentile_min=percentile_min,percentile_max=percentile_max)
      if clim is not None:
         print('\nSetting variable range:')
         print(clim)
         print('\n')
         kwargs.update({'clim':clim})
   # =====================================================

   ############################################################
   # loop over time records
   for i in range(N):
      pobj,bmap   = fobj.make_png_pair(var_opts1,var_opts2,time_index=i,pobj=pobj,\
            HYCOMreg=HYCOMreg,figdir=figdir,**kwargs)

      if i==0:
         # Fix axes position to stop it moving round
         pobj  = pobj.renew(axpos=pobj.ax.get_position())

      pobj.ax.cla()
      if pobj.cbar is not None:
         pobj.cbar.ax.clear()   # cbar.ax.clear()

      print('\n'+str(i+1)+' records done out of '+str(N))
   ############################################################

   plt.close(pobj.fig)
   return
###########################################################


###########################################################
class file_list:
   def __init__(self,directory,pattern,extension,**kwargs):
      import os

      self.object_type  = 'file_list'
      self.directory    = directory
      self.extension    = extension

      lst         = os.listdir(directory)
      file_list   = []
      for fil in lst:

         fname,fext  = os.path.splitext(fil)

         # check pattern
         # print(fil)
         if (fext==extension) and (pattern in fname):
            file_list.append(fil)

      self.number_of_files          = len(file_list)
      if self.number_of_files==0:
         return

      # print(file_list)
      # print(directory)
      wsn            = '/work/shared/nersc/msc/ModelInput'
      gridpath_lut   = {'FR1':wsn+'/FramStrait_Hyc2.2.12/FR1a0.03-clean//topo',\
                        'BS1':wsn+'/BS1a0.045-clean/topo',\
                        'TP4':wsn+'/../REANALYSIS/topo'}


      HB = False
      if extension=='.a':
         self.filetype     = 'HYCOM_binary'
         self.getinfo      = HYCOM_binary_info
         HB                = True
         self.HYCOM_region = file_list[0][:3]
         if 'gridpath' not in kwargs:
            kwargs.update({'gridpath':gridpath_lut[self.HYCOM_region]})
         
      elif extension=='.nc':
         self.filetype     = 'netcdf'
         self.getinfo      = nc_getinfo
         self.HYCOM_region = 'TP4' #TODO pass in HYCOMreg?

      # get main objects
      objects     = []
      datetimes   = []
      for fil in file_list:
         obj   = self.getinfo(self.directory+'/'+fil,**kwargs)
         objects.append(obj)
         datetimes.append(obj.datetimes[0])

      # sort the time values (1st record is earliest)
      # - also reorder objects, datetimes, file_list
      lst               = sorted([(e,i) for i,e in enumerate(datetimes)])
      self.filetimes    = np.array( [e  for e,i in lst])
      self.objects      = [objects  [i] for e,i in lst]
      self.file_list    = [file_list[i] for e,i in lst]

      self.datetimes    = []
      self.file_indices = []
      for i,fobj in enumerate(self.objects):
         lst   = fobj.datetimes
         self.datetimes.extend(lst)
         for j,dto in enumerate(lst):
            self.file_indices.append((i,j))

      self.number_of_time_records   = len(self.datetimes)
      self.reference_date           = min(self.datetimes)
      self.timeunits                = 'hours'
      timevalues                    = [(dt-self.reference_date).total_seconds()/3600. for dt in datetimes]


      #set some extra variables to work with eg make_png_all
      self.basename        = self.objects[0].basename
      self.variables       = self.objects[0].variables    
      self.variables3d     = self.objects[0].variables3d  
      self.all_variables   = self.objects[0].all_variables

      self.date_strings = []
      self.time_strings = []
      for dtm in self.datetimes:
         self.date_strings.append(dtm.strftime('%Y%m%d'))
         self.time_strings.append(dtm.strftime('%H%M%S'))

      # add some methods inherited from individual objects
      self.get_lonlat   = self.objects[0].get_lonlat
      if HB:
         self.get_corners        = self.objects[0].get_corners
         self.get_depths         = self.objects[0].get_depths
         self.create_ESMF_grid   = self.objects[0].create_ESMF_grid

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
   def get_var(self,varname,time_index=0,**kwargs):
      file_index,sub_index = self.file_indices[time_index]
      return self.objects[file_index].get_var(varname,time_index=sub_index,**kwargs)
   ###########################################################


   #######################################################################
   def imshow(self,var_opts,time_index=0,**kwargs):
      """
      pobj   = self.imshow(var_opts,time_index=0,pobj=None,\
           clim=None,add_cbar=True,clabel=None,show=True,\
           test_ijs=None)
      """
      file_index,sub_index = self.file_indices[time_index]
      return imshow(self.objects[file_index],\
            var_opts,time_index=sub_index,**kwargs)
   #######################################################################


   ###########################################################
   def plot_var(self,var_opts,time_index=0,**kwargs):
      """
      pobj,bmap = self.plot_var(var_opts,time_index=0,\
         pobj=None,bmap=None,HYCOMreg='TP4',\
         clim=None,add_cbar=True,clabel=None,show=True,\
         test_lonlats=None,date_label=0):
      """
      file_index,sub_index = self.file_indices[time_index]
      return plot_var(self.objects[file_index],var_opts,\
                  time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def plot_var_pair(self,var_opts1,var_opts2,time_index=0,**kwargs):
      """
      pobj,bmap=self.plot_var_pair(var_opts1,var_opts2,pobj=None,bmap=None,**kwargs)
      """
      file_index,sub_index = self.file_indices[time_index]
      return plot_var_pair(self.objects[file_index],var_opts1,var_opts2,\
            time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def make_png(self,var_opts,time_index=0,**kwargs):
      """
      pobj,bmap=self.make_png(var_opts,pobj=None,bmap=None,figdir='.',time_index=0,date_label=2,**kwargs)
      """
      file_index,sub_index = self.file_indices[time_index]
      return make_png(self.objects[file_index],var_opts,\
                  time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def make_png_pair(self,var_opts1,var_opts2,time_index=0,**kwargs):
      """
      pobj,bmap = self.make_png_pair(var_opts1,var_opts2,\
         pobj=None,bmap=None,figdir='.',date_label=2,**kwargs)
      """
      file_index,sub_index = self.file_indices[time_index]
      return make_png_pair(self.objects[file_index],var_opts1,var_opts2,\
            time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def make_png_all(self,var_opts,**kwargs):
      """
      self.make_png_all(var_opts,HYCOMreg=None,figdir='.')
      """
      out   = make_png_all(self,var_opts,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def make_png_pair_all(self,var_opts1,var_opts2,**kwargs):
      """
      self.make_png_pair_all(var_opts1,var_opts2,HYCOMreg=None,figdir='.')
      """
      out   = make_png_pair_all(self,var_opts1,var_opts2,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def interp2points(self,varname,target_lonlats,time_index=0,mapping=None,**kwargs):
      """
      data = self.compare_ice_edge_obs(varname,target_lonlats,\
            time_index=0,mapping=None,**kwargs)
      INPUTS:
      target_lonlats = (tlon,tlat) = 2-tuple of numpy arrays
      mapping = basemap or pyproj.Proj instance
      kwargs: mask=None = mask to apply to interpolated data
      """
      return interp2points(self,varname,target_lonlats,time_index=0,mapping=None,**kwargs)
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs(self,**kwargs):
      """
      pobj,bmap,obsfil = self.compare_ice_edge_obs(pobj=None,bmap=None,time_index=0,\
         obs_type='OSISAF',date_label=1,figname=None,**kwargs)
      """
      file_index,sub_index = self.file_indices[time_index]
      return compare_ice_edge_obs(self.objects[file_index],time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def compare_ice_edge_obs_all(self,**kwargs):
      return compare_ice_edge_obs_all(self,**kwargs)
   ###########################################################


   ###########################################################
   def areas_of_disagreement(self,time_index=0,**kwargs):
      """
      self.areas_of_disagreement(time_index=0,\
         obs_type='OSISAF',obs_path=None,obs_option='multi',\
         vertices=None,regions=None,\
         do_sort=True,EastOnly=True,\
         plotting=True,**kwargs)

         **kwargs to be passed onto MIZchar.get_MIZ_poly:
            outdir='.' - where to save results
            mask_corners=None - can mask out smaller region
                              - corners = [lons,lats], where lons=[DL,DR,UR,UL]
      """
      file_index,sub_index = self.file_indices[time_index]
      return self.objects[file_index].areas_of_disagreement(time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def AODs_all(self,**kwargs):
      out   = AODs_all(self,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def MIZmap(self,time_index=0,**kwargs):
      file_index,sub_index = self.file_indices[time_index]
      return self.objects[file_index].MIZmap(time_index=sub_index,**kwargs)
   ###########################################################


   ###########################################################
   def MIZmap_all(self,**kwargs):
      out   = MIZmap_all(self,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def time_average(self,varname,**kwargs):
      return time_average(self,varname,**kwargs)
   ###########################################################

   # TODO for list of multi-record files:
   # time_average
   # AODs_all
   # MIZmap_all
   # make_png_all

######################################################################
class polygon_file_list:
   def __init__(self,directory,prefix='TP4archv_wav.',suffix='_dmax.txt',date_format='%Y%m%dT%H%M%SZ',add_day=False):
      import os

      self.object_type  = 'polygon_file_list'
      self.directory    = directory

      lst         = os.listdir(directory)
      file_list   = []
      datetimes   = []
      for fil in lst:

         if not ((suffix in fil) and (prefix in fil)):
            continue

         cdate = fil.strip(suffix).strip(prefix)
         print(cdate)
         dto   = datetime.strptime(cdate,date_format)


         if add_day:
            # model has julian day starting at 0
            # - may sts need to add 1 day
            dto  += timedelta(1)

         file_list.append(fil)
         datetimes.append(dto)

      self.number_of_time_records   = len(file_list)
      if self.number_of_time_records==0:
         print('file list empty')
         return

      self.reference_date  = min(datetimes)
      self.timeunits       = 'hours'
      timevalues           = [(dt-self.reference_date).total_seconds()/3600. for dt in datetimes]

      # sort the time values (1st record is earliest)
      # - also reorder objects, datetimes, file_list
      TV                   = sorted([(e,i) for i,e in enumerate(timevalues)])
      self.timevalues,ii   = np.array(TV).transpose()
      self.datetimes       = [datetimes[int(i)] for i in ii]
      self.file_list       = [self.directory+'/'+file_list[int(i)] for i in ii]

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
   def read_file(self, time_index=0,vertices=None):
      """
      pil = self.read_file(time_index=0)
      pil is a poly_info_list object
      """
      import MIZchar as mc
      tfil  = self.file_list[time_index]
      pil   = mc.single_file(tfil)

      if vertices is None:
         return pil
      else:
         return pil.reduce_area(vertices)
   ###########################################################


   ###########################################################
   def get_solutions(self, time_index=0,vertices=None,**kwargs):
      """
      mil = self.get_solutions(time_index=0,vertices=None,**kwargs)
      mil is a MIZ_info_list object
      **kwargs are for MIZchar.poly_info_list.get_solutions
      """
      import MIZchar as mc
      tfil  = self.file_list[time_index]
      pil   = mc.single_file(tfil)

      if vertices is not None:
         pil   = pil.reduce_area(vertices)
      return pil.get_solutions(**kwargs)
   ###########################################################
