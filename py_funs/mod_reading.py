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


###############################################
class time_series:

   ############################################
   def __init__(self,dates,data,units=None,filename=None):
      self.dates     = dates
      self.data      = data
      self.units     = units
      self.filename  = filename
      self.variables = data.keys()
      return
   ############################################

   ############################################
   def plot(self,var_name,refdate=None,timeunits='days',yscaling=1.,pobj=None,**kwargs):
      if pobj is None:
         pobj  = plot_object()

      if timeunits=='days':
         xfac  = 24*3600. # seconds in 1 day
      elif timeunits=='hours':
         xfac  = 3600. # seconds in 1h
      elif timeunits=='minutes':
         xfac  = 60. # seconds in 1min
      else:
         xfac  = 1. # seconds in 1min
      if refdate is None:
         refdate  = self.dates[0]

      x  = np.array([(dt-refdate).total_seconds()/xfac for dt in self.dates])
      y  = self.data[var_name]*yscaling

      lin  ,= pobj.ax.plot(x,y,**kwargs)
      return pobj,lin
   ############################################
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
      ss = lin.split()
      dates.append(datetime.strptime(ss[0],'%Y%m%dT%H%M%SZ'))
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
      

      self.filename     = tfil
      self.datetime     = None
      self.time_in_days = None
      self.datetime_ref = datetime(1901,1,1)

      if (cdate is not None) and (ctime is not None):
         self.datetime  = datetime.strptime(cdate+' '+ctime,'%Y%m%d %H%M%S')
      elif (cdate is not None):
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
   Z1d            = 1*Z1.data
   Z1d[Z1.mask]   = np.nan

   X1vec = X1.reshape(X1.size)
   Y1vec = Y1.reshape(Y1.size)
   Z1vec = Z1.reshape(Z1.size)
   Z1vec = Z1d.reshape(Z1.size)
   C = [X1vec,Y1vec]
   C = np.array(C)
   C = C.T # input coords from X1,Y1; Z1 is array to interp; X2,Y2 are output matrices

   # Interpolation can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
   Z2    = grd(C,Z1vec,(X2,Y2),method='linear')
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
                  'concentration','sea_ice_concentration'])

   # ice thick alt names
   lists.append(['hicem','hice','ice_thick','icetk',\
                  'sea_ice_thickness','thickness','sea_ice_concentration'])

   # floe size alt names
   lists.append(['dfloe','dmax','Dfloe','Dmax'])

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


#######################################################################
def plot_var(fobj,var_opts,time_index=0,\
      pobj=None,bmap=None,HYCOMreg=None,\
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

   if 'show' in kwargs:
      show           = kwargs['show']
      kwargs['show'] = False
      pobj,bmap      = fobj.plot_var(var_opts,pobj=pobj,bmap=bmap,**kwargs)
   else:
      show        = False
      pobj,bmap   = fobj.plot_var(var_opts,pobj=pobj,bmap=bmap,show=False,**kwargs)

   dtmo     = fobj.datetimes[time_index]
   datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')

   if fobj.HYCOM_region=='TP4':
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

   dtmo     = fobj.datetimes[time_index]
   datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
   if fobj.HYCOM_region=='TP4':
      xyann = (0.05,.925)
   else:
      xyann = (0.4,.925)

   if date_label==1:
      tlabel   = dtmo.strftime('%d %b %Y')
      pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
   elif date_label==2:
      tlabel   = dtmo.strftime('%d %b %Y %H:%M')
      pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)


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
def MIZmap(fobj,var_name='dmax',time_index=0,vertices=None,\
      do_sort=False,EastOnly=True,plotting=True,**kwargs):
   """
   Call  : fobj.MIZmap(var_name='dmax',vertices=None,\
               do_sort=False,EastOnly=True,plotting=True,**kwargs)
   Inputs:
      var_name is variable to find MIZ from
      **kwargs to be passed onto MIZchar.get_MIZ_poly:
         outdir='.' - where to save results
         mask_corners=None - can mask out smaller region
                           - corners = [lons,lats], where lons=[DL,DR,UR,UL]
   Returns: MIZchar.MIZpoly object
   """

   import MIZchar as mc

   vname = check_names(var_name,fobj.variables)
   if var_name == 'dmax':
      # FSD MIZarray(1-
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,300]# for plotting
      lower_limit = .1     # for plotting
   elif var_name == 'fice':
      # conc MIZ
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,1]  # for plotting
      lower_limit = .15    # for plotting
   elif var_name == 'hice':
      # thin ice areas
      Arr         = GetVar(fobj,vname,time_index=time_index)
      clim        = [0,2.] # for plotting
      lower_limit = .01    # for plotting
   else:
      raise ValueError('Wrong selection variable for MIZmap')

   print("MIZchar.get_MIZ_poly\n")
   lon,lat  = fobj.get_lonlat()
   MPdict   = {}
   tfiles   = {}

   if vertices is not None:
      do_sort  = False
   
   if do_sort:
      # possible regions are:
      regions  = ['gre','bar','beau','lab','balt','les','can']

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

      for reg in regions:
         mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,region=reg)
         MPdict.update({reg:mp})

         fname0   = fobj.basename+'_'+var_name +'_'+reg
         tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
         if 'all' in tfile.keys():
            tfiles.update({reg:tfile['all']})

      if 0:
         MPdict['gre'].show_maps()
         return MPdict

   else:
      reg   = 'all'
      mp    = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,vertices=vertices)
      MPdict.update({reg:mp})
      #
      fname0   = fobj.basename+'_'+var_name
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


      if vertices is None:
         ##########################################################
         if do_sort:
            mapreg   = reg
         else:
            mapreg   = fobj.HYCOM_region
         ##########################################################


         ##########################################################
         # process each text file to get MIZ width etc
         print("MIZchar.single_file: "+tfil+"\n")
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
            bmap.plot(loncnr,latcnr,latlon=True,ax=ax,color='g',linewidth=2.5)

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
def areas_of_disagreement(fobj,time_index=0,\
      obs_type='OSISAF',obs_path=None,obs_option='multi',\
      do_sort=True,EastOnly=True,plotting=True,**kwargs):
   # kwargs: outdir='.',do_sort=True

   import MIZchar as mc
   PRINT_INFO  = 1

   dtmo  = fobj.datetimes[time_index]
   if obs_type == 'OSISAF':
      var_name = 'fice'
      bmap     = basemap_OSISAF()
      if obs_path is None:
      	 obs_path   = '/work/shared/nersc/msc/OSI-SAF/'+\
            dtmo.strftime('%Y')+'_nh_polstere/'
      if obs_option is None:
	 obs_option='multi'
      obsfil   = obs_path+\
            '/ice_conc_nh_polstere-100_'+obs_option+'_'+\
            dtmo.strftime('%Y%m%d')+'1200.nc'
   else:
      raise ValueError('Wrong selection variable for areas_of_disagreement')

   # observation grid & compared quantity
   if PRINT_INFO:
      print(obsfil)

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
      Xint,Yint,Zint    = Xmod,Ymod,Zmod.values          # to be interped from model grid onto obs grid;  Zint is np.ma.array

   # add the mask for the ref to output (Arr)
   if PRINT_INFO:
      print('Reprojecting model...')
   Arr   = reproj_mod2obs(Xint,Yint,Zint,Xref,Yref,mask=1*Zref.mask)

   # add the mask for Arr to Zref
   Zref  = np.ma.array(Zref.data,mask=Arr.mask)

   MPdict   = {'Over':{},'Under':{}}
   tfiles   = {'Over':{},'Under':{}}

   if 0:
      # test interpolation and matching of masks
      fig   = plt.figure()
      ax1   = fig.add_subplot(1,2,1)
      ax1.imshow(Arr.transpose(),origin='upper')
      ax2   = fig.add_subplot(1,2,2)
      ax2.imshow(Zref.transpose(),origin='upper')
      plt.show(fig)
      return Xint,Yint,Zint,Xref,Yref,Zref

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
      else:
         regions.remove('balt' )

      # for reg in ['bar']:
      for reg in regions:

         # Arr,Zref are np.ma.array objects
         if PRINT_INFO:
            print('Locating AODs for region '+reg+'...')
         Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref,region=reg)
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
      reg         = 'all'
      if PRINT_INFO:
         print('Locating AODs for all regions...')
      Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref)
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
         ##########################################################


         ##########################################################
         if do_sort:
            mapreg   = reg
         else:
            mapreg   = fobj.HYCOM_region
         ##########################################################


         ##########################################################
         # process each text file to get MIZ width etc
         print("MIZchar.single_file: "+tfil+"\n")
         bmap     = Fplt.start_HYCOM_map(mapreg)
         tfo      = mc.single_file(tfil)
         Psolns   = tfo.get_solutions(METH=5)
         Pdict[OU].update({reg:Psolns})
         
         # Save summary & shapefile
         mc.save_summary  (Psolns,sumname)
         mc.save_shapefile(Psolns,filename=shpname)
         ##########################################################

         
         if plotting:
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
def make_png_all(fobj,var_opts,HYCOMreg=None,figdir='.',**kwargs):

   pobj        = plot_object()
   fig,ax,cbar = pobj.get()

   if HYCOMreg is None:
      HYCOMreg = fobj.HYCOM_region
      if fobj.HYCOM_region is None:
         HYCOMreg = 'TP4'

   bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
   N     = fobj.number_of_time_records

   ############################################################
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
   ############################################################

   plt.close(pobj.fig)
   return
###########################################################


###########################################################
def make_png_pair_all(fobj,var_opts1,var_opts2,\
      HYCOMreg=None,figdir='.',**kwargs):

   pobj        = plot_object()
   fig,ax,cbar = pobj.get()

   if HYCOMreg is None:
      HYCOMreg = fobj.HYCOM_region
      if fobj.HYCOM_region is None:
         HYCOMreg = 'TP4'

   bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
   N     = fobj.number_of_time_records

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


class file_list:
   def __init__(self,directory,pattern,extension,**kwargs):
      import os

      self.directory = directory
      self.extension = extension

      lst         = os.listdir(directory)
      file_list   = []
      for fil in lst:

         fname,fext  = os.path.splitext(fil)

         # check pattern
         if (fext==extension) and (pattern in fname):
            file_list.append(fil)

      # print(file_list)
      # print(directory)
      wsn            = '/work/shared/nersc/msc/ModelInput'
      gridpath_lut   = {'FR1':wsn+'/FramStrait_Hyc2.2.12/FR1a0.03-clean//topo',\
                        'BS1':wsn+'/BS1a0.045-clean/topo',\
                        'TP4':wsn+'/../REANALYSIS/topo'}
      HB = False
      if extension=='.a':
         self.getinfo      = HYCOM_binary_info
         HB                = True
         self.HYCOM_region = file_list[0][:3]
         if 'gridpath' not in kwargs:
            kwargs.update({'gridpath':gridpath_lut[self.HYCOM_region]})
         
      elif extension=='.nc':
         self.getinfo      = nc_getinfo
         self.HYCOM_region = 'TP4' #TODO pass in HYCOMreg?

      # get main objects
      objects     = []
      datetimes   = []
      for fil in file_list:
         obj   = self.getinfo(self.directory+'/'+fil,**kwargs)
         objects.append(obj)
         datetimes.append(obj.datetimes[0])

      self.reference_date  = min(datetimes)
      self.timeunits       = 'hours'
      timevalues           = [(dt-self.reference_date).total_seconds()/3600. for dt in datetimes]

      # sort the time values (1st record is earliest)
      # - also reorder objects, datetimes, file_list
      TV                   = sorted([(e,i) for i,e in enumerate(timevalues)])
      self.timevalues,ii   = np.array(TV).transpose()
      self.objects         = [objects  [int(i)] for i in ii]
      self.datetimes       = [datetimes[int(i)] for i in ii]
      self.file_list       = [file_list[int(i)] for i in ii]

      #set some extra variables to work with eg make_png_all
      self.number_of_time_records   = len(ii)
      self.basename                 = self.objects[0].basename
      self.variables                = self.objects[0].variables    
      self.variables3d              = self.objects[0].variables3d  
      self.all_variables            = self.objects[0].all_variables

      self.date_strings = []
      self.time_strings = []
      for dtm in self.datetimes:
         self.date_strings.append(dtm.strftime('%Y%m%d'))
         self.time_strings.append(dtm.strftime('%H%M%S'))

      # add some methods inherited from individual objects
      self.get_lonlat   = self.objects[0].get_lonlat
      if HB:
         self.get_depths   = self.objects[0].get_depths

      return
   ###########################################################


   ###########################################################
   def plot_var(self,var_opts,time_index=0,**kwargs):
      out   = plot_var(self.objects[time_index],\
                  var_opts,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def plot_var_pair(self,var_opts1,var_opts2,time_index=0,**kwargs):
      out   = plot_var_pair(self.objects[time_index],\
                  var_opts1,var_opts2,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def make_png(self,var_opts,time_index=0,**kwargs):
      out   = make_png(self.objects[time_index],\
                  var_opts,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def make_png_pair(self,var_opts1,var_opts2,time_index=0,**kwargs):
      out   = make_png_pair(self.objects[time_index],\
                  var_opts1,var_opts2,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def make_png_all(self,var_opts,**kwargs):
      out   = make_png_all(self,var_opts,**kwargs)
      return out
   ###########################################################


   ###########################################################
   def make_png_pair_all(self,var_opts1,var_opts2,**kwargs):
      out   = make_png_pair_all(self,var_opts1,var_opts2,**kwargs)
      return out
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

      if not os.path.exists(figdir):
         os.mkdir(figdir)


      N  = len(self.objects)
      for i,obj in enumerate(self.objects):

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

######################################################################
# compare_ice_edge_obs_all
