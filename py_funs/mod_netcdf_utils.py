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
   vbl   = MR.var_object(vals,extra_atts=[attlist,attvals])

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

      self.basename     = ss[-1].strip('.nc')
      self.filetype     = 'netcdf'

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
      self.time_dim     = ('time' in self.dimensions)

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

      self.variable_list   = vkeys
      self.variables       = vkeys
      # self.variables3d     = None  #TODO enable treatment of 3d fields
      self.all_variables   = vkeys
      ########################################################

      nc.close()
      return
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
   def areas_of_disagreement(self,obs_type='OSISAF',time_index=0,do_sort=True,EastOnly=True,\
         plotting=True,HYCOMreg='Arctic',**kwargs):
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

      if '%' in Zobs.units:
         conv_fac = .01
      else:
         conv_fac = 1

      if 1:
         #Zref,Zint should be np.ma.array
         lon_ref,lat_ref   = lon2,lat2
         Xref,Yref,Zref    = Xobs,Yobs,conv_fac*Zobs.values # obs grid is reference;                 
         Xint,Yint,Zint    = Xmod,Ymod,Zmod.values          # to be interped from model grid onto obs grid;  Zint is np.ma.array

      # *interpolate
      # *add the mask for the ref to Arr
      Arr   = reproj_mod2obs(Xint,Yint,Zint,Xref,Yref,mask=1*Zref.mask)

      # add the mask for Arr to Zref
      Zref  = np.ma.array(Zref.data,mask=Arr.mask)

      if 0:
         # test interpolation and matching of masks
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,2,1)
         im1   = ax1.imshow(Arr)
         fig.colorbar(im1)
         #
         ax2   = fig.add_subplot(1,2,2)
         im2   = ax2.imshow(Zref)
         fig.colorbar(im2)
         fig.show()
         return Xint,Yint,Zint,Xref,Yref,Zref

      MPdict   = {'Over':{},'Under':{}}
      tfiles   = {'Over':{},'Under':{}}

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
         Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref)
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
               HYCOMreg   = reg
            ##########################################################


            ##########################################################
            # process each text file to get MIZ width etc
            print("MIZchar.single_file: "+tfil+"\n")
            bmap     = Fplt.start_HYCOM_map(HYCOMreg)
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
               fig.clf()
               # finished region
               ##########################################################

         if PLOTTING:
            plt.close(fig)

      return MPdict,tfiles,Pdict
   ###########################################################

###########################################################
