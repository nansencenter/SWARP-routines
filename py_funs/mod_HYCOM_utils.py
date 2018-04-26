import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime,timedelta
import fns_plotting as Fplt
import os,sys
import mod_reading as MR

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
def get_array_from_HYCOM_binary(infile,recno,dims=None,grid_dir='.',mask_land=True):
   """
   CALL: A=get_array_from_HYCOM_binary(infile,recno,dims=None,grid_dir='.',mask_land=True):
   routine to get the array from HYCOM binary files (*.a or *ICE.uf)
   * infile is path to file
   * recno=1 is 1st record 
   * dims=[nx,ny] - shape of stored arrays
   * grid_dir: location of regional.grid.b and regional.depth.a file
     > regional.grid.b is used to get nx,ny if dims is not given
     > regional.depth.a is used to get land mask if infile is a .uf file
       & mask_land is True
   """

   import struct

   ######################################################################
   # get size of grid
   if dims is None:
      if 'regional.grid' in infile:
         # infile is a grid file (regional.grid.a)
         # - check .b file for size of grid
         bfile = infile[:-2]+'.b'
      else:
         # check regional.grid.b file for size of grid
         bfile = grid_dir+'/regional.grid.b'
         if not os.path.exists(bfile):
            raise ValueError('Grid file not present: '+bfile)

      bid   = open(bfile,'r')
      nx    = int( bid.readline().split()[0] )
      ny    = int( bid.readline().split()[0] )
      bid.close()
   else:
      nx = dims[0]
      ny = dims[1]
   ######################################################################


   ######################################################################
   nrec        = nx*ny
   fname,fext  = os.path.splitext(infile)
   if fext=='.a':
      # usual HYCOM binary file
      fmt_size = 4      # files are single precision
      n0       = 4096   # stores records in multiples of 4096
   elif fext=='.uf':
      # HYCOM ICE restart binary file
      fmt_size = 8      # files are double precision
      n0       = nrec   # stores records in multiples of nrec (no padding like in .a files)
   else:
      raise ValueError('Unknown file extension '+fext)
   ######################################################################


   ######################################################################
   # set record size, skip to record number
   if fmt_size==4:
      fmt_py   = 'f' # python string for single
   else:
      fmt_py   = 'd' # python string for double

   if nrec%n0==0:
      Nhyc  = nrec
   else:
      Nhyc  = (1+nrec/n0)*n0
   rec_size = Nhyc*fmt_size
   #
   fid   = open(infile,'rb')
   for n in range(1,recno):
      fid.seek(rec_size,1) # seek in bytes (1: reference is current position)
   ######################################################################


   ######################################################################
   # read data and close file
   data  = fid.read(rec_size)
   fid.close()
   ######################################################################


   ######################################################################
   # rearrange into correctly sized array
   fld   = struct.unpack('>'+Nhyc*fmt_py,data)  # NB BIG-ENDIAN so need '>'
   fld   = np.array(fld[0:nx*ny])               # select the 1st nx,ny - rest of the Nhyc record is rubbish
   fld   = fld.reshape((ny,nx)).transpose()     # need to transpose because of differences between
   ######################################################################
                                                # python/c and fortran/matlab 

   ######################################################################
   # handle land cells:
   land_thresh = 1.e30# on land 1.2677e30 
   if fext=='.a':
      if mask_land:
         fld[fld>land_thresh] = np.nan
      else:
         fld[fld>land_thresh] = 0.
   elif mask_land:
      # need to read in depths file
      # (fields are 0 on land)
      depfil   = grid_dir+'/regional.depth.a'
      deps     = get_array_from_HYCOM_binary(depfil,1,dims=[nx,ny])

      # depth is set to NaN on land
      fld[np.isnan(deps)]  = np.nan
   ######################################################################

   return fld
##############################################################


##############################################################
def get_record_numbers_HYCOM(bfile):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)


   bid   = open(bfile,'r')
   if 'restart' not in bfile:
      word  = bid.readline().split()[0] # 1st word in line 
      while word!='field':
         word  = bid.readline().split()[0] # 1st word in line 
   else:
      # restart file: just read 1st 2 lines
      for i in range(2):
         bid.readline()

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
      xmin  = float(lin.split()[-2]) # min where defined
      xmax  = float(lin.split()[-1]) # max where defined

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
class HYCOM_grid_info:

   ###################################################################
   def __init__(self,gridpath='.'):
      """
      call hgi = mod_reading.HYCOM_grid_info(gridpath='.')
      methods:
      hgi.get_corners    (q-points)
      hgi.get_centres    (p-points)    - option to only have p-points that are surrounded by q-points
      hgi.get_grid_sizes (scux,scvy)   - option to only have p-points that are surrounded by q-points
      hgi.get_areas      (scux*scvy)   - option to only have p-points that are surrounded by q-points
      """

      self.gridpath  = gridpath
      flist          = [gridpath+'/regional.grid.a',
                        gridpath+'/regional.grid.b']

      # check if grid files are present
      for fil in flist:
         if not os.path.exists(fil):
            raise ValueError(fil +' not present - change "gridpath" argument')
      self.afile,self.bfile   = flist

      self.afile_depth  = gridpath+'/regional.depth.a'
      self.bfile_depth  = gridpath+'/regional.depth.b'

      # get grid size & record numbers
      self._read_bfile()

      return
   ###################################################################


   ###################################################################
   def _read_bfile(self):
      bid   = open(self.bfile,'r')

      # # ==============================================
      # # get grid size
      # line  = bid.readline()
      # while 'idm' not in line:
      #    line  = bid.readline()
      # nx    = int( line.split()[0] )

      # line  = bid.readline()
      # while 'jdm' not in line:
      #    line  = bid.readline()
      # ny    = int( line.split()[0] )
      # bid.close()
      # # ==============================================

      lines = bid.readlines()
      bid.close()

      # get grid size
      self.Nx  = int(lines[0].split()[0])
      self.Ny  = int(lines[1].split()[0])

      R  = {}
      for i,lin in enumerate(lines[3:]):
         word  = lin.split()[0].strip(':') # 1st word in line 
         R.update({word:i+1})

      self.record_numbers  = R
      return
   ###################################################################


   ###################################################################
   def get_array(self,vbl):
      recno = self.record_numbers[vbl]
      return get_array_from_HYCOM_binary(self.afile,recno)
   ###################################################################


   ###################################################################
   def get_corners(self):
      """
      call self.get_corners()
      returns qlon,qlat - corners of grid cells if only using inner points
      """
      qlon  = self.get_array('qlon')
      qlat  = self.get_array('qlat')
      return qlon,qlat
   ###################################################################


   ###################################################################
   def get_centres(self,inner_points=False):
      """
      call self.get_centres(inner_points=False)

      grid is arranged (Arakawa C-grid)
      i=-0.5   : q-v-q-v-q-...v-q-v-q-v
      i=0      : u-p-u-p-u-...p-u-p-u-p
      i=0.5    : q-v-q-v-q-...v-q-v-q-v
      i=2      : u-p-u-p-u-...p-u-p-u-p
      ...
      i=nx-1.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx-1   : u-p-u-p-u-...p-u-p-u-p
      i=nx-0.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx     : u-p-u-p-u-...p-u-p-u-p

      if inner_points:
         return plon[:-1,:-1],plat[:-1,:-1]  # ie p-points not surrounded by q points
         # NB 1st/last rows/cols are "land" (no valid data)
      else:
         return plon,plat                    # ie all p-points
      """

      plon  = self.get_array('plon')
      plat  = self.get_array('plat')
      if inner_points:
         return plon[:-1,:-1],plat[:-1,:-1]
      else:
         return plon,plat
   ###################################################################


   ###################################################################
   def get_grid_sizes(self,loc='p',inner_points=False):
      """
      call self.get_centres(inner_points=False)

      grid is arranged (Arakawa C-grid)
      i=-0.5   : q-v-q-v-q-...v-q-v-q-v
      i=0      : u-p-u-p-u-...p-u-p-u-p
      i=0.5    : q-v-q-v-q-...v-q-v-q-v
      i=2      : u-p-u-p-u-...p-u-p-u-p
      ...
      i=nx-1.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx-1   : u-p-u-p-u-...p-u-p-u-p
      i=nx-0.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx     : u-p-u-p-u-...p-u-p-u-p

      if inner_points:
         return scux[:-1,:-1],scvy[:-1,:-1]  # ie grid sizes corresponding to p-points not surrounded by q points
         # NB 1st/last rows/cols are "land" (no valid data)
      else:
         return scux,scvy                    # ie grid sizes corresponding to all p-points
      """
      if loc not in ['p','q','u','v']:
         raise ValueError('Unkwown value for "loc": '+loc)

      dx = self.get_array('sc'+loc+'x')
      dy = self.get_array('sc'+loc+'y')
      if inner_points:
         return dx[:-1,:-1],dy[:-1,:-1]
      else:
         return dx,dy
   ###################################################################

   ###################################################################
   def get_areas(self,**kwargs):
      """
      call self.get_areas(inner_points=False)

      grid is arranged (Arakawa C-grid)
      i=-0.5   : q-v-q-v-q-...v-q-v-q-v
      i=0      : u-p-u-p-u-...p-u-p-u-p
      i=0.5    : q-v-q-v-q-...v-q-v-q-v
      i=2      : u-p-u-p-u-...p-u-p-u-p
      ...
      i=nx-1.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx-1   : u-p-u-p-u-...p-u-p-u-p
      i=nx-0.5 : q-v-q-v-q-...v-q-v-q-v
      i=nx     : u-p-u-p-u-...p-u-p-u-p

      if inner_points:
         return scux[:-1,:-1],scvy[:-1,:-1]  # ie grid sizes corresponding to p-points not surrounded by q points
         # NB 1st/last rows/cols are "land" (no valid data)
      else:
         return scux,scvy                    # ie grid sizes corresponding to all p-points
      """
      dx,dy   = self.get_grid_sizes(loc='p',**kwargs)
      return dx*dy
   ###################################################################


   ###################################################################
   def get_depths(self,inner_points=False):
      """
      self.get_depths(inner_points=False)
      *can remove outer points if inner_points=True
      *returns masked array
      """
      dep = get_array_from_HYCOM_binary(self.afile_depth,1,grid_dir=self.gridpath)
      if inner_points:
         dep   =  dep[:-1,:-1]

      mask  = np.isnan(dep)
      return np.ma.array(dep,mask=mask)
   ###################################################################


   ###################################################################
   def land_mask(self,**kwargs):
      """
      self.land_mask(**kwargs)
      **kwargs for self.get_depth: inner_points=True/False
      returns:
      Boolean matrix size of full/"inner grid"
      """
      dep   = self.get_depths(**kwargs)
      return (dep==0.) # mask land out
   ###################################################################


   #########################################################################
   def create_ESMF_grid(self,do_mask=True):
      """
      create_ESMF_grid(centres,corners,mask=None,\
         coordTypeKind='r8',coordSys=0,numPeriDims=0
      inputs:
      *coordSys=0:lon,lat (degrees); 1:lon,lat (radians); 2:x,y,z
      *coordTypeKind eg 'r8':real double
      *numPeriDims = number of periodic dimensions
      """
      import ESMF_utils as EU

      corners  = self.get_corners()                  #(qlon,qlat)
      centres  = self.get_centres(inner_points=True) #(plon,plat)
      areas    = self.get_areas  (inner_points=True) # scvx*scuy

      # ========================================================
      mask  = None
      if do_mask:
        # set up the grid mask
        mask = self.land_mask(inner_points=True)
      # ========================================================
      
      Egrid = EU.create_ESMF_grid(centres,corners,AREA=areas,MASK=mask)
      return Egrid
      ###################################################################

   #########################################################################
   def plot_depth(self,figname=None,show=False,imshow=False):
      """
      plot the depth
      """
      from matplotlib import pyplot as plt

      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      depth = self.get_depths()

      mindep      = depth.min()
      maxdep      = depth.max()
      print("\nMin depth (m) = "+str(mindep))
      print("Max depth (m) = "+str(maxdep)+"\n")

      if imshow:
         IM = ax.imshow(depth.transpose(),origin='lower')
      else:
         # need lon,lat
         qlon,qlat   = self.get_corners()

         # make basemap for plotting
         import fns_plotting as Fplt
         bbox  = [qlon.min(),qlat.min(),qlon.max(),qlat.max()]
         bmap  = Fplt.start_map(bbox)

         IM    = bmap.pcolor(qlon,qlat,depth,latlon=True,ax=ax,vmin=mindep,vmax=maxdep)
         Fplt.finish_map(bmap,do_fill=False)

      # colorbar
      cbar  = fig.colorbar(IM)
      cbar.set_label("Depth, m",rotation=270,labelpad=20,fontsize=16)

      if figname is not None:
         print("\nSaving "+figname+"\n")
         fig.savefig(figname)
      else:
         show  = True

      if show:
         plt.show(fig)

      ax.cla()
      plt.close(fig)
      return
      ###################################################################

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

      basename             = fname.split('/')[-1]
      self.basename        = basename[:-2]
      self.HYCOM_region    = basename[:3]
      self.filetype        = 'HYCOM_binary'


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

      self.HYCOM_grid_info = HYCOM_grid_info(self.gridpath)
      nx,ny = self.HYCOM_grid_info.Nx,self.HYCOM_grid_info.Ny
      # # get grid size
      # bid   = open(self.gridpath+'/regional.grid.b','r')
      # line  = bid.readline()
      # while 'idm' not in line:
      #    line  = bid.readline()
      # nx    = int( line.split()[0] )

      # line  = bid.readline()
      # while 'jdm' not in line:
      #    line  = bid.readline()
      # ny    = int( line.split()[0] )
      # bid.close()

      self.dims   = [nx,ny]
      self.Nx     = nx
      self.Ny     = ny
      #######################################################################


      #######################################################################
      # date
      bid   = open(self.bfile,'r')
      line  = bid.readline()
      EOF   = (line=='')
      if 'restart' not in self.bfile:
         while 'model day' not in line and not EOF:
            line  = bid.readline()
            EOF   = (line=='')

         line              = bid.readline()
         self.time_value   = float(line.split()[-5]) # model time (days)
         bid.close()
      else:
         line              = bid.readline()
         self.time_value   = float(line.split()[-2])
         bid.close()


      self.reference_date  = datetime(1900,12,31)
      if (self.time_value==0) and ('archv_wav' in self.basename):
         # bug in output of archv_wav at first time step
         # - determine time from basename
         idx      = self.basename.find('.')+1
         yr       = int(self.basename[idx:idx+4])
         jday     = int(self.basename[idx+5:idx+8])
         refdate  = datetime(yr,1,1)+timedelta(jday)
         cdate    = refdate.strftime('%Y%m%dT')+self.basename[idx+9:idx+15]+'Z'
         #
         self.datetime     = datetime.strptime(cdate,'%Y%m%dT%H%M%SZ')
         self.time_value   = (self.datetime-self.reference_date).total_seconds()/(24.*3600.)

      if 'DAILY' in self.basename:
         # file written at end of day - change to midday
         self.time_value   = self.time_value -.5 # model time (days)

      self.datetime        = self.reference_date+timedelta(self.time_value)
      self.datetimes       = [self.datetime]
      self.time_values     = [self.time_value]
      self.time_units      = 'days'
      #######################################################################

      gi = self.HYCOM_grid_info
      self.create_ESMF_grid   = gi.create_ESMF_grid
      self.get_lonlat         = gi.get_centres
      self.get_corners        = gi.get_corners
      self.get_depths         = gi.get_depths
      # self.get_fixed_lonlat   = gi.get_corners

      return # __init__
   #######################################################################


   #######################################################################
   def get_fixed_lonlat(self,bmap):

      if 1:
         # get qlon,qlat
         lon,lat  = self.HYCOM_grid_info.get_corners()

      elif 1:
         print('hay')
         gfil     = self.gridpath+'/regional.grid.a'
         # get ulon,vlat
         lon   = get_array_from_HYCOM_binary(gfil,5,\
                     dims=self.dims,grid_dir=self.gridpath)
         lat   = get_array_from_HYCOM_binary(gfil,8,\
                     dims=self.dims,grid_dir=self.gridpath)

      else:
         gfil     = self.gridpath+'/regional.grid.a'

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
   def get_var(self,vname,time_index=None,inner_points=False):
      """
      vbl=get_var(vname,time_index=None) - vname is a string of 2d variable name (surface variable name,layer=0)
      vbl=get_var([vname,layer],time_index=None) - vname is a string of 2d or 3d variable name (layer=0 is surface, layer=1 to kdm are ocean layers)
      *time_index is not used - only a place-holder for some routines which handle netcdf as well
      *vbl is a mod_reading.var_opt class (vbl.values is type np.ma.array - eg land is masked)
      """

      if type(vname)==type("hi"):
         # 2d var
         layer = 0
         if vname!="Volume":
            vname = MR.check_names(vname,self.variables)
            recno = self.record_numbers[vname]
            xmin  = self.minvals2d[vname]
            xmax  = self.maxvals2d[vname]
            #
            vbl   = get_array_from_HYCOM_binary(self.afile,recno,\
                        dims=self.dims)
            mask  = np.array(1-np.isfinite(vbl),dtype='bool')
         else:
            vc    = MR.check_names('fice',self.variables)
            crec  = self.record_numbers[vc]
            vh    = MR.check_names('hice',self.variables)
            hrec  = self.record_numbers[vh]
            #xmin  = self.minvals2d[vc]
            #xmax  = self.maxvals2d[vc]
            #
            fice  = get_array_from_HYCOM_binary(self.afile,crec,\
                        dims=self.dims)
            hice  = get_array_from_HYCOM_binary(self.afile,hrec,\
                        dims=self.dims)

            vbl   = fice*hice
            mask  = np.array(1-np.isfinite(fice),dtype='bool')

      elif type(vname)==type([]):
         # 3d var if list
         vname,layer = vname
         vname       = MR.check_names(vname,self.variables3d)
         recno       = self.record_numbers3d[vname][layer]
         xmin        = self.minvals3d       [vname][layer]
         xmax        = self.maxvals3d       [vname][layer]
         vbl         = get_array_from_HYCOM_binary(self.afile,recno,\
                           dims=self.dims)

      extra_atts  = [['dimensions'],['i','j']]
      vbl         = MR.var_object(vbl,extra_atts=extra_atts)

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


      ########################################################
      if inner_points:
         extra_atts  = [['dimensions'],['i','j']]
         vbl         = MR.var_object(vbl[:-1,:-1],extra_atts=extra_atts)
      ########################################################

      return vbl
   #######################################################################


   #######################################################################
   def imshow(self,var_opts,**kwargs):
      """
      pobj   = self.imshow(var_opts,pobj=None,\
           clim=None,add_cbar=True,clabel=None,show=True,\
           test_ijs=None)
      """

      pobj   = MR.imshow(self,var_opts,**kwargs)
      return pobj
   #######################################################################


   #######################################################################
   def plot_var(self,var_opts,**kwargs):
      """
      pobj,bmap = self.plot_var(var_opts,\
         pobj=None,bmap=None,HYCOMreg='TP4',\
         clim=None,add_cbar=True,clabel=None,show=True,\
         test_lonlats=None,date_label=0)
      """

      return MR.plot_var(self,var_opts,**kwargs)
   #######################################################################


   #######################################################################
   def plot_var_pair(self,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs):
      """
      pobj,bmap=plot_var_pair(self,var_opts1,var_opts2,pobj=None,bmap=None,**kwargs)
      """
      return MR.plot_var_pair(self,var_opts1,var_opts2,**kwargs)
   #######################################################################


   ###########################################################
   def make_png(self,var_opts,**kwargs):
      """
      pobj,bmap=self.make_png(var_opts,pobj=None,bmap=None,figdir='.',date_label=True,time_index=0,**kwargs)
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
   def interp2points(self,varname,target_lonlats,**kwargs):
      """
      self.interp2points(varname,target_lonlats,**kwargs):
      *fobj is a file object eg from netcdf
      *varname is a string (name of variable in file object)
      *target_lonlats = [target_lon,target_lat], target_lon/lat are numpy arrays
      *time_index (integer) - for multi-record files
      *mapping is a pyproj.Proj object to project form lat/lon to x,y (stereographic projection)
      * kwargs for mod_reading.interp2points
      """
      return MR.interp2points(self,varname,target_lonlats,**kwargs)
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


# 
# class file_list:
#    def __init__(self,directory,pattern,extension,**kwargs):
#       import os
# 
#       self.directory = directory
#       self.extension = extension
# 
#       lst         = os.listdir(directory)
#       file_list   = []
#       for fil in lst:
# 
#          fname,fext  = os.path.splitext(fil)
# 
#          # check pattern
#          if (fext==extension) and (pattern in fname):
#             file_list.append(fil)
# 
#       wsn            = '/work/shared/nersc/msc/ModelInput'
#       gridpath_lut   = {'FR1':wsn+'/FramStrait_Hyc2.2.12/FR1a0.03-clean//topo',\
#                         'BS1':wsn+'/BS1a0.045-clean/topo',\
#                         'TP4':wsn+'/../REANALYSIS/topo'}
#       HB = False
#       if extension=='.a':
#          self.getinfo      = HYCOM_binary_info
#          HB                = True
#          self.HYCOM_region = file_list[0][:3]
#          if 'gridpath' not in kwargs:
#             kwargs.update({'gridpath':gridpath_lut[self.HYCOM_region]})
#          
#       elif extension=='.nc':
#          self.getinfo      = nc_getinfo
#          self.HYCOM_region = 'TP4' #TODO pass in HYCOMreg?
# 
#       # get main objects
#       objects     = []
#       datetimes   = []
#       time_values = []
#       for fil in file_list:
#          obj   = self.getinfo(self.directory+'/'+fil,**kwargs)
#          objects.append(obj)
#          datetimes.append(obj.datetimes[0])
#          time_values.append(obj.time_values[0])
# 
#       self.reference_date  = obj.reference_date
#       self.time_units      = obj.time_units
# 
#       # sort the time values (1st record is earliest)
#       # - also reorder objects, datetimes, file_list
#       TV                   = sorted([(e,i) for i,e in enumerate(time_values)])
#       self.time_values,ii  = np.array(TV).transpose()
#       self.objects         = [objects  [int(i)] for i in ii]
#       self.datetimes       = [datetimes[int(i)] for i in ii]
#       self.file_list       = [file_list[int(i)] for i in ii]
# 
#       self.date_strings = []
#       self.time_strings = []
#       for dtm in self.datetimes:
#          self.date_strings.append(dtm.strftime('%Y%m%d'))
#          self.time_strings.append(dtm.strftime('%H%M%S'))
# 
#       # add some methods inherited from individual objects
#       self.get_lonlat   = self.objects[0].get_lonlat
#       if HB:
#          self.get_depths   = self.objects[0].get_depths
# 
#       return
#    ###########################################################
# 
# 
#    ###########################################################
#    def make_png_all(self,var_opts,HYCOMreg=None,figdir='.',**kwargs):
# 
#       pobj        = plot_object()
#       fig,ax,cbar = pobj.get()
# 
#       if HYCOMreg is None:
#          HYCOMreg = self.HYCOM_region
#       bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
# 
#       N  = len(self.objects)
#       for i,obj in enumerate(self.objects):
# 
#          pobj,bmap   = obj.make_png(var_opts,\
#                            bmap=bmap,\
#                            figdir=figdir,show=False,**kwargs)
# 
#          ax.cla()
#          if pobj.cbar is not None:
#             pobj.cbar.ax.clear()   # cbar.ax.clear()
# 
#          print('\n'+str(i+1)+' records done out of '+str(N))
# 
#       plt.close(fig)
#       return
#    ###########################################################
# 
# 
#    ###########################################################
#    def make_png_pair_all(self,var_opts1,var_opts2,HYCOMreg=None,figdir='.',**kwargs):
# 
#       pobj        = plot_object()
#       fig,ax,cbar = pobj.get()
# 
#       if HYCOMreg is None:
#          HYCOMreg = self.HYCOM_region
#       bmap        = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
# 
#       N  = len(self.objects)
#       for i,obj in enumerate(self.objects):
# 
#          print('\n'+str(i)+' records done out of '+str(N))
# 
#          pobj,bmap   = obj.make_png_pair(var_opts1,var_opts2,\
#                         pobj=pobj,bmap=bmap,\
#                         figdir=figdir,show=False,**kwargs)
# 
#          if i==0:
#             # Fix axes position to stop it moving round
#             pobj  = pobj.renew(axpos=pobj.ax.get_position())
# 
#          ax.cla()
#          if pobj.cbar is not None:
#             pobj.cbar.ax.clear()   # cbar.ax.clear()
# 
#       plt.close(fig)
#       return
#    ###########################################################
# 
# 
#    ###########################################################
#    def compare_ice_edge_obs_all(self,HYCOMreg=None,figdir='.',**kwargs):
# 
#       pobj        = plot_object()
#       fig,ax,cbar = pobj.get()
# 
#       if HYCOMreg is None:
#          HYCOMreg = self.HYCOM_region
#       bmap  = Fplt.start_HYCOM_map(HYCOMreg,cres='i')
# 
#       if 'show' in kwargs:
#          kwargs['show'] = False
#       else:
#          kwargs.update({'show':False})
# 
#       if 'obs_type' not in kwargs:
#          obs_type = 'OSISAF'
#          kwargs.update({'obs_type':obs_type})
# 
#       N  = len(self.objects)
#       for i,obj in enumerate(self.objects):
# 
#          dtmo     = obj.datetime
#          datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
#          figname  = figdir+'/'+obj.basename+'_v'+obs_type+'_'+datestr+'.png'
#          pobj,bmap   = obj.compare_ice_edge_obs(pobj=pobj,bmap=bmap,\
#                figname=figname,**kwargs)
# 
#          ax.cla()
#          if pobj.cbar is not None:
#             pobj.cbar.ax.clear()   # cbar.ax.clear()
# 
#          print('\n'+str(i+1)+' records done out of '+str(N))
# 
#       plt.close(fig)
#       return
#    ###########################################################

######################################################################
