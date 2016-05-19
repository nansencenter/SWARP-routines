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

      self.afile        = gridpath+'/regional.grid.a'
      self.bfile        = gridpath+'/regional.grid.b'
      self.afile_depth  = gridpath+'/regional.depths.a'
      self.bfile_depth  = gridpath+'/regional.depths.b'

      # get grid size & record numbers
      self._read_bfile()

      return
   ###################################################################


   ###################################################################
   def _read_bfile(self):
      bid   = open(self.bfile,'r')
      lines = bid.readlines()
      bid.close()

      # get grid size
      self.Nx  = int(lines[0].split()[0])
      self.Ny  = int(lines[1].split()[0])

      R     = {}
      for i,lin in enumerate(lines[3:]):
         word  = lin.split()[0] # 1st word in line 
         R.update({word,i})

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

      grid is arranged (Arakawa C-grid):
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
   def get_grid_sizes(self,inner_points=False):
      """
      call self.get_centres(inner_points=False)

      grid is arranged (Arakawa C-grid):
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
      scux  = self.get_array('scux')
      scvy  = self.get_array('scvy')
      if inner_points:
         return scux[:-1,:-1],scvy[:-1,:-1]
      else:
         return scux,scvy
   ###################################################################

   ###################################################################
   def get_areas(self,inner_points=False):
      """
      call self.get_centres(inner_points=False)

      grid is arranged (Arakawa C-grid):
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
      scux,scvy   = self.get_grid_sizes(inner_points=inner_points)
      return scux*scvy
   ###################################################################


   ###################################################################
   def get_depth(self,inner_points=False):
     dep = get_array_from_HYCOM_binary(self.afile_depth,1)
     if inner_points:
        return dep[:-1,:-1]
     else:
        return dep
   ###################################################################


   ###################################################################
   def land_mask(self,inner_points=False):
      dep   = self.get_depth(self,inner_points=inner_points)
      return (dep==0.) # mask land out
   ###################################################################


   ###################################################################
   def create_ESMF_grid(self,do_mask):
      import ESMF
      maxIndex       = [self.Nx-1,self.Ny-1] # no of centres
      coordTypeKind  = 'r8'                  # real double
      coordSys       = ESMF.CoordSys.SPH_DEG      # lon,lat (degrees) or SPH_DEG or CART
      numPeriDims    = 1                     # no of periodic dimensions (lon)
      staggerlocs    = [ESMF.StaggerLoc.CORNER,ESMF.StaggerLoc.CENTER]

      Egrid = ESMF.Grid(maxIndex, numPeriDims=numPeriDims, coordSys=coordSys,\
                        coordTypeKind=coordTypeKind, staggerlocs=staggerlocs)

      # # VM - needed?
      # vm = ESMF.ESMP_VMGetGlobal()
      # localPet, petCount = ESMF.ESMP_VMGet(vm)
  
      # ========================================================
      # CORNERS
      # get the coordinate pointers and set the coordinates
      [x,y]       = [0, 1]
      gridXCorner = Egrid.get_coords(x, ESMF.StaggerLoc.CORNER)
      gridYCorner = Egrid.get_coords(y, ESMF.StaggerLoc.CORNER)

      qlon,qlat         = self.get_corners()
      gridXCorner[:,:]  = qlon
      gridYCorner[:,:]  = qlat
      # ========================================================
  
  
      # ========================================================
      # CENTERS
      # get the coordinate pointers and set the coordinates
      [x,y]       = [0, 1]
      gridXCenter = Egrid.get_coords(x, ESMF.StaggerLoc.CENTER)
      gridYCenter = Egrid.get_coords(y, ESMF.StaggerLoc.CENTER)
  
      plon,plat         = self.get_centres(inner_points=True)
      gridXCenter[:,:]  = plon
      gridYCenter[:,:]  = plat
      # ========================================================


      # ========================================================
      if do_mask:
        # set up the grid mask
        grid.add_item(ESMF.GridItem.MASK)
        mask      = Egrid.get_item(ESMF.GridItem.MASK) # pointer
        mask[:,:] = self.land_mask()
      # ========================================================
      

      return Egrid
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
      if 'DAILY' in self.basename:
         # file written at end of day - change to midday
         self.time_value   = self.time_value -.5 # model time (days)

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
   def get_var(self,vname,time_index=None):
      """
      vbl=get_var(vname,time_index=None) - vname is a string of 2d variable name (surface variable name,layer=0)
      vbl=get_var([vname,layer],time_index=None) - vname is a string of 2d or 3d variable name (layer=0 is surface, layer=1 to kdm are ocean layers)
      *time_index is not used - only a place-holder for some routines which handle netcdf as well
      *vbl is a mod_reading.var_opt class (vbl.values is type np.ma.array - eg land is masked)
      """

      if type(vname)!=type([]):
         # 2d var
         layer = 0
         vname = MR.check_names(vname,self.variables)
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


#    ###########################################################
#    def MIZmap(self,var_name='dmax',vertices=None,\
#          do_sort=False,EastOnly=True,plotting=True,**kwargs):
#       """
#       Call  : self.MIZmap(var_name='dmax',vertices=None,\
#                   do_sort=False,EastOnly=True,plotting=True,**kwargs):
#       Inputs:
#          var_name is variable to find MIZ from
#          **kwargs to be passed onto MIZchar.get_MIZ_poly:
#             outdir='.' - where to save results
#             mask_corners=None - can mask out smaller region
#                               - corners = [lons,lats], where lons=[DL,DR,UR,UL]
#       Returns: MIZchar.MIZpoly object
#       """
# 
#       import MIZchar as mc
# 
#       vname = check_names(var_name,self.variables)
#       if var_name == 'dmax':
#          # FSD MIZarray(1-
#          Arr         = self.get_var(vname)
#          clim        = [0,300]# for plotting
#          lower_limit = .1     # for plotting
#       elif var_name == 'fice':
#          # conc MIZ
#          Arr         = self.get_var(vname)
#          clim        = [0,1]  # for plotting
#          lower_limit = .15    # for plotting
#       elif var_name == 'hice':
#          # thin ice areas
#          Arr         = self.get_var(vname)
#          clim        = [0,2.] # for plotting
#          lower_limit = .01    # for plotting
#       else:
#          raise ValueError('Wrong selection variable for MIZmap')
# 
#       print("MIZchar.get_MIZ_poly\n")
#       lon,lat  = self.get_lonlat()
#       MPdict   = {}
#       tfiles   = {}
# 
#       if vertices is not None:
#          do_sort  = False
#       
#       if do_sort:
#          # possible regions are:
#          regions  = ['gre','bar','beau','lab','balt','les','can']
# 
#          if EastOnly:
#             # concentrate on the eastern Arctic
#             # (and forget Baltic Sea)
#             regions.remove('balt' )
#             regions.remove('les' )
#             regions.remove('can' )
#             regions.remove('beau')
# 
#          # for reg in ['gre']:
#          for reg in regions:
#             mp = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,region=reg)
#             MPdict.update({reg:mp})
# 
#             fname0   = self.basename+'_'+var_name +'_'+reg
#             tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
#             if 'all' in tfile.keys():
#                tfiles.update({reg:tfile['all']})
# 
#          if 0:
#             MPdict['gre'].show_maps()
#             return MPdict
# 
#       else:
#          reg   = 'all'
#          mp    = mc.get_MIZ_poly(Arr.values,lon,lat,var_name=var_name,vertices=vertices)
#          MPdict.update({reg:mp})
#          #
#          fname0   = self.basename+'_'+var_name
#          tfile    = mp.write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
#          if 'all' in tfile.keys():
#             tfiles.update({reg:tfile['all']})
# 
#       Pdict    = {}
#       PLOTTING = False
#       for reg in tfiles.keys():
# 
#          ##########################################################
#          # filenames
#          tfil     = tfiles[reg]                          # text file with polygon outlines characterized
#          figname  = tfil.replace('.txt','.png')          # plot of polygons
#          shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
#          sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
#          ##########################################################
# 
# 
#          if vertices is None:
#             ##########################################################
#             if do_sort:
#                mapreg   = reg
#             else:
#                mapreg   = self.HYCOM_region
#             ##########################################################
# 
# 
#             ##########################################################
#             # process each text file to get MIZ width etc
#             print("MIZchar.single_file: "+tfil+"\n")
#             bmap     = Fplt.start_HYCOM_map(mapreg)
#          else:
#             vlons,vlats = np.array(vertices).transpose()
#             vx,vy       = GS.polar_stereographic_simple(vlons,vlats,NH=True,inverse=False)
#             xy_verts    = [(vx[i],vy[i]) for i in range(len(vx))]
# 
#             # get approximate centre and limits for basemap
#             vP          = shg.Polygon(xy_verts)
#             width       = 5*np.sqrt(vP.area)
#             height      = 5*np.sqrt(vP.area)
#             xcen,ycen   = vP.centroid.coords.xy
#             lonc,latc   = GS.polar_stereographic_simple(np.array([xcen]),np.array([ycen]),\
#                            NH=True,inverse=True)
#             # make basemap
#             bmap        = Basemap(projection='stere',lon_0=lonc,lat_0=latc,\
#                                     width=width,height=height,\
#                                     resolution='i')
# 
#          Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
#          Pdict.update({reg:Psolns})
#          
#          # Save summary & shapefile
#          mc.save_summary  (Psolns,sumname)
#          mc.save_shapefile(Psolns,filename=shpname)
#          ##########################################################
# 
#          
#          if plotting:
#             ##########################################################
#             # Make plot
#             var_opts = make_plot_options(vname,lower_limit=lower_limit)
#             pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=clim)[0]
#             fig      = pobj.fig
#             ax       = pobj.ax
#             PLOTTING = True
# 
#             for MIZi in Psolns:
#                # plot outlines of polygons
#                lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
#                bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)
# 
#                Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
#                if Wavg>26:
#                   MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)
# 
#                   # add text with mean width
#                   xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
#                   xav                  = (xmin+xmax)/2.
#                   ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
#                      color='k',fontsize=16,horizontalalignment='right',\
#                      verticalalignment='top')
# 
#             Fplt.finish_map(bmap)
#             print('Saving '+figname)
#             fig.savefig(figname)
#             # plt.show(fig)
#             ax.cla()
#             fig.clear()
#             # finished region
#             ##########################################################
# 
#       if PLOTTING:
#          plt.close(fig)
#       return mp,Pdict,tfiles
#    ###########################################################
#    
#    
#    ###########################################################
#    def areas_of_disagreement(self,obs_type='OSISAF',do_sort=True,EastOnly=True,plotting=True,**kwargs):
#       # kwargs: outdir='.',do_sort=True
# 
#       import MIZchar as mc
# 
#       if obs_type == 'OSISAF':
#          var_name = 'fice'
#          bmap     = basemap_OSISAF()
#          obsfil   = '/work/shared/nersc/msc/OSI-SAF/'+\
#                self.datetime.strftime('%Y')+'_nh_polstere/'+\
#                'ice_conc_nh_polstere-100_multi_'+\
#                self.datetime.strftime('%Y%m%d')+'1200.nc'
#       else:
#          raise ValueError('Wrong selection variable for areas_of_disagreement')
# 
#       # observation grid & compared quantity
#       nci         = nc_getinfo(obsfil)
#       lon2,lat2   = nci.get_lonlat()
#       Xobs,Yobs   = bmap(lon2,lat2)
#       Zobs        = nci.get_var(var_name)
# 
#       # model grid & compared quantity
#       Zmod        = self.get_var(var_name)
#       lon,lat     = self.get_lonlat()
#       Xmod,Ymod   = bmap(lon,lat)
# 
#       if '%' in Zobs.units:
#          conv_fac = .01
#       else:
#          conv_fac = 1
# 
#       if 1:
#          #Zref,Zint should be np.ma.array
#          lon_ref,lat_ref   = lon2,lat2
#          Xref,Yref,Zref    = Xobs,Yobs,conv_fac*Zobs.values # obs grid is reference;                 
#          Xint,Yint,Zint    = Xmod,Ymod,Zmod.values          # to be interped from model grid onto obs grid;  Zint is np.ma.array
# 
#       # add the mask for the ref to output (Arr)
#       Arr   = reproj_mod2obs(Xint,Yint,Zint,Xref,Yref,mask=1*Zref.mask)
# 
#       # add the mask for Arr to Zref
#       Zref  = np.ma.array(Zref.data,mask=Arr.mask)
# 
#       MPdict   = {'Over':{},'Under':{}}
#       tfiles   = {'Over':{},'Under':{}}
# 
#       if 0:
#          # test interpolation and matching of masks
#          fig   = plt.figure()
#          ax1   = fig.add_subplot(1,2,1)
#          ax1.imshow(Arr.transpose(),origin='upper')
#          ax2   = fig.add_subplot(1,2,2)
#          ax2.imshow(Zref.transpose(),origin='upper')
#          plt.show(fig)
#          return Xint,Yint,Zint,Xref,Yref,Zref
# 
#       if do_sort:
#          # possible regions are:
#          regions  = ['gre','bar','beau','lab','balt','les','can']
# 
#          if EastOnly:
#             # concentrate on the eastern Arctic
#             # (and forget Baltic Sea)
#             regions.remove('balt' )
#             regions.remove('les' )
#             regions.remove('can' )
#             regions.remove('beau')
# 
#          # for reg in ['bar']:
#          for reg in regions:
# 
#             # Arr,Zref are np.ma.array objects
#             Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref,region=reg)
#             MPdict['Over'] .update({reg:Over})
#             MPdict['Under'].update({reg:Under})
# 
#             for OU in ['Over','Under']:
# 
#                fname0   = self.basename+'_v'+obs_type +'_'+OU+'_'+reg
#                tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
#                if 'all' in tfile.keys():
#                   tfiles[OU].update({reg:tfile['all']})
# 
#          if 0:
#             # Show the 4 maps for last region (for each over/under)
#             MPdict['Over'] [reg].show_maps()
#             MPdict['Under'][reg].show_maps()
#             return MPdict
#       else:
#          reg         = 'all'
#          Over,Under  = mc.get_AOD_polys(Arr,Zref,lon_ref,lat_ref)
#          MPdict['Over'] .update({reg:Over})
#          MPdict['Under'].update({reg:Under})
# 
#          for OU in ['Over','Under']:
# 
#             fname0   = self.basename+'_v'+obs_type+'_'+OU
#             tfile    = MPdict[OU][reg].write_poly_stats(filename_start=fname0,do_sort=False,**kwargs)
#             if 'all' in tfile.keys():
#                tfiles[OU].update({reg:tfile['all']})
# 
#          if 0:
#             # Show the 4 maps (for over/under)
#             # MPdict['Over'] ['all'].show_maps()
#             MPdict['Under']['all'].show_maps()
#             return MPdict
# 
#       print(tfiles)
#       print(MPdict)
#       Pdict = {'Over':{},'Under':{}}
#       for OU in ['Over','Under']:
#          PLOTTING = False
#          for reg in tfiles[OU].keys():
# 
#             ##########################################################
#             # filenames
#             tfil     = tfiles[OU][reg]                      # text file with polygon outlines characterized
#             figname  = tfil.replace('.txt','.png')          # plot of polygons
#             shpname  = tfil.replace('.txt','.shp')          # save polygons to shapefile with characteristics eg MIZ width
#             sumname  = tfil.replace('.txt','_summary.txt')  # save average MIZ width etc to summary file
#             ##########################################################
# 
# 
#             ##########################################################
#             if do_sort:
#                mapreg   = reg
#             else:
#                mapreg   = self.HYCOM_region
#             ##########################################################
# 
# 
#             ##########################################################
#             # process each text file to get MIZ width etc
#             print("MIZchar.single_file: "+tfil+"\n")
#             bmap     = Fplt.start_HYCOM_map(mapreg)
#             Psolns   = mc.single_file(tfil,bmap,MK_PLOT=False,METH=5)
#             Pdict[OU].update({reg:Psolns})
#             
#             # Save summary & shapefile
#             mc.save_summary  (Psolns,sumname)
#             mc.save_shapefile(Psolns,filename=shpname)
#             ##########################################################
# 
#             
#             if plotting:
#                ##########################################################
#                # Make plot
#                var_opts = make_plot_options('fice',ice_mask=True)
#                pobj     = self.plot_var(var_opts,bmap=bmap,show=False,clim=[0,1])[0]
#                fig      = pobj.fig
#                ax       = pobj.ax
#                PLOTTING = True
# 
#                for MIZi in Psolns:
#                   # plot outlines of polygons
#                   lon,lat  = np.array(MIZi.ll_bdy_coords).transpose()
#                   bmap.plot(lon,lat,latlon=True,ax=ax,color='k',linewidth=2.5)
# 
#                   Wavg  = MIZi.record['Width_mean']/1.e3 # mean width in km
#                   if Wavg>26:
#                      MIZi.plot_representative_lines(bmap,ax=ax,color='k',linewidth=1.5)
# 
#                      # add text with mean width
#                      xmin,xmax,ymin,ymax  = MIZi.bbox(bmap)
#                      xav                  = (xmin+xmax)/2.
#                      ax.text(xmax,ymin,'%4.1f km' %(Wavg),\
#                         color='k',fontsize=16,horizontalalignment='right',\
#                         verticalalignment='top')
# 
#                Fplt.finish_map(bmap)
#                print('Saving '+figname)
#                fig.savefig(figname)
#                # plt.show(fig)
#                ax.cla()
#                fig.clear()
#                # finished region
#                ##########################################################
# 
#          if PLOTTING:
#             plt.close(fig)
# 
#       return MPdict,tfiles,Pdict
#    ###########################################################
# 
# 
#    ###########################################################
#    def compare_ice_edge_obs(self,pobj=None,bmap=None,\
#          obs_type='OSISAF',date_label=1,figname=None,**kwargs):
# 
#       var_opts1   = make_plot_options('ficem',\
#          vec_opt=0,conv_fac=1,wave_mask=False,ice_mask=True,dir_from=True)
#       var_opts1   = check_var_opts(var_opts1,self.variables)
# 
#       if 'show' in kwargs:
#          show           = kwargs['show']
#          kwargs['show'] = False
#          pobj,bmap      = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,**kwargs)
#       else:
#          show        = True
#          pobj,bmap   = self.plot_var(var_opts1,pobj=pobj,bmap=bmap,show=False,**kwargs)
# 
#       fig,ax,cbar = pobj.get()
# 
#       dtmo  = self.datetimes[0]
#       if obs_type=='OSISAF':
#          obsfil   = '/work/shared/nersc/msc/OSI-SAF/'+\
#                dtmo.strftime('%Y')+'_nh_polstere/'+\
#                'ice_conc_nh_polstere-100_multi_'+\
#                dtmo.strftime('%Y%m%d')+'1200.nc'
#       else:
#          raise ValueError('invalid value of obs_type: '+obs_type)
# 
#       print(obsfil)
#       obs      = nc_getinfo(obsfil)
#       fice     = obs.get_var('ice_conc')
#       lon,lat  = obs.get_lonlat()
#       bmap.contour(lon,lat,fice.values[:,:],[15],colors='g',\
#             linewidths=2,ax=ax,latlon=True)
# 
#       if 'HYCOMreg' in kwargs:
#          reg   = kwargs['HYCOMreg']
#       else:
#          reg   = self.HYCOM_region
#          kwargs.update({'HYCOMreg':reg})
# 
# 
#       #############################################################
#       if reg=='TP4':
#          xyann = (0.05,.925)
#       else:
#          xyann = (0.4,.925)
# 
# 
#       if date_label==1:
#          # date only
#          tlabel   = dtmo.strftime('%d %b %Y')
#          pobj.ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18)
#       elif date_label==2:
#          # date + time
#          tlabel   = dtmo.strftime('%d %b %Y %H:%M')
#          pobj.ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)
#       elif type(date_label)==type('hi'):
#          # manual label
#          pobj.ax.annotate(date_label,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)
#       #############################################################
# 
# 
#       if figname is not None:
#          print('Saving to '+figname)
#          fig.savefig(figname)
# 
#       if show:
#          # fig.show()
#          plt.show(fig)
# 
#       return pobj,bmap,obsfil
#    ###########################################################
# 
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
