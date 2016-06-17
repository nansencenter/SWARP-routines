import os,sys
if 0:
   from   mpl_toolkits.basemap import pyproj
else:
   import pyproj
import mod_reading          as mr
import numpy                as np
import geometry_sphere      as GS
import ESMF_utils           as EU
import fns_plotting         as FP
import ESMF
from matplotlib import pyplot as plt

ODLsrs = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45'+\
          ' +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

###############################################
def wgs84_ellipsoid_mat():
   a,fi = 6378137,298.257223563
   # f=flattening=(a-b)/a=1/fi
   f    = 1./fi
   b    = (1-f)*a
   ecc  = np.sqrt(1-pow(b/a,2)) # eccentricity
   # print(a,fi,f,b,ecc)
   # print(a/(a-b))
   return a,ecc
###############################################


###############################################
def get_ESMF_inputs(lon,lat,ODLmap):
   X,Y      = ODLmap(lon,lat)

   # inner centres
   centres = lon[1:-1,1:-1],lat[1:-1,1:-1]
   nx,ny   = centres[0].shape # shape of centres 
   
   # corners
   # - take mean in stereographic projection
   cX  = np.zeros((nx+1,ny+1))
   cY  = np.zeros((nx+1,ny+1))
   X,Y = ODLmap(lon,lat)
   
   # mean X
   X1  = X[0:-1,0:-1]
   X2  = X[0:-1,1:]
   X3  = X[1:,1:]
   X4  = X[1:,0:-1]
   cX  = .25*(X1+X2+X3+X4)
   
   # mean Y
   Y1  = Y[0:-1,0:-1]
   Y2  = Y[0:-1,1:]
   Y3  = Y[1:,1:]
   Y4  = Y[1:,0:-1]
   cY  = .25*(Y1+Y2+Y3+Y4)
   
   # reverse to get lon,lat
   corners = ODLmap(cX,cY,inverse=True)
   ell_mat = wgs84_ellipsoid_mat()

   areas   = np.zeros(centres[0].shape)
   latC    = np.zeros((4,))
   lonC    = np.zeros((4,))
   for i in range(nx):
      for j in range(ny):
         latC[0] = lat[i,j]
         latC[1] = lat[i,j+1]
         latC[2] = lat[i+1,j+1]
         latC[3] = lat[i+1,j]
         lonC[0] = lon[i,j]
         lonC[1] = lon[i,j+1]
         lonC[2] = lon[i+1,j+1]
         lonC[3] = lon[i+1,j]
         areas[i,j] = GS.area_polygon_ellipsoid(lonC,latC,ell_mat)
         print(lonC)
         print(latC)
         print(ell_mat)
         sys.exit()

   return centres,corners,areas
###############################################


###############################################
def create_field(grid, name, mask_array=None):
  '''
  PRECONDITIONS: An Grid has been created, and 'name' is a string that 
                 will be used to initialize the name of a new Field.\n
  POSTCONDITIONS: An Field has been created.\n
  RETURN VALUES: \n Field :: field \n
  '''
  # defaults to center staggerloc
  field = ESMF.Field(grid, name)

  if mask_array is not None:
     field.data[:,:]   = mask_array


  return field
###############################################


###############################################
def run_regridding(srcfield, dstfield, srcfracfield, dstfracfield):
  '''
  PRECONDITIONS: Two Fields have been created and a regridding operation 
                 is desired from 'srcfield' to 'dstfield'.  The 'srcfracfield'
                 and 'dstfractfield' are Fields created to hold
                 the fractions of the source and destination fields which 
                 contribute to the regridding operation.\n
  POSTCONDITIONS: A regridding operation has set the data on 'dstfield', 
                  'srcfracfield', and 'dstfracfield'.\n
  RETURN VALUES: \n Field :: dstfield \n Field :: srcfracfield \n
                 Field :: dstfracfield \n
  '''
  # call the regridding functions
  print('hi')
  regridSrc2Dst = ESMF.Regrid(srcfield, dstfield, \
                              src_mask_values=np.array([1], dtype=np.int32), \
                              dst_mask_values=np.array([1], dtype=np.int32), \
                              regrid_method=ESMF.RegridMethod.CONSERVE, \
                              unmapped_action=ESMF.UnmappedAction.ERROR, \
                              src_frac_field=srcfracfield, \
                              dst_frac_field=dstfracfield)
  print('ho')
  dstfield = regridSrc2Dst(srcfield, dstfield)
  
  return dstfield, srcfracfield, dstfracfield
###############################################


###############################################
def plot_field(grid,field,bmap):

    from matplotlib import pyplot

    # get corners - better for pcolor plots
    lon = grid.get_coords(0, ESMF.StaggerLoc.CORNER)
    lat = grid.get_coords(1, ESMF.StaggerLoc.CORNER)
    
    # data
    V   = field.data[:,:]

    #plot
    pobj    = FP.plot_object()
    PC      = bmap.pcolor(lon,lat,V,latlon=True,ax=pobj.ax)
    pobj.fig.colorbar(PC)
    FP.finish_map(bmap)

    plt.show(pobj.fig)
    return pobj
###############################################


# start up ESMF
manager = ESMF.Manager(debug = True)  

# inquire for rank and proc from ESMF Virtual Machine
localPet = manager.local_pet
petCount = manager.pet_count

ODLmap  = pyproj.Proj(ODLsrs)   #same as srs="+init=EPSG:3413"
bmap    = FP.start_HYCOM_map('Arctic')

if 0:
   # ease2 projection (North)
   # cf https://nsidc.org/data/ease/versions.html
   ease2 = pyproj.Proj("+init=EPSG:6931")
   # ease1 = pyproj.Proj("+init=EPSG:3408") # works

# directories
if 0:
    #hexagon
    odir  = '/work/shared/nersc/msc/cersat/' # path to observations
    mdir  = '/work/timill/RealTime_Models/TP4a0.12/expt_01.5/data' # path to model data
    flist = mr.file_list(mdir,'DAILY','.a')
else:
    mdir    = '/mnt/sda1/work/Model-outputs/thickness_comp_ifremer/TP4/2015_060'
    odir    = '/mnt/sda1/work/Model-outputs/thickness_comp_ifremer/cersat'
    flist = mr.file_list(mdir,'DAILY','.a',gridpath=mdir+'/../topo')


if 1:
    olist=['cs2_smos_ice_thickness_20150302_20150308.nc']
else:
    olist = os.listdir(odir)

# mlon,mlat   = flist.get_corners()
Mgrid       = flist.create_ESMF_grid()
Mfld        = create_field(Mgrid,'hice')

ofil    = olist[-1]
nci     = mr.nc_getinfo(odir+'/'+ofil)
lon,lat = nci.get_lonlat()

vhobs   = 'analysis_thickness' 
hobs    = nci.get_var(vhobs)

centres,corners,areas   = get_ESMF_inputs(lon,lat,ODLmap)
# print(centres[0].shape)
# print(corners[0].shape)
# print(areas.shape)
# print(mask.shape)
Ogrid       = EU.create_ESMF_grid(centres,corners,AREA=areas,MASK=hobs.values.mask[1:-1,1:-1])
Ofld        = create_field(Ogrid,'hice')
Mfrac_fld   = create_field(Mgrid,'src_frac') # frac of model grid that contributes to regridding
Ofrac_fld   = create_field(Ogrid,'dst_frac') # frac of observation grid that contributes to regridding


if 0:
   # check grid
   xc    = nci.get_var('xc')
   yc    = nci.get_var('yc')
   XC,YC = np.meshgrid(xc,yc)

   lon0,lat0   = ease2(XC[0,0],YC[0,0],inverse=True)
   print(lon0,lat0)

for ofil in olist:
   cdate1   = ofil[-20:-12]
   cdate2   = ofil[-11:-3]
   hav_m    = flist.time_average('hice',\
                start_date=cdate1+'T120000Z',\
                end_date=cdate2+'T120000Z',\
                inner_points=True)

   Mfld.data[:,:]           = 1*hav_m.data
   Mfld.data[hav_m.mask]    = 0.
   if 1:
      plot_field(Mgrid,Mfld,bmap)
   Ofld,Mfrac_fld,Ofrac_fld = run_regridding(Mfld,Ofld,Mfrac_fld,Ofrac_fld)
   plot_field(Ogrid,Ofld,bmap)
