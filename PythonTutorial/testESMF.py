def create_grid(bounds, do_mask):
  '''
  PRECONDITIONS: ESMPy has been initialized, 'bounds' contains the number of
                 indices required for the first two dimensions of a Grid.
                 'do_mask' is a boolean value that gives the option to put
                 a mask on this Grid.\n
  POSTCONDITIONS: An Grid has been created.\n
  RETURN VALUES: \n Grid :: grid \n
  '''

  nx = float(bounds[0])
  ny = float(bounds[1])

  dx = 360.0/nx
  dy = 180.0/ny

  DEG2RAD = 3.141592653589793/180.0

  maxIndex = np.array([nx,ny], dtype=np.int32)

  staggerLocs = [ESMF.StaggerLoc.CORNER, ESMF.StaggerLoc.CENTER]
  grid = ESMF.Grid(maxIndex, numPeriDims=1, staggerlocs=staggerLocs)
  # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_6_3_0rp1/esmpy_doc/html/grid.html
  # ESMF.Grid(maxIndex, numPeriDims=0, coordSys=None, coordTypeKind=None, staggerlocs=None)
  # >numPeriDims   = number of periodic dimensions
  # >coordSys      = ESMF.CoordSys.CART,SPH_DEG,SPH_RAD
  # *CART: 3d cartesian?
  # *SPH_DEG: lon,lat (degrees) --> DEFAULT
  # *SPH_RAD: lon,lat (radians)
  # >coordTypeKind = ESMF.TypeKind.I4,I8,R4,R8
  # *default=R8
  # >staggerlocs   = list of places where you want grid to be defined

  # VM
  vm = ESMF.ESMP_VMGetGlobal()
  localPet, petCount = ESMF.ESMP_VMGet(vm)
  
 # get the coordinate pointers and set the coordinates
  [x,y] = [0, 1]
  gridXCorner = grid.get_coords(x, ESMF.StaggerLoc.CORNER)
  gridYCorner = grid.get_coords(y, ESMF.StaggerLoc.CORNER)
  
  for i in xrange(gridXCorner.shape[x]):
    gridXCorner[i, :] = float(i)*dx - 180.0
 
  for j in xrange(gridYCorner.shape[y]):
    gridYCorner[:, j] = float(j)*dy - 90.0
  
  ##   CENTERS

  # get the coordinate pointers and set the coordinates
  [x,y] = [0, 1]
  gridXCenter = grid.get_coords(x, ESMF.StaggerLoc.CENTER)
  gridYCenter = grid.get_coords(y, ESMF.StaggerLoc.CENTER)
  
  for i in xrange(gridXCenter.shape[x]):
    gridXCenter[i, :] = float(i)*dx + 0.5*dx - 180.0

  for j in xrange(gridYCenter.shape[y]):
    y = (float(j)*dy - 90.0)
    yp1 = (float(j+1)*dy - 90.0)
    gridYCenter[:, j] = (y+yp1)/2.0

  mask = 0
  if do_mask:
    # set up the grid mask
    grid.add_item(ESMF.GridItem.MASK)
    mask = grid.get_item(ESMF.GridItem.MASK)
    
    [x,y] = [0, 1]
    for i in range(mask.shape[x]):
      if (i == 2.0):
        mask[i, :] = 1
      else:
        mask[i, :] = 0

  return grid

import fns_plotting as Fplt
import mod_reading  as mr
from matplotlib import pyplot as plt

gdir  = '/work/shared/nersc/msc/ModelInput/../REANALYSIS/topo'
afil  = gdir+'/regional.grid.a'

plon  = mr.get_array_from_HYCOM_binary(afil,1)
plat  = mr.get_array_from_HYCOM_binary(afil,2)
qlon  = mr.get_array_from_HYCOM_binary(afil,3)
qlat  = mr.get_array_from_HYCOM_binary(afil,4)
ulon  = mr.get_array_from_HYCOM_binary(afil,5)
ulat  = mr.get_array_from_HYCOM_binary(afil,6)
vlon  = mr.get_array_from_HYCOM_binary(afil,7)
vlat  = mr.get_array_from_HYCOM_binary(afil,8)

if 0:
   bmap  = Fplt.start_HYCOM_map('TP4')
   bmap.plot(plon[:-1,:-1],plat[:-1,:-1],'.k',latlon=True)
   bmap.plot(qlon,qlat,'or',latlon=True)
   bmap.plot(ulon,ulat,'^g',latlon=True)
   bmap.plot(vlon,vlat,'vb',latlon=True)
   plt.show()
