import ESMF
import numpy as np

#########################################################################
def create_ESMF_grid(centres,corners,AREA=None,MASK=None,\
      coordTypeKind='r8',coordsys=0,numPeriDims=0):
   """
   create_ESMF_grid(centres,corners,mask=None,\
      coordTypeKind='r8',coordSys=0,numPeriDims=0)
   inputs:
   *coordSys=0:lon,lat (degrees); 1:lon,lat (radians); 2:x,y,z
   *coordTypeKind eg 'r8':real double
   *numPeriDims = number of periodic dimensions
   """
   if coordsys==0:
      coordSys=ESMF.CoordSys.SPH_DEG #lon,lat (degrees)
   elif coordsys==1:
      coordSys=ESMF.CoordSys.SPH_RAD #lon,lat (radians)
   elif coordsys==2:
      coordSys=ESMF.CoordSys.CART #x,y,z   
      

   maxIndex       = np.array(centres[0].shape) # no of centres
   Ndim           = len(maxIndex)
   staggerlocs    = [ESMF.StaggerLoc.CORNER,ESMF.StaggerLoc.CENTER]

   # # VM - needed?
   # vm = ESMF.ESMP_VMGetGlobal()
   # localPet, petCount = ESMF.ESMP_VMGet(vm)

   print('\ncreate_ESMF_grid\n')
   Egrid = ESMF.Grid(maxIndex,\
           num_peri_dims=numPeriDims,\
           coord_sys=coordSys,\
           staggerloc=staggerlocs)#,\
           # TypeKind=coordTypeKind)


   # ========================================================
   # CENTERS
   # get the coordinate pointers and set the coordinates
   for x in range(Ndim):
      gridCenter        = Egrid.get_coords(x, ESMF.StaggerLoc.CENTER)
      # print('centre array: '+str(x))
      # print(centres[x].shape)
      # print(gridCenter[:,:].shape)
      gridCenter[:,:]   = centres[x]
   # ========================================================


   # ========================================================
   # CORNERS
   # get the coordinate pointers and set the coordinates
   for x in range(Ndim):
      gridCorner        = Egrid.get_coords(x, ESMF.StaggerLoc.CORNER)
      # print('corner array: '+str(x))
      # print(corners[x].shape)
      # print(gridCorner[:,:].shape)
      gridCorner[:,:]   = corners[x]
   # ========================================================


   # ========================================================
   if AREA is not None:
     # add areas of grid cells
     Egrid.add_item(ESMF.GridItem.AREA)
     area      = Egrid.get_item(ESMF.GridItem.AREA) # pointer
     area[:,:] = AREA
   # ========================================================


   # ========================================================
   if MASK is not None:
     # set up the grid mask
     Egrid.add_item(ESMF.GridItem.MASK)
     mask      = Egrid.get_item(ESMF.GridItem.MASK) # pointer
     mask[:,:] = MASK
   # ========================================================
   
   return Egrid
###################################################################
