import ESMF

#########################################################################
def create_ESMF_grid(centres,corners,areas=None,mask=None,\
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
      

   maxIndex       = list(centres[0].shape) # no of centres
   Ndim           = len(maxIndex)
   staggerlocs    = [ESMF.StaggerLoc.CORNER,ESMF.StaggerLoc.CENTER]

   # # VM - needed?
   # vm = ESMF.ESMP_VMGetGlobal()
   # localPet, petCount = ESMF.ESMP_VMGet(vm)

   print('\ncreate_ESMF_grid\n')
   Egrid = ESMF.Grid(maxIndex, numPeriDims=numPeriDims, coordSys=coordSys,\
                     coordTypeKind=coordTypeKind, staggerlocs=staggerlocs)


   # ========================================================
   # CORNERS
   # get the coordinate pointers and set the coordinates
   for x in range(Ndim):
      gridCorner        = Egrid.get_coords(x, ESMF.StaggerLoc.CORNER)
      gridCorner[:,:]   = corners[x]
   # ========================================================


   # ========================================================
   # CENTERS
   # get the coordinate pointers and set the coordinates
   for x in range(Ndim):
      gridCenter        = Egrid.get_coords(x, ESMF.StaggerLoc.CENTER)
      gridCenter[:,:]   = centres[x]
   # ========================================================


   # ========================================================
   if areas is not None:
     # add areas of grid cells
     grid.add_item(ESMF.GridItem.AREA)
     area      = Egrid.get_item(ESMF.GridItem.AREA) # pointer
     area[:,:] = self.land_mask(inner_points=True)
   # ========================================================


   # ========================================================
   if mask is not None:
     # set up the grid mask
     grid.add_item(ESMF.GridItem.MASK)
     mask      = Egrid.get_item(ESMF.GridItem.MASK) # pointer
     mask[:,:] = self.land_mask(inner_points=True)
   # ========================================================
   
   return Egrid
###################################################################
