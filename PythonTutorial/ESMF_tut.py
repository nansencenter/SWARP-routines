import os,sys
from   mpl_toolkits.basemap import pyproj
import mod_reading          as mr
import numpy                as np
import geometry_sphere      as GS
import ESMF_utils           as EU
import fns_plotting         as FP
import ESMF

def create_grid(bounds, domask):
  '''
  PRECONDITIONS: ESMPy has been initialized, 'bounds' contains the number of
                 indices required for the first two dimensions of a Grid.
                 'domask' is a boolean value that gives the option to put
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
  grid = ESMF.Grid(maxIndex, num_peri_dims=1, staggerloc=staggerLocs)

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
  if domask:
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

def create_field(grid, name):
  '''
  PRECONDITIONS: An Grid has been created, and 'name' is a string that 
                 will be used to initialize the name of a new Field.\n
  POSTCONDITIONS: An Field has been created.\n
  RETURN VALUES: \n Field :: field \n
  '''
  # defaults to center staggerloc
  field = ESMF.Field(grid, name)

  return field



def build_analyticfield(field, grid):
  '''
  PRECONDITIONS: An Field has been created as 'field'.  'grid' has
                 been created and coordinates have been set on both the
                 center and corner stagger locations. \n
  POSTCONDITIONS: The 'field' has been initialized to an analytic field.\n
  RETURN VALUES: \n Field :: field \n
  '''
  DEG2RAD = 3.141592653589793/180.0

  # get the grid bounds and coordinate pointers
  exLB = grid.lower_bounds[ESMF.StaggerLoc.CENTER]
  exUB = grid.upper_bounds[ESMF.StaggerLoc.CENTER]

  # get the coordinate pointers and set the coordinates
  [x,y] = [0, 1]
  gridXCoord = grid.get_coords(x, ESMF.StaggerLoc.CENTER)
  gridYCoord = grid.get_coords(y, ESMF.StaggerLoc.CENTER)
  
  for i in range(exLB[x], exUB[x]): 
    for j in range(exLB[y], exUB[y]): 
      theta = DEG2RAD*gridXCoord[i, j]
      phi = DEG2RAD*(90.0 - gridYCoord[i, j])
      field.data[i, j] = 2.0 + np.cos(theta)**2 * np.cos(2.0*phi)

  return field

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
  regridSrc2Dst = ESMF.Regrid(srcfield, dstfield, \
                              src_mask_values=np.array([1], dtype=np.int32), \
                              dst_mask_values=np.array([1], dtype=np.int32), \
                              regrid_method=ESMF.RegridMethod.CONSERVE, \
                              unmapped_action=ESMF.UnmappedAction.ERROR, \
                              src_frac_field=srcfracfield, \
                              dst_frac_field=dstfracfield)
  dstfield = regridSrc2Dst(srcfield, dstfield)
  
  return dstfield, srcfracfield, dstfracfield

def compute_mass(valuefield, areafield, fracfield, dofrac):
  '''
  PRECONDITIONS: Two Fields have been created and initialized.  'valuefield'
                 contains data values of a field built on the cells of a grid, 
                 'areafield' contains the areas associated with the grid 
                 cells, and 'fracfield' contains the fractions of each cell
                 which contributed to a regridding operation involving
                 'valuefield.  'dofrac' is a boolean value that gives the
                 option to not use the 'fracfield'.\n
  POSTCONDITIONS: The mass of the data field is computed and returned.\n
  RETURN VALUES: integer :: mass \n
  '''

  [x, y] = [0, 1]
  areafield.get_area()
  mass = 0.0
  frac = 0
  for i in range(valuefield.data.shape[x]):
    for j in range(valuefield.data.shape[y]):
      if dofrac:
        mass += areafield.data[i, j] * valuefield.data[i, j] * \
                fracfield.data[i, j]
      else:
        mass += areafield.data[i, j] * valuefield.data[i, j]
  
  return mass

def compare_fields(interp_field, exact_field, dstfracfield, srcmass, dstmass, \
                   parallel):
  '''
  PRECONDITIONS: 'interp_field' is a Field that holds the values resulting
                  from a regridding operation, 'exact_field' is a Field
                  containing the values of the exact solution that is 
                  expected, and 'dstfracfield' is a Field containing the
                  fractions of the 'interp_field' which contributed to the
                  regridding product.  'srcmass' and 'dstmass' are the 
                  mass values for the source and destination data fields.
                  'parallel' is an internal variable used to determine if this
                  run should be done in serial or parallel.\n
  POSTCONDITIONS: The interpolation accuracy of a regridding operation is
                  determined by comparing the 'interp_field' to the 
                  'exact_field'.  The mass conservation is validated by
                  comparing the 'srcmass' to the 'dstmass'.\n
  RETURN VALUES: None \n
  '''
  
  vm = ESMF.ESMP_VMGetGlobal()
  localPet, _ = ESMF.ESMP_VMGet(vm)
  
  if (interp_field.data.shape != exact_field.data.shape):
    raise TypeError('compare_fields: Fields must be the same size!')

  # initialize to True, and check for False point values
  [x, y] = [0, 1]
  total_error = 0.0
  max_error = 0.0
  min_error = 1000000.0
  for i in range(interp_field.data.shape[x]):
    for j in range(interp_field.data.shape[y]):
      if (exact_field.data[i, j] != 0.0):
        err = abs(interp_field.data[i, j] / dstfracfield.data[i, j] - \
              exact_field.data[i, j]) / abs(exact_field.data[i, j])
      else:
        err = abs(interp_field.data[i, j] / dstfracfield.data[i, j] - \
              exact_field.data[i, j])
      total_error = total_error + err
      if (err > max_error): 
        max_error = err
      if (err < min_error): 
        min_error = err

  if parallel:
    # use mpi4py to collect values
    from mpi4py import MPI
 
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
 
    max_error_global = comm.reduce(max_error, op=MPI.MAX)
    min_error_global = comm.reduce(min_error, op=MPI.MIN)
    total_error_global = comm.reduce(total_error, op=MPI.SUM)
    srcmass_global = comm.reduce(srcmass, op=MPI.SUM)
    dstmass_global = comm.reduce(dstmass, op=MPI.SUM)
 
    if rank == 0:
      # check the mass
      csrv = False
      csrv_error = abs(dstmass_global - srcmass_global)/srcmass_global
      if (csrv_error < 10e-12):
        csrv = True
      itrp = False
      if (max_error_global < 10E-2):
        itrp = True
    
      if (itrp and csrv):
        print "PASS"
      else:
        print "FAIL"
      print "       Total error = "+str(total_error_global)
      print "       Max error   = "+str(max_error_global)
      print "       Min error   = "+str(min_error_global)
      print "       Csrv error  = "+str(csrv_error)
      print "       srcmass     = "+str(srcmass_global)
      print "       dstmass     = "+str(dstmass_global)
  else:
    # check the mass
    csrv = False
    csrv_error = abs(dstmass - srcmass)/srcmass
    if (csrv_error < 10e-12):
      csrv = True
    itrp = False
    if (max_error < 10E-2):
      itrp = True
  
    if (itrp and csrv):
      print "PASS"
    else:
      print "FAIL"
    print "       Total error = "+str(total_error)
    print "       Max error   = "+str(max_error)
    print "       Min error   = "+str(min_error)
    print "       Csrv error  = "+str(csrv_error)
    print "       srcmass     = "+str(srcmass)
    print "       dstmass     = "+str(dstmass)
  
  return

# opening remarks
print "\nWelcome to the Field regridding tutorial!\n"

# start up ESMF
manager = ESMF.Manager(debug = True)  

# inquire for rank and proc from ESMF Virtual Machine
localPet = manager.local_pet
petCount = manager.pet_count

parallel = False
if petCount > 1:
    if petCount != 4:
      raise NameError('PET count must be 4 in parallel mode!')
parallel = True

# create two unique Grids
srcgrid = create_grid([18,10], True)
dstgrid = create_grid([10,8], False)

# srcgrid.write("srcgridSphPer")
# dstgrid.write("dstgridSphPer")

# create the Fields
srcfield = create_field(srcgrid, 'srcfield')
dstfield = create_field(dstgrid, 'dstfield')
exact_field = create_field(dstgrid, 'dstfield_exact')

# create the area fields
srcareafield = create_field(srcgrid, 'srcfracfield')
dstareafield = create_field(dstgrid, 'dstfracfield')

# create the fraction fields
srcfracfield = create_field(srcgrid, 'srcfracfield')
dstfracfield = create_field(dstgrid, 'dstfracfield')

# initialize the Fields to an analytic function
srcfield = build_analyticfield(srcfield, srcgrid)
exact_field = build_analyticfield(exact_field, dstgrid)

# run the ESMF regridding
dstfield, srcfracfield, dstfracfield = run_regridding(srcfield, dstfield, 
                                            srcfracfield, dstfracfield)

# compute the mass
srcmass = compute_mass(srcfield, srcareafield, srcfracfield, True)
dstmass = compute_mass(dstfield, dstareafield, 0, False)

# compare results and output PASS or FAIL
compare_fields(dstfield, exact_field, dstfracfield, srcmass, dstmass, parallel)

# closing remarks
print '\nThanks for using the Field regridding tutorial, have a nice day.\n'
