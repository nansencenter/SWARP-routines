
# The function 'ncVarSizes' returns the sizes of the dimensions
#  of a variable from a netCDF file, ordered "Fortran style"

# Input:
#  fileID  - ID for netCDF file,
#            returned from a call to e.g. open.nc
#  varName - name of variable of interest

# Example:
# --------
# Assuming
#  sshID    <- open.nc("ssh.nc", write=FALSE)
#  varSizes <- ncVarSizes(sshID,"ssh")

# ...and the netCDF file contains
# ncdump -h ssh.nc 
#  netcdf ssh {
#  dimensions:
#	x = 800 ;
#	y = 880 ;
#	time = UNLIMITED ; // (6 currently)
#  :
#	short ssh(time, y, x) ;

# ...this function will return the vector
#  [1] 800 880 6


ncVarSizes <- function(fileID, varName) {


  varInq <- var.inq.nc(fileID,varName)

  varSizes <- c(varInq$ndims)

  for (dimNo in 1:varInq$ndims) {
#    sz <- varInq$dimids
#    print(sz)
#    print(sz[])
    dimInq          <- dim.inq.nc(fileID,varInq$dimids[dimNo])
    varSizes[dimNo] <- dimInq$length
  }

  return(varSizes)
}