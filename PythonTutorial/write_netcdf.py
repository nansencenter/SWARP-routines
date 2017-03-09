from netCDF4 import Dataset as ncopen
import numpy as np
from datetime import datetime

"""
fake data to be written
"""
Ntot     = 100
uice     = 1.+np.zeros(Ntot)
vice     = np.zeros(Ntot)
cno      = np.zeros(Ntot,dtype='int')
cno[50:] = 1
nc = ncopen('sample.nc','w', format='NETCDF4') #'w' stands for write

"""
The above line creates a netCDF file called "sample.nc" in the current folder. nc is a netCDF Dataset object that provides methods for storing data to the
file. nc also doubles as the root group. A netCDF group is basically a directory or folder within the netCDF dataset. This allows you to organize data as you would in a unix file system. Let's create a group for the heck of it:
"""

if 0:
   tempgrp = nc.createGroup('Temp_data')
else:
   tempgrp = nc

"""
Specifying dimensions

The next step is to specify the dimensions of the data. If you plan to save a multidimensional array of data, each dimension of that array needs to be given a name and a length:
"""

tempgrp.createDimension('point_index', Ntot)
tempgrp.createDimension('time', 1)

"""
Building variables

This step essentially pre-allocates NetCDF variables for data storage. NetCDF variables are very similar to numpy arrays. To construct them, you use the createVariable method:
"""
idx   = tempgrp.createVariable('point_index'          ,'i4' ,'point_index') # name, data type (i4 = 32-bit integer), dimension it depends on
Cno   = tempgrp.createVariable('contour_number'       ,'i4' ,'point_index') # name, data type (i4 = 32-bit integer), dimension it depends on
Uice  = tempgrp.createVariable('ice_velocity_east'    ,'f4' ,'point_index') # name, data type (f4 = 32-bit float  ), dimension it depends on
Vice  = tempgrp.createVariable('ice_velocity_north'   ,'f4' ,'point_index') # name, data type (f4 = 32-bit float  ), dimension it depends on
Time  = tempgrp.createVariable('time'                 ,'f4' ,'time')        # name, data type (f4 = 32-bit float  ), dimension it depends on

"""
Passing data into variables

Here, you simply pass your data into the variables you just created:
"""
Uice[:]  = uice #The "[:]" at the end of the variable instance is necessary
Vice[:]  = vice
idx[:]   = np.arange(Ntot,dtype='int')
Cno[:]   = cno

"""
Adding attributes

NetCDF attributes can be used to provide additional information about the dataset (i.e. metadata). You can add attributes to variables, groups and the dataset itself. This is optional but considered good practice:

#Add global attributes
f.description = "Example dataset containing one group"
f.history = "Created " + today.strftime("%d/%m/%y")

#Add local attributes to variable instances
longitude.units = 'degrees east'
latitude.units = 'degrees north'
time.units = 'days since Jan 01, 0001'
temp.units = 'Kelvin'
levels.units = 'meters'
temp.warning = 'This data is not real!'

You can add attributes any way you see fit, but you should be aware of the different attribute conventions that already exist. Most notable are the COARDS and Climate Forecast (CF) conventions. Even if you choose not to conform to any existing standard, I highly recommend creating a convenient and consistent naming system for yourself. For example, I ensure that all my variables have units and long_name attributes.
"""

#get time in days since Jan 01,01
today = datetime.today()
time_num = today.toordinal()
Time[0] = time_num

nc.close()
