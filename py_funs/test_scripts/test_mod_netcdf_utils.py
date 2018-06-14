import numpy as np
import mod_netcdf_utils as mnu

tdep_ll = False
obsdir = os.getenv('INPUT_OBS_DATA_DIR')
if 1:
    f = obsdir+'/AMSR2_ice_conc/Arc_20180520_res3.125_pyres.nc'
    vname = 'fice'
elif 1:
    f = obsdir+'/OSISAF_ice_drift/2018/06/ice_drift_nh_polstere-625_multi-oi_201806031200-201806051200.nc'
    vname = 'lat1'
    tdep_ll = True

nci = mnu.nc_getinfo(f, timedep_lonlat=tdep_ll)
lon, lat = nci.get_lonlat()
V = nci.get_var(vname).values

I, J = np.where(~V.mask)
i0 = I[0]
j0 = J[0]
ij_range = i0, i0+2, j0, j0+2
i0, i1, j0, j1 = ij_range

lon1 = lon[i0:i1, j0:j1]
lat1 = lat[i0:i1, j0:j1]
V1 = V[i0:i1, j0:j1]
print(V1)

lon2, lat2 = nci.get_lonlat(ij_range=ij_range)
V2 = nci.get_var(vname, ij_range=ij_range).values
print(V2)

print('test reduced lon: %s' %(str((lon2==lon1).all())) )
print('test reduced lat: %s' %(str((lat2==lat1).all())) )
print('test reduced var: %s' %(str((V2.data==V1.data).all())) )
