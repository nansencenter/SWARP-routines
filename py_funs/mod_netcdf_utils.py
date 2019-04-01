import numpy as np
import os,sys

from datetime import datetime,timedelta
from netCDF4 import Dataset as ncopen
try:
    from netCDF4 import netcdftime as NCT
except:
    import netcdftime as NCT

from osgeo import gdal, ogr, osr
import mod_reading as MR

def lonlat_names(ncfil):
    nc = ncopen(ncfil)
    for vbl in nc.variables:
        if 'lon' in vbl[:3].lower():
            lonname = vbl
        if 'lat' in vbl[:3].lower():
            latname = vbl
    nc.close()
    return lonname, latname


def get_time_name(nc):
    """
    NEMO outputs call time "time_counter"
    CS2-SMOS thickness files use 'tc' for time dimension, but
    'time_bnds' for time variable
    """
    time_name = None
    for tname in [
            'time',
            'time_counter',
            'time_bnds',
            ]:
        if tname in nc.dimensions:
            time_name = tname
            break

    return time_name


def get_time_converter(time):
    """
    reftime,time_converter = get_time_converter(time)
    *input:
    time = nc.variables['time'],
    where nc is a netCDF4.Dataset object
    *outputs:
    reftime - datetime.datetime object
    time_converter netCDF4.netcdftime.utime object
    """

    tu        = time.units
    time_info = tu.split()

    # determine time format of reference point:
    Unit = time_info[0]# 1st word gives units
    Unit = Unit.strip('s')
    if Unit=='econd':
        Unit = 'second'

    time_info.remove(time_info[0])
    time_info.remove(time_info[0]) # 'since'

    # rest is date and time of reference point
    if len(time_info)>=2:
        cdate, ctime = time_info[0:2]
    else:
        if ('T' in time_info[0]) and ('Z' in time_info[0]):
            # cdate+T...Z format for time
            split1 = tu.split('T')
            ctime  = split1[1].strip('Z')
            cdate  = split1[0].split()[2]

    # reformat cdate to YYYYMMDD
    if '-' in cdate:
        # remove '-'
        # - otherwise assume YYYYMMDD format
        split2 = cdate.split('-')
        for loop_i in range(1,3):
            if len(split2[loop_i])==1:
                split2[loop_i] = '0'+split2[loop_i]
        cdate = split2[0]+split2[1]+split2[2] # should be YYYYMMDD now
    if len(cdate)<8:
        cdate = (8-len(cdate))*'0'+cdate

    # reformat ctime to HHMMSS
    if ':' in ctime:
        # remove ':'
        # - otherwise assume HHMMSS format
        split2 = ctime.split(':')
        for loop_i in range(0,3):
            if (split2[loop_i])==1:
                split2[loop_i] = '0'+split2[loop_i]
        ctime = split2[0]+split2[1]+split2[2] # should be HHMMSS now


    # now can make new string where format is known
    # - this is to pass into netcdftime.utime
    # - NB can't always use strftime/strptime since it only works after 1900
    cyear0 = cdate[:4]
    cmon0  = cdate[4:6]
    cday0  = cdate[6:8]
    chr0   = ctime[:2]
    cmin0  = ctime[2:4]
    csec0  = ctime[4:]
    #
    year0   = int(cyear0)
    mon0    = int(cmon0)
    day0    = int(cday0)
    hr0     = int(chr0)
    min0    = int(cmin0)
    sec0    = int(float(csec0))
    reftime = datetime(year0,mon0,day0,hr0,min0,sec0)

    init_string = (Unit+'s since '+
            cyear0+'-'+cmon0+'-'+cday0+' '+
            chr0+'-'+cmin0+'-'+csec0[:2])

    if 'calendar' in time.ncattrs():
        time_converter = NCT.utime(init_string,calendar=time.calendar)
    else:
        time_converter = NCT.utime(init_string)

    return time_converter


def nc_get_var(ncfil, vblname, time_index=None,
        depth_index=0, ij_range=None, **kwargs):
    """
    vbl=nc_get_var(ncfil, vblname, time_index=None)
    *ncfil is string (filename)
    *vname is string (variable name)
    *time_index is record number to get
    *depth_index is horizon number to get
    *vbl is a mod_reading.var_object instance
    """
    # NB kwargs is not used, but is there as a dummy to avoid having to sort kwargs
    # before calling this function

    with ncopen(ncfil) as nc:
        vbl0 = nc.variables[vblname]

        # get the netcdf attributes
        attlist = vbl0.ncattrs()
        attvals = []
        for att in attlist:
            attval = getattr(vbl0, att)
            attvals.append(attval)

        dims  = vbl0.dimensions
        shape = vbl0.shape

        # do we want to limit the range
        if ij_range is not None:
            i0, i1, j0, j1 = ij_range

        # some attributes that depend on rank
        if vbl0.ndim==1:
            vals = vbl0[:]

        elif vbl0.ndim==2:
            if ij_range is not None:
                vals = vbl0[i0:i1, j0:j1]
            else:
                vals = vbl0[:, :]

        elif vbl0.ndim==3:
            if time_index is None:
                if shape[0]==1:
                    time_index = 0
            if time_index is None:
                if ij_range is not None:
                    vals = vbl0[:, i0:i1, j0:j1]
                else:
                    vals = vbl0[:, :, :]
            else:
                if ij_range is not None:
                    vals = vbl0[time_index, i0:i1, j0:j1]
                else:
                    vals = vbl0[time_index, :, :]
                dims = dims[1:]

        elif vbl0.ndim==4:
            if time_index is None:
                if shape[0]==1:
                    time_index = 0

            if time_index is None:
                if ij_range is not None:
                    vals = vbl0[:, depth_index, i0:i1, j0:j1]
                else:
                    vals = vbl0[:, depth_index, :, :]
                dims = (dims[0], dims[2], dims[3])
            else:
                if ij_range is not None:
                    vals = vbl0[time_index, depth_index, i0:i1, j0:j1]
                else:
                    vals = vbl0[time_index, depth_index, :, :]
                dims = dims[2:]

    attlist.append('dimensions')
    attvals.append(dims)
    return MR.var_object(vals, extra_atts=[attlist, attvals])


def nc_get_time(ncfil, time_name='time'):
    """
    vbl=nc_get_time(ncfil, time_name)
    *ncfil is string (filename)
    *vname is string (variable name)
    *time_index is record number to get
    *vbl is a mod_reading.var_object instance
    """
    with ncopen(ncfil) as nc:
        time = nc.variables[time_name][:]
    return time


def nc_get_var_atts(ncfil,vblname):
    """
    vbl=nc_get_var_atts(ncfil, vblname, time_index=None)
    Parameters:
    *ncfil is string (filename)
    *vname is string (variable name)

    Returns:
    *vbl is a mod_reading.var_object instance
    """
    nc   = ncopen(ncfil)
    vbl0 = nc.variables[vblname]

    # get the netcdf attributes
    attlist = vbl0.ncattrs()
    atts    = {}
    for att in attlist:
        attval = getattr(vbl0,att)
        atts.update({att:attval})

    nc.close()
    return atts


def nc_get_dim(ncfil,vblname):
    """
    vbl=nc_get_var(ncfil,vblname,time_index=None)
    *ncfil is string (filename)
    *vname is string (variable name)
    *vbl    is a mod_reading.var_object instance
    """

    nc = ncopen(ncfil)

    if vblname in nc.variables:
        vbl0 = nc.variables[vblname]

        # get the netcdf attributes
        attlist = vbl0.ncattrs()
        attvals = []
        for att in attlist:
            attval = getattr(vbl0,att)
            attvals.append(attval)
        Xatts = [attlist,attvals]

        vals = vbl0[:]
        nc.close()

        return MR.var_object(vals, extra_atts=Xatts)
    else:
        raise ValueError(vbl_name+' not given as an array')


class nc_getinfo:
    """
    mod_netcdf_utils.nc_getinfo(ncfil, lonlat_file=None)
    """

    def __init__(self,
              ncfil,
              time_index=None,
              lonlat_file=None,
              timedep_lonlat=False):
        # NB time_index not needed here, but adding it as a dummy index
        # means this object can be used along with other similar objects
        # in some functions in mod_reading.py

        self.filename  = ncfil
        if ncfil[0]=='/':
            bn           = os.path.basename(ncfil)
            self.basedir = ncfil.strip(bn)
        else:
            bn           = ncfil
            self.basedir = os.getcwd()+'/'

        self.basename    = os.path.splitext(bn)[0]
        self.filetype    = 'netcdf'
        self.object_type = 'netcdf'

        # things to work with plotting stuff in mod_reading.py
        self.HYCOM_region     = None
        self.get_fixed_lonlat = self.get_lonlat

        # for drifters
        self.timedep_lonlat = timedep_lonlat

        # added here manually
        # - TODO could possibly be determined
        #    from netcdf metadata though
        # - could also be an input
        self.reftime_sig = 'start of forecast'


        # open the file
        nc     = ncopen(ncfil)
        dkeys = list(nc.dimensions.keys())
        vkeys = list(nc.variables.keys())

        # remove dimensions from variables
        self.dimensions = dkeys
        for key in dkeys:
            if key in vkeys:
                vkeys.remove(key)
        Nkeys = len(vkeys)

        # time info:
        self._set_time_info(nc)

        # get global netcdf attributes
        class ncatts:
            def __init__(self,nc):
                for att in nc.ncattrs():
                    attval    = getattr(nc,att)
                    setattr(self,att,attval)
                return
            def atts2list(self):
                return vars(self).keys()
        self.ncattrs = ncatts(nc)



        # grid info:
        if lonlat_file is None:
            lonlat_file = ncfil
        self.lonlat_file = lonlat_file

        # are lon,lat dimensions?
        self.lonname, self.latname = lonlat_names(self.lonlat_file)
        self.lonlat_dim            = (self.lonname in self.dimensions)

        # basic lon-lat info
        with ncopen(self.lonlat_file) as nc2:
            lon = nc2.variables[self.lonname]
            lat = nc2.variables[self.latname]
            if self.lonlat_dim:
                self.lon0 = lon[0]
                self.lat0 = lat[0]
                self.lonlat_corners = np.array([
                        (lon[0], lat[0]),
                        (lon[-1], lat[0]),
                        (lon[-1], lat[-1]),
                        (lon[0], lat[-1]),
                        (lon[0], lat[0]),
                        ]).T

                # get example variable:
                # - for eg plotting, need to make the lon/lat matrices
                # (converted from vectors)
                # have the same shape as the variables
                for vkey in vkeys:
                    if self.lonname in nc.variables[vkey].dimensions:
                        ncvar = nc.variables[vkey]
                        break
                for dkey in ncvar.dimensions:
                    if dkey==self.lonname:
                        self.lon_first = True
                        self.shape     = (len(lon), len(lat))
                        break
                    elif dkey==self.latname:
                        self.lon_first = False
                        self.shape     = (len(lat), len(lon))
                        break
            elif self.timedep_lonlat:
                self.lon0  = None
                self.lat0  = None
                self.lonlat_corners = None
                self.shape = lon[0,:,:].shape
            else:
                self.lon0  = lon[0, 0]
                self.lat0  = lat[0, 0]
                self.lonlat_corners = np.array([
                        (lon[0, 0], lat[0, 0]),
                        (lon[-1, 0], lat[-1, 0]),
                        (lon[-1, -1], lat[-1, -1]),
                        (lon[0, -1], lat[0, -1]),
                        (lon[0, 0], lat[0, 0]),
                        ]).T
                self.lon1  = lon[-1, -1]
                self.lat1  = lat[-1, -1]
                self.shape = lon.shape

        ny, nx      = self.shape
        self.Npts_x = nx    # No of points in x dirn
        self.Npts_y = ny    # No of points in y dirn
        self.Npts   = nx*ny # Total no of points


        # projection info:
        proj_list = [
                'stereographic',
                'projection_3',
                'Polar_Stereographic_Grid',
                'Lambert_Azimuthal_Grid',
                ]

        # could also have mercator or regular lon-lat
        HAVE_PROJ    = 0    # if 0 assume HYCOM native grid
        keys_to_remove = []
        for proj_name in proj_list:
            if proj_name in vkeys:
                proj      = nc.variables[proj_name]
                att_list  = proj.ncattrs()
                HAVE_PROJ = 1
                keys_to_remove.append(proj_name) # don't want to keep projection variable with the other variables
                break

        if HAVE_PROJ:
            # object with the netcdf attributes of projection variable
            # + some extra proj-dependent info
            att_list_full = [att_list[i] for i in range(len(att_list))]
            att_vals_full = []
            for att in att_list:
                att_val = proj.getncattr(att)
                att_vals_full.append(att_val)

            # specific to stereographic
            if proj.getncattr('grid_mapping_name') in [
                    'lambert_azimuthal_equal_area',
                    'polar_stereographic',
                    ]:
                # add x, y resolution to ncinfo.proj_info
                att_list_full.extend(['x_resolution', 'y_resolution'])

                if 'x' in nc.dimensions:
                    xname, yname = 'x', 'y'
                elif 'xc' in nc.dimensions:
                    xname, yname = 'xc', 'yc'
                xx = nc.variables[xname][0:2]
                yy = nc.variables[yname][0:2]
                dx = xx[1]-xx[0]
                dy = yy[1]-yy[0]

                #convert to m
                xunits = nc.variables[xname].units.split()
                fac    = 1.
                if len(xunits)==2:
                    fac = float(xunits[0])
                    xunits.remove(xunits[0])

                if xunits[0]=='km':
                    fac = fac*1.e3
                #
                att_vals_full.extend([dx*fac, dy*fac])

            self.proj_info = MR.proj_obj(att_list_full, att_vals_full)
        else:
            self.proj_info = None

        # variable list
        # - remove some other variables from vkeys
        # - eg projection,lon,lat
        # - time_bnds
        # - TODO model_depth?
        keys_to_remove.append("time_bnds")
        # keys_to_remove.append('model_depth')
        if not self.timedep_lonlat:
            keys_to_remove.extend([self.lonname, self.latname])

        # other variables to remove
        for key in keys_to_remove:
            if key in vkeys:
                vkeys.remove(key)

        self.variables     = vkeys
        self.variables3d   = None  #TODO enable treatment of 3d fields
        self.all_variables = vkeys

        nc.close()
        return


    def _set_time_info(self, nc):
        """
        * sets self.time_name  = name of time variable
        * sets self.time_dim = True or False - is time is a dimension
        * sets self.time_converter = function to convert time value to datetime
        * sets datetimes
        """

        self.time_name = get_time_name(nc)
        self.time_dim  = (self.time_name is not None)

        if not self.time_dim:
            self.datetimes = None
            return

        time = nc.variables[self.time_name]
        fmt  = '%Y-%m-%d %H:%M:%S'

        self.time_converter = get_time_converter(time)

        arr             = time[:] #time values
        self.datetimes  = []
        self.timevalues = []

        Unit = self.time_converter.units.lower()

        for i, tval in enumerate(arr):
            if isinstance(tval, np.int32):
                # can be problems if int32 format
                tval  = int(tval)
            try:
                cdate = self.time_converter.num2date(tval).strftime(fmt)
            except ValueError:
                # might get errors if close to end/start of month
                # eg CS2-SMOS
                tval=round(float(tval))
                cdate = self.time_converter.num2date(tval).strftime(fmt)
            dto    = datetime.strptime(cdate, fmt)         # now a proper datetime object
            self.datetimes.append(dto)

            if i==0:
                self.reftime  = dto

            tdiff = (dto-self.reftime).total_seconds()
            if Unit=='seconds':
                self.timevalues.append(tdiff/3600.)         # convert to hours for readability
                self.timeunits = 'hour'
            elif Unit=='hours':
                self.timevalues.append(tdiff/3600.)         # keep as hours
                self.timeunits = 'hour'
            elif Unit=='days':
                self.timevalues.append(tdiff/3600./24.)    # keep as days
                self.timeunits = 'day'

        self.number_of_time_records    = len(self.datetimes)


    def nearestDate(self, pivot):
        """
        dto,time_index = self.nearestDate(dto0)
        dto0  = datetime.datetime object
        dto    = datetime.datetime object - nearest value in self.datetimes to dto0
        time_index: dto=self.datetimes[time_index]
        """
        dto        = min(self.datetimes, key=lambda x: abs(x - pivot))
        time_index = self.datetimes.index(dto)
        return dto, time_index


    def timeval_to_datetime(self, timeval):

        # check format of time
        i32 = np.array([0],dtype='int32')
        if type(i32[0])==type(timeval):
            timeval  = int(timeval)

        if self.timeunits=='second':
            dt = self.reftime +timedelta(seconds=timeval)
        elif self.timeunits=='hour':
            dt = self.reftime +timedelta(hours=timeval)
        elif self.timeunits=='day':
            dt = self.reftime +timedelta(timeval) #NB works for fraction of days also
        return dt


    def get_lonlat(self, vec2mat=True, **kwargs):
        """
        Parameters:
        * vec2mat (bool): if lon, lat are vectors (eg if they are dimensions),
            convert to matrices
        * ij_range: list of integers to reduce the size of arrays
            - can be: None or [i0, i1, j0, j1]
            - let lon_all, lat_all = self.get_lonlat(vec2mat=True)
            - let lon, lat = self.get_lonlat(vec2mat=True, ij_range=[i0, i1, j0, j1])
            - lon = lon_all[i0:i1, j0:j1], lat = lat_all[i0:i1, j0:j1]

        Returns:
        * lon, lat (np.array)
        """

        if self.timedep_lonlat:
            # return lon,lat using get_var
            lon = self.get_var(self.lonname, **kwargs)
            lat = self.get_var(self.latname, **kwargs)
            return (lon.values.filled(np.nan),
                      lat.values.filled(np.nan))

        ij_range = kwargs.get('ij_range', None)

        with ncopen(self.lonlat_file) as nc:
            lono = nc.variables[self.lonname]
            lato = nc.variables[self.latname]

            if ij_range is None:
                if lono.ndim==2:
                    lon = lono[:, :]
                    lat = lato[:, :]
                else:
                    lon = lono[:]
                    lat = lato[:]
                    if vec2mat:
                        if self.lon_first:
                            # lon in cols, lat in rows
                            lon, lat  = np.meshgrid(lon, lat, indexing='ij')
                        else:
                            # lon in rows, lat in cols
                            lon, lat  = np.meshgrid(lon, lat, indexing='xy')
            else:
                print('ij_range: ', ij_range)
                i0, i1, j0, j1 = ij_range
                if lono.ndim==2:
                    lon = lono[i0:i1, j0:j1]
                    lat = lato[i0:i1, j0:j1]
                else:
                    if self.lon_first:
                        # lon in cols, lat in rows
                        lon = lono[j0:j1]
                        lat = lato[i0:i1]
                    else:
                        # lon in rows, lat in cols
                        lon = lono[i0:i1]
                        lat = lato[j0:j1]
                    if vec2mat:
                        if self.lon_first:
                            # lon in cols, lat in rows
                            lon, lat  = np.meshgrid(lon, lat, indexing='ij')
                        else:
                            # lon in rows, lat in cols
                            lon, lat  = np.meshgrid(lon, lat, indexing='xy')

        return lon, lat

    def get_bbox(self, mapping, **kwargs):
        """
        Parameters:
        * mapping: pyproj mapping
        * kwargs for self.get_lonlat()

        Returns:
        * bbox = [xmin, xmax, ymin, ymax], where x,y are coordinates specified by mapping
        """

        lon, lat = self.get_lonlat(**kwargs)
        x, y = mapping(lon, lat)
        return [x.min(), x.max(), y.min(), y.max()]


    def get_var(self, vname, **kwargs):
        """
        Call: self.get_var(vname, time_index=None)
        Inputs:
        vname = string (name of variable)
        time_index = integer: if time_index
        Returns: mod_reading.var_object instance
        """
        vname = MR.check_names(vname, self.variables)
        return nc_get_var(self.filename, vname, **kwargs)


    def get_var_atts(self, vname):
        """
        Call: self.get_var(vname)
        Inputs:
        vname = string (name of variable)
        Returns: dictionary of attributes
        """

        if 'time' in vname:
            vname = self.time_name
        elif 'lon' in vname:
            vname = self.lonname
        elif 'lat' in vname:
            vname = self.latname
        else:
            vname = MR.check_names(vname, self.variables)

        return nc_get_var_atts(self.filename, vname)


    def get_global_atts(self):
        """
        Call: self.get_global_atts()
        Returns: dictionary of attributes
        """

        return vars(self.ncattrs)


    def get_time(self):
        """
        vbl=nc_get_var(ncfil, vblname, time_index=None)
        *ncfil is string (filename)
        *vname is string (variable name)
        *time_index is record number to get
        *vbl is a mod_reading.var_object instance
        """
        return nc_get_time(self.filename, time_name=self.time_name)


    def get_dim(self, dname):
        """
        call: self.get_dim(dname)
        inputs:
        dname = string (name of dimension)
        returns: mod_reading.var_object instance
        """

        vname = MR.check_names(dname, self.dimensions)
        vbl   = nc_get_dim(self.filename, dname)

        return vbl


    def imshow(self, var_opts, **kwargs):
        """
        pobj    = self.imshow(var_opts, time_index=0, pobj=None,
              clim=None, add_cbar=True, clabel=None, show=True,
              test_ijs=None)
        """
        return MR.imshow(self, var_opts, **kwargs)


    def plot_var(self, var_opts, **kwargs):
        """
        pobj,bmap = self.plot_var(var_opts, time_index=0,
            pobj=None, bmap=None, HYCOMreg='TP4',
            clim=None, add_cbar=True, clabel=None, show=True,
            test_lonlats=None, date_label=0):
        """
        return MR.plot_var(self,var_opts,**kwargs)


    def plot_var_pair(self, var_opts1, var_opts2, pobj=None, bmap=None, **kwargs):
        """
        pobj, bmap=self.plot_var_pair(var_opts1, var_opts2, pobj=None, bmap=None, **kwargs)
        """
        return MR.plot_var_pair(self,var_opts1,var_opts2,**kwargs)


    def make_png(self, var_opts, **kwargs):
        """
        pobj, bmap=self.make_png(var_opts, pobj=None, bmap=None, figdir='.', time_index=0, date_label=2, **kwargs)
        """
        return MR.make_png(self, var_opts, **kwargs)


    def make_png_pair(self, var_opts1, var_opts2, **kwargs):
        """
        pobj,bmap = self.make_png_pair(var_opts1, var_opts2,
            pobj=None, bmap=None, figdir='.', date_label=2, **kwargs)
        """
        return MR.make_png_pair(self,var_opts1,var_opts2,**kwargs)


    def compare_ice_edge_obs(self, **kwargs):
        """
        pobj,bmap,obsfil = self.compare_ice_edge_obs(pobj=None, bmap=None, time_index=0,
            obs_type='OSISAF', date_label=1, figname=None, **kwargs)
        """
        return MR.compare_ice_edge_obs(self, **kwargs)


    def interp2points(self, varname, target_lonlats, **kwargs):
        """
        self.interp2points(varname, target_lonlats, **kwargs)
        * fobj is a file object eg from netcdf
        * varname is a string (name of variable in file object)
        * target_lonlats = [target_lon, target_lat], target_lon/lat are numpy arrays
        * kwargs for mod_reading.interp2points()
        """
        return MR.interp2points(self, varname, target_lonlats, **kwargs)



    def get_external_data(self, ncfil, vname, dto_in=None, time_index=None,
		lonlat_file=None, mapping=None):
         target_lonlats = self.get_lonlat()

         # initialise mnu.nc_getinfo object
         nci = mnu.nc_getinfo(ncfil, lonlat_file=lonlat_file)

         tind = time_index
         if time_index is None:
             tind  = 0
             if type(dto_in) is not type(None):
                 dto,tind = nci.nearestDate(dto_in)
         nci = nc_getinfo(ncfil, lonlat_file=lonlat_file)

         if time_index is not None:
             tind = time_index
         elif type(dto_in) is not type(None):
             dto, tind = self.nearestDate(dto_in)
         else:
             tind = 0

         vout = nci.interp2points(vname, target_lonlats,
                 time_index=tind,
                 mapping=mapping)
         return vout


    def get_area_euclidean(self, pyproj_map, **kwargs):
        """
        Calculates element area from netcdf file
        Assumes regular grid in the given projection
        area  = self.get_area_euclidean(pyproj_map, **kwargs)

        Parameters:
        -----------
        pyproj_map : pyproj.Proj
        kwargs for self.get_lonlat

        Returns:
        --------
        * area (float)
        """
        lon, lat = self.get_lonlat(**kwargs)
        pp = pyproj_map
        x, y = pp(lon, lat)
        dy = np.max([np.abs(np.mean(y[:, 2]-y[:, 1])), np.abs(np.mean(y[1, :]-y[0, :]))])
        dx = np.max([np.abs(np.mean(x[:, 2]-x[:, 1])), np.abs(np.mean(x[1, :]-x[0, :]))])
        area = np.abs(dx*dy)
        return area



    def MIZmap(self, **kwargs):
        """
        Call  : self.MIZmap(var_name='dmax',do_sort=False,EastOnly=True,plotting=True,**kwargs)
        Inputs:
            var_name is variable to find MIZ from
            **kwargs to be passed onto MIZchar.get_MIZ_poly:
                outdir='.',do_sort=True
        Returns: MIZchar.MIZpoly object
        """
        return MR.MIZmap(self, **kwargs)


    def areas_of_disagreement(self, **kwargs):
        """
        MPdict,tfiles,Pdict = self.areas_of_disagreement(obs_type='OSISAF',
                time_index=0,do_sort=True,EastOnly=True,
                plotting=True,HYCOMreg='Arctic',**kwargs)
        kwargs: outdir='.',do_sort=True
        """
        return MR.areas_of_disagreement(self, **kwargs)


    def make_png_all(self, var_opts, **kwargs):
        """
        self.make_png_all(var_opts, HYCOMreg=None, figdir='.')
        """
        MR.make_png_all(self, var_opts, **kwargs)
        return


    def make_png_pair_all(self, var_opts1, var_opts2, **kwargs):
        """
        self.make_png_pair_all(var_opts1,var_opts2,HYCOMreg=None,figdir='.')
        """
        MR.make_png_pair_all(self, var_opts1, var_opts2, **kwargs)
        return

    def compare_ice_edge_obs_all(self, **kwargs):
        return MR.compare_ice_edge_obs(self, **kwargs)

def get_amsr2_gdal_dataset(filename):
    ''' Return geocoded GDAL Dataset matching AMSR2 Arc_*_res3.125_pyres.nc

    Parameters of the projection (min/max of x/y and proj4 string) are hardcoded
    based on experiments using Nansat.

    Parameters
    ----------
        filename : str
            input file name
    Returns
    -------
        dst_ds : GDALDataset
            destination dataset in memory
    '''

    # check that file is correct
    ds = gdal.Open(filename)
    title = ds.GetMetadata()['NC_GLOBAL#title']
    if (not 'Daily averaged Arctic sea ice concentration derived from AMSR2' in
         title):
        raise Exception('Not correct inpu file %s' % filename)

    # hardcode resolution and min/max of X/Y coordinates in meters
    grid_resolution = 3125
    min_x = -3800000
    max_x = 3800000
    off_x = grid_resolution * -16
    min_x = min_x + off_x
    max_x = max_x + off_x

    min_y = -5600000
    max_y = 5600000
    off_y = grid_resolution * 80
    min_y = min_y + off_y
    max_y = max_y + off_y

    # hardcode projection
    srs_proj4 = '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lat_ts=70 +lon_0=-45 +no_defs'
    srs = osr.SpatialReference()
    srs.ImportFromProj4(str(srs_proj4))
    srs_wkt = srs.ExportToWkt()

    # create dataset
    subds0 = gdal.Open(ds.GetSubDatasets()[0][0])
    dst_ds = gdal.GetDriverByName('MEM').Create('tmp', subds0.RasterXSize,
                                                                        subds0.RasterYSize,
                                                                        1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform((min_x, grid_resolution, 0, min_y, 0, grid_resolution))
    dst_ds.SetProjection(srs_wkt)

    # set no_data_value for the band
    band = dst_ds.GetRasterBand(1)
    NoData_value = -999999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()

    return dst_ds
