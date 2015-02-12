# script to extract some variables from a netcdf file

#########################################################################
if [ $# -ne 1 ]; then
   echo "Usage:"
   echo "merge_TP4archv.sh [directory name]"
   echo "- merges netcdf files in directory into one file with ncrcat"
   exit
fi
dir0=$1
#########################################################################

# TODO extract start date/time of forecast
# and put it into output filename
module load nco
odir=`pwd`  # output file will be placed in current directory
cd $dir0

for f in *.nc
do
   echo $f
   #get start of forecast
   start_date=${f:18:8}
   start_time=${f:27:6}
   break
done

#set name of output file
ofil=SWARPwavesice_forecast_start${start_date}T${start_time}Z.nc
echo "Making $ofil"
ncrcat *.nc $ofil

#########################################################################
# rename dimensions
# ncrename -d dim00800,index1 $ofil
# ncrename -d dim00880,index2 $ofil
# ncrename -d record,time     $ofil

# change variable names to standard abbreviations in GRIB2 tables
ncrename -v fice,icec   $ofil #ice cover
ncrename -v hice,icetk  $ofil #ice thickness

# ###########################################################################################
# # now add variable attributes:
# # icec:
# ncatted -O -a long_name,icec,c,c,"Sea ice concentration" $ofil
# ncatted -O -a units,icec,c,c," "                         $ofil
# ncatted -O -a add_offset,icec,c,f,0.                     $ofil
# ncatted -O -a scale_factor,icec,c,f,1.                   $ofil
# ncatted -O -a _FillValue,icec,c,f,-32767.                $ofil

# # icetk:
# ncatted -O -a long_name,icetk,c,c,"Sea ice thickness"     $ofil
# ncatted -O -a units,icetk,c,c,"m"                         $ofil
# ncatted -O -a add_offset,icetk,c,f,0.                     $ofil
# ncatted -O -a scale_factor,icetk,c,f,1.                   $ofil
# ncatted -O -a _FillValue,icetk,c,f,-32767.                $ofil

# # dmax:
# ncatted -O -a long_name,dmax,c,c,"Sea ice maximum floe size"   $ofil
# ncatted -O -a units,dmax,c,c,"m"                               $ofil
# ncatted -O -a add_offset,dmax,c,f,0.                           $ofil
# ncatted -O -a scale_factor,dmax,c,f,1.                         $ofil
# ncatted -O -a _FillValue,dmax,c,f,-32767.                      $ofil

# # Hs:
# ncatted -O -a long_name,Hs,c,c,"Significant wave height" $ofil
# ncatted -O -a units,Hs,c,c,"m"                           $ofil
# ncatted -O -a add_offset,Hs,c,f,0.                       $ofil
# ncatted -O -a scale_factor,Hs,c,f,1.                     $ofil
# ncatted -O -a _FillValue,Hs,c,f,-32767.                  $ofil

# # Tw:
# ncatted -O -a long_name,Tw,c,c,"Mean wave period (T_m_02)"  $ofil
# ncatted -O -a units,Tw,c,c,"s"                              $ofil
# ncatted -O -a add_offset,Tw,c,f,0.                          $ofil
# ncatted -O -a scale_factor,Tw,c,f,1.                        $ofil
# ncatted -O -a _FillValue,Tw,c,f,-32767.                     $ofil
# ###########################################################################################

# ncrename -v icef,sea_ice_concentration    $ofil
# ncrename -v iceh,sea_ice_thickness        $ofil
# ncrename -v dmax,sea_ice_max_floe_size    $ofil
# ncrename -v Hs,significant_wave_height    $ofil
# ncrename -v Tw,mean_wave_period           $ofil

###########################################################################################
# Global attributes for  files
# (o or c = overwrite/create)
ncatted -O -h -a software_version,global,c,c,"NERSC-HYCOM (TOPAZ)"            $ofil
ncatted -O -h -a references,global,c,c,"www.nersc.no"                         $ofil
ncatted -O -h -a comment,global,c,c," "                                       $ofil
ncatted -O -h -a area,global,c,c,"TP4 (12.5km)"                               $ofil
ncatted -O -h -a field_type,global,c,c,"3-hourly"                             $ofil
ncatted -O -h -a forecast_range,global,c,c,"3 day forecast"                   $ofil
# ncatted -O -h -a forecast_type,global,c,c,"forecast"                        $ofil
ncatted -O -h -a institution,global,c,c,"NERSC"                               $ofil
ncatted -O -h -a institution_references,global,c,c,"http://www.nersc.no/"     $ofil
ncatted -O -h -a data_centre,global,c,c,"NERSC"                               $ofil
ncatted -O -h -a data_centre_references,global,c,c,"www.nersc.no"             $ofil
ncatted -O -h -a contact,global,c,c,"timothy.williams@nersc.no"               $ofil
ncatted -O -h -a project,global,c,c,"SWARP"                                   $ofil
ncatted -O -h -a project_references,global,c,c,"swarp.nersc.no"               $ofil
ncatted -O -h -a distribution_statement,global,c,c,"No restrictions"          $ofil
ncatted -O -h -a operational_status,global,c,c,"test"                         $ofil
#
ncatted -O -h -a title,global,o,c,"SWARP Waves-in-ice forecast"               $ofil
# ncatted -O -h -a history,global,o,c,"NERSC-HYCOM output->hyc2proj->ncrcat"    $ofil

# delete old attributes
ncatted -a field_date,global,d,,                                              $ofil

###########################################################################################

mv $ofil $odir
