#!/bin/bash
# script to extract some variables from a netcdf file

# NO NEED FOR A MAIL ALLERT. IF CONVERT WORKS THEN THIS WILL AS WELL. IF CONVERT FAILS THEN IT WILL SEND A MESSAGE.

#########################################################################
tday=$1
firstday=`date --date="$2" +%Y-%m-%d` # convert YYYYMMDD to YYYY-MM-DD
dir0=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$tday/netcdf
#########################################################################

#extract start date/time of forecast
# and put it into output filename
odir=$(pwd)  # output file will be placed in current directory
cd $dir0
mkdir -p tmp

# LOG
log=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/logs/merge_log.txt
if [ -f "$log" ]
then
   rm $log
fi
touch $log

echo " "                                              >> $log
echo "Unpacking files (ncpdq -U)..."                  >> $log

flist=(*.nc)
f=${flist[0]}
# #get start of forecast
# start_date=${f:9:8}        # check date
# start_time=${f:18:2}:00:00 # HH:MM:SS

for f in TP4archv_*.nc
do
   cf=${f:9:8}
   if [ "$cf" -ge "$tday" ]
   then
      echo $f
      # unpack files to make sure that scale factors are same in each file
      ncpdq -U $f tmp/$f
   else
      echo "*.nc older than actual date"              >> $log
   fi
done

#combine unpacked files
echo " "                                              >> $log
echo "Combining unpacked files (ncrcat)..."           >> $log
ncrcat tmp/*.nc tmp.nc

#set name of output file
ofil=SWARPiceonly_forecast_start${tday}T000000Z.nc

# make final file by repacking tmp.nc
echo " "                                              >> $log
echo "Making $ofil (ncpdq)..."                        >> $log
ncpdq tmp.nc $ofil

rm -r tmp tmp.nc

#########################################################################
# rename dimensions
# ncrename -d dim00800,index1 $ofil
# ncrename -d dim00880,index2 $ofil
# ncrename -d record,time     $ofil

# change variable names to standard abbreviations in GRIB2 tables
ncrename -v fice,icec   $ofil #ice cover
ncrename -v hice,icetk  $ofil #ice thickness

# change time:units
ncatted -O -h -a units,time,m,c,"seconds since 1970-01-01T00:00:00Z"  $ofil

###########################################################################################
# most variable attributes set in hyc2proj:
# - just add shape of earth to "int stereographic"
ncatted -O -a "semi_major_axis",stereographic,c,f,6378273.     $ofil
ncatted -O -a "semi_minor_axis",stereographic,c,f,6378273.     $ofil
###########################################################################################

###########################################################################################
# Global attributes for  files
# (o or c = overwrite/create)
ncatted -O -h -a software_version,global,c,c,"NERSC-HYCOM (TOPAZ)"            $ofil
ncatted -O -h -a references,global,c,c,"www.nersc.no"                         $ofil
ncatted -O -h -a comment,global,c,c," "                                       $ofil
ncatted -O -h -a area,global,c,c,"TP4 (12.5km)"                               $ofil
ncatted -O -h -a field_type,global,c,c,"3-hourly"                             $ofil
ncatted -O -h -a forecast_start_date,global,c,c,"${firstday}T00:00:00Z"       $ofil
ncatted -O -h -a forecast_range,global,c,c,"5 day forecast"                   $ofil
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
ncatted -O -h -a title,global,o,c,"SWARP sea ice forecast"               $ofil # o=overwrite/create, c=format (also f=float)
# ncatted -O -h -a history,global,o,c,"NERSC-HYCOM output->hyc2proj->ncrcat"    $ofil

# Restart file date
reformat=0
if [ $reformat -eq 1 ]
then
   # reformat bulletin date
   # eg 2015-05-09T00:00:00Z -> 2015-05-09 UTC 00:00:00
   bdate_info=$(ncinfo $ofil | grep "bulletin_date")
   split=(${bdate_info//\ / })   # make array with delimiter space "\ "
   bdate_time=${split[1]}        # eg 2015-05-09T00:00:00Z
   split=(${bdate_time//T/ })    # make array with delimiter "T"
   bdate=${split[0]}             # eg 2015-05-09
   btime=${split[1]}             # eg 00:00:00Z
   btime=${btime:0:8}            # remove "Z"
   ncatted -O -h -a bulletin_date,global,o,c,"$bdate UTC $btime" $ofil   # add new attribute "restart_file_date"
fi

# Rename bulletin_date to restart_file_date (clearer)
ncrename -a bulletin_date,restart_file_date $ofil

# delete old attribute(s)
ncatted -a field_date,global,d,,            $ofil
ncatted -a history,global,d,,               $ofil
###########################################################################################

mkdir -p /work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$tday/final_output
mv $ofil ../final_output/
