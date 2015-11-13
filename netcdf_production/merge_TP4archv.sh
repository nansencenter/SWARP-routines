#!/bin/bash
# script to extract some variables from a netcdf file

# NO NEED FOR A MAIL ALERT. IF CONVERT WORKS THEN THIS WILL AS WELL. IF CONVERT FAILS THEN IT WILL SEND A MESSAGE.

source $SWARP_ROUTINES/source_files/hex_vars.src
email=`cat $SWARP_ROUTINES/forecast_scripts/fc_alert_email.txt`

#########################################################################
# inputs are start date of forecast
if [ $# -eq 0 ]
then
   # do a test
   tday=`date +%Y%m%d`                       # start date of forecast YYYYMMDD
   echo \${tday}=$tday
   tday_long=`date --date=$tday +%Y-%m-%d`   # start date of forecast YYYY-MM-DD
   lognm=merge_log_test.txt
else
   tday=$1                                   # start date of forecast YYYYMMDD
   tday_long=`date --date=$tday +%Y-%m-%d`   # start date of forecast YYYY-MM-DD
   # tday_long=`date --date="$2" +%Y-%m-%d`  # start date of forecast YYYY-MM-DD
   lognm=merge_log.txt
fi
dir0=$TP4_REALTIME_RES/ice_only/work/$tday/netcdf
odir=$TP4_REALTIME_RES/ice_only/work/$tday/final_output
#########################################################################

cd $dir0
rm -r -f  tmp
mkdir tmp

# LOG
log=$SWARP_ROUTINES/forecast_scripts/logs/$lognm
if [ -f "$log" ]
then
   rm $log
fi
touch $log

echo " "                                              >> $log
echo "Unpacking files (ncpdq -U)..."                  >> $log

flist=(*.nc)
f=${flist[0]}

Nfiles=0
for f in TP4archv_*.nc
do
   cf=${f:9:8}
   if [ "$cf" -ge "$tday" ]
   then
      Nfiles=$((Nfiles+1))
      echo Unpacking   $f
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
# change variable names to standard abbreviations in GRIB2 tables
ncrename -v fice,icec   $ofil #ice cover
ncrename -v hice,icetk  $ofil #ice thickness

# change time:units
ncatted -O -h -a units,time,m,c,"seconds since 1970-01-01T00:00:00Z"  $ofil

###########################################################################################
# most variable attributes set in hyc2proj:
# - just add shape of earth to "int stereographic"
#   (radius from hyc2proj - mod_toproj.F90)
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
ncatted -O -h -a forecast_start_date,global,c,c,"${tday_long}T00:00:00Z"      $ofil
# ncatted -O -h -a forecast_range,global,c,c,"5 day forecast"                   $ofil
ncatted -O -h -a forecast_range,global,c,c,"6 day forecast"                   $ofil
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
   ncatted -O -h -a bulletin_date,global,o,c,"$bdate UTC $btime" $ofil
fi

# Rename bulletin_date to restart_file_date (clearer)
ncrename -a global@bulletin_date,restart_file_date    $ofil

# delete old attribute(s)
ncatted -a field_date,global,d,,                      $ofil
ncatted -a history,global,d,,                         $ofil
###########################################################################################

###########################################################################################
# move to final location
mkdir -p $odir
mv $ofil $odir
###########################################################################################

###########################################################################################
# check number of records (use $Nfiles)
# 5 day, 3-hourly forecast -> 41
# 6 day, 3-hourly forecast -> 49
Ncorrect=49
if [ $Nfiles -ne $Ncorrect ]
then
   efil=swarp_tmp.txt
   echo Warning: merge_TP4archv.sh           >  $efil
   echo Wrong number of records in $ofil     >> $efil
   echo "($Nfiles -  should be $Ncorrect)"   >> $efil
   mail -s "WARNING: ice-only final product faulty" $email < $efil
   rm $efil
else
   efil=swarp_tmp.txt
   echo Confirmation: merge_TP4archv.sh                  >  $efil
   echo Correct number of records ($Ncorrect) in $ofil   >> $efil
   mail -s "Ice-only final product OK" $email < $efil
   rm $efil
fi
###########################################################################################
