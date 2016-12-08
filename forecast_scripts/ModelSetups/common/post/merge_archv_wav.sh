#!/bin/bash
# script to extract some variables from a netcdf file

source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

# ===================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ===================================================================================

#########################################################################
# inputs are start date of forecast
if [ $# -eq 0 ]
then
   # do a test
   tday=`date +%Y%m%d`                       # start date of forecast YYYYMMDD
   echo \${tday}=$tday
   tday_long=`date --date=$tday +%Y-%m-%d`   # start date of forecast YYYY-MM-DD
   lognm=merge_wav_log_test.txt
else
   tday=$2                                   # start date of forecast YYYYMMDD
   tday_long=`date --date=$tday +%Y-%m-%d`   # start date of forecast YYYY-MM-DD
   # tday_long=`date --date="$2" +%Y-%m-%d`  # start date of forecast YYYY-MM-DD
   lognm=merge_wav_log.txt
fi
#########################################################################

#########################################################################
dir0=$THISFC2/$tday/netcdf
odir=$THISFC2/$tday/final_output  
#########################################################################

#extract start date/time of forecast
# and put it into output filename

cd $dir0
mkdir -p tmp

# LOG
logdir=$THISFC/logs
log=$logdir/$lognm
if [ -f "$log" ]
then
   rm $log
fi
touch $log

echo $date                                            >> $log
echo " "                                              >> $log
echo "Unpacking files (ncpdq -U)..."                  >> $log

Nfiles=0
for f in ${rungen}archv_wav*.nc
do
   end=${f#*dump}
   fdate=${end:0:8}
   if [ $fdate -ge $tday ]
   then
      Nfiles=$((Nfiles+1))
      echo Unpacking   $f
      ncpdq -U $f tmp/$f
   fi
done

#combine unpacked files
echo " "                                              >> $log
echo "Combining unpacked files (ncrcat)..."           >> $log
ncrcat tmp/*.nc tmp.nc

# get 3d variables
# - fix missing and fill values
ncdump -h tmp.nc |grep "(time," > tmp.txt
cat tmp.txt | while read line
do
   Line=($line)
   v=${Line[1]}
   v=${v%(*}
   ncatted -O -a _FillValue,$v,o,s,-32767      tmp.nc
   ncatted -O -a missing_value,$v,o,s,-32767   tmp.nc
done
rm tmp.txt

ncatted -O -a _FillValue,model_depth,o,s,-32767 tmp.nc


#set name of output file
ofil=${FC_OUTPUT}_start${tday}T000000Z.nc

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
ncatted -O -h -a units,time,m,c,"seconds since 1970-01-01T00:00:00Z" $ofil
ncatted -O -h -a calendar,time,c,c,"gregorian"                        $ofil
#########################################################################

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

if [ $rungen == "TP4" ]
then
   ncatted -O -h -a area_name,global,c,c,"TP4a0.12"                                 $ofil
   ncatted -O -h -a area_resolution,global,c,c,"12.5km"                             $ofil
   ncatted -O -h -a area_description,global,c,c,"Arctic and North Atlantic Oceans"  $ofil
elif [ $rungen == "BS1" ]
then
   ncatted -O -h -a area_name,global,c,c,"BS1a0.045"                    $ofil
   ncatted -O -h -a area_resolution,global,c,c,"4.5km"                  $ofil
   ncatted -O -h -a area_description,global,c,c,"Barents and Kara Sea"  $ofil
elif [ $rungen == "FR1" ]
then
   ncatted -O -h -a area_name,global,c,c,"FR1a0.03"            $ofil
   ncatted -O -h -a area_resolution,global,c,c,"3.5km"         $ofil
   ncatted -O -h -a area_description,global,c,c,"Fram Strait"  $ofil
fi

if [ $FCtype == "wavesice" ]
then

   ncatted -O -h -a title,global,o,c,"SWARP waves-in-ice forecast"            $ofil # o=overwrite/create, c=format (also f=float)
   ncatted -O -h -a field_type,global,c,c,"6-hourly"                          $ofil
   ncatted -O -h -a forecast_range,global,c,c,"2.5 day forecast"              $ofil
   ncatted -O -h -a wave_forcing,global,c,c,"WAM North Sea Arctic (met.no)"   $ofil
   ncatted -O -h -a wave_forcing_contact,global,c,c,"bruce.hackett@met.no"    $ofil

elif [ $FCtype == "wavesice_ww3arctic" ]
then

   ncatted -O -h -a title,global,o,c,"SWARP waves-in-ice forecast"            $ofil # o=overwrite/create, c=format (also f=float)
   ncatted -O -h -a field_type,global,c,c,"3-hourly"                          $ofil
   ncatted -O -h -a forecast_range,global,c,c,"6 day forecast"                $ofil
   ncatted -O -h -a wave_forcing,global,c,c,"WW3 Arctic"                      $ofil
   #
   wref="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3/iowaga/HINDCAST/ARCTIC"
   ncatted -O -h -a wave_forcing_reference,global,c,c,$wref                   $ofil

elif [ $FCtype == "ice_only" ]
then

   ncatted -O -h -a title,global,o,c,"SWARP ice-only forecast"                $ofil # o=overwrite/create, c=format (also f=float)
   ncatted -O -h -a field_type,global,c,c,"3-hourly"                          $ofil
   ncatted -O -h -a forecast_range,global,c,c,"6 day forecast"                $ofil
   ncatted -O -h -a wave_forcing,global,c,c,"None"                            $ofil

fi

ncatted -O -h -a forecast_start_date,global,c,c,"${tday_long}T00:00:00Z"      $ofil
ncatted -O -h -a institution,global,c,c,"NERSC"                               $ofil
ncatted -O -h -a institution_references,global,c,c,"http://www.nersc.no/"     $ofil
ncatted -O -h -a data_centre,global,c,c,"NERSC"                               $ofil
ncatted -O -h -a data_centre_references,global,c,c,"www.nersc.no"             $ofil
ncatted -O -h -a contact,global,c,c,"timothy.williams@nersc.no"               $ofil
ncatted -O -h -a project,global,c,c,"SWARP"                                   $ofil
ncatted -O -h -a project_references,global,c,c,"swarp.nersc.no"               $ofil
ncatted -O -h -a distribution_statement,global,c,c,"No restrictions"          $ofil
ncatted -O -h -a operational_status,global,c,c,"test"                         $ofil

# wave forcing info
#
# ncatted -O -h -a history,global,o,c,"NERSC-HYCOM output->hyc2proj->ncrcat"    $ofil

# Restart file date
reformat=0
if [ $reformat -eq 1 ]
then
   # reformat bulletin date
   # eg 2015-05-09T00:00:00Z -> 2015-05-09 UTC 00:00:00
   bdate_info=$(ncinfo $ofil | grep "bulletin_date")
   split=(${bdate_info//\ / }) # make array with delimiter space "\ "
   bdate_time=${split[1]}       # eg 2015-05-09T00:00:00Z
   split=(${bdate_time//T/ }) # make array with delimiter "T" 
   bdate=${split[0]}               # eg 2015-05-09
   btime=${split[1]}               # eg 00:00:00Z
   btime=${btime:0:8} # remove "Z"
   ncatted -O -h -a bulletin_date,global,o,c,"$bdate UTC $btime" $ofil   # add new attribute "restart_file_date"
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
# 2.5 day,  6-hourly forecast -> 11
# 5 day,    3-hourly forecast -> 41
# 6 day,    3-hourly forecast -> 49
Ncorrect=$FCrecords
if [ $Nfiles -ne $Ncorrect ]
then
   efil=swarp_tmp.txt
   echo Warning: merge_archv_wav.sh                            >  $efil
   echo Wrong number of records in $ofil                       >> $efil
   echo "($Nfiles -  should be $Ncorrect)"                     >> $efil
   mail -s "WARNING: $rungen $FCtype_long final product faulty" $email <  $efil
   rm $efil
else
   efil=swarp_tmp.txt
   echo Confirmation: merge_archv_wav.sh                       >  $efil
   echo "Correct number of records ($Ncorrect) in $ofil"       >> $efil
   mail -s "$rungen $FCtype_long final product OK" $email              <  $efil
   rm $efil
fi
###########################################################################################
