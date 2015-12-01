#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

# ===================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ===================================================================================

logdir=$THISFC/logs
post=$THISFC/post

datelist=$logdir/datelist.txt
if [ -f $datelist ]
then
   # set var's
   tday=$(cat $datelist | sed '1!d')      # first day in YYYYMMDD format

   # ======================================================================================
   # 1. get all TP4archv.*.[ab] and TP4DAILY*.[ab]
   $post/gather_FCresults.sh  $tday

   # 2. convert all the TP4archv.*.[ab] to netcdf
   $post/convert_TP4archv.sh $tday #$1

   # 3. merge some of the netcdf files into one
   #     (only those later than today),
   #       add necessary attributes,
   #        and copy final output to correct location
   $post/merge_TP4archv.sh   $tday # $day1_long #$1 $2

   # 4. make gifs
   $post/make_gifs.sh $tday

   # 5. save to migrate
   # TODO change to norstore
   $post/backup_FCresults.sh $tday # $day1_long #$1 $2
   # ======================================================================================

   # finish up - add datelist.txt to info
   cp $datelist $THISFC2/$tday/info
else
   # send alert email
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
