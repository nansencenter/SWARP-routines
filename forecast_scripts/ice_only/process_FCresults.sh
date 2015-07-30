#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

# EMAIL ADDRESS
# email="user1@domain.com,user2@domain.com,etc..."
# ===================================================================================
email="gcmdnt90@gmail.com"
# ===================================================================================

datelist=$SWARP_ROUTINES/forecast_scripts/ice_only/datelist.txt
if [ -f $datelist ]
then
   # set var's
   tday=$(cat $datelist | sed '1!d')      # first day in YYYYMMDD format
   # tday_long=$(cat $datelist | sed '2!d') # first day in YYYY-MM-DD format TODO shouldn't need this (just convert $tday)

   # run scripts
   $SWARP_ROUTINES/forecast_scripts/ice_only/gather_FCresults.sh  $tday # $day1_long #$1 $2 
   $SWARP_ROUTINES/netcdf_production/convert_TP4archv.sh $tday #$1
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv.sh   $tday # $day1_long #$1 $2
   $SWARP_ROUTINES/forecast_scripts/ice_only/collect_FCresults.sh $tday # $day1_long #$1 $2

   # add datelist.txt to info
   cp $datelist /work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$tday/info/
else
   # send alert email
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
