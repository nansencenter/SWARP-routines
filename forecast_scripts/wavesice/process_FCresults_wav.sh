#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

# EMAIL ADDRESS
# email="user1@domain.com,user2@domain.com,etc..."
# ===================================================================================
email="gcmdnt90@gmail.com"
# ===================================================================================

datelist=$SWARP_ROUTINES/forecast_scripts/datelist.txt
if [ -f $datelist ]
then
   # set vbl's
   tday=$(cat $datelist | sed '1!d')
   tday_long=`date --date=$tday +%Y-%m-%d`

   # run scripts
   $SWARP_ROUTINES/forecast_scripts/wavesice/gather_FCresults_wav.sh  $tday # $(cat $datelist | sed '2!d') #$1 $2 
   $SWARP_ROUTINES/netcdf_production/convert_TP4archv_wav.sh $tday #$1
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv_wav.sh   $tday # $(cat $datelist | sed '2!d') #$1 $2
   $SWARP_ROUTINES/forecast_scripts/wavesice/collect_FCresults_wav.sh $tday # $(cat $datelist | sed '2!d') #$1 $2

   # finish up
   cp $datelist /work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $datelist | sed '1!d')/info/
else
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
