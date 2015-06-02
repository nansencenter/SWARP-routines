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
   $SWARP_ROUTINES/forecast_scripts/ice_only/gather_FCresults.sh  $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2 
   $SWARP_ROUTINES/netcdf_production/convert_TP4archv.sh $(cat $datelist | sed '1!d') #$1
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv.sh   $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
   $SWARP_ROUTINES/forecast_scripts/ice_only/collect_FCresults.sh $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
   cp $datelist /work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$(cat $datelist | sed '1!d')/info/
else
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
