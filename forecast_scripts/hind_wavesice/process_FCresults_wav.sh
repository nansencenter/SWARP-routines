#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

# EMAIL ADDRESS
# email="user1@domain.com,user2@domain.com,etc..."
# ===================================================================================
email="gcmdnt90@gmail.com"
# ===================================================================================

hindi_datelist=$SWARP_ROUTINES/forecast_scripts/hindi_datelist.txt
if [ -f $hindi_datelist ]
then
   $SWARP_ROUTINES/forecast_scripts/wavesice/gather_FCresults_wav.sh  $(cat $hindi_datelist | sed '1!d') $(cat $hindi_datelist | sed '2!d') #$1 $2 
   $SWARP_ROUTINES/netcdf_production/convert_TP4archv_wav.sh $(cat $hindi_datelist | sed '1!d') #$1
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv_wav.sh   $(cat $hindi_datelist | sed '1!d') $(cat $hindi_datelist | sed '2!d') #$1 $2
   $SWARP_ROUTINES/forecast_scripts/wavesice/collect_FCresults_wav.sh $(cat $hindi_datelist | sed '1!d') $(cat $hindi_datelist | sed '2!d') #$1 $2
   cp $hindi_datelist /work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $hindi_datelist | sed '1!d')/info/
else
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
