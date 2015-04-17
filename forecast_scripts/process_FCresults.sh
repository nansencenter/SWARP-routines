#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model
datelist=$SWARP_ROUTINES/forecast_scripts/datelist
if [ -f $datelist ]
then
   $SWARP_ROUTINES/netcdf_production/convert_TP4archv.sh $(cat $datelist | sed '1!d') #$1
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv.sh   $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
   $SWARP_ROUTINES/forecast_scripts/gather_FCresults.sh  $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2 
   $SWARP_ROUTINES/forecast_scripts/collect_FCresults.sh $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
else
   mail -s "process forecast problems" gcmdnt90@gmail.com < echo "datelist file not found, check ASAP"
fi
