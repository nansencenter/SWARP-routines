#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model
datelist=$SWARP_ROUTINES/forecast_scripts/datalist

$SWARP_ROUTINES/netcdf_production/convert_TP4archv.sh $(cat $datelist | sed '1!d') #$1
$SWARP_ROUTINES/netcdf_production/merge_TP4archv.sh   $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
$SWARP_ROUTINES/forecast_scripts/gather_FCresults.sh  $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2 
$SWARP_ROUTINES/forecast_scripts/collect_FCresults.sh $(cat $datelist | sed '1!d') $(cat $datelist | sed '2!d') #$1 $2
rm datelist
