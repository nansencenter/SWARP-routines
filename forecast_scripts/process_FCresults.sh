#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

$SWARP_ROUTINES/netcdf_production/convert_TP4archv.sh
$SWARP_ROUTINES/netcdf_production/merge_TP4archv.sh
$SWARP_ROUTINES/forecast_scripts/gather_FCresults.sh
$SWARP_ROUTINES/forecast_scripts/collect_FCresults.sh

