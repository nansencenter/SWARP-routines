#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

# AS LONG AS THE OTHER SCRIPTS WORK I DON'T SEE HOW THIS SIMPLE ONE COULD FAIL.
#TODO FOR FUTURE ALLERTS?

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC
# ===================================================================================

tday=$1
tday_long=`date --date=$tday +%Y-%m-%d`
cyear=`date --date=$tday +%Y`
echo "Collecting data produced in date $tday_long"


cd $THISFC2

echo "The archive file name will be "
tfil=SWARP_ice_only_forecast_$tday.tar.gz
echo " $tfil "
tar -zcvf $tfil -C $THISFC2 $tday
chmod 666 $tfil

mkdir -p $THISFC3/$cyear
mv $tfil $THISFC3/$cyear
echo "SWARP products of $tday"
echo "stored in $THISFC3/$cyear"

