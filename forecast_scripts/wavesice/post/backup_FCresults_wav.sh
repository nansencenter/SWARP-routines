#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

# AS LONG AS THE OTHER SCRIPTS WORK I DON'T SEE HOW THIS SIMPLE ONE COULD FAIL.
#TODO FOR FUTURE ALERTS?

source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

tday=$1
tday_long=`date --date=$tday +%Y-%m-%d`
cyear=${tday:0:4}
echo "Collecting data produced in date $tday_long"


cd $THISFC2
tfil=${FC_OUTPUT}_$tday.tar.gz

echo "The archive file name will be "
echo " $tfil "
touch $tfil
tar -zcvf $tfil -C $THISFC2 $tday
chmod 777 $tfil
mkdir -p $THISFC3/$cyear
mv $tfil $THISFC3/$cyear/
echo "SWARP products of $tday"
echo "stored in $THISFC3/$cyear"
