#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

echo "Collecting data produced in date `date +%d/%m/%Y`"

RTM=/work/timill/RealTime_Models
WDIR=$RTM/results/TP4a0.12/ice_only/work
SDIR=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/ice_only

cyear=`date +%Y`
tday=`date +%Y%m%d`

cd $WDIR

echo "The archive file name will be "
tfil=SWARP_forecast_$tday.tar.gz
echo " $tfil "
touch $tfil
tar -zcvf $tfil -C $WDIR $tday
mkdir -p $SDIR/$cyear
mv $tfil $SDIR/$cyear
echo "SWARP products of $tday"
echo "stored in $SDIR/$cyear"
