#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

echo "Collecting data produced in date `date +%d/%m/%Y`"

RTM=/work/timill/RealTime_Models
WDIR=$RTM/results/TP4a.0.12/ice_only/work
SDIR=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecast/ice_only

cd $WDIR

cyear=`date +%Y`
tday=`date +%d_%m_%Y`
TDIR=$WDIR/$tday

echo "The archive file name will be "
tfil=SWARP_forecast_$tday.tar.gz
echo " "
touch $tfil
tar -zcvf $tfil $TDIR
mkdir -p $SDIR/$cyear
mv $tfil $SDIR/$cyear
echo "SWARP products of $tday"
echo "stored in $SDIR/$cyear"

