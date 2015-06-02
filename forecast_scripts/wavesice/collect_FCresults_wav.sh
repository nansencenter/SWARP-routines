#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

# AS LONG AS THE OTHER SCRIPTS WORK I DON'T SEE HOW THIS SIMPLE ONE COULD FAIL.
#TODO FOR FUTURE ALLERTS?

echo "Collecting data produced in date $2"

RTM=/work/timill/RealTime_Models
WDIR=$RTM/results/TP4a0.12/wavesice/work
SDIR=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/wavesice

tday=$1
cyear=${tday::4}

cd $WDIR

echo "The archive file name will be "
tfil=SWARP_wavesice_forecast_$tday.tar.gz
echo " $tfil "
touch $tfil
tar -zcvf $tfil -C $WDIR $tday
mkdir -p $SDIR/$cyear
mv $tfil $SDIR/$cyear/
echo "SWARP products of $tday"
echo "stored in $SDIR/$cyear"

