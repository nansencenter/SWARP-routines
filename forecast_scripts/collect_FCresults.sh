#!/bin/bash
# This script will collect and archive the results of the local TP4 model

echo "Collecting data produced in date `date +%d/%m/%Y`"

RTM=/work/timill/RealTime_Models
DFDIR=$RTM/TP4a0.12/expt_01.1/data
NCDIR=$RTM/post_processing/archv_netcdf
STDIR=$RTM/results/TP4a.0.12/ice_only

cd $WDIR

tpday=( TP4DAILY*  )
tparchv=( TP4archv.* )
iy=${tpday[0]:9:4}
id=${tpday[0]:14:3}
echo "The archive file name will be "
tfil=TP4DAILYarchv.$iy_$fd.tar.gz
echo " "
touch $tfil
tar -zcvf $tfil $tpday $tparchv
mkdir -p $SDIR/${yy}
mv $tfil $SDIR/${yy}
echo "Daily and archv files from $fd/$fy"
echo "stored in $SDIR"

