#!/bin/bash
# This script will collect and archive the results of the local TP4 model

echo "Collecting data produced in date `date +%d/%m/%Y`"

# defining all the dir that will be used
RTM=/work/timill/RealTime_Models
DFDIR=$RTM/TP4a0.12/expt_01.1/data
NCDIR=$RTM/post_processing/archv_netcdf
FODIR=$RTM/results/TP4a0.12/ice_only/final_output
WKDIR=$RTM/results/TP4a0.12/ice_only/work

tday=`date +%d_%m_%Y`
TDIR=$WKDIR/$tday
mkdir -p $WKDIR/$tday/bin
mkdir -p $WKDIR/$tday/netcdf
mkdir -p $WKDIR/$tday/final_output
mkdir -p $WKDIR/$tday/info

#moving TP4archv & TP4DAILY
mv $DFDIR/TP4archv* $TDIR/bin
mv $DFDIR/TP4DAILY* $TDIR/bin

#moving *.nc files
mv $NCDIR/TP4archv*.nc $TDIR/netcdf
mv $NCDIR/SWARP*.nc $TDIR/final_output

#moving the info files
TODO


