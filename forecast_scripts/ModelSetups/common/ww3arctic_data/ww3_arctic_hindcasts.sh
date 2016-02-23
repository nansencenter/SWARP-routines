#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC="$FORECAST/ModelSetups/common/ww3arctic"
fget="$THISFC/ww3_arctic_download.sh"
fsort="$THISFC/ww3_arctic_sort.sh"
ww3a=$wmsc/WAVES_INPUT/WW3_ARCTIC
cycle=3 # get hindcasts

#set manually:
date0=20141231 #refdate
date1=20150101 #start date
date2=20151107 #final date

J1=`date --date="$date1" +%j`
J2=`date --date="$date2" +%j`
# for jday in `seq $J1 $J2`
for jday in `seq $J1 $J1`
do
   dtj=`date --date="$date0 +${jday}days" +%Y%m%d`
   echo ""
   echo $fget  $dtj $cycle
   $fget  $dtj $cycle

   # echo ""
   # echo $fsort $dtj $cycle
   # $fsort $dtj $cycle
done
