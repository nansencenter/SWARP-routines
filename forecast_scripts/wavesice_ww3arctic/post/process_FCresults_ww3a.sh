#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_ww3arctic
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

# ===================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ===================================================================================

logdir=$THISFC/logs
post=$THISFC/post

datelist=$logdir/datelist.txt
if [ -f $datelist ]
then
   # set vbl's
   tday=$(cat $datelist | sed '1!d')

   # run scripts
   $post/gather_FCresults_ww3a.sh   $tday
   $post/convert_TP4archv_ww3a.sh   $tday
   $post/merge_TP4archv_ww3a.sh     $tday
   $post/backup_FCresults_ww3a.sh   $tday

   # finish up
   cp $datelist $THISFC2/$tday/info
else
   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt
fi
