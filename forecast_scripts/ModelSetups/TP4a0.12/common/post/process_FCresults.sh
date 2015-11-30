#!/bin/bash
# This script will execute subscipts that will manage the products of the SWARP model

source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

# ===================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ===================================================================================

logdir=$THISFC/logs
post=$HREGcommon/post

datelist=$logdir/datelist.txt
if [ -f $datelist ]
then
   # set vbl's
   tday=$(cat $datelist | sed '1!d')

#  # run scripts
   $post/gather_FCresults.sh       $THIS_SRC $tday

   if [ $archv_opt -eq 1 ]
   then
      $post/convert_archv.sh       $THIS_SRC $tday
      $post/merge_archv.sh         $THIS_SRC $tday
   fi

   if [ $archv_wav_opt -eq 1 ]
   then
      $post/convert_archv_wav.sh    $THIS_SRC $tday
      $post/merge_archv_wav.sh      $THIS_SRC $tday
   fi

   $post/make_gifs.sh               $THIS_SRC   $tday

   # finish up
   cp $datelist $THISFC2/$tday/info
   cp $THIS_SRC $THISFC2/$tday/info

   $post/backup_FCresults.sh    $THIS_SRC $tday

else

   touch log.txt
   echo "DATELIST NOT FOUND" >> log.txt
   mail -s "Process forecast problems" $email < log.txt
   rm log.txt

fi
