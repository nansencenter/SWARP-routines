#!/bin/bash
# This script will collect and archive the results of the local TP4 model

source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

# =============================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# =============================================================================

# defining all the dir that will be used
xdir=$RTmods/$HYCOMreg/expt_01.$Xno

# LOG
log=$THISFC/logs/gather_log_wav.txt
if [ -f $log ]
then
   rm $log
fi
touch $log

tday=$2
tday_long=`date --date=$tday +%Y-%m-%d`
cyear=${tday:0:4}
echo "Collecting data produced in date $tday_long"

TDIR=$THISFC2/$tday
mkdir -p $TDIR/bin
mkdir -p $TDIR/netcdf
mkdir -p $TDIR/final_output
mkdir -p $TDIR/info

#moving ${rungen}archv_wav
rungen=${HYCOMreg:0:3}
if [ $archv_wav_opt -eq 1 ]
then
   echo "Moving the ${rungen}archv_wav*.[ab] and ${rungen}DAILY*.[ab]"      >> $log
   mv $xdir/data/${rungen}archv_wav* $TDIR/bin
   if [ $? -eq 0 ]
   then
      echo "${rungen}archv_wav* files present"        >> $log
   else
      echo "${rungen}archv_wav* files NOT present"    >> $log
      mail -s "$0 FAILED" $email  < $log
   fi
fi

#moving TP4archv
if [ $archv_opt -eq 1 ]
then
   echo "Moving the ${rungen}archv*.[ab]"       >> $log
   mv $xdir/data/${rungen}archv* $TDIR/bin
   if [ $? -eq 0 ]
   then
      echo "${rungen}archv* files present"      >> $log
   else
      echo "${rungen}archv* files NOT present"  >> $log
      mail -s "$0 FAILED" $email  < $log
   fi
fi

echo "Moving ${rungen}DAILY"
mv $xdir/data/${rungen}DAILY* $TDIR/bin
if [ $? -eq 0 ]
then
   echo "DAILY* files present"   >> $log
# else
#    echo "DAILY* files NOT present"                       >> $log
#    mail -s "gather_FCresults_wav FAILED" $email  < $log
fi

#moving the info files
echo "Moving the info files"     >> $log
cp $xdir/log/mpijob.out    $TDIR/info
cp $THISFC/inputs/flags    $TDIR/info
cp $THISFC/logs/*log.txt   $TDIR/info

echo "Transfer complete"         >> $log
