#!/bin/bash
# This script will collect and archive the results of the local TP4 model

source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_ww3arctic
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

# =============================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# =============================================================================

# defining all the dir that will be used
DFDIR=$TP4_REALTIME/expt_01.$Xno

# LOG
log=$THISFC/logs/gather_log_ww3a.txt
if [ -f $log ]
then
   rm $log
fi
touch $log

tday=$1
tday_long=`date --date=$tday +%Y-%m-%d`
cyear=${tday:0:4}
echo "Collecting data produced in date $tday_long"

TDIR=$THISFC2/$tday
mkdir -p $TDIR/bin
mkdir -p $TDIR/netcdf
mkdir -p $TDIR/final_output
mkdir -p $TDIR/info

# #moving TP4restart
# echo "Moving the TP4restarts"                            >> $log
# cp $DFDIR/data/TP4restart* $TDIR/bin
# if [ $? -eq 0 ]
# then
#    echo "Restart* files present"                           >> $log
# else
#    echo "Restart* files NOT present"                       >> $log
#    mail -s "gather_FCresults_wav FAILED" $email  < $log
# fi

#moving TP4archv & TP4DAILY
echo "Moving the TP4archv*.[ab] and TP4DAILY*.[ab]"      >> $log
mv $DFDIR/data/TP4archv_wav* $TDIR/bin
if [ $? -eq 0 ]
then
   echo "Archv_wav* files present"                           >> $log
else
   echo "Archv_wav* files NOT present"                       >> $log
   mail -s "gather_FCresults_wav FAILED" $email  < $log
fi

mv $DFDIR/data/TP4DAILY* $TDIR/bin
if [ $? -eq 0 ]
then
   echo "DAILY* files present"                           >> $log
#else
#   echo "DAILY* files NOT present"                       >> $log
#   mail -s "gather_FCresults_wav FAILED" $email  < $log
fi

#moving the info files
echo "Moving the info files"                             >> $log
cp $DFDIR/log/mpijob.out $TDIR/info
if [ $? -eq 0 ]
then
   echo "mpijob file present"                           >> $log
else
   echo "mpijob file NOT present"                       >> $log
   mail -s "gather_FCresults_wav FAILED" $email  < $log
fi
cp $TP4_REALTIME/Build_V2.2.12_X01.$Xno/flags $TDIR/info
if [ $? -eq 0 ]
then
   echo "flags file present"                           >> $log
else
   echo "flags file NOT present"                       >> $log
   mail -s "gather_FCresults_ww3a FAILED" $email  < $log
fi
echo "Transfer complete"                               >> $log