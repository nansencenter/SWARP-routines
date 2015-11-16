#!/bin/bash
# This script will collect and archive the results of the local TP4 model

# =============================================================================
# EMAIL ADDRESS
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
email=$(cat $address)
# =============================================================================

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC
# ===================================================================================

# defining all the dir that will be used
source $SWARP_ROUTINES/source_files/hex_vars.src
FCDIR=$SWARP_ROUTINES/forecast_scripts
DFDIR=$TP4_REALTIME/expt_01.$Xno

# LOG
log=$FCDIR/logs/gather_log.txt
if [ -f $log ]
then
   rm $log
fi
touch $log

tday=$1
tday_long=`date --date=$tday +%Y-%m-%d`
echo "Collecting data produced in date $tday_long"       >> $log

TDIR=$THISFC2/$tday
mkdir -p $THISFC2/$tday/bin
mkdir -p $THISFC2/$tday/netcdf
mkdir -p $THISFC2/$tday/final_output
mkdir -p $THISFC2/$tday/info

#moving TP4archv & TP4DAILY
echo "Moving the TP4archv*.[ab] and TP4DAILY*.[ab]"      >> $log
mv $DFDIR/data/TP4archv* $TDIR/bin
if [ $? -eq 0 ]
then
   echo "Archv* files present"                           >> $log
else
   echo "Archv* files NOT present"                       >> $log
   mail -s "gather_FCresults FAILED" $email  < $log
fi
mv $DFDIR/data/TP4DAILY* $TDIR/bin
if [ $? -eq 0 ]
then
   echo "DAILY* files present"                           >> $log
else
   echo "DAILY* files NOT present"                       >> $log
   mail -s "gather_FCresults FAILED" $email  < $log
fi

#moving the info files
echo "Moving the info files"                             >> $log
cp $DFDIR/log/mpijob.out $TDIR/info
if [ $? -eq 0 ]
then
   echo "mpijob file present"                           >> $log
else
   echo "mpijob file NOT present"                       >> $log
   mail -s "gather_FCresults FAILED" $email  < $log
fi
cp $TP4_REALTIME/Build_V2.2.12_X01.$Xno/flags $TDIR/info
if [ $? -eq 0 ]
then
   echo "flags file present"                           >> $log
else
   echo "flags file NOT present"                       >> $log
   mail -s "gather_FCresults FAILED" $email  < $log
fi
echo "Transfer complete"                               >> $log

