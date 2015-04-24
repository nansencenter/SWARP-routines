#!/bin/bash
# This script will collect and archive the results of the local TP4 model

# EMAIL ADDRESS
# email="user1@domain.com,user2@domain.com,etc..."
# =============================================================================
email="gcmdnt90@gmail.com"
# =============================================================================

# defining all the dir that will be used
RTM=/work/timill/RealTime_Models
DFDIR=$RTM/TP4a0.12/expt_01.1
WKDIR=$RTM/results/TP4a0.12/ice_only/work
FCDIR=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts

# LOG
log=$FCDIR/logs/gather_log.txt
if [ -f $log ]
then
   rm $log
fi
touch $log

echo "Collecting data produced in date $2"               >> $log

tday=$1

TDIR=$WKDIR/$tday
mkdir -p $WKDIR/$tday/bin
mkdir -p $WKDIR/$tday/netcdf
mkdir -p $WKDIR/$tday/final_output
mkdir -p $WKDIR/$tday/info

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
cp $TP4_REALTIME/Build_V2.2.12_X01.1/flags $TDIR/info
if [ $? -eq 0 ]
then
   echo "flags file present"                           >> $log
else
   echo "flags file NOT present"                       >> $log
   mail -s "gather_FCresults FAILED" $email  < $log
fi
echo "Transfer complete"                               >> $log
