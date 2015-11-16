#!/bin/bash
#Get latest restart from the internal repo to the working dir

source $SWARP_ROUTINES/source_files/hex_vars.src

# ====================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ====================================================================================

# DIRECTORIES AND DATELIST
datelist=$FORECAST/wavesice/datelist.txt
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=$TP4_REALTIME/expt_01.1/data # where ice_only restart is
wdir=$TP4_REALTIME/expt_01.2/data # where to put restart
logdir=$FORECAST/logs
mkdir -p $logdir

# TEXTFILE AND LOG
in_restart=$FORECAST/ice_only/last_restart.txt
out_restart=$FORECAST/wavesice/last_restart.txt

# will just use the TOPAZ restart file used by ice_only
cat $in_restart > $out_restart
base_restart=`cat $in_restart`
TP4restart=$ddir/${base_restart}


log=$logdir/tp_get_wav_log.txt

if [ $(date +%A) == "Monday" ]
then
   # start a new log
   echo "$date ---------------------------------"  >  $log
   echo ""                                         >> $log
else
   # append to old log
   echo ""                                         >> $log
   echo "$date ---------------------------------"  >> $log
   echo ""                                         >> $log
fi

echo "Restart name: $base_restart"                          >> $log
echo "Looking for: $TP4restart"                             >> $log


if [ -f $wdir/${base_restart}.a ] && [ -f $wdir/${base_restart}.b ] && [ -f $wdir/${base_restart}ICE.uf ]
then
   # check waves-ice dir
   echo ""                                                  >> $log
   echo "Daily restarts already present"                    >> $log
elif [ -f $TP4restart.a ] && [ -f $TP4restart.b ] && [ -f ${TP4restart}ICE.uf ]
then
   # check ice-only dir
   cp ${TP4restart}* $wdir/
else
   echo "Daily restarts NOT found, check ice_only run"      >> $log
   mail -s "WARNING - daily ice_only restarts NOT found" $email < $log
fi

echo ""                                                     >> $log

# CREATING DAILY INFO DIR
idir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $datelist | sed '1!d')/info
mkdir -p $idir
cp $out_restart $idir/
