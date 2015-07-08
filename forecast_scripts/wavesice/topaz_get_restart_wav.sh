#!/bin/bash
#Get latest restart from the internal repo to the working dir

# EMAIL ADDRESS
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
# ====================================================================================
email=$(cat $address)
# ====================================================================================

# DIRECTORIES AND DATELIST
fcdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts
datelist=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/datelist.txt
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data # where ice_only forecast will be done
wdir=/work/timill/RealTime_Models/TP4a0.12/expt_01.2/data # where waves_in_ice forecast will be done
logdir=$fcdir/logs
mkdir -p $logdir

# TEXTFILE AND LOG
out_restart=$fcdir/last_restart.txt
log=$logdir/tp_get_wav_log.txt

if [ $(date +%A) == "Monday" ]
then
   rm $log
fi

touch $log
echo $date  >> $log
echo ""     >> $log

cyear=$(cat $datelist | sed '3!d')			# current year
cmon=$(cat $datelist | sed '4!d')			# current month
cday=$(cat $datelist | sed '6!d')                       # current day julian (1 Jan = 1)
pyear=$(expr $cyear - 1)		              	# previous year
pday=$(expr $cday - 1)
pday=`printf '%3.3d' $pday`

base_restart=TP4restart${cyear}_${cday}_00
TP4restart=$ddir/${base_restart}

echo "Restart name: $base_restart"                          >> $log
echo "Looking for: $TP4restart"                             >> $log

if [ -f $TP4restart.a ] && [ -f $TP4restart.b ] && [ -f ${TP4restart}ICE.uf ]
then
   if [ $cday == 000 ]
   then
      rm $wdir/TP4restart${pyear}*
   else
      rm $wdir/TP4restart${cyear}_${pday}*
   fi
   mv ${TP4restart}* $wdir/
elif [ -f $wdir/${base_restart}.a ] && [ -f $wdir/${base_restart}.b ] && [ -f $wdir/${base_restart}ICE.uf ]
then
   echo ""                                                  >> $log
   echo "Daily restarts already present"                    >> $log
else
   echo "Daily restarts NOT found, check ice_only run"      >> $log
   if [ ! -f /work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/${cday}/final_output/*.nc ]
   then
      mail -s "WARNING - daily ice_only restarts NOT found" $email < $log
   fi
fi

echo ""                                                     >> $log

echo $base_restart > $out_restart

# CREATING DAILY INFO DIR
idir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $datelist | sed '1!d')/info
mkdir -p $idir
mv $out_restart $idir/

