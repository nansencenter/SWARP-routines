#!/usr/bin/ksh
# ======================================================================
# Six arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (-1 at 00:00UTC)
# 4. to_day (+0)
# 5. to_day (+2 at 00:00 and 12:00UTC)
# Other infile parameters must be changed in infile.mal.outer (TP4) or infile.mal
# ======================================================================

# EMAIL ADDRESS
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
# ======================================================================
email=$(cat $address)
# ======================================================================

# CREATING DIRECTORIES
SWARP_ROUTINES=$HOME/GITHUB-REPOSITORIES/SWARP-routines
TP4_REALTIME=/work/timill/RealTime_Models/TP4a0.12
MAINDIR=$SWARP_ROUTINES/forecast_scripts
xdir=/work/timill/RealTime_Models/TP4a0.12/expt_01.1

# CREATING LOG
logdir=$MAINDIR/logs
mkdir -p $logdir
log=$logdir/mk_infile_log.txt
touch $log

if [ $# -ne 5 ]
then
  echo "Usage: $0 <5 inputs> (see script)"                  >> $log
  mail -s "Make_infile4forecast.sh FAILED" $email < $log
  exit
else
  echo "-------------------------------------------------------------------------"  >> $log 
  echo " make_infile4forecast.sh $1"                                                >> $log
  echo "-------------------------------------------------------------------------"  >> $log
  rungen=$1
  thr=`printf '%3.3d' $3`	#three digits modification
  fou=`printf '%3.3d' $4`	#three digits modification
  fiv=`printf '%3.3d' $5`	#three digits modification
  cd ${MAINDIR}/inputs
  if [ -z "${rungen}" -o "${rungen}" == "TP4" ]
  then
     file=infile.mal.outer
  else
     file=infile.mal
  fi
  
  if [ -f $file ]
  then
    echo "Run version:    $1"                                                       >> $log
    echo "Reference year: $2"                                                       >> $log
    echo "From day:       $thr"                                                     >> $log 
    echo "To day 1:       $fou"                                                     >> $log 
    echo "To day 2:       $fiv"                                                     >> $log 
    if [ -f $xdir/infile.in ]
    then
      rm -f $xdir/infile.in
    fi

    cat $file | sed \
    -e "s/nd1/$thr/g" \
    -e "s/nd2/$fou/g" \
    -e "s/nd3/$fiv/g" \
    -e "s/ref_year/$2/g" \
    -e "s/run_gen/$1/g" \
    > $xdir/infile.in

    echo " $xdir/infile.in created"                                                 >> $log
  else
    echo " ERROR : infile.mal does not exist"                                       >> $log
    mail -s "Make_infile4forecast.sh FAILED" $email < $log
    exit
  fi
fi

# CLEANING THE LOG EVERY WEEK
if [ $(date +%A) == "Monday" ]
then
   rm $log
fi


