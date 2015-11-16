#!/usr/bin/ksh
# ======================================================================
# Six arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (0)
# 4. to_day (+2.5)
# Other infile parameters must be changed in infile.mal.outer (TP4) or infile.mal
# ======================================================================

# ======================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ======================================================================

# CREATING DIRECTORIES
SWARP_ROUTINES=$HOME/GITHUB-REPOSITORIES/SWARP-routines
TP4_REALTIME=/work/timill/RealTime_Models/TP4a0.12
MAINDIR=$SWARP_ROUTINES/forecast_scripts
xdir=/work/timill/RealTime_Models/TP4a0.12/expt_01.2

# CREATING LOG
logdir=$MAINDIR/logs
mkdir -p $logdir
log=$logdir/mk_infile_log.txt
touch $log

if [ $# -ne 4 ]
then
  echo "Usage: $0 <4 inputs> (see script)"                  >> $log
  mail -s "Make_infile4forecast.sh FAILED" $email < $log
  exit
else
  echo "-------------------------------------------------------------------------"  >> $log 
  echo " make_infile4forecast.sh $1"                                                >> $log
  echo "-------------------------------------------------------------------------"  >> $log
  rungen=$1
  thr=`printf '%3.3d' $3`	#three digits modification
  fou=`printf '%3.3d' $4`	#three digits modification
  cd ${MAINDIR}/inputs/wavesice
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
    if [ -f $xdir/infile.in ]
    then
      rm -f $xdir/infile.in
    fi

    cat $file | sed \
    -e "s/nd1/$thr/g" \
    -e "s/nd2/$fou/g" \
    -e "s/ref_year/$2/g" \
    -e "s/run_gen/$1/g" \
    > $xdir/infile.in

    echo " $xdir/infile.in created"                                                 >> $log
  else
    echo " ERROR : infile.mal does not exist"                                       >> $log
    mail -s "Make_infile4forecast_wav.sh FAILED" $email < $log
    exit
  fi
fi

# CLEANING THE LOG EVERY WEEK
if [ $(date +%A) == "Monday" ]
then
   rm $log
fi

