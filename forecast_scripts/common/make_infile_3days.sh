#!/usr/bin/ksh
# ======================================================================
# Four arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (last Monday at 00:00UTC)
# 5. to_day (+6 at 00:00UTC)
# Other infile parameters must be changed in infile.mal.outer (TP4) or infile.mal
# ======================================================================

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

# other inputs
reg=$2
ryear=$3
day1=$4
day2=$5
day3=$6
hour=$7
#
nest_outer="F"
if [ $# -ge 8 ]
then
   nest_outer="$8"
fi
#
nest_inner="F"
if [ $# -ge 9 ]
then
   nest_inner="$9"
fi
# ===================================================================================

# ======================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ======================================================================


# CREATING DIRECTORIES
xdir=$TP4_REALTIME/expt_01.$Xno

# CREATING LOG
logdir=$THISFC/logs
mkdir -p $logdir
log=$logdir/mk_infile_log.txt
touch $log

if [ $# -lt 7 ]
then
  echo "Usage: $0 <4 inputs> (see script)"          >> $log
  mail -s "make_infile_3days.sh FAILED"      $email <  $log
  exit
else
  echo "-------------------------------------------------------------------------"  >> $log 
  echo " make_infile_3days.sh ($reg)"                                               >> $log
  echo "-------------------------------------------------------------------------"  >> $log
  file=$FCcommon/inputs/infile.3days
  
  if [ -f $file ]
  then

    echo "Run version               : $reg"           >> $log
    echo "Reference year            : $ryear"         >> $log
    echo "From day 1                : $day1"          >> $log 
    echo "to day 2 (make restart)   : $day2"          >> $log 
    echo "To day 3                  : ${day3}_$hour"  >> $log 
    rm -f $xdir/infile.in

    cat $file | sed \
    -e "s/run_gen/$reg/g" \
    -e "s/ref_year/$ryear/g" \
    -e "s/run_gen/$reg/g" \
    -e "s/nd1/$day1/g" \
    -e "s/nd2/$day2/g" \
    -e "s/nd3/$day3/g" \
    -e "s/HH/$hour/g" \
    -e "s/TF/$nest_outer /g" \
    -e "s/FT/$nest_inner /g" \
    > $xdir/infile.in

    echo " $xdir/infile.in created"                >> $log
  else
     echo " ERROR : $file does not exist"          >> $log
    mail -s "Make_infile_3days.sh FAILED" $email  <  $log
    exit
  fi
fi

# CLEANING THE LOG EVERY WEEK
if [ $(date +%A) == "Monday" ]
then
   rm $log
fi
