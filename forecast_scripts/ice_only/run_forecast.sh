#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh

# get variables
# ($SWARP_ROUTINES defined in crontab or .bash_profile)
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

test_pre=0
print_info=1 # print info to screen (or email in crontab)


# ================================================================================================
# 1. First check if forecast is already running
# NB Do this before datelist.txt is changed 
msg=`$qstat | grep TP4x01${Xno}fc`
if [ ${#msg} -ne 0 ]
then
   if [ $print_info -eq 1 ]
   then
      echo "Ice-only FC is already running - stopping"
      echo "pbs job message:"
      echo $msg
   fi
   exit
fi
# ================================================================================================


manual=0
if [ $# -eq 1 ]
then
   echo ""
   echo "Running manually (no checks)"
   echo ""
   manual=1
fi


# ================================================================================================
# FORECAST DAYS
fc_days=6         # to be consistent with WW3 forecast            
final_hour="00"   # finishes at this hour on the final day
# ================================================================================================


# CREATING THE LOG
logdir=$THISFC/logs
mkdir -p $logdir
log=$logdir/run_forecast_log.txt
if [ -f "$log" ]
then
   rm $log
fi

#=======================================================================
# DATE INFO - TODAY
if [ 1 -eq 1 ]
then
   cday=$(date +%Y%m%d)
else
   # test last Tuesday (restart and dump day are the same)
   cday=$(date --date="last Tuesday" +%Y%m%d)
fi
cyear=$(date --date=$cday +%Y)
jday0=$(date --date=$cday +%j)
jday_today0=$(( $jday0 - 1))              # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0) # 3 digits

rundir=$THISFC2/$cday                     # where the results will end up
mkdir -p $rundir
mkdir -p $rundir/info

# CREATING THE DATELIST
# NB don't make one in $FORECAST since both
# ice-only and waves-ice FC's use it
# - if they run at the same time it gets messed up
datelist=$rundir/info/datelist.txt
echo $cday                             >  $datelist
echo $(date --date=$cday +%Y-%m-%d)    >> $datelist
echo $cyear                            >> $datelist
echo $(date --date=$cday +%m)          >> $datelist
echo $(date --date=$cday +%d)          >> $datelist
echo $jday_today                       >> $datelist
cp $datelist $logdir

if [ $print_info -eq 1 ]
then
   echo Datelist:
   echo $datelist
   echo Contents:
   cat $datelist
   echo " "
fi
#=======================================================================


##############################################################
if [ $manual -eq 0 ]
then
   # Checks before run:
   # 2. check if forecast has already run
   if [ -f $rundir/final_output/SWARP*.nc ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "Ice-only FC has already run - stopping"
      fi
      exit
   fi

fi
###################################################################


###################################################################
# RUNNING TOPAZ_GET_RESTART
# get latest restart file from reanalysis
echo "Launching topaz_get_restart_reanalysis @ $(date)"   > $log
$FCcommon/topaz_get_restart_reanalysis.sh                   $THIS_SRC

#=======================================================================
# DATE INFO - LAST RESTART (REANALYSIS) 
cd $rundir  # just in case we've changed dir in script
out_restart=info/last_restart.txt
rname=$(cat $out_restart)
rgen=${rname:0:3}    # eg TP4
ryear=${rname:10:4}  # year of restart file
rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)
rdate=`date --date="$ryear-01-01 +${rday}days" +%Y%m%d`
r_ndays=`date --date="$ryear-12-31" +%j`
#=======================================================================

restart_OPT=$TP4restart_OPT # temporary variable
if [ $restart_OPT -eq 2 ]
then
   # ======================================================================
   # DATE INFO - YESTERDAY 
   # - need to dump restart for yesterday
   dump_day=$(date --date="$cday -1day" +%Y%m%d)
   dump_year=$(date --date=$dump_day +%Y)
   dump_day_j=$(date --date=$dump_day +%j)
   dump_day_j0=$((dump_day_j-1))
   # ======================================================================

   if [ $ryear -eq $dump_year ]
   then
      if [ $dump_day_j0 -eq $rday ]
      then
         # yesterday is also day of restart from reanalysis,
         # so don't need to dump extra restart
         restart_OPT=1
      else
         dump_day_j1=$(printf '%3.3d' $dump_day_j0) # 3 digits
      fi
   else
      dump_day_j1=$((r_ndays+dump_day_j0))
   fi
fi
###################################################################


# =========================================================================
# DATE INFO - LAST DAY OF FORECAST
# $final_day=julian day relative to $ryear (NB 3 digits)
fin_day=`date --date="$cday +${fc_days}days" "+%Y%m%d"`
fin_year=$(date --date=$fin_day +%Y)
fin_day_j=$(date --date=$fin_day +%j)
fin_day_j0=$((fin_day_j-1))
if [ $ryear -eq $fin_year ]
then
   final_day=$(printf '%3.3d' $fin_day_j0) # 3 digits
else
   final_day=$((r_ndays+fin_day_j0))
fi
# =========================================================================


#################################################################
# MAKE INFILE 
echo "Launching make_infile4forecast @ $(date)"          >> $log

# print to screen
echo "Restart date   : $rdate"
echo "Current date   : $cday"
echo "Final date     : $fin_day"
echo "Final hour     : $final_hour"
echo "Restart files of $rname"
echo "Forecast final day ${fin_year}_$(printf '%3.3d' $fin_day_j0)_$final_hour (${ryear}_${final_day}_$final_hour)"

if [ $restart_OPT -eq 1 ]
then
   echo "$FCcommon/make_infile_2days.sh $THIS_SRC $rgen $ryear $rday $final_day $final_hour"
   $FCcommon/make_infile_2days.sh       $THIS_SRC $rgen $ryear $rday $final_day $final_hour
else
   echo "$FCcommon/make_infile_3days.sh $THIS_SRC $rgen $ryear $rday $dump_day_j1 $final_day $final_hour"
   $FCcommon/make_infile_3days.sh       $THIS_SRC $rgen $ryear $rday $dump_day_j1 $final_day $final_hour
fi

xdir=$TP4_REALTIME/expt_01.$Xno
infile=$xdir/infile.in
if [ $print_info -eq 1 ]
then
   echo ""
   echo "Infile: $infile"
   echo ""
   cat $infile
   echo ""
fi
cp $infile $rundir/info
#################################################################


#################################################################
# Get other inputs
echo "Launching pbsjob @ $(date)"                  >> $log
cd $xdir

# want to save archive files (TP4archv*.[ab]) every 3 hours
cp $THISFC/inputs/blkdat.input.archv_3h  blkdat.input
cp $THISFC/inputs/pbsjob.sh              .
cp $THISFC/inputs/preprocess.sh          .

#choose variables to extract (-DARCHIVE_SELECT)
cp $THISFC/inputs/archv.extract data
#################################################################


#################################################################
# clean data directory before run
if [ -f "./data/TP4DAILY*" ]
then
   rm data/TP4DAILY*
fi
if [ -f "./data/TP4archv*" ]
then
   rm data/TP4archv*
fi

# clean log file - else mpijob.out gets too big
rm -f log/*

if [ ! $test_pre -eq 1 ]
then
   # launch job
   $qsub pbsjob.sh
fi
#################################################################
