#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh

# get variables
# ($SWARP_ROUTINES defined in crontab or .bash_profile)
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_ww3arctic
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC
xdir=$TP4_REALTIME/expt_01.$Xno
restart_OPT=$TP4restart_OPT

test_pre=1
print_info=1 # print info to screen (or email in crontab)
manual=0
if [ $# -eq 1 ]
then
   echo ""
   echo "Running manually (over-riding checks)"
   echo ""
   manual=1
fi

if [ $print_info -eq 1 ]
then
   echo "Expt no     : $Xno"
   echo "Results to  : $THISFC2"
   echo "Email       : $FCemail"
   echo ""
fi

# ================================================================================================
# FORECAST DAYS
fc_days=2
final_hour="12"   # finishes at this hour on the final day
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
   # test last Tuesday (restart and dump day are the same
   # for $TP4restart_OPT=2)
   cday=$(date --date="last Tuesday" +%Y%m%d)
   echo " "
   echo "******************************************************"
   echo "WARNING: testing artificial start date:" $cday
   echo "******************************************************"
   echo " "
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

if [ $restart_OPT -eq 2 ]
then
   # ======================================================================
   # DATE INFO - YESTERDAY 
   # - need to get restart for yesterday
   rgen="TP4"
   rdate=$(date --date="$cday -1day" +%Y%m%d)
   ryear=$(date --date=$rdate +%Y)
   rwday=$(date --date=$rdate +%A)
   rday_j=$(date --date=$rdate +%j)
   rday=$((rday_j-1))
   rdir=$TP4_REALTIME/expt_01.1/data
   rname=${rgen}restart${ryear}_${rday}_00
   if [ $rwday == 'Monday' ]
   then
      # yesterday's restart file is topaz reanalysis file
      restart_OPT=1
   fi
   # ======================================================================
fi
###################################################################


##############################################################
if [ $manual -eq 0 ]
then
   # Checks before run:
   # 1. check if forecast has already run
   if [ -f $rundir/final_output/SWARP*.nc ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "Waves-ice (WW3) FC has already run - stopping"
      fi
      exit
   fi

   # 2. check if forecast is already running
   msg=`$qstat | grep TP4x01${Xno}fc`
   if [ ${#msg} -ne 0 ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "Waves-ice (WW3) FC is already running - stopping"
         echo "pbs job message:"
         echo $msg
      fi
      exit
   fi

   if [ $restart_OPT -eq 2 ]
   then
      # 3. check if ice-only forecast is still running
      if [ ! -f $rundir/../ice_only/final_output/SWARP*.nc ]
      then
         if [ $print_info -eq 1 ]
         then
            echo "Ice-only FC has not yet run - stopping"
         fi
         exit
      fi
   fi
fi
###################################################################


#################################################################
if [ $restart_OPT -eq 1 ]
then
   # get reanalysis restart
   echo "Launching topaz_get_restart_reanalysis.sh @ $(date)"   > $log
   $FCcommon/topaz_get_restart_reanalysis.sh $THIS_SRC            # get latest restart file

   # GETTING INFO FROM LAST_RESTART
   cd $rundir  # just in case we've changed dir in script
   out_restart=info/last_restart.txt
   rname=$(cat $out_restart)
   rgen=${rname:0:3}    # eg TP4
   ryear=${rname:10:4}  # year of restart file
   rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)
   rdate=`date --date="$ryear-01-01 +${rday}days" +%Y%m%d`
else
   afil=$rdir/$rname.a
   bfil=$rdir/$rname.b
   ufil=$rdir/${rname}ICE.uf
   if [ -f $afil ] && [ -f $bfil ] && [ -f $ufil ]
   then
      echo "cp $afil $xdir/data" > $log
      echo "cp $bfil $xdir/data" > $log
      echo "cp $ufil $xdir/data" > $log
      cp $afil $xdir/data
      cp $bfil $xdir/data
      cp $ufil $xdir/data
      #
      out_restart=$rundir/info/last_restart.txt
      echo $rname > $out_restart
      cp $out_restart $logdir
   else
      echo "Yesterday's restart files from ice-only run not present:"
      echo "$afil"
      echo "$bfil"
      echo "$ufil"
      exit
   fi
fi
#################################################################


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
echo "Launching make_infile_2days.sh @ $(date)"                  >> $log

# print to screen
echo "Restart date   : $rdate"
echo "Current date   : $cday"
echo "Final date     : $fin_day"
echo "Final hour     : $final_hour"
echo "Restart files of $rname"
echo "Forecast final day ${fin_year}_$(printf '%3.3d' $fin_day_j0)_$final_hour (${ryear}_${final_day}_$final_hour)"

echo "$FCcommon/make_infile_2days.sh $THIS_SRC $rgen $ryear $rday $final_day $final_hour"
$FCcommon/make_infile_2days.sh       $THIS_SRC $rgen $ryear $rday $final_day $final_hour

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

# Don't want to save archive files (TP4archv*.[ab])
cp $THISFC/inputs/blkdat.input   $xdir
cp $THISFC/inputs/pbsjob.sh      $xdir
cp $THISFC/inputs/preprocess.sh  $xdir

#choose variables to extract (-DARCHIVE_SELECT)
# cp $THISFC/inputs/archv.extract $xdir/data

if [ $print_info -eq 1 ]
then
   echo "cp $THISFC/inputs/blkdat.input   $xdir"
   echo "cp $THISFC/inputs/pbsjob.sh      $xdir"
   echo "cp $THISFC/inputs/preprocess.sh  $xdir"
   # echo "cp $THISFC/inputs/archv.extract $xdir/data"
fi
#################################################################


#################################################################
# clean data directory before run
rm -f $xdir/data/TP4DAILY*
rm -f $xdir/data/TP4archv*

# clean log file - else mpijob.out gets too big
rm -f $xdir/log/*

if [ ! $test_pre -eq 1 ]
then
   # launch job
   echo "Launching pbsjob @ $(date)"   >> $log
   cd $xdir
   $qsub pbsjob.sh
fi
#################################################################
