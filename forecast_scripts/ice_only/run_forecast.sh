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
print_info=0 # print info to screen (or email in crontab)
manual=0
if [ $# -eq 1 ]
then
   echo ""
   echo "Running manually (no checks)"
   echo ""
   manual=1
fi

# FORECAST DAYS
# ================================================================================================
# fc_days=5
fc_days=6 # to be consistent with WW3 forecast            
# thour=`date "+%H"`
# if [ $thour -le 1 ]
# then
#    # too early (winds not there)
#    exit
# fi
# ================================================================================================


# CREATING THE LOG
logdir=$THISFC/logs
mkdir -p $logdir
log=$logdir/run_forecast_log.txt
if [ -f "$log" ]
then
   rm $log
fi

cday=$(date +%Y%m%d)
cyear=$(date +%Y)
jday0=$(date +%j)
jday_today0=$(expr $jday0 - 1)            # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0) # 3 digits
rundir=$THISFC2/$cday                     # where the results will end up

mkdir -p $rundir
mkdir -p $rundir/info

# CREATING THE DATELIST
# NB don't make one in $FORECAST since both
# ice-only and waves-ice FC's use it
# - if they run at the same time it gets messed up
datelist=$rundir/info/datelist.txt
echo $cday              >  $datelist
echo $(date +%Y-%m-%d)  >> $datelist
echo $cyear             >> $datelist
echo $(date +%m)        >> $datelist
echo $(date +%d)        >> $datelist
echo $jday_today        >> $datelist
cp $datelist $logdir

if [ $print_info -eq 1 ]
then
   echo Datelist:
   echo $datelist
   echo Contents:
   cat $datelist
   echo " "
fi

##############################################################
if [ $manual -eq 0 ]
then
   # Checks before run:
   # 1. check if forecast has already run
   if [ -f $rundir/final_output/SWARP*.nc ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "Ice-only FC has already run - stopping"
      fi
      exit
   fi

   # 2. check if forecast is already running
   msg=`$qstat | grep TP4x011fc`
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
fi
###################################################################

# RUNNING TOPAZ_GET_RESTART
echo "Launching topaz_get_restart @ $(date)"   > $log
$THISFC/pre/topaz_get_restart.sh $THIS_SRC         # get latest restart file

# GETTING INFO FROM LAST_RESTART
cd $rundir  # just in case we've changed dir in script
out_restart=info/last_restart.txt
rname=$(cat $out_restart)
rgen=${rname:0:3}    # eg TP4
ryear=${rname:10:4}  # year of restart file
rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)

#################################################################

# MAKE INFILE 
echo "Launching make_infile4forecast @ $(date)"                  >> $log

# if last restart was in different year to this year:
if [ $print_info -eq 1 ]
then
   echo current year: $cyear
   echo restart year: $ryear
   echo restart name: $rname
fi

if [ $cyear -ne $ryear ]
then
   jday_today=$(expr $jday_today + $rday + 1) # integer
   jday_today=$(printf '%3.3d' $jday_today)   # 3 digits (compare to rday)
fi

final_day=$(expr $jday_today + $fc_days)

# print to screen - work out if last day of forecast is in a different year to current year
ndays=$(date --date="${cyear}-12-31" +%j)                 # days in current year
if [ $final_day -gt $((ndays-1)) ]
then
   fc_final_day=`expr $final_day - $ndays`
   fc_year=`expr $cyear + 1`
else
   fc_year=$cyear
fi
echo "Restart files of ${ryear}_$rday"
echo "Forecast final day ${fc_year}_${final_day}"

# $THISFC/pre/make_infile4forecast.sh $rgen $ryear $rday $jday_today $final_day
echo "$THISFC/pre/make_infile4forecast.sh $rgen $ryear $rday $final_day"
$THISFC/pre/make_infile4forecast.sh $rgen $ryear $rday $final_day

xdir=$TP4_REALTIME/expt_01.$Xno
infile=$xdir/infile.in
if [ $print_info -eq 1 ]
then
   cat $infile
fi
# if [ $rday -eq $jday_today ]
# then
#    # delete "today" line
#    echo "restart day is from today"
#    echo "- editing infile.in"
#    sed '17d' $infile >> infile.in.replace
#    mv infile.in.replace $infile
# fi
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
