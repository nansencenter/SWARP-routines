#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh
#  

# FORECAST DAYS
# ================================================================================================
fc_days=5            
# thour=`date "+%H"`
# if [ $thour -le 1 ]
# then
#    # too early (winds not there)
#    exit
# fi
# ================================================================================================

# get variables
# ($SWARP_ROUTINES defined in crontab or .bash_profile)
source $SWARP_ROUTINES/source_files/hex_vars.src

# CREATING THE LOG
logdir=$SWARP_ROUTINES/forecast_scripts/logs
mkdir -p $logdir
log=$logdir/run_forecast_log.txt
if [ -f "$log" ]
then
   rm $log
fi
touch $log

# CREATING THE DATELIST 
datelist=$SWARP_ROUTINES/forecast_scripts/datelist.txt

if [ -f "$datelist" ]
then
   rm $datelist
fi

tday=$(date +%Y%m%d)
touch $datelist
echo $tday                 >> $datelist
echo $(date +%Y-%m-%d)     >> $datelist
echo $(date +%Y)           >> $datelist
echo $(date +%m)           >> $datelist
echo $(date +%d)           >> $datelist
jday0=$(date +%j)                             # julian day of today (1=1st Jan => need to change)
jday_today0=$(expr $jday0 - 1)                                 # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0)
echo $jday_today           >> $datelist

rundir=$TP4_REALTIME/../results/TP4a0.12/ice_only/work/$tday # where the last_restart.txt will end up
mkdir -p $rundir
cd $rundir
mkdir -p ./info
cp $datelist ./info

##############################################################
# Checks before run:

# 1. check if forecast has already run
if [ -f $rundir/final_product/SWARP*.nc ]
then
   # echo "Ice-only FC has already run - stopping"
   exit
fi

# 2. check if forecast is already running
msg=`$qstat | grep TP4x011fc`
if [ ${#msg} -ne 0 ]
then
   # echo "forecast is already running - stopping"
   # echo "pbs job message:"
   # echo $msg
   exit
fi
###################################################################

# RUNNING TOPAZ_GET_RESTART
echo "Launching topaz_get_restart @ $date"                     >> $log
$SWARP_ROUTINES/forecast_scripts/ice_only/topaz_get_restart.sh          # get latest restart file

# GETTING INFOS FROM LAST_RESTART
cd $rundir  # just in case we've changed dir in script
out_restart=info/last_restart.txt

year_today=$(cat $datelist | sed '3!d')
rname=$(cat $out_restart)

rgen=${rname:0:3}    # eg TP4
ryear=${rname:10:4}  # year of restart file
rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)

#################################################################

# MAKE INFILE 
echo "Launching make_infile4forecast @ $date"                  >> $log


# if last restart was in different year to this year:
if [ $year_today -ne $ryear ]
then
   jday_today=$(expr $jday_today + $rday + 1) # integer
   jday_today=$(printf '%3.3d' $jday_today)   # 3 digits (compare to rday)
fi

final_day=$(expr $jday_today + $fc_days)

# print to screen - work out if last day of forecast is in a different year to current year
ndays=$(date --date="${year_today}-12-31" +%j)                 # days in current year
if [ $final_day -gt $((ndays-1)) ]
then
   fc_final_day=`expr $final_day - $ndays`
   fc_year=`expr $year_today + 1`
else
   fc_year=$year_today
fi
echo "Restart files of ${ryear}_$rday"
echo "Forecast final day ${fc_year}_${final_day}"

$SWARP_ROUTINES/forecast_scripts/ice_only/make_infile4forecast.sh $rgen $ryear $rday $jday_today $final_day
xdir=$TP4_REALTIME/expt_01.1
infile=$xdir/infile.in
if [ $rday -eq $jday_today ]
then
   # delete "today" line
   echo "restart day is from today"
   echo "- editing infile.in"
   sed '17d' $infile >> infile.in.replace
   mv infile.in.replace $infile
fi
#################################################################

# Launch job
echo "Launching pbsjob @ $(date)"                  >> $log
cd $xdir

# want to save archive files (TP4archv*.[ab]) every 3 hours
cp $SWARP_ROUTINES/forecast_scripts/inputs/ice_only/blkdat.input.archv_3h blkdat.input
cp $SWARP_ROUTINES/forecast_scripts/inputs/ice_only/pbsjob.sh pbsjob.sh

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
rm log/*

# launch job
echo $qsub
$qsub pbsjob.sh
#################################################################

echo "pbsjob done @ $(date)"                  >> $log
