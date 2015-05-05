#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh
#  

# FORECAST DAYS
# ================================================================================================
fc_days=2.5            
# ================================================================================================

# DIRECTORIES
SWARP_ROUTINES=$HOME/GITHUB-REPOSITORIES/SWARP-routines
TP4_REALTIME=/work/timill/RealTime_Models/TP4a0.12

# CREATING THE LOG
logdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/logs
mkdir -p $logdir
log=$logdir/run_forecast_wav_log.txt
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

touch $datelist
echo $(date +%Y%m%d)       >> $datelist
echo $(date +%Y-%m-%d)     >> $datelist
echo $(date +%Y)           >> $datelist
echo $(date +%m)           >> $datelist
echo $(date +%d)           >> $datelist
echo $(date +%j)           >> $datelist

rundir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $datelist | sed '1!d') # where the last_restart.txt will end up
cd $rundir
mkdir -p ./info
cp $datelist ./info


# RUNNING TOPAZ_GET_RESTART
#TODO get today's restart - should be output by ice_only forecast
# echo "Launching topaz_get_restart @ $date"                     >> $log
# $SWARP_ROUTINES/forecast_scripts/topaz_get_restart.sh          # get latest restart file
# cd $rundir                                                     # just in case we've changed dir in script
# 
# # GETTING INFOS FROM LAST_RESTART
# out_restart=info/last_restart.txt
# 
# cyear=$(cat $datelist | sed '3!d')
# rname=$(cat $out_restart)
# 
# rgen=${rname:0:3}                                              # eg TP4
# ryear=${rname:10:4}                                            # year of restart file
# rday=${rname:15:3}                                             # julian day of restart file (1 Jan = 0)

#################################################################

# MAKE INFILE 
# TODO SHOULD BE READY, W8TING FOR THE RESTARTS 
# echo "Launching make_infile4forecast @ $date"                  >> $log
# jday0=$(cat $datelist | sed '6!d')                             # julian day of today (1=1st Jan => need to change)
# jday_today0=$(expr $jday0 - 1)                                 # julian day of TOPAZ (0=1st Jan)
# jday_today=$(printf '%3.3d' $jday_today0)
# final_day=$(expr $jday_today + $fc_days)
# 
# # print to screen - work out if last day of forecast is in a different year to current year
# ndays=$(date --date="${cyear}-12-31" +%j)                 # days in current year
# if [ $final_day -gt $((ndays-1)) ]
# then
# 	fc_final_day=$(expr $final_day - $ndays)
# 	fc_year=$(expr $cyear + 1)
# else
# 	fc_year=$cyear
# fi
# echo "Restart files of ${cyear}_${jday_today}"                  >> $log
# echo "Forecast final day ${fc_year}_${fc_final_day}"            >> $log
# 
# $SWARP_ROUTINES/forecast_scripts/make_infile4forecast.sh $rgen $ryear $jday_today $final_day
# xdir=$TP4_REALTIME/expt_01.1
# infile=$xdir/infile.in
# if [ $rday -eq $jday_today ]
# then
#    # delete "today" line
#    echo "restart day is from today"
#    echo "- editing infile.in"
#    sed '17d' $infile >> infile.in.replace
#    mv infile.in.replace $infile
# fi
#################################################################

# Launch job
echo "Launching pbsjob @ $(date)"                  >> $log
cd $xdir

# want to save archive files (TP4archv*.[ab]) every 3 hours
cp $SWARP_ROUTINES/forecast_scripts/inputs/wavesice/blkdat.input .
cp $SWARP_ROUTINES/forecast_scripts/inputs/wavesice/pbsjob.sh pbsjob.sh

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
qsub=/opt/torque/2.5.13pre-up/bin/qsub #get full path from which qsub
$qsub pbsjob.sh
#################################################################
   
echo "pbsjob done @ $(date)"                  >> $log

