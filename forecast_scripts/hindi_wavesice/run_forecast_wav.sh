#!/bin/bash
# This script will: 
# 1) Create a hindi_datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh
#  

# FORECAST DAYS
# ================================================================================================

fc_days=2.5           

# INSIDE INFILE.MAL.OUTER THE HALF DAY IS ALREADY CALCULATED

fc_day_mi4f=$(printf '%.f' $(echo "$fc_days-0.5" | bc))

# ================================================================================================

# DIRECTORIES
SWARP_ROUTINES=$HOME/GITHUB-REPOSITORIES/SWARP-routines
fcdir=$SWARP_ROUTINES/forecast_scripts
TP4_REALTIME=/work/timill/RealTime_Models/TP4a0.12

# CREATING THE LOG
logdir=$fcdir/logs
mkdir -p $logdir
log=$logdir/run_forecast_wav_log.txt
if [ -f "$log" ]
then
   rm $log
fi
touch $log

# CREATING THE DATELIST 
hindi_datelist=$fcdir/hindi_datelist.txt

if [ -f "$hindi_datelist" ]
then
   rm $hindi_datelist
fi

# We will need to recreate the hindi_datelist based on the input given by the user
echo 'Insert the date (YYYYMMDD):   '
echo ''
read dadate
dayear=${dadate:0:4}
damonth=${dadate:4:2}
daday=${dadate:6:2}

touch $hindi_datelist
echo $dadate >> $hindi_datelist
echo $dayear-$damonth-$daday >> $hindi_datelist
echo $dayear >> $hindi_datelist
echo $damonth >> $hindi_datelist
echo $daday >> $hindi_datelist
jday0=$(date -d "$dadate" +%j)
jday_today0=$(expr $jday0 - 1)                                 # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0)
echo $jday_today                             >> $hindi_datelist
final_day=$(echo "$jday_today + $fc_days" | bc)    
final_day_mi4f=$(expr $jday_today + $fc_day_mi4f)
cyear=$(cat $hindi_datelist | sed '3!d')

rundir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$(cat $hindi_datelist | sed '1!d') # where the last_restart.txt will end up
mkdir -p $rundir
cd $rundir
mkdir -p ./info
cp $hindi_datelist ./info

cday=$(cat $hindi_datelist | sed '1!d')
if [ -f $rundir/$cday/final_product/* ]
then
   exit
else
   # RUNNING TOPAZ_GET_RESTART
   echo "Launching topaz_get_restart_wav @ $date"                     >> $log
   $fcdir/wavesice/topaz_get_restart_wav.sh                                # get latest restart file
   cd $rundir                                                     # just in case we've changed dir in script

   # GETTING INFOS FROM LAST_RESTART
   out_restart=$rundir/info/last_restart.txt

   rname=$(cat $out_restart)

   rgen=${rname:0:3}                                              # eg TP4
   ryear=${rname:10:4}                                            # year of restart file
   rday=${rname:15:3}                                             # julian day of restart file (1 Jan = 0)

   #################################################################

   # MAKE INFILE 
   echo "Launching make_infile4forecast_wav @ $date"                  >> $log



   # print to screen - work out if last day of forecast is in a different year to current year
   ndays=$(date --date="${cyear}-12-31" +%j)                 # days in current year
   if [ $final_day_mi4f -gt $(($ndays-1)) ]
   then
           fc_final_day=$(expr $final_day - $ndays)
           fc_year=$(expr $cyear + 1)
   else
           fc_year=$cyear
   fi
   echo "Restart files of ${cyear}_${jday_today}"                  >> $log
   echo "Forecast final day ${fc_year}_${fc_final_day}"            >> $log

   $SWARP_ROUTINES/forecast_scripts/wavesice/make_infile4forecast_wav.sh $rgen $ryear $rday $final_day_mi4f
   xdir=$TP4_REALTIME/expt_01.2
   infile=$xdir/infile.in

   ###############################################################

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
fi