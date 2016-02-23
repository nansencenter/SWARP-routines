#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ModelSetups/common/wamnsea_data
fget="$THISFC/wamnsea_download.sh"

# ==============================================================================
# EMAIL 
email=$(cat $FCEmail)
# ==============================================================================

# TIMEOUT OF THE SCRIPT
# We want this script to run everyhour except from 00 to 04 in the morning since
# wam_nsea files are uploaded around 4am

# THE PREVIOUS ASSUMPTION IS BASED ON A VISUAL CHECK MADE ON DATE 22/04/2015
# NB IF YOU SUSPECT THAT THE UPLOAD DATE HAS CHANGED PLEASE MODIFY THE FOLLOWING

# ==============================================================================
timeout_warning=08
# ==============================================================================

time_now=$(date +%H)
tday=$(date +%Y%m%d)
if  [ $time_now -le 4 ]
then
   exit
fi

# check year today, tomorrow, the day after tomorrow
year=$(date +%Y)
dayname=$(date +%A)

# creating (if needed) all the storing directories
mkdir -p $wamnsea/${year}/analysis
mkdir -p $wamnsea/${year}/forecasts

# moving to the working folder
cd $wamnsea/$year

log=$wamnsea/$year/wamn_log.txt
rm $log
touch $log

echo "** Update WAMNSEA wave data from myocean.met.no **"   >> $log
echo " Today is $tday"                                      >> $log
echo " Today year is:.........$year"                        >> $log
echo ""                                                     >> $log

# Loop over previous 30 days
# - check if (an) wave file is there already
# - if not download it
# If latest fc isn't present continue with the script
# else send a warning_log
AN_DOWNLOADS=0
check_days=7
for n in `seq 0 $check_days`
do
   ddate=$(date --date="-$n days" +%Y%m%d) # download date
   dyear=${ddate:0:4}
   ddir1=$wamnsea/${dyear}/analysis
   cd $ddir1
   dfil1=wam_nsea.an.$ddate.nc

   # an file
   if [ ! -f $ddir1/$dfil1 ]
   then
      $fget $ddate an
      AN_DOWNLOADS=$((AN_DOWNLOADS+1))
   # else
   #    echo $dfil1 exists already
   fi
done
#################################################

#################################################
# old fc files are deleted every day
# - just look for today's
ddate=`date +%Y%m%d`
dyear=${ddate:0:4}
ddir2=$wamnsea/${dyear}/forecasts
dfil2=wam_nsea.fc.$ddate.nc
cd $ddir2

FC_DOWNLOADS=0
if [ ! -f $ddir2/$dfil2 ]
then
   $fget $ddate fc
   FC_DOWNLOADS=$((FC_DOWNLOADS+1))
# else
#    echo $dfil2 exists already
fi
#################################################

if [ -f $ddir2/$dfil2 ]
then
   #update the link to the latest forecast file
   cd $wamnsea
   fclat=wam_nsea.fc.latest.nc
   rm -f $fclat
   ln -s $ddir2/$dfil2 $fclat
fi

if  [ ! -f "$wamnsea/$year/forecasts/wam_nsea.fc.$tday.nc" ]
then
   if [ $time_now -gt $timeout_warning ]
   then
      # We want to give the script time till 0900 before sending a warning mail
      echo " The file wam_nsea.fc.$tday.nc doesnt exist"             >> $log
      echo " Please check either myocean ftp server status "         >> $log
      echo " IF server is online and updated check:   "              >> $log
      echo " ~/SWARP-routines/forecast_scripts/wamnsea_update.sh"    >> $log 
      mail -s "WAMNSEA forecast not found for $tday" $email < $log
   else
      exit
   fi
else
   echo "The file wam_nsea.an.$tday.nc exists - continue"   >> $log
fi
