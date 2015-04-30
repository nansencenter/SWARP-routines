#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================

# EMAIL 
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
# ==============================================================================
email=$(cat $address)
# ==============================================================================

# TIMEOUT OF THE SCRIPT
# We want this script to run everyhour except from 00 to 04 in the morning since
# wam_nsea files are uploaded around 3.04 am

# THE PREVIOUS ASSUMPTION IS BASED ON A VISUAL CHECK MADE ON DATE 22/04/2015
# NB IF YOU SUSPECT THAT THE UPLOAD DATE HAS CHANGED PLEASE MODIFY THE FOLLOWING

# ==============================================================================

timeout_start=00
timeout_end=04
timeout_warning=08

# ==============================================================================

time_now=$(date +%H)

if  [ $time_now -gt $timeout_start ] && [ $time_now -lt $timeout_end ]
then
   exit
fi

SHOST=/Wave/METNO-WAM-NORDIC/
LDIR=/work/shared/nersc/msc/WAMNSEA

tday=$(date +%Y%m%d)

# check year today, tomorrow, the day after tomorrow
year=$(date +%Y)
dayname=$(date +%A)
yearplus1=$(expr $year + 1)
year1d=$(date --date="1 day" +"%Y")
year2d=$(date --date="2 day" +"%Y")

# creating (if needed) all the storing directories
mkdir -p $LDIR/${year}/analysis
mkdir -p $LDIR/${year}/forecasts
mkdir -p $LDIR/${yearplus1}/analysis
mkdir -p $LDIR/${yearplus1}/forecasts

# moving to the working folder
cd $LDIR/$year

log=$LDIR/$year/wamn_log.txt
rm $log
touch $log

echo "** Update WAMNSEA wave data from myocean.met.no **"   >> $log
echo " Today is $tday"                                      >> $log
echo " Check if year change next 2 days: "                  >> $log
echo " Today year is:.........$year"                        >> $log
echo " Next year is:..........$yearplus1"                   >> $log
echo " Today +1 day year is:..$year1d"                      >> $log
echo " Today +2 day year is:..$year2d"                      >> $log
echo ""                                                     >> $log

# Downloading EVERYTHING from myocean
# Checking the presence of the latest forecast
# if present continue with the script
# else send a warning_log


  cat > ncftp.in<<EOF
open myocean
cd $SHOST
get wam_nsea.an.*.nc
get wam_nsea.fc.*.nc
set confirm-close no
bye
EOF
   ncftp < ncftp.in
   rm ncftp.in

# Moving every download in the proper dir
mv $LDIR/$year/wam_nsea.an.* $LDIR/$year/analysis/
mv $LDIR/$year/wam_nsea.fc.* $LDIR/$year/forecasts/
   
if  ! [ -f "$LDIR/$year/forecasts/wam_nsea.fc.$tday.nc" ]
then
    if [ $time_now -gt $timeout_warning ]                # We want to give the script time till 0900 before sending a warning mail
    then
        echo " The file wam_nsea.fc.$tday.nc doesnt exist"                         >> $log
        echo " Please check either myocean ftp server status "                     >> $log
        echo " IF server is online and updated check:   "                          >> $log
        echo " ~/SWARP-routines/forecast_scripts/wamnsea_update.sh"                >> $log 
        mail -s "WAMNSEA forecast not found for $tday" $email < $log
     else
        exit
else
  echo "The file wam_nsea.an.$tday.nc exist - continue"                      >> $log
fi

#######################################################################
# OLD SCRIPT (NB USING SLEEP COMMAND WITH CRONTAB MAY CAUSE FAILURE)
#######################################################################

# Loop, try download WAMNSEA file, if not exist try again after 30 min
## after 4.5h exit and write Error message
#count=0
#while  ! [ -f wam_nsea.an.${tday}.nc ]
#do
#
#  cat > ncftp.in<<EOF
#open myocean
#cd $SHOST
#get wam_nsea.an.*.nc
#get wam_nsea.fc.*.nc
#set confirm-close no
#bye
#EOF
#   ncftp < ncftp.in
#   rm ncftp.in
#
# if  ! [ -f wam_nsea.an.${tday}.nc ]
#  then
#   count=`expr $count + 1`
#   echo " The file wam_nsea.an.${tday}.nc doesnt exist"
#   echo ""
#    if [  $count -ge 9 ]
#     then 
#      echo " *************************************************** "
#      echo " ***** Download WAMNSEA 10km data fails - exit ***** "
#      echo " *************************************************** "
#      exit
#    fi
#   echo " $count :count - We wait 30 min, try to download agin"
#   echo ""
#   sleep 1800
# else
#   echo "The file wam_nsea.an.${tday}.nc exist - continue"
#   echo ""
# fi
#
#done

######################################################################

# load cdo on hexagon with direct path
# module load cdo
# /work/apps/cdo/1.5.4-cray/bin/cdo

# check year, save in correct year file

# Use cdo functions to add last +1 day (.an. file)
#  to "keep" file (wam_nsea_2012.nc)
# Remove and recreate forecast file (wam_nsea_fc_${year}.nc)
#  that also include the +1 and +2 forecast days

echo " Add data for this day (.an. file) at end of year file"     >> $log
echo " Todays inputfile: wam_nsea.an.${tday}.nc"                  >> $log
echo " Cat into year file: wam_nsea_${year}.nc"                   >> $log
echo " Then rm wam_nsea_fc_${year}.nc and "                       >> $log
echo " cdo copy wam_nsea_${year}.nc wam_nsea.fc.${tday}.nc "      >> $log
echo " into new version of wam_nsea_fc_${year}.nc"                >> $log
echo ""                                                           >> $log

# Check if .an. file exist, if not exit
if [ -f "$LDIR/$year/analysis/wam_nsea.an.${tday}.nc" ]
 then
  /work/apps/cdo/1.5.4-cray/bin/cdo cat analysis/wam_nsea.an.${tday}.nc $LDIR/${year}/wam_nsea_${year}.nc
else
 echo " PROBLEM! wam_nsea.an.${tday}.nc doesn't exist"                              >> $log
 echo " check ~/SWARP-routines/forecast_scripts/wamnsea_update.sh"                  >> $log
 mail -s "WAMNSEA forecast not found for $tday" $email < $log
 exit
fi

echo "Remove forecast year file wam_nsea_fc_${year}.nc if exist"  >> $log
echo ""                                                           >> $log
if [ -f "$LDIR/$year/wam_nsea_fc_${year}.nc" ]
 then
  rm $LDIR/$year/wam_nsea_fc_${year}.nc
fi

echo "Check for year today, today +1, today +2"                   >> $log
if [ $year -eq $year1d -a $year -eq $year2d ]
 then
  echo "Same year all 3 days"                                     >> $log # merge fc in the file of the year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy wam_nsea_${year}.nc $LDIR/$year/forecasts/wam_nsea.fc.${tday}.nc wam_nsea_fc_${year}.nc
elif [ $year -eq $year1d -a $yearplus1 -eq $year2d ]
 then
  echo "Same year today and today +1 - new year today +2"         >> $log # split the fc file, day+1 in this year, day+2 in the next year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy wam_nsea_${year}.nc wam_nsea_fc_${year}.nc
  /work/apps/cdo/1.5.4-cray/bin/cdo splityear forecasts/wam_nsea.fc.${tday}.nc wam_nsea_fc_split_
  /work/apps/cdo/1.5.4-cray/bin/cdo cat wam_nsea_fc_split_${year}.nc wam_nsea_fc_${year}.nc
  /work/apps/cdo/1.5.4-cray/bin/cdo cat wam_nsea_fc_split_${yearplus1}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  rm wam_nsea_fc_split_*
  echo ""                                                         >> $log
  echo "NEW YEAR IS A MAGICAL TIME FOR BUGS! BE CAREFUL!"  >> $log
  mail -s "WAMNSEA UPDATE - 2 day to new year" $email < $log
elif [ $yearplus1 -eq $year1d -a  $yearplus1 -eq $year2d ]
 then
  echo "New year today +1 and today +2"                           >> $log # the whole fc file is in the next year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy forecasts/wam_nsea.fc.${tday}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  echo "Happy new year user - remember to check if the next upload is correct"   >> $log
  echo "NEW YEAR IS A SUPER MAGICAL TIME FOR BUGS! BE EXTRA CAREFUL!"            >> $log
  mail -s "WAMNSEA UPDATE - 1 day to new year" $email < $log
else
  echo " !!WARNING!! "                                                              >> $log
  echo " Check ~/SWARP-routines/forecast_scripts/wamnsea_update.sh"                 >> $log
  mail -s "WAMNSEA merging didn't work for $tday" $email < $log  
fi

