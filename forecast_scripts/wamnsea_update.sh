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
# wam_nsea files are uploaded around 4am

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

# Loop over previous 30 days
# - check if (an) wave file is there already
# - if not download it
# If latest fc isn't present continue with the script
# else send a warning_log
AN_DOWNLOADS=0
for n in `seq 0 30`
do
   ddate=$(date --date="-$n days" +%Y%m%d) # download date
   dyear=${ddate:0:4}
   ddir1=$LDIR/${dyear}/analysis
   dfil1=wam_nsea.an.$ddate.nc

   # an file
   if [ ! -f $ddir1/$dfil1 ]
   then
      echo getting $dfil1
      cat > ncftp.in<<EOF
open myocean
cd $SHOST
get $dfil1
set confirm-close no
bye
EOF
      ncftp < ncftp.in
      rm ncftp.in
      mv $dfil1 analysis

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
ddir2=$LDIR/${dyear}/forecasts
dfil2=wam_nsea.fc.$ddate.nc

FC_DOWNLOADS=0
if [ ! -f $ddir2/$dfil2 ]
then
   echo getting $dfil2
   cat > ncftp.in<<EOF
open myocean
cd $SHOST
get $dfil2
set confirm-close no
bye
EOF
   ncftp < ncftp.in
   rm ncftp.in
   mv $dfil2 forecasts
   FC_DOWNLOADS=$((FC_DOWNLOADS+1))
# else
#    echo $dfil2 exists already
fi
#################################################


if  [ ! -f "$LDIR/$year/forecasts/wam_nsea.fc.$tday.nc" ]
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
     fi
else
  echo "The file wam_nsea.an.$tday.nc exists - continue"                      >> $log
fi


# check year, save in correct year file

# Use cdo functions to add last +1 day (.an. file)
#  to "keep" file (wam_nsea_${tyear}.nc)
# Remove and recreate forecast file (wam_nsea_fc_${year}.nc)
#  that also include the +1 and +2 forecast days

echo " Add data for this day (.an. file) at end of year file"     >> $log
echo " Todays inputfile: wam_nsea.an.${tday}.nc"                  >> $log
echo " Cat into year file: wam_nsea_${year}.nc"                   >> $log
echo " Then rm wam_nsea_fc_${year}.nc and "                       >> $log
echo " cdo copy wam_nsea_${year}.nc wam_nsea.fc.${tday}.nc "      >> $log
echo " into new version of wam_nsea_fc_${year}.nc"                >> $log
echo ""                                                           >> $log

if [ $AN_DOWNLOADS -gt 0 ]
then
   # TODO change this in case any other new an files appear
   $cdo cat analysis/wam_nsea.an.${tday}.nc $LDIR/${year}/wam_nsea_${year}.nc
fi

if [ $FC_DOWNLOADS -eq 0 ]
then
   # no more to do
   exit
fi

###############################################################################################
# if there's a new fc file then continue to combine them
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
  $cdo copy wam_nsea_${year}.nc $LDIR/$year/forecasts/wam_nsea.fc.${tday}.nc wam_nsea_fc_${year}.nc
elif [ $year -eq $year1d -a $yearplus1 -eq $year2d ]
 then
  echo "Same year today and today +1 - new year today +2"         >> $log # split the fc file, day+1 in this year, day+2 in the next year
  $cdo copy wam_nsea_${year}.nc wam_nsea_fc_${year}.nc
  $cdo splityear forecasts/wam_nsea.fc.${tday}.nc wam_nsea_fc_split_
  $cdo cat wam_nsea_fc_split_${year}.nc wam_nsea_fc_${year}.nc
  $cdo cat wam_nsea_fc_split_${yearplus1}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  rm wam_nsea_fc_split_*
  echo ""                                                         >> $log
  echo "NEW YEAR IS A MAGICAL TIME FOR BUGS! BE CAREFUL!"  >> $log
  mail -s "WAMNSEA UPDATE - 2 day to new year" $email < $log
elif [ $yearplus1 -eq $year1d -a  $yearplus1 -eq $year2d ]
 then
  echo "New year today +1 and today +2"                           >> $log # the whole fc file is in the next year
  $cdo copy forecasts/wam_nsea.fc.${tday}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  echo "Happy new year user - remember to check if the next upload is correct"   >> $log
  echo "NEW YEAR IS A SUPER MAGICAL TIME FOR BUGS! BE EXTRA CAREFUL!"            >> $log
  mail -s "WAMNSEA UPDATE - 1 day to new year" $email < $log
else
  echo " !!WARNING!! "                                                              >> $log
  echo " Check ~/SWARP-routines/forecast_scripts/wamnsea_update.sh"                 >> $log
  mail -s "WAMNSEA merging didn't work for $tday" $email < $log  
fi
###############################################################################################
