#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC="$FORECAST/ModelSetups/common/ww3arctic_data"
fget="$THISFC/ww3_arctic_download.sh"
fsort="$THISFC//ww3_arctic_sort.sh"
ww3a=$wmsc/WAVES_INPUT/WW3_ARCTIC
cycle=2 # get 05:20 update from ifremer server
print_info=0

# ==============================================================================
# EMAIL 
email=$(cat $FCemail)
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
yday=`date --date="yesterday" "+%Y%m%d"`

if  [ $time_now -lt 5 ]
then
   exit
fi

# check year today, tomorrow, the day after tomorrow
year=$(date +%Y)

# creating (if needed) all the storing directories
mkdir -p $ww3a/${year}

# moving to the working folder
cd $ww3a/$year
log=$ww3a/$year/ww3a_log.txt

# Loop over previous 30 days
# - check if eg wave file is there already
# - if not download all materials for that day
DLS=0 # no of downloads
dfirst=1
for n in `seq 1 3`
do
   ddate=$(date --date="-$n days" +%Y%m%d) # download date
   dyear=${ddate:0:4}
   ddir1=$ww3a/${dyear}/originals/${ddate}00
   dfil1=SWARP_WW3_ARCTIC-12K_$ddate.nc

   if [ $print_info -eq 1 ]
   then
      echo "Checking $ddate"
   fi

   if [ $n -eq 1 ]
   then
      ddir2=$ww3a/${dyear}/forecast
      dfil2=SWARP_WW3_ARCTIC-12K_${ddate}.fc.nc
      #
      ddir3=$ww3a/${dyear}_ef/forecast
      dfil3=SWARP_WW3_ARCTIC-12K_${ddate}_ef.fc.nc
   fi

   if [ ! -f $ddir1/$dfil1 ]
   then
      # download, sort and convert files
      if [ $dfirst -eq 1 ]
      then
         echo "** Update WW3_ARCTIC wave data from ftp.ifremer.fr **"   >  $log
         echo " Today is $tday $(date "+%H:%M:%S")"                     >> $log
         echo ""                                                        >> $log
         dfirst=0
      fi
      echo "Downloading WW3 Arctic files for $ddate"  >> $log
      $fget    $ddate $cycle

      echo "Sorting WW3 Arctic files for $ddate"  >> $log
      $fsort   $ddate $cycle
      DLS=$((DLS+1))
   # else
   #    echo "WW3 Arctic files present for $ddate"      >> $log
   fi

   if [ ! -f $ddir1/$dfil1 ]
   then
      echo "> WW3 Arctic files still not present for $ddate"  >> $log
      echo " "                                                >> $log
   fi
done
#################################################

# ===============================================
# WARNINGS
# ===============================================
if [ $DLS -gt 0 ]
then
   echo " "                         >> $log
   echo "Number of downloads: $DLS" >> $log
fi


if [ -f $ddir2/$dfil2 ]
then

   fclat=$ww3a/SWARP_WW3_ARCTIC-12K.fc.latest.nc
   fctarg=$ddir2/$dfil2

   if [ ! -L $fclat ]
   then
      echo ln -s $fctarg $fclat
      ln -s $fctarg $fclat
   fi
   linkloc=`readlink $fclat`

   if [ ! $linkloc == $fctarg ]
   then
      echo " "
      echo "Link location     : $linkloc"
      echo "Latest forecast   : $fctarg"
      echo "Setting new link..."
      echo ln -s $fctarg $fclat
      rm $fclat
      ln -s $fctarg $fclat
      linkloc=`readlink $fclat`
      echo "New link location : $linkloc"
      echo " "

   elif [ $print_info -eq 1 ]
   then

      echo " "
      echo "Link name         : $fclat"
      echo "Link location     : $linkloc"
      echo "Latest forecast   : $fctarg"
      echo "Link already points to latest forecast"
      echo " "
   fi

else

   if [ $time_now -gt $timeout_warning ]
   then
      # We want to give the script time till 0900 before sending a warning mail
      echo " The file $dfil2 doesn't exist"                       >> $log
      echo " Please check either ifremer ftp server status "      >> $log
      echo " IF server is online and updated check:   "           >> $log
      echo " `readlink -f $0`"                                    >> $log 
      echo " and $ww3a/$year"                                     >> $log 
      #
      mail -s "WW3 Arctic forecast not found for $tday" $email < $log
   fi

fi

if [ -f $ddir3/$dfil3 ]
then

   # update the link to the latest forecast file
   # - ef type
   fclat=$ww3a/SWARP_WW3_ARCTIC-12K_ef.fc.latest.nc
   fctarg=$ddir3/$dfil3

   if [ ! -L $fclat ]
   then
      echo ln -s $fctarg $fclat
      ln -s $fctarg $fclat
   fi
   linkloc=`readlink $fclat`


   if [ ! $linkloc == $fctarg ]
   then

      echo " "
      echo "Link location     : $linkloc"
      echo "Latest forecast   : $fctarg"
      echo "Setting new link..."
      echo ln -s $fctarg $fclat
      rm $fclat
      ln -s $fctarg $fclat
      linkloc=`readlink $fclat`
      echo "New link location : $linkloc"
      echo " "

   elif [ $print_info -eq 1 ]
   then

      echo " "
      echo "Link name         : $fclat"
      echo "Link location     : $linkloc"
      echo "Latest forecast   : $fctarg"
      echo "Link already points to latest forecast"
      echo " "

   fi

else

   if [ $time_now -gt $timeout_warning ]
   then
      # We want to give the script time till 0900 before sending a warning mail
      echo " The file $dfil3 doesn't exist"                    >> $log
      echo " Please check either ifremer ftp server status "   >> $log
      echo " IF server is online and updated check:   "        >> $log
      echo " `readlink -f $0`"                                 >> $log 
      echo " and $ww3a/${year}_ef"                             >> $log 
      #
      mail -s "WW3 Arctic forecast (Ef) not found for $tday" $email < $log
   fi

fi

# echo " "
# echo "ww3_arctic_update.sh finished"
# echo "log in $log"
# echo " "
# 
# echo " " >> $log
# cat $log
# 
# echo " "
# cat $ww3a/latest_download.txt
