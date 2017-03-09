#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
HERE="$FORECAST/ModelSetups/common/ww3arctic_data"
fget="$HERE/ww3_arctic_download.sh"
fsort="$HERE/ww3_arctic_sort.sh"
fget2="$HERE/ww3_arctic_get_best_est.sh"
ww3a=$wmsc/WAVES_INPUT/WW3_ARCTIC
cycle=2 # get 05:20 update from ifremer server
print_info=0

# ==============================================================================
# EMAIL 
email=$(cat $FCemail)
# ==============================================================================


# ==============================================================================
# TIMEOUT OF THE SCRIPT
# We want this script to run everyhour except from 00 to 05 in the morning since
# ww3arctic files are uploaded around 5.20am
# - after 8am if files are still not there, send email
timeout_warning=08
# ==============================================================================


time_now=$(date +%H)
tday=$(date +%Y%m%d)
yday=`date --date="yesterday" "+%Y%m%d"`
yyear=${yday:0:4}

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

# only today's FC is kept on the ftp server
# (previously  this was indexed according to yesterday's date,
#  since the run began then)
ddir1=$ww3a/${yday:0:4}/originals/${yday}00
dfil1=SWARP_WW3_ARCTIC-12K_$yday.nc

if [ $print_info -eq 1 ]
then
   echo "Checking $yday"
fi

# 2d fields
ddir2=$ww3a/${yyear}/forecast
dfil2=SWARP_WW3_ARCTIC-12K_${yday}.fc.nc

# full spectrum
ddir3=$ww3a/${yyear}_ef/forecast
dfil3=SWARP_WW3_ARCTIC-12K_${yday}_ef.fc.nc
DLS=0

if [ ! -f $ddir1/$dfil1 ]
then
   # download, sort and convert files
   echo "** Update WW3_ARCTIC wave data from ftp.ifremer.fr **"   >  $log
   echo " Today is $tday $(date "+%H:%M:%S")"                     >> $log
   echo ""                                                        >> $log

   echo "Downloading WW3 Arctic files for $yday"  >> $log
   $fget $yday

   echo "Sorting WW3 Arctic files for $ddate"  >> $log
   $fsort   $ddate  1
   DLS=1
elif [ $print_info -eq 1 ]
then
   echo "WW3 Arctic files present for $ddate"      >> $log
fi

if [ ! -f $ddir1/$dfil1 ]
then
   echo "> WW3 Arctic files still not present for $ddate"  >> $log
   echo " "                                                >> $log
fi

if [ $DLS -eq 1 ]
then
   # if doing download, also update best est for 2 days ago:
   ddate=`date --date="$tday -2days" "+%Y%m%d"`
   $fget2 $ddate
fi


for n in `seq -30 -3`
do
   # also get best est's for last few days if not present
   ddate=`date --date="$tday ${n}days" "+%Y%m%d"`
   ddir=$ww3a/${ddate:0:4}/analysis_m1
   dfil=SWARP_WW3_ARCTIC-12K_${ddate}.nc
   if [ ! -f $ddir/$dfil ]
   then
      $fget2 $ddate
   fi
done


# ===============================================
# WARNINGS
# ===============================================
if [ $DLS -gt 0 ]
then
   echo " "                         >> $log
   echo "Number of downloads: $DLS" >> $log
fi

echo $ddir2/$dfil2
if [ -f $ddir2/$dfil2 ]
then

   fclat=$ww3a/SWARP_WW3_ARCTIC-12K.fc.latest.nc
   fctarg=$ddir2/$dfil2
   echo $fctarg

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
