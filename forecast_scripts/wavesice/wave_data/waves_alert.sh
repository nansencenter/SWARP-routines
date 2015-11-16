#!/bin/bash
# Sending a mail with the content of check_wamnsea.py
# $1 = YYYYMMDD

source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -eq 0 ]
then
   echo Usage:
   echo "./waves_alert.sh [date]"
   echo "date in YYYYMMDD format"
   exit
fi

#######################################################################################################
# EMAIL ADRESS
email=$(cat $FCemail)
#######################################################################################################

wdir=/work/timill/RealTime_Models/check_wamnsea/$1
echo "waves_alert.sh called by check_wamnsea.py:"  >  tmp.txt
echo " "                                           >> tmp.txt
cat $wdir/lst/*.txt                                >> tmp.txt

if [ 1 -eq 0 ]
then
   # just email person from file
   mutt -s "WAM forecast for $1" -a $wdir/img/*.png -a $wdir/lst/*.txt -- $email < tmp.txt
else
   # add Giacomo for his interest
   mutt -s "WAM forecast for $1" -a $wdir/img/*.png -a $wdir/lst/*.txt -- $email,gcmdnt90@gmail.com < tmp.txt
fi
rm tmp.txt
