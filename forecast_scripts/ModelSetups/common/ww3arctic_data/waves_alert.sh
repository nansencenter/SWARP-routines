#!/bin/bash
# Sending a mail with the content of check_wamnsea.py
# $1 = YYYYMMDD

source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -lt 2 ]
then
   echo Usage:
   echo "./waves_alert.sh [date] [region]"
   echo "date in YYYYMMDD format"
   echo "region = Beaufort or Barents"
   exit
fi
fcday=`date --date="$1 +5days" "+%Y%m%d"`

#######################################################################################################
# EMAIL ADRESS
email=$(cat $FCemail)
if [ 1 -eq 1 ]
then
   # add Giacomo for his interest
   email=$email,gcmdnt90@gmail.com
fi
#######################################################################################################

wdir=/work/timill/RealTime_Models/check_ww3arctic/$1
echo "waves_alert.sh called by check_ww3arctic.py:"   >  tmp.txt
echo " "                                              >> tmp.txt

# send email
if [  -f $wdir/lst/$fcday*$2*.txt ]
then
   cat $wdir/lst/$fcday*$2*.txt                                                                             >> tmp.txt
   mutt -s "WW3a forecast for $fcday - $2" -a $wdir/img/$fcday*$2.png -a $wdir/lst/$fcday*$2*.txt -- $email <  tmp.txt
else
   echo "No large waves close to ice"                                                                       >> tmp.txt
   mutt -s "WW3a forecast for $fcday - $2" -a $wdir/img/$fcday*$2.png -- $email                             <  tmp.txt
fi

rm tmp.txt
