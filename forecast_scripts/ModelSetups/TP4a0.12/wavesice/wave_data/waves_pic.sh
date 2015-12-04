#!/bin/bash
# Sending a mail with the content of check_wamnsea.py
# $1 = YYYYMMDD

source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -eq 0 ]
then
   echo Usage:
   echo "./waves_pic.sh [date]"
   echo "date in YYYYMMDD format"
   exit
fi

#######################################################################################################
# EMAIL ADRESS
email=$(cat $FCemail)
#######################################################################################################

# wdir=$SWARP_ROUTINES/forecast_scripts/alert_waves/out
wdir=$RTmods/check_wamnsea/$1
echo "waves_pic.sh called by check_wamnsea.py:" >  tmp.txt
echo " "                                        >> tmp.txt
echo "No large waves close to ice"              >> tmp.txt
# mutt -s "WAM forecast for $1" -a tmp.txt -a $wdir/img/*.png -- $email < /dev/null
mutt -s "WAM forecast for $1" -a $wdir/img/*.png -- $email < tmp.txt
rm tmp.txt
