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

# EMAIL ADRESS
#######################################################################################################
fc_email=$SWARP_ROUTINES/forecast_scripts/fc_alert_email.txt
email=$(cat $fc_email)
#######################################################################################################

wdir=/work/timill/RealTime_Models/check_ww3arctic/$1
echo "waves_alert.sh called by check_ww3arctic.py:"   >  tmp.txt
echo " "                                              >> tmp.txt
cat $wdir/lst/*.txt                                   >> tmp.txt

# send email
mutt -s "WW3a forecast for $1" -a $wdir/img/*.png -a $wdir/lst/*.txt -- $email < tmp.txt
rm tmp.txt
