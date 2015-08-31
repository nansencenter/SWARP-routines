#!/bin/bash
# Sending a mail with the content of check_wamnsea.py
# $1 = YYYYMMDD

source $SWARP_ROUTINES/source_files/hex_vars.src

# EMAIL ADRESS
#######################################################################################################
fc_email=$SWARP_ROUTINES/forecast_scripts/fc_alert_email.txt
email=$(cat $fc_email)
#######################################################################################################

# wdir=$SWARP_ROUTINES/forecast_scripts/alert_waves/out
wdir=/work/timill/RealTime_Models/check_wamnsea/$1
mutt -s "Wave Breaking events for $1" -a $wdir/lst/$1_list.txt -a $wdir/img/$1.png -- $email < /dev/null
