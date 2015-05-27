#!/bin/bash
# Sending a mail with the content of check_wamnsea.py
# $1 = YYYYMMDD

# EMAIL ADRESS
#######################################################################################################
fc_email=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
email=$(cat $fc_email)
#######################################################################################################

wdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/alert_waves/out

mutt -s "Wave Breaking events for $1" -a $wdir/lst/$1_list.txt -a $wdir/lst/$1_threshold_list.txt -a $wdir/img/$1.png -a $wdir/img/$1_threshold.png -- $email < /dev/null
