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
fcday=`date --date="$1 +2days" "+%Y%m%d"`

#######################################################################################################
# EMAIL ADRESS
email=$(cat $FCemail)
if [ 1 -eq 1 ]
then
   # add Giacomo for his interest
   email=$email,gcmdnt90@gmail.com
fi
#######################################################################################################

wdir=$RTmods/check_wamnsea/$1
echo "waves_alert.sh called by check_wamnsea.py:"  >  tmp.txt
echo " "                                           >> tmp.txt
cat $wdir/lst/*.txt                                >> tmp.txt

if [  -f $wdir/lst/$fcday*.txt ]
then
   cat $wdir/lst/$fcday*.txt                                                                       >> tmp.txt
   mutt -s "WAM forecast for $fcday" -a $wdir/img/$fcday*.png -a $wdir/lst/$fcday*.txt -- $email  <  tmp.txt
else
   echo "No large waves close to ice"                                                              >> tmp.txt
   mutt -s "WAM forecast for $fcday" -a $wdir/img/$fcday*.png -- $email                           <  tmp.txt
fi

rm tmp.txt
