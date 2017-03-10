#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#18-03-15 VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder
me=`readlink -f $0`
here=`dirname $me`

source $SWARP_ROUTINES/source_files/hex_vars.src
mkdir -p /work/$USER/tmp
cd       /work/$USER/tmp

# =========================================================================================================
# EMAIL ADDRESS FOR THE WEEKLY UPDATE
email=$(cat $FCemail)
# =========================================================================================================

print_info=1

# DIRECTORIES AND TIME DEFINITION
rdir=/work/fanf/TOPAZ_RT 
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
fcldir=$FCcommon/logs
mkdir -p $fcldir

cyear=$(date +%Y) 
pyear=$(expr $cyear - 1)
tday=$(date +%Y%m%d-%A)
cday=$(date +%Y%m%d)
tdd=$(date +%d)
jday=10#$(date +%j)
if [ "$(date +%A)" == "Monday" ]
then
   #latest restart is today
   JDAY=$((jday-1))
   JDAY=`printf %3.3d $JDAY`
   Rlatest=TP4restart${cyear}_${JDAY}_00
   dt=$cday
   dty=$cyear
   dtj=$jday
else
   #latest restart is last Monday
   dt=`date --date="last Monday" "+%Y%m%d"`
   dty=`date --date="last Monday" "+%Y"`
   dtj=10#`date --date="last Monday" "+%j"`
   JDAY=$((dtj-1))
   JDAY=`printf %3.3d $JDAY`
   Rlatest=TP4restart${dty}_${JDAY}_00
fi

$here/Retrieve_TOPAZ_RT_restart.sh $Rlatest
