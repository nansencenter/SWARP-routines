#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

# AS LONG AS THE OTHER SCRIPTS WORK I DON'T SEE HOW THIS SIMPLE ONE COULD FAIL.
#TODO FOR FUTURE ALERTS?

source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

tday=$2
tday_long=`date --date=$tday +%Y-%m-%d`
cyear=${tday:0:4}
echo "Collecting data produced in date $tday_long"


cd $THISFC2/$tday
tfil=${FC_OUTPUT}_$tday.tar.gz

echo "The archive file name will be "
echo " $tfil "
touch $tfil

list="bin info final_output"
list="$list netcdf"
if [ -d figures ]
then
   list="$list figures"
fi
tar -zcvf $tfil $list
chmod 777 $tfil

# TODO change to scp & Norstore
# ncftp protocol?
if [ 1 -eq 1 ]
then
   # migrate
   mkdir -p $THISFC3/$cyear
   mv $tfil $THISFC3/$cyear/
   echo "SWARP products of $tday"
   echo "stored in $THISFC3/$cyear"
else
   # Norstore
   source $hidden/ssh_info_Nor.src
   dest=/scratch/$huser/SWARP_FC
   scp -i $keyname $tfil $huser@login3.norstore.uio.no:$dest
fi
