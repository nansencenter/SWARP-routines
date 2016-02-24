#!/bin/bash
# This script will collect and archive the whole directory with the final products of the SWARP model

# AS LONG AS THE OTHER SCRIPTS WORK I DON'T SEE HOW THIS SIMPLE ONE COULD FAIL.
#TODO FOR FUTURE ALERTS?

# TODO get this? user on NorStore
huser=timill
keyname=/home/nersc/timill/.ssh/key_hex2store

source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC

#tday=$2
tday=`date +%Y%m%d`
tday_long=`date --date=$tday +%Y-%m-%d`
tnr=`date --date=$tday +%u`
cyear=${tday:0:4}
jday0=10#$(date --date=$tday +%j)
jday_today0=$(( $jday0 - 1))              # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0) # 3 digits

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

# check if there is a monday then tar and backup on Norstore 
# 1) check if monday and if nesting then in $THIS_SRC ?
# if [[ $tnr == 1 ]]; then
#  echo "It is monday so save restart and backup on Norstore" 
#  TP4restart2016_017_00.a
#  tar ${rungen}restart${cyear}_YD.tar.gz ${rungen}restart${cyear}_$TYD_00.a ${rungen}restart${cyear}_$TYD_00.b
#  ${rungen}restart${cyear}_$TYD_00ICE.uf 
# fi
# 2) check if it is TP4 and if nesting then save nesting files into separate tar

# TODO change to only scp to Norstore
#if [ 1 -eq 1 ]
# then
# Norstore
   dest=/projects/NS2993K/SWARP_FC/results/$HYCOMreg/$FCtype
   dest2=/norstore_osl/home/timill/SWARP-routines/forecast_scripts/ModelSetups/common/NorStore
   #  make sure that the "results" directories .../results .../results/$HYCOMreg and .../results/$HYCOMreg/$FCtype exist?
   command1=`ssh -i $keyname $huser@login3.norstore.uio.no "sh $dest2/swarp_fc_mkdir.sh $HYCOMreg $FCtypei results"`
   echo "$command1"
   # Move file to /project/... folder $dest
   echo "scp $tfil to $huser@login3.norstore.uio.no:$dest"
   scp -i $keyname $tfil $huser@login3.norstore.uio.no:$dest
   # change permission
   command2=`ssh -i $keyname $huser@login3.norstore.uio.no "chmod 644 $dest/$tfil"`
   echo "$command2"
# else
   # migrate
   mkdir -p $THISFC3/$cyear
   mv $tfil $THISFC3/$cyear/
   echo "SWARP products of $tday"
   echo "stored in $THISFC3/$cyear"
# fi
