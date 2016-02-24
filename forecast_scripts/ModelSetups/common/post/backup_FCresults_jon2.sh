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
datelist=$THISFC/logs/datelist.txt

tday=$2
#tday=`date +%Y%m%d`
tday_long=`date --date=$tday +%Y-%m-%d`
tnr=`date --date=$tday +%u`
cyear=${tday:0:4}
dm2l=`date -d "$tday - 2 days" +%Y%m%d`
dm8l=`date -d "$tday - 8 days" +%Y%m%d`
# TOPAZ year day start at 000 (so -1 extra day)
TYDm2=`date -d "$tday - 3 days" +%j`
TYDm8=`date -d "$tday - 9 days" +%j`
#echo "tday $tday, tday_long $tday_long, tnr $tnr"
#echo "dm2 $TYDm2 $dm2l, dm8 $TYDm8 $dm8l"
#exit

echo "Collecting data produced in date $tday_long"

cd $THISFC2/$tday
tfil=${FC_OUTPUT}_$tday.tar.gz

#echo "The archive file name will be "
#echo " $tfil "
#touch $tfil

list="bin info final_output"
list="$list netcdf"
if [ -d figures ]
then
   list="$list figures"
fi
##tar -zcvf $tfil $list
##chmod 777 $tfil
echo "tar -zcvf $tfil $list"
echo "chmod 777 $tfil"

# TODO move this to gather_FCresults or separate script?
# check if there is a Wednesday then tar and backup restart files on Norstore 
echo "nesting_inner= $nesting_inner, FCtype= $FCtype" 
if [[ $tnr -eq 3 ]] && [[ $nesting_inner == "T" ]] && [[ $FCtype == "ice_only" ]]; then
   # 1) Restart from Monday!
   echo "It is Wednesday and ice_only, so backup monday TYDm2 (-2 day) Topaz restart file on Norstore" 
   resfilebase=${rungen}restart${cyear}_$TYDm2
   resfile=${resfilebase}.tar.gz
   echo "Save all restart into: $resfile"
   resdir=$RTmods/$HYCOMreg/expt_01.$Xno/data
   #  TP4restart2016_017_00.a ...
   tar zcvfh $resfile $resdir/${resfilebase}_00.a $resdir/${resfilebase}_00.b $resdir/${resfilebase}_00ICE.uf
   chmod 777 $resfile
   echo "tar zcvfh $resfile $resdir/${resfilebase}_00.a $resdir/${resfilebase}_00.b $resdir/${resfilebase}_00ICE.uf"
   echo "chmod 777 $resfile"
   # 2) then save nesting files from prevous Monday (-9) to last Sun (-3) into separate tar
   nest=$RTmods/$HYCOMreg/expt_01.$Xno/SCRATCH/nest
   nestout=FC_${HYCOMreg}_nest7days_${TYDm8}.tar.gz
   nestlist=()
   for nr in $(seq -w $TYDm8 $TYDm2)
   do
      echo "from $TYDm8 to $TYDm2, now nr: $nr"
      addlist=()
      addlist=$(ls ${nest}/nest_${cyear}_${nr}*)
      echo "addlist: ${addlist[*]})"
      if [ -n "$addlist" ]; then
         nestlist+=("${addlist[*]}")
      else
         echo "!! Nesting backup, no files for topaz day $nr exists, though continue !!"
      fi
   done
   if [ -n "$nestlist" ]; then
      tar zcvfh $nestout  ${nestlist[*]} 
   #echo "NESTLIST ALL: ${nestlist[*]}"
   else
      echo "!! No nesting files at all for $TYDm9 to $TYDm3 !!" 
   fi
else
   echo "IS NOT TRUE EXIT"
   exit
fi

## test above

echo "files are saved in: $THISFC2/$tday"

exit

# TODO change to only scp to Norstore
#if [ 1 -eq 1 ]
# then
# Norstore
   # Results
   dest1=/projects/NS2993K/SWARP_FC/results/$HYCOMreg/$FCtype
   dest=/norstore_osl/home/timill/SWARP-routines/forecast_scripts/ModelSetups/common/NorStore
   #  make sure that the results directories /results .../results/$HYCOMreg and .../results/$HYCOMreg/$FCtype exists?
   command1=`ssh -i $keyname $huser@login3.norstore.uio.no "sh $dest/swarp_fc_mkdir.sh $HYCOMreg $FCtype results"`
   echo "$command1"
   # scp file to /project/... folder $dest
   echo "scp $tfil to $huser@login3.norstore.uio.no:$dest1"
   scp -i $keyname $tfil $huser@login3.norstore.uio.no:$dest1
   # change permission
   command2=`ssh -i $keyname $huser@login3.norstore.uio.no "chmod 644 $dest1/$tfil"`
   echo "$command2"
   # Restart
   if [ -s $resfile ]; then
      dest3=/projects/NS2993K/SWARP_FC/restart/$HYCOMreg/$FCtype
      command3=`ssh -i $keyname $huser@login3.norstore.uio.no "sh $dest/swarp_fc_mkdir.sh $HYCOMreg $FCtype restarts"`
      echo "$command3"
      echo "scp $tfil to $huser@login3.norstore.uio.no:$dest3"
    scp -i $keyname $resfile $huser@login3.norstore.uio.no:$dest3
   fi
   # Nesting
   if [ -s $nestout ]; then
      dest4=/projects/NS2993K/SWARP_FC/nesting/$HYCOMreg/$FCtype
      command4=`ssh -i $keyname $huser@login3.norstore.uio.no "sh $dest/swarp_fc_mkdir.sh $HYCOMreg $FCtype nesting"`
      echo "$command4"
      echo "scp $nestout to $huser@login3.norstore.uio.no:$dest4"
      scp -i $keyname $nestout $huser@login3.norstore.uio.no:$dest4
   fi
   # else
   # migrate
   mkdir -p $THISFC3/$cyear
   mv $tfil $THISFC3/$cyear/
   echo "SWARP products of $tday"
   echo "stored in $THISFC3/$cyear"
# fi
