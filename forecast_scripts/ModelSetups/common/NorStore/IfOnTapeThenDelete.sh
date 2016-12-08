#!/bin/bash
# check weekly (day in week?) if week .tar.gz files are written to Tape and size is ok
# get current date define timespan to check for files!
# Read list of files in WriteToTape_list/yyyymmdd_week_filelist.txt
# Check if files are written to tape
# then clear all project/../daily and scratch/.../weekly files

# INPUT
# 1) HYCOMreg read in expt HYCOM region, to sort and get right paths 
# 2) forecast version, same as folder name, tex "ice_only"

proj=NS2993K
thisdate=`date`
me=`readlink -f $0`
here=`dirname $me`

if [ $# -ne 2 ]
then
   echo Usage: "`basename $me` [region] [FCtype]"
   echo region: TP4, BS1, FR1
   echo FCtype: ice_only, wavesice_ww3arctic
   exit
fi

echo "******************************************"
echo " Running IfOnTapeThenDelete.sh $1 $2"
echo " Time: $thisdate"
echo "******************************************"

modset="$SWARP_ROUTINES/forecast_scripts/ModelSetups/"
HYCOMreg=$1
FCtype=$2

# Set up paths 
if [ $HYCOMreg == TP4 ]
then
   hycomreg=TP4a0.12
elif [ $HYCOMreg == FR1 ]
then
   hycomreg=FR1a0.03
elif [ $HYCOMreg == BS1 ]
then
   hycomreg=BS1a0.045
fi

source $modset/$hycomreg/$FCtype/inputs/THISFC.src
path1=/projects/$proj/SWARP_FC/results/$hycomreg/$FCtype
path2=/scratch/timill/SWARP_FC/results/$hycomreg/$FCtype
path3=SWARP_FC/results/$hycomreg/$FCtype

mkdir -p $path1/CheckedWriteToTape_list
Cdir=$path1/CheckedWriteToTape_list
Wdir=$path1/WriteToTape_list
# n_list=$path1/NOT_WriteToTape_list #not used

#cdate=`date -d "$cdate" +"%Y%m%d"`
#year=`date -d "$cdate" +"%Y"`
#date1=`date +"%Y%m%d" -d "$cdate-6 days"`
#
#weekdir=${date1}_week

# Loop over files in WriteToTape_list/
# check if file is written to tape
# if not, send error mess (e-mail?)
# if true, move file from $Wdir to $Cdir
# example filename 20160118_week_filelist.txt

# Files written to tape:
Tdir=/tape/NS2993K/SWARP_FC/results/$hycomreg/$FCtype/
t_files=(`lst $Tdir`)
if [ 1 -eq 0 ]
then
   echo "files on tape:"
   for f in ${t_files[@]}
   do
      lst -lh $Tdir/$f
   done
   exit
fi

for w_file in "$Wdir"/*.txt
do
   echo ""
   echo "Weekly file to check: $w_file"
   echo "date of weekly file: $w_date"
   fname=`basename $w_file`
   w_date=${fname:0:8}
   
   # Need to compare date in files
   found=0
   for t_file in ${t_files[@]}
   do
      tfil=`basename $t_file`
      tdate=${tfil:0:8}
      if [ $tdate == $w_date ]
      then
         info=(`lst -l $Tdir/$t_file`)
         echo "found $w_date"
         if [ ${info[4]} -gt 0 ]
         then
            echo file not empty
            found=1
         else
            echo file empty
         fi
         break
      fi
   done


   if [ $found -eq 0 ]
   then
      # don't delete files
      echo file not on tape
      echo - continuing without deleting
      continue
   else
      echo file on tape
      echo deleting from disk
   fi


   # remove daily files:
   d_list=(`cat $w_file`)
   for df in ${d_list[@]}
   do
      echo rm -f $path1/$df
      rm -f $path1/$df
   done

   # move the weekly file list to Checked...
   echo mv $w_file $Cdir
   mv $w_file $Cdir

done
