#!/bin/bash
# Backup SWARP forecast from /project/ to /tape area ones a week
# Run on Sundays, check all dates for the week, if not send error message
# Untar all files and make a new "weekly" tar file
# ASSUME input file type as:
#  SWARP${FC_OUTPUT}_${idate}.tar.gz (SWARPiceonly_forecast_20150101.tar.gz)
# Use Norstore backup commands:
# https://www.norstore.no/services/tape-storage
#  WriteToTape <filename> <projectname> [<relative_path_on_tape>] [--replicate]
# INPUT
# 1) HYCOMreg read in expt HYCOM region, to sort and get right paths 
# 2) forecast version, same as folder name, tex "ice_only"
# 3) fictive current date set to YYYYMMDD  (if not def. the current date is assumed)
# OUTPUT
# save files to NorStore tape, file name of FIRST DATE IN WEEK!!
# outfile=YYYYMMDD_weekly.tar.gz
# write weekly text file in $path1/WriteToTape_list, list of existing daily files
# write weekly text file in $path1/NOT_WriteToTape_list, list of NOT existing daily files
# TODO do we need a replicate of the backup on a different tape??

thisdate=`date`
echo "******************************************"
echo " Running: swarp_fc_tape_backup.sh $1 $2 $3"
echo " Time: $thisdate"
echo "******************************************"
proj=NS2993K
modset="$SWARP_ROUTINES/forecast_scripts/ModelSetups/"
if [ -z "$1" ]; then
    errormessage="EXIT NEED at least input 1) and 2) about HYCOMreg & FCtype"
    echo $errormessage 1>&2
    echo "$errormessage"
    exit
fi
HYCOMreg=$1
if [ -z "$2" ]; then
   errormessage="EXIT NEED at least input 1) and 2) about HYCOMreg & FCtype" 
   echo $errormessage 1>&2
   echo "$errormessage"
   exit
fi
FCtype=$2
# get FC_OUTPUT from THISFC.src
source $modset/$HYCOMreg/$FCtype/inputs/THISFC.src
if [ -z "$3" ]; then
   cdate=`date`
else
   cdate=$3
fi
# Set up paths to 
path1=/projects/$proj/SWARP_FC/results/$HYCOMreg/$FCtype
path0=/scratch/timill/SWARP_FC/results/$HYCOMreg
path2=/scratch/timill/SWARP_FC/results/$HYCOMreg/$FCtype
path3=SWARP_FC/results/$HYCOMreg/$FCtype
if [ ! -d $path1 ]; then
 errormessage=" The specified folder do not exist, so EXIT!" 
 echo "$errormessage" 1>&2
 exit
fi
# make new direction in scratch for untar day files and tar week file
[ ! -d $path0 ] && mkdir -p $path0
[ ! -d $path2 ] && mkdir -p $path2

cdate=`date -d "$cdate" +"%Y%m%d"`
year=`date -d "$cdate" +"%Y"`
date1=`date +"%Y%m%d" -d "$cdate-6 days"`
weekdir=${date1}_week

# check for files the last week (assume Sunday, so last 7 days)
filelist=""
datelist=""
notfilelist=""
for i in {0..6}
do
   idate=`date +"%Y%m%d" -d "$cdate-$i days"`
   file=${FC_OUTPUT}_${idate}.tar.gz
   if [ -s $path1/$file ]
   then
      echo "file for $idate exist:"
      #echo "   $file"
      filelist="$filelist $file"
      datelist="$datelist $idate"
   else
      errormessage="file for $idate DOES NOT exist:"
      echo $errormessage 1>&2
      notfilelist="$notfilelist $file"
      #echo "   $file"
   fi
done

datelist=($datelist)
filelist=($filelist)
notfilelist=($notfilelist)

echo "DATELIST:"
echo "${datelist[@]}"
echo "FILELIST:"
echo "${filelist[@]}"
echo "NOTFILELIST:"
echo "${notfilelist[@]}"

# if some daily files is missing then write a list with the week name 
if [[ ${#notfilelist} -gt 0 ]]
then
   mkdir -p $path1/NOT_WriteToTape_list
   printf "%s\n" "${notfilelist[@]}" > $path1/NOT_WriteToTape_list/${weekdir}_notfilelist.txt
fi

# collect files into week folder and sub date folders, untar day files and then tar into week file 
if [[ ${#filelist} -gt 0 ]]
then
   mkdir -p $path1/WriteToTape_list
   printf "%s\n" "${filelist[@]}" > $path1/WriteToTape_list/${weekdir}_filelist.txt
   mkdir -p $path2/$weekdir
   nr=${#filelist[@]}
   ind=`expr $nr - 1`
   echo "Number of files that exist: $nr" 
   # loop over all dates/files  
   for i in $(eval echo "{0..$ind}")
   do 
      datedir=${datelist[$i]}
      file="${filelist[$i]}" 
      if [ -s "$path1/$file" ]
      then
         mkdir -p $path2/$weekdir/$datedir
         echo "untar $path1/$file to $path2/$weekdir/$datedir"
         tar -xzf $path1/$file -C $path2/$weekdir/$datedir > /dev/null
         # rm files from project/..? rm $path1/$file => done in IfOnTapeThenDelete.sh
      fi
   done
   outfile=${weekdir}.tar.gz
   echo "tar $path2/$weekdir into $path2/$outfile"
   tar -czf $path2/$outfile -C / $path2/$weekdir > /dev/null
else
   errormessage=" The filelist is emty for week: $date1 to $cdate"
   echo $errormessage 2>&1
   exit
fi

# if new week.tar.gz exist rm weekdir
if [ -f $path2/$outfile ]; then
   echo "rm folder after compressed week folder: $path2/$weekdir"
   rm -r $path2/$weekdir
   chmod 644 $path2/$outfile
fi

echo ""
echo "Collecting $FCtype daily forecasts from $path1"
echo "Compress into new week file: $path2/$outfile"
echo "Save to Tape: /tape/$proj/$path3"

# write to tape (need replica? then add --replicate)
# TODO! DO NOT remove outfile before it is saved to tape, need to check this with lst? 
yes | /norstore_osl/diverse/bin/WriteToTape $path2/$outfile $proj $path3
if [[ "$?" == "0" ]]; then
 echo "Writing to tape is succesful"
elif [ "$?" == "1" ]; then
 errormessage="ERROR - problems with WriteToTape - EXIT"
 echo "$errormessage" 2>&1
 exit
fi

## echo "**** WRITE TO TAPE ****"
## echo "yes | WriteToTape $path2/$outfile $proj $path3" 

echo "THE END - EXIT" 
