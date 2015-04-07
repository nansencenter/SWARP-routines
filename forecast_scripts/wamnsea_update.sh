#!/bin/bash
# Run from cron job regularely
# Updating the WAMNSEA data in /work/shared/nersc/msc/WAMNSEA/
# Cron settings can be changed by editing 'mycronfile' and updating it
# using 'crontab mycronfile'.
# ==============================================================================
# 1. Update wam_nsea_fc_YYYY.nc file with forecast data from met.no WAMNSEA 10km 
# ==============================================================================
MAINDIR=/home/nersc/bergh/Realtime
SHOST=/Wave/METNO-WAM-NORDIC/
LDIR=/work/shared/nersc/msc/WAMNSEA

thedate=`date +"%Y%m%d"`
#thedate='20150319'
# check year today, tomorrow, the day after tomorrow
year=`date +"%Y"`
dayname=`date +%A`
yearplus1=`expr $year + 1`
year1d=$(date --date="1 day" +"%Y")
year2d=$(date --date="2 day" +"%Y")
mkdir -p $LDIR/${year}
mkdir -p $LDIR/${yearplus1}

echo '** Update WAMNSEA wave data from myocean.met.no **'
echo "Check if year change next 2 days:"
echo " Today year is:.........$year"
echo " Next year is:..........$yearplus1"
echo " Today +1 day year is:..$year1d"
echo " Today +2 day year is:..$year2d"
echo ""

cd $LDIR

# Loop, try download WAMNSEA file, if not exist try again after 30 min
# after 4.5h exit and write Error message
count=0
while  ! [ -f wam_nsea.an.${thedate}.nc ]
do

  cat > ncftp.in<<EOF
open myocean
cd $SHOST
get wam_nsea.an.${thedate}.nc
get wam_nsea.fc.${thedate}.nc
set confirm-close no
bye
EOF
   ncftp < ncftp.in
   rm ncftp.in

 if  ! [ -f wam_nsea.an.${thedate}.nc ]
  then
   count=`expr $count + 1`
   echo " The file wam_nsea.an.${thedate}.nc doesnt exist"
   echo ""
    if [  $count -ge 9 ]
     then 
      echo " *************************************************** "
      echo " ***** Download WAMNSEA 10km data fails - exit ***** "
      echo " *************************************************** "
      exit
    fi
   echo " $count :count - We wait 30 min, try to download agin"
   echo ""
   sleep 1800
 else
   echo "The file wam_nsea.an.${thedate}.nc exist - continue"
   echo ""
 fi

done

# end of download while loop

# load cdo on hexagon with direct path
# module load cdo
# /work/apps/cdo/1.5.4-cray/bin/cdo

# check year, save in correct year file

# Use cdo functions to add last +1 day (.an. file)
#  to "keep" file (wam_nsea_2012.nc)
# Remove and recreate forecast file (wam_nsea_fc_${year}.nc)
#  that also include the +1 and +2 forecast days

echo " Add data for this day (.an. file) at end of year file"
echo " Todays inputfile: wam_nsea.an.${thedate}.nc"
echo " Cat into year file: wam_nsea_${year}.nc"
echo " Then rm wam_nsea_fc_${year}.nc and "
echo " cdo copy wam_nsea_${year}.nc wam_nsea.fc.${thedate}.nc "
echo "  into new version of wam_nsea_fc_${year}.nc"
echo ""

# Check if .an. file exist, if not exit
if [ -f wam_nsea.an.${thedate}.nc ]
 then
  /work/apps/cdo/1.5.4-cray/bin/cdo cat wam_nsea.an.${thedate}.nc $LDIR/${year}/wam_nsea_${year}.nc
else
 echo "PROBLEM! wam_nsea.an.${thedate}.nc doesn't exist, no update for this day, exit"
 exit
fi

cd $LDIR/${year}

echo "Remove forecast year file wam_nsea_fc_${year}.nc if exist"
echo ""
if [ -f wam_nsea_fc_${year}.nc ]
 then
  rm wam_nsea_fc_${year}.nc
fi

echo "Check for year today, today +1, today +2"
if [ $year -eq $year1d -a $year -eq $year2d ]
 then
  echo "Same year all 3 days" #merge fc in the file of the year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy wam_nsea_${year}.nc $LDIR/wam_nsea.fc.${thedate}.nc wam_nsea_fc_${year}.nc
elif [ $year -eq $year1d -a $yearplus1 -eq $year2d ]
 then
  echo "Same year today and today +1 - new year today +2" #split the fc file, day+1 in this year, day+2 in the next year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy wam_nsea_${year}.nc wam_nsea_fc_${year}.nc
  /work/apps/cdo/1.5.4-cray/bin/cdo splityear wam_nsea.fc.${thedate}.nc wam_nsea_fc_split_
  /work/apps/cdo/1.5.4-cray/bin/cdo cat wam_nsea_fc_split_${year}.nc wam_nsea_fc_${year}.nc
  /work/apps/cdo/1.5.4-cray/bin/cdo cat wam_nsea_fc_split_${yearplus1}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  rm wam_nsea_fc_split_*
elif [ $yearplus1 -eq $year1d -a  $yearplus1 -eq $year2d ]
 then
  echo "New year today +1 and today +2" #the whole fc file is in the next year
  /work/apps/cdo/1.5.4-cray/bin/cdo copy wam_nsea.fc.${thedate}.nc $LDIR/${yearplus1}/wam_nsea_fc_${yearplus1}.nc
  echo "Happy new year user - remember to check if the next upload is correct"
fi
echo ""

#check for missed files. the get function in ncftp will download only if the file DOESN'T exist already in the folder

#NOTE if the file exists but it's different this script will overwrite said file
#every DAY the script will check and merge every file missed in the past week

cd $LDIR
if [ $dayname == $dayname ]
then
	cat > ncftp.in<<EOF
open myocean
cd $SHOST
get -f wam_nsea.*.nc
set confirm-close no
bye
EOF
	ncftp < ncftp.in
	rm ncftp.in
	/work/apps/cdo/1.5.4-cray/bin/cdo -O mergetime *.an.${year}*.nc ./${year}/wam_nsea_${year}.nc
fi

cd $MAINDIR 
