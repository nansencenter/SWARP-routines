#!/bin/bash
# each day this script gets yesterday's AMSR2 concentration
# - last update is around 2000, so we run at 2100

# get yesterday's date
ydate=`date -d "yesterday" '+%Y%m%d'`

# check the last 5 days
for day in `seq 0 5`
do
   cdate=`date -d "$ydate -${day}days" +%Y%m%d`
   cyear=`date -d "$cdate" '+%Y'`
   cmon=`date -d "$cdate" '+%m'`
   cday=`date -d "$cdate" '+%d'`
   gfil=Arc_${cdate}_res3.125_pyres.nc
   zfil=$gfil.gz
   ftpdir="ftp://ftp-projects.zmaw.de/seaice/AMSR2/3.125km/"

   # where the conc files are stored
   hex_dir=/work/shared/nersc/msc/AMSR2_3125
   cd $hex_dir
   mkdir -p Arc_$cyear
   cd Arc_$cyear

   # if file not present, download
   if [ ! -f $gfil ]
   then
      mkdir -p tmp
      cd tmp
      wget $ftpdir/$zfil
      gunzip $zfil
      mv *.nc ..
      cd ..
      rm -r tmp
   fi
done
