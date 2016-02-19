#!/bin/bash
# each day this script gets yesterday's AMSR2 concentration
# - last update is around 2000, so we run at 2100

# get yesterday's date
yyear=`date -d "yesterday" '+%Y'`
ymon=`date -d "yesterday" '+%m'`
yday=`date -d "yesterday" '+%d'`
ydate=$yyear$ymon$yday
gfil=Arc_${ydate}_res3.125_pyres.nc.gz
ftpdir="ftp://ftp-projects.zmaw.de/seaice/AMSR2/3.125km/"

# where the conc files are stored
hex_dir=/work/shared/nersc/msc/AMSR2_3125
cd $hex_dir
mkdir -p Arc_$yyear
cd Arc_$yyear

mkdir -p tmp
cd tmp
wget $ftpdir/$gfil
gunzip $gfil
mv *.nc ..
cd ..
rm -r tmp
