#!/bin/bash
# each day this script gets yesterday's OSISAF concentration from myocean

cd /work/timill/DMI_tmp
rm -r *

# download .imz files
ftpdir=ftp://ftp.dmi.dk/icemar/wa/
wget -r -nH -e robots=off --no-parent --cut-dirs=5 -R "index.html*" $ftpdir
lst=*.imz

mkdir -p tmp
cd tmp

# unpack and cp to storage directory
DMIdir=/work/shared/nersc/msc/DMI_icecharts/sigrid
for g in ../*.imz
do
   f=`basename $g`
   echo $f
   cyear=${f:0:4}
   unzip ../$f
   echo mv * $DMIdir/$cyear
   mv * $DMIdir/$cyear
done

# #########################################################################
# # Launch validation script
# # - compares today's observation to relevant forecasts
# $Vdir/ice_edge_OSISAF_1obs.sh $ydate
# #########################################################################
