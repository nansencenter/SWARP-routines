#!/bin/bash
# each day this script gets yesterday's OSISAF concentration from myocean
ME=`readlink -f $0`
Vdir=`dirname $ME`

tmpdir=/work/timill/DMI_tmp
cd $tmpdir

if [ 1 -eq 1 ]
then
   # download .imz files
   rm -r *
   ftpdir=ftp://ftp.dmi.dk/icemar/wa/
   wget -r -nH -e robots=off --no-parent --cut-dirs=5 -R "index.html*" $ftpdir
fi
lst=*.imz

Tmpdir=$tmpdir/tmp
mkdir -p $Tmpdir
cd $Tmpdir

# unpack and cp to storage directory
DMIdir=/work/shared/nersc/msc/DMI_icecharts/sigrid
for g in ../*.imz
do
   pwd
   f=`basename $g`
   echo $f
   cyear=${f:0:4}
   unzip ../$f

   sroot=${f%*.imz}

   odir=$DMIdir/$cyear
   mkdir -p $odir

   Odir=$odir/$sroot
   mkdir -p $Odir
   echo mv * $Odir
   mv * $Odir

   #########################################################################
   # Launch script to extract MIZ & calc MIZ width/area
   # - compares today's observation to relevant forecasts
   $Vdir/DMI_process.sh $cyear
   #########################################################################

   #clean (only needed if testing)
   rm -f $Tmpdir/*

done
