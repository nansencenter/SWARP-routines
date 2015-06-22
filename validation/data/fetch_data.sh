#!/bin/bash
# FETCH DAILY BINARY FILES

# IF NOT INPUT TAKE ICE_ONLY ELSE WAVESICE
echo "Type the ""ice"" for ICE_ONLY or ""waves"" for WAVESICE followed by [ENTER]:  "
read typo
echo "Type the date of the run (YYYYMMDD), followed by [ENTER]:  "
read fdate
ryear=${fdate:0:4}
rmonth=${fdate:4:2}
rday=${fdate:6:2}
jdayp1=$(date -d "$fdate" +%j)
jdayp=$(expr $jdayp1 - 1)
jday=`printf '%3.3d' $jdayp`
echo "Correspond to the following hycom julian day:   $jday"
echo ""

bdir=`pwd`

mkdir -p tmp

REPOCHECK=0

if [ $typo == "ice" ]
then
   bak=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/ice_only/${ryear}/*${fdate}.tar.gz
   edir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work
   daily=TP4DAILY*
   wdir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$fdate/bin
   for file in $wdir/$daily
   do
      if [ -z "`find $wdir -type f`" ] 
      then
         echo 'WARNING - No files detected, [1] to check the repositories, [Enter] to exit'
         echo ''
         read REPOCHECK
      fi
      echo ${file##*/}
   done
   if [ -z $REPOCHECK ]
   then
      exit
   elif [ $REPOCHECK == 1 ]
   then
      cp -v $bak $edir/
      cd $edir
      tar -zxvf *${fdate}.tar.gz
      rm *.tar.gz
      cd $bdir 
      for file in $wdir/$daily
      do
         echo ${file##*/}
      done
   fi
   echo "Enter hycom julian day:  "
   read steve
   thedate=$(date -d "`date +%Y`-01-01 +$(( ${steve} ))days" +%Y%m%d)
   echo "You chose to copy products from: $thedate"
   echo "[ENTER] to confirm"
   read ok
   echo "please wait..."
   old=*.nc
   rm $old
   cp $wdir/${daily}*${steve}* ./tmp/
elif [ $typo == "waves" ]
then
   bak=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/wavesice/${ryear}/*${fdate}.tar.gz
   edir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work
   daily=TP4archv*
   wdir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$fdate/bin
   for file in $wdir/$daily
   do
      if [ -z "`find $wdir -type f`" ] 
      then
         echo 'WARNING - No files detected, [1] to check the repositories, [Enter] to exit'
         echo ''
         read REPOCHECK
      fi
      echo ${file##*/}
   done
   if [ -z $REPOCHECK ]
   then
      exit
   elif [ $REPOCHECK == 1 ]
   then
      cp -v $bak $edir/
      cd $edir
      tar -zxvf *${fdate}.tar.gz
      rm *.tar.gz
      cd $bdir 
      for file in $wdir/$daily
      do
         echo ${file##*/}
      done
   fi
   echo "Enter hycom julian day:  "
   read steve
   thedate=$(date -d "`date +%Y`-01-01 +$(( ${steve} ))days" +%Y%m%d)
   echo "You chose to copy products from: $thedate"
   echo "[ENTER] to confirm"
   read ok
   echo "please wait..."
   old=*.nc
   rm $old   
   cp $wdir/${daily}_*_*_*${steve}* ./tmp/
else
   echo "Please enter either ""ice"" or "" waves"" "
   exit
fi


# MAKE .nc FILE
h2p_in="$GIT_REPOS/SWARP-routines/netcdf_production/Input"
ddir=`pwd`

cd tmp 

# Info for hyc2proj
ln -s $h2p_in/proj.in .
ln -s $h2p_in/extract.archv . 
ln -s $h2p_in/extract.daily . 
ln -s $h2p_in/depthlevels.in .

# Info about grid from topo dir
tp4_input=/work/shared/nersc/msc/ModelInput/TOPAZ4/TP4a0.12/topo
ln -s $tp4_input/grid.info .
ln -s $tp4_input/regional.grid.a .
ln -s $tp4_input/regional.grid.b .
ln -s $tp4_input/regional.depth.a .
ln -s $tp4_input/regional.depth.b .

# CONVERT
f=$ddir/TP4*.a
len=${#ddir}
g=${f:$((len+1))}
len=${#g}
base=${g:0:$((len-2))}
ln -s $f $ddir/$base.b .

echo "********************************************************"
echo "Running hyc2proj on $g..."
echo "********************************************************"
echo " "

hyc2proj $g
cp *.nc $ddir
cd $ddir
rm -r tmp

# REMOVE BINARY FILES (TOO BIG FOR THIS DISK)
#rm *.a *.b

# FETCH OSI-SAF
year=${fdate:0:4}
osidir=/work/shared/nersc/msc/OSI-SAF/${year}_nh_polstere
osifil=*${thedate}*.nc
cp $osidir/$osifil ./
