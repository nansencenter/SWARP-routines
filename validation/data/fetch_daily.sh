#!/bin/bash
# FETCH DAILY BINARY FILES

# IF NOT INPUT TAKE ICE_ONLY ELSE WAVESICE
#echo "Type the ""ice"" for ICE_ONLY or ""waves"" for WAVESICE followed by [ENTER]:  "
#read typo
#echo "Type the date of the run (YYYYMMDD), followed by [ENTER]:  "
#read fdate

fdate=$1
typo=$3
jdayp1=$(date -d "$fdate" +%j)
jdayp=$(expr $jdayp1 - 1)
steve=$(($jdayp-$2))
jday=`printf '%3.3d' $steve`
steve=$(($jday))
#echo "Correspond to the following hycom julian day:   $jday"
#echo ""

if [ $typo == "ice" ]
then
   daily=TP4DAILY*
   wdir=$ICE_RESULTS/$fdate/bin
    for file in $wdir/$daily
   do
      echo ${file##*/}
   done
   #echo "Enter last 3 numbers (hycom julian day):  "
   #read steve
   thedate=$(date -d "`date +%Y`-01-01 +$(( ${steve} ))days" +%Y%m%d)
   #echo "You chose to copy products from: $thedate"
   #echo "[ENTER] to confirm"
   #read ok
   #echo "please wait..."
	 rm ./*.nc
   cp $wdir/${daily}*${steve}* ./
elif [ $typo == "waves" ]
then
   daily=TP4DAILY*
   wdir=$WAVES_RESULTS/$fdate/bin
   for file in $wdir/$daily
   do
      echo ${file##*/}
   done
   #echo "Enter last 3 numbers (hycom julian day):  "
   #read steve
   thedate=$(date -d "`date +%Y`-01-01 +$(( ${steve} ))days" +%Y%m%d)
   #echo "You chose to copy products from: $thedate"
   #echo "[ENTER] to confirm"
   #read ok
   #echo "please wait..."
	 rm ./*.nc
   cp $wdir/${daily}_*_*_*${steve}* ./
else
   echo "Please enter either ""ice"" or "" waves"" "
   exit
fi


# MAKE .nc FILE
h2p_in="$GIT_REPOS/SWARP-routines/netcdf_production/Input"
ddir=`pwd`

mkdir -p tmp 
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

#echo "********************************************************"
#echo "Running hyc2proj on $g..."
#echo "********************************************************"
#echo " "

hyc2proj $g
cp *.nc $ddir
cd $ddir
rm -r tmp

# REMOVE BINARY FILES (TOO BIG FOR THIS DISK)
rm *.a *.b

# FETCH OSI-SAF
year=${fdate:0:4}
osidir=/work/shared/nersc/msc/OSI-SAF/${year}_nh_polstere
osifil=*${thedate}*.nc
cp $osidir/$osifil ./
