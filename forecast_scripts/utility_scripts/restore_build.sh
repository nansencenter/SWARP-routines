# script to restore forecast Build directories
# if it has been cleaned

if [ $# -eq 0 ]
then
   echo Usage:
   echo restore_build.sh Xno
   echo "where Xno = 1 (ice_only), 2 (waves-in-ice - WAM), 3 (waves-in-ice - WW3)"
   exit
else
   X=01.$1
fi

echo " "
cd $TP4_REALTIME/Build_V2.2.12_X$X
pwd

if [ ! -f ../expt_$X/blkdat.input ]
then
   # need some files from expt directory
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_expt.sh $1
fi

echo cp svn_Build/* .
cp svn_Build/* .

echo ./setup_patch.sh 133
./setuppatch.sh 133

if [ $1 -eq 1 ]
then
   THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only
elif [ $1 -eq 2 ]
then
   THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice
elif [ $1 -eq 3 ]
then
   THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_ww3arctic
fi

echo "cp $THISFC/inputs/flags ."
echo " "
cp $THISFC/inputs/flags .
