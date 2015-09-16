# script to restore forecast Build directories
# if it has been cleaned

if [ $# -eq 0 ]
then
   echo Usage:
   echo restore_build.sh Xno
   echo "where Xno = 1 (ice_only), 2 (wavesice)"
   exit
else
   X=01.$1
fi

echo " "
cd $TP4_REALTIME/Build_V2.2.12_X$X
pwd

echo cp svn_Build/* .
cp svn_Build/* .

echo ./setup_patch.sh 133
./setuppatch.sh 133

if [ $1 -eq 1 ]
then
   echo cp $SWARP_ROUTINES/forecast_scripts/inputs/ice_only/flags .
   cp $SWARP_ROUTINES/forecast_scripts/inputs/ice_only/flags .
else
   echo cp $SWARP_ROUTINES/forecast_scripts/inputs/wavesice/flags .
   cp $SWARP_ROUTINES/forecast_scripts/inputs/wavesice/flags .
fi

echo " "
